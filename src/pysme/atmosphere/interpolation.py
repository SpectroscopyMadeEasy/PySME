# -*- coding: utf-8 -*-
""" Handles reading and interpolation of atmopshere (grid) data """
import itertools
import logging

import numpy as np
from astropy import constants as const
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.interpolate import RBFInterpolator

from ..abund import Abund
from ..large_file_storage import setup_atmo
from .atmosphere import Atmosphere as Atmo
from .atmosphere import AtmosphereError, AtmosphereGrid
from .savfile import SavFile

logger = logging.getLogger(__name__)

# Radius and Surface Gravity of the Sun
# Required for spherical models
R_sun = const.R_sun.to_value("cm")
g_sun = (const.G * const.M_sun / const.R_sun ** 2).to_value("cm/s**2")
logg_sun = np.log10(g_sun)


class RbfGrid:
    """
    Precomputed data needed to interpolate model atmospheres from a grid via
    radial basis functions (RBF), as an alternative to pysme's default
    corner-model interpolation. Built once per atmosphere grid and reused for
    every subsequent RBF interpolation query against that grid.
    """

    # Relative weights used to rescale Teff/logg/[M/H] onto a shared step size before
    # nearest-neighbor search, so that neighbors are preferentially chosen close in Teff,
    # next in [M/H], and least strictly in logg when grid spacing would otherwise tie them.
    _TEFF_WEIGHT = 9
    _LOGG_WEIGHT = 11
    _MET_WEIGHT = 10

    # Number of nearest grid points used to build the local RBF interpolant.
    # WARNING: this was tuned empirically (trial and error against real atmosphere grids),
    # balancing interpolation accuracy against runtime cost. Do not change it casually:
    # too few neighbors can make the interpolated atmosphere inaccurate or discontinuous,
    # while too many can make RBFInterpolator's linear solve extremely slow (its cost grows
    # steeply with neighbor count) -- capable of hanging for hours on realistic grids.
    # Re-validate against a known grid before changing this value.
    _RBF_NEIGHBORS = 64

    def __init__(self, sav_grid, geometry):
        """
        Precompute everything needed to interpolate model atmospheres from
        ``sav_grid`` via radial basis functions: rescaled grid coordinates for
        nearest-neighbor lookup, which fields to interpolate, the common depth
        length, and the flattened array of grid values the RBF interpolant is
        fit to.

        Parameters
        ----------
        sav_grid : AtmosphereGrid
            Grid of precomputed model atmospheres to interpolate between.
        geometry : {"PP", "SPH"}
            Geometry of the atmosphere models. Only used to decide whether
            "height" and "radius" participate in the interpolation (spherical
            models only).

        Raises
        ------
        AtmosphereError
            If ``sav_grid`` is not an AtmosphereGrid.
        """
        # Set by self.initialize_gridpoints():
        self.teffs = None   # Scaled to have standard step size, used for finding neighbors in interpolation step
        self.loggs = None
        self.mets = None
        self.grid_teffs = None # Not scaled, same as in sav_grid
        self.grid_loggs = None
        self.grid_mets = None
        self.teff_step_scale = 1.0 # to convert between teffs and grid_teffs
        self.logg_step_scale = 1.0
        self.met_step_scale = 1.0
        self.gridpoints = None
        self._rbf_degree = None # Set by self.initialize_gridpoints()
        # set by self.initialize_tags():
        self.vtags = None
        self.stags = None
        # Set by self.initialize_depth():
        self.depthpoints = None
        # Set by self.read_in_grid():
        self.array_for_interpol = None # Combines all vectors and scalars for interpolation
        # Set by self.build_interpolator():
        self.interpolator = None # Cached RBFInterpolator, built once and reused across queries

        if not isinstance(sav_grid, AtmosphereGrid):
            raise AtmosphereError("Rbf initialization was requested but no sav-grid was provided")

        # Kept so interpolate_RBF can look up the nearest grid point's
        # administrative fields (opflag, wlstd, abundance pattern) that aren't
        # part of the RBF fit itself.
        self.sav_grid = sav_grid
        self.abund_format = getattr(sav_grid, "abund_format", "sme")

        self.initialize_gridpoints(sav_grid) # read in and scale teff, logg, and metallicity of the grid
        self.vtags, self.stags = self.initialize_tags(sav_grid, geometry) # Check if e.g. height and radius are used
        self.depthpoints = self.initialize_depth(sav_grid) # Set self.depthpoints
        self.read_in_grid(sav_grid) # Stores all data in self.array_for_interpol
        self.interpolator = self.build_interpolator()

    @staticmethod
    def _min_step(sorted_unique_vals, fallback=1.0):
        """
        Smallest gap between consecutive sorted unique values, or ``fallback``
        if there's only one value. A single-value axis (e.g. a spherical grid
        with only one metallicity) has no gaps to measure, but its scale
        constant is irrelevant anyway since every grid point shares that
        coordinate and it can't affect relative distances between them.
        """
        diffs = np.abs(sorted_unique_vals[1:] - sorted_unique_vals[:-1])
        return np.min(diffs) if diffs.size > 0 else fallback

    def check_within_bounds(self, teff, logg, monh, interpolation_policy="allow"):
        """
        Check a query point against the grid's Teff/logg/[M/H] range.
        ``RBFInterpolator`` (in local-neighbors mode) extrapolates with no
        warning for points outside the grid, so this is the only thing that
        stops that from happening silently.

        Mirrors the naming/semantics of the (currently unmerged)
        ``interpolation_policy`` feature on ``upstream/codex/interpolation-policy``:
        ``"allow"`` (default) logs and proceeds (extrapolating); anything else
        (e.g. ``"error"``) raises.

        Parameters
        ----------
        interpolation_policy : {"allow", "error"}, optional
            What to do about an out-of-range query (default: "allow").

        Raises
        ------
        AtmosphereError
            If the query is out of range and ``interpolation_policy`` isn't "allow".
        """
        for name, value, (lo, hi) in [
            ("Teff", teff, self.teff_bounds),
            ("log(g)", logg, self.logg_bounds),
            ("[M/H]", monh, self.met_bounds),
        ]:
            if value < lo or value > hi:
                if interpolation_policy == "allow":
                    logger.info(
                        "RBF interpolation: requested %s=%.3f outside grid range [%.3f, %.3f]; extrapolating.",
                        name, value, lo, hi,
                    )
                else:
                    raise AtmosphereError(
                        f"RBF interpolation: requested {name}={value} outside grid range [{lo}, {hi}]."
                    )

    def build_interpolator(self):
        """Build the RBFInterpolator over the whole grid, once, for reuse across queries."""
        grid_points = np.vstack((self.teffs, self.loggs, self.mets)).T
        kwargs = {} if self._rbf_degree is None else {"degree": self._rbf_degree}
        return RBFInterpolator(
            grid_points, self.array_for_interpol, neighbors=self._RBF_NEIGHBORS, **kwargs
        )

    def to_interp_space_vector(self, vtag, values, atmo):
        """
        Map a grid vector quantity (e.g. temp, rhox, height) into the space the
        RBF interpolant operates in. Everything is interpolated in log space;
        "height" is additionally combined with "radius" first (see
        ``from_interp_space_vector`` for the inverse transform).
        """
        values = np.asarray(values)
        if vtag == "height":
            return np.log10(values+atmo["radius"]) # combine height and radius for the interpolation
        return np.log10(values)

    def to_interp_space_scalar(self, stag, value):
        """Map a grid scalar quantity (e.g. teff, logg, radius) into the space the RBF interpolant operates in."""
        if stag == "radius":
            return np.log10(value)
        return value


    def initialize_gridpoints(self, sav_grid):
        """
        Read the (teff, logg, [M/H]) coordinate of every grid point in
        ``sav_grid`` and rescale each axis (via ``_TEFF_WEIGHT``,
        ``_LOGG_WEIGHT``, ``_MET_WEIGHT``) onto a common step size, so nearest
        neighbors can be found with a plain Euclidean distance in
        ``build_interpolator``. Sets ``self.grid_teffs``/``grid_loggs``/``grid_mets``
        (unscaled) and ``self.teffs``/``loggs``/``mets`` (scaled).
        """
        # Pull out all gridpoint parameters from sav_grid
        temp_teffs = np.array(sorted(np.unique(sav_grid['teff'])))
        temp_loggs = np.array(sorted(np.unique(sav_grid['logg'])))
        grid_teffs, grid_loggs, grid_mets = [], [], [] 
        for t in temp_teffs:
            for l in temp_loggs:
                preselection = sav_grid[(sav_grid['teff'] == t) & (sav_grid['logg'] == l)]
                available_mets = np.unique(preselection['monh'])
                for m in available_mets:
                    grid_teffs.append(t)
                    grid_loggs.append(l)
                    grid_mets.append(m)
        self.grid_teffs = np.array(grid_teffs)
        self.grid_loggs = np.array(grid_loggs)
        self.grid_mets = np.array(grid_mets)
        self.teff_bounds = (self.grid_teffs.min(), self.grid_teffs.max())
        self.logg_bounds = (self.grid_loggs.min(), self.grid_loggs.max())
        self.met_bounds = (self.grid_mets.min(), self.grid_mets.max())

        # Scale s.t. all step sizes match, except for a slight scaling that will prefer teff over met over logg steps when all would otherwise be equidistant
        temp_mets = np.array(sorted(np.unique(grid_mets)))
        teff_step_original = self._min_step(temp_teffs)
        logg_step_original = self._min_step(temp_loggs)
        met_step_original = self._min_step(temp_mets)

        # RBFInterpolator's default (degree=1) polynomial includes a linear
        # term per axis; if any axis has a single value grid-wide (e.g. a
        # spherical grid with only one metallicity), that axis's column is
        # constant and collinear with the constant term, making the
        # polynomial matrix singular. Dropping to degree=0 (constant term
        # only) sidesteps this - only relevant when an axis is degenerate,
        # so normal (3-axis-varying) grids keep the default degree.
        if len(temp_teffs) < 2 or len(temp_loggs) < 2 or len(temp_mets) < 2:
            self._rbf_degree = 0

        self.teff_step_scale = round(1/teff_step_original * self._TEFF_WEIGHT,4) # weighs Teff as more important / "closer" when selecting neighbors
        self.logg_step_scale = round(1/logg_step_original * self._LOGG_WEIGHT,4) # weighs logg as less important / further away when selecting neighbors
        self.met_step_scale = round(1/met_step_original * self._MET_WEIGHT,4)
        self.teffs = np.around(np.array(grid_teffs) * self.teff_step_scale, 0).astype(np.int32) # the attributes are scaled; these are used for selecting neighbors
        self.loggs = np.around(np.array(grid_loggs)*self.logg_step_scale, 0).astype(np.int32)
        self.mets = np.around(np.array(grid_mets) * self.met_step_scale, 0).astype(np.int32)
        self.gridpoints = len(self.teffs)
    
    def initialize_tags(self, sav_grid, geometry):
        """
        Determine which vector fields (``vtags``) and scalar fields (``stags``)
        of an atmosphere are present in ``sav_grid`` and should be
        interpolated. "height" and "radius" are only included for spherical
        ("SPH") geometry.

        Returns
        -------
        vtags : list of str
            Names of the per-depth-point (vector) fields to interpolate.
        stags : list of str
            Names of the scalar fields to interpolate.
        """
        relevant_atmo = sav_grid[(sav_grid['teff'] == self.grid_teffs[0]) & (sav_grid['logg'] == self.grid_loggs[0]) & (sav_grid['monh'] == self.grid_mets[0])]

        vtags = ["temp", "xne", "xna", "rho"]
        if "tau" in relevant_atmo.dtype.names:
            vtags += ["tau"]
        if "rhox" in relevant_atmo.dtype.names:
            vtags += ["rhox"]
        if "height" in relevant_atmo.dtype.names and geometry == "SPH":
            vtags += ["height"] # Height will store height+radius during interpol

        stags = ["teff", "logg", "monh", "vturb", "lonh"]
        if "radius" in relevant_atmo.dtype.names and geometry == "SPH":
            stags += ["radius"]

        return vtags, stags

    def initialize_depth(self, sav_grid):
        """
        Determine the common number of depth points (``self.depthpoints``) to
        use for every vector field, based on the first grid atmosphere. Falls
        back to the count of finite, positive "tau" values if "tau" has
        trailing NaNs or zero-padding, since those do not represent real depth
        points.
        """
        # Use first atmosphere to define sizes
        relevant_atmo = sav_grid[(sav_grid['teff'] == self.grid_teffs[0]) & (sav_grid['logg'] == self.grid_loggs[0]) & (sav_grid['monh'] == self.grid_mets[0])]
        temp = [len(relevant_atmo[vtag]) for vtag in self.vtags]
        # check for irregularities in tau
        if np.isnan(relevant_atmo["tau"][-1]): # Nan at the end
            temp_notnan = relevant_atmo["tau"][~np.isnan(relevant_atmo["tau"])]
            return len(temp_notnan)
        elif relevant_atmo["tau"][-1] == 0: # zeroes at the end
            temp_notnan = np.array(relevant_atmo["tau"])[np.array(relevant_atmo["tau"])>0.0]
            return len(temp_notnan)
        return np.min(temp) #return unaltered array length
        

    def read_in_grid(self, sav_grid):
        """
        Read every grid atmosphere's ``vtags``/``stags`` fields, convert them
        via ``to_interp_space_vector``/``to_interp_space_scalar``, and stack
        them into ``self.array_for_interpol`` (one flattened row per grid
        point: all vector fields concatenated, then the scalars, then the
        fully-realized H=12 abundance pattern). This is the array the RBF
        interpolant is fit to.

        The abundance pattern is appended as a plain (untransformed) block
        rather than a vtag: it's already logarithmic and isn't
        ``depthpoints`` long. Since ``RBFInterpolator`` applies the same
        local weights to every output column, interpolating the realized
        pattern directly keeps it automatically consistent with the
        separately-interpolated "monh" stag - no extra rescaling needed.
        """
        self.n_abund = sav_grid["abund"].shape[-1]
        abund_offset = len(self.vtags) * self.depthpoints + len(self.stags)
        width = abund_offset + self.n_abund
        self.array_for_interpol = np.zeros(dtype='float64', shape=(self.gridpoints, width))
        for i in range(self.gridpoints): # Loop over all atmospheres
            relevant_atmo = sav_grid[(sav_grid['teff'] == self.grid_teffs[i]) & (sav_grid['logg'] == self.grid_loggs[i]) & (sav_grid['monh'] == self.grid_mets[i])]
            for n, vtag in enumerate(self.vtags): # Read in vectors
                self.array_for_interpol[i,n*self.depthpoints:(n+1)*self.depthpoints] = self.to_interp_space_vector(vtag, relevant_atmo[vtag][:self.depthpoints], relevant_atmo)
            for s, stag in enumerate(self.stags): # Read in scalars
                self.array_for_interpol[i,len(self.vtags)*self.depthpoints + s] = self.to_interp_space_scalar(stag, relevant_atmo[stag])
            self.array_for_interpol[i, abund_offset:] = self._realized_abund_h12(relevant_atmo)

    def _realized_abund_h12(self, relevant_atmo):
        """
        Fully-realized (metallicity-baked-in) abundance pattern for one grid
        atmosphere, converted to the H=12 scale. ``relevant_atmo["abund"]``
        is already an ``Abund`` instance (Atmosphere builds it from the raw
        grid field on access), whose internal ``_pattern`` is exactly this
        realized H=12 pattern (verified against real MARCS grid data: metal
        abundances shift dex-for-dex with [M/H]) - no metallicity needs to
        be added or removed here.
        """
        return np.array(relevant_atmo["abund"]._pattern, dtype=float, copy=True)


class AtmosphereInterpolator:
    def __init__(self, depth=None, interp=None, geom=None, lfs_atmo=None, verbose=0):
        self.depth = depth
        self.interp = interp
        self.geom = geom
        if lfs_atmo is None:
            lfs_atmo = setup_atmo()
        self.lfs_atmo = lfs_atmo

        self.source = None
        self.atmo_grid = None
        self.rbf_grid = None
        self.rbf_grid_key = None
        self.verbose = verbose

    def interp_atmo_grid(self, atmo_grid, teff, logg, monh, interpolation_policy="allow"):
        """
        General routine to interpolate in 3D grid of model atmospheres

        Parameters
        -----
        Teff : float
            effective temperature of desired model (K).
        logg : float
            logarithmic gravity of desired model (log cm/s/s).
        MonH : float
            metalicity of desired model.
        atmo_in : Atmosphere
            Input atmosphere
        verbose : {0, 1}, optional
            how much information to plot/print (default: 0)
        plot : bool, optional
            wether to plot debug information (default: False)
        reload : bool
            wether to reload atmosphere information from disk (default: False)
        interpolation_policy : {"allow", "error"}, optional
            Only used by the RBF path (``interp="RBF"``): what to do when the
            requested Teff/logg/[M/H] falls outside the grid's range.
            "allow" (default) logs and extrapolates; "error" raises
            AtmosphereError. Matches the naming/semantics of the
            (currently unmerged) ``interpolation_policy`` feature.

        Returns
        -------
        atmo : Atmosphere
            interpolated atmosphere data
        """

        # Internal parameters.
        if not isinstance(atmo_grid, AtmosphereGrid):
            if self.atmo_grid is None or self.source != atmo_grid:
                atmo_file = self.lfs_atmo.get(atmo_grid)
                self.source = atmo_grid
                self.atmo_grid = SavFile(
                    atmo_file, source=self.source, lfs=self.lfs_atmo
                )
        else:
            self.atmo_grid = atmo_grid
            self.source = atmo_grid.source

        if self.geom is not None:
            if self.geom == "PP":
                atmo_grid = self.atmo_grid[self.atmo_grid.radius <= 1]
            elif self.geom == "SPH":
                atmo_grid = self.atmo_grid[self.atmo_grid.radius > 1]
        else:
            atmo_grid = self.atmo_grid

        assert len(atmo_grid) > 0, "No atmospheres for this geometry"

        # Get field names in ATMO and ATMO_GRID structures.
        depth = self.determine_depth_scale(self.depth, atmo_grid)
        interp = self.determine_interpolation_scale(self.interp, atmo_grid)



        #### Option 1: if atmo_interp is rbf, call alternative method
        if interp == "RBF":
            rbf_key = (self.source, self.geom)
            if not isinstance(self.rbf_grid, RbfGrid) or self.rbf_grid_key != rbf_key:
                self.rbf_grid = RbfGrid(atmo_grid, self.geom)
                self.rbf_grid_key = rbf_key
            rbf_grid = self.rbf_grid
            atmo = self.interpolate_RBF(teff, logg, monh, rbf_grid, interpolation_policy=interpolation_policy)
            if self.geom == "SPH":
                geom = "SPH"
                radius = atmo.radius
            else:
                geom = "PP"
                radius = None

        
        #### Option 2 (default): If atmo_interp is Rhox or Tau or sth similar, do the usual pysme interpolation by pairs
        else:
            # Find the corner models bracketing the given values
            icor = self.find_corner_models(teff, logg, monh, atmo_grid)

            
            # Interpolate the corner models
            atmo = self.interpolate_corner_models(
                teff, logg, monh, icor, atmo_grid, interp=interp
            )

            # TODO: Or should we only consider spherical models for interpolation of spherical modesl are requested?
            geom, radius = self.spherical_model_correction(atmo_grid, icor, logg)

        # Create ATMO.GEOM, if necessary, and set value.
        if self.geom is not None and self.geom != geom:
            if self.geom == "SPH":
                raise AtmosphereError(
                    "Input ATMO.GEOM='%s' was requested but the model only supports PP (at this point)."
                    % self.geom
                )
            else:
                logger.info(
                    "Input ATMO.GEOM='%s' overrides '%s' from grid.",
                    self.geom,
                    geom,
                )

        # Add standard ATMO input fields, if they are missing from ATMO_IN.
        atmo.depth = depth
        atmo.interp = interp
        atmo.geom = geom
        if radius is not None:
            atmo.radius = radius
        atmo.source = self.source
        atmo.method = "grid"

        return atmo

    def interpolate_RBF(self, teff, logg, monh, rbf_grid, interpolation_policy="allow"):
        """
        Interpolate a model atmosphere at (teff, logg, monh) using radial
        basis function interpolation over ``rbf_grid``, as an alternative to
        the corner-model interpolation used by ``interpolate_corner_models``.

        Parameters
        ----------
        teff : float
            effective temperature of desired model (K).
        logg : float
            logarithmic gravity of desired model (log cm/s/s).
        monh : float
            metallicity of desired model.
        rbf_grid : RbfGrid
            Precomputed RBF grid (rescaled coordinates, cached interpolator)
            to interpolate from.
        interpolation_policy : {"allow", "error"}, optional
            What to do when (teff, logg, monh) falls outside the grid's
            range (default: "allow", i.e. log and extrapolate).

        Returns
        -------
        atmo : Atmosphere
            interpolated atmosphere data

        Raises
        ------
        AtmosphereError
            If the query is out of range and ``interpolation_policy`` isn't "allow".
        """
        rbf_grid.check_within_bounds(teff, logg, monh, interpolation_policy)
        target_params = np.array([round(teff * rbf_grid.teff_step_scale,3), round(logg * rbf_grid.logg_step_scale,3), round(monh * rbf_grid.met_step_scale,3)]) #scaled so all three params have same step-size
        atmo = Atmo(interp="rhox") #purely to initliaize the structure

        # INTERPOLATION (uses the interpolator cached on rbf_grid, built once for the whole grid)
        target_vector = rbf_grid.interpolator([target_params])[0]

        # READ-OUT
        def from_interp_space_vector(vtag, values, atmo):
            values = np.asarray(values)
            if vtag == "height":
                return 10.0 ** (values) - atmo["radius"] # separate height and radius
            return 10.0 ** values
        def from_interp_space_scalar(stag, value):
            if stag == "radius":
                return 10.0 ** value
            return value

        for s, stag in enumerate(rbf_grid.stags): # Read out scalars
            atmo[stag] = from_interp_space_scalar(stag, target_vector[len(rbf_grid.vtags)*rbf_grid.depthpoints + s])
        for n, vtag in enumerate(rbf_grid.vtags): # Read out vectors
            atmo[vtag] = from_interp_space_vector(vtag, target_vector[n * rbf_grid.depthpoints : (n+1) * rbf_grid.depthpoints], atmo)

        # opflag and wlstd are administrative fields, not smoothly varying
        # physical structure, so copy them from the nearest grid point
        # instead - the same way the standard corner-model path copies them
        # verbatim rather than interpolating them.
        nearest = np.argmin(
            (rbf_grid.teffs - target_params[0]) ** 2
            + (rbf_grid.loggs - target_params[1]) ** 2
            + (rbf_grid.mets - target_params[2]) ** 2
        )
        atmo.opflag = rbf_grid.sav_grid["opflag"][nearest]
        atmo.wlstd = rbf_grid.sav_grid["wlstd"][nearest]

        # The abundance pattern, unlike opflag/wlstd, does vary smoothly (it's
        # metallicity- and enhancement-dependent), so it's part of the RBF fit
        # itself (see read_in_grid/_realized_abund_h12) rather than copied
        # from the nearest point - interpolating the realized pattern with
        # the same local weights used everywhere else keeps it automatically
        # consistent with the interpolated "monh" above.
        abund_offset = len(rbf_grid.vtags) * rbf_grid.depthpoints + len(rbf_grid.stags)
        interpolated_pattern = target_vector[abund_offset:]
        # atmo.monh was already set above (it's one of rbf_grid.stags) and,
        # since Atmosphere.monh is a passthrough to atmo.abund.monh, it would
        # otherwise be silently reset to 0 by replacing atmo.abund below -
        # capture it first and carry it over to the new Abund. This does not
        # reintroduce double-counting: the pattern is already fully realized,
        # and get_pattern_abundance() (unlike Abund.__call__/__repr__) never
        # applies .monh on top of it - matching how a real grid atmosphere's
        # Abund object already stores both its own realized pattern and its
        # own (redundant, non-additive) monh value.
        interpolated_monh = atmo.monh
        atmo.abund = Abund(monh=interpolated_monh, pattern=interpolated_pattern, type="H=12")

        return atmo

    def interp_atmo_pair(self, atmo1, atmo2, frac, interpvar="RHOX", itop=0):
        """
        Interpolate between two model atmospheres, accounting for shifts in
        the mass column density or optical depth scale.

        How it works:

        1. The second atmosphere is fitted onto the first, individually for
        each of the four atmospheric quantitites: T, xne, xna, rho.
        The fitting uses a linear shift in both the (log) depth parameter and
        in the (log) atmospheric quantity. For T, the midpoint of the two
        atmospheres is aligned for the initial guess. The result of this fit
        is used as initial guess for the other quantities. A penalty function
        is applied to each fit, to avoid excessively large shifts on the
        depth scale.
        2. The mean of the horizontal shift in each parameter is used to
        construct the common output depth scale.
        3. Each atmospheric quantity is interpolated after shifting the two
        corner models by the amount determined in step 1), rescaled by the
        interpolation fraction (frac).

        Parameters
        ------
        atmo1 : Atmosphere
            first atmosphere to interpolate
        atmo2 : Atmosphere
            second atmosphere to interpolate
        frac : float
            interpolation fraction: 0.0 -> atmo1 and 1.0 -> atmo2
        interpvar : {"RHOX", "TAU"}, optional
            atmosphere interpolation variable (default:"RHOX").
        itop : int, optional
            index of top point in the atmosphere to use. default
            is to use all points (itop=0). use itop=1 to clip top depth point.
        atmop : array[5, ndep], optional
            interpolated atmosphere prediction (for plots)
            Not needed if atmospheres are provided as structures. (default: None)
        verbose : {0, 5}, optional
            diagnostic print level (default 0: no printing)
        plot : {-1, 0, 1}, optional
            diagnostic plot level. Larger absolute
            value yields more plots. Negative values cause a wait for keypress
            after each plot. (default: 0, no plots)
        old : bool, optional
            also plot result from the old interpkrz2 algorithm. (default: False)

        Returns
        ------
        atmo : Atmosphere
            interpolated atmosphere
            .RHOX (vector[ndep]) mass column density (g/cm^2)
            .TAU  (vector[ndep]) reference optical depth (at 5000 Å)
            .temp (vector[ndep]) temperature (K)
            .xne  (vector[ndep]) electron number density (1/cm^3)
            .xna  (vector[ndep]) atomic number density (1/cm^3)
            .rho  (vector[ndep]) mass density (g/cm^3)
        """
        # """
        # History
        # -------
        # 2004-Apr-15 Valenti
        #     Initial coding.
        # MB
        #     interpolation on TAU scale
        # 2012-Oct-30 TN
        #     Rewritten to use either column mass (RHOX) or
        #     reference optical depth (TAU) as vertical scale. Shift-interpolation
        #     algorithms have been improved for stability in cool dwarfs (<=3500 K).
        #     The reference optical depth scale is preferred in terms of interpolation
        #     accuracy across most of parameter space, with significant improvement for
        #     both cool models (where depth vs temperature is rather flat) and hot
        #     models (where depth vs temperature exhibits steep transitions).
        #     Column mass depth is used by default for backward compatibility.
        # 2013-May-17 Valenti
        #     Use frac to weight the two shifted depth scales,
        #     rather than simply averaging them. This change fixes discontinuities
        #     when crossing grid nodes.
        # 2013-Sep-10 Valenti
        #     Now returns an atmosphere structure instead of a
        #     [5,NDEP] atmosphere array. This was necessary to support interpolation
        #     using one variable (e.g., TAU) and radiative transfer using a different
        #     variable (e.g. RHOX). The atmosphere array could only store one depth
        #     variable, meaning the radiative transfer variable had to be the same
        #     as the interpolation variable. Returns atmo.RHOX if available and also
        #     atmo.TAU if available. Since both depth variables are returned, if
        #     available, this routine no longer needs to know which depth variable
        #     will be used for radiative transfer. Only the interpolation variable
        #     is important. Thus, the interpvar= keyword argument replaces the
        #     type= keyword argument. Very similar code blocks for each atmospheric
        #     quantity have been unified into a single code block inside a loop over
        #     atmospheric quantities.
        # 2013-Sep-21 Valenti
        #     Fixed an indexing bug that affected the output depth
        #     scale but not other atmosphere vectors. Itop clipping was not being
        #     applied to the depth scale ('RHOX' or 'TAU'). Bug fixed by adding
        #     interpvar to vtags. Now atmospheres interpolated with interp_atmo_grid
        #     match output from revision 398. Revisions back to 399 were development
        #     only, so no users should be affected.
        # 2014-Mar-05 Piskunov
        #     Replicated the removal of the bad top layers in models
        #     for interpvar eq 'TAU'
        # """

        # Internal program parameters.
        min_drhox = min_dtau = 0.01  # minimum fractional step in RHOX
        # min_dtau = 0.01  # minimum fractional step in TAU
        interpvar = interpvar.lower()
        ##
        ## Select interpolation variable (RHOX vs. TAU)
        ##

        def _field(record, name, default=None):
            try:
                return getattr(record, name)
            except AttributeError:
                try:
                    return record[name]
                except (KeyError, IndexError, ValueError):
                    return default

        def _as_abund(record):
            abund = _field(record, "abund")
            if isinstance(abund, Abund):
                return abund
            abund_format = getattr(
                self.atmo_grid, "abund_format", _field(record, "abund_format", "sme")
            )
            return Abund(
                monh=_field(record, "monh", 0),
                pattern=np.array(abund, dtype=float, copy=True),
                type=abund_format,
            )

        # Check which depth scales are available in both input atmospheres.
        tags1 = atmo1.dtype.names
        tags2 = atmo2.dtype.names
        ok_tau = "tau" in tags1 and "tau" in tags2
        ok_rhox = "rhox" in tags1 and "rhox" in tags2
        if not ok_tau and not ok_rhox:
            raise AtmosphereError(
                "atmo1 and atmo2 structures must both contain RHOX or TAU"
            )

        # Set interpolation variable, if not specified by keyword argument.
        if interpvar is None:
            interpvar = "tau" if ok_tau else "rhox"
        if interpvar != "tau" and interpvar != "rhox":
            raise AtmosphereError("interpvar must be 'TAU' (default) or 'RHOX'")

        has_spherical_height = (
            _field(atmo1, "radius", 0) > 1
            and _field(atmo2, "radius", 0) > 1
            and "height" in tags1
            and "height" in tags2
        )

        def to_interp_space(vtag, values, atm):
            values = np.asarray(values)
            if vtag == "height":
                return np.log10(values + atm["radius"]) # combine height and radius for interpolation
            return np.log10(values)

        def from_interp_space(vtag, values):
            values = np.asarray(values)
            return 10.0 ** values

        ##
        ## Define depth scale for both atmospheres
        ##

        # Define depth scale for atmosphere #1
        mask1 = np.full(len(atmo1[interpvar]), True)
        mask1[:itop] = False

        itop1 = itop
        while (atmo1[interpvar][itop1 + 1] / atmo1[interpvar][itop1] - 1) <= min_drhox:
            mask1[itop1] = False
            itop1 += 1

        mask1 &= atmo1[interpvar] != 0
        mask1 &= np.isfinite(atmo1[interpvar])
        ndep1 = np.count_nonzero(mask1)

        depth1 = np.log10(atmo1[interpvar][mask1])

        # Define depth scale for atmosphere #2
        mask2 = np.full(len(atmo2[interpvar]), True)
        mask2[:itop] = False
        itop2 = itop
        while (atmo2[interpvar][itop2 + 1] / atmo2[interpvar][itop2] - 1) <= min_drhox:
            mask2[itop2] = False
            itop2 += 1

        mask2 &= atmo2[interpvar] != 0
        mask2 &= np.isfinite(atmo2[interpvar])
        ndep2 = np.count_nonzero(mask2)

        depth2 = np.log10(atmo2[interpvar][mask2])

        ##
        ## Prepare to find best shift parameters for each atmosphere vector.
        ##

        # List of atmosphere vectors that need to be shifted.
        # The code below assumes 'TEMP' is the first vtag in the list.
        vtags = ["temp", "xne", "xna", "rho", interpvar]
        if interpvar == "rhox" and ok_tau:
            vtags += ["tau"]
        if interpvar == "tau" and ok_rhox:
            vtags += ["rhox"]
        if has_spherical_height:
            vtags += ["height"]
        nvtag = len(vtags)

        # Adopt arbitrary uncertainties for shift determinations.
        err1 = np.full(ndep1, 0.05)

        # Initial guess for TEMP shift parameters.
        # Put depth and TEMP midpoints for atmo1 and atmo2 on top of one another.
        npar = 4
        ipar = np.zeros(npar, dtype="f4")
        temp1 = np.log10(_field(atmo1, "temp")[mask1])
        temp2 = np.log10(_field(atmo2, "temp")[mask2])
        mid1 = np.argmin(np.abs(temp1 - 0.5 * (temp1[1] + temp1[-2])))
        mid2 = np.argmin(np.abs(temp2 - 0.5 * (temp2[1] + temp2[-2])))
        ipar[0] = depth1[mid1] - depth2[mid2]  # horizontal
        ipar[1] = temp1[mid1] - temp2[mid2]  # vertical

        # Apply a constraint on the fit, to avoid runaway for cool models, where
        # the temperature structure is nearly linear with both TAU and RHOX.
        constraints = np.zeros(npar)
        constraints[0] = 0.5  # weakly constrain the horizontal shift

        # For first pass ('TEMP'), use all available depth points.
        igd = np.isfinite(depth1)
        ngd = igd.size

        ##
        ## Find best shift parameters for each atmosphere vector.
        ##

        # Loop through atmosphere vectors.
        pars = np.zeros((nvtag, npar))
        for ivtag, vtag in enumerate(vtags):

            # Find vector in each structure.
            if vtag not in tags1:
                raise AtmosphereError("atmo1 does not contain " + vtag)
            if vtag not in tags2:
                raise AtmosphereError("atmo2 does not contain " + vtag)

            vect1 = to_interp_space(vtag, atmo1[vtag][mask1], atmo1)
            vect2 = to_interp_space(vtag, atmo2[vtag][mask2], atmo2)

            # Fit the second atmosphere onto the first by finding the best horizontal
            # shift in depth2 and the best vertical shift in vect2.
            pars[ivtag], _ = self.interp_atmo_constrained(
                depth1[igd],
                vect1[igd],
                err1[igd],
                ipar,
                x2=depth2,
                y2=vect2,
                y1=vect1,
                ndep=ngd,
                constraints=constraints,
            )

            # After first pass ('TEMP'), adjust initial guess and restrict depth points.
            if ivtag == 0:
                ipar = [pars[0, 0], 0.0, 0.0, 0.0]
                try:
                    igd = np.where(
                        (depth1 >= min(depth2[igd]) + ipar[0])
                        & (depth1 <= max(depth2[igd]) + ipar[0])
                    )[0]
                except IndexError:
                    pass
                if igd.size < 2:
                    raise AtmosphereError("unstable shift in temperature")

        ##
        ## Use mean shift to construct output depth scale.
        ##

        # Calculate the mean depth2 shift for all atmosphere vectors.
        xsh = np.sum(pars[:, 0]) / nvtag

        # Base the output depth scale on the input scale with the fewest depth points.
        # Combine the two input scales, if they have the same number of depth points.
        depth1f = depth1 - xsh * frac
        depth2f = depth2 + xsh * (1 - frac)
        if ndep1 > ndep2:
            depth = depth2f
        elif ndep1 == ndep2:
            depth = depth1f * (1 - frac) + depth2f * frac
        elif ndep1 < ndep2:
            depth = depth1f
        ndep = len(depth)

        ##
        ## Interpolate input atmosphere vectors onto output depth scale.
        ##

        # Loop through atmosphere vectors.
        vects = np.zeros((nvtag, ndep))
        for ivtag, (vtag, par) in enumerate(zip(vtags, pars)):

            # Extract data
            vect1 = to_interp_space(vtag, atmo1[vtag][mask1], atmo1)
            vect2 = to_interp_space(vtag, atmo2[vtag][mask2], atmo2)

            # Identify output depth points that require extrapolation of atmosphere vector.
            depth1f = depth1 - par[0] * frac
            depth2f = depth2 + par[0] * (1 - frac)
            x1max = np.max(depth1f)
            x2max = np.max(depth2f)
            iup = (depth > x1max) | (depth > x2max)
            nup = np.count_nonzero(iup)
            checkup = (nup >= 1) and abs(frac - 0.5) <= 0.5 and ndep1 == ndep2

            # Combine shifted vect1 and vect2 structures to get output vect.
            vect1f = self.interp_atmo_func(depth, -frac * par, x2=depth1, y2=vect1)
            vect2f = self.interp_atmo_func(depth, (1 - frac) * par, x2=depth2, y2=vect2)
            vect = (1 - frac) * vect1f + frac * vect2f
            ends = [vect1[ndep1 - 1], vect[ndep - 1], vect2[ndep2 - 1]]
            if (
                checkup
                and np.median(ends) != vect[ndep - 1]
                and (
                    abs(vect1[ndep1 - 1] - 4.2) < 0.1
                    or abs(vect2[ndep2 - 1] - 4.2) < 0.1
                )
            ):
                vect[iup] = vect2f[iup] if x1max < x2max else vect1f[iup]
            vects[ivtag] = vect

        ##
        ## Construct output structure
        ##

        # Construct output structure with interpolated atmosphere.
        # Might be wise to interpolate abundances, in case those ever change.
        atmo = Atmo(interp=interpvar)
        stags = ["teff", "logg", "monh", "vturb", "lonh", "abund"]
        ndep_orig = len(_field(atmo1, "temp"))
        if has_spherical_height:
            pair_logg = (1 - frac) * atmo1["logg"] + frac * atmo2["logg"]
            mass1 = atmo1["logg"] - logg_sun - 2 * np.log10(R_sun / atmo1["radius"])
            mass2 = atmo2["logg"] - logg_sun - 2 * np.log10(R_sun / atmo2["radius"])
            mass = np.mean([mass1,mass2])
            pair_radius = R_sun * 10 ** ((logg_sun - pair_logg + mass) * 0.5)
            
        else:
            pair_radius = 0
            pair_logg = 0

        for tag in tags1:

            # Default is to copy value from atmo1. Trim vectors.
            value = atmo1[tag]
            if np.size(value) == ndep_orig and tag != "abund":
                value = value[:ndep]

            # Scalar quantities that should be interpolated using frac.
            if tag in stags:
                if tag == "abund":
                    value = (1 - frac) * _as_abund(atmo1) + frac * _as_abund(atmo2)
                elif tag in tags2:
                    value = (1 - frac) * atmo1[tag] + frac * atmo2[tag]
                else:
                    value = atmo1[tag]
                
            # Vector quantities that have already been interpolated.
            if tag in vtags:
                ivtag = [i for i in range(nvtag) if tag == vtags[i]][0]

                value = from_interp_space(tag, vects[ivtag])
                if tag == "height":
                    if not has_spherical_height:
                        raise AtmosphereError(
                            "Cannot interpolate 'height': both atmospheres must be "
                            "spherical (radius > 1) to recover height from the "
                            "combined height+radius quantity used for interpolation."
                        )
                    value = value - pair_radius

            # Store the pair-consistent radius, so a subsequent interpolation stage
            # recombines 'height' with the same radius that was just subtracted from it.
            if tag == "radius" and has_spherical_height:
                value = pair_radius

            # Remaining cases.
            if tag == "ndep":
                value = ndep

            # Create or add to output structure.
            atmo[tag] = value
        return atmo

    def determine_depth_scale(self, depth, atmo_grid):
        """
        Determine ATMO.DEPTH radiative transfer depth variable. Order of precedence:

        1. Value of ATMO_IN.DEPTH, if it exists and is set
        2. Value of ATMO_GRID[0].DEPTH, if it exists and is set
        3. 'RHOX', if ATMO_GRID.RHOX exists (preferred over 'TAU' for depth)
        4. 'TAU', if ATMO_GRID.TAU exists

        Check that INTERP is valid and the corresponding field exists in ATMO.

        Parameters
        ----------
        depth : {"RHOX", "TAU", None}
            requested value, or None for autoselection based in available grid
        atmo_grid : AtmosphereGrid
            input atmosphere grid to interpolate on

        Returns
        -------
        depth : {"TAU", "RHOX"}
            The chosen depth scale

        Raises
        ------
        AtmosphereError
            If an invalid value was set in atmo_in/atmo_grid
        """

        gtags = atmo_grid.dtype.names

        if depth is not None:
            depth = depth
        elif "depth" in gtags and atmo_grid.depth is not None:
            depth = atmo_grid.depth
        elif "rhox" in gtags:
            depth = "RHOX"
        elif "tau" in gtags:
            depth = "TAU"
        else:
            raise AtmosphereError("no value for ATMO.DEPTH")
        if depth != "TAU" and depth != "RHOX":
            raise AtmosphereError(
                "ATMO.DEPTH must be 'TAU' or 'RHOX', not '%s'" % depth
            )
        if depth.lower() not in gtags:
            raise AtmosphereError(
                "ATMO.DEPTH='{}', but ATMO. {} does not exist".format(depth, depth)
            )
        return depth

    def determine_interpolation_scale(self, interp, atmo_grid):
        """
        Determine ATMO.INTERP interpolation variable. Order of precedence:

        1. Value of ATMO_IN.INTERP, if it exists and is set
        2. Value of ATMO_GRID[0].INTERP, if it exists and is set
        3. 'TAU', if ATMO_GRID.TAU exists (preferred over 'RHOX' for interpolation)
        4. 'RHOX', if ATMO_GRID.RHOX exists

        Check that INTERP is valid and the corresponding field exists in ATMO.

        Parameters
        ----------
        interp : {"RHOX", "TAU", "RBF", None}
            requested interpolation axis, or None for autoselect
        atmo_grid : AtmosphereGrid
            Atmosphere grid for interpolation

        Returns
        -------
        interp: {"RHOX", "TAU", "RBF"}
            the interpolation axis

        Raises
        ------
        AtmosphereError
            if a non valid parameter is set in atmo_in/atmo_grid
        """

        gtags = atmo_grid.dtype.names

        if interp is not None:
            interp = interp
        elif atmo_grid.interp is not None:
            interp = atmo_grid.interp
        elif "tau" in gtags:
            interp = "TAU"
        elif "rhox" in gtags:
            interp = "RHOX"
        else:
            raise AtmosphereError("no value for ATMO.INTERP")
        if interp not in ["TAU", "RHOX", "RBF"]:
            raise AtmosphereError(
                "ATMO.INTERP must be 'TAU', 'RHOX', or 'RBF', not '%s'" % interp
            )
        if interp != "RBF" and interp.lower() not in gtags:
            raise AtmosphereError(
                "ATMO.INTERP='{}', but ATMO. {} does not exist".format(interp, interp)
            )
        return interp

    def find_corner_models(self, teff, logg, monh, atmo_grid):
        """
        Find the models in the grid that bracket the given stellar parameters

        The purpose of the first major set of code blocks is to find values
        of [M/H] in the grid that bracket the requested [M/H]. Then in this
        subset of models, find values of log(g) in the subgrid that bracket
        the requested log(g). Then in this subset of models, find values of
        Teff in the subgrid that bracket the requested Teff. The net result
        is a set of 8 models in the grid that bracket the requested stellar
        parameters. Only these 8 "corner" models will be used in the
        interpolation that constitutes the remainder of the program.

        Parameters
        ----------
        teff : float
            effective temperature
        logg : float
            surface gravity
        monh : float
            metallicity
        atmo_grid : AtmosphereGrid
            atmosphere grid to search models in
        verbose : int, optional
            determines how many debugging messages to print, by default 0
        """

        nb = 2  # number of bracket points

        # *** DETERMINATION OF METALICITY BRACKET ***
        # Find unique set of [M/H] values in grid.
        mlist = np.unique(atmo_grid.monh)  # list of unique [M/H]

        # Test whether requested metalicity is in grid.
        mmin = np.min(mlist)  # range of [M/H] in grid
        mmax = np.max(mlist)
        if monh > mmax:  # true: [M/H] too large
            logger.info(
                "interp_atmo_grid: requested [M/H] (%.3f) larger than max grid value (%.3f). extrapolating.",
                monh,
                mmax,
            )
        if monh < mmin:  # true: logg too small
            raise AtmosphereError(
                "interp_atmo_grid: requested [M/H] (%.3f) smaller than min grid value (%.3f). returning."
                % (monh, mmin)
            )

        # Find closest two [M/H] values in grid that bracket requested [M/H].
        if monh <= mmax:
            mlo = np.max(mlist[mlist <= monh])
            mup = np.min(mlist[mlist >= monh])
        else:
            mup = mmax
            mlo = np.max(mlist[mlist < mup])
        mb = [mlo, mup]  # bounding [M/H] values

        # Trace diagnostics.
        if self.verbose >= 5:
            logger.info("[M/H]: %.3f < %.3f < %.3f", mlo, monh, mup)

        # *** DETERMINATION OF LOG(G) BRACKETS AT [M/H] BRACKET VALUES ***
        # Set up for loop through [M/H] bounds.
        gb = np.zeros((nb, nb))  # bounding gravities
        for i in range(nb):
            # Find unique set of gravities at boundary below [M/H] value.
            im = atmo_grid.monh == mb[i]  # models on [M/H] boundary
            glist = np.unique(atmo_grid[im].logg)  # list of unique gravities

            # Test whether requested logarithmic gravity is in grid.
            gmin = np.min(glist)  # range of gravities in grid
            gmax = np.max(glist)
            if logg > gmax:  # true: logg too large
                logger.info(
                    "interp_atmo_grid: requested log(g) (%.3f) larger than max grid value (%.3f). extrapolating.",
                    logg,
                    gmax,
                )

            if logg < gmin:  # true: logg too small
                raise AtmosphereError(
                    "interp_atmo_grid: requested log(g) (%.3f) smaller than min grid value (%.3f). returning."
                    % (logg, gmin)
                )

            # Find closest two gravities in Mlo subgrid that bracket requested gravity.
            if logg <= gmax:
                glo = np.max(glist[glist <= logg])
                gup = np.min(glist[glist >= logg])
            else:
                gup = gmax
                glo = np.max(glist[glist < gup])
            gb[i] = [glo, gup]  # store boundary values.

            # Trace diagnostics.
            if self.verbose >= 5:
                logger.info(
                    "log(g) at [M/H]=%.3f: %.3f < %.3f < %.3f",
                    mb[i],
                    glo,
                    logg,
                    gup,
                )

        # End of loop through [M/H] bracket values.
        # *** DETERMINATION OF TEFF BRACKETS AT [M/H] and LOG(G) BRACKET VALUES ***
        # Set up for loop through [M/H] and log(g) bounds.
        tb = np.zeros((nb, nb, nb))  # bounding temperatures
        for ig in range(nb):
            for im in range(nb):
                # Find unique set of gravities at boundary below [M/H] value.
                it = (atmo_grid.monh == mb[im]) & (
                    atmo_grid.logg == gb[im, ig]
                )  # models on joint boundary
                tlist = np.unique(atmo_grid[it].teff)  # list of unique temperatures

                # Test whether requested temperature is in grid.
                tmin = np.min(tlist)  # range of temperatures in grid
                tmax = np.max(tlist)
                if teff > tmax:  # true: Teff too large
                    raise AtmosphereError(
                        "interp_atmo_grid: requested Teff (%i) larger than max grid value (%i). returning."
                        % (teff, tmax)
                    )
                if teff < tmin:  # true: logg too small
                    logger.info(
                        "interp_atmo_grid: requested Teff (%i) smaller than min grid value (%i). extrapolating.",
                        teff,
                        tmin,
                    )

                # Find closest two temperatures in subgrid that bracket requested Teff.
                if teff > tmin:
                    tlo = np.max(tlist[tlist <= teff])
                    tup = np.min(tlist[tlist >= teff])
                else:
                    tlo = tmin
                    tup = np.min(tlist[tlist > tlo])
                tb[im, ig, :] = [tlo, tup]  # store boundary values.

                # Trace diagnostics.
                if self.verbose >= 5:
                    logger.info(
                        "Teff at log(g)=%.3f and [M/H]=%.3f: %i < %i < %i",
                        gb[im, ig],
                        mb[im],
                        tlo,
                        teff,
                        tup,
                    )

        # End of loop through log(g) and [M/H] bracket values.

        # Find and save atmo_grid indices for the 8 corner models.
        icor = np.zeros((nb, nb, nb), dtype=int)
        for it, ig, im in itertools.product(range(nb), repeat=3):
            iwhr = np.where(
                (atmo_grid.teff == tb[im, ig, it])
                & (atmo_grid.logg == gb[im, ig])
                & (atmo_grid.monh == mb[im])
            )[0]
            nwhr = iwhr.size
            if nwhr != 1:
                logger.info(
                    "interp_atmo_grid: %i models in grid with [M/H]=%.1f, log(g)=%.1f, and Teff=%i",
                    nwhr,
                    mb[im],
                    gb[im, ig],
                    tb[im, ig, it],
                )
            icor[im, ig, it] = iwhr[0]

        # Trace diagnostics.
        def _field(model, name):
            try:
                return getattr(model, name)
            except AttributeError:
                return model[name]

        if self.verbose >= 1:
            logger.info("Teff=%i,  log(g)=%.3f,  [M/H]=%.3f:", teff, logg, monh)
            logger.info("indx  M/H  g   Teff     indx  M/H  g   Teff")
            for im in range(nb):
                for ig in range(nb):
                    i0 = icor[im, ig, 0]
                    i1 = icor[im, ig, 1]
                    logger.info(
                        i0,
                        _field(atmo_grid[i0], "monh"),
                        _field(atmo_grid[i0], "logg"),
                        _field(atmo_grid[i0], "teff"),
                        i1,
                        _field(atmo_grid[i1], "monh"),
                        _field(atmo_grid[i1], "logg"),
                        _field(atmo_grid[i1], "teff"),
                    )
        return icor

    def interpolate_corner_models(
        self, teff, logg, monh, icor, atmo_grid, interp="RHOX", itop=1
    ):
        """
        Interpolate over the 8 corner models in the cube around the stellar parameters

        The code below interpolates between 8 corner models to produce
        the output atmosphere. In the first step, pairs of models at each
        of the 4 combinations of log(g) and Teff are interpolated to the
        desired value of [M/H]. These 4 new models are then interpolated
        to the desired value of log(g), yielding 2 models at the requested
        [M/H] and log(g). Finally, this pair of models is interpolated
        to the desired Teff, producing a single output model.

        Interpolation is done on the logarithm of all quantities to improve
        linearity of the fitted data. Kurucz models sometimes have very small
        fractional steps in mass column at the top of the atmosphere. These
        cause wild oscillations in splines fitted to facilitate interpolation
        onto a common depth scale. To circumvent this problem, we ignore the
        top point in the atmosphere by setting itop=1.

        Parameters
        ----------
        teff : float
            effective temperature
        logg : float
            surface gravity
        monh : float
            metallicity
        icor : array_like
            indices of the corner models
        atmo_grid : AtmosphereGrid
            atmosphere grid with input models
        interp : str, optional
            interpolation axis, by default "RHOX"
        itop : int, optional
            index of the topmost point in the RHOX scale, by default 1

        Returns
        -------
        atmo3: Atmosphere
            Interpolated atmosphere
        """

        # We do this for every pair of atmosphere models
        def interpolate(m0, m1, p, param, **kwargs):
            try:
                p0 = getattr(m0, param)
            except AttributeError:
                p0 = m0[param]
            try:
                p1 = getattr(m1, param)
            except AttributeError:
                p1 = m1[param]
            pfrac = (p - p0) / (p1 - p0) if p0 != p1 else 0
            return self.interp_atmo_pair(m0, m1, pfrac, interpvar=interp, **kwargs)

        # Interpolate 8 corner models to create 4 models at the desired [M/H].
        atmo = [[None, None], [None, None]]
        for (i, j) in np.ndindex(2, 2):
            m0 = atmo_grid[icor[0, j, i]]
            m1 = atmo_grid[icor[1, j, i]]
            atmo[i][j] = interpolate(m0, m1, monh, "monh", itop=itop)

        # Interpolate 4 models at the desired [M/H] to create 2 models at desired
        # [M/H] and log(g).
        atmo2 = [None, None]
        for k in range(2):
            atmo2[k] = interpolate(atmo[k][0], atmo[k][1], logg, "logg")

        # Interpolate the 2 models at desired [M/H] and log(g) to create final
        # model at desired [M/H], log(g), and Teff
        atmo3 = interpolate(atmo2[0], atmo2[1], teff, "teff")

        return atmo3

    def spherical_model_correction(self, atmo_grid, icor, logg):
        """
        Correct the radius of the model grids, for the stellar parameters if necessary

        If all interpolated models were spherical, the interpolated model should
        also be reported as spherical. This enables spherical-symmetric radiative
        transfer in the spectral synthesis.

        Formulae for mass and radius at the corners of interpolation cube:
            log(M/Msol) = log g - log g_sol - 2*log(R_sol / R)
            2 * log(R / R_sol) = log g_sol - log g + log(M / M_sol)

        Parameters
        ----------
        atmo_grid : AtmosphereGrid
            complete atmosphere grid
        icor : array_like
            indices for the corner models that were interpolated from the
            atmosphere grid

        Returns
        -------
        geom : {"PP", "SPH"}
            The geometry of the interpolated models
        radius : None, array
            The corrected radius of the models, if geom == "SPH". Otherwise None.
        """

        gtags = atmo_grid.dtype.names
        selection = atmo_grid[icor]
        if "radius" in gtags and "height" in gtags and np.min(selection["radius"]) > 1:
            mass_cor = (
                selection["logg"]
                - logg_sun
                - 2 * np.log10(R_sun / selection["radius"])
            )
            mass = np.mean(mass_cor)
            radius = R_sun * 10 ** ((logg_sun - logg + mass) * 0.5)
            geom = "SPH"
        else:
            geom = "PP"
            radius = None

        return geom, radius

    def interp_atmo_constrained(
        self, x, y, err, par, x2, y2, constraints=None, **kwargs
    ):
        """
        Apply a constraint on each parameter, to have it approach zero

        Parameters
        -------
        x : array[n]
            x data
        y : array[n]
            y data
        err : array[n]
            errors on y data
        par : list[4]
            initial guess for fit parameters - see interp_atmo_func
        x2 : array
            independent variable for tabulated input function
        y2 : array
            dependent variable for tabulated input function
        ** kwargs : dict
            passes keyword arguments to interp_atmo_func

            ndep : int
                number of depth points in supplied quantities
            constraints : array[nconstraint], optional
                error vector for constrained parameters.
                Use errors of 0 for unconstrained parameters.

        Returns
        --------
        ret : list of floats
            best fit parameters
        yfit : array of size (n,)
            best fit to data
        """

        # Evaluate with fixed paramters 3, 4
        # ret = mpfitfun("interp_atmo_func", x, y, err, par, extra=extra, yfit=yfit)
        # TODO: what does the constrain do?
        func = lambda x, *p: self.interp_atmo_func(x, [*p, *par[2:]], x2, y2, **kwargs)
        popt, _ = curve_fit(func, x, y, sigma=err, p0=par[:2])

        ret = [*popt, *par[2:]]
        yfit = func(x, *popt)
        return ret, yfit

    def interp_atmo_func(self, x1, par, x2, y2, ndep=None, y1=None):
        """
        Apply a horizontal shift to x2.
        Apply a vertical shift to y2.
        Interpolate onto x1 the shifted y2 as a function of shifted x2.

        Parameters
        ---------
        x1 : array[ndep1]
            independent variable for output function
        par : array[3]
            shift parameters
            par[0] - horizontal shift for x2
            par[1] - vertical shift for y2
            par[2] - vertical scale factor for y2
        x2 : array[ndep2]
            independent variable for tabulated input function
        y2 : array[ndep2]
            dependent variable for tabulated input function
        ndep : int, optional
            number of depth points in atmospheric structure (default is use all)
        y1 : array[ndep2], optional
            data values being fitted

        Note
        -------
        Only pass y1 if you want to restrict the y-values of extrapolated
        data points. This is useful when solving for the shifts, but should
        not be used when calculating shifted functions for actual use, since
        this restriction can cause discontinuities.
        """
        # Constrained fits may append non-atmospheric quantities to the end of
        # input vector.
        # Extract the output depth scale:
        if ndep is None:
            ndep = len(x1)

        # Shift input x-values.
        # Interpolate input y-values onto output x-values.
        # Shift output y-values.
        x2sh = x2 + par[0]
        y2sh = y2 + par[1]
        y = np.zeros_like(x1)
        y[:ndep] = interp1d(x2sh, y2sh, kind="linear", fill_value="extrapolate")(
            x1[:ndep]
        )  # Note, this implicitly extrapolates

        # Scale output y-values about output y-center.
        ymin = np.min(y[:ndep])
        ymax = np.max(y[:ndep])
        ycen = 0.5 * (ymin + ymax)
        y[:ndep] = ycen + (1.0 + par[2]) * (y[:ndep] - ycen)

        # If extra.y1 was passed, then clip minimum and maximum of output y1.
        if y1 is not None:
            y[:ndep] = np.clip(y[:ndep], np.min(y1), np.max(y1))

        # Set all leftover values to zero
        y[ndep:] = 0
        return y
