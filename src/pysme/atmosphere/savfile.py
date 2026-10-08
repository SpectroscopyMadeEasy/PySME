# -*- coding: utf-8 -*-
import os
from os.path import basename
from tempfile import NamedTemporaryFile

import numpy as np
from scipy.io import readsav

from .atmosphere import AtmosphereGrid


class SavFile(AtmosphereGrid):
    """IDL savefile atmosphere grid"""

    _cache = {}

    def _identify_abund_format(self):
        """Recover legacy He metadata, including old converted numpy caches.

        With H present, log10(N_He/N_all) is negative; a nonnegative He
        fraction belongs to the legacy format. Never reinterpret general
        Abund inputs, or explicitly labelled non-SME grids.
        """
        if self.abund_format != "sme":
            return
        he = self["abund"][:, 1]
        finite = he[np.isfinite(he)]
        if finite.size and np.any(finite >= 0):
            if np.any(finite < 0) or np.any(finite >= 1):
                raise ValueError("Atmosphere grid has inconsistent helium abundance formats")
            self.abund_format = "sme_legacy"

    def __new__(cls, filename, source=None, lfs=None):
        filename = os.fspath(filename)
        # Check if the file is in the cache
        if filename in cls._cache:
            return cls._cache[filename]
        # Try loading the datafile using Numpy which is faster
        # and was generated in a previous iteration of PySME
        try:
            self = cls.load(filename)
        except (OSError, ValueError, KeyError):
            self = None
        if self is not None:
            self._identify_abund_format()
            cls._cache[filename] = self
            return self

        # Otherwise we parse the sav file
        data = readsav(filename)

        npoints = data["atmo_grid_maxdep"]
        ngrids = data["atmo_grid_natmo"]
        self = super(SavFile, cls).__new__(cls, ngrids, npoints)
        if source is None:
            source = basename(filename)
        self.source = source

        try:
            self.info = b" ".join(data["atmo_grid_intro"]).decode()
        except (KeyError, TypeError):
            self.info = ""

        # TODO cover all cases
        if "marcs" in self.source:
            self.citation_info = r"""
                @ARTICLE{2008A&A...486..951G,
                    author = {{Gustafsson}, B. and {Edvardsson}, B. and {Eriksson}, K. and
                    {J{\o}rgensen}, U.~G. and {Nordlund}, {\r{A}}. and {Plez}, B.},
                    title = "{A grid of MARCS model atmospheres for late-type stars. I. Methods and general properties}",
                    journal = {Astronomy and Astrophysics},
                    keywords = {stars: atmospheres, Sun: abundances, stars: fundamental parameters, stars: general, stars: late-type, stars: supergiants, Astrophysics},
                    year = "2008",
                    month = "Aug",
                    volume = {486},
                    number = {3},
                    pages = {951-970},
                    doi = {10.1051/0004-6361:200809724},
                    archivePrefix = {arXiv},
                    eprint = {0805.0554},
                    primaryClass = {astro-ph},
                    adsurl = {https://ui.adsabs.harvard.edu/abs/2008A&A...486..951G},
                    adsnote = {Provided by the SAO/NASA Astrophysics Data System}
                }
            """
        elif "atlas" in self.source:
            self.citation_info = r"""
                @MISC{2017ascl.soft10017K,
                    author = {{Kurucz}, Robert L.},
                    title = "{ATLAS9: Model atmosphere program with opacity distribution functions}",
                    keywords = {Software},
                    year = "2017",
                    month = "Oct",
                    eid = {ascl:1710.017},
                    pages = {ascl:1710.017},
                    archivePrefix = {ascl},
                    eprint = {1710.017},
                    adsurl = {https://ui.adsabs.harvard.edu/abs/2017ascl.soft10017K},
                    adsnote = {Provided by the SAO/NASA Astrophysics Data System}
                }
                @ARTICLE{2002A&A...392..619H,
                    author = {{Heiter}, U. and {Kupka}, F. and {van't Veer-Menneret}, C. and {Barban}, C. and {Weiss}, W.~W. and {Goupil}, M. -J. and {Schmidt}, W. and {Katz}, D. and {Garrido}, R.},
                    title = "{New grids of ATLAS9 atmospheres I: Influence of convection treatments on model structure and on observable quantities}",
                    journal = {\aap},
                    keywords = {stars: atmospheres, stars: fundamental parameters, stars: variables delta Scuti stars, convection, Astrophysics},
                    year = 2002,
                    month = sep,
                    volume = {392},
                    pages = {619-636},
                    doi = {10.1051/0004-6361:20020788},
                    archivePrefix = {arXiv},
                    eprint = {astro-ph/0206156},
                    primaryClass = {astro-ph},
                    adsurl = {https://ui.adsabs.harvard.edu/abs/2002A&A...392..619H},
                    adsnote = {Provided by the SAO/NASA Astrophysics Data System}
                }
            """
        elif "ll" in self.source:
            self.citation_info = r"""
                @ARTICLE{2004A&A...428..993S,
                    author = {{Shulyak}, D. and {Tsymbal}, V. and {Ryabchikova}, T. and {St{\"u}tz}, Ch. and {Weiss}, W.~W.},
                    title = "{Line-by-line opacity stellar model atmospheres}",
                    journal = {\aap},
                    keywords = {stars: atmospheres, stars: abundances, stars: chemically peculiar, stars: fundamental parameters, stars: individual: CU Vir, stars: individual: HD 124224},
                    year = 2004,
                    month = dec,
                    volume = {428},
                    pages = {993-1000},
                    doi = {10.1051/0004-6361:20034169},
                    adsurl = {https://ui.adsabs.harvard.edu/abs/2004A&A...428..993S},
                    adsnote = {Provided by the SAO/NASA Astrophysics Data System}
                }
            """
        else:
            self.citation_info = ""  # ???

        atmo_grid = data["atmo_grid"]

        if "RADIUS" in atmo_grid.dtype.names and "HEIGHT" in atmo_grid.dtype.names:
            self.geom = "SPH"
            self["radius"] = atmo_grid["radius"]
            self["height"] = np.stack(atmo_grid["height"])
            # If the radius is given in absolute values
            # self["radius"] /= np.max(self["radius"])
        else:
            self.geom = "PP"
            self["radius"][:] = 1

        self.abund_format = "sme"

        # Scalar Parameters (one per atmosphere)
        self["teff"] = atmo_grid["teff"]
        self["logg"] = atmo_grid["logg"]
        self["monh"] = atmo_grid["monh"]
        self["vturb"] = atmo_grid["vturb"]
        self["lonh"] = atmo_grid["lonh"]
        self["wlstd"] = atmo_grid["wlstd"]
        # Vector Parameters (one array per atmosphere)
        self["rhox"] = np.stack(atmo_grid["rhox"])
        self["tau"] = np.stack(atmo_grid["tau"])
        self["temp"] = np.stack(atmo_grid["temp"])
        self["rho"] = np.stack(atmo_grid["rho"])
        self["xne"] = np.stack(atmo_grid["xne"])
        self["xna"] = np.stack(atmo_grid["xna"])
        self["abund"] = np.stack(atmo_grid["abund"])
        self["opflag"] = np.stack(atmo_grid["opflag"])
        self._identify_abund_format()

        # Store in cache
        cls._cache[filename] = self
        # And also replace the IDL file with a numpy file in the cache
        # We have to use a try except block, as this will crash with
        # permissions denied on windows, when trying to copy an open file
        # here the temporary file
        # Therefore we close the file, after copying and then delete it manually
        if lfs is not None:
            try:
                with NamedTemporaryFile(delete=False) as named:
                    self.save(named)
                    named.flush()
                lfs.move_to_cache(named.name, key=lfs.get_url(self.source))
            finally:
                try:
                    os.remove(named.name)
                except:
                    pass

        return self
