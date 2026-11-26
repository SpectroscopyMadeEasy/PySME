"""Experimental native physical line-depression diagnostics.

No atmospheric model, EOS or transfer is started on import. Native state calls
must run under SME_DLL's process lock after the ordinary transfer preparation.
"""
import ctypes
import operator
import numpy as np

CHANNELS = ("line_driving", "source_plus_numerical", "continuum_emission",
            "full_emission", "source_correction", "numerical_correction")
PTR = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")


def _vector(values, name):
    a = np.ascontiguousarray(values, dtype=np.float64)
    if a.ndim != 1 or not len(a) or not np.isfinite(a).all():
        raise ValueError(f"{name} must be a nonempty finite vector")
    return a


def _quadrature(value):
    value = operator.index(value)
    if not 4 <= value <= 128:
        raise ValueError("quadrature must be between 4 and 128")
    return value


class LineDepressionResult:
    """Absolute native intensity contributions, channel axis penultimate.

    Path arrays: (6, npath-1); state arrays: (nwave, nray, 6, nradial-1).
    Metadata ends with Ic, If, attenuated continuum/full boundaries, physical
    depression boundary, minimum internal Ic/If; state adds path count/tangent.
    Ordinary native outer-interval omission is preserved.
    """
    channel_names = CHANNELS

    def __init__(self, channels, metadata):
        self.channels = channels
        self.metadata = metadata

    @property
    def total(self):
        return self.channels[..., 0, :] + self.channels[..., 4, :] + self.channels[..., 5, :]

    @property
    def boundary(self):
        return self.metadata[..., 4]

    def density(self, edges):
        """Convert shell integrals using this coordinate's own monotonic widths."""
        edges = _vector(edges, "edges")
        widths = np.diff(edges)
        if len(widths) != self.channels.shape[-1] or not (np.all(widths > 0) or np.all(widths < 0)):
            raise ValueError("edges must be strictly monotonic and match shell count")
        return self.total / abs(widths)

    def disc_flux(self, weights):
        """Return absolute outer-radius flux bins/boundary using global area weights.

        Chunk weights need not sum to one; sum returned chunk fluxes afterwards.
        No rotation, macroturbulence or instrumental convolution is included.
        """
        w = _vector(weights, "weights")
        if self.channels.ndim != 4 or len(w) != self.channels.shape[1] or np.any(w < 0):
            raise ValueError("Nonnegative projected-area weights must match rays")
        return (np.pi*np.einsum("r,wri->wi", w, self.total),
                np.pi*np.einsum("r,wr->w", w, self.boundary))


class NativeLineDepression:
    """Bind diagnostic symbols from a library or the loaded Python extension.

    Prefer SME_DLL.GetLineDepressionContributionfunction for atmospheric state:
    it resolves symbols through the already loaded extension, avoiding a second
    SMElib state. Direct binding is useful for prescribed-opacity benchmarks.
    """
    def __init__(self, library):
        self.library = ctypes.CDLL(str(library))
        try:
            self.path = self.library.SMEPhysicalCFLDPath
            self.state = self.library.SMEPhysicalCFLD
        except AttributeError as exc:
            raise RuntimeError("Loaded SMElib lacks the optional physical CFLD extension") from exc
        self.path.argtypes = [ctypes.c_int] + [PTR]*6 + [ctypes.c_int, ctypes.c_int, PTR, PTR]
        self.state.argtypes = [ctypes.c_int, PTR, ctypes.c_int, PTR, ctypes.c_int, ctypes.c_int, PTR, PTR]
        self.path.restype = self.state.restype = ctypes.c_int

    def prescribed_path(self, mass, continuum_opacity, full_opacity,
                        continuum_source, full_source, planck, tangent=False, quadrature=24):
        arrays = [_vector(v, "path field") for v in (mass, continuum_opacity,
                  full_opacity, continuum_source, full_source, planck)]
        n = len(arrays[0])
        if n < 3 or any(len(v) != n for v in arrays) or np.any(np.diff(arrays[0]) <= 0):
            raise ValueError("At least three increasing path nodes and equal field lengths required")
        if np.any(arrays[1] <= 0) or np.any(arrays[2] <= 0):
            raise ValueError("Extinctions must be positive")
        channels, meta = np.zeros((6, n-1)), np.zeros(7)
        rc = self.path(n, *arrays, int(bool(tangent)), _quadrature(quadrature), channels, meta)
        if rc:
            raise RuntimeError(f"Native physical CFLD status {rc}")
        return LineDepressionResult(channels, meta)

    def spherical_state(self, wave, mu, layers, quadrature=24):
        wave, mu = _vector(wave, "wave"), _vector(mu, "mu")
        layers = operator.index(layers)
        if layers < 3 or np.any(wave <= 0) or np.any(np.diff(wave) <= 0):
            raise ValueError("Positive increasing wavelengths and at least three layers required")
        if np.any(mu <= 0) or np.any(mu > 1):
            raise ValueError("Ray mu must satisfy 0 < mu <= 1")
        channels = np.zeros((len(wave), len(mu), 6, layers-1))
        meta = np.zeros((len(wave), len(mu), 9))
        rc = self.state(len(wave), wave, len(mu), mu, layers, _quadrature(quadrature), channels, meta)
        if rc:
            raise RuntimeError(f"Native physical CFLD status {rc}; check prepared spherical state, ray limit and dimensions")
        return LineDepressionResult(channels, meta)
