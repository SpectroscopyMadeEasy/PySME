"""Minimal prescribed-field CFLD demo; no stellar atmosphere, EOS or line list.

Requires the companion diagnostic-enabled SMElib and this PySME branch.
Run: python examples/line_depression.py
Use --library to select an isolated diagnostic build explicitly.
"""
import argparse
import json
import os
from pathlib import Path

import numpy as np
from pysme.line_depression import NativeLineDepression


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--library', default=os.environ.get('CFLD_LIBRARY'),
                        help='Diagnostic native library; defaults to the linked PySME library')
    parser.add_argument('--output', type=Path,
                        help='Optional new NPZ file; existing files are refused')
    args = parser.parse_args()
    if args.output and args.output.exists():
        parser.error('Output already exists; choose a new filename')

    library = args.library
    if not library:
        from pysme.sme_synth import ensure_smelib_ready
        library = ensure_smelib_ready().__file__
    native = NativeLineDepression(library)
    # Controlled dimensionless optical-depth problem, not a MARCS atmosphere.
    # Set continuum extinction to unity: m is also the reference continuum tau.
    tau = np.geomspace(1e-7, 50.0, 401)
    continuum_opacity = np.ones_like(tau)
    full_opacity = 3.0 * continuum_opacity  # constant line/continuum ratio = 2
    planck = 1.0 + 0.8 * tau
    result = native.prescribed_path(tau, continuum_opacity, full_opacity,
                                    planck, planck, planck, quadrature=24)
    depression = result.metadata[0] - result.metadata[1]
    reconstructed = result.total.sum() + result.boundary
    relative_error = (reconstructed - depression) / depression
    if not np.isfinite(relative_error) or abs(relative_error) > 1e-6:
        raise RuntimeError(f'Closure failed: relative residual {relative_error}')

    no_line = native.prescribed_path(tau, continuum_opacity, continuum_opacity,
                                     planck, planck, planck, quadrature=24)
    np.testing.assert_allclose(no_line.total, 0.0, rtol=0, atol=1e-10)
    np.testing.assert_allclose(no_line.boundary, 0.0, rtol=0, atol=1e-10)
    print(json.dumps(dict(status='PASS_PRESCRIBED_DEMO',
                         native_intensity_depression=float(depression),
                         integrated_cfld_plus_boundary=float(reconstructed),
                         relative_closure=float(relative_error),
                         no_line_zero=True,
                         scope='Single prescribed ray; no stellar/EOS/scattering calculation'),
                     indent=2))
    if args.output:
        # Density coordinate is reference continuum tau, NOT tauRoss/tau5000.
        with args.output.open('xb') as stream:
            np.savez_compressed(stream, tau_reference_edges=tau,
                                logtau_reference_centres=(np.log10(tau[:-1]) +
                                                          np.log10(tau[1:]))/2,
                                channels=result.channels, metadata=result.metadata,
                                complete_interval_contributions=result.total,
                                complete_per_dex=result.density(np.log10(tau)))
        print(f'Saved {args.output}')


if __name__ == '__main__':
    main()
