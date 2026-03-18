# FAQ


## How do I change the default log file?

Call `util.start_logging(filename)`.

```py
from pysme import util
util.start_logging("your_log_file.log")
```

## Get output of `A module that was compiled using NumPy 1.x cannot be run in NumPy 2.2.0 as it may crash.`

This is an ABI mismatch between your NumPy and the installed `_smelib` binary.

- For `<= v0.6.26`, PySME may attempt a runtime rebuild/download path.
- For `>= v0.6.27`, the recommended fix is to reinstall PySME in the target environment so binaries are rebuilt/reinstalled consistently.

If you see output like:
```
running build_ext
building '_smelib' extension
```
the extension is being rebuilt.

## I get an error "Derivatives in the starting point are not finite"

Make sure your initial stellar parameters are within the
atmosphere grid defined by the atmosphere file set in sme.atmo.source

## I get an error "lnGAS: DGESVX failed to solved for corrections the partial pressures."

The most possible reason would be the abvundance of the element in error
is too low or nan, thus the EOS code cannot compute its EOS. 

## I upgraded PySME, but behavior/version output still looks old

This is usually caused by a stale compiled extension (`_smelib`) left in the environment.
Python may load the stale binary instead of the expected one from the current wheel.

First, check what is actually loaded:

```bash
python - <<'PY'
import pysme
from pysme.smelib import _smelib
print("PySME:", pysme.__version__)
print("PySME file:", pysme.__file__)
print("_smelib file:", _smelib.__file__)
PY
```

If `_smelib` points to an unexpected legacy file, do a clean reinstall in the same environment:

```bash
python -m pip uninstall -y pysme-astro
python - <<'PY'
import site, shutil
from pathlib import Path
for sp in site.getsitepackages():
    p = Path(sp)
    shutil.rmtree(p / "pysme", ignore_errors=True)
    for di in p.glob("pysme_astro-*.dist-info"):
        shutil.rmtree(di, ignore_errors=True)
print("cleaned")
PY
python -m pip install --no-cache-dir pysme-astro==X.Y.Z  # example version
```

Replace `X.Y.Z` with the target release version you want to test.
