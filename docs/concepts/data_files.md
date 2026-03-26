Data Files You Need
===================

PySME uses local cache data under `~/.sme/`.

The supported way to change data/cache locations is to edit `~/.sme/config.json`.
PySME reads directory settings from that file.

Typical files:

- Atmosphere grids
- NLTE grids
- Optional custom resources

## Configure data locations

Edit `~/.sme/config.json` and set the data paths you want:

```json
{
  "data.file_server": "https://sme.astro.uu.se/atmos",
  "data.file_servers": [
    "https://sme.astro.uu.se/atmos"
  ],
  "data.hlineprof": "/path/to/pysme/hlineprof",
  "data.atmospheres": "/path/to/pysme/atmospheres",
  "data.nlte_grids": "/path/to/pysme/nlte_grids",
  "data.cache.atmospheres": "/path/to/pysme/atmospheres/cache",
  "data.cache.nlte_grids": "/path/to/pysme/nlte_grids/cache",
  "data.pointers.atmospheres": "datafiles_atmospheres.json",
  "data.pointers.nlte_grids": "datafiles_nlte.json"
}
```

Notes:

- `data.file_servers` is optional and defines an ordered list of path-compatible mirrors.
- If `data.file_servers` is present, PySME expands relative pointer paths against those servers and falls back automatically.
- If `data.file_servers` is absent, PySME uses legacy `data.file_server` behavior.
- Pointer values in `datafiles_*.json` can be:
  - a single relative path (expanded against each path-compatible mirror server)
  - a list of paths/URLs (tried in order)
  - a full URL (`https://...` or `file://...`)
- Absolute URLs are tried as written and are not joined with `data.file_servers`.
- This allows mixed-source fallbacks such as `NADC -> Uppsala -> Zenodo` in a single pointer list.
- In the current PySME defaults, Uppsala is configured as the path-compatible mirror in `data.file_servers`, while NADC and Zenodo are usually referenced as explicit URLs in the pointer files.
- The NADC packaged data mirror page is <https://nadc.china-vo.org/res/r101793/>.
- `~` is supported in path values and will be expanded to your home directory.
- Changing these paths does not migrate old files automatically; move existing files manually if needed.

```{admonition} Accessing data files
The atmosphere and nlte data files should be downloaded from the server automatically when used, so network connection is required when using PySME (not only during installation).
These files can also be downloaded manually by:
- Download data files as part of IDL SME from [here](http://www.stsci.edu/~valenti/sme.html).
- copy them into their respective storage locations in ~/.sme/atmospheres and ~/.sme/nlte_grids
  - atmospheres: everything from SME/atmospheres
- Download the nlte_grids in [zenodo](https://doi.org/10.5281/zenodo.3888393).
- The files (mainly atmosphere models and NLTE departure coefficient grids) required by PySME will be saved inside `~/.sme/`. These files can be large thus if your home directory is small, we recommend to create a softlink for `~/.sme`.
```
