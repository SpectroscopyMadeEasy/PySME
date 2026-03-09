# For dev


```{warning}
This is the page mainly for developers of PySME and the note on its function, which may have inaccurate information.
```

## PySME components

Since SME is the C++/Fortran library for the spectral synthesis part, PySME is divided into a few components to make the code work.
1. python package. This is the main interface which the users interact with.
2. `SMR_DLL` class in `sme_synth.py`. This class perform some necessary manipulation for the input and pass it to the functions inside `_smelib`
3. `_smelib` compiled in `src/pysme/smelib`. This is the C-extension which attach the python side data and C++ side data.
4. SMElib core shared library (`sme_synth.*`) which performs the actual synthesis with C++/Fortran.

## What does PySME do from installation to the first use?

Behavior differs by version.

### `<= v0.6.26` (legacy runtime-fetch flow)

1. `pip install pysme-astro`: python package is installed.
2. `import pysme`: PySME may download SMElib binaries at runtime from the SMElib release page.
3. If `_smelib` is missing/incompatible, PySME may try to compile the interface at runtime.

### `>= v0.6.27` (bundled-wheel flow)

1. `pip install pysme-astro`: wheel already contains SMElib shared library and `_smelib`.
2. `import pysme`: PySME loads bundled binaries from the installed package.
3. No runtime download/compile is expected for normal wheel installs.

### Source/editable install (`>= v0.6.27`)

Binary components are built during installation (`pip install -e .` or wheel build), not on first import.

## How to rebuild  smelib

### `>= v0.6.27` (recommended)

From the project root:

```bash
pip install -e . --no-build-isolation
```

This rebuilds/reinstalls package binaries (including SMElib and `_smelib`) for the active environment.

### `<= v0.6.26` (legacy manual workflow)

```bash
cd src/pysme/smelib
python3 setup.py build_ext --inplace
cp _smelib.cpython-310-x86_64-linux-gnu.so /path/to/site-packages/pysme/smelib/
```

## Versions in SME/PySME

There are a few versions in SME/PySME:
- SMElib version
    - The version of the SME library, defined in `SMElib/src/sme/sme_synth_faster.h`, varialbe `VERSION`.
    - This is defined manually.
- SMElib release version
    - Since the SME library here is mainly for the use of PySME, slight modification may be required for the library compilation.
    - The release version is desigend for such small (not effection the synthesize result) function.
    - It will keep the same MAJOR and MINOR versioning, and use PATCH version to indicate the updates.
    - In the ideal situation, SMELIb release version will be the same as SMElib version. 
- PySME version
    - The version of the PySME, which is independent of the versions mentioned above.
    - We will try our best to use the MINOR version to refer to the SMELib version, but large release may take a leap on this.

### Version matching between SMElib and PySME

For the convenience of development, older PySME versions used runtime SMElib download on first import.
From `v0.7.0`, PySME uses strict one-to-one version matching with SMElib at release time.
Each PySME release is pinned to one specific SMElib release, and the mapping is recorded explicitly in the table below.
This avoids accidental runtime mismatch with a newer upstream SMElib release.

For other versions, we will not solve the future compatibility issue.
You can download the corresponding SMElib, and override the default library manually.
PySME versions not listed in the table follow the legacy behavior of resolving to the latest available SMElib release, which is not recommended for reproducible use.

The follwing table shows the version matching between SMElib and PySME.

|PySME version|SMElib release version|SMElib version|
|:--:|:--:|:--:|
|v0.7.0|v6.13.16|6.13 (June 2025)|
|v0.6.23|v6.13.12|6.13 (June 2025)|
|v0.4.199|v6.0.6|6.03 (July 2019)|
|v0.4.167-v0.4.198|v6.0.6|6.03 (July 2019)|

The PySME version range indicate the versions which manually matchting of the SMELib is required, and the single PySME version indicates the one with freezed SMElib (thus only pip install is required before using).
Note that:
- Small bug fix on SMElib will be labeled with the PATCH version, thus the actual SMElib used would be 6.13.x. The change of SMElib MAJOR and MINOR version will follow that in SME. 
- The support for Apple Silicon Mac is not complete for PySME v0.4. You need to clone the SMELib from the source code, compile and replace the library files by yourself.

### How to update PySME/SMElib version

The versions are controlled by `git tag`. 
For PySME, setting a new tag will update the version to the latest tag.
For SMElib, setting a new release will based on a tag. 

### Docs versioning policy (RTD)

PySME docs use a tag-driven versioning strategy.

1. Create release tags in `vX.Y.Z` format (for example `v0.6.23`).
2. `docs/conf.py` resolves `release` from `src/pysme/_version.py` (versioneer), so docs title/version follow git tags.
3. `.readthedocs.yaml` installs both docs requirements and the project itself, so tag metadata is available during RTD build.

RTD project settings must also be configured once (in the RTD web UI):

- Keep both `latest` (branch docs) and `stable` (latest release tag) enabled.
- Set default docs version to `stable`.
- Add an automation rule to activate SemVer tags (`v*`) as RTD versions.
- Hide or deactivate outdated branch versions if needed.

## NLTE departure coefficients

What happens when the NLTE grid is added into PySME?

1. `sme.nlte.set_nlte()`, then nothing happens.
2. Inside `synthesize_spectrum()`, `sme.nlte.update_coefficients(sme, dll, self.lfs_nlte)` will be triggered if there is nlte grids in `sme`.

### `nlte.update_coefficients`

1. All the $b$s will be reset by `dll.ResetDepartureCoefficients()`.
2. Format of line list will be check, and NLTE calculation will not be performed if the format is not `long`.
3. Get the NLTE grid using `self.get_grid(sme, elem, lfs_nlte)`
4. If no lines are found for this element, remove it from NLTE list.
5. Get the $b$ matrix using `grid.get(sme.abund, sme.teff, sme.logg, sme.monh, sme.atmo)`. 
    - This is also the core part of `update_coefficients`.
6. Input the $b$s using `dll.InputDepartureCoefficients(bmat[:, lr], li)`.

#### `self.get_grid`

This function is to get the $b$ grid for a specific element.
Note that the whole $b$ grid is not readed, but only the levels related with the input line list.
Thus the configuration or level matching is done inside this function.

#### configuration/level matching

This is done inside `Grid.__init__`. The level matching can be done with levels or energy, and the default one is energy.

#### `grid.get`

1. Find the relative abundance input
2. `self.read_grid(rabund, teff, logg, monh)`
3. `return self.interpolate(rabund, teff, logg, monh, atmo)`.

## `apply_continuum`

Apply continuum to the spectra according to `cscale`.
    
For each segment i: 
    - If `cscale_type` contains 'spline', then `cscale` should be an array which with the same length of sme.wave[i], and its value will be multiplied to sme.synth[i].
        - If the length of `cscale` is not the same as the length of `sme.wave[i]`, then it will be interpolated according to `cwave`, which for now is the same as `sme.wave[i]` thus this situation is not allowed.
        - We may just remove the `spline` option, if not stated in the paper. 
    - If `cscale_type` does not contain 'spline', then the code will calculate `np.polyval` according to `cscale`, and apply to `sme.synth[i]`.

## `radial_velocity_mode`

This variable determines how the cscale and vrad being assigned, and can be specfied manually.

## NLTE grid

If you create your own NLTE departure coefficient grid following the current public grid to a `DirectAccess` file and you use the `nlte.DirectAccess.wrtie` function, note that the `wrtie` and `read` code changes the shape of b-grid, tau and rhox. 
Some extra test is needed to clear the situation. 

## SME

### `Trasf` function

1. `AutoIonization`
2. Calculate Line center opacity using `LINEOPAC`.
3. Calculate Line contribution limits using `OPMTRX`.
    - Step 2 and 3 go through all the input lines

## IDLSME

### Installation

TBD.
