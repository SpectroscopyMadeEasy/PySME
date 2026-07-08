Profile-Based NLTE Corrections (Experimental)
==============================================

PySME also contains an experimental **profile-based NLTE correction** layer.
This is separate from the standard `sme.nlte` departure-coefficient workflow.

- `sme.nlte` controls traditional NLTE calculations based on element-specific NLTE grids and departure coefficients.
- `sme.profile_nlte` controls profile corrections that are applied to selected wavelength regions after synthesis.

At the moment this interface is experimental and intentionally conservative.

## Current Scope

The current implementation supports:

- one **element provider** at a time
- a provider may internally cover multiple lines / multiple wavelength windows for that element
- a default provider for hydrogen:
  - `pysme_h_3dnlte_rbf`

The current hydrogen provider is treated as one bundle. It is **not** a per-line toggle.
If the provider is enabled, it operates on its supported Balmer-line windows as one provider.

Multiple profile-NLTE element providers at the same time are currently **not supported**.

## Enable The Feature

Use the `sme.profile_nlte` object:

```python
from pysme.sme import SME_Structure

sme = SME_Structure()
sme.profile_nlte.enabled = True
sme.profile_nlte.element = "H"
```

You may optionally set the provider explicitly:

```python
sme.profile_nlte.provider = "pysme_h_3dnlte_rbf"
```

If `provider` is left as `None`, PySME uses the default provider for the selected element.

## Correction Construction

The hydrogen profile-NLTE correction layer now also exposes how the correction
profile itself is constructed:

```python
sme.profile_nlte.correction_construction = "separate"
```

Supported values are:

- `None`
  - keep the historical default behavior
- `"separate"`
  - construct the correction using the current legacy path
- `"ratio"`
  - construct the correction from the profile ratio on the internal grid

`None` is the default and is intentionally conservative.
At runtime, PySME resolves this setting in the following order:

1. `sme.profile_nlte.correction_construction`
2. environment variable `PYSME_PROFILE_NLTE_H_CORR_MODE`
3. built-in default `separate`

This means existing code keeps the old result unless you opt in explicitly.

## Current Hydrogen Provider

The built-in hydrogen provider currently supports these air-wavelength windows:

- `4335.0 - 4345.0 A`
- `4855.0 - 4868.0 A`
- `6550.0 - 6575.0 A`

These correspond to the current bundled Balmer-line correction windows.

## Runtime Summary

After synthesis, PySME writes the provider status to:

```python
sme.profile_nlte.summary
```

The summary currently includes:

- `requested`
- `applied`
- `element`
- `provider`
- `species`
- `profile_kind`
- `data_key`
- `data_source`
- `parameter_axes`
- `reference_label`
- `citation_info`
- `supported_windows_air`
- `applied_windows_air`
- `fallback`
- `fallback_reason`
- `correction_construction`

Example:

```python
{
    "requested": True,
    "applied": True,
    "element": "H",
    "provider": "pysme_h_3dnlte_rbf",
    "species": "H 1",
    "profile_kind": "intensity_ratio",
    "data_key": "data.hlineprof",
    "data_source": "lineprof.dat",
    "parameter_axes": ["teff", "logg", "monh", "mu"],
    "reference_label": "PySME bundled hydrogen profile dataset",
    "citation_info": None,
    "supported_windows_air": [
        [4335.0, 4345.0],
        [4855.0, 4868.0],
        [6550.0, 6575.0],
    ],
    "applied_windows_air": [
        [6550.0, 6575.0],
    ],
    "fallback": False,
    "fallback_reason": None,
    "correction_construction": "separate",
}
```

This is the recommended way to determine whether a synthesis actually used the profile-based correction.

## Provider Metadata

Each provider is described internally by a manifest-like metadata object.
This provider metadata is what drives:

- the default provider mapping for each element
- supported wavelength windows
- the species check
- summary metadata shown to the user
- future reference / citation bookkeeping

For the current hydrogen provider, PySME records:

- the provider name
- the element and species
- the supported windows
- the profile kind
- the data location key
- the bundled data-source name
- a short reference label
- an optional citation payload

For the current hydrogen provider, the scientific reference is:

- **Amarsi et al. (2018, A&A, 615, A139)**  
  *Effective temperature determinations of late-type stars based on 3D non-LTE Balmer line formation*

The provider summary therefore includes both:

- `reference_label`
- `citation_info`

so the citation can be surfaced directly to users and downstream tooling.

## Fallback Behavior

If the provider is requested but cannot be applied, PySME leaves the synthesis on the standard path and records the reason in `sme.profile_nlte.summary`.

Current fallback reasons include:

- `no_wavelength_overlap`
  - the current synthesis window does not overlap any provider-supported wavelength region
- `no_matching_species_in_linelist`
  - the current linelist does not contain the provider species
- `outside_profile_grid`
  - the stellar parameters fall outside the supported provider grid

When one of these cases occurs, PySME also emits a `WARNING`-level log message.

## Relationship To The Legacy Flag

The old hydrogen-only flag

```python
sme.tdnlte_H = True
```

is still accepted for backward compatibility.

Internally, it now maps onto the new profile-NLTE request path and fills `sme.profile_nlte.summary`.
New code should prefer `sme.profile_nlte`.

## Important Limitations

- This is currently an experimental feature.
- It is not a replacement for `sme.nlte`.
- It is not yet a general multi-element profile-correction framework.
- Only one profile-NLTE element provider can be active at a time.
- The current hydrogen provider is handled as one element-level bundle, not as independently switchable individual lines.
- `get_H_3dnlte_correction_rbf()` currently has its own correction-construction switch.
  It is separate from the normalized-spectrum resampling switch described in the advanced synthesis controls.
