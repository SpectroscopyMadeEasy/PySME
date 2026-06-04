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
- `supported_windows_air`
- `applied_windows_air`
- `fallback`
- `fallback_reason`

Example:

```python
{
    "requested": True,
    "applied": True,
    "element": "H",
    "provider": "pysme_h_3dnlte_rbf",
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
}
```

This is the recommended way to determine whether a synthesis actually used the profile-based correction.

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
