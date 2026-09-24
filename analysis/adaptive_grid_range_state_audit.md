# Legacy adaptive-grid range-state residual audit

## Conclusion

The small solar and giant legacy-versus-batched residuals are fully explained by a separate legacy RKINTS range-state mutation. They are not caused by generation batching, wavelength interpolation, or second weak-line pruning.

With second pruning disabled in both paths:

- solar legacy RKINTS shortens 33 precomputed line ranges and differs from batched transfer by `1.8470e-5` maximum;
- giant legacy RKINTS shortens three precomputed line ranges and differs by `1.9150e-7` maximum;
- preventing only those range writes makes the legacy and batched wavelength arrays and normalized fluxes exactly equal at every node for both models.

The production batched path should therefore keep the precomputed ALMAX/CDR physical ranges immutable. The historical RKINTS path remains unchanged for compatibility.

## Source mechanism

During the legacy seed scan, RKINTS inserts a midpoint immediately before each accepted line centre. After evaluating that midpoint it executes range-discovery-era logic equivalent to:

```cpp
if (Wlim_left[line] < midpoint && WLCENT[line] > midpoint &&
    ALMAX[line] < accrt)
    Wlim_left[line] = midpoint;
```

This was meaningful when RKINTS participated in discovering a line's range. It is not valid after `ALMAXRange` has already supplied a physical range: the midpoint happens to lie in a weak part of the profile, so the already-certified blue support is overwritten by a grid-layout-dependent value.

The mutation changes opacity at later refinement nodes but does not necessarily change the chosen grid. That is exactly what is observed here: all wavelength arrays remain identical, while flux changes at the affected line wings.

## Controlled experiment

Three calculations were made with the same ALMAX mask and physical ranges:

1. production generation-batched transfer, whose active mask and ranges are immutable;
2. legacy RKINTS with second pruning disabled and historical range writes enabled;
3. legacy RKINTS with second pruning disabled and only the seed-scan range writes disabled through a temporary diagnostic switch.

All other opacity, source-function, RT, seed-grid, minimum-spacing, and `accwi` logic remained unchanged.

| Model | Batched points | Legacy points | Mutated ranges | Legacy mutable − batch max | Legacy immutable − batch max |
|---|---:|---:|---:|---:|---:|
| Solar dwarf | 2,447 | 2,447 | 33 | `1.8470e-5` | **0** |
| Giant | 9,933 | 9,933 | 3 | `1.9150e-7` | **0** |

For both models, the grids have zero Hausdorff distance in both legacy variants. With immutable ranges, the spectra also have zero maximum, RMS, mean, and p99 difference.

Example solar mutations include:

- W I 5192.706 Å: `[5192.406, 5193.006]` becomes `[5192.658, 5193.006]`;
- TiO 5194.0364 Å: `[5193.7364, 5194.3364]` becomes `[5194.02795, 5194.3364]`;
- TiO 5195.2098 Å: `[5194.9098, 5195.5098]` becomes `[5195.1695, 5195.5098]`.

All three giant changes are TiO blue-range truncations of the same form.

## Runtime context

The production batched calls took 0.157 s for the Sun and 4.35 s for the giant. The two legacy variants took approximately 1.51 s and 103 s respectively. Disabling the range writes does not improve legacy runtime; the speedup remains attributable to evaluating sorted generations through the interval index.

## Production decision

- Do not reproduce this mutation in the new batched path.
- Do not remove it from legacy RKINTS in the same change; that path remains an explicit compatibility/reference implementation.
- Treat precomputed ALMAX/CDR ranges as immutable synthesis input.
- If legacy internal line selection is later migrated to batching, decide separately whether it still needs a true range-discovery phase before transfer.

## Artifacts

- `adaptive_grid_range_state_audit.py`
- `adaptive_grid_range_state_audit.json`
- the temporary `SME_RKINTS_RANGE_STATE` diagnostic switch was removed after
  this audit; the production batched path enforces immutable ranges directly
