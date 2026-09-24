# Fixed-grid line-centric opacity prototype

## Decision

Park the line-centric rewrite and move the next performance effort to prepared
synthesis / atmosphere-EOS state reuse.

The range-preserving loop transpose is numerically equivalent and accelerates
the native fixed-grid `Transf` by about 1.4--2.3x, but it improves complete
cool-star synthesis wall time by only 1.12x at 200 A and 1.18x at 800 A.  A
more aggressive prototype that discovers the range directly on the output
grid reaches only 1.28x at 10 A for the cool star and changes normalized flux
by up to 8.27e-4.  It therefore fails both the requested wall-time threshold
and the requirement that this be a loop reorder rather than a new accuracy
trade-off.

All experimental C++ paths were removed after measurement.  The production
SMElib path is unchanged by this prototype.

## Inputs and metric

- line list: `merged_3700-9500_hfs.lin`
- center: 5200 A
- fixed-grid sampling: 0.05 A
- widths: 10, 200, and 800 A
- line-list padding: +/-150 A
- input lines: 378,154; 592,757; and 1,258,596 respectively
- models: solar dwarf and 4250 K, log(g)=4.5, [M/H]=+0.3 cool dwarf
- `accrt=1e-4`
- continuum grid: adaptive, 1 A nominal spacing, `rtol=1e-3`
- primary metric: complete `synthesize_spectrum` wall time after VALD parsing

## Existing-path scaling

| model | width | complete wall | native Transf | OPMTRX | OPMTRX / Transf |
|---|---:|---:|---:|---:|---:|
| solar | 10 A | 1.89 s | 0.0084 s | 0.0058 s | 69% |
| solar | 200 A | 3.04 s | 0.1216 s | 0.0958 s | 79% |
| solar | 800 A | 6.97 s | 0.6276 s | 0.5229 s | 83% |
| cool metal-rich | 10 A | 2.81 s | 0.0931 s | 0.0767 s | 82% |
| cool metal-rich | 200 A | 5.62 s | 1.3786 s | 1.3015 s | 94% |
| cool metal-rich | 800 A | 14.06 s | 5.3854 s | 5.1117 s | 95% |

The wide-band opacity kernel does become dominant inside `Transf`, so the
wide-band benchmark was necessary.  However, ALMAX/range work, atmosphere,
EOS, and other setup still dominate enough of the full call that a large
kernel-only speedup is diluted.

## Range-preserving line-centric prototype

The first prototype retained the existing `Wlim_left/right` ranges and only
transposed ordinary-line accumulation from wavelength-depth-line to
line-wavelength-depth.  H, table-driven He, and autoionization remained on the
existing `OPMTRX` path.

| cool model width | reference wall | prototype wall | complete speedup | max abs flux delta |
|---|---:|---:|---:|---:|
| 200 A | 5.771 s | 5.141 s | 1.123x | 7.82e-12 |
| 800 A | 14.781 s | 12.577 s | 1.175x | 1.13e-11 |

The native `Transf` speedup was about 1.7x for the cool model, but the complete
wall-time gain did not reach the 1.3--1.5x go threshold.

## Fully fused on-grid range prototype

The second prototype skipped the old off-grid geometric wing scan for normal
lines.  Starting at each line center, it walked the output grid in both
directions, accumulated opacity, and stopped when the maximum depth-wise
line/continuum ratio fell below `accrt`.

On the 10 A test:

| model | reference wall | fused wall | speedup | max abs flux delta | RMS delta |
|---|---:|---:|---:|---:|---:|
| solar | 1.965 s | 2.016 s | 0.975x | 3.27e-4 | 1.08e-4 |
| cool metal-rich | 3.074 s | 2.406 s | 1.278x | 8.27e-4 | 2.17e-4 |

The flux is systematically higher because many lines individually below the
local cutoff still contribute collectively.  The old range scan brackets a
line with an off-grid probe and retains all fixed-grid samples inside that
range; the fused prototype instead removed each below-threshold on-grid
sample.  Those two operations are not numerically equivalent even with the
same nominal `accrt`.

Preserving the old ranges fixes the spectrum to roundoff, but then the range
probes and the fixed-grid profile evaluations occur at different wavelengths
and cannot actually be reused.  What remains is the range-preserving loop
transpose measured above, whose complete wall-time gain is too small.

## Go / no-go outcome

- cool-star complete speedup >=1.5x: **no**
- wide-band complete speedup >=1.3--1.5x: **no** for the numerically equivalent path
- numerical equivalence at the same `accrt`: **yes** only for the 1.12--1.18x path
- no duplicated special-physics engine: the prototype kept special lines on
  the reference path, but formalizing even this hybrid would add a second
  accumulation engine for a sub-threshold overall gain

The appropriate next target is prepared synthesis / state reuse, where the
measured atmosphere and EOS fractions offer a larger end-to-end opportunity.
