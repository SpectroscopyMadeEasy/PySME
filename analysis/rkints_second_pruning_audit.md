# Audit of the legacy RKINTS second weak-line pruning

## Executive conclusion

**Recommendation: remove the sequential `MARK=2` weak-line deactivation from the new generation-batched RKINTS path.** Keep the historical RKINTS path temporarily as a compatibility/reference implementation.

The second pruning is a historical performance approximation, not an accuracy safeguard. In the decisive metal-poor tests it removes lines with substantial independent contributions, moves the result away from a minimally pruned reference, and can make the final active set depend on the order of the pruning decisions. Its measured production-style benefit on the same batched grid is 6.1 ms, or 13.4%, in the one affected main case; it has zero semantic benefit in the other three main models. That saving is small compared with the already demonstrated 5.8–24.2× gain from generation batching.

No production default was changed during the audit itself. Its temporary C++
environment controls and logging hooks were removed after the production
decision was implemented.

## 1. Problem

Legacy `RKINTS` first receives an active line set selected by `ALMAX`, then scans line centres from blue to red. At each eligible centre it performs a disk-centre radiative-transfer calculation and applies:

```text
depression = 1 - I(mu=1, line centre) / I_cont(mu=1, line centre)
depression < accwi  ->  MARK[line] = 2
```

For this audit, `accrt=1e-4` and `accwi=0.003`. The question is whether this second, emergent-intensity-based deletion should survive in the batched implementation.

## 2. Three line-selection/truncation concepts

These mechanisms have different responsibilities and should not be conflated:

1. **Physical/validity range (`accrt`)** determines the wavelengths over which an already-active line still has a sufficient opacity contribution. It truncates support, not line membership.
2. **ALMAX preselection** deterministically decides which lines enter the synthesis active set from a precomputed strength proxy.
3. **RKINTS second pruning (`accwi`)** evaluates the *total blended disk-centre intensity* at a line centre and permanently changes that line's `MARK` state for subsequent work.

This audit concerns only (3). Physical line ranges remain required.

## 3. Legacy `MARK=2` semantics

The centre calculation is in `RKINTS`: `OPMTRX` builds total and continuum quantities, `TBINTG` obtains the disk-centre intensity, and a depression below `EPS2/accwi` sets `MARK[line]=2`.

`OPMTRX` subsequently rejects a marked line before profile, opacity, or NLTE processing:

```cpp
if (MARK[LINE] || WAVE <= Wlim_left[LINE] || WAVE >= Wlim_right[LINE])
    continue;
```

Therefore `MARK=2` affects:

- line opacity and, for NLTE lines, line emissivity/source contribution;
- all later emergent intensities, including later line-centre pruning decisions;
- adaptive-grid error estimates and refinement, because the spectrum being sampled has changed;
- later range/grid bookkeeping and the final returned spectrum.

The measured quantity is neither isolated line depth nor a line's marginal contribution. It is the total blended spectrum depression in the *current mutable active state*, divided by an exact centre continuum intensity from the same call. Since earlier decisions change that state, the rule is structurally order-dependent.

## 4. Benchmark setup

The LTE audit used the supplied merged line list, extracted to the existing 10 Å plus approximately ±150 Å working asset, and four atmospheres:

| Model | Teff / log g / [M/H] / vmic | Main window |
|---|---:|---:|
| Solar dwarf | 5777 / 4.44 / 0.0 / 1.0 | 5195–5205 Å |
| Metal-poor dwarf | 6000 / 4.0 / -2.0 / 1.2 | 5195–5205 Å |
| Cool metal-rich dwarf | 4250 / 4.5 / +0.3 / 1.0 | 5195–5205 Å |
| Giant | 4500 / 2.0 / 0.0 / 1.5 | 5195–5205 Å |

Two additional metal-poor windows were run: Mg b/strong-line 5165–5175 Å and a relatively clean atomic window at 5300–5310 Å.

Four semantics were evaluated:

- **L historical:** ALMAX plus sequential RKINTS pruning on the legacy adaptive grid.
- **L final mask:** the final historical `MARK` mask replayed on the same fine fixed grid as the other modes; this isolates pruning from grid-generation differences.
- **A:** production ALMAX at `1e-4`, no second pruning.
- **P:** permissive ALMAX at `1e-5`, no second pruning.
- **F:** minimally pruned reference at `1e-6`, no second pruning.

P and F use a common physical range floor of `1e-6`; thus their difference from A changes line membership without relaxing the reference support calculation. All controlled accuracy comparisons use a 0.000625 Å fixed grid. For the decisive metal-poor case, tightening F from `1e-6` (14,505 lines) to `1e-7` (16,207 lines) changed normalized flux by only `3.16e-6` maximum and `1.51e-6` RMS, so F is converged well below the effects being judged.

## 5. How many lines are additionally removed?

Counts are global active-line counts; overlap and tested-centre counts refer to the padded transfer window.

| Model/window | Legal input | After ALMAX | Overlapping ALMAX | Centres tested | Additional `MARK=2` | Of overlap | Of tested |
|---|---:|---:|---:|---:|---:|---:|---:|
| Solar 5195–5205 | 378,154 | 15,332 | 886 | 475 | 0 | 0% | 0% |
| Metal-poor 5195–5205 | 378,154 | 6,853 | 406 | 243 | 179 | 44.1% | 73.7% |
| Cool metal-rich 5195–5205 | 378,154 | 377,405 | 21,034 | 2,512 | 0 | 0% | 0% |
| Giant 5195–5205 | 378,154 | 346,712 | 19,463 | 2,483 | 0 | 0% | 0% |
| Metal-poor Mg b | 378,154 | 6,853 | 420 | 297 | 209 | 49.8% | 70.4% |
| Metal-poor clean atomic | 378,154 | 6,853 | 310 | 202 | 171 | 55.2% | 84.7% |

The effect is highly atmosphere-dependent. It does nothing in three main models, but removes roughly half of the locally relevant ALMAX lines in all three tested metal-poor windows.

## 6. Properties of removed lines

In the main metal-poor window, the 179 removed lines comprise 100 Fe-group lines, 28 MgH lines, 50 other atomic lines, and one other molecular line. Their ALMAX distribution is:

| ALMAX bin | Removed lines |
|---|---:|
| `1e-4–1e-3` | 109 |
| `1e-3–1e-2` | 56 |
| `1e-2–1e-1` | 14 |

Thus 70/179 are above `1e-3`, and 14 are above `1e-2`; these are not exclusively threshold-borderline lines. The largest removed ALMAX is `0.0873` for Nd III 5203.9236 Å, whose legacy blended centre depression is only `5.73e-4`. Ti I 5201.0814 Å has ALMAX `0.0173` and a legacy centre depression `1.63e-3`, yet removing it alone changes normalized flux by `1.84e-3`.

The clean atomic window is still more striking: one removed Nd II line has ALMAX `0.205`. This proves that a shallow *blended emergent depression* is not a safe proxy for negligible line opacity or negligible marginal spectral effect.

Complete rows in the JSON include species, wavelength, ALMAX, `[Wlim_left, Wlim_right]`, centre continuum intensity, blended disk-centre intensity/depression, compact line index, and decision sequence.

## 7. Order-dependence experiment

A diagnostic scheduler changed only the order of the `MARK=2` decisions; `OPMTRX` continued to accumulate opacity in the original line-list order. This avoids violating SMElib's wavelength-order assumptions and removes floating-point accumulation order as an explanation.

| Window | Forward pruned | Reverse pruned | Changed membership | Max spectral difference | RMS | p99 |
|---|---:|---:|---:|---:|---:|---:|
| 5195–5205 | 179 | 179 | 0 | 0 | 0 | 0 |
| Mg b | 209 | 209 | 0 | 0 | 0 | 0 |
| Clean atomic | 171 | 169 | 2 | `3.22e-3` | `1.66e-4` | `4.21e-5` |

The two forward-only deletions are Fe I 5300.4024 Å (ALMAX `0.01798`) and Co I 5301.0229 Å (ALMAX `6.27e-4`). Their forward blended depressions are `2.985e-3` and `2.896e-3`, immediately below `accwi=0.003`. Changing earlier deletions moves them across the discontinuous threshold.

Therefore the legacy rule is not merely theoretically order-dependent; a physically unchanged line set produces a different final active set and a maximum `3.22e-3` flux difference solely from decision order.

## 8. Legacy versus ALMAX-only

On the historical adaptive grid, metal-poor L minus generation-batched A reaches `1.450e-3` intrinsic, `2.785e-4` at R=20,000, and `5.993e-4` at R=60,000. The difference is positive at every sampled nonzero point: pruning produces higher, shallower flux, as expected for removed absorption opacity.

For the other main models, historical L versus batched A is zero for the cool model, `1.92e-7` maximum for the giant, and `1.85e-5` for the Sun. Since none of these models pruned any line, those small residuals cannot be attributed to the second pruning. They arise from other historical adaptive-grid/range-state semantics and should be audited separately.

That separate audit is now complete: legacy RKINTS shortened 33 precomputed
line ranges in the solar case and three in the giant case during its seed scan.
Disabling only those writes makes legacy and batched wavelengths and fluxes
exactly identical. See `adaptive_grid_range_state_audit.md`.

## 9. High-accuracy/full reference

The scientifically controlled comparison is the same fixed fine grid, where only the active mask differs:

| Metal-poor 5195–5205, relative to F | Max | RMS | R=20k max | R=60k max |
|---|---:|---:|---:|---:|
| L final mask | `1.897e-3` | `1.883e-4` | `5.203e-4` | `1.097e-3` |
| A, ALMAX-only | `1.416e-4` | `5.142e-5` | `8.143e-5` | `1.091e-4` |
| P, `1e-5` | `4.883e-6` | `6.304e-7` | `1.358e-6` | `2.622e-6` |

A is closer than L by factors of 13.4 in maximum error, 3.7 in RMS, 6.4 at R=20k, and 10.1 at R=60k.

The result repeats in both extra windows:

| Window | L max vs F | A max vs F | L R=20k / R=60k | A R=20k / R=60k |
|---|---:|---:|---:|---:|
| Mg b | `8.47e-4` | `7.93e-5` | `3.47e-4 / 4.67e-4` | `5.77e-5 / 7.08e-5` |
| Clean atomic | `3.29e-3` | `8.09e-5` | `6.96e-4 / 1.75e-3` | `3.00e-5 / 5.27e-5` |

P is within `9.05e-6` of F in Mg b and `6.53e-6` in the clean window. There is no tested accuracy case in which the second pruning improves agreement with the minimally pruned reference.

## 10. Flux-difference direction

For the main metal-poor fixed-grid comparison, `F_L - F_A` is nonnegative at all points (strictly positive at 99.58% of points in the returned window), with maximum `1.846e-3`. Historical adaptive L minus batched A is positive at 100% of nonzero comparison points, with mean `3.83e-5`.

Similarly, L-final-minus-F and A-minus-F are positive at all points in the controlled fixed-grid comparison: both selection approximations omit absorption, but L omits materially more. No opposite-sign pruning effect was found in these LTE windows. The wavelength residuals are plotted in `rkints_second_pruning_residuals.png`.

## 11. Remove-one-line / cumulative removed-line test

The largest main-window L–A residual occurs at 5201.08139 Å. The responsible Ti I 5201.0814 Å line has:

```text
ALMAX                         1.7265e-2
legacy blended depression    1.6296e-3  (< accwi)
remove-one max delta F        1.8404e-3
```

The entire local L–A peak is `1.8459e-3`; that one line explains essentially all of it. The next-largest tested marginal is MgH 5201.0022 Å at `2.83e-4`. Removing the top five tested lines from A gives `1.8456e-3`; removing all 179 exactly reproduces the final legacy mask spectrum. This closes the causal chain: the metal-poor discrepancy is the opacity of lines removed by the second pruning, not an unrelated cached-state artifact.

## 12. Metal-poor discrepancy decomposition

Two effects must be kept separate:

1. **Second pruning:** on one common fine grid, L-final versus A reaches `1.846e-3` in the returned window and is explained by the removed lines, dominated locally by Ti I 5201.0814 Å.
2. **Historical adaptive-grid/range state:** historical L versus F and batched A versus F both contain additional grid/interpolation differences (about `4.95e-3` maximum). These can partly cancel after convolution and should not be used to judge the pruning rule.

The common-grid mask experiment is therefore the primary causal comparison. Under that comparison, removing the historical approximation moves the spectrum decisively toward F.

## 13. Performance cost of retaining all ALMAX-selected lines

The following measurements use the exact same batched grid, continuum cache, interval index, and transfer machinery. Only the final active mask changes. Times are medians of three runs.

| Model | L-mask total | A total | A overhead | L/A candidates | A candidates | OPMTRX L / A |
|---|---:|---:|---:|---:|---:|---:|
| Solar | 0.1307 s | 0.1317 s | +0.8% noise | 137,758 | 137,758 | 0.0596 / 0.0594 s |
| Metal-poor | 0.0453 s | 0.0514 s | **+13.4%** | 16,120 | 27,122 | 0.00954 / 0.01487 s |
| Cool metal-rich | 4.9218 s | 5.0216 s | +2.0% noise | 8,588,328 | 8,588,328 | 3.376 / 3.373 s |
| Giant | 4.3504 s | 4.3874 s | +0.9% noise | 7,640,032 | 7,640,032 | 2.937 / 2.928 s |

Only the metal-poor case has a real mask difference. There, retaining the 179 lines adds 11,002 candidate visits (+68.2%) and 5.3 ms of OPMTRX time (+55.9%), but just 6.1 ms or 13.4% total transfer time. In the other models the masks and candidate counts are identical, so the small timing differences are run noise.

## 14. NLTE implications

`OPMTRX` tests `MARK` before entering its NLTE branch. A pruned NLTE line therefore loses both its NLTE extinction and its line emissivity/source contribution. A disk-centre blended-depth heuristic is especially difficult to justify as a gate for both terms.

The existing solar H-alpha departure-coefficient test (6561–6564.2 Å) contained 983 supported lines, 136 ALMAX-selected lines, and one NLTE H line. It tested 88 centres and pruned none, including no NLTE line, so L and A were identical in this sanity case. A differs from F by `1.52e-5` intrinsic (`1.03e-5` at R=20k, `1.33e-5` at R=60k), while P differs by only `2.39e-7`.

This case does not demonstrate an actual erroneous NLTE deletion, but the source semantics establish a real architectural risk. Removing second pruning makes NLTE behavior easier to reason about and test.

## 15. Turbospectrum comparison

The local Turbospectrum `bsyn.f` synthesis path uses a line-centric wavelength traversal. It stops extending a line profile after a local support criterion such as `kappa_line/kappa_cont < eps` is met. It does **not** perform the RKINTS sequence “admit line, evaluate blended emergent centre depth, permanently deactivate line, alter all redder wavelengths.”

This is a useful precedent for separating local profile support from global active-line membership, but is not by itself the reason for the PySME recommendation.

## 16. Recommended production semantics

The responsibilities should be explicit:

```text
ALMAX     -> whether a line participates in synthesis
accrt     -> how far in wavelength that line contributes
accwi     -> transfer-wavelength interpolation/refinement accuracy only
```

In the new batched RKINTS path, the active mask should be immutable during grid construction and transfer. This gives:

- better correctness against the minimally pruned reference;
- deterministic, order-independent line membership;
- one consistent opacity/source state for batched grid generations;
- cleaner NLTE semantics;
- compatibility with interval indexing, threads, and future parallel evaluation;
- an API in which selection, support, and sampling tolerances do not have hidden cross-effects.

The historical sequential path may retain `MARK=2` behind an explicit compatibility/reference mode while migration tests are completed. The small non-pruning solar/giant legacy-versus-batched residuals should be tracked separately; they are not a reason to preserve second pruning.

## 17. Go / no-go conclusion

This audit meets the strong-case criteria to remove second pruning:

1. **How many lines?** Zero in solar/cool/giant main cases; 179/209/171 in the three metal-poor windows, corresponding to 44–55% of overlapping ALMAX lines.
2. **Negligible?** No. Removed ALMAX reaches 0.087 in the main window and 0.205 in the clean window; one removed Ti I line changes flux by `1.84e-3`.
3. **Order-dependent?** Yes. Reversing only decision order changes two lines and the spectrum by up to `3.22e-3`.
4. **Origin of the metal-poor `1.46e-3`?** It is the opacity of second-pruned lines; Ti I 5201.0814 Å dominates the largest feature.
5. **Which is closer to F?** ALMAX-only, by 3.7–13.4× on the main intrinsic metrics and by 6.4–10.1× at R=20k/60k.
6. **Performance cost?** 13.4% (6.1 ms) in the affected metal-poor main case; effectively zero when no extra lines are pruned.
7. **Observable impact?** Main metal-poor L versus F reaches `5.20e-4` at R=20k and `1.10e-3` at R=60k, versus `8.14e-5` and `1.09e-4` for A.
8. **NLTE risk?** Yes: `MARK=2` suppresses both NLTE extinction and emissivity/source terms, although the H-alpha sanity case pruned none.
9. **Should `accwi` only control sampling?** Yes.
10. **Remove from batched RKINTS?** **Yes.** Treat it as a historical performance approximation, not required physics or accuracy logic.

## Reproducibility artifacts

- `rkints_second_pruning_audit.py`: LTE audit and benchmarks
- `rkints_second_pruning_nlte_audit.py`: NLTE sanity case
- `rkints_second_pruning_reference_convergence.py`: F-reference convergence
- `rkints_second_pruning_plot.py`: residual plot generation
- `rkints_second_pruning_audit.json`: metrics and per-line records
- `rkints_second_pruning_audit.npz`: spectra, masks, and wavelengths
- `rkints_second_pruning_logs/`: raw diagnostic decision logs
- `rkints_second_pruning_residuals.png`: principal residuals
