# Abundance

`sme.abund` stores the elemental abundance pattern, while `sme.monh` stores the
global metallicity offset.

In PySME, the abundance used in synthesis is:

```text
A(X) = pattern(X) + monh
```

for elements heavier than He. Hydrogen and helium are not shifted by `monh`.

## Recommended explicit abundance views

PySME provides explicit abundance-scale views so you can choose which quantity
you want to read or write:

```py
sme.abund.A["Ti"]          # absolute abundance used in synthesis, H=12 scale
sme.abund.xh["Ti"]         # [Ti/H]
sme.abund.xm["Ti"]         # [Ti/M]_PySME = [Ti/H] - sme.monh
sme.abund.pattern["Ti"]    # current internal SME pattern before adding monh
sme.abund.reference["Ti"]  # fixed abundance reference pattern
sme.abund.solar["Ti"]      # alias to reference for built-in solar patterns
```

Here `M` in `xm` means PySME's global metallicity parameter `sme.monh`, not a
recomputed total metal mass fraction `Z`.

## Legacy direct assignment

Legacy direct assignment is still supported:

```py
sme.abund["Ti"] = 5.105
```

but this is ambiguous because it writes the **internal SME pattern**, not the
abundance used in synthesis. New code should prefer:

```py
sme.abund.A["Ti"] = 3.909
```

or, if you intentionally want to set the internal pattern directly:

```py
sme.abund.pattern["Ti"] = 5.105
```

## Practical examples

```py
from pysme.abund import Abund

sme.abund = Abund(pattern="asplund2021", monh=-0.20)
```

### 1. Set a target absolute abundance used in synthesis

If you want final `A(Mg) = 7.40`:

```py
sme.abund.A["Mg"] = 7.40
```

### 2. Set a scaled-solar abundance offset relative to H

If you want `[Mg/H] = -0.15` relative to the reference pattern:

```py
sme.abund.xh["Mg"] = -0.15
```

### 3. Set an element enhancement at fixed metallicity

If you want `[Mg/M]_PySME = +0.30` at the current `monh`:

```py
sme.abund.xm["Mg"] = 0.30
```

For the common case where Fe follows `monh`, this is equivalent to the usual
idea of `[Mg/Fe]`. PySME does not add a dedicated `xfe` view in this first API.

### 4. Work directly with the internal SME pattern

If you want to inspect or set the underlying pattern before `monh` is applied:

```py
current_pattern = dict(sme.abund.pattern)
sme.abund.pattern["Mg"] = 7.60
```

### 5. Inspect the fixed abundance reference

The reference pattern remains fixed even after you modify the current pattern:

```py
ref_ti = sme.abund.reference["Ti"]
```

If the abundance object was initialized from a built-in solar pattern, the same
value is also available through:

```py
sme.abund.solar["Ti"]
```

For custom user-provided abundance patterns, use `reference` instead of
`solar`.
