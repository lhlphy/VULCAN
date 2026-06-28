# VULCAN Kinetics Format

## Contents

- Ordinary thermal reactions
- Pressure-dependent reactions
- Unit conventions
- Formula conversions
- Reverse reactions and photolysis
- Verification checklist

## Ordinary Thermal Reactions

VULCAN reads:

```text
id [ reactants -> products ] A B C Reference Tmin-Tmax
```

and evaluates:

```text
k(T) = A * T**B * exp(-C/T)
```

Example:

```text
115 [ OH + CH4 -> H2O + CH3 ] 1.68E-18 2.180 1231.0 Reference 195-2025
```

The parser requires the first three tokens after `]` to be numeric. Reference and temperature text are metadata for humans; VULCAN does not enforce the stated validity range.

## Pressure-Dependent Reactions

Under `# 3-body and Disscoiation Reactions`, VULCAN reads:

```text
id [ reaction ] A0 B0 C0 Ainf Binf Cinf Reference Temperature-ranges
```

with:

```text
k0(T)   = A0   * T**B0   * exp(-C0/T)
kinf(T) = Ainf * T**Binf * exp(-Cinf/T)
```

The current code applies its own simple falloff expression. It does not generally parse arbitrary Troe, PLOG, Chebyshev, or collider-efficiency fields. A special hard-coded reaction is not a template for adding new reactions.

Expected kinetic units depend on the limiting molecularity:

- `A + B + M -> products + M`: `k0` in `cm6 molecule-2 s-1`; `kinf` in `cm3 molecule-1 s-1`.
- `AB + M -> products + M`: `k0` in `cm3 molecule-1 s-1`; `kinf` in `s-1`.

Under `# 3-body reactions without high-pressure rates`, only `A0 B0 C0` are read.

## Unit Conventions

VULCAN evolves number density in `molecule cm-3`. For a rate coefficient expressed per mole, convert by kinetic order `m`:

```text
A_molecule = A_mol / N_A**(m - 1)
```

Examples:

- First order: no Avogadro conversion.
- Bimolecular: divide by `N_A`.
- Termolecular: divide by `N_A**2`.

Use the kinetic order of the specific limit. A low-pressure association and its high-pressure limit have different orders.

Do not confuse:

- `cal/mol`, `kcal/mol`, `J/mol`, and `kJ/mol`
- `cm3/mol/s` and `m3/mol/s`
- `bar`, `atm`, `Torr`, and number density
- activation energy `Ea` and activation temperature `C = Ea/R`

## Formula Conversions

### NIST modified Arrhenius

Source:

```text
k = Aref * (T/Tref)**n * exp(-Ea/(R*T))
```

VULCAN:

```text
A = Aref / Tref**n
B = n
C = Ea / R
```

Convert `Aref` to molecule-based units before or after the temperature normalization.

Example NIST record for `OH + CH4 -> H2O + CH3`:

```text
Aref = 3.38E-13 cm3 molecule-1 s-1
Tref = 298 K
n    = 2.23
Ea   = 9840 J mol-1
```

Conversion:

```text
A = 1.02663E-18
B = 2.23
C = 1183.48 K
```

### KIDA/UMIST Kooij

Source:

```text
k = alpha * (T/Tref)**beta * exp(-gamma/T)
```

VULCAN:

```text
A = alpha / Tref**beta
B = beta
C = gamma
```

`Tref` is commonly 300 K, but read the source definition.

### CHEMKIN/RMG Arrhenius

Source:

```text
k = A * T**b * exp(-Ea/(R*T))
```

VULCAN:

```text
A_vulcan = A converted to cm/molecule units
B_vulcan = b
C_vulcan = Ea/R
```

Read the mechanism header for length, quantity, time, and activation-energy units. Never infer units from the numeric magnitude alone.

## Reverse Reactions and Photolysis

- VULCAN generates reverse rates for ordinary forward reactions before `# reverse stops` using thermodynamic data.
- Confirm all species exist in `thermo/all_compose.txt` and `thermo/NASA9/`.
- Do not add both directions unless the network design explicitly requires it; duplicate checks may fail and detailed balance may be broken.
- Reactions after `# reverse stops` are not standard automatically reversible thermal entries.
- Photolysis entries have no `A B C`. They map a species and branch index to cross sections and quantum yields.

## Verification Checklist

1. Confirm reaction direction and channel.
2. Confirm source formula and units.
3. Confirm `Tmin/Tmax` cover the model range of interest.
4. Confirm pressure regime, bath gas, and falloff form.
5. Convert units and formula independently.
6. Compare source and VULCAN `k(T)` at `Tmin`, a midpoint, and `Tmax`.
7. Check values against at least one alternative source.
8. Record extrapolation warnings explicitly.
9. Run VULCAN's network regeneration without `-n` after an approved network edit.

