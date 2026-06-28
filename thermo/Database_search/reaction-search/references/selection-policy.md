# Selecting One Kinetic Record

## Contents

- Hard filters
- Evidence ranking
- Tie-breaking
- Pressure-dependent records
- Required audit trail

## Hard Filters

Reject a record before scoring when any of these fail:

1. Exact reactants and products, including electronic state, charge, and isomer.
2. Correct reaction direction or a justified detailed-balance conversion.
3. Correct gas phase and reaction type.
4. Target temperature coverage, unless the user explicitly accepts extrapolation.
5. Target pressure/falloff regime and bath gas where relevant.
6. Compatible units and sufficient formula metadata.

Do not choose a wide temperature range for the wrong product channel.

## Evidence Ranking

Use this default ranking after hard filters:

1. Current evaluated recommendation based on multiple experiments.
2. Recent direct experiment with explicit uncertainty and broad target coverage.
3. Older direct experiment covering the target conditions.
4. High-level theoretical kinetics with an appropriate treatment of tunneling, variational effects, and pressure dependence.
5. Database estimate, rate rule, analog, or structure-activity relationship.
6. Untraceable mechanism value.

Label imported database values by their original method. A value appearing in KIDA or UMIST is not automatically an experiment.

## Tie-Breaking

Among records at the same evidence level:

1. Prefer complete coverage of the requested range.
2. Prefer the wider experimentally supported range.
3. Prefer newer evaluations or measurements when quality is comparable.
4. Prefer explicit uncertainty and pressure-dependence statements.
5. Prefer direct absolute measurements over relative-rate measurements.
6. Prefer data with accessible primary references and reproducible units.

Do not rank by publication year alone. A recent estimate does not outrank a slightly older evaluated experiment.

If no single record covers the full target range, do not silently splice fits. Report the gap and either:

- choose the record covering the scientifically important range,
- retain a documented piecewise expression if VULCAN is extended to support it,
- fit a new expression with residual/error documentation, or
- request a narrower model range.

## Pressure-Dependent Records

- Separate low-pressure, high-pressure, and falloff records.
- Match collider and collider efficiencies.
- Do not compare `k0` numerically with `kinf` without units and number density.
- Preserve Troe/PLOG/Chebyshev/master-equation data in the evidence record even if current VULCAN cannot represent it directly.
- Require an explicit approximation decision before reducing a pressure-dependent model.

## Required Audit Trail

For the selected record, keep:

- database and record identifier
- full citation and URL
- retrieval date
- original formula and units
- method category
- uncertainty
- `Tmin/Tmax` and `Pmin/Pmax` or pressure statement
- conversion equations
- converted VULCAN parameters
- numerical spot checks
- rejected alternatives and rejection reasons
- raw download path

