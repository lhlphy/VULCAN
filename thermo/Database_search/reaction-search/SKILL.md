---
name: reaction-search
description: Manually search, compare, select, and convert gas-phase chemical reaction kinetics for VULCAN. Use when Codex must find forward reaction parameters from NIST, KIDA, IUPAC, UMIST, RMG, JPL, GRI-Mech, or primary literature; resolve multiple rate records; check temperature and pressure applicability; convert rate expressions and units into VULCAN network format; or document the evidence for adding or updating a VULCAN reaction.
---

# Reaction Search

Find reaction kinetics manually, preserve the source evidence, and produce one defensible VULCAN-ready recommendation. Do not treat the database search as a fully automatic authority.

## Paths

- Treat the parent directory of this skill as the working root: `thermo/Database_search`.
- Store downloaded HTML, text, JSON, CSV, YAML, and PDFs in `../data_temp/`.
- Do not modify a VULCAN network unless the user explicitly requests it.
- Keep final source URLs, citations, retrieval date, and raw rate expression in the result.

## Required Workflow

1. Normalize the requested reaction.
   - Record reactants, products, direction, phase, electronic states, isomers, charge, target temperature range, target pressure range, and bath gas.
   - Treat `O`, `O_1`, `CH2`, `CH2_1`, cyclic/linear isomers, and ions as distinct species.
2. Check the active VULCAN network for an existing reaction and reference.
3. Search NIST first for ordinary gas-phase thermal kinetics.
   - Use the NIST web form manually or run `scripts/nist_search_example.py` to save the raw search page.
4. Search secondary sources according to the reaction domain.
   - Use KIDA or UMIST for astrochemical and planetary reactions.
   - Use IUPAC or JPL for atmospheric and photochemical evaluations.
   - Use RMG and combustion mechanisms for high-temperature combustion chemistry and estimates.
   - Search primary papers when evaluated databases are missing, inconsistent, or outside the target range.
5. Save raw query results under `../data_temp/`; never save temporary downloads inside the skill folder.
6. Apply hard applicability filters before ranking records. Read `references/selection-policy.md`.
7. Select one final record.
   - Prefer experimental or evaluated kinetics that cover the requested temperature and pressure.
   - Among otherwise comparable experimental records, prefer wider applicable temperature coverage, newer evaluations or measurements, explicit uncertainty, and documented pressure behavior.
   - Never select a wider range that does not match the reaction channel, bath gas, or pressure regime.
8. Convert the original expression and units to VULCAN form. Read `references/vulcan-format.md` and use `scripts/convert_to_vulcan.py` for arithmetic.
9. Verify the conversion by evaluating the source and converted formulas at the range endpoints and at one interior temperature.
10. Report the chosen record, rejected alternatives, selection reason, valid range, pressure behavior, uncertainty, original expression, converted VULCAN line, and warnings.

## Source Guidance

Read `references/databases.md` before querying an unfamiliar source. Database websites do not share a stable public API or output schema. Use a source-specific workflow and preserve the raw response.

Default order for thermal gas-phase reactions:

1. NIST evaluated/review record
2. NIST direct experiments covering the target conditions
3. KIDA/IUPAC/JPL recommendation appropriate to the domain
4. Recent primary experiment
5. High-level theory with tunneling/master-equation treatment when relevant
6. RMG, UMIST, analog, or structure-activity estimate

Do not interpret this order mechanically. A 200-300 K atmospheric evaluation is not suitable for a 2000 K VULCAN calculation.

## VULCAN Guardrails

- Use `k = A * T**B * exp(-C/T)` for ordinary thermal reactions.
- Match number-density units used by VULCAN; convert molar CHEMKIN/RMG coefficients by reaction order.
- Preserve pressure-dependent low- and high-pressure limits. Do not collapse Troe, PLOG, Chebyshev, or master-equation data to one `A B C` fit without an explicit justified approximation.
- Confirm that every species has composition and NASA9 data before relying on VULCAN's automatic reverse rate.
- Put an irreversible reaction after `# reverse stops` only with explicit physical justification.
- Do not assign Arrhenius parameters to photolysis. Photolysis requires cross sections, quantum yields, branch mapping, and stellar flux.
- Keep negative activation temperatures when the source genuinely reports them; do not silently take an absolute value.
- Treat source temperature limits as hard validity metadata, not comments to discard.

## Bundled Scripts

- `scripts/nist_search_example.py`: submit one NIST reaction-form search, save raw HTML, and print matched reaction groups.
- `scripts/kida_search_example.py`: submit one KIDA species search, fetch reaction-class HTML fragments, and save all raw pages.
- `scripts/convert_to_vulcan.py`: convert NIST, Kooij, or CHEMKIN Arrhenius parameters to VULCAN `A B C` and evaluate sample temperatures.

These are source-specific examples, not a unified database API. Inspect the saved raw pages whenever a parser returns no records or the website changes.

## Final Result Format

Return a compact evidence block:

```text
Reaction:
Target conditions:
Chosen source and reference:
Method: evaluation / experiment / theory / estimate
Original rate expression and units:
Temperature range:
Pressure range or pressure dependence:
VULCAN A, B, C:
VULCAN network line:
Verification temperatures and k(T):
Alternatives rejected:
Warnings:
Raw data path:
```

