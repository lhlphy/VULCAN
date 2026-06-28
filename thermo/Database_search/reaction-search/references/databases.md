# Database Reference

## Contents

- NIST Chemical Kinetics Database
- KIDA
- IUPAC Atmospheric Kinetics
- UMIST Database for Astrochemistry
- RMG Database
- NASA/JPL Evaluation
- GRI-Mech and mechanism files
- Access and parsing rules

No database has one universal temperature or pressure range. Treat `Tmin`, `Tmax`, pressure, collider, phase, and method as properties of each reaction record.

## Verified Access Snapshot

Observed from the local `exoctm` Python/OpenSSL environment on 2026-06-20. Recheck at query time.

| Source | Access result | Preferred route |
|---|---|---|
| NIST | Search and detail HTML returned successfully | Browser or form POST with a `User-Agent` and cookies |
| KIDA | Search, AJAX fragments, and detail pages returned successfully | Website search or the bundled example script |
| IUPAC | Home and catalogue returned successfully | Catalogue and evaluated datasheets |
| UMIST | New Rate22 search and download pages returned successfully | `umistdatabase.uk`, not the obsolete hostname |
| RMG | Database, kinetics search, and library pages returned successfully | Website or versioned RMG-database download |
| JPL | Direct site had an SSL EOF in this environment | Retry; otherwise use NASA NTRS or a traceable evaluated mirror |
| GRI-Mech | Historical website links/certificate were unreliable | Use a maintained Cantera/GitHub mechanism mirror |

## Reference Reaction Range Example

For `OH + CH4 -> H2O + CH3`, records encountered during validation illustrate why ranges must be stored per entry:

| Source/record | Expression | Temperature | Pressure statement | Method note |
|---|---|---:|---|---|
| NIST 2021 review | `3.38E-13*(T/298)^2.23*exp(-9840/(R*T))` | 200-2025 K | pressure independent | extensive literature review |
| KIDA/IUPAC evaluation | `1.85E-12*exp(-1690/T)` | 200-300 K | no finite range listed for this bimolecular entry | evaluated recommendation |
| JPL 2011 as recorded by KIDA | `2.45E-12*exp(-1780/T)` | 200-500 K | no finite range listed for this bimolecular entry | review/evaluation |
| UMIST as recorded by KIDA | `3.77E-13*(T/300)^2.42*exp(-1160/T)` | 10-800 K | no finite range listed | estimate; lower confidence |
| GRI-Mech 3.0 | CHEMKIN `A=1.0E8, b=1.6, Ea=3120 cal/mol` | no strict per-record validity range in the mechanism file | elementary expression has no explicit pressure dependence | mechanism value |

Do not transfer these ranges to another reaction. They are examples of source metadata and evidence quality.

## NIST Chemical Kinetics Database

- URL: `https://kinetics.nist.gov/kinetics/`
- Best use: default source for gas-phase elementary kinetics and literature history.
- Access: HTML form POST plus session cookies; detail pages are HTML.
- Data: many experimental, theoretical, review, and evaluated records for one reaction.
- Typical expression:

  `k = Aref * (T/Tref)^n * exp(-Ea/(R*T))`

- Range metadata: temperature per record; pressure units/preferences and pressure-dependence notes may be present.
- Important: scripted requests may require a browser-like `User-Agent`. There is no stable public JSON API. A result page layout change can break HTML extraction.
- Selection: prefer an evaluated/review record covering the target conditions, then direct experiments. Do not average all records.

## KIDA

- URL: `https://kida.astrochem-tools.org/`
- Best use: astrochemistry, planetary atmospheres, neutral/ion reactions, low-temperature kinetics, and expert recommendations.
- Access: HTML search page plus internal AJAX POST endpoints returning HTML fragments.
- Data: reaction channels, `alpha`, `beta`, `gamma`, temperature limits, uncertainty factors, method, expert status, and references.
- Typical Kooij expression:

  `k = alpha * (T/300)^beta * exp(-gamma/T)`

- Range metadata: temperature per entry. Pressure is often absent for simple bimolecular entries; termolecular records require separate interpretation.
- Important: internal AJAX routes are website implementation details, not a promised API. Save both the initial and fragment HTML.
- Selection: distinguish recommended/expert-reviewed values from imported or estimated entries.

## IUPAC Atmospheric Kinetics

- URL: `https://iupac.aeris-data.fr/catalogue/`
- Best use: evaluated atmospheric gas-phase, heterogeneous, aqueous, and photolysis data.
- Access: web catalogue backed by a single-page application; datasheets and catalogue records may use different formats.
- Data: recommended expressions, uncertainty, temperature range, pressure/falloff parameters, and evaluation discussion.
- Range metadata: reaction-specific. Atmospheric recommendations are often narrower than combustion or exoplanet ranges.
- Important: use the evaluated datasheet as evidence. Do not assume an internal frontend JSON route is a stable public API.

## UMIST Database for Astrochemistry

- URL: `https://umistdatabase.uk/`
- Current release family: Rate22; use the current site, not the obsolete `udfa.ajmarkwick.net` hostname.
- Best use: large astrochemical networks and low-temperature reaction coverage.
- Access: HTML search and versioned downloadable database files.
- Data: reaction type, reactants/products, Kooij parameters, `Tmin`, `Tmax`, accuracy metadata, and reference codes.
- Typical expression:

  `k = alpha * (T/300)^beta * exp(-gamma/T)`

- Range metadata: temperature per row; pressure is usually represented by reaction type rather than a general pressure interval.
- Important: prefer the versioned download format over scraping HTML. Many values are estimates or network choices rather than evaluated laboratory recommendations.

## RMG Database

- URL: `https://rmg.mit.edu/database/kinetics/search/`
- Download: use the linked versioned RMG-database archive or repository.
- Best use: combustion/high-temperature chemistry, reaction families, rate rules, estimates, pressure-dependent kinetics, and mechanism generation.
- Access: web search and structured Python database files.
- Data models include Arrhenius, MultiArrhenius, ArrheniusEP, Troe, Lindemann, PDepArrhenius, and Chebyshev.
- Range metadata: record-specific; pressure-dependent models may include explicit temperature and pressure bounds.
- Important: use library/depository records before family estimates. Preserve the model type and uncertainty. Do not force Chebyshev/PLOG/Troe data into one Arrhenius fit.

## NASA/JPL Evaluation

- URL: `https://jpldataeval.jpl.nasa.gov/`
- Best use: evaluated atmospheric kinetics, photolysis cross sections, quantum yields, and falloff data.
- Access: evaluation reports and PDFs, not a reaction-by-reaction public API.
- Data: recommended expressions, uncertainty functions, temperature ranges, pressure dependence, and photochemical tables.
- Important: cite the report number and page/table. If the site is inaccessible, use NASA NTRS or a database record that explicitly cites the JPL recommendation, then retain that provenance.

## GRI-Mech and Mechanism Files

- Cantera mirror: `https://github.com/Cantera/cantera/blob/main/data/gri30.yaml`
- Best use: a fixed natural-gas combustion mechanism and comparison source, not a general kinetics database.
- Access: CHEMKIN or Cantera YAML.
- Data: mechanism-specific Arrhenius and falloff records with declared unit systems.
- Range metadata: a mechanism may not state a strict validity range for every individual reaction. Species thermo ranges are not reaction-rate validity ranges.
- Important: the historical GRI-Mech website may have certificate or stale-link problems. Prefer a maintained, versioned mirror and cite the mechanism version.

## Access and Parsing Rules

1. Check HTTP access at query time; access notes become stale.
2. Use a descriptive `User-Agent` and low request rate.
3. Respect site terms, download licenses, and robots policies.
4. Cache raw responses under `data_temp` with source and retrieval date.
5. Never infer zero from a blank field; preserve it as missing.
6. Detect units from the record, file header, or user preferences rather than assuming them.
7. Fail loudly if expected headings or formula fields disappear.
8. Keep one parser per source and version. Normalize only after source parsing succeeds.
