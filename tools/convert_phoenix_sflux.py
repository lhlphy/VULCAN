#!/usr/bin/env python3
"""Convert a PHOENIX surface spectrum to VULCAN's stellar-flux format.

Expected input columns:
    Wavelength(m)    Flux(W/m^3)

VULCAN expects:
    WL(nm)           Flux(erg/cm^2/s/nm)

For spectral flux density F_lambda,

    wavelength_nm = wavelength_m * 1e9
    flux_cgs_per_nm = flux_si_per_m * 1e7 / 1e4 * 1e-9
                    = flux_si_per_m * 1e-6

The input may contain ``*BEGIN_DATA`` and ``*END`` marker lines.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Iterable


M_TO_NM = 1.0e9
W_M3_TO_ERG_CM2_S_NM = 1.0e-6


def read_phoenix_spectrum(path: Path) -> list[tuple[float, float]]:
    """Read the first two numeric columns from a PHOENIX text spectrum."""
    rows: list[tuple[float, float]] = []
    in_data = False
    saw_begin_marker = False

    with path.open("r", encoding="utf-8-sig") as stream:
        for line_number, raw_line in enumerate(stream, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            if line.upper() == "*BEGIN_DATA":
                saw_begin_marker = True
                in_data = True
                continue
            if line.upper().startswith("*END"):
                break
            if saw_begin_marker and not in_data:
                continue

            fields = line.replace(",", " ").split()
            if len(fields) < 2:
                continue
            try:
                wavelength_m = float(fields[0])
                flux_w_m3 = float(fields[1])
            except ValueError:
                # Permit a plain-text column header before the data.
                if not rows:
                    continue
                raise ValueError(
                    f"{path}:{line_number}: expected two numeric columns"
                ) from None

            if not math.isfinite(wavelength_m) or wavelength_m <= 0.0:
                raise ValueError(
                    f"{path}:{line_number}: wavelength must be finite and positive"
                )
            if not math.isfinite(flux_w_m3) or flux_w_m3 < 0.0:
                raise ValueError(
                    f"{path}:{line_number}: flux must be finite and non-negative"
                )
            if rows and wavelength_m <= rows[-1][0]:
                raise ValueError(
                    f"{path}:{line_number}: wavelengths must be strictly increasing"
                )
            rows.append((wavelength_m, flux_w_m3))

    if len(rows) < 2:
        raise ValueError(f"{path}: found fewer than two spectral data rows")
    return rows


def convert_rows(
    rows: Iterable[tuple[float, float]],
) -> list[tuple[float, float]]:
    """Convert m and W/m^3 to nm and erg/cm^2/s/nm."""
    return [
        (wavelength_m * M_TO_NM, flux_w_m3 * W_M3_TO_ERG_CM2_S_NM)
        for wavelength_m, flux_w_m3 in rows
    ]


def trapezoid_integral(rows: list[tuple[float, float]]) -> float:
    """Integrate the second column over the first with the trapezoid rule."""
    return math.fsum(
        0.5 * (left[1] + right[1]) * (right[0] - left[0])
        for left, right in zip(rows, rows[1:])
    )


def write_vulcan_spectrum(
    path: Path, rows: list[tuple[float, float]], overwrite: bool
) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(f"{path} already exists; pass --overwrite to replace it")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="\n") as stream:
        stream.write("# WL(nm)    Flux(ergs/cm**2/s/nm)\n")
        for wavelength_nm, flux_cgs_nm in rows:
            stream.write(f"{wavelength_nm:.10e}    {flux_cgs_nm:.10e}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert Wavelength(m), Flux(W/m^3) PHOENIX spectra to the "
            "VULCAN sflux format."
        )
    )
    parser.add_argument("input", type=Path, help="input PHOENIX text spectrum")
    parser.add_argument(
        "output",
        type=Path,
        nargs="?",
        help="output path (default: sflux-<input-name>.txt beside the input)",
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="replace an existing output file"
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_path = args.input.resolve()
    output_path = args.output
    if output_path is None:
        output_path = input_path.with_name(f"sflux-{input_path.name}.txt")
    else:
        output_path = output_path.resolve()

    input_rows = read_phoenix_spectrum(input_path)
    output_rows = convert_rows(input_rows)
    write_vulcan_spectrum(output_path, output_rows, args.overwrite)

    input_bolometric_w_m2 = trapezoid_integral(input_rows)
    output_bolometric_cgs = trapezoid_integral(output_rows)
    expected_output_cgs = input_bolometric_w_m2 * 1.0e3
    relative_error = abs(output_bolometric_cgs / expected_output_cgs - 1.0)

    print(f"Converted {len(output_rows)} rows")
    print(f"Wavelength range: {output_rows[0][0]:.6g}-{output_rows[-1][0]:.6g} nm")
    print(f"Input integral:  {input_bolometric_w_m2:.10e} W/m^2")
    print(f"Output integral: {output_bolometric_cgs:.10e} erg/cm^2/s")
    print(f"Integral conversion relative error: {relative_error:.3e}")
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
