#!/usr/bin/env python3
"""
Generate a hybrid-pressure coefficient sidecar for VULCAN GCM preprocessing.

The coefficient construction follows the user-supplied ``uniform_spacing()``
layout: one pure-pressure top layer followed by uniformly spaced sigma layers.
The output stores ``pk_Pa`` and ``bk`` on interfaces plus the corresponding
reference half-level pressure profile ``reference_phalf_Pa``.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import xarray as xr


SUPPORTED_PRESSURE_UNITS = {
    "pa": 1.0,
    "pascal": 1.0,
    "pascals": 1.0,
    "mb": 100.0,
    "mbar": 100.0,
    "millibar": 100.0,
    "millibars": 100.0,
    "bar": 1.0e5,
    "bars": 1.0e5,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a uniform-spacing hybrid coefficient NetCDF sidecar.",
    )
    parser.add_argument(
        "template_nc",
        type=Path,
        help="Template NetCDF used to infer the vertical size and default output name.",
    )
    parser.add_argument(
        "--output-path",
        type=Path,
        default=None,
        help="Output NetCDF path. Defaults to <template_stem>_hybrid_coeffs.nc next to the template file.",
    )
    parser.add_argument(
        "--ps0-pa",
        type=float,
        required=True,
        help="Reference surface pressure ps0 in Pa used by the hybrid constructor.",
    )
    parser.add_argument(
        "--ptop-pa",
        type=float,
        default=100.0,
        help="Top interface pressure in Pa. Default: 100 Pa (1 mb).",
    )
    parser.add_argument(
        "--pint-factor",
        type=float,
        default=0.4,
        help="Factor used in pint = ptop + pint_factor * ps0 / km. Default: 0.4.",
    )
    parser.add_argument(
        "--km",
        type=int,
        default=None,
        help="Number of full levels. Defaults to the template pfull length.",
    )
    parser.add_argument(
        "--reference-pressure-nc",
        type=Path,
        default=None,
        help="Optional NetCDF whose phalf axis is used for hard validation.",
    )
    parser.add_argument(
        "--reference-pressure-var",
        default="phalf",
        help="Pressure variable used in --reference-pressure-nc. Default: phalf.",
    )
    return parser.parse_args()


def load_dataset(path: Path) -> xr.Dataset:
    return xr.open_dataset(path, decode_times=False)


def pressure_to_pa(da: xr.DataArray, label: str) -> np.ndarray:
    units = str(da.attrs.get("units", "")).strip().lower()
    if units not in SUPPORTED_PRESSURE_UNITS:
        raise ValueError(
            f"Unsupported pressure units for {label}: {units or '<missing>'}. "
            f"Supported units: {', '.join(sorted(SUPPORTED_PRESSURE_UNITS))}."
        )
    return np.asarray(da.values, dtype=float) * SUPPORTED_PRESSURE_UNITS[units]


def infer_km(template_nc: Path, km_override: int | None) -> int:
    if km_override is not None:
        if km_override < 1:
            raise ValueError("--km must be positive.")
        return int(km_override)

    with load_dataset(template_nc) as ds:
        if "pfull" in ds:
            return int(ds["pfull"].size)
        if "phalf" in ds:
            return int(ds["phalf"].size - 1)
    raise ValueError("Unable to infer km from the template file; provide --km explicitly.")


def build_uniform_coeffs(km: int, ps0_pa: float, ptop_pa: float, pint_factor: float) -> tuple[np.ndarray, np.ndarray]:
    if km < 2:
        raise ValueError("km must be at least 2 for the hybrid grid construction.")
    if ps0_pa <= 0.0 or ptop_pa <= 0.0:
        raise ValueError("ps0_pa and ptop_pa must stay positive.")

    pint_pa = ptop_pa + pint_factor * ps0_pa / km
    pk_pa = np.zeros(km + 1, dtype=float)
    bk = np.zeros(km + 1, dtype=float)

    pk_pa[0] = ptop_pa
    pk_pa[1] = pint_pa

    for idx in range(2, km + 1):
        bk[idx] = (idx - 1) / (km - 1)
        pk_pa[idx] = pint_pa - bk[idx] * pint_pa

    pk_pa[-1] = 0.0
    bk[-1] = 1.0
    return pk_pa, bk


def validate_reference_profile(
    reference_nc: Path,
    reference_var: str,
    expected_phalf_pa: np.ndarray,
) -> None:
    with load_dataset(reference_nc) as ds:
        if reference_var not in ds:
            raise ValueError(f"Validation file {reference_nc} does not contain {reference_var}.")
        reference_phalf_pa = pressure_to_pa(ds[reference_var], f"{reference_nc}:{reference_var}")

    if reference_phalf_pa.shape != expected_phalf_pa.shape:
        raise ValueError(
            f"Reference phalf shape {reference_phalf_pa.shape} does not match "
            f"expected {expected_phalf_pa.shape}."
        )
    if not np.allclose(reference_phalf_pa, expected_phalf_pa, rtol=1.0e-10, atol=1.0e-8):
        max_diff = float(np.max(np.abs(reference_phalf_pa - expected_phalf_pa)))
        raise ValueError(
            "The generated reference_phalf_Pa does not match the provided reference pressure axis. "
            f"Max abs diff = {max_diff:.6e} Pa."
        )


def build_dataset(
    template_nc: Path,
    output_path: Path,
    km: int,
    ps0_pa: float,
    ptop_pa: float,
    pint_factor: float,
    pk_pa: np.ndarray,
    bk: np.ndarray,
) -> xr.Dataset:
    reference_phalf_pa = pk_pa + bk * ps0_pa
    return xr.Dataset(
        data_vars={
            "pk_Pa": xr.DataArray(pk_pa, dims=("interface",), attrs={"units": "Pa"}),
            "bk": xr.DataArray(bk, dims=("interface",), attrs={"units": "1"}),
            "reference_phalf_Pa": xr.DataArray(
                reference_phalf_pa,
                dims=("interface",),
                attrs={"units": "Pa", "description": "Reference half-level pressure used for hybrid-grid validation."},
            ),
        },
        coords={"interface": xr.DataArray(np.arange(km + 1, dtype=int), dims=("interface",))},
        attrs={
            "template_file": str(template_nc),
            "output_file": str(output_path),
            "ps0_Pa": float(ps0_pa),
            "ptop_Pa": float(ptop_pa),
            "pint_factor": float(pint_factor),
            "km": int(km),
            "ks": 1,
            "generator_formula": "pk[0]=ptop; pk[1]=ptop+pint_factor*ps0/km; bk[2:]=(i-1)/(km-1); pk[i]=pint-bk[i]*pint; reference_phalf=pk+bk*ps0",
        },
    )


def main() -> None:
    args = parse_args()
    output_path = args.output_path or args.template_nc.with_name(f"{args.template_nc.stem}_hybrid_coeffs.nc")
    km = infer_km(args.template_nc, args.km)
    pk_pa, bk = build_uniform_coeffs(km, args.ps0_pa, args.ptop_pa, args.pint_factor)
    reference_phalf_pa = pk_pa + bk * args.ps0_pa

    if args.reference_pressure_nc is not None:
        validate_reference_profile(args.reference_pressure_nc, args.reference_pressure_var, reference_phalf_pa)

    dataset = build_dataset(
        template_nc=args.template_nc,
        output_path=output_path,
        km=km,
        ps0_pa=args.ps0_pa,
        ptop_pa=args.ptop_pa,
        pint_factor=args.pint_factor,
        pk_pa=pk_pa,
        bk=bk,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    dataset.to_netcdf(output_path)
    print(f"Wrote hybrid coefficient sidecar to {output_path}")


if __name__ == "__main__":
    main()
