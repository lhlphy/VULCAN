#!/usr/bin/env python3
"""
Extract one GCM column and convert it into VULCAN-friendly driver files.

This script reads a 3D GCM NetCDF output, selects the nearest column to a
requested latitude/longitude, converts pressure from mb to dyne/cm^2, corrects
the background pressure into the total pressure implied by SiO and SiO2 mole
numbers, and optionally interpolates each time-dependent profile onto a VULCAN
pressure grid.

Notes
-----
- The input file uses a non-standard calendar, so time is read with
  ``decode_times=False`` and kept as raw hours.
- ``sphum`` and ``cld_amt`` are interpreted as SiO and SiO2 mass loadings.
- A background mean molecular weight is required to recover the SiO/SiO2 mole
  numbers and therefore the total pressure.
- For pressures above the modeled GCM top (P < P_min), the interpolation uses
  the topmost GCM value instead of extrapolating.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
from typing import Iterable

import numpy as np
import xarray as xr


MBAR_TO_DYN_CM2 = 1.0e3
DEFAULT_SIO_MW = 44.0845
DEFAULT_SIO2_MW = 60.0835


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract a single GCM column and prepare VULCAN driver files."
    )
    parser.add_argument("input_nc", type=Path, help="Input GCM NetCDF file.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("VULCAN/atm/gcm_columns"),
        help="Directory for generated files.",
    )
    parser.add_argument(
        "--output-prefix",
        default=None,
        help="Prefix for generated files. Defaults to the NetCDF stem.",
    )
    parser.add_argument("--target-lat", type=float, default=0.0, help="Requested latitude in degrees.")
    parser.add_argument("--target-lon", type=float, default=0.0, help="Requested longitude in degrees.")
    parser.add_argument(
        "--lat-index",
        type=int,
        default=None,
        help="Use an explicit latitude index instead of nearest-neighbor lookup.",
    )
    parser.add_argument(
        "--lon-index",
        type=int,
        default=None,
        help="Use an explicit longitude index instead of nearest-neighbor lookup.",
    )
    parser.add_argument(
        "--vulcan-cfg",
        type=Path,
        default=None,
        help="Optional vulcan_cfg.py used to build the target VULCAN pressure grid.",
    )
    parser.add_argument(
        "--target-pressure-file",
        type=Path,
        default=None,
        help="Optional text file whose first column is the target pressure grid in dyne/cm^2.",
    )
    parser.add_argument(
        "--write-tp-snapshot-index",
        type=int,
        default=None,
        help="If set, also write one VULCAN atm_file text snapshot for the selected time index.",
    )
    parser.add_argument(
        "--background-mean-mol-weight",
        "--mean-mol-weight",
        dest="background_mean_mol_weight",
        type=float,
        required=True,
        help="Background-gas mean molecular weight used to convert SiO/SiO2 mass loadings into mole numbers.",
    )
    parser.add_argument(
        "--humidity-basis",
        choices=("total", "background"),
        default="total",
        help=(
            "Interpret sphum/cld_amt as mass fractions relative to total column mass "
            "('total') or relative to the background-gas mass only ('background')."
        ),
    )
    parser.add_argument(
        "--clip-negative",
        action="store_true",
        default=True,
        help="Clip negative sphum/cld_amt values to zero. Enabled by default.",
    )
    parser.add_argument(
        "--no-clip-negative",
        dest="clip_negative",
        action="store_false",
        help="Keep negative sphum/cld_amt values unchanged.",
    )
    return parser.parse_args()


def load_dataset(path: Path) -> xr.Dataset:
    return xr.open_dataset(path, decode_times=False)


def resolve_lon_delta(lon_values: np.ndarray, target_lon: float) -> np.ndarray:
    wrapped_target = np.mod(target_lon, 360.0)
    delta = np.abs(np.mod(lon_values - wrapped_target + 180.0, 360.0) - 180.0)
    return delta


def pick_column_indices(
    ds: xr.Dataset,
    target_lat: float,
    target_lon: float,
    lat_index: int | None,
    lon_index: int | None,
) -> tuple[int, int]:
    lat_values = ds["grid_yt"].values
    lon_values = ds["grid_xt"].values

    if lat_index is None:
        lat_index = int(np.argmin(np.abs(lat_values - target_lat)))
    if lon_index is None:
        lon_index = int(np.argmin(resolve_lon_delta(lon_values, target_lon)))

    return lat_index, lon_index


def clip_nonphysical(values: np.ndarray) -> np.ndarray:
    return np.maximum(values, 0.0)


def mass_loading_to_moles_per_background_mole(
    mass_loading: np.ndarray,
    background_mean_mol_weight: float,
    species_mol_weight: float,
    humidity_basis: str,
    total_condensible_mass_loading: np.ndarray,
) -> np.ndarray:
    if humidity_basis == "background":
        return mass_loading * background_mean_mol_weight / species_mol_weight

    background_mass_fraction = 1.0 - total_condensible_mass_loading
    if np.any(background_mass_fraction <= 0.0):
        raise ValueError(
            "sphum + cld_amt must stay below 1 when humidity-basis='total'."
        )
    return (mass_loading / background_mass_fraction) * background_mean_mol_weight / species_mol_weight


def derive_total_pressure_and_mole_fractions(
    background_pressure_dyn_cm2: np.ndarray,
    sio_mass_loading: np.ndarray,
    sio2_mass_loading: np.ndarray,
    background_mean_mol_weight: float,
    humidity_basis: str,
) -> dict[str, np.ndarray]:
    total_condensible_mass_loading = sio_mass_loading + sio2_mass_loading
    sio_per_bg_mole = mass_loading_to_moles_per_background_mole(
        sio_mass_loading,
        background_mean_mol_weight,
        DEFAULT_SIO_MW,
        humidity_basis,
        total_condensible_mass_loading,
    )
    sio2_per_bg_mole = mass_loading_to_moles_per_background_mole(
        sio2_mass_loading,
        background_mean_mol_weight,
        DEFAULT_SIO2_MW,
        humidity_basis,
        total_condensible_mass_loading,
    )

    mole_ratio_total_to_background = 1.0 + sio_per_bg_mole + sio2_per_bg_mole
    total_pressure_dyn_cm2 = background_pressure_dyn_cm2 * mole_ratio_total_to_background

    sio_mole_fraction = sio_per_bg_mole / mole_ratio_total_to_background
    sio2_mole_fraction = sio2_per_bg_mole / mole_ratio_total_to_background
    background_mole_fraction = 1.0 / mole_ratio_total_to_background

    if humidity_basis == "background":
        total_mass = 1.0 + total_condensible_mass_loading
    else:
        total_mass = np.ones_like(total_condensible_mass_loading)
    background_moles = np.ones_like(total_condensible_mass_loading) / background_mean_mol_weight
    total_moles = background_moles * mole_ratio_total_to_background
    mixture_mean_mol_weight = total_mass / total_moles

    return {
        "total_pressure_dyn_cm2": total_pressure_dyn_cm2,
        "sio_mole_fraction": sio_mole_fraction,
        "sio2_mole_fraction": sio2_mole_fraction,
        "background_mole_fraction": background_mole_fraction,
        "mixture_mean_mol_weight": mixture_mean_mol_weight,
    }


def log_pressure_interp(
    p_source_dyn_cm2: np.ndarray,
    values: np.ndarray,
    p_target_dyn_cm2: np.ndarray,
) -> np.ndarray:
    mask = np.isfinite(p_source_dyn_cm2) & np.isfinite(values) & (p_source_dyn_cm2 > 0.0)
    if mask.sum() == 0:
        return np.full_like(p_target_dyn_cm2, np.nan, dtype=float)

    p_valid = p_source_dyn_cm2[mask]
    v_valid = values[mask]
    order = np.argsort(p_valid)
    p_sorted = p_valid[order]
    v_sorted = v_valid[order]

    logp_source = np.log10(p_sorted)
    logp_target = np.log10(np.asarray(p_target_dyn_cm2, dtype=float))

    # left: P < P_min, i.e. above the GCM top, use the topmost value.
    # right: P > P_max, use the deepest available value.
    return np.interp(logp_target, logp_source, v_sorted, left=v_sorted[0], right=v_sorted[-1])


def interpolate_time_series(
    p_source_dyn_cm2: np.ndarray,
    values_2d: np.ndarray,
    p_target_dyn_cm2: np.ndarray,
) -> np.ndarray:
    out = np.empty((values_2d.shape[0], len(p_target_dyn_cm2)), dtype=float)
    for i in range(values_2d.shape[0]):
        if np.ndim(p_source_dyn_cm2) == 1:
            p_source_row = p_source_dyn_cm2
        else:
            p_source_row = p_source_dyn_cm2[i]
        out[i] = log_pressure_interp(p_source_row, values_2d[i], p_target_dyn_cm2)
    return out


def load_vulcan_cfg_pressure_grid(vulcan_cfg_path: Path) -> np.ndarray:
    spec = importlib.util.spec_from_file_location("vulcan_cfg_runtime", vulcan_cfg_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load VULCAN config from {vulcan_cfg_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return np.logspace(np.log10(module.P_b), np.log10(module.P_t), module.nz)


def load_pressure_grid_from_text(path: Path) -> np.ndarray:
    table = np.genfromtxt(path, names=True, dtype=None, encoding=None, skip_header=1)
    if getattr(table, "dtype", None) is not None and "Pressure" in table.dtype.names:
        return np.asarray(table["Pressure"], dtype=float)

    raw = np.genfromtxt(path, comments="#", dtype=float)
    raw = np.atleast_2d(raw)
    return np.asarray(raw[:, 0], dtype=float)


def build_native_dataset(
    column: xr.Dataset,
    background_pressure_dyn_cm2: np.ndarray,
    actual_lat: float,
    actual_lon: float,
    background_mean_mol_weight: float,
    humidity_basis: str,
) -> xr.Dataset:
    sio_mass_fraction = column["sphum"].values.astype(float)
    sio2_mass_fraction = column["cld_amt"].values.astype(float)
    thermo = derive_total_pressure_and_mole_fractions(
        background_pressure_dyn_cm2[np.newaxis, :],
        sio_mass_fraction,
        sio2_mass_fraction,
        background_mean_mol_weight,
        humidity_basis,
    )
    background_pressure_2d = np.broadcast_to(background_pressure_dyn_cm2, sio_mass_fraction.shape)
    data_vars: dict[str, tuple[tuple[str, str], np.ndarray, dict[str, str]]] = {
        "background_pressure_dyn_cm2": (
            ("time", "level"),
            background_pressure_2d,
            {"description": "Background-only pressure from GCM pfull", "units": "dyne/cm^2"},
        ),
        "total_pressure_dyn_cm2": (
            ("time", "level"),
            thermo["total_pressure_dyn_cm2"],
            {"description": "Total pressure including background gas, SiO, and SiO2 treated as gas", "units": "dyne/cm^2"},
        ),
        "temp_K": (("time", "level"), column["temp"].values.astype(float), {"units": "K"}),
        "sio_mass_fraction": (
            ("time", "level"),
            sio_mass_fraction,
            {"description": "Raw GCM sphum mapped to SiO vapor", "units": "1"},
        ),
        "sio2_mass_fraction": (
            ("time", "level"),
            sio2_mass_fraction,
            {"description": "Raw GCM cld_amt mapped to SiO2 condensate", "units": "1"},
        ),
        "sio_mole_fraction": (
            ("time", "level"),
            thermo["sio_mole_fraction"],
            {"description": "SiO mole fraction implied by the supplied mass loading", "units": "1"},
        ),
        "sio2_mole_fraction": (
            ("time", "level"),
            thermo["sio2_mole_fraction"],
            {"description": "SiO2 mole fraction implied by the supplied mass loading", "units": "1"},
        ),
        "background_mole_fraction": (
            ("time", "level"),
            thermo["background_mole_fraction"],
            {"description": "Background-gas mole fraction after adding SiO and SiO2", "units": "1"},
        ),
        "mixture_mean_mol_weight": (
            ("time", "level"),
            thermo["mixture_mean_mol_weight"],
            {"description": "Mean molecular weight of the full mixture", "units": "g/mol"},
        ),
    }
    return xr.Dataset(
        data_vars={name: xr.DataArray(values, dims=dims, attrs=attrs) for name, (dims, values, attrs) in data_vars.items()},
        coords={
            "time": xr.DataArray(column["time"].values.astype(float), dims=("time",), attrs={"units": "hours"}),
            "level": xr.DataArray(np.arange(background_pressure_dyn_cm2.size), dims=("level",)),
            "background_pressure_level_dyn_cm2": xr.DataArray(
                background_pressure_dyn_cm2.astype(float),
                dims=("level",),
                attrs={"units": "dyne/cm^2", "description": "Static background pressure grid from GCM pfull"},
            ),
        },
        attrs={
            "source_file": str(column.attrs.get("source_file", "")),
            "selected_lat_deg": float(actual_lat),
            "selected_lon_deg": float(actual_lon),
            "pressure_source_units": "mb",
            "pressure_output_units": "dyne/cm^2",
            "pressure_definition": "total pressure including background gas, SiO vapor, and SiO2 treated as gas",
            "humidity_basis": humidity_basis,
            "background_mean_mol_weight_g_per_mol": float(background_mean_mol_weight),
        },
    )


def build_interpolated_dataset(
    native_ds: xr.Dataset,
    p_target_dyn_cm2: np.ndarray,
) -> xr.Dataset:
    p_source = native_ds["total_pressure_dyn_cm2"].values
    temp_interp = interpolate_time_series(p_source, native_ds["temp_K"].values, p_target_dyn_cm2)
    sio_interp = interpolate_time_series(p_source, native_ds["sio_mass_fraction"].values, p_target_dyn_cm2)
    sio2_interp = interpolate_time_series(p_source, native_ds["sio2_mass_fraction"].values, p_target_dyn_cm2)
    sio_mole_interp = interpolate_time_series(p_source, native_ds["sio_mole_fraction"].values, p_target_dyn_cm2)
    sio2_mole_interp = interpolate_time_series(p_source, native_ds["sio2_mole_fraction"].values, p_target_dyn_cm2)
    background_mole_interp = interpolate_time_series(
        p_source, native_ds["background_mole_fraction"].values, p_target_dyn_cm2
    )
    mixture_mmw_interp = interpolate_time_series(
        p_source, native_ds["mixture_mean_mol_weight"].values, p_target_dyn_cm2
    )
    background_pressure_interp = interpolate_time_series(
        p_source, native_ds["background_pressure_dyn_cm2"].values, p_target_dyn_cm2
    )

    data_vars: dict[str, tuple[tuple[str, str], np.ndarray, dict[str, str]]] = {
        "temp_K": (("time", "pressure"), temp_interp, {"units": "K"}),
        "background_pressure_dyn_cm2": (
            ("time", "pressure"),
            background_pressure_interp,
            {"description": "Background-only pressure interpolated onto the target total-pressure grid", "units": "dyne/cm^2"},
        ),
        "sio_mass_fraction": (
            ("time", "pressure"),
            sio_interp,
            {"description": "SiO mass fraction interpolated in log-pressure space", "units": "1"},
        ),
        "sio2_mass_fraction": (
            ("time", "pressure"),
            sio2_interp,
            {"description": "SiO2 mass fraction interpolated in log-pressure space", "units": "1"},
        ),
        "sio_mole_fraction": (
            ("time", "pressure"),
            sio_mole_interp,
            {"description": "SiO mole fraction interpolated in log-pressure space", "units": "1"},
        ),
        "sio2_mole_fraction": (
            ("time", "pressure"),
            sio2_mole_interp,
            {"description": "SiO2 mole fraction interpolated in log-pressure space", "units": "1"},
        ),
        "background_mole_fraction": (
            ("time", "pressure"),
            background_mole_interp,
            {"description": "Background-gas mole fraction interpolated in log-pressure space", "units": "1"},
        ),
        "mixture_mean_mol_weight": (
            ("time", "pressure"),
            mixture_mmw_interp,
            {"description": "Mixture mean molecular weight interpolated in log-pressure space", "units": "g/mol"},
        ),
    }

    return xr.Dataset(
        data_vars={name: xr.DataArray(values, dims=dims, attrs=attrs) for name, (dims, values, attrs) in data_vars.items()},
        coords={
            "time": native_ds["time"],
            "pressure": xr.DataArray(p_target_dyn_cm2.astype(float), dims=("pressure",), attrs={"units": "dyne/cm^2"}),
        },
        attrs={
            **native_ds.attrs,
            "interpolation": "log10(P) interpolation with constant extension above the GCM top and below the deepest layer",
        },
    )


def write_tp_snapshot(
    dataset: xr.Dataset,
    time_index: int,
    output_path: Path,
) -> None:
    if "pressure" in dataset.coords:
        pressure = dataset["pressure"].values
    else:
        pressure = dataset["total_pressure_dyn_cm2"].isel(time=time_index).values
    temp = dataset["temp_K"].isel(time=time_index).values

    header = "#(dyne/cm2) (K)\nPressure\tTemp\n"
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write(header)
        for p_value, t_value in zip(pressure, temp):
            handle.write(f"{p_value:.6E}\t{t_value:.6f}\n")


def write_summary_json(
    output_path: Path,
    source_file: Path,
    lat_index: int,
    lon_index: int,
    actual_lat: float,
    actual_lon: float,
    background_pressure_top_dyn_cm2: float,
    background_pressure_bottom_dyn_cm2: float,
    total_pressure_top_dyn_cm2_min: float,
    total_pressure_top_dyn_cm2_max: float,
    total_pressure_bottom_dyn_cm2_min: float,
    total_pressure_bottom_dyn_cm2_max: float,
    time_values: Iterable[float],
) -> None:
    time_values = list(time_values)
    payload = {
        "source_file": str(source_file),
        "lat_index": lat_index,
        "lon_index": lon_index,
        "selected_lat_deg": actual_lat,
        "selected_lon_deg": actual_lon,
        "background_pressure_top_dyn_cm2": background_pressure_top_dyn_cm2,
        "background_pressure_bottom_dyn_cm2": background_pressure_bottom_dyn_cm2,
        "total_pressure_top_dyn_cm2_min": total_pressure_top_dyn_cm2_min,
        "total_pressure_top_dyn_cm2_max": total_pressure_top_dyn_cm2_max,
        "total_pressure_bottom_dyn_cm2_min": total_pressure_bottom_dyn_cm2_min,
        "total_pressure_bottom_dyn_cm2_max": total_pressure_bottom_dyn_cm2_max,
        "time_start_hr": float(np.min(time_values)),
        "time_end_hr": float(np.max(time_values)),
        "n_time": int(len(time_values)),
    }
    output_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    prefix = args.output_prefix or args.input_nc.stem
    ds = load_dataset(args.input_nc)

    lat_idx, lon_idx = pick_column_indices(ds, args.target_lat, args.target_lon, args.lat_index, args.lon_index)
    column = ds[["temp", "sphum", "cld_amt"]].isel(grid_yt=lat_idx, grid_xt=lon_idx).load()
    column.attrs["source_file"] = str(args.input_nc)

    if args.clip_negative:
        column["sphum"].values = clip_nonphysical(column["sphum"].values)
        column["cld_amt"].values = clip_nonphysical(column["cld_amt"].values)

    background_pressure_dyn_cm2 = ds["pfull"].values.astype(float) * MBAR_TO_DYN_CM2
    actual_lat = float(ds["grid_yt"].values[lat_idx])
    actual_lon = float(ds["grid_xt"].values[lon_idx])
    thermo = derive_total_pressure_and_mole_fractions(
        background_pressure_dyn_cm2[np.newaxis, :],
        column["sphum"].values.astype(float),
        column["cld_amt"].values.astype(float),
        args.background_mean_mol_weight,
        args.humidity_basis,
    )
    native_ds = build_native_dataset(
        column,
        background_pressure_dyn_cm2,
        actual_lat,
        actual_lon,
        args.background_mean_mol_weight,
        args.humidity_basis,
    )
    native_path = args.output_dir / f"{prefix}_native_column.nc"
    native_ds.to_netcdf(native_path)

    summary_path = args.output_dir / f"{prefix}_selection_summary.json"
    write_summary_json(
        summary_path,
        args.input_nc,
        lat_idx,
        lon_idx,
        actual_lat,
        actual_lon,
        float(np.min(background_pressure_dyn_cm2)),
        float(np.max(background_pressure_dyn_cm2)),
        float(np.min(thermo["total_pressure_dyn_cm2"][:, 0])),
        float(np.max(thermo["total_pressure_dyn_cm2"][:, 0])),
        float(np.min(thermo["total_pressure_dyn_cm2"][:, -1])),
        float(np.max(thermo["total_pressure_dyn_cm2"][:, -1])),
        native_ds["time"].values,
    )

    target_pressure = None
    if args.vulcan_cfg is not None:
        target_pressure = load_vulcan_cfg_pressure_grid(args.vulcan_cfg)
    elif args.target_pressure_file is not None:
        target_pressure = load_pressure_grid_from_text(args.target_pressure_file)

    target_ds = None
    if target_pressure is not None:
        target_ds = build_interpolated_dataset(native_ds, np.asarray(target_pressure, dtype=float))
        target_path = args.output_dir / f"{prefix}_on_vulcan_pressure.nc"
        target_ds.to_netcdf(target_path)

    if args.write_tp_snapshot_index is not None:
        source_ds = target_ds if target_ds is not None else native_ds
        snapshot_path = args.output_dir / f"{prefix}_tp_t{args.write_tp_snapshot_index:03d}.txt"
        write_tp_snapshot(source_ds, args.write_tp_snapshot_index, snapshot_path)

    ds.close()


if __name__ == "__main__":
    main()
