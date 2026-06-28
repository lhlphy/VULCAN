#!/usr/bin/env python3
"""
Shared helpers for GCM-nudged VULCAN post-processing.
"""

from __future__ import annotations

import pickle
import sys
from pathlib import Path

import numpy as np
import xarray as xr


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import vulcan_cfg


BAR_TO_DYN_CM2 = 1.0e6


def resolve_repo_path(raw_path: str | Path) -> Path:
    path = Path(raw_path)
    if path.is_absolute():
        return path
    return (REPO_ROOT / path).resolve()


def default_vulcan_output_path() -> Path:
    return resolve_repo_path(Path(vulcan_cfg.output_dir) / vulcan_cfg.out_name)


def default_driver_path() -> Path:
    return resolve_repo_path(vulcan_cfg.gcm_nudge_file)


def default_plot_dir() -> Path:
    return resolve_repo_path(vulcan_cfg.plot_dir)


def load_vulcan_output(vul_path: str | Path) -> dict:
    with resolve_repo_path(vul_path).open("rb") as handle:
        return pickle.load(handle)


def require_history(vul_data: dict) -> tuple[np.ndarray, np.ndarray]:
    variable = vul_data["variable"]
    if "y_time" not in variable or "t_time" not in variable:
        raise ValueError(
            "The selected .vul file does not contain y_time/t_time. "
            "Run VULCAN with save_evolution=True first."
        )

    y_time = np.asarray(variable["y_time"], dtype=float)
    t_time = np.asarray(variable["t_time"], dtype=float)
    if y_time.ndim != 3:
        raise ValueError(f"Unexpected y_time shape {y_time.shape}; expected (nt, nz, nsp).")
    if t_time.ndim != 1 or t_time.shape[0] != y_time.shape[0]:
        raise ValueError("t_time and y_time lengths do not match.")
    return y_time, t_time


def history_to_ymix(vul_data: dict) -> tuple[np.ndarray, np.ndarray]:
    y_time, t_time = require_history(vul_data)
    n0 = np.asarray(vul_data["atm"]["n_0"], dtype=float)
    ymix_time = y_time / n0[np.newaxis, :, np.newaxis]
    return ymix_time, t_time


def load_driver(driver_path: str | Path) -> dict[str, np.ndarray]:
    resolved = resolve_repo_path(driver_path)
    with xr.open_dataset(resolved, decode_times=False) as ds:
        return {
            "path": str(resolved),
            "time_hours": np.asarray(ds["time"].values, dtype=float),
            "pressure_dyn_cm2": np.asarray(ds["pressure"].values, dtype=float),
            "temp_K": np.asarray(ds["temp_K"].values, dtype=float),
            "sio_mass_fraction": np.asarray(ds["sio_mass_fraction"].values, dtype=float),
            "sio2_mass_fraction": np.asarray(ds["sio2_mass_fraction"].values, dtype=float),
        }


def driver_time_hours_relative(driver: dict[str, np.ndarray]) -> np.ndarray:
    time_hours = np.asarray(driver["time_hours"], dtype=float)
    return time_hours - time_hours[0]


def interpolate_with_hold(
    source_time_hours: np.ndarray,
    source_values: np.ndarray,
    target_time_hours: np.ndarray,
    hold_last_days: float | None = None,
) -> np.ndarray:
    source_time_hours = np.asarray(source_time_hours, dtype=float)
    source_values = np.asarray(source_values, dtype=float)
    target_time_hours = np.asarray(target_time_hours, dtype=float)
    hold_last_days = vulcan_cfg.gcm_nudge_hold_last_days if hold_last_days is None else hold_last_days

    out = np.full(target_time_hours.shape, np.nan, dtype=float)
    if source_time_hours.size == 0:
        return out

    in_gcm = target_time_hours <= source_time_hours[-1]
    if np.any(in_gcm):
        out[in_gcm] = np.interp(target_time_hours[in_gcm], source_time_hours, source_values)

    hold_end_hours = source_time_hours[-1] + hold_last_days * 24.0
    in_hold = (target_time_hours > source_time_hours[-1]) & (target_time_hours <= hold_end_hours)
    if np.any(in_hold):
        out[in_hold] = source_values[-1]

    return out


def nearest_pressure_indices(
    pressure_dyn_cm2: np.ndarray,
    levels_bar: list[float] | tuple[float, ...] | np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    pressure_dyn_cm2 = np.asarray(pressure_dyn_cm2, dtype=float)
    target_dyn_cm2 = np.asarray(levels_bar, dtype=float) * BAR_TO_DYN_CM2
    indices = np.array(
        [int(np.argmin(np.abs(pressure_dyn_cm2 - target))) for target in target_dyn_cm2],
        dtype=int,
    )
    actual_bar = pressure_dyn_cm2[indices] / BAR_TO_DYN_CM2
    return indices, actual_bar


def pick_species(
    species: list[str],
    requested_species: list[str] | None,
    all_species: bool,
    ymix_time: np.ndarray,
    min_mix: float,
) -> list[str]:
    if requested_species:
        missing = [sp for sp in requested_species if sp not in species]
        if missing:
            raise ValueError("Requested species are not in the output: " + ", ".join(missing))
        return requested_species

    if all_species:
        species_arr = np.asarray(species)
        active_mask = np.max(np.where(ymix_time > 0.0, ymix_time, 0.0), axis=(0, 1)) >= min_mix
        return species_arr[active_mask].tolist()

    default_species = [sp for sp in vulcan_cfg.plot_spec if sp in species]
    if default_species:
        return default_species

    species_arr = np.asarray(species)
    active_mask = np.max(np.where(ymix_time > 0.0, ymix_time, 0.0), axis=(0, 1)) >= min_mix
    return species_arr[active_mask].tolist()


def sanitize_for_log(values: np.ndarray, floor: float) -> np.ndarray:
    arr = np.asarray(values, dtype=float).copy()
    arr[arr < floor] = np.nan
    return arr


def dedupe_legend_entries(handles: list, labels: list[str]) -> tuple[list, list[str]]:
    seen: set[str] = set()
    out_handles = []
    out_labels = []
    for handle, label in zip(handles, labels):
        if not label or label in seen:
            continue
        seen.add(label)
        out_handles.append(handle)
        out_labels.append(label)
    return out_handles, out_labels


def format_pressure_label(requested_bar: float, actual_bar: float) -> str:
    return f"target {requested_bar:.1e} bar\nactual {actual_bar:.2e} bar"
