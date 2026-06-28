"""
Runtime helpers for the dedicated GCM SiO/SiO2 nudging mode.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import xarray as xr

import vulcan_cfg


REQUIRED_DRIVER_VARS = (
    "temp_K",
    "background_pressure_dyn_cm2",
    "sio_mole_fraction",
    "sio2_mole_fraction",
    "background_mole_fraction",
)

_DRIVER_CACHE: dict[Path, dict[str, np.ndarray | float | str]] = {}


def _resolve_driver_path(path: str | Path | None = None) -> Path:
    raw_path = Path(path or vulcan_cfg.gcm_nudge_file)
    if raw_path.is_absolute():
        return raw_path
    return (Path(__file__).resolve().parent / raw_path).resolve()


def _validate_required_coords(ds: xr.Dataset) -> None:
    missing = [name for name in ("time", "pressure") if name not in ds.coords]
    if missing:
        raise IOError("Missing required GCM nudging coordinates: " + ", ".join(missing))


def _extract_driver_array(ds: xr.Dataset, name: str) -> np.ndarray:
    if name not in ds:
        raise IOError(f"Missing required GCM nudging variable: {name}")
    return np.asarray(ds[name].transpose("time", "pressure").values, dtype=float)


def _validate_driver_shapes(driver: dict[str, np.ndarray | float | str]) -> None:
    time_size = driver["time_raw_hours"].shape[0]
    pressure_size = driver["pressure"].shape[0]
    for name in REQUIRED_DRIVER_VARS:
        values = driver[name]
        if values.shape != (time_size, pressure_size):
            raise IOError(
                f"Driver variable {name} has shape {values.shape}, "
                f"expected {(time_size, pressure_size)}."
            )


def _validate_monotonic(driver: dict[str, np.ndarray | float | str]) -> None:
    time_sec = driver["time_sec"]
    pressure = driver["pressure"]
    if time_sec.ndim != 1 or pressure.ndim != 1:
        raise IOError("Driver time and pressure coordinates must be one-dimensional.")
    if np.any(np.diff(time_sec) <= 0.0):
        raise IOError("Driver time coordinate must be strictly increasing.")
    if np.any(pressure <= 0.0):
        raise IOError("Driver pressure coordinate must stay positive.")
    if not np.all(np.diff(pressure) < 0.0):
        raise IOError("Driver pressure coordinate must be strictly decreasing from bottom to top.")


def load_gcm_nudge_driver(path: str | Path | None = None) -> dict[str, np.ndarray | float | str]:
    resolved_path = _resolve_driver_path(path)
    if resolved_path in _DRIVER_CACHE:
        return _DRIVER_CACHE[resolved_path]

    print(f"Loading GCM nudging driver once from {resolved_path}")
    with xr.open_dataset(resolved_path, decode_times=False) as ds:
        _validate_required_coords(ds)
        driver = {
            "source_file": str(resolved_path),
            "time_raw_hours": np.asarray(ds["time"].values, dtype=float),
            "pressure": np.asarray(ds["pressure"].values, dtype=float),
        }
        for name in REQUIRED_DRIVER_VARS:
            driver[name] = _extract_driver_array(ds, name)

    driver["time_sec"] = (driver["time_raw_hours"] - driver["time_raw_hours"][0]) * 3600.0
    driver["frame_dt_sec"] = np.diff(driver["time_sec"])

    _validate_driver_shapes(driver)
    _validate_monotonic(driver)

    if driver["frame_dt_sec"].size == 0:
        raise IOError("GCM nudging driver needs at least two time samples.")

    _DRIVER_CACHE[resolved_path] = driver
    return driver


def validate_pressure_grid(expected_pressure: np.ndarray, path: str | Path | None = None) -> dict[str, np.ndarray | float | str]:
    driver = load_gcm_nudge_driver(path)
    expected_pressure = np.asarray(expected_pressure, dtype=float)
    if expected_pressure.shape != driver["pressure"].shape:
        raise IOError(
            f"Driver pressure grid shape {driver['pressure'].shape} does not match "
            f"VULCAN grid shape {expected_pressure.shape}."
        )
    if not np.allclose(expected_pressure, driver["pressure"], rtol=1.0e-10, atol=1.0e-20):
        max_diff = np.max(np.abs(expected_pressure - driver["pressure"]))
        raise IOError(
            "Driver pressure grid does not match the configured VULCAN grid. "
            f"Max abs diff = {max_diff:.3e} dyne/cm^2."
        )
    return driver


def get_nudge_setup(
    tau_fraction: float | None = None,
    hold_last_days: float | None = None,
    free_run_days: float | None = None,
    path: str | Path | None = None,
) -> dict[str, np.ndarray | float | str]:
    driver = load_gcm_nudge_driver(path)
    tau_fraction = vulcan_cfg.gcm_nudge_tau_fraction if tau_fraction is None else tau_fraction
    hold_last_days = vulcan_cfg.gcm_nudge_hold_last_days if hold_last_days is None else hold_last_days
    free_run_days = vulcan_cfg.gcm_nudge_free_run_days if free_run_days is None else free_run_days

    last_frame_time = float(driver["time_sec"][-1])
    hold_last_end = last_frame_time + hold_last_days * 86400.0
    timed_run_end = hold_last_end + free_run_days * 86400.0

    return {
        **driver,
        "tau_fraction": float(tau_fraction),
        "hold_last_end_sec": float(hold_last_end),
        "timed_run_end_sec": float(timed_run_end),
        "min_frame_dt_sec": float(np.min(driver["frame_dt_sec"])),
        "last_frame_dt_sec": float(driver["frame_dt_sec"][-1]),
        "last_frame_time_sec": float(last_frame_time),
    }


def interpolate_nudge_targets(
    t_eval_sec: float,
    tau_fraction: float | None = None,
    hold_last_days: float | None = None,
    free_run_days: float | None = None,
    path: str | Path | None = None,
) -> dict[str, np.ndarray | float | bool | str]:
    setup = get_nudge_setup(
        tau_fraction=tau_fraction,
        hold_last_days=hold_last_days,
        free_run_days=free_run_days,
        path=path,
    )
    time_sec = setup["time_sec"]
    frame_dt_sec = setup["frame_dt_sec"]
    tau_fraction = setup["tau_fraction"]

    if t_eval_sec <= 0.0:
        return {
            "active": True,
            "phase": "gcm",
            "tau_sec": frame_dt_sec[0] * tau_fraction,
            "sio_target": setup["sio_mole_fraction"][0].copy(),
            "sio2_target": setup["sio2_mole_fraction"][0].copy(),
        }

    if t_eval_sec <= setup["last_frame_time_sec"]:
        upper = int(np.searchsorted(time_sec, t_eval_sec, side="right"))
        lower = max(upper - 1, 0)
        upper = min(upper, time_sec.size - 1)

        if lower == upper:
            lower = max(upper - 1, 0)

        dt_frame = time_sec[upper] - time_sec[lower]
        weight = 0.0 if dt_frame == 0.0 else (t_eval_sec - time_sec[lower]) / dt_frame
        sio_target = (1.0 - weight) * setup["sio_mole_fraction"][lower] + weight * setup["sio_mole_fraction"][upper]
        sio2_target = (1.0 - weight) * setup["sio2_mole_fraction"][lower] + weight * setup["sio2_mole_fraction"][upper]
        return {
            "active": True,
            "phase": "gcm",
            "tau_sec": dt_frame * tau_fraction,
            "sio_target": sio_target,
            "sio2_target": sio2_target,
        }

    if t_eval_sec <= setup["hold_last_end_sec"]:
        return {
            "active": True,
            "phase": "hold_last",
            "tau_sec": setup["last_frame_dt_sec"] * tau_fraction,
            "sio_target": setup["sio_mole_fraction"][-1].copy(),
            "sio2_target": setup["sio2_mole_fraction"][-1].copy(),
        }

    return {
        "active": False,
        "phase": "free_run",
        "tau_sec": np.nan,
        "sio_target": None,
        "sio2_target": None,
    }
