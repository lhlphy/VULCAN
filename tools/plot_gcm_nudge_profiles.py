#!/usr/bin/env python3
"""
Plot chemical vertical profiles at selected times for a GCM-nudged VULCAN run.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

try:
    from tools.gcm_nudge_postprocess_common import (
        default_driver_path,
        default_plot_dir,
        default_vulcan_output_path,
        history_to_ymix,
        load_driver,
        load_vulcan_output,
        pick_species,
        resolve_repo_path,
        sanitize_for_log,
    )
except ImportError:
    from gcm_nudge_postprocess_common import (
        default_driver_path,
        default_plot_dir,
        default_vulcan_output_path,
        history_to_ymix,
        load_driver,
        load_vulcan_output,
        pick_species,
        resolve_repo_path,
        sanitize_for_log,
    )
import vulcan_cfg


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot vertical chemical profiles at selected times for a GCM-nudged VULCAN run.",
    )
    parser.add_argument(
        "--vul-file",
        type=Path,
        default=default_vulcan_output_path(),
        help="Path to the VULCAN .vul output.",
    )
    parser.add_argument(
        "--driver-file",
        type=Path,
        default=default_driver_path(),
        help="Path to the preprocessed *_on_vulcan_pressure.nc driver.",
    )
    parser.add_argument(
        "--species",
        nargs="+",
        default=None,
        help="Species to plot. Defaults to vulcan_cfg.plot_spec for compact figures.",
    )
    parser.add_argument(
        "--all-species",
        action="store_true",
        help="Plot every species above --min-mix instead of only the requested/default list.",
    )
    parser.add_argument(
        "--times-hr",
        type=float,
        nargs="+",
        default=None,
        help="Requested times in hours since the first GCM frame.",
    )
    parser.add_argument(
        "--min-mix",
        type=float,
        default=1.0e-20,
        help="Minimum mixing ratio used for species filtering and log plotting.",
    )
    parser.add_argument(
        "--ncols",
        type=int,
        default=2,
        help="Number of subplot columns.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=default_plot_dir() / "gcm_nudge_profiles.png",
        help="Output figure path.",
    )
    return parser.parse_args()


def default_times_hours(driver: dict[str, np.ndarray], final_hour: float) -> list[float]:
    driver_time_rel = driver["time_hours"] - driver["time_hours"][0]
    gcm_last = float(driver_time_rel[-1])
    hold_last_end = gcm_last + 24.0 * vulcan_cfg.gcm_nudge_hold_last_days
    candidates = [0.0, 0.25 * gcm_last, 0.5 * gcm_last, 0.75 * gcm_last, gcm_last, hold_last_end, final_hour]

    out = []
    for value in candidates:
        clipped = min(max(value, 0.0), final_hour)
        if not out or abs(clipped - out[-1]) > 1.0e-8:
            out.append(clipped)
    return out


def plot_profiles(
    vul_file: str | Path,
    driver_file: str | Path,
    output: str | Path,
    species_to_plot: list[str] | None = None,
    all_species: bool = False,
    times_hr: list[float] | tuple[float, ...] | None = None,
    min_mix: float = 1.0e-20,
    ncols: int = 2,
) -> Path:
    vul_data = load_vulcan_output(vul_file)
    driver = load_driver(driver_file)
    ymix_time, t_time_sec = history_to_ymix(vul_data)

    species = list(vul_data["variable"]["species"])
    selected_species = pick_species(species, species_to_plot, all_species, ymix_time, min_mix)
    if not selected_species:
        raise ValueError("No species passed the selection criteria.")

    t_time_hours = np.asarray(t_time_sec, dtype=float) / 3600.0
    requested_times = list(times_hr) if times_hr is not None else default_times_hours(driver, float(t_time_hours[-1]))
    time_indices = np.array([int(np.argmin(np.abs(t_time_hours - req))) for req in requested_times], dtype=int)
    actual_times = t_time_hours[time_indices]

    pressure_bar = np.asarray(vul_data["atm"]["pco"], dtype=float) / 1.0e6

    n_panels = len(selected_species)
    ncols = max(1, ncols)
    nrows = int(np.ceil(n_panels / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.0 * ncols, 4.2 * nrows), sharey=True)
    axes = np.atleast_1d(axes).ravel()

    time_colors = plt.cm.viridis(np.linspace(0.0, 1.0, len(actual_times)))

    legend_handles = []
    legend_labels = []
    for panel_idx, (ax, sp) in enumerate(zip(axes, selected_species)):
        sp_idx = species.index(sp)
        for color, time_idx, actual_time in zip(time_colors, time_indices, actual_times):
            line, = ax.plot(
                sanitize_for_log(ymix_time[time_idx, :, sp_idx], min_mix),
                pressure_bar,
                color=color,
                linewidth=1.4,
                label=f"{actual_time:.1f} hr" if panel_idx == 0 else None,
            )
            if panel_idx == 0:
                legend_handles.append(line)
                legend_labels.append(f"{actual_time:.1f} hr")

        ax.set_title(sp)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.invert_yaxis()
        ax.set_xlim(min_mix, 1.2)
        ax.grid(True, which="both", alpha=0.2)
        ax.set_xlabel("Mixing ratio")
        if panel_idx % ncols == 0:
            ax.set_ylabel("Pressure (bar)")

    for ax in axes[n_panels:]:
        ax.remove()

    fig.legend(legend_handles, legend_labels, loc="upper center", ncol=min(6, len(legend_labels)), frameon=False)
    fig.suptitle("Chemical profiles at selected times", y=0.995)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.93))

    output_path = resolve_repo_path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    print(f"Saved {output_path}")
    return output_path


def main() -> None:
    args = parse_args()
    plot_profiles(
        vul_file=args.vul_file,
        driver_file=args.driver_file,
        output=args.output,
        species_to_plot=args.species,
        all_species=args.all_species,
        times_hr=args.times_hr,
        min_mix=args.min_mix,
        ncols=args.ncols,
    )


if __name__ == "__main__":
    main()
