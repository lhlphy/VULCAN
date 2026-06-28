#!/usr/bin/env python3
"""
Plot time-dependent mixing ratios at selected pressure levels for a GCM-nudged VULCAN run.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

try:
    from tools.gcm_nudge_postprocess_common import (
        dedupe_legend_entries,
        default_driver_path,
        default_plot_dir,
        default_vulcan_output_path,
        driver_time_hours_relative,
        format_pressure_label,
        history_to_ymix,
        interpolate_with_hold,
        load_driver,
        load_vulcan_output,
        nearest_pressure_indices,
        pick_species,
        resolve_repo_path,
        sanitize_for_log,
    )
except ImportError:
    from gcm_nudge_postprocess_common import (
        dedupe_legend_entries,
        default_driver_path,
        default_plot_dir,
        default_vulcan_output_path,
        driver_time_hours_relative,
        format_pressure_label,
        history_to_ymix,
        interpolate_with_hold,
        load_driver,
        load_vulcan_output,
        nearest_pressure_indices,
        pick_species,
        resolve_repo_path,
        sanitize_for_log,
    )
import vulcan_cfg


DEFAULT_LEVELS_BAR = [10.0, 1.0, 1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot VULCAN species evolution at selected pressure levels and overlay GCM sphum/cld_amt/temp.",
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
        "--levels-bar",
        type=float,
        nargs="+",
        default=DEFAULT_LEVELS_BAR,
        help="Pressure levels in bar to plot.",
    )
    parser.add_argument(
        "--species",
        nargs="+",
        default=None,
        help="Species to plot. Defaults to vulcan_cfg.plot_spec for manageable figures.",
    )
    parser.add_argument(
        "--all-species",
        action="store_true",
        help="Plot every species above --min-mix instead of only the requested/default list.",
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
        default=4,
        help="Number of subplot columns.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=default_plot_dir() / "gcm_nudge_time_series.png",
        help="Output figure path.",
    )
    return parser.parse_args()


def plot_time_series(
    vul_file: str | Path,
    driver_file: str | Path,
    output: str | Path,
    levels_bar: list[float] | tuple[float, ...] | None = None,
    species_to_plot: list[str] | None = None,
    all_species: bool = False,
    min_mix: float = 1.0e-20,
    ncols: int = 4,
) -> Path:
    levels_bar = list(DEFAULT_LEVELS_BAR if levels_bar is None else levels_bar)

    vul_data = load_vulcan_output(vul_file)
    driver = load_driver(driver_file)
    ymix_time, t_time_sec = history_to_ymix(vul_data)

    species = list(vul_data["variable"]["species"])
    selected_species = pick_species(species, species_to_plot, all_species, ymix_time, min_mix)
    if not selected_species:
        raise ValueError("No species passed the selection criteria.")

    pressure_dyn_cm2 = np.asarray(vul_data["atm"]["pco"], dtype=float)
    level_indices, actual_levels_bar = nearest_pressure_indices(pressure_dyn_cm2, levels_bar)

    t_time_hours = np.asarray(t_time_sec, dtype=float) / 3600.0
    driver_time_hours = driver_time_hours_relative(driver)
    gcm_last_hour = float(driver_time_hours[-1])
    hold_last_end_hour = gcm_last_hour + 24.0 * vulcan_cfg.gcm_nudge_hold_last_days

    n_panels = len(levels_bar)
    ncols = max(1, ncols)
    nrows = int(np.ceil(n_panels / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.7 * ncols, 3.7 * nrows), sharex=True)
    axes = np.atleast_1d(axes).ravel()

    species_colors = plt.cm.tab20(np.linspace(0.0, 1.0, max(len(selected_species), 2)))
    driver_styles = {
        "sphum": {"color": "black", "linestyle": "--", "linewidth": 1.5},
        "cld_amt": {"color": "dimgray", "linestyle": "--", "linewidth": 1.5},
        "temp": {"color": "firebrick", "linestyle": "--", "linewidth": 1.5},
    }

    legend_handles = []
    legend_labels = []

    for panel_idx, (ax, target_bar, actual_bar, level_idx) in enumerate(
        zip(axes, levels_bar, actual_levels_bar, level_indices)
    ):
        for color, sp in zip(species_colors, selected_species):
            sp_idx = species.index(sp)
            line, = ax.plot(
                t_time_hours,
                sanitize_for_log(ymix_time[:, level_idx, sp_idx], min_mix),
                color=color,
                linewidth=1.2,
                label=sp if panel_idx == 0 else None,
            )
            if panel_idx == 0:
                legend_handles.append(line)
                legend_labels.append(sp)

        sio_series = interpolate_with_hold(
            driver_time_hours,
            driver["sio_mass_fraction"][:, level_idx],
            t_time_hours,
        )
        sio2_series = interpolate_with_hold(
            driver_time_hours,
            driver["sio2_mass_fraction"][:, level_idx],
            t_time_hours,
        )

        line_sio, = ax.plot(
            t_time_hours,
            sanitize_for_log(sio_series, min_mix),
            label="sphum (GCM)" if panel_idx == 0 else None,
            **driver_styles["sphum"],
        )
        line_sio2, = ax.plot(
            t_time_hours,
            sanitize_for_log(sio2_series, min_mix),
            label="cld_amt (GCM)" if panel_idx == 0 else None,
            **driver_styles["cld_amt"],
        )

        ax_temp = ax.twinx()
        temp_series = interpolate_with_hold(
            driver_time_hours,
            driver["temp_K"][:, level_idx],
            t_time_hours,
        )
        line_temp, = ax_temp.plot(
            t_time_hours,
            temp_series,
            label="temp (GCM)" if panel_idx == 0 else None,
            **driver_styles["temp"],
        )

        if panel_idx == 0:
            legend_handles.extend([line_sio, line_sio2, line_temp])
            legend_labels.extend(["sphum (GCM)", "cld_amt (GCM)", "temp (GCM)"])

        ax.axvline(gcm_last_hour, color="0.65", linestyle=":", linewidth=1.0)
        ax.axvline(hold_last_end_hour, color="0.75", linestyle=":", linewidth=1.0)

        ax.set_title(format_pressure_label(target_bar, actual_bar))
        ax.set_yscale("log")
        ax.set_ylim(min_mix, 1.2)
        ax.set_xlim(0.0, max(float(t_time_hours[-1]), hold_last_end_hour))
        ax.grid(True, which="both", alpha=0.2)

        if panel_idx % ncols == 0:
            ax.set_ylabel("Mixing ratio / mass fraction")
        if panel_idx >= (nrows - 1) * ncols:
            ax.set_xlabel("Time since first GCM frame (hr)")
        if panel_idx == ncols - 1:
            ax_temp.set_ylabel("Temperature (K)")

    for ax in axes[n_panels:]:
        ax.remove()

    legend_handles, legend_labels = dedupe_legend_entries(legend_handles, legend_labels)
    fig.legend(legend_handles, legend_labels, loc="upper center", ncol=min(6, len(legend_labels)), frameon=False)
    fig.suptitle(
        "GCM-nudged time evolution by pressure level\n"
        "Vertical dotted lines: GCM end and end of 1-day hold-last forcing",
        y=0.995,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.93))

    output_path = resolve_repo_path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    print(f"Saved {output_path}")
    return output_path


def main() -> None:
    args = parse_args()
    plot_time_series(
        vul_file=args.vul_file,
        driver_file=args.driver_file,
        output=args.output,
        levels_bar=args.levels_bar,
        species_to_plot=args.species,
        all_species=args.all_species,
        min_mix=args.min_mix,
        ncols=args.ncols,
    )


if __name__ == "__main__":
    main()
