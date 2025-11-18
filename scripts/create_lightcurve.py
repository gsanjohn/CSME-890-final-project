import numpy as np
import pandas as pd
from scripts.find_blocks import upper_limit_gehrels, lower_limit_gehrels
from scripts.get_binned_corr_factor import get_binned_corr_factor


def generate_lightcurve(
    events_df: pd.DataFrame,
    gti_df: pd.DataFrame,
    binsize=100,
    energy_range=(3.0, 79.0),
    gti_average=True,
):
    """
    Generate a regularly binned light curve from filtered event data.

    Parameters:
        events_df (pd.DataFrame): Filtered events with columns ['TIME', 'Energy', 'Exposure'].
        gti_df (pd.DataFrame): GTI intervals with columns ['START', 'STOP'].
        binsize (int): Time bin size in seconds (default: 100).
        energy_range (tuple): (Emin, Emax) energy range in keV (default: (3.0, 79.0)).
        gti_average (bool): Use GTI-based binning if True (default: True).

    Returns:
        (pd.DataFrame): Regularly binned light curve with columns:
                      ['bin_start', 'bin_end', 'count_rate', 'upper_limit', 'lower_limit'].
    """
    # Validate required columns
    required_columns = ["TIME", "Energy"]
    for col in required_columns:
        if col not in events_df.columns:
            raise ValueError(f"Missing required column: {col}")

    # Filter events by energy range
    filtered_events = events_df[
        (events_df["Energy"] >= energy_range[0])
        & (events_df["Energy"] <= energy_range[1])
    ]
    if filtered_events.empty:
        raise ValueError("No events found in the specified energy range.")

    # Define time bins based on GTIs and events
    if gti_average:
        bins = []
        for _, gti in gti_df.iterrows():
            bins += list(np.arange(gti["START"], gti["STOP"], binsize))
        bins = sorted(set(bins))
    else:
        observation_start = gti_df["START"].min()
        observation_stop = max(
            gti_df["STOP"].max(), events_df["TIME"].max()
        )  # Extend to max event time
        bins = np.arange(observation_start, observation_stop + binsize, binsize)

    # Bin the events
    bin_edges = np.array(bins)
    event_counts, _ = np.histogram(filtered_events["TIME"], bins=bin_edges)

    # Calculate bin durations
    bin_durations = np.diff(bin_edges)
    count_rates = event_counts / bin_durations

    # Confidence intervals using Gehrels' method
    upper_limits = upper_limit_gehrels(0.8413, event_counts) / bin_durations
    lower_limits = lower_limit_gehrels(0.8413, event_counts) / bin_durations

    # Set count_rate, upper_limit, and lower_limit to NaN for rows where count_rate = 0.0
    mask_zero_counts = event_counts == 0
    count_rates[mask_zero_counts] = np.nan
    upper_limits[mask_zero_counts] = np.nan
    lower_limits[mask_zero_counts] = np.nan

    # Remove bins outside GTIs (allow partial overlaps)
    valid_bins = []
    for start, stop in zip(bin_edges[:-1], bin_edges[1:]):
        if any((gti_df["START"] < stop) & (gti_df["STOP"] > start)):
            valid_bins.append((start, stop))

    # Correct for Livetime, vignetting, and PSF completeness for each light curve bin
    correction_factor = get_binned_corr_factor(events_df, valid_bins)

    count_rates[: len(valid_bins)] /= correction_factor
    upper_limits[: len(valid_bins)] /= correction_factor
    lower_limits[: len(valid_bins)] /= correction_factor

    # Format as DataFrame
    valid_bins = np.array(valid_bins)
    lightcurve_df = pd.DataFrame(
        {
            "bin_start": valid_bins[:, 0],
            "bin_end": valid_bins[:, 1],
            "count_rate": count_rates[: len(valid_bins)],
            "upper_limit": upper_limits[: len(valid_bins)],
            "lower_limit": lower_limits[: len(valid_bins)],
        }
    )

    # Filter out bins longer than the specified binsize - trim in between GTI bins
    lightcurve_df = lightcurve_df[
        (lightcurve_df["bin_end"] - lightcurve_df["bin_start"]) <= binsize
    ]

    print(f"    Light Curve Time Range: {bin_edges[0]} - {bin_edges[-1]}")
    print(
        f"    Event Time Range: {events_df['TIME'].min()} - {events_df['TIME'].max()}"
    )
    print(f"    Removed {len(valid_bins) - len(lightcurve_df)} oversized bins.")

    return lightcurve_df


def plot_prep(gti_df, bb_df):
    plot_frame = gti_df.copy()
    # Filter events in the current block
    if not bb_df.empty:
        for a, row in gti_df.iterrows():
            row_start = row["start"]
            row_stop = row["stop"]

            for _, block in bb_df.iterrows():
                if row_stop <= block["stop"]:
                    if row_start >= block["start"]:
                        plot_frame["rate"][a] = block["rate"]
                        plot_frame["upperlim"][a] = block["upperlim"]
                        plot_frame["lowerlim"][a] = block["lowerlim"]

    return plot_frame
