import pandas as pd
import numpy as np


def calculate_average_rate(
    events: pd.DataFrame,
    gti: pd.DataFrame,
    flare_gti=None,
    calculate_average_count_rate=False,
):
    """
    Calculate the average count rate during GTI intervals, optionally filtered by flare GTI.

    Parameters:
        events (pd.DataFrame): Unified photon event DataFrame with 'TIME' and 'Exposure' columns.
        gti (pd.DataFrame): Common GTI DataFrame with 'START' and 'STOP' columns.
        flare_gti (pd.DataFrame, optional): Flare GTI DataFrame with 'START' and 'STOP' columns.
        calculate_average_count_rate (bool): Whether to calculate average count rates.

    Returns:
        (pd.DataFrame): Summary DataFrame with 'START', 'STOP', and 'Count Rate'.
    """
    if not calculate_average_count_rate:
        print("Skipping average count rate calculation as per configuration.")
        return None

    try:
        # Validate inputs
        if events.empty or gti.empty:
            raise ValueError("Input events or GTI DataFrame is empty.")
        if not all(col in events.columns for col in ["TIME", "Exposure"]):
            raise ValueError(
                "Events DataFrame must contain 'TIME' and 'Exposure' columns."
            )
        if not all(col in gti.columns for col in ["START", "STOP"]):
            raise ValueError("GTI DataFrame must contain 'START' and 'STOP' columns.")

        # Use flare GTI if provided, else default to common GTI
        intervals = flare_gti if flare_gti is not None else gti

        # Results list to collect average count rates for each interval
        results = []

        for _, interval in intervals.iterrows():
            start, stop = interval["START"], interval["STOP"]
            # Filter events within the interval
            events_in_interval = events[
                (events["TIME"] >= start) & (events["TIME"] < stop)
            ]

            if not events_in_interval.empty:
                # Total corrected exposure is the sum of 1 / Exposure
                total_corrected_exposure = (1.0 / events_in_interval["Exposure"]).sum()
                # Calculate interval duration
                interval_duration = stop - start
                # Calculate average count rate
                count_rate = total_corrected_exposure / interval_duration
            else:
                count_rate = np.nan  # No events in this interval

            # Append the result
            results.append({"START": start, "STOP": stop, "Count Rate": count_rate})

        # Convert results to a DataFrame
        results_df = pd.DataFrame(results)

        # Log output
        print(f"    Processed {len(results_df)} intervals.")
        print(
            f"    Count rate range: {results_df['Count Rate'].min()} - {results_df['Count Rate'].max()}"
        )
        return results_df

    except Exception as e:
        raise RuntimeError(f"Error in calculating average rate: {e}")
