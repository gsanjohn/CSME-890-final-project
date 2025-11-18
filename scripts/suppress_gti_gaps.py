import numpy as np


def suppress_gti_gaps(event_df, gti_df, original_tt_stop):
    """
    Remove gaps between GTIs, adjusting photon event times to ensure a continuous timeline.

    Parameters:
        event_df (pd.DataFrame): DataFrame containing photon events with a 'TIME' column.
        gti_df (pd.DataFrame): DataFrame with GTI intervals (START, STOP).
        original_tt_stop (float): Original end time of the observation.

    Returns:
       (tuple): Updated event_df, updated tt_stop, cumulative_gap_times
    """
    try:
        # Validate inputs
        if event_df.empty or gti_df.empty:
            raise ValueError("Event or GTI DataFrame is empty.")
        if not all(col in gti_df.columns for col in ["START", "STOP"]):
            raise ValueError("GTI DataFrame must contain 'START' and 'STOP' columns.")
        if "TIME" not in event_df.columns:
            raise ValueError("Event DataFrame must contain a 'TIME' column.")

        # Calculate gaps between GTIs
        gap_durations = gti_df["START"][1:].values - gti_df["STOP"][:-1].values
        cumulative_gap_times = np.cumsum(gap_durations)
        cumulative_gap_times = np.insert(
            cumulative_gap_times, 0, 0
        )  # Include a zero for the first GTI

        # Adjust event times
        adjusted_times = event_df["TIME"].copy()
        for i in range(1, len(gti_df)):
            # Identify events in each GTI and adjust their times
            mask = (event_df["TIME"] >= gti_df["START"].iloc[i]) & (
                event_df["TIME"] < gti_df["STOP"].iloc[i]
            )
            adjusted_times[mask] -= cumulative_gap_times[i]

        # Update event DataFrame
        updated_event_df = event_df.copy()
        updated_event_df["TIME"] = adjusted_times

        # Adjust tt_stop
        last_event_time = updated_event_df["TIME"].iloc[-1]
        updated_tt_stop = last_event_time

        # Log cumulative gaps for debugging
        print(f"    Total gap time removed: {cumulative_gap_times[-1]}")
        print(f"    New observation stop time: {updated_tt_stop}")

        return updated_event_df, updated_tt_stop, cumulative_gap_times

    except Exception as e:
        raise RuntimeError(f"Error in suppressing GTI gaps: {e}")
