import pandas as pd
import numpy as np


def clean_gti(
    gti: pd.DataFrame, events: pd.DataFrame, threshold=30, starttrim=10, stoptrim=10
):
    """
    Clean GTIs by trimming and filtering based on duration and events.

    Parameters:
        gti (pd.DataFrame): GTI DataFrame with 'START' and 'STOP' columns.
        events (pd.DataFrame): Events DataFrame with 'TIME' column.
        threshold (float): Minimum duration for a valid GTI (seconds).
        starttrim (float): Seconds to trim from GTI start times.
        stoptrim (float): Seconds to trim from GTI stop times.

    Returns:
        (pd.DataFrame): Cleaned GTI DataFrame.
    """
    try:
        # Trim GTI intervals
        gti["START"] += starttrim
        gti["STOP"] -= stoptrim

        # Calculate GTI durations
        gti["DURATION"] = gti["STOP"] - gti["START"]

        # Filter GTIs by duration
        valid_gti = gti[gti["DURATION"] > threshold]

        # Count events within each GTI
        valid_gti = valid_gti.assign(
            EVENT_COUNT=valid_gti.apply(
                lambda row: (
                    (events["TIME"] >= row["START"]) & (events["TIME"] <= row["STOP"])
                ).sum(),
                axis=1,
            )
        )

        # Retain GTIs with at least one event
        final_gti = valid_gti[valid_gti["EVENT_COUNT"] > 0]

        # Log changes for debugging
        print(f"    GTIs before cleaning: {len(gti)}")
        print(f"    GTIs after duration filter: {len(valid_gti)}")
        print(f"    GTIs after event filter: {len(final_gti)}")

        # Drop intermediate columns and return
        return final_gti  # for debugging use
        # return final_gti.drop(columns=['DURATION', 'EVENT_COUNT'])
    except Exception as e:
        raise RuntimeError(f"Error during GTI cleaning: {e}")
