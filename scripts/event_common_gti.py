import pandas as pd


def filter_events_with_common_gti(events, gti):
    """
    Filter events to retain those within the merged (common) GTIs.

    Parameters:
        events (pd.DataFrame): Event DataFrame with 'TIME' column.
        gti (pd.DataFrame): Merged GTI DataFrame with 'START' and 'STOP' columns.

    Returns:
        (pd.DataFrame): Filtered event DataFrame with valid observation times.
    """
    try:
        filtered_events = []

        for _, interval in gti.iterrows():
            start, stop = interval["START"], interval["STOP"]
            events_in_gti = events[(events["TIME"] >= start) & (events["TIME"] <= stop)]
            filtered_events.append(events_in_gti)

        filtered_events = pd.concat(filtered_events, ignore_index=True)

        print(f"    Original events: {len(events)}")
        print(f"    Filtered events: {len(filtered_events)}")

        return filtered_events
    except Exception as e:
        raise RuntimeError(f"Error during event filtering: {e}")
