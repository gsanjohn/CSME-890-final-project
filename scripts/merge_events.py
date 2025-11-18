import pandas as pd


def merge_events(eventsA, eventsB, exposureA, exposureB):
    """
    Merge photon event data from modules A and B into a unified dataset.

    Parameters:
        eventsA (pd.DataFrame): Filtered event data for module A, including TIME and Energy.
        eventsB (pd.DataFrame): Filtered event data for module B, including TIME and Energy.
        exposureA (np.ndarray): Correction factors for module A events.
        exposureB (np.ndarray): Correction factors for module B events.

    Returns:
        (pd.DataFrame): Unified photon event DataFrame with TIME, Energy, Exposure, and Module columns.
    """
    try:
        # Add correction factors and module labels to each DataFrame
        eventsA["Exposure"] = exposureA
        eventsA["Module"] = "A"

        eventsB["Exposure"] = exposureB
        eventsB["Module"] = "B"

        # Concatenate the DataFrames
        events_combined = pd.concat([eventsA, eventsB], ignore_index=True)

        # Sort by TIME
        events_combined = events_combined.sort_values(by="TIME").reset_index(drop=True)

        # Log statistics for the merged DataFrame
        print(f"    Number of events Module A: {len(eventsA)}")
        print(f"    Number of events Module B: {len(eventsB)}")
        print(f"    Total events after merging: {len(events_combined)}")
        print(
            f"    Time range: {events_combined['TIME'].min()} - {events_combined['TIME'].max()}"
        )
        print(
            f"    Energy range: {events_combined['Energy'].min()} - {events_combined['Energy'].max()}"
        )

        return events_combined

    except Exception as e:
        raise RuntimeError(f"Error in merging events: {e}")


def merge_event(events, exposure):
    """
    Merge photon event data from one module into a unified dataset.

    Parameters:
        events (pd.DataFrame): Filtered event data for one module, including TIME and Energy..
        exposure (np.ndarray): Correction factors for one module events.

    Returns:
        (pd.DataFrame): Unified photon event DataFrame with TIME, Energy, Exposure, and Module columns.
    """
    try:
        # Add correction factors and module labels to each DataFrame
        events["Exposure"] = exposure

        # Debugging:
        # print("\n1")
        # # Check the lengths of events and exposures
        # print(f"Length of eventsA: {len(eventsA)}, Length of exposureA: {len(exposureA)}")
        # if len(eventsA) != len(exposureA):
        #     print("Warning: Mismatch detected between eventsA and exposureA!")

        # print(f"Length of eventsB: {len(eventsB)}, Length of exposureB: {len(exposureB)}")
        # if len(eventsB) != len(exposureB):
        #     print("Warning: Mismatch detected between eventsB and exposureB!")

        # Concatenate the DataFrames
        events_combined = events

        # Debugging:
        # print("\n2")
        # print(f"Rows in Module A before concatenation: {len(eventsA)}")
        # print(f"Rows in Module B before concatenation: {len(eventsB)}")
        # print(f"Total rows after concatenation: {len(events_combined)}")

        # Sort by TIME
        events_combined = events_combined.sort_values(by="TIME").reset_index(drop=True)

        # Debugging:
        # print("\n3")
        # print(f"Number of rows before sorting: {len(events_combined)}")
        # print(f"Number of NaN TIME values: {events_combined['TIME'].isna().sum()}\n")

        # Log statistics for the merged DataFrame
        print(f"    Number of events Module: {len(events)}")
        # print(f"    Number of events Module B: {len(eventsB)}")
        print(f"    Total events after merging: {len(events_combined)}")
        print(
            f"    Time range: {events_combined['TIME'].min()} - {events_combined['TIME'].max()}"
        )
        print(
            f"    Energy range: {events_combined['Energy'].min()} - {events_combined['Energy'].max()}"
        )

        return events_combined

    except Exception as e:
        raise RuntimeError(f"Error in merging events: {e}")
