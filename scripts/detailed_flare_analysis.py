from scripts.expo_events import bayesian_blocks
import pandas as pd
# import matplotlib.pyplot as plt


def detailed_flare_analysis(
    event_df, bayesian_blocks_df, p0=0.1, flare_analysis_flag=False
):
    """
    Perform a detailed Bayesian Block analysis on the flaring activity block.

    Parameters:
        event_df (pd.DataFrame): Event DataFrame with 'TIME' column.
        bayesian_blocks_df (pd.DataFrame): Results from Step 10 (Bayesian Block Analysis).
        p0 (float): False positive rate for detailed segmentation.
        flare_analysis_flag (bool): Whether to perform detailed flare analysis.

    Returns:
        (pd.DataFrame): Refined Bayesian Block DataFrame for the flaring activity.
    """
    if not flare_analysis_flag:
        print("Skipping detailed flare analysis as per configuration.")
        return None

    try:
        # Identify the flaring block
        flare_block = bayesian_blocks_df.loc[
            bayesian_blocks_df["rate"].idxmax()
        ]  # Assuming flaring block is where the event rate is maximum

        flaring_start, flaring_stop = flare_block["start"], flare_block["stop"]
        print(
            f"Flaring activity identified between {flaring_start} and {flaring_stop}."
        )

        # Extract events in the flaring block
        flaring_events = event_df[
            (event_df["TIME"] >= flaring_start) & (event_df["TIME"] < flaring_stop)
        ]
        if flaring_events.empty:
            print(
                "No events found in the flaring interval. Skipping detailed analysis."
            )
            return None

        # Apply Bayesian segmentation with a maybe higher p0
        flaring_times = flaring_events["TIME"].values
        detailed_change_points = bayesian_blocks(flaring_times, fitness="events", p0=p0)

        # Calculate detailed block statistics
        block_start = detailed_change_points[:-1]
        block_stop = detailed_change_points[1:]
        durations = block_stop - block_start

        # Count events in each interval
        counts = (
            pd.cut(flaring_times, bins=detailed_change_points)
            .value_counts()
            .sort_index()
            .values
        )
        rates = counts / durations

        # Create detailed Bayesian block DataFrame
        detailed_blocks_df = pd.DataFrame(
            {
                "start": block_start,
                "stop": block_stop,
                "duration": durations,
                "counts": counts,
                "rate": rates,
                "fp_rate": p0,
            }
        )

        # Log results
        # print(f"Refined Bayesian Blocks for Flaring Activity:")
        # print(detailed_blocks_df)

        # Plot the detailed segmentation
        # plt.figure(figsize=(10, 6))
        # plt.hist(flaring_times, bins=100, histtype='step', color='gray', label='Photon Events')
        # for cp in detailed_change_points:
        #     plt.axvline(cp, color='r', linestyle='--', label='Detailed Change Point' if 'Detailed Change Point' not in plt.gca().get_legend_handles_labels()[1] else None)
        # plt.title("Detailed Bayesian Segmentation of Flaring Activity")
        # plt.xlabel("Time")
        # plt.ylabel("Event Counts")
        # plt.legend()
        # plt.show()

        return detailed_blocks_df

    except Exception as e:
        raise RuntimeError(f"Error in detailed flare analysis: {e}")
