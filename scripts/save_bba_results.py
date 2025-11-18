import pandas as pd
from datetime import datetime, timedelta
import numpy as np
from tabulate import tabulate


def convert_nustar_to_utc(nustar_time):
    """
    Convert NuSTAR mission time to UTC.

    Parameters:
        nustar_time (float): NuSTAR mission time in seconds since the reference epoch.

    Returns:
        (str): UTC timestamp in ISO format (e.g., '2018-04-24T03:06:32').
    """
    nustar_epoch = datetime(2010, 1, 1)  # NuSTAR reference epoch
    utc_time = nustar_epoch + timedelta(seconds=nustar_time)
    return utc_time  # .isoformat()


def calculate_event_counts_and_rates(bba_df, event_df, exposure_col="Exposure"):
    """
    Calculate the number of events (`Counts`), the length of each block (`Length (s)`),
    and the event rate (`Rate (cts/s)`) for each Bayesian Block interval.

    Parameters:
        bba_df (pd.DataFrame): DataFrame containing BBA block intervals with 'start' and 'stop' columns.
        event_df (pd.DataFrame): DataFrame containing events with 'TIME' and optional `Exposure` columns.
        exposure_col (str): Name of the column in `event_df` containing the exposure correction factor.

    Returns:
        (pd.DataFrame): Updated BBA DataFrame with 'Counts', 'Length (s)', and 'Rate (cts/s)' columns.
    """
    counts = []
    lengths = []
    rates = []

    # block_frame=bba_df.copy()

    # initialize tuple to store info
    block_list = [
        {
            "start": None,
            "stop": None,
            "exposure": 0,
            "NuSTAR duration": 0,
            "counts": 0,
            "total exposure": 0,
            "rate": 0,
            "upperlim": 0,
            "lowerlim": 0,
        }
        for i in range(int(bba_df["block_label"].iloc[-1]) + 1)
    ]

    for _, block in bba_df.iterrows():
        block_start = block["start"]  # old version
        block_stop = block["stop"]  # old version
        block_label = int(block["block_label"])

        # start times
        if block_list[block_label]["start"] is None:
            block_list[block_label]["start"] = block_start

        # stop times
        block_list[block_label]["stop"] = block_stop

        # Calculate block length
        block_length = block_stop - block_start
        lengths.append(block_length)  # old version
        block_list[block_label]["exposure"] += block_length

        # Filter events in the current block
        block_events = event_df[
            (event_df["TIME"] >= block_start) & (event_df["TIME"] < block_stop)
        ]

        # Calculate number of events
        count = len(block_events)
        counts.append(count)
        block_list[block_label]["counts"] += count

        # get correction fraction (PSF, Vignetting, ect)
        correction_factor = np.mean(block_events["CORRECTION_FACTOR"])

        # Calculate event rate if doing exposure correction
        if not block_events.empty:
            total_exposure = block_events[exposure_col].sum()
            block_list[block_label]["total exposure"] += total_exposure
            rate = count / block_length
            rate /= correction_factor
        else:
            rate = 0.0
        rates.append(rate)

        # get the count rates of the chunks
        block_list[block_label]["rate"] = (
            block_list[block_label]["counts"] / block_list[block_label]["exposure"]
        ) / correction_factor

        # get livetime duration
        block_list[block_label]["NuSTAR duration"] = (
            block_list[block_label]["stop"] - block_list[block_label]["start"]
        )

        # also compute upper and lower count rate limits (for the block version I am working with)
        # counts = block_list[block_label]['counts']
        # length = block_list[block_label]['duration']

        # Gehrels approximation
        if block_list[block_label]["counts"] > 0:
            upper_limit = (
                block_list[block_label]["counts"]
                + np.sqrt(block_list[block_label]["counts"] + 0.75)
            ) / block_list[block_label]["exposure"]
            lower_limit = (
                (
                    block_list[block_label]["counts"]
                    - np.sqrt(block_list[block_label]["counts"] - 0.25)
                )
                / block_list[block_label]["exposure"]
                if block_list[block_label]["counts"] > 1
                else 0.0
            )

        else:
            upper_limit = 0.0
            lower_limit = 0.0

        block_list[block_label]["upperlim"] = upper_limit / correction_factor
        block_list[block_label]["lowerlim"] = lower_limit / correction_factor

        block.start_list = [
            block_list[i]["start"] for i in range(int(bba_df["block_label"].iloc[-1]))
        ]
        # upper_limits.append(upper_limit)
        # lower_limits.append(lower_limit)

    # Add new columns to the BBA DataFrame
    bba_df["Counts"] = counts
    bba_df["Length (s)"] = lengths
    bba_df["Rate (cts/s)"] = rates

    # construct new flare df
    flare_df = pd.DataFrame(block_list)

    return bba_df, flare_df


def calculate_confidence_limits(bba_df):
    """
    Calculate the 1-sigma confidence limits (upper and lower) for each block's event rate.

    Parameters:
        bba_df (pd.DataFrame): DataFrame containing BBA results with 'Counts' and 'Length (s)' columns.

    Returns:
        (pd.DataFrame): Updated BBA DataFrame with '1sig upper lim' and '1sig lower lim' columns.
    """
    upper_limits = []
    lower_limits = []

    for _, row in bba_df.iterrows():
        counts = row["Counts"]
        length = row["Length (s)"]

        # Gehrels approximation
        if counts > 0:
            upper_limit = (counts + np.sqrt(counts + 0.75)) / length
            lower_limit = (
                (counts - np.sqrt(counts - 0.25)) / length if counts > 1 else 0.0
            )
        else:
            upper_limit = 0.0
            lower_limit = 0.0

        upper_limits.append(upper_limit)
        lower_limits.append(lower_limit)

    # Add the calculated columns to the DataFrame
    bba_df["1sig upper lim"] = upper_limits
    bba_df["1sig lower lim"] = lower_limits

    return bba_df


def save_bba_results_txt(bba_df, analysis_params, output_file, save_results_note=None):
    """
    Save the Bayesian Block Analysis results to a text file with clean formatting.

    Parameters:
        bba_df (pd.DataFrame): DataFrame containing BBA results with columns like 'start', 'stop', 'Counts', 'length', 'rate', etc.
        analysis_params (dict): Dictionary containing analysis parameters (e.g., ncp_prior, fp_rate).
        output_file (str): Path to save the results text file.
        save_results_note (list): Optional list of strings to populate the "Notes" column.

    """
    # Add UT start and stop columns to BBA DataFrame
    bba_df["UT_start"] = bba_df["start"].apply(convert_nustar_to_utc).astype(str)
    bba_df["UT_stop"] = bba_df["stop"].apply(convert_nustar_to_utc).astype(str)

    # Add Notes column or populate it with save_results_note
    if save_results_note is not None:
        if len(save_results_note) != len(bba_df):
            raise ValueError(
                "Length of save_results_note must match the number of rows in bba_df."
            )
        bba_df["Notes"] = save_results_note
    elif "Notes" not in bba_df.columns:
        bba_df["Notes"] = ""  # Default to empty string if Notes not present

    # Create rows for the table
    table_data = []
    for _, row in bba_df.iterrows():
        table_data.append(
            [
                row["UT_start"].split(".")[0],  # Format UT start (precision: seconds)
                row["UT_stop"].split(".")[0],  # Format UT stop (precision: seconds)
                f"{row['start']:.1f}",
                f"{row['stop']:.1f}",
                # f"{int(row['start'])}",
                # f"{int(row['stop'])}",
                f"{int(row['Counts'])}",
                f"{row['Length (s)']:.1f}",
                f"{row['Rate (cts/s)']:.3E}",
                f"{row['1sig upper lim']:.3E}",
                f"{row['1sig lower lim']:.3E}",
                row["Notes"],
            ]
        )

    # Table headers
    headers = [
        "UT start",
        "UT stop",
        "NuSTAR start",
        "NuSTAR stop",
        "Counts",
        "Length (s)",
        "Rate (cts/s)",
        "1sig upper lim",
        "1sig lower lim",
        "Notes",
    ]

    # Center align columns
    align = "center"

    # Generate formatted table
    formatted_table = tabulate(
        table_data,
        headers=headers,
        tablefmt="plain",
        numalign=align,
        stralign=align,
        floatfmt=(".1f"),
    )

    with open(output_file, "w") as f:
        # Write analysis parameters
        f.write(f"Number of blocks: {len(bba_df)}\n")
        f.write(f"ncp_prior: {analysis_params['ncp_prior']}\n")
        f.write(f"fp_rate: {analysis_params['fp_rate']:.1E}\n")
        f.write(f"do_iter: {analysis_params['do_iter']}\n\n")

        # Write formatted table
        f.write(formatted_table)

    print("    Can customize notes in the bba_results.txt file.")


def save_flare_results_txt(
    bba_df, analysis_params, output_file, save_results_note=None
):
    """
    Save the flare Bayesian Block Analysis results to a text file with clean formatting.

    Parameters:
        bba_df (pd.DataFrame): DataFrame containing BBA results with columns like 'start', 'stop', 'Counts', 'length', 'rate', etc.
        analysis_params (dict): Dictionary containing analysis parameters (e.g., ncp_prior, fp_rate).
        output_file (str): Path to save the results text file.
        save_results_note (list): Optional list of strings to populate the "Notes" column.
    """
    # drop the total exposure column
    bba_df.drop(["total exposure"], axis=1)
    # Add UT start and stop columns to BBA DataFrame
    bba_df["UT_start"] = bba_df["start"].apply(convert_nustar_to_utc).astype(str)
    bba_df["UT_stop"] = bba_df["stop"].apply(convert_nustar_to_utc).astype(str)

    # Add Notes column or populate it with save_results_note
    if save_results_note is not None:
        if len(save_results_note) != len(bba_df):
            raise ValueError(
                "Length of save_results_note must match the number of rows in bba_df."
            )
        bba_df["Notes"] = save_results_note
    elif "Notes" not in bba_df.columns:
        bba_df["Notes"] = ""  # Default to empty string if Notes not present

    # Create rows for the table
    table_data = []
    for _, row in bba_df.iterrows():
        table_data.append(
            [
                row["UT_start"].split(".")[0],  # Format UT start (precision: seconds)
                row["UT_stop"].split(".")[0],  # Format UT stop (precision: seconds)
                f"{row['start']:.1f}",
                f"{row['stop']:.1f}",
                # f"{int(row['start'])}",
                # f"{int(row['stop'])}",
                f"{row['exposure']:.1f}",
                f"{row['NuSTAR duration']:.1f}",
                f"{int(row['counts'])}",
                f"{row['rate']:.3E}",
                f"{row['upperlim']:.3E}",
                f"{row['lowerlim']:.3E}",
                row["Notes"],
            ]
        )

    # Table headers
    headers = [
        "UT start",
        "UT stop",
        "NuSTAR start",
        "NuSTAR stop",
        "Length (s)",
        "Livetime (s)",
        "Counts",
        "Rate (cts/s)",
        "1sig upper lim",
        "1sig lower lim",
        "Notes",
    ]

    # Center align columns
    align = "center"

    # Generate formatted table
    formatted_table = tabulate(
        table_data,
        headers=headers,
        tablefmt="plain",
        numalign=align,
        stralign=align,
        floatfmt=(
            ".1f",
            ".1f",
            ".1f",
            ".1f",
            ".1f",
            ".1f",
            ".1f",
            ".4f",
            ".4f",
            ".4f",
            ".1f",
        ),
    )

    with open(output_file, "w") as f:
        # Write analysis parameters
        f.write(f"Number of blocks: {len(bba_df)}\n")
        f.write(f"ncp_prior: {analysis_params['ncp_prior']}\n")
        f.write(f"fp_rate: {analysis_params['fp_rate']:.1E}\n")
        f.write(f"do_iter: {analysis_params['do_iter']}\n\n")

        # Write formatted table
        f.write(formatted_table)

    # print(f"    Can customize notes in the bba_results.txt file.")
