import numpy as np
import pandas as pd


def find_blocks(time, exposure, fp_rate=0.05, ncp_prior=None, do_iter=False):
    """
    Identify change points for Bayesian Blocks segmentation.

    Parameters:
        time (np.ndarray): Photon arrival times.
        exposure (np.ndarray): Exposure correction factors.
        fp_rate (float): False positive rate for change points.
        ncp_prior (float, optional): Prior number of change points. Defaults to None.
        do_iter (bool): Iterative refinement flag.

    Returns:
        (dict): Bayesian block structure with change points, counts, rates, etc.
    """
    # Step 1: Sort time and exposure
    sorted_indices = np.argsort(time)
    time = time[sorted_indices]
    exposure = exposure[sorted_indices]

    # Step 2: Calculate durations and event counts per block
    block_start = [time[0]]
    block_stop = []
    event_counts = []
    cumulative_exposure = 0
    count = 0

    for i, t in enumerate(time):
        count += 1
        cumulative_exposure += 1.0 / exposure[i]
        if ncp_prior and cumulative_exposure > ncp_prior:  # Trigger segmentation
            block_stop.append(t)
            block_start.append(t)
            event_counts.append(count)
            count = 0
            cumulative_exposure = 0

    block_stop.append(time[-1])
    event_counts.append(count)

    # Step 3: Calculate rates
    block_durations = np.array(block_stop) - np.array(block_start)
    rates = np.array(event_counts) / block_durations

    # print(ncp_prior)

    return {
        "change_points": block_start,
        "start": block_start,
        "stop": block_stop,
        "duration": block_durations,
        "counts": event_counts,
        "rate": rates,
    }


def upper_limit_gehrels(confidence_level, counts):
    """
    Compute the upper confidence limit for count rates.

    Parameters:
        confidence_level (float): Confidence level (e.g., 0.8413 for 1-sigma).
        counts (array-like): Event counts.

    Returns:
        (np.ndarray): Upper confidence limits.
    """
    counts = np.asarray(counts)  # Ensure counts is a NumPy array
    return counts + np.sqrt(counts + 0.75)


def lower_limit_gehrels(confidence_level, counts):
    """
    Compute the lower confidence limit for count rates.

    Parameters:
        confidence_level (float): Confidence level (e.g., 0.8413 for 1-sigma).
        counts (array-like): Event counts.

    Returns:
        (np.ndarray): Lower confidence limits.
    """
    counts = np.asarray(counts)  # Ensure counts is a NumPy array
    return np.maximum(counts - np.sqrt(counts - 0.25), 0)  # Avoid negative values


def format_bayesian_block_output(results, fp_rate, ncp_prior, do_iter):
    """
    Format Bayesian Block results into a DataFrame.

    Parameters:
        results (pd.DataFrame):
        fp_rate(float): false positive rate used in bayesian block analysis
        ncp_prior(float): number of change point prior used in bayesian block analysis
        do_iter(boolean): iterate over blocks

    Returns:
        (pd.DataFrame): Formated block information with upper and lower count rate limits
    """
    # print('right before printing', ncp_prior)
    df = pd.DataFrame(
        {
            "change_points": results["change_points"],
            "start": results["start"],
            "stop": results["stop"],
            "duration": results["duration"],
            "counts": results["counts"],
            "rate": results["rate"],
            "upperlim": upper_limit_gehrels(0.8413, results["counts"])
            / results["duration"],
            "lowerlim": lower_limit_gehrels(0.8413, results["counts"])
            / results["duration"],
            "ncp_prior": ncp_prior,
            "fp_rate": fp_rate,
            "do_iter": do_iter,
        }
    )
    return df
