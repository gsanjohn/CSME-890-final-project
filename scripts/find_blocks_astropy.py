# from astropy.stats import bayesian_blocks
from scripts.expo_events import bayesian_blocks
import pandas as pd
import numpy as np
from scripts.find_blocks import upper_limit_gehrels, lower_limit_gehrels


def bba_astropy(time, ncp_prior, fp_rate=0.05, x_list=None):
    """
    Perform Bayesian Block segmentation using Astropy.

    Parameters:
        time (np.ndarray): Photon arrival times.
        ncp_prior (float): (optional) Number of change point prior
        fp_rate (float): False positive rate for change points.

    Returns:
        (pd.DataFrame): DataFrame with Bayesian Block intervals and statistics.
    """
    # make x values (correct for exposure)
    exp_list = x_list
    # change_points = bayesian_blocks(time, fitness="events", p0=fp_rate)
    change_points = bayesian_blocks(
        time, ex=exp_list, fitness="events", ncp_prior=ncp_prior
    )  # exposure version
    # change_points = bayesian_blocks(time, fitness="events", ncp_prior=ncp_prior)
    # Tuning Hyperparameter p0:
    # Perform Bayesian Blocks segmentation with different p0 values
    # for p0 in [0.5, 0.1, 0.05, 0.01, 0.005]:
    #     change_points = bayesian_blocks(time, fitness="events", p0=p0)
    #     print(f"    p0 = {p0}, Number of change points: {len(change_points) - 1}")
    print(f"    Change points: {change_points}")

    # Calculate intervals, durations, and rates
    block_start = change_points[:-1]
    block_stop = change_points[1:]
    durations = block_stop - block_start
    counts = np.histogram(time, bins=change_points)[0]
    rates = counts / durations
    block_label = np.array(range(0, len(durations)))

    # Format results
    df = pd.DataFrame(
        {
            "start": block_start,
            "stop": block_stop,
            "duration": durations,
            "counts": counts,
            "rate": rates,
            "upperlim": upper_limit_gehrels(0.8413, counts) / durations,
            "lowerlim": lower_limit_gehrels(0.8413, counts) / durations,
            "fp_rate": fp_rate,
            "block_label": block_label,
        }
    )
    return df
