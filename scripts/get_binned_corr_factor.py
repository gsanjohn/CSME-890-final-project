import numpy as np


#### Need to do it twice for A and B for IDL
def get_binned_corr_factor(merged_events, bin_times):
    """
    Calculate correction factors for photon events based on the correction factor FITS file.

    Parameters:
        merged_events (pd.DataFrame): data frame containing photon arrival times and corrections ("exposure").
        bin_times (array-like): start and stop times of data bins.

    Returns:
        (np.ndarray): Array same length as bin_times.
    """

    fraction = merged_events["CORRECTION_FACTOR"]
    tstart = merged_events["TIME"]
    corr_factors = np.zeros(len(bin_times))

    for i in range(len(bin_times)):
        ilo = np.argmin(np.abs(bin_times[i][0] - tstart))
        ihi = np.argmin(np.abs(bin_times[i][1] - tstart))
        corr_factors[i] = np.mean(fraction[ilo:ihi])

    return corr_factors
