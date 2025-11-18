from astropy.io import fits
import pandas as pd


def barycorr(
    corr_file: str,
    event_file_a: str,
    eventsA: pd.DataFrame,
    eventsB: pd.DataFrame,
    gtiA: pd.DataFrame,
    gtiB: pd.DataFrame,
):
    """
    Load barycenter corrected event file and xselected event file, get time correction, apply time correction.

    Args:
        corr_file (str): Path to the corrected event FITS file
        event_file_a (str): path to event FITS file
        eventsA (pd.DataFrame): events dataframe
        eventsB (pd.DataFrame): events dataframe
        gtiA (pd.DataFrame): gti dataframe
        gtiB (pd.DataFrame): gti dataframe

    Returns:
        (pd.DataFrame): Event DataFrame with 'TIME' corrected, gti DataFrame with 'GTI' corrected.

    """

    # Read TSTART from barycenter corrected event
    with fits.open(corr_file) as hdul:
        header = hdul[0].header
        tstartBC = header.get("TSTART")

    # Read TSTART from events file
    with fits.open(event_file_a) as hdul:
        header = hdul[0].header
        tstart = header.get("TSTART")

    # get delta t
    deltaT = tstartBC - tstart

    # add correction to events
    eventsA["TIME"] += deltaT
    eventsB["TIME"] += deltaT

    # add correction to GTIs
    gtiA += deltaT
    gtiB += deltaT

    return eventsA, eventsB, gtiA, gtiB
