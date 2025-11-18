# get_event_corr_factor.py
import numpy as np
import pandas as pd
from astropy.io import fits


def get_event_corr_factor(filename, tt):
    """
    Calculate correction factors for photon events based on the correction factor FITS file.

    Parameters:
        filename (str): Path to the correction factor FITS file.
        tt (array-like): Photon arrival times (e.g., DataFrame['TIME']).

    Returns:
        (np.ndarray): Array of correction factors corresponding to photon times.
    """
    try:
        # Read the correction factor FITS file
        with fits.open(filename) as hdul:
            data = hdul[1].data

            # Convert relevant columns to a DataFrame with proper byte order
            df_corr = pd.DataFrame(
                {
                    "TSTART": data["TSTART"].byteswap().newbyteorder(),
                    "TSTOP": data["TSTOP"].byteswap().newbyteorder(),
                    "FRACTION": data["FRACTION"].byteswap().newbyteorder(),
                }
            )

        # Debugging: Print the first few rows of the correction factor DataFrame
        # print(f"Correction Factor DataFrame (Converted):\n{df_corr.head()}")

        # Initialize the correction factor array with a default value (e.g., 1.0)
        correction_factors = np.full(len(tt), 1.0)

        # Match photon times to correction intervals
        for i, t in enumerate(tt):
            match = df_corr[(df_corr["TSTART"] <= t) & (df_corr["TSTOP"] > t)]
            if not match.empty:
                correction_factors[i] = match.iloc[0]["FRACTION"]
            else:
                # Log warning if no match is found
                # pass
                print(
                    f"Warning: No match found for photon time {t}. Using default correction factor 1.0."
                )

        return correction_factors

    except Exception as e:
        raise RuntimeError(f"Error in GetEventcorrfactor: {e}")
