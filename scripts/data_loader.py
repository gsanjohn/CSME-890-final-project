import pandas as pd
import numpy as np
from astropy.io import fits


def duplicate_fits(event_file: str, lccorrfile: str, output_dir: str):
    """
    Takes an event file and the corresponding light curve correction file and makes a duplicate of them in the output directory.

    Parameters:
        event_file (str): FITS file to be duplicated
        lccorrfile (str): light curve correction file to be duplicated
        output_dir (str): output directory for duplicated files

    """
    input_file = [event_file, lccorrfile]
    output_file = ["event_file_C", "lccorrfile_C"]

    for i in range(len(input_file)):
        # Open original FITS
        with fits.open(input_file[i]) as hdul:
            new_hdus = []

            for hdu in hdul:
                if isinstance(hdu, fits.PrimaryHDU):
                    # Primary HDU (usually image data) → keep header, empty data
                    new_hdu = fits.PrimaryHDU(data=None, header=hdu.header)

                elif isinstance(hdu, fits.ImageHDU):
                    # Image extension → keep header, drop data
                    new_hdu = fits.ImageHDU(data=None, header=hdu.header)

                elif isinstance(hdu, fits.BinTableHDU) or isinstance(
                    hdu, fits.TableHDU
                ):
                    # Table extensions → preserve column structure but no rows
                    cols = hdu.columns
                    new_hdu = hdu.__class__.from_columns(cols, nrows=0)
                    # Copy extra header keywords
                    for key, val in hdu.header.items():
                        if key not in new_hdu.header:
                            new_hdu.header[key] = val

                else:
                    # For any other HDU type, just copy header and strip data
                    new_hdu = hdu.copy()
                    new_hdu.data = None

                new_hdus.append(new_hdu)

            # Write new FITS file with same structure, empty data
            hdul_new = fits.HDUList(new_hdus)
            hdul_new.writeto(output_dir + output_file[i], overwrite=True)

    print(
        f"Copy of FITS with same structure (all extensions) saved to {output_dir, output_file}"
    )


def empty_df(array):
    df_copy = pd.DataFrame(np.nan, index=array.index, columns=array.columns)
    return df_copy


def validate_gti_columns(columns: list):
    """
    Validate that the GTI file contains the required columns.

    Parameters:
        columns (list): List of column names.

    Raises:
        ValueError: If required columns are missing.
    """
    required_columns = {"START", "STOP"}
    if not required_columns.issubset(set(columns)):
        raise ValueError(
            f"GTI file is missing required columns: {required_columns - set(columns)}"
        )


def load_gti_file(file_path):
    """
    Load a GTI file using Astropy and convert it to a DataFrame.

    Parameters:
        file_path (str): Path to the GTI FITS file.

    Returns:
        (pd.DataFrame): GTI DataFrame with 'START' and 'STOP' columns.
    """
    try:
        with fits.open(file_path) as hdul:
            gti_data = hdul[2].data
            validate_gti_columns(gti_data.names)
            return pd.DataFrame({"START": gti_data["START"], "STOP": gti_data["STOP"]})
    except Exception as e:
        raise RuntimeError(f"Failed to load GTI from file {file_path}: {e}")


def load_event_file(file_path):
    """
    Load an event file using Astropy and convert it to a DataFrame.

    Parameters:
        file_path (str): Path to the event FITS file.

    Returns:
        (pd.DataFrame): Event DataFrame with 'TIME', 'PI', and 'Energy' columns.
    """
    try:
        with fits.open(file_path) as hdul:
            event_data = hdul[1].data

            # Convert columns to native byte order
            time_data = event_data["TIME"].byteswap().newbyteorder()
            pi_data = event_data["PI"].byteswap().newbyteorder()

            # Create DataFrame
            df = pd.DataFrame({"TIME": time_data, "PI": pi_data})

            # Calculate Energy
            df["Energy"] = df["PI"] * 0.04 + 1.6
            return df
    except Exception as e:
        raise RuntimeError(f"Failed to load event file {file_path}: {e}")
