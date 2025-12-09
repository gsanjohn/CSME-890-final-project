# 1. make like 20 or so event times (use random number generators)
# 2. make list of 20 or so event energies (use random number generators)
# 3. change header keywords
# 4. change GTIs (like a few minutes of good times) -- same for both files


# 5. make barycorr file


from astropy.table import Table
from astropy.io import fits
# /Users/gracesanger-johnson/Documents/grad school/Fall 2025/Computation/CSME-890-final-project/scripts/data_loader.py

# Duplicated files now need to edit to make it a dummy dataset


# def ModifyEvents(file, HeaderArray):
# fits_file = fits.open(file)[1]  # 1 for event data
# table = Table(fits_file)
# table.keep_columns(HeaderArray)
# table.write(file + "_test", format="fits", overwrite=True)


# ModifyEvents("testing/event_file_A", ["TIME", "PI"])

###############

from astropy.io import fits
import numpy as np

rng = np.random.default_rng()


def CreateBC(eventfile):
    # 1. Create data for the primary HDU
    primary_hdu = fits.PrimaryHDU()
    primary_hdu.header["TSTART"] = "10"
    # 4. Create an HDUList and append the HDUs
    hdul = fits.HDUList([primary_hdu])

    # 5. Write the HDUList to a FITS file
    filename = eventfile
    hdul.writeto(filename, overwrite=True)

    print(f"Multi-extension FITS file '{filename}' created successfully.")


CreateBC("./testing/test_eventsBC.fits")
