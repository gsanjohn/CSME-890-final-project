import pytest
from astropy.io import fits
import pandas as pd

from scripts.data_loader import (
    duplicate_fits,
    load_gti_file,
    load_event_file,
    empty_df,
)
from scripts.barycenter_corr import barycorr
from scripts.event_filter import filter_events_by_energy
from scripts.clean_gti import clean_gti
from scripts.merge_gtis import merge_gtis
from scripts.event_common_gti import filter_events_with_common_gti
from scripts.get_event_corr_factor import get_event_corr_factor
from scripts.merge_events import merge_events
from scripts.calculate_average_rate import calculate_average_rate
from scripts.suppress_gti_gaps import suppress_gti_gaps
from scripts.find_blocks import find_blocks, format_bayesian_block_output
from scripts.find_blocks_astropy import bba_astropy
from scripts.detailed_flare_analysis import detailed_flare_analysis
from scripts.insert_gaps import insert_gti_gaps
from scripts.save_bba_results import (
    convert_nustar_to_utc,
    save_bba_results_txt,
    save_flare_results_txt,
    calculate_event_counts_and_rates,
    calculate_confidence_limits,
)
from scripts.create_lightcurve import generate_lightcurve, plot_prep
from scripts.plot_lc import plot_lightcurve
from scripts.make_readme import write_readme


def fits_diff(file1, file2):
    diff = fits.FITSDiff(file1, file2)
    return diff


@pytest.mark.parametrize(
    "event_file, lccorrfile, output_dir",
    [
        (
            "./tests/data/test_eventsA.fits",
            "./tests/data/test_LCcorrA.fits",
            "./tests/test_output/",
        ),
        (
            "./tests/data/test_eventsB.fits",
            "./tests/data/test_LCcorrB.fits",
            "./tests/test_output/",
        ),
    ],
)
def test_duplicate_fits(event_file, lccorrfile, output_dir):
    duplicate_fits(event_file, lccorrfile, output_dir)
    diff1 = fits_diff(event_file, "./tests/test_output/event_file_C")
    assert diff1.identical, diff1.report()
    diff2 = fits_diff(lccorrfile, "./tests/test_output/lccorrfile_C")
    assert diff2.identical, diff2.report()


@pytest.mark.parametrize("file_path", [("./tests/data/test_eventsA.fits")])
def test_load_event_fileA(file_path):
    output = load_event_file(file_path)

    expected_df_A = pd.DataFrame(
        {
            "TIME": [
                9.0,
                32.0,
                130.0,
                166.0,
                174.0,
                182.0,
                221.0,
                265.0,
                278.0,
                296.0,
                313.0,
                326.0,
                334.0,
                369.0,
                378.0,
                422.0,
                571.0,
                586.0,
                595.0,
                596.0,
            ],
            "PI": [
                1226.0,
                1574.0,
                1347.0,
                768.0,
                1026.0,
                34.0,
                1080.0,
                829.0,
                1479.0,
                599.0,
                584.0,
                737.0,
                135.0,
                1380.0,
                836.0,
                1511.0,
                141.0,
                1316.0,
                844.0,
                1064.0,
            ],
        }
    )
    # Calculate Energy
    expected_df_A["Energy"] = expected_df_A["PI"] * 0.04 + 1.6

    pd.testing.assert_frame_equal(output, expected_df_A)
