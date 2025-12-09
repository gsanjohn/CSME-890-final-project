import pytest
from astropy.io import fits

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
