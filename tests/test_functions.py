import pytest
from astropy.io import fits
import pandas as pd
import numpy as np

from scripts.data_loader import duplicate_fits, load_gti_file, load_event_file
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


eventsA_df = pd.DataFrame(
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

eventsA_df["Energy"] = eventsA_df["PI"] * 0.04 + 1.6

gtiA_df = pd.DataFrame(
    {"START": [0.0, 100.0, 450.0, 504.0], "STOP": [10.0, 200.0, 500.0, 506.0]}
)

eventsB_df = pd.DataFrame(
    {
        "TIME": [
            5.0,
            67.0,
            81.0,
            113.0,
            115.0,
            133.0,
            142.0,
            148.0,
            167.0,
            211.0,
            216.0,
            219.0,
            249.0,
            307.0,
            312.0,
            325.0,
            386.0,
            393.0,
            454.0,
            518.0,
        ],
        "PI": [
            971.0,
            1728.0,
            173.0,
            1875.0,
            1228.0,
            94.0,
            52.0,
            985.0,
            1190.0,
            105.0,
            401.0,
            287.0,
            222.0,
            730.0,
            190.0,
            871.0,
            257.0,
            34.0,
            1690.0,
            424.0,
        ],
    }
)
eventsB_df["Energy"] = eventsB_df["PI"] * 0.04 + 1.6
gtiB_df = pd.DataFrame(
    {"START": [0.0, 100.0, 450.0, 504.0], "STOP": [10.0, 200.0, 500.0, 506.0]}
)
lccorrA_df = pd.DataFrame(
    {
        "TSTART": [0.0, 100.0, 200.0, 300.0, 400.0, 500.0],
        "TSTOP": [100.0, 200.0, 300.0, 400.0, 500.0, 600.0],
        "FRACTION": [
            0.715315770313184,
            0.587911189664267,
            0.549962835941608,
            0.407291029188348,
            0.456721739757098,
            0.932465828637084,
        ],
    }
)
lccorrB_df = pd.DataFrame(
    {
        "TSTART": [0.0, 100.0, 200.0, 300.0, 400.0, 500.0],
        "TSTOP": [100.0, 200.0, 300.0, 400.0, 500.0, 600.0],
        "FRACTION": [
            0.478264015916526,
            0.103395804052955,
            0.612026275541879,
            0.588533709372726,
            0.845335600199995,
            0.610604541451639,
        ],
    }
)


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
def test_load_event_file(file_path):
    output = load_event_file(file_path)
    pd.testing.assert_frame_equal(output, eventsA_df)


@pytest.mark.parametrize("file_path", [("./tests/data/test_eventsA.fits")])
def test_load_gti_file(file_path):
    output = load_gti_file(file_path)

    pd.testing.assert_frame_equal(output, gtiA_df, check_dtype=False)


@pytest.mark.parametrize(
    "corr_file, event_file_a, eventsA, eventsB, gtiA, gtiB, time",
    [
        (
            "./tests/data/test_eventsBC.fits",
            "./tests/data/test_eventsA.fits",
            eventsA_df,
            eventsB_df,
            gtiA_df,
            gtiB_df,
            10,
        )
    ],
)
def test_barycorr(corr_file, event_file_a, eventsA, eventsB, gtiA, gtiB, time):
    out_eventsA, out_eventsB, out_gtiA, out_gtiB = barycorr(
        corr_file, event_file_a, eventsA, eventsB, gtiA, gtiB
    )
    assert (out_eventsA["TIME"] == eventsA_df["TIME"]).all()
    assert (out_eventsB["TIME"] == eventsB_df["TIME"]).all()
    assert out_gtiA.equals(gtiA_df)
    assert out_gtiB.equals(gtiB_df)


@pytest.mark.parametrize("data, energy_min, energy_max", [(eventsA_df, 0, 20)])
def test_filter_events_by_energy(data, energy_min, energy_max):
    output = filter_events_by_energy(data, energy_min, energy_max)
    expected_energy = [2.96, 7.00, 7.24]
    assert (output["Energy"].values == expected_energy).all()
