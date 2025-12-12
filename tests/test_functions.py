import pytest
from astropy.io import fits
import pandas as pd
from pandas.testing import assert_frame_equal, assert_series_equal
import numpy as np
from datetime import datetime


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
    "corr_file, event_file_a, event_file_b, time",
    [
        (
            "./tests/data/test_eventsBC.fits",
            "./tests/data/test_eventsA.fits",
            "./tests/data/test_eventsA.fits",
            10,
        )
    ],
)
def test_barycorr(corr_file, event_file_a, event_file_b, time):
    gtiA = load_gti_file(event_file_a)
    gtiB = load_gti_file(event_file_b)
    eventsA = load_event_file(event_file_a)
    eventsB = load_event_file(event_file_b)

    out_eventsA, out_eventsB, out_gtiA, out_gtiB = barycorr(
        corr_file, event_file_a, eventsA, eventsB, gtiA, gtiB
    )

    assert (out_eventsA["TIME"] == eventsA["TIME"]).all()
    assert (out_eventsB["TIME"] == eventsB["TIME"]).all()
    assert out_gtiA.equals(gtiA_df + 10)
    assert out_gtiB.equals(gtiB_df + 10)


@pytest.mark.parametrize("data, energy_min, energy_max", [(eventsA_df, 0, 20)])
def test_filter_events_by_energy(data, energy_min, energy_max):
    output = filter_events_by_energy(data, energy_min, energy_max)
    expected_energy = [2.96, 7.00, 7.24]
    assert (output["Energy"].values == expected_energy).all()


@pytest.mark.parametrize(
    "eventsfile, threshold, starttrim, stoptrim",
    [("./tests/data/test_eventsA.fits", 60, 1, 1)],
)
def test_clean_gti(eventsfile, threshold, starttrim, stoptrim):
    gti = load_gti_file(eventsfile)
    events = load_event_file(eventsfile)
    output = clean_gti(gti, events, threshold, starttrim, stoptrim)
    expected_output = pd.DataFrame(
        {"START": 101.0, "STOP": 199.0, "DURATION": 98.0, "EVENT_COUNT": 4}, index=[1]
    )
    assert output.equals(expected_output)


@pytest.mark.parametrize(
    "eventsfileA,eventsfileB",
    [("./tests/data/test_eventsA.fits", "./tests/data/test_eventsB.fits")],
)
def test_merge_gtis(eventsfileA, eventsfileB):
    gtiA = load_gti_file(eventsfileA)
    gtiB = load_gti_file(eventsfileB)
    output = merge_gtis(gtiA, gtiB)

    assert (output["START"] == gtiA["START"]).all()
    assert (output["STOP"] == gtiA["STOP"]).all()


@pytest.mark.parametrize("data", [("./tests/data/test_eventsA.fits")])
def test_filter_events_with_common_gti(data):
    events = load_event_file(data)
    gti = load_gti_file(data)
    output = filter_events_with_common_gti(events, gti)
    expected = pd.DataFrame({"TIME": [9.0, 130.0, 166.0, 174.0, 182.0]})
    assert (output["TIME"] == expected["TIME"]).all()


@pytest.mark.parametrize(
    "lccorr,data",
    [("./tests/data/test_LCcorrA.fits", "./tests/data/test_eventsA.fits")],
)
def test_get_event_corr_factor(lccorr, data):
    events = load_event_file(data)
    output = get_event_corr_factor(lccorr, events["TIME"])
    expected = np.array(
        [
            0.71531577,
            0.71531577,
            0.58791119,
            0.58791119,
            0.58791119,
            0.58791119,
            0.54996284,
            0.54996284,
            0.54996284,
            0.54996284,
            0.40729103,
            0.40729103,
            0.40729103,
            0.40729103,
            0.40729103,
            0.45672174,
            0.93246583,
            0.93246583,
            0.93246583,
            0.93246583,
        ]
    )
    # print(all([a == b for a, b in zip(output, expected)]))
    # assert len(output) == len(expected)
    # assert all([a == b for a, b in zip(output, expected)])

    assert output == pytest.approx(expected)


@pytest.mark.parametrize(
    "eventA, eventB,lcA,lcB",
    [
        (
            "./tests/data/test_eventsA.fits",
            "./tests/data/test_eventsB.fits",
            "./tests/data/test_LCcorrA.fits",
            "./tests/data/test_LCcorrB.fits",
        )
    ],
)
def test_merge_events(eventA, eventB, lcA, lcB):
    eventsA = load_event_file(eventA)[0:5]
    eventsB = load_event_file(eventB)[0:5]
    exposureA = get_event_corr_factor(lcA, eventsA["TIME"])
    exposureB = get_event_corr_factor(lcB, eventsB["TIME"])

    # print(eventsA)
    # print(eventsB)
    # print(exposureA)
    # print(eventsA)

    output = merge_events(eventsA, eventsB, exposureA, exposureB)
    expected = pd.DataFrame(
        {
            "TIME": [5.0, 9.0, 32.0, 67.0, 81.0, 113.0, 115.0, 130.0, 166.0, 174.0],
            "PI": [
                971.0,
                1226.0,
                1574.0,
                1728.0,
                173.0,
                1875.0,
                1228.0,
                1347.0,
                768.0,
                1026.0,
            ],
            "Energy": [
                40.44,
                50.64,
                64.56,
                70.72,
                8.52,
                76.60,
                50.72,
                55.48,
                32.32,
                42.64,
            ],
            "Exposure": [
                0.478264,
                0.715316,
                0.715316,
                0.478264,
                0.478264,
                0.103396,
                0.103396,
                0.587911,
                0.587911,
                0.587911,
            ],
            "Module": [
                "B",
                "A",
                "A",
                "B",
                "B",
                "B",
                "B",
                "A",
                "A",
                "A",
            ],
        }
    )
    print(output)
    # assert output.equals(expected)
    assert_frame_equal(output, expected, check_dtype=False)


@pytest.mark.parametrize(
    "lccorr,data",
    [("./tests/data/test_LCcorrA.fits", "./tests/data/test_eventsA.fits")],
)
def test_calculate_average_rate(lccorr, data):
    get_events = load_event_file(data)
    format_events = get_event_corr_factor(lccorr, get_events["TIME"])
    gti = load_gti_file(data)
    output = calculate_average_rate(
        format_events, gti, flare_gti=None, calculate_average_count_rate=False
    )
    assert output is None


@pytest.mark.parametrize(
    "data,original_tt_stop",
    [("./tests/data/test_eventsA.fits", 600)],
)
def test_suppress_gti_gaps(data, original_tt_stop):
    events = load_event_file(data)[0:5]
    gti = load_gti_file(data)
    updated_event_df, updated_tt_stop, cumulative_gap_times = suppress_gti_gaps(
        events, gti, original_tt_stop
    )
    expected_df = pd.DataFrame(
        {
            "TIME": [
                9.0,
                32.0,
                40.0,
                76.0,
                84.0,
            ],
            "PI": [
                1226.0,
                1574.0,
                1347.0,
                768.0,
                1026.0,
            ],
            "Energy": [50.64, 64.56, 55.48, 32.32, 42.64],
        }
    )
    expected_tt_stop = 84.0
    expected_gap_times = np.array([0.0, 90.0, 340.0, 344.0])

    assert_frame_equal(updated_event_df, expected_df, check_dtype=False)
    assert updated_tt_stop == expected_tt_stop
    assert (cumulative_gap_times == expected_gap_times).all


@pytest.mark.parametrize(
    "eventsfileA, eventsfileB, BBA_df",
    [
        (
            "./tests/data/test_eventsA.fits",
            "./tests/data/test_eventsB.fits",
            pd.DataFrame({"start": [0], "stop": [40]}),
        )
    ],
)
def test_insert_gti_gaps(eventsfileA, eventsfileB, BBA_df):
    gtiA = load_gti_file(eventsfileA)
    gtiB = load_gti_file(eventsfileB)
    merged_gtis = merge_gtis(gtiA, gtiB)
    output_corrected_blocks, output_gti_gaps = insert_gti_gaps(BBA_df, merged_gtis)

    expected_corrected_blocks = pd.DataFrame({"start": [0, 100], "stop": [10, 130]})
    expected_gti_gaps = pd.DataFrame(
        {
            "GAP_START": [10, 200, 500],
            "GAP_STOP": [100, 450, 504],
            "GAP_DURATION": [90, 250, 4],
        }
    )

    assert_frame_equal(
        output_corrected_blocks, expected_corrected_blocks, check_dtype=False
    )
    assert_frame_equal(output_gti_gaps, expected_gti_gaps, check_dtype=False)


@pytest.mark.parametrize(
    "nustar_time,utc_time",
    [(0, datetime(2010, 1, 1, 0, 0)), (50000, datetime(2010, 1, 1, 13, 53, 20))],
)
def test_convert_nustar_to_utc(nustar_time, utc_time):
    output = convert_nustar_to_utc(nustar_time)

    assert output == utc_time


@pytest.mark.parametrize(
    "data,lccorr,expected_counts,expected_rates",
    [("./tests/data/test_eventsA.fits", "./tests/data/test_LCcorrA.fits", 500, 5)],
)
def test_calculate_event_counts_and_rates(
    data, lccorr, expected_counts, expected_rates
):
    events = load_event_file(data)
    full_events = get_event_corr_factor(lccorr, events["TIME"])
    bba_frame = pd.DataFrame(
        {
            "start": [0],
            "stop": [100],
            "duration": [100],
            "counts": [500],
            "rate": [500],
            "block_label": [0],
        }
    )
    output = calculate_event_counts_and_rates(
        bba_frame, full_events, exposure_col="Exposure"
    )

    assert output == expected_counts
    assert output == expected_rates
    # assert (output["COUNTS"] == expected_counts).all()
    # assert (output["COUNTS"] == expected_rates).all()


@pytest.mark.parametrize(
    "data, expected_data",
    [
        (
            pd.DataFrame({"Counts": [10, 20, 30], "Length (s)": [100, 200, 300]}),
            pd.DataFrame(
                {
                    "1sig upper lim": [0.132787, 0.122776, 0.118484],
                    "1sig lower lim": [0.068775, 0.077780, 0.081819],
                }
            ),
        ),
        (
            pd.DataFrame({"Counts": [40, 50, 60], "Length (s)": [400, 500, 600]}),
            pd.DataFrame(
                {
                    "1sig upper lim": [0.115959, 0.114248, 0.112990],
                    "1sig lower lim": [0.084238, 0.085893, 0.087117],
                }
            ),
        ),
    ],
)
def test_calculate_confidence_limits(data, expected_data):
    output = calculate_confidence_limits(data)

    assert_series_equal(
        output["1sig upper lim"], expected_data["1sig upper lim"], check_dtype=False
    )

    assert_series_equal(
        output["1sig lower lim"], expected_data["1sig lower lim"], check_dtype=False
    )


@pytest.mark.parametrize(
    "data,lc_data,binsize,energy_range,gti_average, expected",
    [
        (
            "./tests/data/test_eventsA.fits",
            "./tests/data/test_LCcorrA.fits",
            100,
            (3, 79),
            True,
            pd.DataFrame(
                {
                    "bin_start": [0.0, 450.0],
                    "bin_end": [
                        100.0,
                        504.0,
                    ],
                    "count_rate": [0.02796, np.NaN],
                    "upper_limit": [0.051143, np.NaN],
                    "lower_limit": [0.009466, np.NaN],
                },
                index=[0, 2],
            ),
        ),
        (
            "./tests/data/test_eventsA.fits",
            "./tests/data/test_LCcorrA.fits",
            100,
            (0, 100),
            True,
            pd.DataFrame(
                {
                    "bin_start": [0.0, 450.0],
                    "bin_end": [100.0, 504.0],
                    "count_rate": [0.02796, np.NaN],
                    "upper_limit": [0.051143, np.NaN],
                    "lower_limit": [0.009466, np.NaN],
                },
                index=[0, 2],
            ),
        ),
    ],
)
def test_generate_lightcurve(
    data, lc_data, binsize, energy_range, gti_average, expected
):
    events = load_event_file(data)
    gti = load_gti_file(data)
    events["CORRECTION_FACTOR"] = get_event_corr_factor(lc_data, events["TIME"])
    output = generate_lightcurve(
        events_df=events,  # Use the events_merged DataFrame from Step 7
        gti_df=gti,  # Use the GTI DataFrame from Step 4
        binsize=binsize,
        energy_range=energy_range,
        gti_average=True,
    )  # Set to True if GTI-based binning

    output.reset_index(drop=True)
    print("expected", expected)
    print("output", output)
    assert_frame_equal(output, expected, check_dtype=False)
    # assert_series_equal(output["bin_end"], expected["bin_end"], check_dtype=False)
    # assert_series_equal(output["count_rate"], expected["count_rate"], check_dtype=False)
