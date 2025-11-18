import pandas as pd
import numpy as np
import os

from scripts.data_loader import duplicate_fits, load_gti_file, load_event_file, empty_df
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


def main():
    """Bayesian Block alorithm for use with NuSTAR data.

    3 pairs of files are necessary as input:
        - lightcurve correction factor file (produced by nuproducts)
        - event list extracted from region of interest (produced by xselect)
        - barycorr file to do barycenter correction

    Has capability to handle only one module. To do that, simply comment out path names for the missing module



    Flow of the code:
        1. define filenames and paths
        2. load events and filter by energy
        3. Load and clean GTIs
        4. Merge GTIs
        5. Filter events outside of common GTIs
        6. Get exposure correction factor
        7. Merge events from FPM A and B
        8. Get average count rates in GTIs
        9. Remove time gaps from events (in preperation for BB analysis)
        10. Perform Bayesian Block Routine
        11. Re-insert GTI gaps into data
        12. Print out txt file of BBA results
        13. Bin events into light curve
        14. Make light curve plot

    NOTE: parameters do not appear in the function signature, but instead are editable parameters of the code.

    Parameters:
        year (str): year of the observation being analyzed (used for leap second correction)
        obsID (str): observation ID of observation being analyzed
        output_dir (str): path to output directory
        event_file_a (str): path to xselected event file for module A (optional)
        event_file_b (str): path to xselected event file for module B (optional)
        barycorr_event (str): path to barycenter corrected event file for module A or B (optional)
        lccorrfileA (str): path to light curve correction file for module A (optional)
        lccorrfileB (str): path to light curve correction file for module B (optional)
        energy_min (int): energy range in keV for bayesian block analysis. Default 3.0
        energy_max (int): energy range in keV for bayesian block analysis. Default 30.0
        fp_rate (int): false positive rate used for bayesian block detection. Default 0.01
        binsize (int): bin size in seconds for light curve plotting. Default 100
        energy_range (list): energy range in keV for plotted light curve. Default: (3.0, 79.0)




    Returns:

        12_bba_results.txt: text file with times, count rates, upper and lower limits
        14_lightcurve_plot.png: Light curve plot with 1 sigma upper and lower count rates per block
        ReadMe.txt: text file inputs and various script settings





    Grace Sanger-Johnson 9/23/2025 based on runBB.pro IDL package created by Nicolas Barriere
    """
    year = "2016"
    obsID = "40202001002"
    path = (
        "/Users/gracesanger-johnson/xray_astro/SgrA/"
        + year
        + "/"
        + obsID
        + "/event_cl2/BBA_results"
    )

    # Define output directory
    output_dir = path + "/3-30_BB_output/"

    # make sure directory exists
    os.makedirs(output_dir, exist_ok=True)
    # output_dir = path+"/testing/fp_testing/"

    # Set default file paths (DO NOT EDIT):
    event_file_a = None
    event_file_b = None

    # Define input file paths
    barycorr_event = path + "/nu" + obsID + "A01_cl_barycorr.evt"
    event_file_a = path + "/nu" + obsID + "A01_xselected.evt"
    event_file_b = path + "/nu" + obsID + "B01_xselected.evt"
    lccorrfileA = path + "/SgrA_correct_50ac_fpmA_lcsrccorrfile.fits"
    lccorrfileB = path + "/SgrA_correct_50ac_fpmB_lcsrccorrfile.fits"

    # --- Define Parameters / General Used Ones ---
    energy_min = 3.0
    # energy_max = 79.0
    energy_max = 30.0

    # --- For Step 3 / GTI Cleaning ---

    print("\n------ Start Data Processing Pipeline ------\n")

    ######################## Step 1: Load Data #######################
    # copy other fits files if one module is missing
    one_mod = False
    if event_file_a is None:
        one_mod = True
        present_mod = "B"
        # gti_file_a, event_file_a, lccorrfileA= one_module(gti_file_b, event_file_b, lccorrfileB,output_dir)
        # duplicate_fits(gti_file_b, event_file_b, lccorrfileB, output_dir)
        duplicate_fits(event_file_b, lccorrfileB, output_dir)
        # gti_file_a = output_dir + "gti_file_C"
        event_file_a = output_dir + "event_file_C"
        lccorrfileA = output_dir + "lccorrfile_C"

    elif event_file_b is None:
        one_mod = True
        present_mod = "A"
        # gti_file_b, event_file_b, lccorrfileB = one_module(gti_file_a, event_file_a, lccorrfileA,output_dir)
        # duplicate_fits(gti_file_a, event_file_a, lccorrfileA, output_dir)
        duplicate_fits(event_file_a, lccorrfileA, output_dir)
        # gti_file_b = output_dir + "gti_file_C"
        event_file_b = output_dir + "event_file_C"
        lccorrfileB = output_dir + "lccorrfile_C"

    print("Step 1: Loading GTI and event data...")
    gtiA = load_gti_file(event_file_a)
    gtiB = load_gti_file(event_file_b)
    eventsA = load_event_file(event_file_a)
    eventsB = load_event_file(event_file_b)
    print(f"    Module A: {len(gtiA)} GTI intervals, {len(eventsA)} events.")
    print(f"    Module B: {len(gtiB)} GTI intervals, {len(eventsB)} events.")

    print("Step 1.5: apply barycenter correction....")
    barycorr(barycorr_event, event_file_a, eventsA, eventsB, gtiA, gtiB)
    print("Step 1 Complete: Loaded GTI and event data for modules A and B.\n")

    # Debugging: Print energy ranges
    # print(f"Energy range in Module A: min={eventsA['Energy'].min()}, max={eventsA['Energy'].max()}")
    # print(f"Energy range in Module B: min={eventsB['Energy'].min()}, max={eventsB['Energy'].max()}")

    # print("\nInspecting the first 10 rows of Module A event data:")
    # print(eventsA.head(10))
    # print("\nInspecting the first 10 rows of Module B event data:")
    # print(eventsB.head(10))

    ################## Step 2: Filter Events by Energy ##################
    print("\nStep 2: Filtering events by energy...")
    filtered_eventsA = filter_events_by_energy(eventsA, energy_min, energy_max)
    filtered_eventsB = filter_events_by_energy(eventsB, energy_min, energy_max)

    print(f"    Module A: {len(filtered_eventsA)} events remaining.")
    print(f"    Module B: {len(filtered_eventsB)} events remaining.")
    print(f"    Filtered events for energy range [{energy_min}, {energy_max}] keV.")

    # Save or pass the filtered data for the next steps
    filtered_eventsA.to_csv(output_dir + "2_filtered_events_moduleA.csv", index=False)
    filtered_eventsB.to_csv(output_dir + "2_filtered_events_moduleB.csv", index=False)

    print("Step 2 Complete: Filtered event data saved for both modules.\n")

    ######################## Step 3: Clean GTIs ########################
    # Parameters for GTI cleaning
    gti_threshold = 30  # Minimum duration (seconds)
    gti_starttrim = 15  # Trim seconds from start
    gti_stoptrim = 15  # Trim seconds from stop

    print("\nStep 3: Cleaning GTIs...")
    print("    Module A: ")
    gtiA_cleaned = clean_gti(
        gtiA,
        filtered_eventsA,
        threshold=gti_threshold,
        starttrim=gti_starttrim,
        stoptrim=gti_stoptrim,
    )
    print("    Module B: ")
    gtiB_cleaned = clean_gti(
        gtiB,
        filtered_eventsB,
        threshold=gti_threshold,
        starttrim=gti_starttrim,
        stoptrim=gti_stoptrim,
    )

    # Save cleaned GTIs for debugging or further steps
    gtiA_cleaned.to_csv(output_dir + "3_gtiA_cleaned.csv", index=False)
    gtiB_cleaned.to_csv(output_dir + "3_gtiB_cleaned.csv", index=False)

    print("Step 3 Complete: Cleaned GTIs saved.\n")

    ######################## Step 4: Merge GTIs ########################
    print("\nStep 4: Merging GTIs...")
    merged_gti = merge_gtis(gtiA_cleaned, gtiB_cleaned)

    # Save merged GTIs for debugging or further steps
    merged_gti.to_csv(output_dir + "4_merged_gti.csv", index=False)

    print("Step 4 Complete: Merged GTIs saved to 'merged_gti.csv'.\n")

    ########### Step 5: Filtering Events Using Common GTIs #############
    print("\nStep 5: Filtering events using common GTIs...")
    print("    Module A: ")
    filtered_eventsA_with_gti = filter_events_with_common_gti(
        filtered_eventsA, merged_gti
    )
    print("    Module B: ")
    filtered_eventsB_with_gti = filter_events_with_common_gti(
        filtered_eventsB, merged_gti
    )

    # Save the filtered events for debugging or further analysis
    filtered_eventsA_with_gti.to_csv(
        output_dir + "5_filtered_eventsA_with_gti.csv", index=False
    )
    filtered_eventsB_with_gti.to_csv(
        output_dir + "5_filtered_eventsB_with_gti.csv", index=False
    )

    print("Step 5 Complete: Filtered events saved.\n")

    ########### Step 6: Get Correction Factors for Rvents #############
    print("\nStep 6: Getting correction factors for photon events...")

    # Apply the correction factor module for Module A
    print("    Module A: ")
    filtered_eventsA_with_gti["CORRECTION_FACTOR"] = get_event_corr_factor(
        lccorrfileA, filtered_eventsA_with_gti["TIME"]
    )

    # Apply the correction factor module for Module B
    print("    Module B: ")
    filtered_eventsB_with_gti["CORRECTION_FACTOR"] = get_event_corr_factor(
        lccorrfileB, filtered_eventsB_with_gti["TIME"]
    )

    # Save the updated DataFrames
    filtered_eventsA_with_gti.to_csv(
        output_dir + "6_filtered_eventsA_with_corr.csv", index=False
    )
    filtered_eventsB_with_gti.to_csv(
        output_dir + "6_filtered_eventsB_with_corr.csv", index=False
    )

    print("Step 6 Complete: Correction factors saved.\n")

    #################### Step 7: Merge Events ######################
    print("\nStep 7: Merging photon events...")

    if one_mod is True:
        print("making the missing module dataframe empty...")
        if present_mod == "A":
            filtered_eventsB_with_gti = empty_df(filtered_eventsB_with_gti)
        elif present_mod == "B":
            filtered_eventsA_with_gti = empty_df(filtered_eventsA_with_gti)

    # Merge photon events with exposure factors and module labels
    events_merged = merge_events(
        filtered_eventsA_with_gti,
        filtered_eventsB_with_gti,
        filtered_eventsA_with_gti["CORRECTION_FACTOR"].values,
        filtered_eventsB_with_gti["CORRECTION_FACTOR"].values,
    )

    # delete NANs if present
    if one_mod is True:
        events_merged = events_merged.dropna()

    # Save the merged DataFrame
    events_merged.to_csv(f"{output_dir}7_events_merged.csv", index=False)

    print("Step 7 Complete: Merged events saved.\n")

    ################# Step 8: Calculate Average Count Rate ###################

    # Configuration flag for count rate calculation
    calculate_average_count_rate = True  # Set to False to skip this step

    # Optional flare GTI DataFrame (None if not provided)
    flare_gti = None  # Replace with the flare GTI DataFrame if available

    print("\nStep 8: Calculating the average count rate during GTIs...")

    # Calculate average count rates
    average_rates = calculate_average_rate(
        events=events_merged,
        gti=merged_gti,
        flare_gti=flare_gti,
        calculate_average_count_rate=calculate_average_count_rate,
    )

    if average_rates is not None:
        # Save the results to a CSV
        average_rates.to_csv(f"{output_dir}8_average_count_rates.csv", index=False)
        print("Step 8 Complete: Average count rates saved.\n")

    ################# Step 9: Remove Gaps between GTIs ###################

    print("\nStep 9: Removing gaps between GTIs...")

    # Original observation end time
    original_tt_stop = merged_gti["STOP"].max()
    print("    Original observation stop time:", original_tt_stop)

    # Suppress gaps in event times
    events_no_gaps, updated_tt_stop, cumulative_gaps = suppress_gti_gaps(
        event_df=events_merged, gti_df=merged_gti, original_tt_stop=original_tt_stop
    )

    # Save the updated event DataFrame
    events_no_gaps.to_csv(f"{output_dir}9_events_no_gaps.csv", index=False)

    # Optionally save cumulative gaps for debugging
    cumulative_gaps_df = pd.DataFrame({"Cumulative Gap Time": cumulative_gaps})
    cumulative_gaps_df.to_csv(f"{output_dir}9_cumulative_gaps.csv", index=False)

    print("Step 9 Complete: Events with suppressed gaps saved.\n")

    ################# Step 10: Bayesian Block Analysis ###################

    # Config Parameters
    plan_a = False  # Plan A (Hardcode -- does not really exist -- DO NOT USE)
    plan_b = True  # Plan B (Custom Astropy)
    fp_rate = 0.01  # False positive rate
    # ncp_prior = None  # Prior number of change points
    ncp_prior = 4 - np.log10(
        fp_rate / (0.0136 * (len(events_no_gaps) ** 0.478))
    )  # As Shuo did
    do_iter = False  # Iterative refinement flag
    x_list = events_no_gaps["Exposure"].values
    # x_list=np.random.uniform(low=0.4, high=0.68, size=len(events_no_gaps['TIME'].values)) #for testing purposes
    print("\nStep 10: Bayesian Block Analysis...")

    if plan_a:
        print("    Using Custom Bayesian Block Implementation...")
        results = find_blocks(
            events_no_gaps["TIME"].values,
            events_no_gaps["Exposure"].values,
            fp_rate,
            ncp_prior,
            do_iter,
        )
        bayesian_blocks_df = format_bayesian_block_output(
            results, fp_rate, ncp_prior, do_iter
        )
        # Save results
        bayesian_blocks_df.to_csv(output_dir + "10_bayesian_blocks.csv", index=False)

    else:
        print("    Using Astropy Bayesian Block Implementation...")
        bayesian_blocks_df = bba_astropy(
            events_no_gaps["TIME"].values, ncp_prior, fp_rate, x_list=x_list
        )

        # Save results
        bayesian_blocks_df.to_csv(
            output_dir + "10_bayesian_blocks_astropy.csv", index=False
        )

    print("Step 10 Complete: Bayesian Block Analysis results saved.\n")

    ############ Step 10.1: Detailed Analysis of the Flaring Activity Block ############

    # Configuration for detailed flare analysis
    do_detailed_flare_analysis = True  # Do detailed analysis
    do_detailed_flare_analysis = False  # Not do detailed analysis
    flare_p0 = 0.5
    flaring_start = 262239084.32836443  # Configured flaring start time
    flaring_stop = 262241049.05592683  # Configured flaring stop time

    print("\nStep 10.1: Detailed Analysis of Flaring Activity Block...")

    if do_detailed_flare_analysis:
        # Create the flaring block
        flaring_block = events_no_gaps[
            (events_no_gaps["TIME"] >= flaring_start)
            & (events_no_gaps["TIME"] < flaring_stop)
        ]

        if flaring_block.empty:
            print(
                f"    No events found in the flaring interval [{flaring_start}, {flaring_stop}]. Skipping detailed analysis.\n"
            )
        else:
            print(
                f"    Flaring block extracted with {len(flaring_block)} events in the interval [{flaring_start}, {flaring_stop}]."
            )

            # Perform detailed flare analysis
            detailed_blocks_df = detailed_flare_analysis(
                event_df=flaring_block,  # Pass only the flaring block
                bayesian_blocks_df=bayesian_blocks_df,
                p0=flare_p0,
                flare_analysis_flag=do_detailed_flare_analysis,
            )

            if detailed_blocks_df is not None:
                # Save the refined flare analysis results
                detailed_blocks_df.to_csv(
                    output_dir + "10.1_detailed_flare_blocks.csv", index=False
                )
                print("Step 10.1 Complete: Detailed flare blocks saved.\n")
    else:
        print("Skipping Step 10.1 as per configuration.\n")

    print("\nStep 10.5: Put change points back into real time")

    ################ Step 11: Restore GTI Gaps to Blocks ######################

    print("\nStep 11: Insert Data Gaps into Bayesian Blocks...")

    try:
        # Adjust BBA blocks by reintroducing GTI gaps
        corrected_bba_blocks, gti_gaps_df = insert_gti_gaps(
            bba_df=bayesian_blocks_df,
            gti_df=merged_gti,
        )

        # Save the output files
        corrected_bba_csv = f"{output_dir}11_corrected_bba_blocks.csv"
        corrected_bba_blocks.to_csv(corrected_bba_csv, index=False)
        print("    Corrected BBA Blocks saved.")

        gti_gaps_csv = f"{output_dir}11_gti_gaps.csv"
        gti_gaps_df.to_csv(gti_gaps_csv, index=False)
        print("    GTI Gaps saved.")

        # Validate corrected BBA blocks against final event time
        final_bba_stop = corrected_bba_blocks["stop"].max()
        final_event_time = events_merged["TIME"].max()  # From Step 7

        print(f"    Final BBA Stop Time: {final_bba_stop}")
        print(f"    Final Event Time: {final_event_time}")

        if abs(final_bba_stop - final_event_time) > 1e-6:
            print(
                "Validation Failed: Corrected BBA blocks do not align with final event time."
            )
        else:
            print(
                "Validation Passed: Corrected BBA blocks align with final event time."
            )

    except Exception as e:
        print(f"Error in Step 11: {e}")

    print("Step 11 Complete: Gaps reintroduced into Bayesian Blocks.\n")

    ####################### Step 12: Save BBA Results #########################

    print("\nStep 12: Save Bayesian Block Analysis Results...")

    try:
        # Calculate event counts, block lengths, and rates for each BBA block
        corrected_bba_blocks, flare_blocks = calculate_event_counts_and_rates(
            corrected_bba_blocks, events_merged
        )

        # Calculate confidence limits for each block
        corrected_bba_blocks = calculate_confidence_limits(corrected_bba_blocks)

        # Observation metadata
        analysis_params = {
            "ncp_prior": ncp_prior,
            "fp_rate": fp_rate,
            "do_iter": do_iter,
        }

        # Define output file path
        output_file = f"{output_dir}12_gti_results.txt"

        # Custom optional notes for each block
        # save_results_note = [
        #     "Block 1 Note", "Block 2 Note", "Block 3 Note",
        #     "Block 4 Note", "Block 5 Note", "Block 6 Note",
        #     "Block 7 Note", "Block 8 Note", "Block 9 Note",
        #     "Block 10 Note", "Block 11 Note"
        # ]

        # Save BBA results as a text file
        save_bba_results_txt(
            bba_df=corrected_bba_blocks,
            analysis_params=analysis_params,
            output_file=output_file,
            # save_results_note=save_results_note
        )

        # save the flare version
        save_flare_results_txt(
            bba_df=flare_blocks,
            analysis_params=analysis_params,
            output_file=f"{output_dir}12_bba_results.txt",
            # save_results_note=save_results_note
        )

        print("Step 12 Complete: BBA results saved.\n")

    except Exception as e:
        print(f"Error in Step 12: {e}")

    ################# Step 13: Create Binned Light Curve ####################

    print("\nStep 13: Generate Regularly Binned Light Curve...")

    # try:
    binsize = 100  # Example bin size in seconds
    energy_range = (3.0, 79.0)  # Example energy range in keV

    # Generate the light curve
    lightcurve_df = generate_lightcurve(
        events_df=events_merged,  # Use the events_merged DataFrame from Step 7
        gti_df=merged_gti,  # Use the GTI DataFrame from Step 4
        binsize=binsize,
        energy_range=energy_range,
        gti_average=True,  # Set to True if GTI-based binning
    )

    # Save the light curve to a CSV file
    lightcurve_output_path = output_dir + f"13_LC_{binsize}.csv"
    lightcurve_df.to_csv(lightcurve_output_path, index=False)

    print("Step 13 Complete: Light curve saved.\n")

    # except Exception as e:
    #     print(f"Error in Step 13: {e}")

    ####################### Step 14: Plot Light Curve with Bayesian Blocks #########################

    print("\nStep 14: Plot Light Curve with Bayesian Blocks...")

    # Define file paths
    lc_plot_path = output_dir + "14_lightcurve_plot.png"
    lc_csv_path = output_dir + "14_lightcurve.csv"
    bb_csv_path = output_dir + "14_bayesian_blocks.csv"

    # prep for plotting by getting bb count rates in gti intervals
    plot_frame = plot_prep(corrected_bba_blocks, flare_blocks)
    plotframe_csv_path = output_dir + "14_plot_frame.csv"
    plot_frame.to_csv(plotframe_csv_path, index=False)

    # Plot and save light curve
    plot_lightcurve(
        lightcurve_df=lightcurve_df,  # Regularly binned light curve
        # bb_df=corrected_bba_blocks, #old call (basically gti defined)
        bb_df=plot_frame,  # Bayesian Blocks
        line_times=flare_blocks,
        output_path=lc_plot_path,
        lc_csv_path=lc_csv_path,
        bb_csv_path=bb_csv_path,
        convert_nustar_to_utc=convert_nustar_to_utc,
        plot_lines=True,
    )

    print("Step 14 Complete: Light Curve and Bayesian Blocks saved.\n")

    print("\n------ Creating ReadMe with run information ------\n")

    flags_dict = {
        "Barycenter correction": "True",
        "BBA min energy": energy_min,
        "BBA max energy": energy_max,
        "plotting energy range": energy_range,
        "input event file A": event_file_a,
        "input event file B": event_file_b,
        "input lccorr file A": lccorrfileA,
        "input lccorr file B": lccorrfileB,
    }

    write_readme(output_dir, flags_dict)

    print("\n------ End of Data Processing Pipeline ------\n")


main()
