import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import pandas as pd


def plot_lightcurve(
    lightcurve_df,
    bb_df,
    output_path,
    lc_csv_path,
    bb_csv_path,
    line_times,
    convert_nustar_to_utc,
    plot_lines=False,
    title="Light Curve",
    xlabel="Time (UTC)",
    ylabel="Count Rate (cts/s)",
):
    """
    Plot the regularly binned light curve with Bayesian Block overlay, convert times to UTC, and save final data.

        - Saves the light curve plot as an image file.
        - Saves the processed Light Curve DataFrame (`lightcurve_df`) as a `.csv` file.
        - Saves the processed Bayesian Block DataFrame (`bb_df`) as a `.csv` file.

    Parameters:
        lightcurve_df (pd.DataFrame):
        bb_df (pd.DataFrame):
        output_path (str):
        lc_csv_path (str): output path for lightcurve csv
        bb_csv_path (str):output path for bayesian block csv
        line_times (array): UTC times for vertical lines to plot
        convert_nustar_to_utc (function): Convert NuSTAR mission time to UTC.
        plot_lines (Boolean): (default False)
        title (str): (default "Light Curve")
        xlabel (str): (default "Time (UTC))
        ylabel (str): (default "Count Rate (cts/s)")
    """

    if lightcurve_df.empty:
        raise ValueError("Error: Light curve DataFrame is empty. Cannot plot.")

    # ---- 1️ Compute Bin Centers and Error Bars in Mission Time (FLOATS) ----

    bin_centers = (lightcurve_df["bin_start"] + lightcurve_df["bin_end"]) / 2
    bin_centers = mdates.date2num(bin_centers.apply(convert_nustar_to_utc))
    bin_half_widths = (lightcurve_df["bin_end"] - lightcurve_df["bin_start"]) / 2
    bin_half_widths = bin_half_widths / 86400

    vertical_errors = [
        lightcurve_df["count_rate"] - lightcurve_df["lower_limit"],  # Lower error
        lightcurve_df["upper_limit"] - lightcurve_df["count_rate"],  # Upper error
    ]

    # ---- 2️ Convert Time to UTC (AFTER Calculations) ----

    lightcurve_df["UTC_bin_start"] = lightcurve_df["bin_start"].apply(
        convert_nustar_to_utc
    )
    lightcurve_df["UTC_bin_end"] = lightcurve_df["bin_end"].apply(convert_nustar_to_utc)

    # force conversion to datetime
    # lightcurve_df['UTC_bin_start'] = mdates.date2num(lightcurve_df['UTC_bin_start'])
    # lightcurve_df['UTC_bin_end'] = mdates.date2num(lightcurve_df['UTC_bin_end'])

    # Convert bin centers **AFTER computing as float**
    # lightcurve_df['UTC_bin_center'] = pd.Series(bin_centers).apply(convert_nustar_to_utc)
    lightcurve_df["UTC_bin_center"] = pd.Series(bin_centers)

    # Ensure only UTC columns are strings
    lightcurve_df["UTC_bin_start"] = lightcurve_df["UTC_bin_start"].astype(str)
    lightcurve_df["UTC_bin_end"] = lightcurve_df["UTC_bin_end"].astype(str)
    lightcurve_df["UTC_bin_center"] = lightcurve_df["UTC_bin_center"].astype(str)

    # ---- 3️ Convert Bayesian Blocks Time to UTC ----
    if not bb_df.empty:
        bb_df["UTC_start"] = bb_df["start"].apply(convert_nustar_to_utc)  # .astype(str)
        bb_df["UTC_stop"] = bb_df["stop"].apply(convert_nustar_to_utc)  # .astype(str)

        # force datetime conversion
        bb_df["UTC_start"] = mdates.date2num(bb_df["UTC_start"])
        bb_df["UTC_stop"] = mdates.date2num(bb_df["UTC_stop"])

        line_times["UTC_start"] = line_times["start"].apply(
            convert_nustar_to_utc
        )  # .astype(str)
        line_times["UTC_stop"] = line_times["stop"].apply(
            convert_nustar_to_utc
        )  # .astype(str)

        line_times["UTC_start"] = mdates.date2num(line_times["UTC_start"])
        line_times["UTC_stop"] = mdates.date2num(line_times["UTC_stop"])

    # Extract YYYY-MM-DD for title
    observation_date = lightcurve_df["UTC_bin_start"].iloc[0].split("T")[0]
    # observation_date = lightcurve_df["UTC_bin_start"][0] #.astype(str)

    # ---- 4️ Plot the Light Curve ----
    plt.figure(figsize=(12, 6))

    # Plot regularly binned light curve
    plt.errorbar(
        bin_centers,  # Use FLOAT values, NOT UTC
        lightcurve_df["count_rate"],
        xerr=bin_half_widths,
        yerr=vertical_errors,
        fmt="o",
        color="black",
        ecolor="silver",
        capsize=1,
        elinewidth=1,
        markersize=0.3,
        label="Binned LC",
    )

    # Plot Bayesian Block horizontal lines (red)
    if not bb_df.empty:
        for _, row in bb_df.iterrows():
            plt.plot(
                [row["UTC_start"], row["UTC_stop"]],
                [row["rate"], row["rate"]],
                color="red",
                linewidth=0.7,
                label="Bayesian Blocks" if _ == 0 else None,
            )
            plt.plot(
                [row["UTC_start"], row["UTC_stop"]],
                [row["upperlim"], row["upperlim"]],
                color="red",
                linestyle=(0, (5, 5)),
                linewidth=0.7,
            )  # Dashed
            plt.plot(
                [row["UTC_start"], row["UTC_stop"]],
                [row["lowerlim"], row["lowerlim"]],
                color="red",
                linestyle=(0, (5, 5)),
                linewidth=0.7,
            )  # Dashed

            # if plot_lines is True:
            #     plt.axvline(
            #         x=row["UTC_start"], color="green", linestyle="--", linewidth=0.7
            #     )
            #     plt.axvline(
            #         x=row["UTC_stop"], color="red", linestyle="--", linewidth=0.7
            #     )

    # if not line_times.empty:
    #     for _, row in line_times.iterrows():
    #         plt.axvline(
    #             x=row["UTC_start"], color="green", linestyle="--", linewidth=0.7
    #         )
    #         plt.axvline(x=row["UTC_stop"], color="red", linestyle="--", linewidth=0.7)

    # Labels, title, and legend
    plt.title(f"{title} ({observation_date})")  # Include observation date in title
    plt.xlabel(xlabel)

    # plt.ylim(0, 2.6)

    # xticks = np.arange(bb_df["UTC_start"].iloc[0], bb_df["UTC_stop"].iloc[-1], 0.015)
    # xlabels = [f'\\${x:1.2f}' for x in xticks]
    # plt.xticks(xticks, labels=xlabels)
    # plt.xlim(17267.303013,17267.369417 )

    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    plt.gcf().autofmt_xdate()
    plt.ylabel(ylabel)
    plt.legend()

    # Save and close plot
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print("    Light Curve plot saved.")

    # ---- 5️ Save Light Curve and Bayesian Blocks as CSV ----
    lightcurve_df.to_csv(lc_csv_path, index=False)
    print("    Light Curve CSV saved.")

    if not bb_df.empty:
        bb_df.to_csv(bb_csv_path, index=False)
        print("    Bayesian Blocks CSV saved.")
