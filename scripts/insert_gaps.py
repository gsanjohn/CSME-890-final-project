import pandas as pd


def insert_gti_gaps(bba_df, gti_df):
    """
    Adjust BBA blocks to account for GTI gaps, splitting blocks as necessary.

    Parameters:
        bba_df (pd.DataFrame): Bayesian Blocks DataFrame with 'start' and 'stop' columns.
        gti_df (pd.DataFrame): GTI DataFrame with 'START' and 'STOP' columns.

    Returns:
        (pd.DataFrame): Corrected BBA DataFrame with gaps reintroduced and blocks split if needed.
        (pd.DataFrame): GTI gaps DataFrame with calculated gap durations.
    """
    # print(f"BBA Blocks Time Range: {bba_df['start'].min()} - {bba_df['stop'].max()}")
    # print(f"GTI Time Range: {gti_df['START'].min()} - {gti_df['STOP'].max()}")

    # Step 1: Identify Gaps Between GTIs
    gti_gaps = []
    for i in range(len(gti_df) - 1):
        gap_start = gti_df.iloc[i]["STOP"]
        gap_stop = gti_df.iloc[i + 1]["START"]
        gap_duration = gap_stop - gap_start

        if gap_duration > 0:
            gti_gaps.append(
                {
                    "GAP_START": gap_start,
                    "GAP_STOP": gap_stop,
                    "GAP_DURATION": gap_duration,
                }
            )
            # print(f"Identified Gap: START={gap_start}, STOP={gap_stop}, DURATION={gap_duration}")

    gti_gaps_df = pd.DataFrame(gti_gaps)

    # Step 2: Adjust BBA Blocks Iteratively
    cumulative_shift = 0
    corrected_blocks = bba_df.copy()

    for _, gap in gti_gaps_df.iterrows():
        gap_start = gap["GAP_START"]
        gap_stop = gap["GAP_STOP"]
        gap_duration = gap["GAP_DURATION"]
        print(
            f"Processing Gap: START={gap_start}, STOP={gap_stop}, DURATION={gap_duration}"
        )

        new_corrected_blocks = []

        for _, block in corrected_blocks.iterrows():
            block_start = block["start"]
            block_stop = block["stop"]

            if block_stop <= gap_start:
                new_corrected_blocks.append(
                    {**block.to_dict(), "start": block_start, "stop": block_stop}
                )
                print(
                    f"Block Unaffected by Gap: start={block_start}, stop={block_stop}"
                )
            elif block_start >= gap_start:
                new_corrected_blocks.append(
                    {
                        **block.to_dict(),
                        "start": block_start + gap_duration,
                        "stop": block_stop + gap_duration,
                    }
                )
                print(
                    f"Block Shifted After Gap: start={block_start + gap_duration}, stop={block_stop + gap_duration}"
                )
            else:
                print(f"Block Overlaps Gap: start={block_start}, stop={block_stop}")
                if block_start < gap_start:
                    new_corrected_blocks.append(
                        {**block.to_dict(), "start": block_start, "stop": gap_start}
                    )
                    print(f"  Sub-block 1: start={block_start}, stop={gap_start}")
                new_corrected_blocks.append(
                    {
                        **block.to_dict(),
                        "start": gap_stop,
                        "stop": block_stop + gap_duration,
                    }
                )
                print(
                    f"  Sub-block 2: start={gap_stop}, stop={block_stop + gap_duration}"
                )

        corrected_blocks = pd.DataFrame(new_corrected_blocks)
        cumulative_shift += gap_duration
        print(f"Cumulative Shift Updated: {cumulative_shift} seconds")

    return corrected_blocks, gti_gaps_df
