import pandas as pd


def merge_gtis(gtiA, gtiB):
    """
    Merge GTIs from modules A and B by finding overlapping intervals.

    Parameters:
        gtiA (pd.DataFrame): Cleaned GTI DataFrame from module A with 'START' and 'STOP'.
        gtiB (pd.DataFrame): Cleaned GTI DataFrame from module B with 'START' and 'STOP'.

    Returns:
        (pd.DataFrame): Merged GTI DataFrame with 'START' and 'STOP' columns.
    """
    try:
        merged_gti = []

        # Iterate through all GTIs in A and B
        for _, rowA in gtiA.iterrows():
            for _, rowB in gtiB.iterrows():
                # Find the overlap between the two intervals
                overlap_start = max(rowA["START"], rowB["START"])
                overlap_stop = min(rowA["STOP"], rowB["STOP"])

                # Add the interval if there is a valid overlap
                if overlap_start < overlap_stop:
                    merged_gti.append({"START": overlap_start, "STOP": overlap_stop})

        # Convert to DataFrame
        merged_gti = pd.DataFrame(merged_gti)

        print(f"    GTIs in module A: {len(gtiA)}")
        print(f"    GTIs in module B: {len(gtiB)}")
        print(f"    GTIs after merging: {len(merged_gti)}")

        return merged_gti
    except Exception as e:
        raise RuntimeError(f"Error during GTI merging: {e}")
