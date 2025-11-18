import os
from datetime import datetime


def write_readme(output_dir, flags_dict, filename="README.txt", mode="w"):
    """
    Write a README text file with timestamp and user-defined flags

    Parameters:
        output_dir (str): Directory where README file will be saved.
        flags_dict (dict): Dictionary of flags/inputs to record.
        filename(str): Name of the README file (default "README.txt").
        mode (str): File mode: "a" to append, "w" to overwrite (default "w").

    """

    os.makedirs(output_dir, exist_ok=True)  # Ensure directory exists
    filepath = os.path.join(output_dir, filename)

    with open(filepath, mode) as f:
        f.write("=== Run Information ===\n")
        f.write(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("Flags and Inputs:\n")
        for key, value in flags_dict.items():
            f.write(f"  {key}: {value}\n")
        f.write("\n")  # Blank line for readability

    print(f"README written to {filepath}")
