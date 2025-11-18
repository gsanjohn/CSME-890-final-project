def filter_events_by_energy(data, energy_min, energy_max):
    """
    Filter events based on energy range.

    Parameters:
        data (pd.DataFrame): DataFrame containing event data.
        energy_min (float): Minimum energy threshold.
        energy_max (float): Maximum energy threshold.

    Returns:
        (pd.DataFrame): Filtered DataFrame with events within the specified energy range.
    """
    try:
        filtered_data = data[
            (data["Energy"] >= energy_min) & (data["Energy"] <= energy_max)
        ]
        return filtered_data
    except KeyError:
        raise KeyError("The input data does not contain the 'Energy' column.")
