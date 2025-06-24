import pandas as pd
import yaml
import os
import gc

# Load the configuration file
with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)

# Generate a list of slide names from the specified directory
slides = os.listdir(config['path_to_data'])

for slide in slides:
    # Load the peaks and pixels data for the specified slide
    print(f"Loading data for slide: {slide}")
    peaks = pd.read_csv(f"{config['path_to_data']}/{slide}/results/mse_peaks.csv", index_col=0)
    pixels = pd.read_feather(f"{config['path_to_data']}/{slide}/results/mse_pixels.feather")

    # Reset the index
    peaks.reset_index(drop=True, inplace=True)
    pixels.reset_index(drop=True, inplace=True)

    print("Cleaning peaks and pixels data...")

    # Exclude the pixels outside the lesion
    peaks = peaks[pixels["Density_Lesion"] > 0.5]
    pixels = pixels[pixels["Density_Lesion"] > 0.5]

    # Exclude the pixels with defects
    peaks = peaks[pixels["Density_Defects"] < 0.1]
    pixels = pixels[pixels["Density_Defects"] < 0.1]

    # Reset the index
    peaks.reset_index(drop=True, inplace=True)
    pixels.reset_index(drop=True, inplace=True)

    # Create a DataFrame to store the peaks information
    print("Calculating peaks statistics...")
    peaks_mz = pd.DataFrame({
        'mass': peaks.columns,
        'count': (peaks > 0).sum(axis=0).values,
        'frequency': (peaks > 0).mean(axis=0).values,
        'max_intensity': peaks.max(axis=0).values,
        'mean_intensity': (peaks[peaks != 0]).mean(axis=0).values,
        'std_intensity': (peaks[peaks != 0]).std(axis=0).values,
        'median_intensity': (peaks[peaks != 0]).median(axis=0).values})

    # Compute the correlation between the detected peaks and IHC pixels Density_Collagen, Density_CD8, Density_Tumor, Density_Stroma
    print("Calculating correlations...")
    correlation = pd.DataFrame({**{f'Pearson_{density}': peaks[peaks != 0][pixels["Density_Tumor"] < 0.1].corrwith(pixels[density][pixels["Density_Tumor"] < 0.1], method='pearson') for density in ['Density_CD8', 'Density_Collagen']},
                                **{f'Spearman_{density}': peaks[peaks != 0][pixels["Density_Tumor"] < 0.1].corrwith(pixels[density][pixels["Density_Tumor"] < 0.1], method='spearman') for density in ['Density_CD8', 'Density_Collagen']},
                                **{f'Pearson_{density}': peaks[peaks != 0].corrwith(pixels[density], method='pearson') for density in ['Density_Tumor', 'Density_Stroma']},
                                **{f'Spearman_{density}': peaks[peaks != 0].corrwith(pixels[density], method='spearman') for density in ['Density_Tumor', 'Density_Stroma']}})

    # Reset the index
    correlation.reset_index(inplace=True, drop=True)

    # Concatenate the correlation results with the peaks_mz DataFrame
    print("Concatenating correlation results with the statistics DataFrame...")
    peaks_mz = pd.concat([peaks_mz, correlation], axis=1)

    # Save the peaks_mz DataFrame to a CSV file
    print("Saving the peaks statistics and correlation results to CSV...")
    peaks_mz.to_csv(f"{config['path_to_data']}/{slide}/results/peaks_mz.csv", index=False)

    # Clear memory
    del peaks, pixels, peaks_mz, correlation
    gc.collect()