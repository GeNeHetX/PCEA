import pandas as pd
import os
import gseapy as gp
import yaml

from utils import create_protein_peptide_map, match_protein_peptide_map

# Load the configuration file
with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)

# Define paths
path_peptides = "data/MALDI_IHC/results/peptides_mass_matrisome"
path_pcea = "data/MALDI_IHC/results/peptide_composition_enrichment_matrisome"

# Create the directory for results if it does not exist
os.makedirs(path_pcea, exist_ok=True)

slides = os.listdir(config['path_to_data'])

for slide in slides:
    # Load the peaks statistics
    print(f"Processing slide: {slide}")
    print("Loading peaks statistics...")
    peaks_mz = pd.read_csv(f"{config['path_to_data']}/{slide}/results/peaks_mz.csv")

    # Create a protein peptide map from the peptide path
    print("Creating protein peptide map...")
    protein_peptide_map = create_protein_peptide_map(path_peptides, type=str)

    # Match the protein peptide map with the peaks
    print("Matching protein peptide map with peaks...")
    protein_peptide_map = match_protein_peptide_map(protein_peptide_map, peaks_mz['mass'], tolerance=0.01)

    # Create a ranked list DataFrame with the correlation
    print("Creating ranked list DataFrame...")
    # ranked_list_df = peaks_mz[['mass', 'Pearson_CD8']].copy()
    ranked_list_df = peaks_mz[['mass', 'Spearman_Density_CD8']].copy()

    # Remove the na values
    ranked_list_df = ranked_list_df.dropna()

    # Change the mass to a string type
    ranked_list_df['mass'] = ranked_list_df['mass'].astype(str)

    # Reset the index
    ranked_list_df.reset_index(drop=True, inplace=True)

    # Perform PCEA prerank analysis
    print("Performing prerank analysis...")
    prerank_results = gp.prerank(
        rnk=ranked_list_df,         # Your ranked list DataFrame
        gene_sets=protein_peptide_map, # Your protein-peptide mapping dictionary
        min_size=3,                # Minimum number of peptides from a protein present in the ranked list to consider the protein (adjust as needed)
        max_size=500,               # Maximum number of peptides (adjust as needed)
        permutation_num=10**5,       # Number of permutations for significance testing (1000 is standard)
        threads=40,
        outdir=None,               # Set to a directory path string to save results (e.g., './gsea_results')
        seed=42,                   # For reproducibility
        verbose=True               # Print progress
    )

    # Get the main results DataFrame
    print("Arranging results DataFrame...")
    results_df = prerank_results.res2d

    # Rename the columns for clarity
    results_df.rename(columns={'Name': 'Slide', 'Term': 'Protein', 'Gene %': 'Peptide %', 'Lead_genes': 'Lead_peptides'}, inplace=True)

    # Add the slide name to the results
    results_df['Slide'] = slide

    # Save the results to a CSV file
    print(f"Saving results for slide: {slide}")
    results_df.to_csv(f"{path_pcea}/{slide}.csv", index=False)