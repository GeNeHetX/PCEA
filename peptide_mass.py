import pandas as pd
import os

from utils import get_expasy_peptide_mass, get_uniprotkbid


# Peptides mass path
path_peptides = "data/MALDI_IHC/results/peptides_mass_matrisome"

# Create the directory if it does not exist
os.makedirs(path_peptides, exist_ok=True)

# Load the matrisome genes
matrisome_genes = pd.read_csv("Hs_Matrisome_Masterlist_2012.csv")

# Get the gene names for the top and bottom genes and remove duplicates
genes = matrisome_genes['Gene Symbol']

# Get the peptide mass for each gene
for gene in genes:
    print(f"\nGene: {gene}")

    # Get a list of protein names for the gene
    protein_names = get_uniprotkbid(gene)
    
    # Get the peptide mass for each protein name
    for protein_name in protein_names:
        print(f"Protein: {protein_name}")

        # Check if the file already exists
        if not os.path.exists(f"{path_peptides}/{protein_name}_peptide_mass.csv"):

            # Get the peptide mass
            peptide_mass = get_expasy_peptide_mass(protein_name)
            print(f"{peptide_mass.shape[0]} peptides")

            # Save the peptide mass to a csv file
            peptide_mass.to_csv(f"{path_peptides}/{protein_name}_peptide_mass.csv", index=False)