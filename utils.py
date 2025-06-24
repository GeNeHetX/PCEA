import numpy as np
import pandas as pd
import requests
import sys
import os

from bs4 import BeautifulSoup


def get_url(url, **kwargs):
    """
    Send a GET request to a given URL and return the response.
    
    Args:
        url (str): URL to send the GET request to.
        **kwargs: Additional keyword arguments to pass to requests.get.
        
    Returns:
        requests.models.Response: Response object.
    """
    response = requests.get(url, **kwargs);

    if not response.ok:
        print(response.text)
        response.raise_for_status()
        sys.exit()

    return response


def get_uniprotkbid(gene_name: str, taxonomy_id: int = 9606, reviewed: bool = True) -> list[str]:
    """
    Get a list of UniProtKB IDs for a given gene name.

    Args:
        gene_name (str): Gene name.
        taxonomy_id (int, optional): Taxonomy ID. Defaults to 9606 Homo sapiens (Human).
        reviewed (bool, optional): Consider only reviewed entries. Defaults to True.

    Returns:
        list[str]: List of UniProtKB IDs.
    """
    # Query parameters
    website = "https://rest.uniprot.org/uniprotkb"
    reviewed = " AND (reviewed:true)" if reviewed else ""

    # Run a search query for a protein
    response = get_url(f"{website}/search?query=(gene:{gene_name}) AND (taxonomy_id:{taxonomy_id}){reviewed}&fields=id&size=100")

    # # Store the protein names
    protein_names = [entry["uniProtkbId"] for entry in response.json()["results"]]

    return protein_names


def get_expasy_peptide_mass(protein_name: str, min_mass: float = 500, max_mass: float = 3500, enzyme: str = "Trypsin") -> pd.DataFrame:
    """
    Fetch peptide mass data from Expasy PeptideMass service.

    Args:
        protein_name (str): Name of the protein.
        min_mass (float, optional): Minimum mass of peptides. Defaults to 500.
        max_mass (float, optional): Maximum mass of peptides. Defaults to 3500.
        enzyme (str, optional): Enzyme used for digestion. Defaults to "Trypsin".

    Returns:
        pd.DataFrame: DataFrame containing peptide mass, position and sequence.
    """
    # Define the URL for the Expasy PeptideMass service
    url = "https://web.expasy.org/cgi-bin/peptide_mass/peptide-mass.pl"

    # Define the parameters for the PeptideMass service
    params = {
        "protein": protein_name,
        "reagents": "nothing+(in+reduced+form)",
        "mplus": "mh",
        "masses": "monoisotopic",
        "enzyme": enzyme,
        "MC": "0",
        "minmass": str(min_mass),
        "maxmass": str(max_mass),
        "order": "mass"
    }

    # Send a GET request to the Expasy PeptideMass service
    response = get_url(url, params=params)

    # Parse the HTML response
    soup = BeautifulSoup(response.text, 'html.parser')

    # Find the table with class "type-1"
    table = soup.find('table', {'class': 'type-1'})

    # Extract the table rows
    rows = table.find_all('tr')

    # Extract the table headers
    headers = [header.text for header in rows[0].find_all('th')]

    # Extract the table data
    data = []
    for row in rows[1:]:
        cells = row.find_all('td')
        data.append([cell.text for cell in cells])

    # Create a DataFrame from the data
    df = pd.DataFrame(data=data, columns=headers)

    return df

def proteins_share(path_peptides: str, mass_arr: np.ndarray, tolerance: float=0.2) -> dict:
    """Extracts the proteins shared masses with the mass array.

    Args:
        path_peptides (str): path to the proteins peptides masses csv files
        mass_arr (np.ndarray): an array with the masses
        tolerance (float): The tolerance to consider a mass shared between the proteins. (Default is 0.2)

    Returns:
        dict: A dictionary with the proteins names as keys and the shared masses as values.
    """
    # Get the names of the proteins
    protein_names = [name.split('_')[0] for name in os.listdir(path_peptides)]

    # Initialize the dictionary 
    proteins_share = {name: [] for name in protein_names}

    # Iterate over the proteins
    for protein_name in protein_names:

        # Load the masses of the protein
        protein_mass = pd.read_csv(f"{path_peptides}/{protein_name}_HUMAN_peptide_mass.csv")['mass'].values

        # Iterate over the masses of the protein
        for mass in protein_mass:

            # Check if the mass is in the mass array with a tolerance
            if np.min(np.abs(mass - mass_arr)) < tolerance:

                # Append the mass to the list of masses
                proteins_share[protein_name].append(str(mass))

    # Order the proteins by the number of masses
    proteins_share = {k: v for k, v in sorted(proteins_share.items(), key=lambda item: len(item[1]), reverse=True)}
    return proteins_share


def create_protein_peptide_map(path_peptides: str, type=float) -> dict:
    """
    Creates a mapping of protein names to their corresponding peptide masses.

    Args:
        path_peptides (str): Path to the directory containing peptide mass CSV files.

    Returns:
        dict: A dictionary where keys are protein names and values are lists of peptide masses.
    """
    # Extract protein names from the filenames in the specified directory
    protein_names = [name.split('_')[0] for name in os.listdir(path_peptides)]
    
    # Create a dictionary mapping each protein name to its list of peptide masses
    protein_peptide_map = {
        protein_name: pd.read_csv(f"{path_peptides}/{protein_name}_HUMAN_peptide_mass.csv")['mass'].astype(type).to_list()
        for protein_name in protein_names
    }
    
    return protein_peptide_map


def match_protein_peptide_map(protein_peptide_map: dict, match_list: np.ndarray, tolerance: float=0.1) -> dict:
    """
    Updates the values in protein_peptide_map to match the mz values in match_list within a given tolerance.

    Args:
        protein_peptide_map (dict): A dictionary where keys are protein names and values are lists of peptide masses.
        match_list (np.ndarray): An array containing mz values to match.
        tolerance (float): The tolerance within which mz values are considered a match.

    Returns:
        dict: Updated protein_peptide_map with values replaced by matching mz values from match_list.
    """
    updated_map = {}

    mz_values = match_list.astype(float)

    for protein in protein_peptide_map.keys():
        updated_map[protein] = protein_peptide_map[protein][:]
        for i, mz in enumerate(updated_map[protein]):
            if mz in mz_values:
                continue
            elif any(abs(float(mz) - float(mz_value)) <= tolerance for mz_value in mz_values):
                updated_map[protein][i] = str(mz_values.values[np.argmin(np.abs(mz_values - float(mz)))])

    return updated_map