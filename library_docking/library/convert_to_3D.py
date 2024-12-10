import os
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import SDMolSupplier

# File paths
sdf_file = "AD194827_Del-1-1.sdf"  # Path to the large SDF file
output_smiles_file = "selected_smiles_output.txt"  # File to save SMILES strings
output_sdf_dir = "3d_sdf"  # Directory to save SDF files from PubChem

# Step 1: Read the SDF file
if not os.path.exists(sdf_file):
    print(f"Error: The file '{sdf_file}' does not exist.")
    exit()

supplier = SDMolSupplier(sdf_file)
os.makedirs(output_sdf_dir, exist_ok=True)  # Create directory for SDF files from PubChem

# Step 2: Extract specific molecules based on indices, get SMILES, and download SDFs
with open(output_smiles_file, 'w') as smiles_file:
    for i, mol in enumerate(supplier, start=1):  # Start from 1 to match the index in SDF
        if mol is None:
            print(f"Warning: Molecule {i} is invalid and will be skipped.")
            continue  # Skip invalid molecules

        smiles = Chem.MolToSmiles(mol)
        smiles_file.write(f"Molecule {i}: {smiles}\n")
        print(f"Extracted SMILES for Molecule {i}: {smiles}")

        # Step 3: Use PubChemPy to search PubChem and download SDF
        try:
            compound = pcp.get_compounds(smiles, 'smiles')
            if compound:
                cid = compound[0].cid
                sdf_output_file = os.path.join(output_sdf_dir, f"molecule_{i}.sdf")
                
                # Fetch the SDF data
                sdf_data = pcp.get_sdf(cid, record_type='3d')
                
                if sdf_data:
                    with open(sdf_output_file, 'w') as sdf_file:
                        sdf_file.write(sdf_data)
                    print(f"Successfully downloaded SDF for SMILES: {smiles}, Index: {i}")
                else:
                    print(f"Failed to retrieve SDF data for CID: {cid}")

        except pcp.PubChemHTTPError as e:
            print(f"Error occurred while searching PubChem for SMILES: {smiles}. Error: {e}")
        except Exception as e:
            print(f"Unexpected error occurred for SMILES: {smiles}. Error: {e}")

print(f"\nSMILES strings successfully written to '{output_smiles_file}' and SDF files downloaded to '{output_sdf_dir}'")