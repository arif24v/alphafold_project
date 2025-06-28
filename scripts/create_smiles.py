import os
import csv
from rdkit import Chem

def convert_sdf_to_smiles(input_directory, output_csv):
    # Prepare the output data
    output_data = []

    # Walk through the directory to find sdf files
    for root, _, files in os.walk(input_directory):
        for file in files:
            if file.endswith('.sdf'):
                file_path = os.path.join(root, file)
                try:
                    # Read the .sdf file
                    suppl = Chem.SDMolSupplier(file_path)
                    
                    for mol in suppl:
                        if mol is None:
                            continue  # Skip invalid molecules
                        smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
                        output_data.append([file, smiles])
                except Exception as e:
                    print(f"Error converting {file}: {e}")
                    output_data.append([file, "Conversion failed"])

    # Write the output CSV
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Filename", "SMILES"])
        writer.writerows(output_data)

    print(f"Conversion complete. SMILES strings saved to {output_csv}")

# Specify the input directory and output CSV file
input_directory = "./../library_docking/library/3d_sdf"
output_csv = "output_smiles.csv"

convert_sdf_to_smiles(input_directory, output_csv)