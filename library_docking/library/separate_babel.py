import os
import subprocess
from rdkit import Chem

# Specify the input SDF file
input_sdf = "AD194827_Del-1-1.sdf"  # Replace with your SDF filename

# Specify the output directory
output_dir = "pdbqt_files"
os.makedirs(output_dir, exist_ok=True)

# Read the SDF file using RDKit
supplier = Chem.SDMolSupplier(input_sdf, removeHs=False)

for idx, mol in enumerate(supplier):
    if mol is None:
        print(f"Warning: Molecule at index {idx} could not be read and will be skipped.")
        continue

    # Generate a unique molecule name or use an existing property
    mol_name = f"molecule_{idx+1}"

    # Write the molecule to a temporary SDF file
    temp_sdf = os.path.join(output_dir, f"{mol_name}.sdf")
    writer = Chem.SDWriter(temp_sdf)
    writer.write(mol)
    writer.close()

    # Define the output PDBQT filename
    output_pdbqt = os.path.join(output_dir, f"{mol_name}.pdbqt")

    # Construct the Open Babel command
    command = [
        'obabel', temp_sdf,
        '-O', output_pdbqt,
        '--partialcharge', 'gasteiger',
        '-xp'  # Assign AutoDock atom types
    ]

    # Run the Open Babel command
    try:
        subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(f"Processed and saved {output_pdbqt}")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while converting {temp_sdf}:")
        print(e.stderr.decode())
    except FileNotFoundError:
        print("Error: The 'obabel' command was not found. Please ensure that Open Babel is installed and added to your system's PATH.")

    # Remove the temporary SDF file
    os.remove(temp_sdf)

print(f"Conversion complete. PDBQT files are saved in '{output_dir}' directory.")
