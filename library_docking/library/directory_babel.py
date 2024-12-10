import os
import subprocess
from rdkit import Chem

# Specify the input directory containing SDF files
input_dir = "3d_sdf"  # Replace with your directory containing SDF files

# Specify the output directory
output_dir = "3d_pdbqt_files"
os.makedirs(output_dir, exist_ok=True)

# Iterate over all SDF files in the input directory
for sdf_file in os.listdir(input_dir):
    if not sdf_file.endswith(".sdf"):
        continue  # Skip non-SDF files

    input_sdf = os.path.join(input_dir, sdf_file)

    # Read the SDF file using RDKit
    supplier = Chem.SDMolSupplier(input_sdf, removeHs=False)

    # Ensure the SDF contains exactly one molecule
    mols = [mol for mol in supplier if mol is not None]
    if len(mols) != 1:
        print(f"Warning: {sdf_file} contains {len(mols)} molecules. Skipping...")
        continue

    mol = mols[0]  # There's exactly one molecule

    # Generate the output PDBQT filename matching the SDF base name
    base_name = os.path.splitext(sdf_file)[0]
    output_pdbqt = os.path.join(output_dir, f"{base_name}.pdbqt")

    # Write the molecule to a temporary SDF file
    temp_sdf = os.path.join(output_dir, f"{base_name}.sdf")
    writer = Chem.SDWriter(temp_sdf)
    writer.write(mol)
    writer.close()

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
