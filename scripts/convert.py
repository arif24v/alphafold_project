from rdkit import Chem
from rdkit.Chem import AllChem
import os

# Function to convert SDF to PDB
def sdf_to_pdb(input_sdf, output_pdb):
    # Read SDF file
    suppl = Chem.SDMolSupplier(input_sdf)
    
    # Iterate over molecules in SDF file
    for mol in suppl:
        if mol is not None:
            # Generate 3D coordinates if not already present
            AllChem.EmbedMolecule(mol)
            AllChem.MMFFOptimizeMolecule(mol)

            # Write molecule to PDB file
            writer = Chem.PDBWriter(output_pdb)
            writer.write(mol)
            writer.close()

# Directory containing SDF files
sdf_directory = 'sdf'

# Output directory for PDB files
pdb_directory = 'pdb'

# Loop through each SDF file in the directory
for filename in os.listdir(sdf_directory):
    if filename.endswith('.sdf'):
        input_sdf = os.path.join(sdf_directory, filename)
        output_pdb = os.path.join(pdb_directory, os.path.splitext(filename)[0] + '.pdb')
        
        # Convert SDF to PDB
        sdf_to_pdb(input_sdf, output_pdb)

print("Conversion complete.")
