from rdkit import Chem
from rdkit.Chem import AllChem
import os
from meeko import MoleculePreparation  # Import Meeko's MoleculePreparation

# Load the SDF file
input_sdf = "AD194827_Del-1-1.sdf"  # Replace with your SDF filename
output_dir = "pdbqt_files_new"
os.makedirs(output_dir, exist_ok=True)

supplier = Chem.SDMolSupplier(input_sdf)
for i, mol in enumerate(supplier):
    if mol is None:
        continue

    # Add hydrogens
    mol = Chem.AddHs(mol)

    # Generate 3D conformer
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.MMFFOptimizeMolecule(mol)

    # Prepare molecule for AutoDock Vina using Meeko
    preparator = MoleculePreparation()
    preparator.prepare(mol)
    pdbqt_string = preparator.write_pdbqt_string()

    # Save PDBQT file
    pdbqt_filename = os.path.join(output_dir, f"molecule_{i}.pdbqt")
    with open(pdbqt_filename, "w") as f:
        f.write(pdbqt_string)

    print(f"Processed and saved {pdbqt_filename}")

print(f"Saved PDBQT files in {output_dir}")
