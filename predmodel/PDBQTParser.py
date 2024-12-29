import os

class PDBQTParser:
    def __init__(self, file_path):
        self.file_path = file_path
        self.results = self.parse_pdbqt()
    
    def parse_pdbqt(self):
        results = []
        vina_result = None
        current_model = None

        with open(self.file_path, 'r') as file:
            for line in file:
                if line.startswith("REMARK VINA RESULT"):
                    vina_result = float(line.split()[3]) 
                elif line.startswith("MODEL"):
                    current_model = {
                        "model_id": int(line.split()[1]),
                        "atoms": []
                    }
                elif line.startswith("ATOM"):
                    parts = line.split()
                    atom_type = parts[2]
                    x = float(parts[6])
                    y = float(parts[7])
                    z = float(parts[8])
                    current_model["atoms"].append([atom_type, x, y, z])
                elif line.startswith("ENDMDL"):
                    molecule_name = os.path.basename(self.file_path)
                    molecule_name = molecule_name.split('_out.pdbqt')[0] 
                    results.append([molecule_name, current_model["model_id"], vina_result, current_model["atoms"]])
                    current_model = None

        return results

    def print(self):
        for record in self.results:
            file, model, vina_result, atoms = record
            print(f"\nFile: {file}, Model: {model}, Vina Result: {vina_result}")
            print(f"{'Atom Type':<10}{'X (Å)':<10}{'Y (Å)':<10}{'Z (Å)':<10}")
            for atom in atoms:
                atom_id, atom_type, x, y, z = atom
                print(f"{atom_type:<10}{x:<10.3f}{y:<10.3f}{z:<10.3f}")


