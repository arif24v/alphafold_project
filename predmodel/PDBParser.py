class PDBParser:
    def __init__(self, file_path):
        self.file_path = file_path
        self.atom_data = []

    def parse(self):
        with open(self.file_path, 'r') as pdb_file:
            for line in pdb_file:
                if line.startswith("ATOM"):
                    element = line[77:].strip()

                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())

                    self.atom_data.append([element, x, y, z])
        return self.atom_data

    def get_atom_data(self):
        return self.atom_data
    def print(self):
        for atom in self.atom_data:
            print(atom)