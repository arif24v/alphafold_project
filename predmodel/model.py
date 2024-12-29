import os
import math
import matplotlib.pyplot as plt
from PDBQTParser import PDBQTParser
from PDBParser import PDBParser
from AtomicCalculator import AtomicCalculator

calculator = AtomicCalculator()
directory_path = 'output_data/breadth_search/breadth_search_poses/'
protein = PDBParser('hCAR1.pdb')
patoms = protein.parse()

gibs_pot_all = []

for filename in os.listdir(directory_path):
    if filename.endswith('.pdbqt'):
        molecule_path = os.path.join(directory_path, filename)
        molecule = PDBQTParser(molecule_path)
        atoms1 = molecule.parse_pdbqt()
        
        sum1 = 0
        gibs_pot = []
        for model in atoms1:
            for atom in model[3]:
                for patom in patoms:
                    sum1 += AtomicCalculator.lennard_jones_potential(
                        self=AtomicCalculator, particle1=patom, particle2=atom
                    )
            gibs_pot_all.append([model[2], sum1])
            sum1 = 0
        

x = [item[0] for item in gibs_pot]
y = [math.log10((item[1])) for item in gibs_pot]

plt.scatter(x, y, color='blue')

plt.xlabel('Model Identifier')
plt.ylabel('Gibbs Potential')
plt.title('Gibbs Potential Points')

plt.grid(True)
plt.show()

