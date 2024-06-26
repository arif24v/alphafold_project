USE THIS COMMAND TO SIMULATE WITH AUTODOCK FROM THE COMMAND LINE:

/Users/arifv/Desktop/GitHub/alphafold_project/vina --config [CONFIG.CONF] --log[OUTPUT.TXT]

CONF FILE SHOULD BE OF THE FORMAT:
receptor = [receptor].pdbqt
ligand = [ligand].pdbqt
center_x = 8.08
center_y = -4.85
center_z = -3.92
size_x = 48.44
size_y = 47.68
size_z = 40.62

Chemical Sketch Tool (download as pdf): https://www.rcsb.org/docs/search-and-browse/advanced-search/chemical-sketch-tool

Convert SDF to PDB: https://www.cheminfo.org/Chemistry/Cheminformatics/FormatConverter/index.html

CIF to PDB: https://rna-tools.online/tools/topdb/66d61e19

Goals: 

Modify hCAR1 (348 aa) minimally to increase binding affinity for TCPOBOP
	Find proper hCAR PDB

Modify CITCO to become a better activator of hCAR
	Use chemical draw to modify functional groups
