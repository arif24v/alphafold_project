import os

# Source and destination directories
source_dir = 'pdb'
dest_dir = 'pdbqt'

# Ensure destination directory exists
os.makedirs(dest_dir, exist_ok=True)

# Iterate through files in the source directory
for filename in os.listdir(source_dir):
    if os.path.isfile(os.path.join(source_dir, filename)):
        # Create empty file in the destination directory with .pdbqt extension
        dest_filename = os.path.join(dest_dir, f"{os.path.splitext(filename)[0]}.pdbqt")
        open(dest_filename, 'a').close()  # Create empty file

print("Empty .pdbqt files created in the destination directory.")
