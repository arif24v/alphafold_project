import os

# Specify the paths to the two text files
file1_path = "pdbqt_list.txt"  # Replace with your file path
file2_path = "sdf_list.txt"  # Replace with your file path

# Function to extract base names without extensions
def get_base_names(file_path):
    with open(file_path, 'r') as f:
        return set(os.path.splitext(line.strip())[0] for line in f if line.strip())

# Get base names from both files
file1_base_names = get_base_names(file1_path)
file2_base_names = get_base_names(file2_path)

# Find overlaps and unique entries
overlap = file1_base_names & file2_base_names
unique_to_file1 = file1_base_names - file2_base_names
unique_to_file2 = file2_base_names - file1_base_names

# Print the results
print(f"Number of overlapping file names (base names): {len(overlap)}")
print(f"Number of unique file names in {file1_path} (base names): {len(unique_to_file1)}")
print(f"Number of unique file names in {file2_path} (base names): {len(unique_to_file2)}")

# Optionally, save the results to output files
with open("overlap.txt", 'w') as overlap_file:
    overlap_file.write("\n".join(sorted(overlap)))

with open("unique_to_pdbqt.txt", 'w') as unique1_file:
    unique1_file.write("\n".join(sorted(unique_to_file1)))

with open("unique_to_sdf.txt", 'w') as unique2_file:
    unique2_file.write("\n".join(sorted(unique_to_file2)))

print("Results saved")
