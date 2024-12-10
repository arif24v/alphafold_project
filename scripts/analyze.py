import os

# Define the directory containing the log files
log_dir = "../library_docking/mass_docking_3d"
output_file = "top_local_search.txt"

# Function to check if a file contains an affinity <= -5.6
def check_affinity(file_path):
    with open(file_path, 'r') as file:
        bool = False
        for line in file:
            if "-----+------------+----------+----------" in line:
                bool = True
                continue
            if bool:
                try:
                    affinity = float(line.split()[1])
                    if affinity <= -6.7:
                        affinities.append(affinity)
                        return True
                except: 
                    continue
                bool = False
    return False

# List to hold the file names that meet the condition
filtered_files = []
affinities = []

# Iterate through all the files in the directory
for file_name in os.listdir(log_dir):    
    if file_name.endswith(('_log.log')): #if file_name.endswith(('_log.log', '_log1.log', '_log2.log')):
        file_path = os.path.join(log_dir, file_name)
        if check_affinity(file_path):
            filtered_files.append(file_name)

# Write the filtered file names to the output file
with open(output_file, 'w') as out_file:
    index = 0
    for file_name in filtered_files:
        out_file.write(f"{file_name}: {affinities[index]} \n")
        index+=1

print(f"Filtered file names have been written to {output_file}")