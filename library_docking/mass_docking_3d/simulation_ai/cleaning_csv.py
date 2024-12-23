import csv
import os
import re

# File paths
input_csv_path = "combined_output.csv"  # CSV with molecule data
log_files_directory = "./../output_data/localized_search/localized_search_scores/"  # Directory containing .log files
output_csv_path = "new_combined_output.csv"  # Output CSV with appended log data

# Function to extract the best affinity score from a .log file
def extract_best_score(log_file_path):
    try:
        with open(log_file_path, "r") as file:
            for line in file:
                # Match the first line in the mode table with an affinity score
                match = re.match(r"\s*1\s+(-?\d+\.\d+)", line)
                if match:
                    return float(match.group(1))  # Return the best score
    except Exception as e:
        print(f"Error reading {log_file_path}: {e}")
    return None  # Return None if no score is found or an error occurs

# Read the input CSV, process, and write the output
with open(input_csv_path, mode="r") as infile, open(output_csv_path, mode="w", newline="") as outfile:
    csv_reader = csv.reader(infile)
    csv_writer = csv.writer(outfile)

    # Write the header with an additional column for the best score
    header = next(csv_reader)
    header.append("Best Log Score (kcal/mol)")
    csv_writer.writerow(header)

    for row in csv_reader:
        # Get molecule identifier
        molecule_name = row[0]
        log_file_name = f"{molecule_name}.pdbqt_log.log"
        log_file_path = os.path.join(log_files_directory, log_file_name)

        # Extract the best score from the log file
        best_score = extract_best_score(log_file_path)

        # Append the best score to the row
        row.append(best_score if best_score is not None else "N/A")
        csv_writer.writerow(row)

print(f"Data with best scores has been written to {output_csv_path}")
