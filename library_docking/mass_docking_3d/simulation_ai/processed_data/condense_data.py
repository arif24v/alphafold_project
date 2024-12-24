import csv

# File paths
input_csv_path = "new_combined_output.csv"  # Input CSV with unwanted columns
output_csv_path = "condensed_output.csv"  # Output CSV after removing columns

# Columns to delete (zero-indexed)
columns_to_delete = [2, 3, 4, 5]

# Read the input CSV, remove specified columns, and write the output
with open(input_csv_path, mode="r") as infile, open(output_csv_path, mode="w", newline="") as outfile:
    csv_reader = csv.reader(infile)
    csv_writer = csv.writer(outfile)

    for row in csv_reader:
        # Keep only the columns not in `columns_to_delete`
        condensed_row = [col for idx, col in enumerate(row) if idx not in columns_to_delete]
        csv_writer.writerow(condensed_row)

print(f"Condensed data has been written to {output_csv_path}")
