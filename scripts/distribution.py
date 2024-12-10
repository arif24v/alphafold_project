import os
import matplotlib.pyplot as plt

# Define the directory containing the log files
log_dir = "../library_docking/mass_docking_3d"

def check_affinity(file_path):
    """Extract the top affinity after the header."""
    with open(file_path, 'r') as file:
        for line in file:
            if "-----+------------+----------+----------" in line:
                try:
                    next_line = next(file)
                    affinity = float(next_line.split()[1])
                    return affinity
                except (StopIteration, ValueError, IndexError):
                    break
    return None

def process_files(log_dir):
    """Process all files in the directory and collect top binding affinities."""
    top_affinities = []

    for file_name in os.listdir(log_dir):
        if file_name.endswith('_log1.log'):
            file_path = os.path.join(log_dir, file_name)
            top_affinity = check_affinity(file_path)
            if top_affinity is not None:
                top_affinities.append(top_affinity)

    return top_affinities

def plot_histogram(data, output_path="breadth_3d_search.png"):
    """Plot a histogram of the top binding affinities."""
    plt.figure(figsize=(10, 6))
    plt.hist(data, bins=30, edgecolor='k', alpha=0.7)
    plt.title('Distribution of Top Binding Affinities')
    plt.xlabel('Binding Affinity')
    plt.ylabel('Frequency')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.savefig(output_path)
    plt.show()

# Main script
if __name__ == "__main__":
    top_affinities = process_files(log_dir)

    # Plot the histogram
    plot_histogram(top_affinities)
