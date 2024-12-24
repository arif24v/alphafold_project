import pandas as pd
import matplotlib.pyplot as plt

# File path to the CSV
csv_path = "condensed_output.csv"

# Load the CSV file into a DataFrame
data = pd.read_csv(csv_path)

# Extract the AI Confidence column
ai_confidences = data["Local Score"]

# Plot the distribution
plt.figure(figsize=(10, 6))
plt.hist(ai_confidences, bins=20, edgecolor='black', alpha=0.7)
plt.title("Distribution of Localized Search Scores", fontsize=14)
plt.xlabel("Binding Affinity (kcal/mol)", fontsize=12)
plt.ylabel("Frequency", fontsize=12)
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.show()
