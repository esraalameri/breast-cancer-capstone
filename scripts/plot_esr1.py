
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

# Load our processed dataset
df = pd.read_csv(
    "data/processed/breast_cancer_expression.csv"
)

# Create a folder for figures
Path("results/figures").mkdir(parents=True, exist_ok=True)

# Define the subtype order
order = ["LumA", "LumB", "Her2", "Basal", "Normal"]

# Create the boxplot
plt.figure(figsize=(9, 6))

sns.boxplot(
    data=df,
    x="PAM50Call_RNAseq",
    y="ESR1",
    order=order
)

plt.title("ESR1 Expression Across Breast Cancer Molecular Subtypes")
plt.xlabel("PAM50 Molecular Subtype")
plt.ylabel("ESR1 Expression (log2 normalised count + 1)")

plt.tight_layout()

# Save the figure
plt.savefig(
    "results/figures/esr1_boxplot.png",
    dpi=300
)

plt.show()

print("ESR1 boxplot saved successfully!")
