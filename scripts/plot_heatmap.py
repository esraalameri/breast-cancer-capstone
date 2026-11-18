
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

# Load the processed dataset
df = pd.read_csv("data/processed/breast_cancer_expression.csv")

# Select genes and subtype order
genes = ["ESR1", "PGR", "ERBB2", "MKI67",
         "CCNB1", "BRCA1", "BRCA2", "TP53"]

subtypes = ["LumA", "LumB", "Her2", "Basal", "Normal"]

# Calculate median expression for each gene in each subtype
median_expression = (
    df.groupby("PAM50Call_RNAseq")[genes]
    .median()
    .reindex(subtypes)
)

# Standardise each gene across subtypes (row-wise z-score)
heatmap_data = median_expression.T
heatmap_data = heatmap_data.sub(
    heatmap_data.mean(axis=1), axis=0
).div(heatmap_data.std(axis=1), axis=0)

# Create the heatmap
plt.figure(figsize=(9, 7))

sns.heatmap(
    heatmap_data,
    annot=True,
    fmt=".2f",
    cmap="RdBu_r",
    center=0,
    linewidths=0.5,
    cbar_kws={"label": "Relative expression (z-score)"}
)

plt.title("Gene Expression Patterns Across PAM50 Subtypes")
plt.xlabel("Breast Cancer Molecular Subtype")
plt.ylabel("Gene")

plt.tight_layout()

# Save the figure
Path("results/figures").mkdir(parents=True, exist_ok=True)
plt.savefig("results/figures/eight_gene_heatmap.png", dpi=300)

plt.show()

print("Eight-gene heatmap saved successfully!")
