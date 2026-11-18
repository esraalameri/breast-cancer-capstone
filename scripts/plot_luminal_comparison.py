
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

# Load processed breast cancer data
df = pd.read_csv(
    "data/processed/breast_cancer_expression.csv"
)

# Keep only Luminal A and Luminal B
luminal = df[
    df["PAM50Call_RNAseq"].isin(["LumA", "LumB"])
].copy()

# Genes to compare
genes = ["ESR1", "MKI67", "CCNB1"]

# Create three side-by-side boxplots
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

for ax, gene in zip(axes, genes):
    sns.boxplot(
        data=luminal,
        x="PAM50Call_RNAseq",
        y=gene,
        order=["LumA", "LumB"],
        ax=ax
    )

    ax.set_title(gene)
    ax.set_xlabel("Molecular subtype")
    ax.set_ylabel("Expression (log2 normalised count + 1)")

fig.suptitle(
    "Gene Expression Differences Between Luminal A and Luminal B",
    fontsize=14
)

plt.tight_layout()

# Save figure
Path("results/figures").mkdir(parents=True, exist_ok=True)

plt.savefig(
    "results/figures/luminal_gene_comparison.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print("Combined luminal comparison figure saved successfully!")
