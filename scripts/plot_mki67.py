
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

# Load processed dataset
df = pd.read_csv(
    "data/processed/breast_cancer_expression.csv"
)

# Create output folder
Path("results/figures").mkdir(parents=True, exist_ok=True)

# Define subtype order
order = ["LumA", "LumB", "Her2", "Basal", "Normal"]

# Plot MKI67 expression
plt.figure(figsize=(9, 6))

sns.boxplot(
    data=df,
    x="PAM50Call_RNAseq",
    y="MKI67",
    order=order
)

plt.title("MKI67 Expression Across Breast Cancer Molecular Subtypes")
plt.xlabel("PAM50 Molecular Subtype")
plt.ylabel("MKI67 Expression (log2 normalised count + 1)")

plt.tight_layout()

plt.savefig(
    "results/figures/mki67_boxplot.png",
    dpi=300
)

plt.show()

# Calculate median expression
print("\nMedian MKI67 expression by subtype:")
print(
    df.groupby("PAM50Call_RNAseq")["MKI67"]
    .median()
    .sort_values(ascending=False)
)
