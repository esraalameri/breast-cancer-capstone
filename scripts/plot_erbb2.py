
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

# Plot ERBB2 expression
plt.figure(figsize=(9, 6))

sns.boxplot(
    data=df,
    x="PAM50Call_RNAseq",
    y="ERBB2",
    order=order
)

plt.title("ERBB2 Expression Across Breast Cancer Molecular Subtypes")
plt.xlabel("PAM50 Molecular Subtype")
plt.ylabel("ERBB2 Expression (log2 normalised count + 1)")

plt.tight_layout()

plt.savefig(
    "results/figures/ERBB2_boxplot.png",
    dpi=300
)

plt.show()

# Calculate median expression
print("\nMedian ERBB2 expression by subtype:")
print(
    df.groupby("PAM50Call_RNAseq")["ERBB2"]
    .median()
    .sort_values(ascending=False)
)
