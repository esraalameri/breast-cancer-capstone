
import pandas as pd
from scipy.stats import mannwhitneyu

# Load processed TCGA-BRCA data
df = pd.read_csv(
    "data/processed/breast_cancer_expression.csv"
)

# Extract MKI67 expression for Luminal A and Luminal B
lumA = df.loc[
    df["PAM50Call_RNAseq"] == "LumA", "MKI67"
].dropna()

lumB = df.loc[
    df["PAM50Call_RNAseq"] == "LumB", "MKI67"
].dropna()

# Perform two-sided Mann-Whitney U test
statistic, p_value = mannwhitneyu(
    lumA,
    lumB,
    alternative="two-sided"
)

print("MKI67: Luminal A vs Luminal B")
print(f"Luminal A median: {lumA.median():.4f}")
print(f"Luminal B median: {lumB.median():.4f}")
print(f"U statistic: {statistic:.4f}")
print(f"P-value: {p_value:.4e}")

