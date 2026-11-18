
import pandas as pd
from scipy.stats import kruskal

# Load processed dataset
df = pd.read_csv(
    "data/processed/breast_cancer_expression.csv"
)

# Define molecular subtypes
subtypes = ["LumA", "LumB", "Her2", "Basal", "Normal"]

# Collect ESR1 expression values for each subtype
groups = [
    df.loc[df["PAM50Call_RNAseq"] == subtype, "ESR1"]
    .dropna()
    for subtype in subtypes
]

# Perform Kruskal-Wallis test
statistic, p_value = kruskal(*groups)

print("Kruskal-Wallis Test: ESR1 Expression")
print(f"H statistic: {statistic:.4f}")
print(f"P-value: {p_value:.4e}")

if p_value < 0.05:
    print("Result: Statistically significant difference detected.")
else:
    print("Result: No statistically significant difference detected.")


from itertools import combinations
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from pathlib import Path

# Compare every pair of PAM50 subtypes
results = []

for subtype1, subtype2 in combinations(subtypes, 2):

    group1 = df.loc[
        df["PAM50Call_RNAseq"] == subtype1, "ESR1"
    ].dropna()

    group2 = df.loc[
        df["PAM50Call_RNAseq"] == subtype2, "ESR1"
    ].dropna()

    statistic, p_value = mannwhitneyu(
        group1,
        group2,
        alternative="two-sided"
    )

    results.append({
        "Subtype_1": subtype1,
        "Subtype_2": subtype2,
        "P_value": p_value
    })

# Convert results into a table
results_df = pd.DataFrame(results)

# Correct for multiple comparisons
results_df["Adjusted_P_value"] = multipletests(
    results_df["P_value"],
    method="holm"
)[1]

# Save the statistical results
Path("results").mkdir(exist_ok=True)

results_df.to_csv(
    "results/esr1_pairwise_tests.csv",
    index=False
)

print("\nPairwise ESR1 comparisons:")
print(results_df.to_string(index=False))
