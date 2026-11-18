import pandas as pd
from pathlib import Path

# Define file locations
raw = Path("data/raw")
output = Path("data/processed")
output.mkdir(parents=True, exist_ok=True)

# Load the clinical dataset
clinical = pd.read_csv(
    raw / "TCGA.BRCA.sampleMap_BRCA_clinicalMatrix",
    sep="\t"
)

# Load the gene expression dataset
expression = pd.read_csv(
    raw / "TCGA.BRCA.sampleMap_HiSeqV2.gz",
    sep="\t",
    compression="gzip"
)

# Select our eight genes of interest
genes = [
    "ESR1", "PGR", "ERBB2", "MKI67",
    "CCNB1", "BRCA1", "BRCA2", "TP53"
]

expression = expression[expression["sample"].isin(genes)]

print("Genes found:", expression["sample"].tolist())

# Reshape so each row represents a sample
expression = expression.set_index("sample").T
expression.index.name = "sampleID"
expression = expression.reset_index()

# Keep only sample IDs and PAM50 subtypes
clinical = clinical[["sampleID", "PAM50Call_RNAseq"]]

# Match expression data with clinical annotations
merged = expression.merge(
    clinical,
    on="sampleID",
    how="inner",
    validate="one_to_one"
)

# Remove samples without valid PAM50 annotations
subtypes = ["LumA", "LumB", "Her2", "Basal", "Normal"]

merged = merged[
    merged["PAM50Call_RNAseq"].isin(subtypes)
].copy()

# Save the processed dataset
merged.to_csv(
    output / "breast_cancer_expression.csv",
    index=False
)

print("\nFinal dataset dimensions:", merged.shape)
print("\nSubtype distribution:")
print(merged["PAM50Call_RNAseq"].value_counts())
print("\nSaved to data/processed/breast_cancer_expression.csv")