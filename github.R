# ==============================================================================
# PROJECT: TRANSCRIPTOMIC ANALYSIS & IMMUNO-ONCOLOGY (NSCLC)
# AUTHOR:  Mohamed Aymen Lahmer ( Bioinformatics, Bordeaux)
# DATE:    December 2025
# 
# DESCRIPTION: 
# Custom pipeline developed by M.A. Lahmer to identify immunotherapy responders
# in Lung Adenocarcinoma (LUAD) using a cytotoxic CD8+ T-cell signature.
# "From Raw Data to Biomarkers."
# ==============================================================================

# --- 1. ENVIRONMENT INITIALIZATION ---
print(">>> [Lahmer-Pipeline] Loading bioinformatics libraries...")

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
})

# --- 2. DATA MINING (GDC Portal) ---
# Querying the Genomic Data Commons for Lung Adenocarcinoma (LUAD)
# Workflow: STAR - Counts (Raw expression data)

query_lahmer <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts"
)

# NOTE (M.A. Lahmer): Prototyping on 20 samples for efficiency.
# Scalable to the full cohort (500+ patients) for production.
query_lahmer$results[[1]] <- query_lahmer$results[[1]][1:20,]

print(">>> [Lahmer-Pipeline] Downloading TCGA data...")
GDCdownload(query_lahmer)

# preparing the main data object
data_aymen <- GDCprepare(query_lahmer)

# --- 3. DATA ENGINEERING & PREPROCESSING ---
print(">>> [Lahmer-Pipeline] Cleaning and Normalizing sequences...")

raw_matrix <- assay(data_aymen)
gene_info <- rowData(data_aymen)

# Mapping Ensembl IDs to Gene Symbols
rownames(raw_matrix) <- make.unique(gene_info$gene_name)

# Log2 Transformation for variance stabilization (Bioinfo Standard)
matrix_mohamed_log <- log2(raw_matrix + 1)

# --- 4. CUSTOM IMMUNE SIGNATURE DEFINITION ---
# Selecting specific markers for Cytotoxic Activity (CD8+ T-Cells).
# Goal: Detect "Hot Tumors" vs "Cold Tumors".

genes_target_lahmer <- c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "IFNG")

# Intersection check
genes_detected <- intersect(genes_target_lahmer, rownames(matrix_mohamed_log))
signature_matrix <- matrix_mohamed_log[genes_detected, ]

# --- VISUALIZATION 1: Heatmap ---
# Visualizing the immune landscape across patients
heatmap(as.matrix(signature_matrix),
        scale = "row",
        col = cm.colors(256),
        main = "Immune Signature - Analysis by M.A. Lahmer", # Branding
        xlab = "Patients (TCGA-LUAD)",
        ylab = "Cytotoxic Markers")

# --- 5. CLINICAL CORRELATION & SURVIVAL ANALYSIS ---

# Calculating the "Lahmer Score": Mean expression of cytotoxic genes
score_lahmer <- colMeans(signature_matrix)

# Retrieving Clinical Data
clinical_data <- colData(data_aymen)

df_aymen <- data.frame(
  Patient = rownames(clinical_data),
  Status = clinical_data$vital_status,
  Immune_Score = score_lahmer
)

# Removing NA values for rigorous statistics
df_aymen <- df_aymen[!is.na(df_aymen$Status), ]

# --- VISUALIZATION 2: Clinical Boxplot ---
boxplot(Immune_Score ~ Status, data = df_aymen,
        col = c("#ff6b6b", "#4ecdc4"),
        main = "Impact of Immunity on Survival (Lahmer Pipeline)",
        ylab = "Infiltration Score (CD8+)",
        xlab = "Vital Status")

# Statistical Testing (Wilcoxon Rank Sum Test)
test_result <- wilcox.test(Immune_Score ~ Status, data = df_aymen)

print(">>> [Lahmer-Pipeline] STATISTICAL TEST RESULTS:")
print(test_result)

if(test_result$p.value < 0.05){
  print(">>> SUCCESS! The signature significantly predicts survival.")
} else {
  print(">>> Trend observed. Increasing sample size is recommended for significance.")
}

print(">>> [Lahmer-Pipeline] Analysis complete. Ready for deployment.")