# 🧬 Lung Cancer Immuno-Profiling & Survival Prediction

**Author:** Mohamed Aymen Lahmer  

## 🎯 Project Overview
Immunotherapy has revolutionized Non-Small Cell Lung Cancer (NSCLC) treatment. However, response rates vary. This project aims to identifying "Hot Tumors" (responsive) vs "Cold Tumors" using a custom transcriptomic signature.

## 🛠️ Methodology (The "Lahmer Pipeline")
1. **Data Mining:** Automated extraction of TCGA-LUAD data using `TCGAbiolinks`.
2. **Preprocessing:** Normalization (Log2) and Gene Symbol mapping.
3. **Feature Selection:** Focus on Cytotoxic T-Cell markers (*CD8A, GZMB, PRF1*).
4. **Clinical Analysis:** Correlation between the "Immune Score" and patient survival.

## 📊 Key Results
*(Analysis performed by M.A. Lahmer)*<img width="664" height="664" alt="impact of immunity on survival" src="https://github.com/user-attachments/assets/34540c63-84e9-4206-bb30-ebdc8887bf4b" />


## 🚀 How to Run
```r
# Install dependencies
BiocManager::install("TCGAbiolinks")
# Run the analysis
source("AYMEN.R")
