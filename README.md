# BioinformHer R Class Project – scRNA-seq Marker Verification

This repository contains Pavithra's submission for the BioinformHer R Class Project, focused on validating lineage-specific marker expression in a simulated single-cell RNA sequencing (scRNA-seq) dataset.

## 📁 Contents

- `BioinformHer_Project.R` – Clean, commented R script answering all 20 questions using base R
- `BioinformHer_Report.pdf` – Detailed report with biological interpretations and reasoning
- `.RData` – Saved R workspace for reproducibility
- `.Rhistory` – R command history (optional)

## 🧬 Project Summary

The dataset includes 12 single cells manually labeled into Tcell, Bcell, and Myeloid clusters. Each cell is quantified for four key markers:

- **CD3D** – T-cell marker  
- **MS4A1** – B-cell marker  
- **LYZ** – Myeloid marker  
- **MKI67** – Proliferation marker  

Using base R functions (`apply`, `tapply`, `sapply`, `order`, logical indexing), the project verifies whether marker expression patterns support the cluster labels. The analysis includes:

- Total expression and activity labeling  
- Cluster-wise mean comparisons  
- Matrix operations and thresholding  
- Visualization and interpretation  
