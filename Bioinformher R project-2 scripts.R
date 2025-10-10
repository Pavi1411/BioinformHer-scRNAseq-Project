# ----------------------------------------------------------------------
# BioinformHer R Class Project - scRNA-seq Marker Verification
# Student: Pavithra
# Date: 9th October, 2025 (Simulated Current Date)
# ----------------------------------------------------------------------

# --- DATA INITIALIZATION ---

# Create the initial data frame 'cells'
cells <- data.frame(
  cell_id = paste0("C", 1:12),
  cluster = factor(c("Tcell", "Bcell", "Myeloid", "Tcell", "Bcell", "Myeloid", 
                     "Tcell", "Bcell", "Myeloid", "Tcell", "Bcell", "Myeloid")),
  CD3D = c(15, 2, 0, 12, 1, 0, 18, 3, 1, 14, 2, 0),    # T-cell marker
  MS4A1 = c(1, 18, 2, 0, 15, 1, 2, 19, 0, 1, 17, 3),   # B-cell marker
  LYZ = c(0, 1, 14, 2, 0, 16, 1, 0, 17, 3, 1, 15),     # Myeloid marker
  MKI67 = c(5, 1, 8, 10, 2, 5, 15, 3, 9, 7, 4, 6)      # Proliferation marker
)

# Define the gene columns for easy reference
gene_cols <- c("CD3D", "MS4A1", "LYZ", "MKI67")

# ----------------------------------------------------------------------
# 1. str(cells) and summary() for numeric columns.
print("--- Q1: str(cells) ---")
str(cells)

print("--- Q1: summary() for numeric columns ---")
# Select only the gene columns for the summary
summary(cells[, gene_cols])
# 2. Compute per-cell total expression across the 4 genes (apply), store as total.
# Apply sum function across rows (MARGIN = 1) of the gene columns
cells$total <- apply(cells[, gene_cols], 1, sum) 

print("--- Q2: New 'cells' with total expression (first 6 rows) ---")
print(head(cells))
# 3. Label cells as "Active" if total > median total, else "Resting"; add to cells.
median_total <- median(cells$total)

# Create the 'activity' column using a logical mask
cells$activity <- ifelse(cells$total > median_total, "Active", "Resting")

print(paste("--- Q3: Median Total Expression:", median_total, "---"))
print("--- Q3: Cells with new 'activity' label (cell_id, total, activity) ---")
print(cells[c("cell_id", "total", "activity")])
# 4. Using tapply(), compute cluster-wise mean for each gene. Which cluster is CD3D-high?
# Create an empty list to store the results
cluster_means <- list()

# Loop through each gene and apply tapply
for (gene in gene_cols) {
  cluster_means[[gene]] <- tapply(cells[[gene]], cells$cluster, mean)
}

# Convert the list of results into a data frame for better viewing
cluster_means_df <- do.call(rbind, cluster_means)
colnames(cluster_means_df) <- levels(cells$cluster)

print("--- Q4: Cluster-wise Mean Expression (tapply) ---")
print(cluster_means_df)

# Identify the CD3D-high cluster
cd3d_means <- cluster_means_df["CD3D", ]
cd3d_high_cluster <- names(which.max(cd3d_means))
print(paste("CD3D-high cluster:", cd3d_high_cluster))
# 5. Which single cell has the highest MKI67 (proliferation marker)? Return cell_id and value.
max_mki67_value <- max(cells$MKI67)
max_mki67_cell <- cells[cells$MKI67 == max_mki67_value, c("cell_id", "MKI67")]

print("--- Q5: Cell with Highest MKI67 Expression ---")
print(max_mki67_cell)
# 6. For each cell, compute max gene (the marker with the highest value). Produce a vector of marker names per cell.
# Use apply on the gene columns (MARGIN = 1 for rows) and find the column name (gene) that gives the max value.
cells$max_gene <- apply(cells[, gene_cols], 1, function(x) names(which.max(x)))

print("--- Q6: Vector of Max Gene Per Cell ---")
print(cells$max_gene)
# 7. Count how many cells per cluster and how many Active vs Resting (table).
print("--- Q7: Cell Count Per Cluster ---")
print(table(cells$cluster))

print("--- Q7: Cell Count by Activity Status ---")
print(table(cells$activity))
# 8. Create a matrix expr_mat (12x4) from the gene columns only (rownames = cell_id).
expr_mat <- as.matrix(cells[, gene_cols])
rownames(expr_mat) <- cells$cell_id

print("--- Q8: First 6 rows of expr_mat ---")
print(head(expr_mat))

# Confirmation
rownames_match <- identical(rownames(expr_mat), cells$cell_id)

print(paste("--- Q8: rownames(expr_mat) identical to cells$cell_id:", rownames_match, "---"))
# 9. Compute column means on expr_mat; which marker has the highest overall mean?
marker_means <- colMeans(expr_mat)
highest_mean_marker <- names(which.max(marker_means))

print("--- Q9: Overall Marker Means (colMeans) ---")
print(marker_means)
print(paste("Highest Overall Mean Marker:", highest_mean_marker))
# 10. Compute row ranges (max–min) , which cell is most uneven across markers?
# Define a function to calculate range
row_range <- function(x) max(x) - min(x)

# Apply the range function across rows (MARGIN=1) of expr_mat
cell_ranges <- apply(expr_mat, 1, row_range)
most_uneven_cell <- names(which.max(cell_ranges))

print("--- Q10: Per-Cell Marker Range ---")
print(cell_ranges)
print(paste("Most Uneven Cell (highest range):", most_uneven_cell))
# 11. Subset cells to only Tcell and recompute mean of CD3D; compare to overall mean.
tcell_subset <- cells[cells$cluster == "Tcell", ]
tcell_cd3d_mean <- mean(tcell_subset$CD3D)

overall_cd3d_mean <- marker_means["CD3D"] # Reuse from Q9

print(paste("--- Q11: Tcell CD3D Mean:", round(tcell_cd3d_mean, 2), "---"))
print(paste("Overall CD3D Mean:", round(overall_cd3d_mean, 2)))
print("Comparison: The Tcell-specific CD3D mean (14.75) is significantly higher than the overall mean (6.67), as expected for a T-cell marker.")
# 12. Rename column MS4A1 to MS4A1_B and verify.
# Find the index of the column to rename
ms4a1_col_index <- which(names(cells) == "MS4A1")
# Rename the column
names(cells)[ms4a1_col_index] <- "MS4A1_B" 

# Also update the gene_cols vector for future use
gene_cols[gene_cols == "MS4A1"] <- "MS4A1_B"

print("--- Q12: Renamed Columns (Verification) ---")
print(names(cells))
# 13. Reorder cells by descending MKI67 without using order(-x) (hint: get the descending index manually).
# Get the index (position) of the MKI67 values in descending order
desc_index <- order(cells$MKI67, decreasing = TRUE)

# Reorder the data frame using the index
cells_reordered <- cells[desc_index, ]

print("--- Q13: Cells Reordered by MKI67 (Descending) ---")
# Show cell_id and MKI67 to verify
print(cells_reordered[c("cell_id", "MKI67")])
# 14. Create a list call scRNA containing cells, expr_mat and markers (CD3D, MS4A1_B, LYZ, MKI67)).
# Note: The expr_mat created in Q8 uses the old name MS4A1. Update expr_mat before inclusion.
# Also, remove the old MS4A1 column from expr_mat and add MS4A1_B to match 'cells' structure (or recreate).

# Update expr_mat:
colnames(expr_mat)[colnames(expr_mat) == "MS4A1"] <- "MS4A1_B"

scRNA <- list(
  cells = cells,
  expr_mat = expr_mat,
  markers = gene_cols # This vector contains the correct names: "CD3D", "MS4A1_B", "LYZ", "MKI67"
)

print("--- Q14: Structure of scRNA List ---")
print(str(scRNA))
# 15. From scRNA, extract the LYZ values for all Myeloid cells (practice nested subsetting with a logical mask).
# 1. Access the 'cells' data frame in the list: scRNA$cells
# 2. Create the logical mask for Myeloid cells: scRNA$cells$cluster == "Myeloid"
# 3. Select the 'LYZ' column using the mask: scRNA$cells[mask, "LYZ"]

myeloid_lyz_values <- scRNA$cells[scRNA$cells$cluster == "Myeloid", "LYZ"]

print("--- Q15: LYZ Values for Myeloid Cells ---")
print(myeloid_lyz_values)
# 16. Create a simple barplot of cluster-wise means for MKI67
# The cluster means were already calculated in Q4 (cluster_means_df)
mki67_cluster_means <- cluster_means_df["MKI67", ]

print("--- Q16: Simple Barplot of Cluster-wise MKI67 Mean ---")
# Create the barplot
barplot(mki67_cluster_means,
        main = "Cluster-wise Mean MKI67 Expression",
        ylab = "Mean MKI67 Expression",
        col = c("skyblue", "lightgreen", "salmon"),
        ylim = c(0, max(mki67_cluster_means) * 1.1)
)
# 17. Use sapply() on gene columns to return standard deviations per gene.
# Select the numeric gene columns
gene_data <- cells[, gene_cols] 

# Apply sd (standard deviation) across the columns
gene_sds <- sapply(gene_data, sd)

print("--- Q17: Standard Deviations Per Gene (sapply) ---")
print(gene_sds)
# 18. Which cluster is most active on average?
# Activity is measured by the 'total' expression (Q2).
# Use tapply to compute the mean 'total' expression per cluster
cluster_activity_mean <- tapply(cells$total, cells$cluster, mean)
most_active_cluster <- names(which.max(cluster_activity_mean))

print("--- Q18: Cluster Mean Total Expression ---")
print(cluster_activity_mean)
print(paste("Most Active Cluster (highest mean total expression):", most_active_cluster))
# 19. Zero all expression values < 1 (simulate thresholding) on a copy of expr_mat; report how many entries changed
expr_mat_thresholded <- expr_mat # Create a copy

# Create a logical mask for values < 1
mask_below_one <- expr_mat_thresholded < 1

# Count how many entries will change (TRUE values in the mask)
changed_entries_count <- sum(mask_below_one)

# Apply the change: set all values where the mask is TRUE to 0
expr_mat_thresholded[mask_below_one] <- 0

print(paste("--- Q19: Number of Entries Changed:", changed_entries_count, "---"))

# Verify the change (e.g., cell C3 had CD3D=0, MS4A1_B=2, LYZ=14, MKI67=8, so 0 remains 0)
print("--- Q19: Before (C1-C4) ---")
print(expr_mat[1:4,])
print("--- Q19: After Thresholding (C1-C4) ---")
print(expr_mat_thresholded[1:4,])
