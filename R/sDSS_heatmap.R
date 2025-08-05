# sDSS_heatmap.R
# Drug sensitivity heatmap using sDSS_asym scores

library(ComplexHeatmap)
library(circlize)

# Define data path
sDSS_path <- file.path("../data", "sDSS_asym.csv")

# Load drug sensitivity matrix
sDSS <- read.csv(sDSS_path, row.names = 1)

# Clean invalid values: keep NA, remove NaN and Inf
sDSS[is.nan(as.matrix(sDSS)) | is.infinite(as.matrix(sDSS))] <- NA

# Convert to matrix
sDSS_mat <- as.matrix(sDSS)

# Define color scale and NA color
heatmap_colors <- colorRamp2(c(-10, 5, 30), c("steelblue4", "white", "tomato3"))

# Draw the heatmap interactively
Heatmap(
  t(sDSS_mat),
  col = heatmap_colors,
  na_col = "grey90",  # visualize NA as grey
  name = "sDSS_asym",
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_gp = grid::gpar(fontsize = 10),
  row_names_gp = grid::gpar(fontsize = 10),
  clustering_distance_columns = "spearman",
  clustering_distance_rows = "pearson",
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2"
)
