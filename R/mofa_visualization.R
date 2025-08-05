# mofa_visualization.R

library(MOFA2)
library(ggplot2)
library(ggpubr)
library(ggrepel)

# === Load the trained MOFA model ===
MOFAobject <- load_model("../model/MOFA_model.hdf5")

# === 1. Data Overview ===
plot_data_overview(MOFAobject) + theme(text = element_text(size = 14))

# === 2. Variance Explained ===

## 2A. Variance explained per factor
plot_variance_explained(MOFAobject, max_r2 = 10) +
  theme(text = element_text(size = 14),
        panel.background = element_rect(fill = "white", color = "gray90")) +
  scale_fill_gradientn(
    colors = c("white", "steelblue4", "steelblue4"),
    values = scales::rescale(c(0, 6, 10)),  # Fully qualified call
    limits = c(0, 10),
    name = "Variance (%)"
  )

## 2B. Total variance explained per view
vp <- plot_variance_explained(MOFAobject, plot_total = TRUE)[[2]]
plot_data <- ggplot_build(vp)$data[[1]]
plot_data$x <- factor(c("Mutations", "CNVs", "RNA", "Fusions", "sDSS_asym"),
                      levels = c("Mutations", "CNVs", "RNA", "Fusions", "sDSS_asym")[order(-plot_data$y)])
ggplot(plot_data, aes(x = x, y = y)) +
  geom_bar(stat = "identity", fill = "dodgerblue4", color = "black", width = 0.9) +
  theme(text = element_text(size = 14)) +
  xlab("") + ylab("Variance explained (%)") +
  ggpubr::theme_pubr()

# === 3. Factor 1 – Dot Plots for Key Views ===

## 3A. RNA – Factor 1
data <- plot_weights(MOFAobject, view = "RNA", factor = 1,
                     nfeatures = 4, scale = TRUE, return_data = TRUE)
data <- data[order(data$value), ]
ggplot(data, aes(x = value, y = reorder(feature, value))) +
  geom_point(size = 0.8, color = "gray55") +
  ggrepel::geom_label_repel(label = ifelse(data$value > 0.9, as.character(data$feature), ''),
                            box.padding = 0.35, point.padding = 0.5,
                            color = "dodgerblue4", segment.color = 'gray29', max.overlaps = Inf) +
  coord_cartesian(clip = 'off') + xlim(-1, 1) +
  ylab("Gene") + xlab("Weight (RNA)") +
  ggpubr::theme_pubr() + ggpubr::rremove("y.text") + ggpubr::rremove("y.ticks")

## 3B. sDSS_asym – Factor 1
data <- plot_weights(MOFAobject, view = "sDSS_asym", factor = 1,
                     nfeatures = 1, scale = TRUE, return_data = TRUE)
data <- data[order(data$value), ]
ggplot(data, aes(x = value, y = reorder(feature, value))) +
  geom_point(size = 0.8, color = "gray55") +
  ggrepel::geom_label_repel(label = ifelse(data$value > 0.9, as.character(data$feature), ''),
                            box.padding = 0.35, point.padding = 0.5,
                            color = "dodgerblue4", segment.color = 'gray29') +
  coord_cartesian(clip = 'off') + xlim(-1, 1) +
  ylab("Drug") + xlab("Weight (sDSS_asym)") +
  ggpubr::theme_pubr() + ggpubr::rremove("y.text") + ggpubr::rremove("y.ticks")

## 3C. Mutations – Factor 1
data <- plot_weights(MOFAobject, view = "Mutations", factor = 1,
                     nfeatures = 1, scale = TRUE, return_data = TRUE)
data <- data[order(data$value), ]
ggplot(data, aes(x = value, y = reorder(feature, value))) +
  geom_point(size = 0.8, color = "gray55") +
  ggrepel::geom_label_repel(label = ifelse(data$value > 0.8, as.character(data$feature), ''),
                            box.padding = 0.35, point.padding = 0.5,
                            color = "dodgerblue4", segment.color = 'gray29') +
  coord_cartesian(clip = 'off') + xlim(-1, 1) +
  ylab("Mutation") + xlab("Weight (Mutations)") +
  ggpubr::theme_pubr() + ggpubr::rremove("y.text") + ggpubr::rremove("y.ticks")

# === 4. Factor 3 – RNA and Drug Dot Plots ===

## 4A. RNA – Factor 3
data <- plot_weights(MOFAobject, view = "RNA", factor = 3,
                     nfeatures = 4, scale = TRUE, return_data = TRUE)
data <- data[order(data$value), ]
ggplot(data, aes(x = value, y = reorder(feature, value))) +
  geom_point(size = 0.8, color = "gray55") +
  ggrepel::geom_label_repel(label = ifelse(data$value > 0.9, as.character(data$feature), ''),
                            box.padding = 0.35, point.padding = 0.5,
                            color = "dodgerblue4", segment.color = 'gray29') +
  coord_cartesian(clip = 'off') + xlim(-1, 1) +
  ylab("Gene") + xlab("Weight (RNA – Factor 3)") +
  ggpubr::theme_pubr() + ggpubr::rremove("y.text") + ggpubr::rremove("y.ticks")

## 4B. sDSS_asym – Factor 3
data <- plot_weights(MOFAobject, view = "sDSS_asym", factor = 3,
                     nfeatures = 1, scale = TRUE, return_data = TRUE)
data <- data[order(data$value), ]
ggplot(data, aes(x = value, y = reorder(feature, value))) +
  geom_point(size = 0.8, color = "gray55") +
  ggrepel::geom_label_repel(label = ifelse(data$value > 0.97, as.character(data$feature), ''),
                            box.padding = 0.35, point.padding = 0.5,
                            color = "dodgerblue4", segment.color = 'gray29') +
  coord_cartesian(clip = 'off') + xlim(-1, 1) +
  ylab("Drug") + xlab("Weight (sDSS – Factor 3)") +
  ggpubr::theme_pubr() + ggpubr::rremove("y.text") + ggpubr::rremove("y.ticks")

# === 5. Fusion Features – Factors 2 & 5 ===

## Combined dot plot for both factors
for (f in c(2, 5)) {
  data <- plot_weights(MOFAobject, view = "Fusions", factor = f,
                       nfeatures = 4, scale = TRUE, return_data = TRUE)
  data <- data[order(data$value), ]
  
  p <- ggplot(data, aes(x = value, y = reorder(feature, value))) +
    geom_point(size = 0.8, color = "gray55") +
    ggrepel::geom_label_repel(label = ifelse(data$value > 0.9, as.character(data$feature), ''),
                              box.padding = 0.35, point.padding = 0.5,
                              color = "dodgerblue4", segment.color = 'gray29',
                              max.overlaps = Inf) +
    coord_cartesian(clip = 'off') + xlim(-1, 1) +
    ylab("Fusion") + xlab(paste0("Weight (Fusions – Factor ", f, ")")) +
    ggpubr::theme_pubr() +
    ggpubr::rremove("y.text") + ggpubr::rremove("y.ticks")
  
  print(p)  # <- Add this to display the plot
}

