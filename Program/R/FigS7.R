topDir <- "~/IMIC/table/flux_sum_analysis/"
IMIC <- read.table(paste0(topDir, '14_MAG_media_mets_flux.csv'),header = T, sep = ',')
coco <- read.table(paste0(topDir, 'coco_14_MAG_media_mets_flux.csv'),header = T, sep = ',')
coco_ab1 <- read.table(paste0(topDir, 'coco_ab1_14_MAG_media_mets_flux.csv'),header = T, sep = ',')
micom <- read.table(paste0(topDir, 'MICOM_14_MAG_media_mets_flux.csv'),header = T, sep = ',')
micom_ab1 <- read.table(paste0(topDir, 'MICOM_ab1_14_MAG_media_mets_flux.csv'),header = T, sep = ',')

measured_mets <- read.table(paste0(topDir, 'measured_metabolites.csv'), header = T, sep = ',')

# Define the function
process_flux_data <- function(IMIC, measured_mets) {
  # Add the compartment info
  measured_mets$ID <- paste0(measured_mets$ID, '[e]')
  
  # Find the common IDs 
  common_ids <- intersect(measured_mets$ID, IMIC$out_transport)
  filtered_IMIC <- IMIC[IMIC$out_transport %in% common_ids, ]
  filtered_IMIC$Metabolite_name <- measured_mets$Metabolite.name[match(filtered_IMIC$out_transport, measured_mets$ID)]
  
  # Remove the column of MetaNetX IDs
  filtered_IMIC <- filtered_IMIC[-1]
  
  # Remove the metabolites with flux sum == 0 at all the time points
  filtered_IMIC <- filtered_IMIC[rowSums(filtered_IMIC[, 1:5], na.rm = TRUE) != 0, ]
  
  # Value transform 
  filtered_IMIC[, 1:5] <- filtered_IMIC[, 1:5]/100
  
  # Organic acids
  sugur_organic_acids <- c('D-Mannose', 'D-Glucose', 'D-Galactose', 'D-Fructose', 'L-Lactate', 'D-Lactate', 'Acetate')
  sugars <- filtered_IMIC[filtered_IMIC$Metabolite_name %in% sugur_organic_acids, ]
  colnames(sugars) <- c('20d','40d','60d','90d','180d', 'Sugars and organic acids')
  sugars$`Sugars and organic acids` <- gsub('^[DL]-', '', sugars$`Sugars and organic acids` )
  sugars <- aggregate(. ~ `Sugars and organic acids` , data = sugars, sum)
  
  # Merge the 'D-Mannose', 'D-Glucose', 'D-Galactose' and 'D-Fructose' in sugars$`Sugar and organic acids`
  sugars$`Sugars and organic acids` <- ifelse(sugars$`Sugars and organic acids` %in% c('Mannose', 'Glucose', 'Galactose', 'Fructose'), 'Free sugars', sugars$`Sugars and organic acids`)
  sugars <- aggregate(. ~ `Sugars and organic acids`, data = sugars, sum)
  
  # Amino acids
  amino_acid_names <- setdiff(filtered_IMIC$Metabolite_name, sugur_organic_acids)
  amino_acid <- filtered_IMIC[filtered_IMIC$Metabolite_name %in% amino_acid_names, ]
  colnames(amino_acid) <- c('20d','40d','60d','90d','180d', 'Amino acids')
  amino_acid$`Amino acids` <- gsub('^[DL]-', '', amino_acid$`Amino acids`)
  amino_acid <- aggregate(. ~ `Amino acids`, data = amino_acid, sum)
  
  # Reshape the data for plotting
  sugars_melted <- reshape2::melt(sugars, id.vars = c('Sugars and organic acids'), variable.name = 'Time', value.name = 'Flux')
  amino_melted <- reshape2::melt(amino_acid, id.vars = c('Amino acids'), variable.name = 'Time', value.name = 'Flux')
  
  return(list(sugars_melted = sugars_melted, amino_melted = amino_melted))
}

# Call the function
# IMIC
IMIC_result <- process_flux_data(IMIC, measured_mets)
IMIC_sugar <- IMIC_result$sugars_melted
IMIC_amino <- IMIC_result$amino_melted

# CoCo-GEM
coco_result <- process_flux_data(coco, measured_mets)
coco_sugar <- coco_result$sugars_melted
coco_amino <- coco_result$amino_melted

# CoCo-GEM ab1
coco_ab1_result <- process_flux_data(coco, measured_mets)
coco_ab1_sugar <- coco_ab1_result$sugars_melted
coco_ab1_amino <- coco_ab1_result$amino_melted

# MICOM
micom_result <- process_flux_data(micom, measured_mets)
micom_sugar <- micom_result$sugars_melted
micom_amino <- micom_result$amino_melted

# MICOM ab1
micom_ab1_result <- process_flux_data(micom_ab1, measured_mets)
micom_ab1_sugar <- micom_ab1_result$sugars_melted
micom_ab1_amino <- micom_ab1_result$amino_melted

library(ggplot2)
library(reshape2)

# Plot the data
# Free sugars
a <- ggplot(IMIC_sugar, aes(x = Time, y = Flux, color = `Sugars and organic acids`, group = `Sugars and organic acids`)) +
  geom_line() +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "#A65628", "#FF7F00", "yellow2")) +
  labs(title = 'a IMIC', y = 'Flux-sum (mmol/g-DCW-hr)') +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        panel.border = element_rect(color = 'black', linewidth = 1, fill = NA),
        text = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        plot.title = element_text(face = "bold", size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1))
a

b <- ggplot(micom_ab1_sugar, aes(x = Time, y = Flux, color = `Sugars and organic acids`, group = `Sugars and organic acids`)) +
  geom_line() +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "#A65628", "#FF7F00", "yellow2")) +
  labs(title = 'b MICOM', y = '') +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        panel.border = element_rect(color = 'black', linewidth = 1, fill = NA),
        text = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        plot.title = element_text(face = "bold", size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1))

b

c <- ggplot(coco_ab1_sugar, aes(x = Time, y = Flux, color = `Sugars and organic acids`, group = `Sugars and organic acids`)) +
  geom_line() +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "#A65628", "#FF7F00", "yellow2")) +
  labs(title = 'c CoCo-GEM', y = '') +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        panel.border = element_rect(color = 'black', linewidth = 1, fill = NA),
        text = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        plot.title = element_text(face = "bold", size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1))

c

d <- ggplot(micom_sugar, aes(x = Time, y = Flux, color = `Sugars and organic acids`, group = `Sugars and organic acids`)) +
  geom_line() +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "#A65628", "#FF7F00", "yellow2")) +
  labs(title = 'd MICOM', y = '') +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        panel.border = element_rect(color = 'black', linewidth = 1, fill = NA),
        text = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        plot.title = element_text(face = "bold", size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1))

d

e <- ggplot(coco_sugar, aes(x = Time, y = Flux, color = `Sugars and organic acids`, group = `Sugars and organic acids`)) +
  geom_line() +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "#A65628", "#FF7F00", "yellow2")) +
  labs(title = 'e CoCo-GEM', y = '') +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        panel.border = element_rect(color = 'black', linewidth = 1, fill = NA),
        text = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        plot.title = element_text(face = "bold", size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1))

e

f <- ggplot(IMIC_amino, aes(x = Time, y = Flux, color = `Amino acids`, group = `Amino acids`)) +
  geom_line() +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "orange4", 
    "yellow4", "tan", "#F781BF", "#999999", "#66C2A5", 
    "orange1")) +
  labs(title = 'f IMIC', y = 'Flux-sum (mmol/g-DCW-hr)') +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        panel.border = element_rect(color = 'black', linewidth = 1, fill = NA),
        text = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        plot.title = element_text(face = "bold", size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1))

f

g <- ggplot(micom_ab1_amino, aes(x = Time, y = Flux, color = `Amino acids`, group = `Amino acids`)) +
  geom_line() +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "orange4", 
    "yellow4", "tan", "#F781BF", "#999999", "#66C2A5", 
    "orange1")) +
  labs(title = 'g MICOM', y = '') +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        panel.border = element_rect(color = 'black', linewidth = 1, fill = NA),
        text = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        plot.title = element_text(face = "bold", size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1))

g

h <- ggplot(coco_ab1_amino, aes(x = Time, y = Flux, color = `Amino acids`, group = `Amino acids`)) +
  geom_line() +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "orange4", 
    "yellow4", "tan", "#F781BF", "#999999", "#66C2A5", 
    "orange1")) +
  labs(title = 'h CoCo-GEM', y = '') +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        panel.border = element_rect(color = 'black', linewidth = 1, fill = NA),
        text = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        plot.title = element_text(face = "bold", size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1))

h

i <- ggplot(micom_amino, aes(x = Time, y = Flux, color = `Amino acids`, group = `Amino acids`)) +
  geom_line() +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "orange4", 
    "yellow4", "tan", "#F781BF", "#999999", "#66C2A5", 
    "orange1")) +
  labs(title = 'i MICOM', y = '') +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        panel.border = element_rect(color = 'black', linewidth = 1, fill = NA),
        text = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        plot.title = element_text(face = "bold", size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1))

i

j <- ggplot(coco_amino, aes(x = Time, y = Flux, color = `Amino acids`, group = `Amino acids`)) +
  geom_line() +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "orange4", 
    "yellow4", "tan", "#F781BF", "#999999", "#66C2A5", 
    "orange1")) +
  labs(title = 'j CoCo-GEM', y = '') +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        panel.border = element_rect(color = 'black', linewidth = 1, fill = NA),
        text = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        plot.title = element_text(face = "bold", size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1))

j

library(patchwork)
sugur_plot <- (a + b + c + d + e) + plot_layout(nrow = 1, guides = "collect") & theme(legend.position = 'bottom')
amino_plot <- (f + g + h + i + j) + plot_layout(nrow = 1, guides = "collect") & theme(legend.position = 'bottom')

combined_plot <- (sugur_plot / amino_plot) + plot_layout(ncol = 1)

# Save the combined plot as an SVG file
ggsave(filename = "~/IMIC/Figure/Fig S7.svg", plot = combined_plot, width = 35, height = 25, units = "cm")

# calculate the correlation between predicted flux and measured concentration
han_data <- read.table(paste0(topDir, 'Measured amino acids from Han et al 2020.csv'), header = T, sep = ',')

# Create the data frame with the first column filled
hab_data1 <- data.frame(`Amino acids` = han_data$Amino.acid)
colnames(hab_data1) <- 'Amino acids'

# List of time points and corresponding columns in `han_data`
time_points <- c("20d", "40d", "60d", "90d", "180d")
han_columns <- list(
  c('G1.20d', 'G2.20d', 'G3.20d'),
  c('G1.40d', 'G2.40d', 'G3.40d'),
  c('G1.70d', 'G2.70d', 'G3.70d'),
  c('G1.90d', 'G2.90d', 'G3.90d'),
  c('G1.180d', 'G2.180d', 'G3.180d')
)

# Compute row means dynamically
for (i in seq_along(time_points)) {
  hab_data1[[time_points[i]]] <- rowMeans(han_data[, han_columns[[i]]], na.rm = TRUE)
}

# Reshape to long format
hab_data1 <- reshape2::melt(hab_data1, id.vars = 'Amino acids', variable.name = 'Time', value.name = 'Concentration')

# Merge `hab_data1` with multiple amino acid datasets in a loop
amino_datasets <- list(IMIC_amino, coco_amino, micom_amino, coco_ab1_amino, micom_ab1_amino)
amino_names <- c("IMIC", "CoCo", "MICOM", "CoCo_ab1", "MICOM_ab1")

merged_amino <- lapply(amino_datasets, function(df) merge(df, hab_data1, by = c("Amino acids", "Time"), all = FALSE))
names(merged_amino) <- amino_names

# Function to perform correlation analysis
corr_analysis <- function(aa_data) {
  aa_data %>%
    group_by(`Amino acids`) %>%
    summarise(
      rho = cor.test(Flux, Concentration, method = "spearman")$estimate,
      p_value = cor.test(Flux, Concentration, method = "spearman")$p.value
    ) %>%
    ungroup() %>%
    mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      # now round for presentation
      rho     = round(rho,     5),
      p_value = round(p_value, 5),
      p_adj   = round(p_adj,   4)
    )
}

# Apply correlation analysis to each merged dataset
cor_results_list <- lapply(merged_amino, corr_analysis)

# Combine all results into a single data frame
cor_results <- do.call(cbind, lapply(cor_results_list, function(df) select(df, rho, p_adj)))
colnames(cor_results) <- paste0(rep(amino_names, each = 2), c("_r", "_p"))

# Set row names
rownames(cor_results) <- cor_results_list[[1]]$`Amino acids`

# Save the correlation results to a CSV file
write.csv(cor_results, file = "~/IMIC/table/flux_sum_analysis/correlation_results.csv")


