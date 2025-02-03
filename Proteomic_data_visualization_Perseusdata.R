
####pre-processing, cleaning and visualization, differential analysis and graphs 
#of protein intensity file obtained from maxQuant software

install.packages("pheatmap")
install.packages("FactoMineR")
install.packages("plotly")
install.packages("limma")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db")  # For human data; use appropriate package for other organisms

library(data.table)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(preprocessCore)
library(FactoMineR)
library(plotly)
library(stats)
library(limma)
library(impute)
library(clusterProfiler)
library(org.Mm.eg.db)

setwd("C:/Users/alkar/Documents/Wang lab/Proteiomics")

data <- fread('LWANG-combined-data-raw-Maxquant-extracted-newsamples-outliered.csv')

###### Reshape data and visualize
intensity_long <- data %>%
  pivot_longer(cols = starts_with("Intensity_"), names_to = "Sample", values_to = "Intensity")

# Plot
intensity_long$Intensity <- as.numeric(intensity_long$Intensity)

ggplot(intensity_long, aes(x = Sample, y = Intensity, fill = Sample)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Protein Intensity Across Samples", x = "Sample", y = "Intensity")

#Box plot
ggplot(intensity_long, aes(x = Sample, y = Intensity, fill = Sample)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Distribution of Protein Intensities", x = "Sample", y = "Intensity")

# Violin Plot
ggplot(intensity_long, aes(x = Sample, y = Intensity, fill = Sample)) +
  geom_violin(trim = FALSE) +
  theme_minimal() +
  labs(title = "Distribution of Protein Intensities (Violin Plot)", x = "Sample", y = "Intensity")

####### Extract intensity columns and process
class(data)
data <- as.data.frame(data) #convert to data frame to use select function
class(data)

intensity_data <- dplyr::select(data, starts_with("Intensity_"))

# Log2-transform the intensity data
intensity_data <- log2(intensity_data + 1)

# Normalize the data (e.g., quantile normalization)
intensity_data <- normalize.quantiles(as.matrix(intensity_data))
colnames(intensity_data) <- colnames(data %>% select(starts_with("Intensity_")))
rownames(intensity_data) <- data$"Protein IDs"  # Ensure row names match protein IDs
intensity_data[intensity_data == 0] <- NA

head (intensity_data) 

# Check NA data
sum(is.na(intensity_data))

#keep rows with less than 70% missing values
keep_rows_less_na <- function(intensity_data, threshold = 0.7) {
  # Calculate the proportion of missing values in each row
  na_prop <- rowMeans(is.na(intensity_data))
  # Keep rows where the proportion of missing values is less than the threshold
  intensity_data[na_prop < threshold, , drop = FALSE]
}

# Apply the function
intensity_data <- keep_rows_less_na(intensity_data, threshold = 0.7)

# Print the filtered matrix
head(intensity_data)


#Impute Missing Values
intensity_data <- impute.knn(intensity_data)$data


# Check imputed data
sum(is.na(intensity_data))  # Should be 0


# Perform PCA
pca_result <- PCA(intensity_data, scale.unit = TRUE, ncp = 5, graph = FALSE)

# Plot PCA
ggplot(data = as.data.frame(pca_result$ind$coord), aes(x = Dim.1, y = Dim.2)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA of Protein Intensities", x = "PC1", y = "PC2")

# Dot Plot
ggplot(intensity_long, aes(x = Sample, y = Intensity, color = Sample)) +
  geom_point(size = 4, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Protein Intensity Across Samples", x = "Sample", y = "Intensity")

# Create an interactive bar plot
plot_ly(intensity_long, x = ~Sample, y = ~Intensity, type = "bar", color = ~Sample) %>%
  layout(title = "Protein Intensity Across Samples", xaxis = list(title = "Sample"), yaxis = list(title = "Intensity"))

###### Create a metadata table for groups
metadata <- data.frame(
  Sample = c("Intensity_BTRT1", "Intensity_BTRT2", "Intensity_BTRT3", 
             "Intensity_BUS1","Intensity_BUS2", "Intensity_BUS3", "Intensity_CTRL1", 
             "Intensity_CTRL2", "Intensity_CTRL3", "Intensity_CTRL4",
             "Intensity_Sham1", "Intensity_Sham2", "Intensity_Sham3", "Intensity_Sham4"),
  Group = c("BTRT", "BTRT", "BTRT", 
            "BUS","BUS", "BUS", "CTRL", "CTRL", "CTRL", "CTRL",
            "Sham", "Sham", "Sham", "Sham")
)

# View metadata
print(metadata)

# Merge with metadata
intensity_long <- merge(intensity_long, metadata, by = "Sample")

# View the merged data
head(intensity_long)


ggplot(intensity_long, aes(x = Group, y = Intensity, fill = Group)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Protein Intensity by Group", x = "Group", y = "Intensity")

ggplot(intensity_long, aes(x = Group, y = Intensity, fill = Group)) +
  geom_violin(trim = FALSE) +
  theme_minimal() +
  labs(title = "Distribution of Protein Intensities by Group", x = "Group", y = "Intensity")

# Plot heatmap with group annotations
annotation_col <- data.frame(Group = metadata$Group)
rownames(annotation_col) <- metadata$Sample

pheatmap(intensity_data, scale = "row", 
         annotation_col = annotation_col, 
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         main = "Protein Intensity Heatmap by Group",
         fontsize_row = 3, 
         fontsize_col = 3,
         cluster_rows = FALSE,
         cluster_cols = FALSE)


dist_matrix <- dist(t(intensity_data), method = "euclidean")
hc <- hclust(dist_matrix, method = "average")
plot(hc, main = "Hierarchical Clustering of Samples", xlab = "", sub = "")


#### Perform PCA Groupwise
pca_result <- prcomp(t(intensity_data), scale. = TRUE)

# View PCA summary
summary(pca_result)

# Extract PC1 and PC2 scores
pca_scores <- as.data.frame(pca_result$x[, 1:2])  # First two principal components
colnames(pca_scores) <- c("PC1", "PC2")

# Add group information from metadata
pca_scores$Group <- metadata$Group
print(pca_scores)

# Create PCA plot
ggplot(pca_scores, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA Plot (PC1 vs. PC2)", x = "PC1", y = "PC2")

# Extract PC4 and PC3 scores
pca_scores2 <- as.data.frame(pca_result$x[, 3:4])  # First two principal components
colnames(pca_scores2) <- c("PC3", "PC4")

# Add group information from metadata
pca_scores2$Group <- metadata$Group
print(pca_scores2)

# Create PCA plot
ggplot(pca_scores2, aes(x = PC4, y = PC3, color = Group)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA Plot (PC1 vs. PC2)", x = "PC1", y = "PC2")

######################differential protein expression test


#check if colnames are same
print(colnames(intensity_data))
print(metadata$Sample)

# Create a design matrix
design <- model.matrix(~ 0 + Group, data = metadata)
colnames(design) <- levels(as.factor(metadata$Group))

# Fit the linear model
fit <- lmFit(intensity_data, design)

# Define contrasts (e.g., Treatment vs. Control)
contrast_matrix <- makeContrasts(Sham - CTRL, levels = design)

# Fit the contrasts
fit2 <- contrasts.fit(fit, contrast_matrix)

# Apply empirical Bayes moderation
fit2 <- eBayes(fit2)

# Extract results
results <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH")

# View the results
head(results)

# Create a volcano plot
ggplot(results, aes(x = logFC, y = -log10(P.Value), color = adj.P.Val < 0.05 & abs(logFC) > 1)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  scale_color_manual(values = c("gray", "red")) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10(p-value)")

# Add a column to indicate upregulated, downregulated, and non-significant proteins
results$Expression <- ifelse(
  results$adj.P.Val < 0.05 & results$logFC > 1, "Upregulated",
  ifelse(results$adj.P.Val < 0.05 & results$logFC < -1, "Downregulated", "Non-significant")
)

# Create a volcano plot
ggplot(results, aes(x = logFC, y = -log10(P.Value), color = Expression)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  scale_color_manual(values = c("Downregulated" = "blue", "Non-significant" = "gray", "Upregulated" = "red")) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10(p-value)")


# Bar plot
protein_counts <- table(results$Expression)
ggplot(as.data.frame(protein_counts), aes(x = Var1, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Number of Upregulated and Downregulated Proteins", x = "Expression", y = "Count")



################upregulated and downregulated proteins

#check if the proteins-ids are in same sequence

protein_ids <- rownames(intensity_data) 

if (length(protein_ids) != nrow(results)) {
  stop("Replacement data and data frame have different row numbers.")
}


# Assign the replacement data to the data frame
results$Protein_ids <- protein_ids

# View the updated data frame
head(results)

print(head(rownames(results)))
print(head(rownames(intensity_data)))

# Reorder intensity_data to match results
intensity_data <- intensity_data[rownames(results), ]

if (!all(rownames(results) %in% rownames(intensity_data))) {
  stop("Row names of results do not match row names of intensity_data.")
}

significant_proteins <- results[results$adj.P.Val < 0.05, ]
print(nrow(significant_proteins))


significant_intensity <- intensity_data[rownames(significant_proteins), ]
head(significant_intensity)

significant_intensity <- intensity_data[rownames(results[results$adj.P.Val < 0.05, ]), ]

pheatmap(significant_intensity, scale = "row", 
         annotation_col = annotation_col, 
         fontsize_row = 5, 
         fontsize_col = 5,
         clustering_distance_rows = "euclidean", 
         clustering_distance_cols = "euclidean", 
         main = "Heatmap of Significant Proteins")


# Filter top 10 upregulated and downregulated proteins
top_upregulated <- significant_proteins[significant_proteins$Expression == "Upregulated", ][1:10, ]
top_downregulated <- significant_proteins[significant_proteins$Expression == "Downregulated", ][1:10, ]

# Check if ProteinID is included
print(top_upregulated$Protein_ids)
print(top_downregulated$Protein_ids)

# Combine the top proteins
top_proteins <- rbind(top_upregulated, top_downregulated)

# Ensure matching row numbers
if (length(top_proteins$Protein_ids) != length(top_proteins$logFC) || 
    length(top_proteins$Protein_ids) != length(top_proteins$Expression)) {
  stop("Columns have different lengths. Check the data.")
}

# Create the data frame
top_proteins <- data.frame(
  ProteinID = as.factor(top_proteins$Protein_ids),
  logFC = as.numeric(top_proteins$logFC),
  Expression = as.factor(top_proteins$Expression)
)

# Remove missing values
top_proteins <- na.omit(top_proteins)

# Create the bar plot
ggplot(top_proteins, aes(x = reorder(ProteinID, logFC), y = logFC, fill = Expression)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Top 10 Upregulated and Downregulated Proteins", x = "Protein", y = "Log2 Fold Change") +
  coord_flip()



###########customise heatmap#############

# Define a custom color palette
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
# Example row and column annotations
row_annotations <- data.frame(Cluster = rep(c("A", "B"), each = 5))
rownames(row_annotations) <- paste0("Row", 1:10)

col_annotations <- data.frame(Group = rep(c("Control", "Treatment"), each = 5))
rownames(col_annotations) <- paste0("Col", 1:10)# Heatmap with annotations
pheatmap(heatmap_data, 
         annotation_row = row_annotations, 
         annotation_col = col_annotations)
# Define annotation colors
annotation_colors <- list(
  Cluster = c(A = "blue", B = "red"),
  Group = c(Control = "green", Treatment = "orange")
)

# Heatmap with custom annotation colors
pheatmap(heatmap_data, 
         annotation_row = row_annotations, 
         annotation_col = col_annotations, 
         annotation_colors = annotation_colors)
# Custom row and column labels
row_labels <- paste0("Protein", 1:10)
col_labels <- paste0("Sample", 1:10)

# Heatmap with custom labels
pheatmap(heatmap_data, 
         labels_row = row_labels, 
         labels_col = col_labels)

# Heatmap with adjusted font size
pheatmap(heatmap_data, 
         fontsize_row = 8, 
         fontsize_col = 8)
# Heatmap with no row clustering
pheatmap(heatmap_data, 
         cluster_rows = FALSE)

# Heatmap with no column clustering
pheatmap(heatmap_data, 
         cluster_cols = FALSE)
# Heatmap with a title
pheatmap(heatmap_data, 
         main = "Protein Expression Heatmap")
# Save heatmap to a PNG file
pheatmap(heatmap_data, 
         filename = "heatmap.png")

# Save heatmap to a PDF file
pheatmap(heatmap_data, 
         filename = "heatmap.pdf")


##########pathway enrichment ###############

# Filter significant proteins
significant_proteins <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, ]

protein_ids <- rownames(significant_proteins)
print(protein_ids)



# Convert UniProt IDs to Entrez IDs
gene_ids <- mapIds(org.Mm.eg.db, keys = protein_ids, column = "ENTREZID", keytype = "UNIPROT")
print(gene_ids)


# Remove NA values
gene_ids <- na.omit(gene_ids)

print(gene_ids)

# Perform GO enrichment analysis
go_results <- enrichGO(gene = gene_ids, 
                       OrgDb = org.Mm.eg.db, 
                       ont = "BP",  # Biological Process
                       pvalueCutoff = 0.1, 
                       qvalueCutoff = 0.1, 
                       readable = TRUE)

# View results
print(go_results)

barplot(go_results, showCategory = 10, title = "GO Biological Process")

# Dot plot
dotplot(go_results, showCategory = 10, title = "GO Biological Process")

# Dot plot with reduced legend font size
dotplot(go_results, showCategory = 10) +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 8),
    axis.title.x = element_text(size = 8),
    legend.text = element_text(size = 8),  # Legend labels
    legend.title = element_text(size = 8) # Legend title
  )

# Gene-concept network
cnetplot(go_results, showCategory = 5) +
  theme(
    legend.title = element_text(size = 8), 
    legend.text = element_text(size = 8)
  )


# Perform KEGG pathway enrichment
kegg_results <- enrichKEGG(gene = gene_ids, 
                           organism = "hsa", 
                           pvalueCutoff = 0.1, 
                           qvalueCutoff = 0.1)

# View results
head(kegg_results)


##############single protein visulatization

# Example: Extract intensity values for one protein (e.g., "ProteinA")
protein_of_interest <- "Q9Z2X2"
protein_data <- intensity_data[protein_of_interest, ]

# Convert to a data frame for easier plotting
protein_data <- data.frame(Sample = names(protein_data), Intensity = as.numeric(protein_data))
print(protein_data)

# Merge with metadata to add group information
protein_data <- merge(protein_data, metadata, by = "Sample")
print(protein_data)

# Calculate mean intensity for each group
mean_intensity <- protein_data %>%
  group_by(Group) %>%
  summarize(MeanIntensity = mean(Intensity, na.rm = TRUE))

# View the mean intensity
print(mean_intensity)

# Create a bar plot
ggplot(mean_intensity, aes(x = Group, y = MeanIntensity, fill = Group)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = paste("Mean Intensity of", protein_of_interest, "Across Groups"), x = "Group", y = "Mean Intensity")


# Create a line plot
ggplot(mean_intensity, aes(x = Group, y = MeanIntensity, group = 1)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "red", size = 3) +
  theme_minimal() +
  labs(title = paste("Mean Intensity of", protein_of_interest, "Across Groups"), x = "Group", y = "Mean Intensity")

# Interactive bar plot with plotly
plot_ly(mean_intensity, x = ~Group, y = ~MeanIntensity, type = "bar", color = ~Group) %>%
  layout(title = paste("Mean Intensity of", protein_of_interest, "Across Groups"), xaxis = list(title = "Group"), yaxis = list(title = "Mean Intensity"))

