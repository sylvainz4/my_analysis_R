#Setting working directory
setwd("/Users/sylvainzemsi/Downloads")
#Installation of R packages
install.packages("readxl")
library(readxl)
library(tibble)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(broom)
library(tidyr)
library(corrplot)
install.packages("circlize")
library(circlize)
# Importation of data
pet <- read_excel("./For_Correlation.xlsx")

#Plot a heatmap with clustering
# Replace non-numeric values with NA
pet[pet == "non-numeric"] <- NA
# Convert all columns to numeric
pet <- as.data.frame(sapply(pet, as.numeric))
# Specify the names of the variables for x-axis and y-axis
x_var_names <- c("A5TZ77", "A5U7L9", "A5TYA9", "A5U2V0", "P9WJM1")

y_var_names <- c("Total_hard_volume", "Kurtosis", "Skewness", 
                 "TLG", "Lobes_affected", "Total_cavity")
# Subset the data
pet_subset <- pet[, c(x_var_names, y_var_names)]
# Extract subsets
x_data <- pet_subset %>% select(all_of(x_var_names))
y_data <- pet_subset %>% select(all_of(y_var_names))

# Function to replace NA with 0
replace_na_with_zero <- function(x) {
  replace(x, is.na(x), 0)
}

# Apply the function to both x_data and y_data
x_data <- as.data.frame(lapply(x_data, replace_na_with_zero))
y_data <- as.data.frame(lapply(y_data, replace_na_with_zero))

# Compute the correlation matrix to avoid repetition on x-axis and y-axis
correlation_matrix <- cor(y_data, x_data, use = "complete.obs")
view(correlation_matrix)
print(correlation_matrix)
# Perform hierarchical clustering
row_dist <- dist(correlation_matrix)
col_dist <- dist(t(correlation_matrix))
row_clustering <- hclust(row_dist)
col_clustering <- hclust(col_dist)
# Reorder the rows and columns based on clustering
row_order <- row_clustering$order
col_order <- col_clustering$order
# Reorder the correlation matrix
ordered_correlation_matrix <- correlation_matrix[row_order, col_order]
# Convert the ordered correlation matrix to a long format data frame for ggplot
correlation_df <- as.data.frame(ordered_correlation_matrix) %>%
  rownames_to_column("y_var") %>%
  pivot_longer(-y_var, names_to = "x_var", values_to = "correlation")
# Ensure x_var and y_var are factors with levels in the clustered order
correlation_df$x_var <- factor(correlation_df$x_var, levels = colnames(ordered_correlation_matrix))
correlation_df$y_var <- factor(correlation_df$y_var, levels = rownames(ordered_correlation_matrix))
# Plot the heatmap of the correlation matrix with clustering
heatmap_plot <- ggplot(correlation_df, aes(x = x_var, y = y_var, fill = correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0,
                       breaks = c(-1, -0.5, 0, 0.5, 1),
                       labels = c("-1", "-0.5", "0", "0.5", "1"),
                       guide = guide_colorbar(title = "Correlation")) +
  labs(x = "Mycobacterium Antigens", y = "PET-CT Scan", title = "Correlation Antigens vs Pet-CT") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

print(heatmap_plot)

#Adding statistical significance to the heatmap plot
#loading required packages
install.packages("ggcorrplot")
install.packages("corrplot")
install.packages("ppcor")
library(ggcorrplot)
library(corrplot)
library(ppcor)
# Assuming x_data and y_data are data frames or matrices (the current plot in the dissertation)
combined_data <- cbind(x_data, y_data)

megadat_correlation <- cor(combined_data, method="kendall", use = "complete.obs")#remember your_matrix is a matrix with only numeric values
megadat_correlation_p <- cor.mtest(combined_data)

corrplot(megadat_correlation, col=COL2(diverging="PRGn"), tl.col="black", tl.cex=0.8, cl.cex=0.8, 
         p.mat=megadat_correlation_p$p, lowCI.mat=megadat_correlation_p$lowCI, 
         uppCI.mat=megadat_correlation_p$uppCI, 
         pch.col="firebrick3", 
         tl.srt=45, 
         diag=FALSE, 
         addgrid.col="grey70", 
         pch.cex=1.5, 
         plotCI="rect", sig.level=0.05, order="original", type="upper")




# Determine the indices for x_data and y_data
x_indices <- 1:ncol(x_data)
y_indices <- (ncol(x_data) + 1):(ncol(x_data) + ncol(y_data))

# Extract the submatrix of correlations between x_data and y_data
megadat_correlation_subset <- megadat_correlation[x_indices, y_indices]

# Define the custom function to compute p-values and confidence intervals for the submatrix
cor.mtest_subset <- function(mat1, mat2, conf.level = 0.95) {
  n <- ncol(mat1)
  m <- ncol(mat2)
  p.mat <- matrix(NA, n, m)
  lowCI.mat <- matrix(NA, n, m)
  uppCI.mat <- matrix(NA, n, m)
  
  for (i in 1:n) {
    for (j in 1:m) {
      tmp <- cor.test(mat1[, i], mat2[, j], conf.level = conf.level)
      p.mat[i, j] <- tmp$p.value
      lowCI.mat[i, j] <- tmp$conf.int[1]
      uppCI.mat[i, j] <- tmp$conf.int[2]
    }
  }
  
  return(list(p.mat = p.mat, lowCI.mat = lowCI.mat, uppCI.mat = uppCI.mat))
}


# Calculate p-values and confidence intervals for the submatrix
megadat_correlation_p <- cor.mtest_subset(x_data, y_data)

# Visualize the submatrix with corrplot
# Ensure row and column names are correctly set
rownames(megadat_correlation_subset) <- colnames(x_data)
colnames(megadat_correlation_subset) <- colnames(y_data)

# Ensure the p-value and CI matrices are correctly labeled and match dimensions
rownames(megadat_correlation_p$p.mat) <- colnames(x_data)
colnames(megadat_correlation_p$p.mat) <- colnames(y_data)

rownames(megadat_correlation_p$lowCI.mat) <- colnames(x_data)
colnames(megadat_correlation_p$lowCI.mat) <- colnames(y_data)

rownames(megadat_correlation_p$uppCI.mat) <- colnames(x_data)
colnames(megadat_correlation_p$uppCI.mat) <- colnames(y_data)

# Visualize the submatrix with corrplot
corrplot(megadat_correlation_subset, 
         col = COL2(diverging = "PRGn"), 
         tl.col = "black", 
         tl.cex = 1.2, 
         cl.cex = 0.8, 
         p.mat = megadat_correlation_p$p.mat, 
         lowCI.mat = megadat_correlation_p$lowCI.mat, 
         uppCI.mat = megadat_correlation_p$uppCI.mat, 
         pch.col = "firebrick3", 
         tl.srt = 45, 
         diag = FALSE, 
         addgrid.col = "grey70", 
         pch.cex = 3, 
         plotCI = "rect", 
         sig.level = 0.05, 
         order = "original", 
         type = "full") 
 #Spirometry correlation
# Convert all columns to numeric
spiro <- as.data.frame(sapply(pet, as.numeric))
# Specify the names of the variables for x-axis and y-axis
x_var_names <- c("A5TZ77", "A5U7L9", "A5TYA9", "A5U2V0", "P9WJM1")

z_var_names <- c("fev1", "fvc", "fef2575", "fef75")
# Subset the data
spiro_subset <- spiro[, c(x_var_names, z_var_names)]

# Replace NA values with 0 in all columns
spiro_subset[is.na(spiro_subset)] <- 0

# Extract subsets using dplyr::select
x_data_s <- spiro_subset %>% dplyr::select(all_of(x_var_names))
z_data_s <- spiro_subset %>% dplyr::select(all_of(z_var_names))


# Assuming x_data and y_data are data frames or matrices (the current plot in the dissertation)
combined_data_s <- cbind(x_data_s, z_data_s)

megadat_correlation <- cor(combined_data_s, method="kendall", use = "complete.obs")#remember your_matrix is a matrix with only numeric values
megadat_correlation_p <- cor.mtest(combined_data_s)

corrplot(megadat_correlation, col=COL2(diverging="PRGn"), tl.col="black", tl.cex=0.8, cl.cex=0.8, 
         p.mat=megadat_correlation_p$p, lowCI.mat=megadat_correlation_p$lowCI, 
         uppCI.mat=megadat_correlation_p$uppCI, 
         pch.col="firebrick3", 
         tl.srt=45, 
         diag=FALSE, 
         addgrid.col="grey70", 
         pch.cex=1.5, 
         plotCI="rect", sig.level=0.05, order="original", type="upper")

# Subset the correlation matrix to include only the relevant correlations
relevant_correlation <- megadat_correlation[1:length(x_var_names), (length(x_var_names)+1):ncol(megadat_correlation)]
relevant_p_values <- megadat_correlation_p$p[1:length(x_var_names), (length(x_var_names)+1):ncol(megadat_correlation)]
relevant_lowCI <- megadat_correlation_p$lowCI[1:length(x_var_names), (length(x_var_names)+1):ncol(megadat_correlation)]
relevant_uppCI <- megadat_correlation_p$uppCI[1:length(x_var_names), (length(x_var_names)+1):ncol(megadat_correlation)]

# Custom labels for x and y axes
x_labels <- colnames(x_data_s)
y_labels <- colnames(z_data_s)

# Plot the correlation matrix with corrplot
corrplot(relevant_correlation, col = COL2(diverging = "PRGn"), 
         tl.col = "black", tl.cex = 0.8, cl.cex = 0.8,
         p.mat = relevant_p_values, lowCI.mat = relevant_lowCI, 
         uppCI.mat = relevant_uppCI, pch.col = "firebrick3", 
         tl.srt = 45, addgrid.col = "grey70", 
         pch.cex = 1.5, plotCI = "rect", sig.level = 0.05, 
         order = "original", method = "circle",
         tl.labels = y_labels, mar = c(1, 1, 1, 1), cl.ratio = 0.5)