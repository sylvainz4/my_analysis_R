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
install.packages("ggcorrplot")
install.packages("corrplot")
install.packages("ppcor")
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
x_var_names <- c("A5TYA9")

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

Antigens <- pet[, c("PID", "A5TZ77", "A5U7L9", "A5TYA9", "A5U2V0", "P9WJM1")]

# Impute NA values with 0
Antigens_imputed <- Antigens %>%
  mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))

# Shapiro-Wilk test for normality
shapiro_test_results <- lapply(Antigens_imputed, shapiro.test)
shapiro_test_results

# Reshape data to long format
Antigens_long <- Antigens_imputed %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value")

# Perform Kruskal-Wallis test
kruskal_result <- kruskal.test(value ~ variable, data = Antigens_long)

# Print Kruskal-Wallis test result
print(kruskal_result)

# Install and load necessary package
install.packages("dunn.test")
library(dunn.test)

# Perform Dunn's test for multiple comparisons
dunn_result <- dunn.test(Antigens_long$value, Antigens_long$variable, method = "bonferroni")

# Print Dunn's test result
print(dunn_result)

# Create boxplot
ggplot(Antigens_long, aes(x = variable, y = value)) +
  geom_point() +
  labs(title = "Boxplot of Antigens", x = "Variable", y = "Value")
# Identify participants with more than one non-zero value
Antigens_long <- Antigens_imputed %>%
  pivot_longer(cols = -PID, names_to = "variable", values_to = "value") %>%
  group_by(PID) %>%
  summarise(non_zero_count = sum(value != 0)) %>%
  filter(non_zero_count > 1)

# Merge with original data to get values
Antigens_filtered <- Antigens_imputed %>%
  filter(PID %in% Antigens_long$PID) %>%
  pivot_longer(cols = -PID, names_to = "variable", values_to = "value")


# Convert 'value' column to numeric (handles "NA" strings as true NA)
Antigens_filtered <- Antigens_filtered %>%
  mutate(value = as.numeric(value))

# Replace NA values with 0
Antigens_filtered <- Antigens_filtered %>%
  mutate(value = ifelse(is.na(value), 0, value))


# Create visualization
ggplot(Antigens_filtered, aes(x = variable, y = factor(PID))) +
  geom_tile(aes(fill = value != 0), color = "white") +
  scale_fill_manual(values = c("FALSE" = "lightgrey", "TRUE" = "blue")) +
  labs(title = "Participants with Non-Zero and Zero Values",
       x = "Variable",
       y = "Participant",
       fill = "Non-Zero Value") +
  theme_minimal()

# Reorder the PID factor levels based on the desired order
Antigens_filtered$PID <- factor(Antigens_filtered$PID, levels = unique(Antigens_filtered$PID))

# Create visualization with gradient fill based on actual value
ggplot(Antigens_filtered, aes(x = variable, y = factor(PID))) +
  geom_tile(aes(fill = value), color = "white") +
  scale_fill_gradient(low = "lightgrey", high = "blue") +
  labs(title = "Plot of Participants and Antigen Identified",
       x = "Mycobacterium Tuberculosis Antigens",
       y = "Participant",
       fill = "Antigen Value") +
  theme_minimal()
