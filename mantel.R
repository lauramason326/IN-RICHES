## KB code, LM edits
library(dplyr)
library(tidyr)
library(tidyverse)
library(vegan)
library(readr)
library(ggplot2)

metadata<-read_delim("~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/soilChem_bothdatasets_matchesOnly.csv")%>%
  filter(selected=="y")%>%
  filter(use=="y")

# Data (or features) file
# This file is tab-delimited,Formatted with features as columns and samples as rows.
#Import data and make OTU table
df <- read_tsv("~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/singleM/profiles_renamed.tsv")%>%
  filter(sample != "sample") %>%  # Remove repeated headers
  mutate(
    coverage = as.numeric(coverage),
    taxonomy = ifelse(is.na(taxonomy) | taxonomy == "", "unclassified", taxonomy)
  ) %>%
  filter(!is.na(coverage)) %>%
  mutate(across(where(is.character), ~ gsub("IN_RICHES", "IN-RICHES", .)))

otu_table <- df %>%
  pivot_wider(
    names_from = taxonomy,
    values_from = coverage,
    values_fill = list(coverage = 0)) %>%
  # filter(str_detect(sample, "IN-RICHES"))
  filter(sample %in% metadata$sample)



otu_rel_cov <- otu_table %>%
  rowwise() %>%
  mutate(total_coverage = sum(c_across(-sample))) %>%
  ungroup() %>%
  mutate(across(-c(sample, total_coverage), ~ .x / total_coverage)) %>%
  dplyr::select(-total_coverage) %>%
  na.omit()


ncol(metadata)
nrow(metadata)

# Get matching metadata
meta <- metadata %>%
  select(1,4:17)%>%
  na.omit(.)

otu_data <- otu_rel_cov  %>% 
  filter(sample %in% meta $sample)%>%
  column_to_rownames("sample")

meta <- meta %>%
filter(sample %in% otu_rel_cov $sample)%>%
column_to_rownames("sample") %>%
  mutate(across(everything(), as.numeric)) 



nrow(meta)

colnames(meta) <- make.names(colnames(meta), unique = TRUE)



# Create distance matrix from OTU table (drop sample column)
comm_dist <- otu_data %>%
  dist(method = "euclidean")


# Initialize dataframe for Mantel test results
  mantel_table <- data.frame(matrix(nrow = ncol(meta) - 1, ncol = 3))  # -1 to exclude Sample_ID column
  colnames(mantel_table) <- c("predictor", "mantelr", "pval")
  
  # Perform Mantel test for each column in meta
 for (i in 2:ncol(meta)) {  # Start from 2 to skip Sample_ID column
    # Calculate distance matrix for the current metadata variable
    metadata_dist <- dist(meta[, i, drop = FALSE], method = "euclidean")
    
    # Perform Mantel test
    mantel_result <- mantel(comm_dist, metadata_dist)
    
    # Store results
    mantel_table[i - 1, "predictor"] <- colnames(meta)[i]
    mantel_table[i - 1, "mantelr"] <- mantel_result$statistic
    mantel_table[i - 1, "pval"] <- mantel_result$signif
 }
  
  mantel_table <- data.frame(matrix(nrow = ncol(meta), ncol = 3))
  colnames(mantel_table) <- c("predictor", "mantelr", "pval")
  
  # Loop over metadata columns
  for (i in 1:ncol(meta)) {
    metadata_dist <- dist(meta[, i, drop = FALSE], method = "euclidean")
    mantel_result <- mantel(comm_dist, metadata_dist)
    mantel_table[i, "predictor"] <- colnames(meta)[i]
    mantel_table[i, "mantelr"] <- mantel_result$statistic
    mantel_table[i, "pval"] <- mantel_result$signif
  }
  
print(mantel_table)

write.csv(mantel_table, "~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/singleM/mantel_results_bothDatasets.csv")







#barplot
#I added the category by hand
mantel_table<-read_delim("~/Documents/Wrighton_lab/Agribiome/Schipanski_CoverCrop_2023/ITS/Mantel_results_summer.csv")

mantel_table$Group <- ifelse(mantel_table$`pval` <= 0.1, "Significant", "Not Significant")

# Ensure 'Significant' is ordered before 'Not Significant'
mantel_table$Group <- factor(mantel_table$Group, levels = c("Significant", "Not Significant"))

# Order the data first by Group, then by Category, and then by mantelr within each group
mantel_table <- mantel_table %>%
  arrange(-mantelr)

#mantel_table <- mantel_table %>%
#arrange(category,-mantelr)

# Create an ordered factor for plotting
mantel_table$name_Ordered <- factor(mantel_table$predictor, levels = mantel_table$predictor)


# Plot with legend and adjusted y-axis ticks
mantel_plot <- ggplot(mantel_table, aes(x = name_Ordered, y = mantelr, fill = Group)) + 
  geom_bar(stat = "identity") + 
  geom_text(aes(label = paste("P =", sprintf("%.3f", pval))), vjust = -0.5, size = 2.5,angle = 45) +
  scale_fill_manual(values = c("Significant" = "#b35806", "Not Significant" = "#D6D6D6"),
                    labels = c("Significant", "Not Significant")) +
  theme_minimal() +
  labs(title = "",
       x = "",
       y = "Rho") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.line.x = element_line(color = "black", size = 1),
        axis.line.y = element_line(color = "black", size = 1),
        legend.position = "none") + # Hide the legend if not needed
  scale_y_continuous(breaks = c(-1,0,1,2,3,4))
