## Maaslin2
## Last updated: 6/11/25

###################
# load libraries
##############

library(Maaslin2)
library(tidyverse) 
library(RColorBrewer)
library(dplyr)

##################
#Metadata file
#Formatted with features as columns and samples as rows.
metadata<-read_delim("~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/metadata_server(producer info).csv")%>%
  filter(selected=="y")%>%
  column_to_rownames(., var = "sample")


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
  filter(sample %in% rownames(metadata))

otu_table<- read_delim("~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/kraken/fungi_gt3_jul14.csv")

otu_rel_cov <- otu_table %>%
  dplyr::rowwise() %>%
  dplyr::mutate(total_coverage = sum(dplyr::c_across(dplyr::where(is.numeric)), na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ .x / total_coverage)) %>%
  dplyr::select(-total_coverage) %>%
  tidyr::drop_na()


# Step 2: Filter out lineages with < 0.05 total abundance
column_sums <- colSums(otu_rel_cov[, -1])  # assumes first column is 'sample'
keep_columns_abundance <- names(column_sums[column_sums >= 0.05])

otu_rel_cov <- otu_rel_cov %>%
  select(sample, all_of(keep_columns_abundance))

##If only completing the first filtering
otu_rel_cov <- otu_rel_cov %>%
  column_to_rownames(var = "sample")

# Step 3: Filter out lineages present in < 10% of samples
otu_data_only <- otu_rel_cov[, -1]  # exclude sample column
prevalence <- colMeans(otu_data_only > 0)
keep_columns_prevalence <- names(prevalence[prevalence >= 0.10])

# Step 4: Apply both filters and finalize
otu_rel_cov <- otu_rel_cov %>%
  select(sample, all_of(keep_columns_abundance), all_of(keep_columns_prevalence)) %>%
  select(sample, all_of(intersect(keep_columns_abundance, keep_columns_prevalence))) %>%
  column_to_rownames(var = "sample")


names(metadata)
# run Maaslin2
Maaslin2(input_data = otu_rel_cov,
         input_metadata = metadata,
         output = "~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/kraken/Minimize_Soil_Disturbance_both_filters_random", 
         fixed_effects = ("Minimize_Soil_Disturbance"),
         random_effects = c("locationName"))

## cat all results, remove all but one header row and save
########volcano plot
df <- read_delim("~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/kraken/all_results_random_bothFilters.tsv")  # replace with your actual file path
str(df)

df <- df %>%
  mutate(
    coef = as.numeric(coef),
    stderr = as.numeric(stderr),
    pval = as.numeric(pval),
    qval = as.numeric(qval)
  )



# Process to add log-transformed q-values and enrichment labels
plot_data <- df %>%
  mutate(
    log_q = -log10(qval),
    Enrichment = case_when(
      qval < 0.1 & coef > 0 ~ paste(metadata, "Enriched (↑)"),
      qval < 0.1 & coef < 0 ~ paste(metadata, "Enriched (↓)"),
      TRUE ~ "Non-significant"
    ),
    feature_label = str_replace_all(feature, "\\.\\.", " > ")
  ) %>%
  filter(!grepl("^dataset Enriched", Enrichment))

#check
plot_data %>% filter(is.na(Enrichment))


plot_data %>%
  dplyr::count(Enrichment) %>%
  arrange(desc(n))

ggplot(plot_data, aes(x = coef, y = log_q)) +
  geom_point(aes(fill = Enrichment),
             shape = 21, size = 3, stroke = 0, color = "black", alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c(
    "Conservation_Cover Enriched (↑)" = "#82C503",
    "Conservation_Cover Enriched (↓)" = "#82C503",
    "Continuous_Living_Root Enriched (↑)" = "#CFFC00",
    "Continuous_Living_Root Enriched (↓)" = "#CFFC00",
    "Cover_Crop Enriched (↑)" = "#E56A54",
    "Cover_Crop Enriched (↓)" = "#E56A54",
    "Irrigated Enriched (↑)" = "#D9782D",
    "Irrigated Enriched (↓)" = "#D9782D",
    "Livestock_Integration Enriched (↑)" = "#FF1493",  # updated yellow
    "Livestock_Integration Enriched (↓)" = "#FF1493",
    "Minimize_Soil_Disturbance Enriched (↑)" = "#7E5475",
    "Minimize_Soil_Disturbance Enriched (↓)" = "#7E5475",
    "Mulching Enriched (↑)" = "#999900",
    "Mulching Enriched (↓)" = "#999900",
    "No_till Enriched (↑)" = "#008FB3",
    "No_till Enriched (↓)" = "#008FB3",
    "Nutrient_Management Enriched (↑)" ="#4afdb4",
    "Nutrient_Management Enriched (↓)" = "#4afdb4",
    "Plant_Diversity Enriched (↑)" = "#105456",
    "Plant_Diversity Enriched (↓)" = "#105456",
    "Prescribed_Grazing Enriched (↑)" = "#B86BFF",  
    "Prescribed_Grazing Enriched (↓)" = "#B86BFF",
    "Reduced_Till Enriched (↑)" = "#FFC038",
    "Reduced_Till Enriched (↓)" = "#FFC038",
    "Soil_Armor Enriched (↑)" = "#C8C372",
    "Soil_Armor Enriched (↓)" = "#C8C372",
    "Soil_Carbon_Amendment Enriched (↑)" = "#FFE135", 
    "Soil_Carbon_Amendment Enriched (↓)" = "#FFE135", 
    "TreeShrub Enriched (↑)" = "#fbb7bb",
    "TreeShrub Enriched (↓)" = "#fbb7bb",
    "Non-significant" = "gray80"
  )) +
  labs(
    x = "Coefficient",
    y = "-log10(q-value)",
    fill = "Enrichment"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

## barplot
data <- read_delim("~/Documents/Wrighton_lab/Agribiome/AG ML project/enriched.csv") 
str(data)
df <- data %>%
  mutate(count = ifelse(enrichment == "no", -count, count),
         label = paste0(variable, " Enriched (", ifelse(enrichment == "yes", "↑", "↓"), ")"))

# Define your custom colors
custom_colors <- c(
  "Irrigated Enriched (↑)" = "#D9782D",
  "Irrigated Enriched (↓)" = "#D9782D",
  "Livestock_Integration Enriched (↑)" = "#FF1493",
  "Livestock_Integration Enriched (↓)" = "#FF1493",
  "Mulching Enriched (↑)" = "#999900",
  "Mulching Enriched (↓)" = "#999900",
  "No_till Enriched (↑)" = "#008FB3",
  "No_till Enriched (↓)" = "#008FB3",
  "Nutrient_Management Enriched (↑)" ="#4afdb4",
  "Nutrient_Management Enriched (↓)" = "#4afdb4",
  "Prescribed_Grazing Enriched (↑)" = "#B86BFF",  
  "Prescribed_Grazing Enriched (↓)" = "#B86BFF",
  "Reduced_Till Enriched (↑)" = "#FFC038",
  "Reduced_Till Enriched (↓)" = "#FFC038",
  "Soil_Armor Enriched (↑)" = "#C8C372",
  "Soil_Armor Enriched (↓)" = "#C8C372",
  "TreeShrub Enriched (↑)" = "#fbb7bb",
  "TreeShrub Enriched (↓)" = "#fbb7bb"
)

# Plot
ggplot(df, aes(x = reorder(variable, abs(count), sum), y = count, fill = label)) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = custom_colors) +
  labs(x = NULL, y = "Count", fill = "Variable & Enrichment") +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 10)
  ) +
  geom_text(aes(label = abs(count)),
            position = position_stack(vjust = 0.5), size = 3, color = "black")
