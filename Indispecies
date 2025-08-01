# run_indispecies.R

library(tidyverse)
library(indicspecies)

# -----------------------------
# 1. Load and clean metadata
# -----------------------------
metadata <- read_delim("~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/metadata_server(producer info).csv") %>%
  filter(selected == "y", seq == "y") %>%
  select(-Treatment.longform, -treatment) %>%
  na.omit() %>%
  as.data.frame()   # convert to data.frame so rownames actually work
# -----------------------------
# 1b. Load and clean soil vars
# -----------------------------
metadata <- read_delim("~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/Soils data/soil_levels.csv") %>%
  filter(selected == "y", seq == "y") %>%
  na.omit() %>%
  as.data.frame()
# Set rownames to sample IDs
rownames(metadata) <- metadata$sample

# -----------------------------
# 2. Load and clean OTU data
# -----------------------------
df <- read_tsv("~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/singleM/profiles_renamed.tsv") %>%
  filter(sample != "sample") %>%  # remove repeated headers
  mutate(
    coverage = as.numeric(coverage),
    taxonomy = ifelse(is.na(taxonomy) | taxonomy == "", "unclassified", taxonomy)
  ) %>%
  filter(!is.na(coverage)) %>%
  mutate(across(where(is.character), ~ gsub("IN_RICHES", "IN-RICHES", .)))

# Pivot to wide OTU table
otu_table <- df %>%
  pivot_wider(
    names_from = taxonomy,
    values_from = coverage,
    values_fill = list(coverage = 0)
  )

# Keep only samples present in metadata
otu_table <- otu_table %>% filter(sample %in% metadata$sample)

# Calculate relative coverage
otu_rel_cov <- otu_table %>%
  rowwise() %>%
  mutate(total_coverage = sum(c_across(-sample))) %>%
  ungroup() %>%
  mutate(across(-c(sample, total_coverage), ~ .x / total_coverage)) %>%
  select(-total_coverage) %>%
  as.data.frame()

# Set rownames to sample IDs
rownames(otu_rel_cov) <- otu_rel_cov$sample

# Remove 'sample' column for numeric data
otu_data_numeric <- otu_rel_cov %>% select(-sample)

# -----------------------------
# 3. Align metadata to OTU data
# -----------------------------
# Keep only common samples and same order
common_samples <- intersect(rownames(metadata), rownames(otu_rel_cov))
metadata <- metadata[common_samples, ]
otu_rel_cov <- otu_rel_cov[common_samples, ]
otu_data_numeric <- otu_data_numeric[common_samples, ]

# Final check
stopifnot(identical(rownames(metadata), rownames(otu_rel_cov)))

# -----------------------------
# 4. Run indispecies
# -----------------------------
# Define grouping variable (replace 'No_till' if needed)
group <- factor(metadata$per_Clay_level)
#names(metadata)
# Remove NAs in grouping variable if needed
if(any(is.na(group))) {
  keep_rows <- !is.na(group)
  group <- group[keep_rows]
  otu_data_numeric <- otu_data_numeric[keep_rows, ]
}

# Run indicator species analysis
set.seed(42)  # for reproducibility
indval_result <- multipatt(otu_data_numeric, group, func = "IndVal.g", control = how(nperm=999))

# Show summary
summary(indval_result)

# -----------------------------
# 5. Extract and adjust p-values
# -----------------------------
res_df <- as.data.frame(indval_result$sign)
res_df$FDR <- p.adjust(res_df$p.value, method = "fdr")

# View significant indicator species
sig_res <- res_df %>% filter(FDR < 0.05)
#print(sig_res)

write.csv(sig_res, "~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/kraken/indspec_Clay.csv")

## Fungi
# Fungi
otu_table_filtered <- read_delim(
  "~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/kraken/fungi_gt3_jul14.csv")
otu_table <- otu_table_filtered %>% filter(sample %in% metadata$sample)

#convert otu to relative abundance (filtering with chem and/or metadata sample names so that everything matches)
otu_rel_cov <- otu_table %>%
  rowwise() %>%
  mutate(total_coverage = sum(c_across(-sample))) %>%
  ungroup() %>%
  mutate(across(-c(sample, total_coverage), ~ .x / total_coverage)) %>%
  select(-total_coverage) %>%
  as.data.frame()



