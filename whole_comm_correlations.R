# OTU correlations
## the purpose of this is to run spearmans correlations on the taxa and soil vars
## output to a csv for sorting
## probably need to run on chem vars one by one....

library(dplyr)
library(tidyverse)
library(purrr)
library(ggplot2)

#if you are using envfit below, add that data and use it to filter the samples. 
chem <- read_delim("~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/soilChem_bothdatasets_matchesOnly.csv")%>%
  select(1,4:17)%>%
  filter(sample %in% otu_data$sample)%>%
  na.omit(.)


# Bacteria/archaea
## Import data and make OTU table
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
    values_fill = list(coverage = 0)) 


otu_rel_cov <- otu_table %>%
  rowwise() %>%
  mutate(total_coverage = sum(c_across(-sample))) %>%
  ungroup() %>%
  mutate(across(-c(sample, total_coverage), ~ .x / total_coverage)) %>%
  dplyr::select(-total_coverage) %>%
  filter(sample %in% chem $sample)%>%
  na.omit()


# Fungi
otu_table_filtered <- read_delim(
  "~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/kraken/fungi_gt3_jul14.csv")


#convert otu to relative abundance (filtering with chem and/or metadata sample names so that everything matches)
otu_rel_cov <- otu_table_filtered %>%
  rowwise() %>%
  mutate(total_coverage = sum(c_across(-sample))) %>%
  ungroup() %>%
  mutate(across(-c(sample, total_coverage), ~ .x / total_coverage)) %>%
  dplyr::select(-total_coverage) %>%
  filter(sample %in% chem $sample)%>%
  na.omit()


### loop through OTU abudances
### I did this one by one
names(chem)

otu_table <- otu_rel_cov  %>%
  left_join(chem %>% select(sample, pH), by = "sample")

otu_list <- colnames(otu_table)[!(colnames(otu_table) %in% c("sample", "pH"))]

results <- map_dfr(otu_list, function(otu) {
  test <- cor.test(otu_table[[otu]], otu_table$pH, method = "spearman")
  tibble(
    OTU = otu,
    rho = test$estimate,
    p_value = test$p.value
  )
})

## adjust p values
results <- results %>%
  mutate(FDR = p.adjust(p_value, method = "fdr"))

write.csv(results, "~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/kraken/ALL_pH_FDR.csv")

##  Filter significant results (FDR < 0.05 ) AND filter out weak correlations (rho >= -0.30 & rho <= 0.30)
significant_results <- results %>%
  filter(FDR < 0.05) %>%
  filter(rho >= -0.30 & rho <= 0.30)

write.csv(significant_results , "~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/kraken/SIG_pH_FDR.csv")


