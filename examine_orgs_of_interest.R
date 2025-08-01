## assessing orgs of interest
## I want abundance box plots of 1) orgs of interest vs low, med, high soil var, and 2) orgs of interest in irrigated vs dryland

library(dplyr)
library(tidyr)
library(FSA)
library(purrr)

# Load metadata
metadata <- read_delim("~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/metadata_server(producer info).csv") %>%
  filter(seq == "y", selected == "y")
metadata$Nutrient_Management

# Load soils data
soil_data <- read_delim("~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/Soils data/soil_levels.csv") %>%
  filter(selected == "y", seq == "y") %>%
  na.omit()

# Load organism data and process
df <- read_tsv("~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/singleM/profiles_renamed.tsv") %>%
  filter(sample != "sample") %>%
  mutate(
    coverage = as.numeric(coverage),
    taxonomy = ifelse(is.na(taxonomy) | taxonomy == "", "unclassified", taxonomy)
  ) %>%
  mutate(across(where(is.character), ~ gsub("IN_RICHES", "IN-RICHES", .)))

otu_table <- df %>%
  pivot_wider(
    names_from = taxonomy,
    values_from = coverage,
    values_fill = list(coverage = 0)
  ) %>%
  filter(sample %in% metadata$sample)

otu_table_filtered <- read_delim("~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/kraken/fungi_gt3_jul14.csv")%>%
 mutate(across(where(is.character), ~ gsub("IN_RICHES", "IN-RICHES", .)))%>%
  filter(sample != "CDA_0524_COAllianceCenterEE_23206")

nrow(otu_table_filtered)

# Calculate relative abundances
otu_rel_cov <-otu_table_filtered %>%
  rowwise() %>%
  mutate(total_coverage = sum(c_across(-sample))) %>%
  ungroup() %>%
  mutate(across(-c(sample, total_coverage), ~ .x / total_coverage)) %>%
  select(-total_coverage) %>%
  na.omit()

nrow(otu_rel_cov)

# 6. Join metadata & soil data
metadata_soil <- metadata %>%
  left_join(soil_data, by = "sample")  
metadata_soil$EC



# 7. Merge rel abundance with metadata
df_wide <- otu_rel_cov %>%
  left_join(metadata_soil, by = "sample")


# 8. Subset to organisms of interest
organism_cols <- c(
  "Root; d__Archaea; p__Thermoproteota; c__Nitrososphaeria; o__Nitrososphaerales; f__Nitrososphaeraceae; g__Nitrososphaera; s__Nitrososphaera sp025929855",
  "Root; d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae; o__Chthoniobacterales; f__UBA10450; g__Udaeobacter",
  "Root; d__Bacteria; p__Actinomycetota; c__Thermoleophilia; o__Gaiellales; f__Gaiellaceae; g__JAIBJP01; s__JAIBJP01 sp020029625",
  "Root; d__Archaea; p__Thermoproteota; c__Nitrososphaeria; o__Nitrososphaerales; f__Nitrososphaeraceae; g__TH5896; s__TH5896 sp014523695")
  


organism_cols <- c("p__Ascomycota;c__Sordariomycetes;o__Hypocreales;f__Nectriaceae;g__Fusarium;s__Fusarium_coffeatum",
 "p__Basidiomycota;c__Microbotryomycetes;o__Sporidiobolales;f__Sporidiobolaceae;g__Rhodotorula;s__Rhodotorula_graminis",
  "p__Ascomycota;c__Eurotiomycetes;o__Chaetothyriales;f__Herpotrichiellaceae;g__Fonsecaea;s__Fonsecaea_erecta",
  "p__Ascomycota;c__Eurotiomycetes;o__Onygenales;f__Ajellomycetaceae;g__Blastomyces;s__Blastomyces_percursus",
  "p__Microsporidia;g__Pseudoloma;s__Pseudoloma_neurophilia;NA;NA;NA",
  "p__Mollusca;c__Bivalvia;o__Arcoida;f__Arcidae;g__Anadara;s__Anadara_broughtonii")

# Keep only these organism columns + sample, totalC_group, watering
df_subset <- df_wide %>%
  select(all_of(c("sample", "EC_level", organism_cols)))%>%
  na.omit(.)


df_subset$Nutrient_Management

nrow(df_subset)

# Pivot to long
df_long <- df_subset %>%
  mutate(EC_level = factor(EC_level, levels=c("low","med","high")))%>%
  #mutate(Mulching = factor(Mulching, levels = c("yes", "no")))%>%
  pivot_longer(
    cols = all_of(organism_cols),
    names_to = "organism",
    values_to = "abundance"
  )


# Step 1: build the 'long' dataframe with two variable types
df_clay <- df_long %>%
  mutate(Variable = "per_Total_N_level", Group = per_Total_N_level)
nrow(df_clay)


df_Irrigated <- df_long %>%
  mutate(Variable = "Irrigated", Group = Irrigated)
nrow(df_watering)


# Combine
df_plot <- bind_rows(df_clay, df_Irrigated)

####run stat
# KW by organism
kruskal_results <- df_long %>%
  group_by(organism) %>%
  summarise(
    p_kw = tryCatch(kruskal.test(abundance ~ EC_level)$p.value, error = function(e) NA),
    .groups = "drop"
  ) %>%
  mutate(p_kw_adj = p.adjust(p_kw, method = "fdr"))

print(kruskal_results)

# Get list of organisms to test
#need to run the dunns test - three levels - it looks like we are adjusting the p value 2x, but we are not - we are generating adjusted p values for two diffent tests
significant_orgs <- kruskal_results %>%
  filter(p_kw < 0.05) %>%
  pull(organism)

## run the dunn test
sub_df <- df_long %>% filter(organism =="p__Microsporidia;g__Pseudoloma;s__Pseudoloma_neurophilia;NA;NA;NA")

dunnTest(abundance ~ EC_level, data = sub_df, method = "bh")

  
##plot
new_labels <- c(
  "Root; d__Archaea; p__Thermoproteota; c__Nitrososphaeria; o__Nitrososphaerales; f__Nitrososphaeraceae; g__Nitrososphaera; s__Nitrososphaera sp025929855" = "Nitrososphaera sp025929855",
  "Root; d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae; o__Chthoniobacterales; f__UBA10450; g__Udaeobacter" = "Udaeobacter",
  "Root; d__Bacteria; p__Actinomycetota; c__Thermoleophilia; o__Gaiellales; f__Gaiellaceae; g__JAIBJP01; s__JAIBJP01 sp020029625" = "JAIBJP01 sp020029625",
  "Root; d__Archaea; p__Thermoproteota; c__Nitrososphaeria; o__Nitrososphaerales; f__Nitrososphaeraceae; g__TH5896; s__TH5896 sp014523695" = "TH5896 sp014523695"
)


new_labels <- c(
  "p__Ascomycota;c__Sordariomycetes;o__Hypocreales;f__Nectriaceae;g__Fusarium;s__Fusarium_coffeatum" = "Fusarium_coffeatum",
  "p__Basidiomycota;c__Microbotryomycetes;o__Sporidiobolales;f__Sporidiobolaceae;g__Rhodotorula;s__Rhodotorula_graminis" = "Rhodotorula_graminis",
  "p__Ascomycota;c__Eurotiomycetes;o__Chaetothyriales;f__Herpotrichiellaceae;g__Fonsecaea;s__Fonsecaea_erecta" = "Fonsecaea_erecta",
  "p__Ascomycota;c__Eurotiomycetes;o__Onygenales;f__Ajellomycetaceae;g__Blastomyces;s__Blastomyces_percursus" = "Blastomyces_percursus",
  "p__Microsporidia;g__Pseudoloma;s__Pseudoloma_neurophilia;NA;NA;NA" = "Pseudoloma_neurophilia",
  "p__Mollusca;c__Bivalvia;o__Arcoida;f__Arcidae;g__Anadara;s__Anadara_broughtonii"= "Anadara_broughtonii"
)


p_totalC <- ggplot(df_long, aes(x = EC_level, y = abundance, fill = EC_level)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.5) +
  facet_wrap(~ organism, scales = "free_y",labeller = labeller(organism = new_labels)) +
  labs(
    x = "EC level",
    y = "Relative abundance",
    title = "Abundance across EC level"
  ) +
  scale_fill_manual(values = c("low" = "#1b9e77", "med" = "#d95f02", "high" = "#7570b3")) +  # optional custom colors
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold"),
    legend.title = element_blank()
  )

p_irrigation <- ggplot(df_long, aes(x = Prescribed_Grazing, y = abundance, fill = Prescribed_Grazing)) +
  geom_boxplot(outlier.shape = NA, alpha=0.7) +
  geom_jitter(width = 0.2, size=0.5, alpha=0.5) +
  facet_wrap(~ organism, scales = "free_y",labeller = labeller(organism = new_labels)) +
  labs(
    x = "Prescribed_Grazing",
    y = "Relative abundance",
    title = "Prescribed_Grazing"
  ) +
  scale_fill_manual(values = c("yes" = "#FFAA33", "no" = "#33CC3F")) +  # optional custom colors
 #scale_fill_manual(values = c("Dryland" = "#FFAA33", "Irrigated" = "#33CCFF")) + 
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold")
  )



