##fungi - beta Div
library(vegan)
library(tidyverse)
library(grid)
library(forcats)
library(dplyr)
library(purrr)
library(tibble)


metadata<-read_delim("~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/metadata_server(producer info).csv")%>%
  #filter(grepl("IN-RICHES", sample))%>%
  #filter(seq == "y") %>%
  filter(selected == "y") %>%
  filter(sample != "CDA_0524_COAllianceCenterEE_23206")

#Kraken beta div
# upstream processing - this is needed for basically everything 
## upload and data wrangling
### upload "longer" format from concatenating the krona otu files from Bracken
data_longer <- read_delim(
  "~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/kraken/Fungi_otu_jul14.txt", 
  col_names = FALSE, delim = "\t"
)
colnames(data_longer) <- c("sample", "Count", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Add joined taxonomy string for convenience
data_longer <- data_longer %>%
  mutate(Taxonomy = paste(Phylum, Class, Order, Family, Genus, Species, sep = ";"))

# Pivot wider: samples as rows, use this for the "otu table"
data_wider <- data_longer %>%
  pivot_wider(
    names_from = Taxonomy,
    values_from = Count,
    values_fill = list(Count = 0))%>%
  select(1,8:1411)%>% #get rid of extra taxonomic info
  group_by(sample) %>% #collapse by sample
  summarize(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)), .groups = "drop")%>%
  mutate(across(where(is.character), ~ gsub("IN_RICHES", "IN-RICHES", .)))%>%
  filter(sample != "CDA_0524_COAllianceCenterEE_23206")

#filter out every taxon that appears less than 3 times 
col_totals <- colSums(select(data_wider, -sample), na.rm = TRUE)
# get names of taxa whose total abundance >= 3
taxa_to_keep <- names(col_totals)[col_totals >= 3]
#  keep only 'sample' + these taxon columns
otu_table_filtered <- data_wider %>%
  select(sample, all_of(taxa_to_keep))




otu_table_filtered <- read_delim("~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/kraken/fungi_gt3_jul14.csv")%>%
#  filter(grepl("IN_RICHES", sample))%>%
  mutate(across(where(is.character), ~ gsub("IN_RICHES", "IN-RICHES", .)))%>%
  filter(sample != "CDA_0524_COAllianceCenterEE_23206")

#convert otu to relative abundance (filtering with chem and/or metadata sample names so that everything matches)
otu_rel_cov <- otu_table_filtered %>%
  rowwise() %>%
  mutate(total_coverage = sum(c_across(-sample))) %>%
  ungroup() %>%
  mutate(across(-c(sample, total_coverage), ~ .x / total_coverage)) %>%
  dplyr::select(-total_coverage) %>%
  na.omit()

meta <- metadata%>% 
  inner_join(otu_rel_cov, by="sample")
nrow(meta)


nrow(otu_rel_cov)



ncol(meta)
nrow(meta)

#write.csv(meta,"~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/singleM/meta_first_pass.csv" )
table <-meta%>% 
  select(34:706) 


data.t <- as.matrix(vegdist(table, method = 'bray'))
names(metadata)

set.seed(1234)
#PERMANOVA - for the first pass, I am keeping this really simple
vars <- c("Organic", "Irrigated", "Conservation_Cover", "Conservation_Crop_Rotation", 
          "No_till", "Soil_Carbon_Amendment", "Cover_Crop", "Reduced_Till", 
          "Mulching", "Forage_Harvest_Management", "Prescribed_Grazing", 
          "Nutrient_Management", "TreeShrub")

# Count number of unique levels in each
sapply(meta[, vars], function(x) length(unique(x)))

priciples <- adonis2(data.t ~ Organic + Irrigated+Conservation_Cover+Conservation_Crop_Rotation+
                        No_till+Soil_Carbon_Amendment+Cover_Crop+Reduced_Till+
                        Mulching+ Forage_Harvest_Management+ Prescribed_Grazing+
                        Nutrient_Management+ TreeShrub,  data = meta,
                        permutations = 999,
                        by = "terms",
                        strata = meta$locationName)




adonis_formulas <- list(
  #dataset_results=data.t~dataset,
#  cropping_system_results = data.t ~ Cropping_System,
  locationName_results = data.t ~ locationName,
  Irrigation_results = data.t ~ Irrigated,
  #organic_nutrient_mangement_results = data.t ~ Organic,
  Soil_Armor_results = data.t ~ `Soil_Armor`,
  Livestock_Integration_results = data.t ~ `Livestock_Integration`,
  Plant_Diversity_results = data.t ~ `Plant_Diversity`,
  Minimize_Soil_Disturbance_results = data.t ~ `Minimize_Soil_Disturbance`,
  Conservation_Cover_results = data.t ~ `Conservation_Cover`,
 # Conservation_Crop_Rotation_results = data.t ~ `Conservation_Crop_Rotation`,
  no_till_results = data.t ~ `No_till`,
  Soil_Carbon_Amendment_results = data.t ~ Soil_Carbon_Amendment,
  Cover_Crop_results = data.t ~ Cover_Crop,
  Reduced_Till_results = data.t ~ `Reduced_Till`,
 # Mulching_results = data.t ~ Mulching,
  #Forage_Harvest_Management_results = data.t ~ Forage_Harvest_Management,
  Prescribed_Grazing_results = data.t ~ `Prescribed_Grazing`
  #Row_Arrangement_results = data.t ~ Row_Arrangement,
  #Nutrient_Management_results = data.t ~ Nutrient_Management,
  #TreeShrub_results = data.t ~ TreeShrub
)

results_with_aic <- map_dfr(names(adonis_formulas), function(var_name) {
  # Run adonis2 with correct strata (actual vector, not string)
  model <- adonis2(
    formula = adonis_formulas[[var_name]], 
    data = meta, 
    permutations = 999, 
    by = "terms",
      strata = meta$locationName
  )
  
  model_df <- as.data.frame(model) %>%
    rownames_to_column("term") %>%
    mutate(variable = var_name)
  
  SSR_resid <- model_df %>% filter(term == "Residual") %>% pull(SumOfSqs)
  df_model <- model_df %>% filter(term != "Residual" & term != "Total") %>% pull(Df) %>% sum()
  n <- nrow(meta)
  
  # Calculate approximate AIC for PERMANOVA
  AIC_like <- n + log(SSR_resid / n) + 2 + df_model
  
  model_df %>%
    filter(term != "Total") %>%
    mutate(SSR = SumOfSqs, AIC_approx = AIC_like) %>%
    select(variable, term, Df, SSR, R2, F, `Pr(>F)`, AIC_approx)
})
#print results to csv
write.csv(results_with_aic,"~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/kraken/both_datasetsaic_approx_1x_wStrata_locationname_actualFungi.csv")


## permanova  - all metrics
all_var_permanova = adonis2(data.t ~ locationName+
                              Irrigated+
                              Reduced_Till+
                              Cropping_System+
                              Nutrient_Management+
                              Conservation_Crop_Rotation+
                              Soil_Armor, data = meta, permutations = 999, by = "terms")



all_permanova2 = adonis2(data.t ~ Cropping_System+locationName+Soil_Armor, data = meta, permutations = 999, by = "terms")

all_permanova3 = adonis2(data.t ~ Cropping_System+locationName, data = meta, permutations = 999, by = "terms")
practices_all <- adonis2(data.t~Irrigated + Conservation_Cover + Conservation_Crop_Rotation + No_till + Soil_Carbon_Amendment+
                           Cover_Crop + Reduced_Till+ Prescribed_Grazing + Cropping_System, data = meta, by = "terms",strata = meta$locationName)


princ_all <- adonis2(data.t ~ Livestock_Integration + Continuous_Living_Root + Plant_Diversity + Minimize_Soil_Disturbance +Soil_Armor, data = meta, by = "terms",strata = meta$locationName)

princ_all <- adonis2(data.t ~ Livestock_Integration + Continuous_Living_Root + Plant_Diversity + Minimize_Soil_Disturbance +Soil_Armor, data = meta, by = "terms",strata = meta$Cropping_System)

## just regen practices and principles
practices_all <- adonis2(data.t~Organic + Irrigated+Conservation_Cover + Conservation_Crop_Rotation + No_till + Soil_Carbon_Amendment+
                           Cover_Crop + Reduced_Till + Mulching + Prescribed_Grazing + TreeShrub+
                           Nutrient_Management, data = meta, by = "terms", strata = meta$locationName)

practices2 <- adonis2(data.t~ Irrigated + No_till+ Reduced_Till + Nutrient_Management, data = meta, by = "terms")

adonis_formulas <- list(
notill_reduced = data.t ~ No_till + Reduced_Till,
notill_res = data.t~No_till,
reduced_res = data.t~Reduced_Till)

# Initialize list to collect results
results_list <- list()

# Loop through each formula by name
for (var_name in names(adonis_formulas)) {
  model <- adonis2(
    formula = adonis_formulas[[var_name]],
    data = meta,
    permutations = 999,
    by = "terms",
    strata = Cropping_System  # add strata if needed
  )
  
  model_df <- as.data.frame(model) %>%
    rownames_to_column("term") %>%
    mutate(variable = var_name)
  
  SSR_resid <- model_df %>% filter(term == "Residual") %>% pull(SumOfSqs)
  df_model <- model_df %>% filter(!term %in% c("Residual", "Total")) %>% pull(Df) %>% sum()
  n <- nrow(meta)
  
  # Calculate approximate AIC for PERMANOVA
  AIC_like <- n + log(SSR_resid / n) + 2 + df_model
  
  model_df <- model_df %>%
    filter(term != "Total") %>%
    mutate(SSR = SumOfSqs, AIC_approx = AIC_like) %>%
    select(variable, term, Df, SSR, R2, F, `Pr(>F)`, AIC_approx)
  
  results_list[[var_name]] <- model_df
}

# Combine all results into one dataframe
results_with_aic <- bind_rows(results_list)

# Write to CSV
write.csv(results_with_aic, "~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/kraken/aic_practices_GGL.csv")


##principles
adonis_formulas<- list(
  soilarmor_livestock = data.t ~ Soil_Armor + Livestock_Integration,
  soilarmor_livingroot = data.t ~ Soil_Armor + Continuous_Living_Root,
  soilarmor_plantdiv = data.t ~ Soil_Armor + Plant_Diversity,
  soilarmor_minsoil = data.t ~ Soil_Armor + Minimize_Soil_Disturbance,
  livestock_livingroot = data.t ~ Livestock_Integration + Continuous_Living_Root,
  livestock_plantdiv = data.t ~ Livestock_Integration + Plant_Diversity,
  livestock_minsoil = data.t ~ Livestock_Integration + Minimize_Soil_Disturbance,
  livingroot_plantdiv = data.t ~ Continuous_Living_Root + Plant_Diversity,
  livingroot_minsoil = data.t ~ Continuous_Living_Root + Minimize_Soil_Disturbance,
  plantdiv_minsoil = data.t ~ Plant_Diversity + Minimize_Soil_Disturbance,
  soilarmor_livestock_livingroot = data.t ~ Soil_Armor + Livestock_Integration + Continuous_Living_Root,
  soilarmor_livestock_plantdiv = data.t ~ Soil_Armor + Livestock_Integration + Plant_Diversity,
  soilarmor_livestock_minsoil = data.t ~ Soil_Armor + Livestock_Integration + Minimize_Soil_Disturbance,
  soilarmor_livingroot_plantdiv = data.t ~ Soil_Armor + Continuous_Living_Root + Plant_Diversity,
  soilarmor_livingroot_minsoil = data.t ~ Soil_Armor + Continuous_Living_Root + Minimize_Soil_Disturbance,
  soilarmor_plantdiv_minsoil = data.t ~ Soil_Armor + Plant_Diversity + Minimize_Soil_Disturbance,
  livestock_livingroot_plantdiv = data.t ~ Livestock_Integration + Continuous_Living_Root + Plant_Diversity,
  livestock_livingroot_minsoil = data.t ~ Livestock_Integration + Continuous_Living_Root + Minimize_Soil_Disturbance,
  livestock_plantdiv_minsoil = data.t ~ Livestock_Integration + Plant_Diversity + Minimize_Soil_Disturbance,
  livingroot_plantdiv_minsoil = data.t ~ Continuous_Living_Root + Plant_Diversity + Minimize_Soil_Disturbance,
  soilarmor_livestock_livingroot_plantdiv = data.t ~ Soil_Armor + Livestock_Integration + Continuous_Living_Root + Plant_Diversity,
  soilarmor_livestock_livingroot_minsoil = data.t ~ Soil_Armor + Livestock_Integration + Continuous_Living_Root + Minimize_Soil_Disturbance,
  soilarmor_livestock_plantdiv_minsoil = data.t ~ Soil_Armor + Livestock_Integration + Plant_Diversity + Minimize_Soil_Disturbance,
  soilarmor_livingroot_plantdiv_minsoil = data.t ~ Soil_Armor + Continuous_Living_Root + Plant_Diversity + Minimize_Soil_Disturbance,
  livestock_livingroot_plantdiv_minsoil = data.t ~ Livestock_Integration + Continuous_Living_Root + Plant_Diversity + Minimize_Soil_Disturbance
)



# Initialize list to collect results
results_list <- list()

# Loop through each formula by name
for (var_name in names(adonis_formulas)) {
  model <- adonis2(
    formula = adonis_formulas[[var_name]],
    data = meta,
    permutations = 999,
    by = "terms"
    # strata = your_strata_vector_here  # add strata if needed
  )
  
  model_df <- as.data.frame(model) %>%
    rownames_to_column("term") %>%
    mutate(variable = var_name)
  
  SSR_resid <- model_df %>% filter(term == "Residual") %>% pull(SumOfSqs)
  df_model <- model_df %>% filter(!term %in% c("Residual", "Total")) %>% pull(Df) %>% sum()
  n <- nrow(meta)
  
  # Calculate approximate AIC for PERMANOVA
  AIC_like <- n + log(SSR_resid / n) + 2 + df_model
  
  model_df <- model_df %>%
    filter(term != "Total") %>%
    mutate(SSR = SumOfSqs, AIC_approx = AIC_like) %>%
    select(variable, term, Df, SSR, R2, F, `Pr(>F)`, AIC_approx)
  
  results_list[[var_name]] <- model_df
}

# Combine all results into one dataframe
results_with_aic <- bind_rows(results_list)

# Write to CSV
write.csv(results_with_aic, "~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/kraken/aic_principles2-4.csv")
