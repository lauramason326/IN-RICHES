##singleM - beta Div
library(vegan)
library(tidyverse)
library(grid)
library(forcats)
library(dplyr)
library(purrr)
library(tibble)


metadata<-read_delim("~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/metadata_server(producer info).csv")%>%
  #filter(grepl("IN-RICHES", sample))%>%
  filter(selected=="y")


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


nrow(otu_rel_cov)

ncol(metadata)
nrow(metadata)
names(metadata)

#wrangle
#there is an na in the conservation crop rotation that I just cannot seem to get rid of
meta <- metadata%>% 
  inner_join(otu_rel_cov, by="sample")

# Samples in otu_table
samples_otu <- otu_table$sample
# Samples in metadata
samples_meta <- meta$sample
# Samples in otu_table but NOT in metadata
missing_in_meta <- setdiff(samples_otu, samples_meta)
print(missing_in_meta)


names(meta)
ncol(meta)
nrow(meta)

#write.csv(meta,"~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/singleM/meta_first_pass.csv" )
table <-meta%>% 
  select(34:4519) 

set.seed(1234)
data.t <- as.matrix(vegdist(table, method = 'bray'))
names(metadata)

### just regen practices / principles
# first make sure all vars have more than 1 level


vars <- c("Organic", "Irrigated", "Conservation_Cover", "Conservation_Crop_Rotation", 
          "No_till", "Soil_Carbon_Amendment", "Cover_Crop", "Reduced_Till", 
          "Mulching", "Forage_Harvest_Management", "Prescribed_Grazing", 
          "Nutrient_Management", "TreeShrub")

# Count number of unique levels in each
sapply(meta[, vars], function(x) length(unique(x)))


practices_all <- adonis2(data.t~Organic+Irrigated + Conservation_Cover + 
                           Conservation_Crop_Rotation + No_till + Soil_Carbon_Amendment+
                           Cover_Crop + Reduced_Till + Prescribed_Grazing + Mulching + 
                           Nutrient_Management + TreeShrub, data = meta, by = "terms",strata = meta$locationName, permutations=999)

practices2 <- adonis2(data.t~ Irrigated + Conservation_Crop_Rotation+ No_till+ Reduced_Till + Nutrient_Management, data = meta, by = "terms")
principles = adonis2(data.t ~ Soil_Armor+Livestock_Integration+Continuous_Living_Root+Plant_Diversity+Minimize_Soil_Disturbance+Cropping_System, data= meta, by ="terms", strata = meta$locationName)

adonis_formulas <- list(
  practices7 = data.t ~ Irrigated + Conservation_Crop_Rotation+ No_till,
  practices8 = data.t ~ Irrigated + No_till+ Reduced_Till,
  practices9 = data.t ~ Irrigated + Reduced_Till+ Nutrient_Management,
  practices10 = data.t ~ Irrigated + Nutrient_Management+ Conservation_Crop_Rotation)

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
  AIC_like <- n * log(SSR_resid / n) + 2 * df_model
  
  model_df <- model_df %>%
    filter(term != "Total") %>%
    mutate(SSR = SumOfSqs, AIC_approx = AIC_like) %>%
    select(variable, term, Df, SSR, R2, F, `Pr(>F)`, AIC_approx)
  
  results_list[[var_name]] <- model_df
}

# Combine all results into one dataframe
results_with_aic <- bind_rows(results_list)

# Write to CSV
write.csv(results_with_aic, "~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/singleM/aic_practices_3x.csv")



#PERMANOVA - for the first pass, I am keeping this really simple

adonis_formulas <- list(
  cropping_system_results = data.t ~ Cropping_System,
 # locationName_results = data.t ~ locationName,
  Irrigation_results = data.t ~ Irrigated,
#organic_nutrient_mangement_results = data.t ~ Organic,
  Soil_Armor_results = data.t ~ `Soil_Armor`,
  Livestock_Integration_results = data.t ~ `Livestock_Integration`,
 Plant_Diversity_results = data.t ~ `Plant_Diversity`,
  Minimize_Soil_Disturbance_results = data.t ~ `Minimize_Soil_Disturbance`,
  Conservation_Cover_results = data.t ~ `Conservation_Cover`,
 Conservation_Crop_Rotation_results = data.t ~ `Conservation_Crop_Rotation`,
  no_till_results = data.t ~ `No_till`,
  Soil_Carbon_Amendment_results = data.t ~ Soil_Carbon_Amendment,
  Cover_Crop_results = data.t ~ Cover_Crop,
  Reduced_Till_results = data.t ~ `Reduced_Till`,
  #Mulching_results = data.t ~ Mulching,
  #Forage_Harvest_Management_results = data.t ~ Forage_Harvest_Management,
  Prescribed_Grazing_results = data.t ~ `Prescribed_Grazing`
  #Row_Arrangement_results = data.t ~ Row_Arrangement,
 #Nutrient_Management_results = data.t ~ Nutrient_Management,
 # TreeShrub_results = data.t ~ TreeShrub
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
  AIC_like <- n * log(SSR_resid / n) + 2 * df_model
  
  model_df %>%
    filter(term != "Total") %>%
    mutate(SSR = SumOfSqs, AIC_approx = AIC_like) %>%
    select(variable, term, Df, SSR, R2, F, `Pr(>F)`, AIC_approx)
})
#print results to csv
write.csv(results_with_aic,"~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/singleM/research_locationName_strata.csv")

## permanova  - all metrics
all_var_permanova = adonis2(data.t ~ Cropping_System+
                            locationName+
                            Irrigated+
                            Prescribed_Grazing+
                            Nutrient_Management+
                            dataset+
                            Reduced_Till+
                            Cover_Crop+
                            Plant_Diversity+
                            No_till+
                            Mulching, data = meta, permutations = 999, by = "terms")

two_var_permanova = adonis2(data.t ~ Cropping_System+locationName, data = meta, permutations = 999, by = "terms")

# Run and extract tidy results

results_with_aic <- map_dfr(names(adonis_formulas), function(var_name) {
  # Run adonis2 with correct strata (actual vector, not string)
  model <- adonis2(
    formula = adonis_formulas[[var_name]], 
    data = meta, 
    permutations = 999, 
    by = "terms",
    strata = meta$dataset
  )
  
  model_df <- as.data.frame(model) %>%
    rownames_to_column("term") %>%
    mutate(variable = var_name)
  
  SSR_resid <- model_df %>% filter(term == "Residual") %>% pull(SumOfSqs)
  df_model <- model_df %>% filter(term != "Residual" & term != "Total") %>% pull(Df) %>% sum()
  n <- nrow(meta)
  
  # Calculate approximate AIC for PERMANOVA
  AIC_like <- n * log(SSR_resid / n) + 2 * df_model
  
  model_df %>%
    filter(term != "Total") %>%
    mutate(SSR = SumOfSqs, AIC_approx = AIC_like) %>%
    select(variable, term, Df, SSR, R2, F, `Pr(>F)`, AIC_approx)
})

write.csv(results_with_aic,"~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/singleM/aic_approx_3x_bothDatasets_wstrata_pt2.csv")


##all of them
allvar1 <- adonis2(data.t ~ 
                  Cropping_System +
                  locationName +
                  Nutrient_Management +
                  Mulching +
                  Prescribed_Grazing +
                  locationName:Prescribed_Grazing +
                  locationName:Mulching, data = meta, by="terms", strata = meta$dataset, permutations = 999)

allvar2 <- adonis2(data.t ~ Cropping_System + locationName + Nutrient_Management + Prescribed_Grazing + locationName:Prescribed_Grazing, data = meta, by="terms", strata = meta$dataset, permutations = 999)

allvar3 <- adonis2(data.t ~ 
                     Cropping_System +
                     locationName +
                     Prescribed_Grazing +
                     locationName:Prescribed_Grazing, data = meta, by="terms", strata = meta$dataset, permutations = 999)


adonis_formulas <- list(multivar_good_score = data.t~ locationName*Prescribed_Grazing +locationName*Mulching) 

results_with_aic <- map_dfr(names(adonis_formulas), function(var_name) {
  # Run adonis2 with correct strata (actual vector, not string)
  model <- adonis2(
    formula = adonis_formulas[[var_name]], 
    data = meta, 
    permutations = 999, 
    by = "terms",
    strata = meta$dataset
  )
  
  model_df <- as.data.frame(model) %>%
    rownames_to_column("term") %>%
    mutate(variable = var_name)
  
  SSR_resid <- model_df %>% filter(term == "Residual") %>% pull(SumOfSqs)
  df_model <- model_df %>% filter(term != "Residual" & term != "Total") %>% pull(Df) %>% sum()
  n <- nrow(meta)
  
  # Calculate approximate AIC for PERMANOVA
  AIC_like <- n * log(SSR_resid / n) + 2 * df_model
  
  model_df %>%
    filter(term != "Total") %>%
    mutate(SSR = SumOfSqs, AIC_approx = AIC_like) %>%
    select(variable, term, Df, SSR, R2, F, `Pr(>F)`, AIC_approx)
})

write.csv(results_with_aic,"~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/singleM/aic_approx_allvar.csv")

### trying with location name as the strata
#these are the models with the best aic scores


new1 = adonis2(data.t ~ Prescribed_Grazing + Mulching,
        data = meta, strata = meta$locationName, by = "terms")

new2 = adonis2(data.t ~ Cropping_System * Prescribed_Grazing,
        data = meta, strata = meta$locationName, by = "terms")

new3 = adonis2(data.t ~ Cropping_System + Nutrient_Management + Prescribed_Grazing,
               data = meta, strata = meta$locationName, by = "terms")

new4 = adonis2(data.t ~ Prescribed_Grazing + Cropping_System,
               data = meta, strata = meta$locationName, by = "terms")

new5 = adonis2(data.t ~ Prescribed_Grazing + Nutrient_Management,
        data = meta, strata = meta$locationName, by = "terms")

new6 = adonis2(data.t ~ Cropping_System * Irrigated + Prescribed_Grazing,
               data = meta, strata = meta$locationName, by = "terms")




adonis_formulas <- list(new1 = data.t ~ Prescribed_Grazing + Mulching,
new2 = data.t ~ Cropping_System * Prescribed_Grazing,
new3 = data.t ~ Cropping_System + Nutrient_Management + Prescribed_Grazing,
new4 = data.t ~ Prescribed_Grazing + Cropping_System,
new5 = data.t ~ Prescribed_Grazing + Nutrient_Management,
new6 = data.t ~ Cropping_System * Irrigated + Prescribed_Grazing) 



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
  AIC_like <- n * log(SSR_resid / n) + 2 * df_model
  
  model_df %>%
    filter(term != "Total") %>%
    mutate(SSR = SumOfSqs, AIC_approx = AIC_like) %>%
    select(variable, term, Df, SSR, R2, F, `Pr(>F)`, AIC_approx)
})

write.csv(results_with_aic,"~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/singleM/aic_approx_3x_bothDatasets_wstrata_pt2.csv")


### just regen practices / principles
practices_all <- adonis2(data.t~Organic + Irrigated + Conservation_Cover + Conservation_Crop_Rotation + No_till + Soil_Carbon_Amendment+
  Cover_Crop + Reduced_Till + Mulching +Forage_Harvest_Management + Prescribed_Grazing +
  Nutrient_Management + TreeShrub+Cropping_System, data = meta, by = "terms",strata = meta$locationName)

practices2 <- adonis2(data.t~ Irrigated + Conservation_Crop_Rotation+ No_till+ Reduced_Till + Nutrient_Management, data = meta, by = "terms")
principles = adonis2(data.t ~ Soil_Armor+Livestock_Integration+Continuous_Living_Root+Plant_Diversity+Minimize_Soil_Disturbance+Cropping_System, data= meta, by ="terms", strata = meta$locationName)

adonis_formulas <- list(
  practices7 = data.t ~ Irrigated + Conservation_Crop_Rotation+ No_till,
  practices8 = data.t ~ Irrigated + No_till+ Reduced_Till,
  practices9 = data.t ~ Irrigated + Reduced_Till+ Nutrient_Management,
  practices10 = data.t ~ Irrigated + Nutrient_Management+ Conservation_Crop_Rotation)

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
  AIC_like <- n * log(SSR_resid / n) + 2 * df_model
  
  model_df <- model_df %>%
    filter(term != "Total") %>%
    mutate(SSR = SumOfSqs, AIC_approx = AIC_like) %>%
    select(variable, term, Df, SSR, R2, F, `Pr(>F)`, AIC_approx)
  
  results_list[[var_name]] <- model_df
}

# Combine all results into one dataframe
results_with_aic <- bind_rows(results_list)

# Write to CSV
write.csv(results_with_aic, "~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/singleM/aic_practices_3x.csv")


##principles
adonis_formulas <- list(princ1 = data.t~Soil_Armor,
                        princ2 = data.t~ Livestock_Integration,
                        princ3 = data.t~Continuous_Living_Root,
                        princ4 = data.t~Plant_Diversity,
                        princ5 = data.t~  Minimize_Soil_Disturbance,
                        princ6 = data.t ~ Soil_Armor*Livestock_Integration*Continuous_Living_Root*Plant_Diversity*Minimize_Soil_Disturbance,
                        princ7 = data.t ~ Soil_Armor+Livestock_Integration+Continuous_Living_Root+Plant_Diversity+Minimize_Soil_Disturbance)
  
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
  AIC_like <- n * log(SSR_resid / n) + 2 * df_model
  
  model_df <- model_df %>%
    filter(term != "Total") %>%
    mutate(SSR = SumOfSqs, AIC_approx = AIC_like) %>%
    select(variable, term, Df, SSR, R2, F, `Pr(>F)`, AIC_approx)
  
  results_list[[var_name]] <- model_df
}

# Combine all results into one dataframe
results_with_aic <- bind_rows(results_list)

# Write to CSV
write.csv(results_with_aic, "~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/singleM/aic_principles_prac.csv")
