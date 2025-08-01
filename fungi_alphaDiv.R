library(tidyverse)
library(vegan)
library(plyr)
library(dplyr)
library(GGally)
library(ecodist)
library(matrixStats)
library(viridis)
library(car)


metadata<-read_delim("~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/metadata_server(producer info).csv")%>%
  filter(selected=="y")%>%
  filter(sample != "CDA_0524_COAllianceCenterEE_23206")

### upload "longer" format from concatenating the krona otu files from Bracken
data_longer <- read_delim(
  "~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/Fungi_otu_jul14.txt", 
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


richness<-  ddply(otu_rel_cov,~sample,function(x) {
  data.frame(Richness=sum(x[-1]>0))
})

evenness <- ddply(otu_rel_cov, ~sample, function(x) {
  simpson_div <- diversity(x[-1], index = "simpson")  # Simpson's diversity
  data.frame(Evenness = simpson_div / log(sum(x[-1] > 0)))  # Evenness
})

shannon<- ddply(otu_rel_cov,~sample,function(x) {
  data.frame(Shannon=diversity(x[-1], index="shannon"))
})


#combine calculations into one 
join <- shannon%>% 
  inner_join(richness) %>% 
  inner_join(evenness) %>% 
  inner_join(metadata) %>% 
  pivot_longer(2:4, names_to = "Type", values_to = "Diversity")

write.csv(join, "~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/fungi_alpha_div.csv")

div<-read_delim("~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/fungi_alpha_div.csv")

# Plotting code

plot <- ggplot(div, aes(x = dataset, y = Diversity, fill = Type)) +
  geom_boxplot() +
  facet_wrap(~Type, scales="free")+
  theme_bw()

plot


#stats
names(div)
div$Type

test<- div%>%
  filter(Type =="Evenness")

shapiro.test(test$Diversity)

summary(richness$Diversity)
summary(evenness$Diversity)
summary(shannon$Diversity)

# Summary for richness
richness %>%
  group_by(Minimize_Soil_Disturbance) %>%
  summarise(
    n = n(),
    mean = mean(Diversity, na.rm = TRUE),
    sd = sd(Diversity, na.rm = TRUE),
    min = min(Diversity, na.rm = TRUE),
    max = max(Diversity, na.rm = TRUE)
  )

# Summary for evenness
evenness %>%
  group_by(TreeShrub) %>%
  summarise(
    n = n(),
    mean = mean(Diversity, na.rm = TRUE),
    sd = sd(Diversity, na.rm = TRUE),
    min = min(Diversity, na.rm = TRUE),
    max = max(Diversity, na.rm = TRUE)
  )%>%
  print(n = Inf)

# Summary for shannon
shannon %>%
  group_by(Minimize_Soil_Disturbance) %>%
  summarise(
    n = n(),
    mean = mean(Diversity, na.rm = TRUE),
    sd = sd(Diversity, na.rm = TRUE),
    min = min(Diversity, na.rm = TRUE),
    max = max(Diversity, na.rm = TRUE)
  )%>%
  print(n = Inf)


leveneTest(evenness$Diversity, evenness$dataset)

#syntax: kruskal.test(PHOS_nmol.h.g ~ Season, data = AU)
names(metadata)


variables <- c( "dataset","locationName","Cropping_System","Organic","Irrigated",               
                "Conservation_Cover","Conservation_Crop_Rotation","No_till","Soil_Carbon_Amendment","Cover_Crop", "Reduced_Till" ,
                "Mulching","Forage_Harvest_Management","Prescribed_Grazing","Row_Arrangement",
                "Nutrient_Management","TreeShrub","Soil_Armor","Livestock_Integration",
                "Continuous_Living_Root","Plant_Diversity","Minimize_Soil_Disturbance" )

# Run Kruskal-Wallis tests and store results
test<- div%>%
  filter(Type =="Shannon")

kruskal_results <- lapply(variables, function(var) {
  formula <- as.formula(paste(var, "~ Diversity"))
  test <- kruskal.test(formula, data = test)
  data.frame(
    Variable = var,
    Statistic = test$statistic,
    P_Value = test$p.value,
    stringsAsFactors = FALSE
  )
})

# Combine all results into one data frame
kruskal_df <- do.call(rbind, kruskal_results)
kruskal_df 


# Write to CSV
write.csv(kruskal_df, "~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/fungal_kruskal_results_shannon_bothDatasets.csv", row.names = FALSE)



##stats with outlier removal



