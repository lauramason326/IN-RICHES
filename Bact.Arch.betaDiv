##singleM - beta Div
library(vegan)
library(tidyverse)
library(grid)
library(forcats)
library(dplyr)

## add metadata so that you can grab the relevant samples first
metadata<-read_delim("~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/metadata_server(producer info).csv")%>%
  filter(selected=="y")%>%
  filter(seq == "y")%>%
  filter(sample != "CDA_0524_LarimerCD_23992")

#if you are using envfit below, add that data and use it to filter the samples. 
chem <- read_delim("~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/soilChem_bothdatasets_matchesOnly.csv")%>%
  select(1,4:17)%>%
  filter(sample != "CDA_0524_LarimerCD_23992")%>%
  filter(sample %in% metadata$sample)%>%
na.omit(.)


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
  filter(sample != "CDA_0524_LarimerCD_23992")

nrow(otu_table)

otu_rel_cov <- otu_table %>%
  rowwise() %>%
  mutate(total_coverage = sum(c_across(-sample))) %>%
  ungroup() %>%
  mutate(across(-c(sample, total_coverage), ~ .x / total_coverage)) %>%
  dplyr::select(-total_coverage) %>%
#  filter(sample %in% chem $sample)%>%
  na.omit()

nrow(otu_rel_cov)
nrow(metadata)

#wrangle
meta <- metadata%>% 
  inner_join(otu_rel_cov, by="sample")%>%
 filter(sample %in% chem$sample)
  
  
nrow(meta)
ncol(meta)
names(meta)

table <-meta%>% 
  dplyr::select(34:4519)

data.t <- as.matrix(vegdist(table, method = 'bray'))


#NMDS

NMDS<-metaMDS(data.t, trymax= 200)
MDS1 = NMDS$points[,1]
MDS2 = NMDS$points[,2]
NMDS1=MDS1
NMDS2=MDS2
NMDS = data.frame(NMDS1,NMDS2)

# Plot
grob <- grobTree(textGrob("Stress0.1673569", x=0.025,  y=0.90, hjust=0,
                         gp=gpar(col="black", fontsize=12)))


NMDS_plot <- ggplot(NMDS) + 
  geom_point(aes(x = NMDS1, y = NMDS2, fill = meta$Irrigated), 
             shape = 21, size = 6, color = "black") +
  labs(x = "NMDS1", y = "NMDS2", fill = "Irrigated") +
  theme_bw() +
  annotation_custom(grob)

NMDS_plot
  
ggsave("~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/both_presribedGrazing.pdf", width=10,height=8,dpi=300)

#ENVFIT
envfit.data <- chem 

# Fit environmental vectors
envfit_perm <- envfit(NMDS, envfit.data, permutations = 999)

# Extract R² and p-values
envfit_r2   <- envfit_perm$vectors$r
envfit_pval <- envfit_perm$vectors$pvals

# Combine into a single dataframe with variable names
env.scores <- data.frame(
  Variable = names(envfit_r2),
  R2       = envfit_r2,
  pval     = envfit_pval
)

write.csv(env.scores,"~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/env_scores_bothDatasets.csv" )

sig.env.scrs <- subset(env.scores, pval<=0.05) #subset data to show variables significant at 0.05

plot_env <- NMDS_plot +
  geom_segment(data = sig.env.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), 
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant env variables
  ggrepel::geom_text_repel(data = sig.env.scrs, aes(x=NMDS1, y=NMDS2, label = env.variables),
                           cex = 4, direction = "both", segment.size = 0.25)#+ #add labels for env variables


plot_env

ggsave("~/Documents/Wrighton_lab/Agribiome/Iowa_Regen/NMDS10_vectors.pdf",width=10,height=8,dpi=300)






