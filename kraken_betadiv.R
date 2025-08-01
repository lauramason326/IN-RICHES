#Kraken beta div
# upstream processing - this is needed for basically everything 
## upload and data wrangling
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
  mutate(across(where(is.character), ~ gsub("IN_RICHES", "IN-RICHES", .)))

head(data_wider)

ncol(data_wider)
nrow(data_wider)
rownames(data_wider)

#filter out every taxon that appears less than 3 times 
col_totals <- colSums(select(data_wider, -sample), na.rm = TRUE)
# get names of taxa whose total abundance >= 3
taxa_to_keep <- names(col_totals)[col_totals >= 3]
#  keep only 'sample' + these taxon columns
otu_table_filtered <- data_wider %>%
  select(sample, all_of(taxa_to_keep))

otu_table_filtered <- otu_table_filtered %>%
  filter(sample != "CDA_0524_COAllianceCenterEE_23206")

write.csv(otu_table_filtered,"~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/fungi_gt3_jul14.csv")
nrow(data_wider)
nrow(otu_table_filtered)
rownames(otu_table_filtered)

## add metadata so that you can grab the relevant samples 
metadata<-read_delim("~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/metadata_server(producer info).csv")%>%
  filter(selected=="y")%>%
  filter(sample != "CDA_0524_COAllianceCenterEE_23206")

#if you are using envfit below, add that data and use it to filter the samples. 
chem <- read_delim("~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/soilChem_bothdatasets_matchesOnly.csv")%>%
  select(1,4:17)%>%
  filter(sample %in% otu_table_filtered$sample)%>%
  filter(sample != "CDA_0524_COAllianceCenterEE_23206")%>%
  na.omit(.)

nrow(chem)

#convert otu to relative abundance (filtering with chem and/or metadata sample names so that everything matches)
otu_rel_cov <- otu_table_filtered %>%
  rowwise() %>%
  mutate(total_coverage = sum(c_across(-sample))) %>%
  ungroup() %>%
  mutate(across(-c(sample, total_coverage), ~ .x / total_coverage)) %>%
  dplyr::select(-total_coverage) %>%
  filter(sample %in% chem $sample)%>%
  na.omit()


nrow(otu_rel_cov)

#convert metadata to meta and set up for NMDS
#wrangle
meta <- metadata%>% 
  inner_join(otu_rel_cov, by="sample")

nrow(meta)
ncol(meta)
names(meta)

table <-meta%>% 
  dplyr::select(34:706)

data.t <- as.matrix(vegdist(table, method = 'bray'))


#NMDS
NMDS<-metaMDS(data.t, trymax= 200)
MDS1 = NMDS$points[,1]
MDS2 = NMDS$points[,2]
NMDS1=MDS1
NMDS2=MDS2
NMDS = data.frame(NMDS1,NMDS2)

# Plot
grob <- grobTree(textGrob("Stress 0.1665636", x=0.025,  y=0.90, hjust=0,
                          gp=gpar(col="black", fontsize=12)))
names(metadata)

NMDS_plot <- ggplot(NMDS, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(color = meta$Cover_Crop), size = 6) + 
 #geom_text(aes(label = meta$sample), vjust = -1, hjust = 1, check_overlap = TRUE) +
  labs(x = "NMDS1", y = "NMDS2", title = "fungi both datasets") +
  theme_bw() +
  annotation_custom(grob)


NMDS_plot

ggsave("~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/fungi_coverCrop.pdf", width=10,height=8,dpi=300)

#ENVFIT
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

# Save all envfit scores to CSV
write.csv(env.scores, "~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/env_scores_fungi.csv", row.names = FALSE)

# Subset to only significant variables (p <= 0.05)
sig.env.scrs <- subset(env.scores, pval <= 0.099)


plot_env <- NMDS_plot +
  geom_segment(data = sig.env.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), 
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant env variables
  ggrepel::geom_text_repel(data = sig.env.scrs, aes(x=NMDS1, y=NMDS2, label = env.variables),
                           cex = 4, direction = "both", segment.size = 0.25)#+ #add labels for env variables


plot_env

ggsave("~/Documents/Wrighton_lab/Agribiome/Iowa_Regen/NMDS10_vectors.pdf",width=10,height=8,dpi=300)



