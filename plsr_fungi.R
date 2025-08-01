###########PLS for fungi/kraken output
#MB script, LM edits
#interpretting results: 

##### Load relevant libraries 
library(ggplot2)
library(pls)
library(tidyverse)
library(dplyr)
source("~/Documents/rscripts/VIP.R") ## This is a custom script (from the internet) to compute VIPs, has to be in the same directory as this script (or this path has to be changed)

metadata<-read_delim("~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/soilChem_bothdatasets_matchesOnly.csv")
######PLS wants the samples as rows####
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


# Get matching metadata
metadata_filtered <- metadata %>%
  filter(selected == "y", sample != "CDA_0524_COAllianceCenterEE_23206") %>%
  select(1,4:17) %>%
  filter(!is.na(sample))


otu_data_filtered <- otu_rel_cov %>%
  filter(sample %in% metadata_filtered$sample)


metadata_filtered <- metadata_filtered %>%
  filter(sample %in% otu_data_filtered$sample) %>%
  na.omit()

otu_data_filtered <- otu_data_filtered %>%
  filter(sample %in% metadata_filtered$sample)

OTU_table <- otu_data_filtered %>% column_to_rownames("sample") %>% as.matrix()
meta <- metadata_filtered %>% column_to_rownames("sample") %>% as.matrix()

#make SURE everything is in the same order and that (in reorganizing) nothing was lost
all(rownames(OTU_table) == rownames(meta))  # Should return TRUE
meta <- meta[rownames(OTU_table), ] #If false, run this and then check

dim(OTU_table)
dim(meta)

#build a Partial Least Squares regression mode for each variable and determine maximum R2 using leave one out cross validation
th_r2<-0.1 # We will only look at the PLS if the correlation is better than 0.1

for (i in 1:ncol(meta)){ # We treat each variable independently
  parameter<-colnames(meta)[i]
  obs_values<-meta[,i] # these are the observed values we'll try to predict
  print(paste("Trying to predict ",parameter," --- ",i,sep=""))
  # We perform the sPLS, trying to predict our specific metabolite vector (metabolite[,i]) using our whole OTU table.
  #Validation is LOO, so "Leave one out", i.e. we train model on n-1 samples and try to predict the value for the remaining one. 
  #the "method" argument is to chose the correct type of sPLS
  pls_result<-plsr(obs_values ~ OTU_table, validation="LOO",method="oscorespls") 
  # Now we check the vector of r2 
  #(sPLS tries to use different numbers of OTUs and provides a correlation between predicted and observed for each of them, 
  #so we get a vectore of r2 and not just one r2 value)
  r2_vector<-R2(pls_result)
  max<-0
  max_comp<--1
  for (j in 1:length(r2_vector$val)){
    if(!(is.na(r2_vector$val[j]))){
      if(r2_vector$val[j]>th_r2){
        if(r2_vector$val[j]>max){
          max<-r2_vector$val[j]
          max_comp<-r2_vector$comp[j]
        }
      }
    }
  }
  print(paste(" the max r2 is ",max," corresponding to comp ",max_comp,sep=""))
  # So here we print the highest r2 across all predictions
}

## So now we can look at the metadata, and see if any can be predicted by the OTUs abundance
## Looks like pH is not bad (r2 0.60) so we set i to the corresponding column (25)

i<-2
# And we regenerate the corresponding results
parameter<-colnames(meta)[i]
obs_values<-meta[,i] # these are the observed values we'll try to predict
pls_result<-plsr(obs_values ~ OTU_table, validation="LOO",method="oscorespls") 
r2_vector<-R2(pls_result)
max<-0
max_comp<--1

for (j in 1:length(r2_vector$val)){
  if(!(is.na(r2_vector$val[j]))){
    if(r2_vector$val[j]>th_r2){
      if(r2_vector$val[j]>max){
        max<-r2_vector$val[j]
        max_comp<-r2_vector$comp[j]
      }
    }
  }
}


# Plotting predicted vs observed
df<-data.frame(x=obs_values,y=pls_result$validation$pred[,,max_comp])
colnames(df)<-c("x","y")
#pdf(paste("measured_vs_predicted_",module,"-vs-",parameter,".pdf"))
ggplot(data=df) + geom_point(aes(x=x,y=y)) + 
  geom_smooth(aes(x=x,y=y),method=lm) + xlab("Measured") + 
  ylab("Predicted") +
  ggtitle(paste("Comparison of ",parameter," measured vs predicted -- r2=",max)) + 
  theme(axis.text=element_text(color="black",size=8),axis.ticks=element_line(color="black"))
#dev.off()

## So if r2 is not "great", but there seems to be some common features to all the "high-sulfate" samples (the two above 300), 
#as well as some of the "medium sulfate" (on the right side of the plot, > 300 but predicted around 100-200),
#while all the "low-sulfate" samples (bottom left, <100) are also well predicted.
# So next we checking the VIP (variable importance in projection), and output a table of the 100 highest values, 
#this will tell us which OTU you need to know the abundance of to correctly predict the feature of interest

output<-paste("VIP_values_with_",parameter,".csv",sep="")
cat("Rank,OTU,VIP\n",file = output,append=FALSE)
vip_result<-VIP(pls_result)
vip_components<-sort(vip_result[max_comp,],decreasing=TRUE)[1:100]
for (k in 1:100){
  cat(paste(k,names(vip_components[k]),vip_components[k],"\n",sep=","),file=output,append=TRUE)
}

write.csv(vip_components, "~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/fungi_EC_predictors.csv")

## Check the correlation between predicted and modeled (should be consistent with what plsr gave us)
cor.test(df$x,df$y)
## Alternatively, we can also use the built-in function "predplot" but I find it less pretty
## Can be good to double check the ggplot2 plot though (should be the same)
predplot(pls_result,ncomp=max_comp)

## Now we can also plot individual OTUs vs X on the VIP list
for (k in 1:5){
  OTU<-unlist(names(vip_components[k]))
  vec_OTU<-OTU_table[,OTU]
  df<-data.frame(x=vec_OTU,y=obs_values)
  print(ggplot(data=df) + geom_point(aes(x=x,y=y)) + xlab(OTU) + ylab(paste("Measured ",parameter,sep="")) 
        + ggtitle(paste("Comparison of ",OTU," vs measured ",parameter)) 
        + theme(axis.text=element_text(color="black",size=10),axis.ticks=element_line(color="black")))
}


#Plot linear regression of  abundance of individual predictive orgs by individual soil chem
#pick taxa and variable
taxon_to_plot <- "p__Ascomycota;c__Eurotiomycetes;o__Onygenales;f__Onygenaceae;g__Amauroascus;s__Amauroascus_niger"
soil_var <- "EC"

# Turn OTU_table into dataframe and keep sample names
otu_df <- as.data.frame(OTU_table) %>%
  rownames_to_column("sample")

# Same for meta
meta_df <- as.data.frame(meta) %>%
  rownames_to_column("sample")
# Join them by sample
plot_df <- left_join(otu_df, meta_df, by = "sample")

# Make a new column with a short clean name
plot_df$taxon_abundance <- plot_df[[ taxon_to_plot ]]

# Fit the model
model <- lm(taxon_abundance ~ plot_df[[soil_var]], data=plot_df)
r2_value <- summary(model)$r.squared
r2_label <- paste0("R² = ", round(r2_value, 2))

# Check summary
summary(model)

# Plot with ggplot, now using the clean column name
ggplot(plot_df, aes_string(x=soil_var, y="taxon_abundance")) +
  geom_point(color="steelblue", alpha=0.7, size=3) +
  geom_smooth(method="lm", color="black", se=FALSE) +
  annotate("text", 
           x=Inf, y=Inf, 
           label=r2_label, 
           hjust=1.1, vjust=1.5, 
           size=5, fontface="bold") +
  labs(
    title=paste("Abundance of Amauroascus_niger vs", soil_var),
    y="Relative abundance",
    x=soil_var
  ) +
  theme_minimal()


# Plot VIP scores of all predictors against individual variables and size points by abundance
predictor <- read_delim("~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/singleM/totalC_predictor.csv")


taxa_of_interest <- predictor$org
mean_abundances <- sapply(taxa_of_interest, function(taxon) {
  # use backticks in case names have special chars
  if(taxon %in% colnames(otu_df)) {
    mean(otu_df[[taxon]], na.rm=TRUE)
  } else {
    NA
  }
})
predictor$mean_abundance <- mean_abundances

correlations <- sapply(taxa_of_interest, function(taxon) {
  if(taxon %in% colnames(otu_df)) {
    # merge with metadata
    temp_df <- data.frame(
      abundance = otu_df[[taxon]],
      soil = meta_df[[soil_var]]
    )
    cor(temp_df$abundance, temp_df$soil, use="complete.obs")
  } else {
    NA
  }
})
predictor$cor_with_soil <- correlations

#Color by genus
library(stringr)

predictor$class <- str_extract(predictor$org, "c__[^;]+")
predictor$class <- str_remove(predictor$genus, "c__")  # remove prefix, keep clean genus name

ggplot(predictor, aes(x=cor_with_soil, y=vip, size=mean_abundance, color=genus)) +
  geom_point(alpha=0.8) +
  labs(
    x = paste("Correlation with", soil_var),
    y = "VIP score",
    size = "Mean abundance",
    color = "Genus",
    title = paste("VIP vs correlation with", soil_var)
  ) +
  theme_minimal() +
  theme(
    legend.position = "right"
  )

#Predictor abubundance x soil var and set point size to VIP score
#orgs of interest
# concatenate the vip values csvs generated for each var
predictors <- read_delim("~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/fungi_EC_predictors.csv")
colnames(predictors) <- c("org", "vip")

predictors <- predictors %>%
  mutate(
    class = str_extract(org, "c__[^;]+") %>%  # Extract family like "f__QHBO01"
      str_remove("^c__")                      # Remove the prefix "f__"
  )

head(predictors)
#OTU table
otu_data <- otu_rel_cov %>%
  filter(sample %in% metadata$sample)
head(otu_data)

#soil var data
meta <- metadata %>%
  select(1,4:17)%>%
  filter(sample %in% otu_data$sample)%>%
  na.omit(.)
head(meta)


# Filter otu_data to just the columns of organisms in predictor
otu_sub <- otu_data %>%
  select(sample, all_of(predictors$org))

# Pivot longer so each row = one sample + organism abundance
otu_long <- otu_sub %>%
  pivot_longer(
    cols = -sample,
    names_to = "org",
    values_to = "abundance"
  )

# Join with VIP scores and filter scores greater than 8
plot_data <- otu_long %>%
  inner_join(predictors, by = "org") %>%
  inner_join(meta, by = "sample")%>%
  filter(vip >= 5)

# Choose soil variable
soil_var <- "EC"

# Plot
ggplot(plot_data, aes(x = .data[[soil_var]], y = abundance)) +
  geom_point(aes(size = vip, color = org), alpha = 0.7) +
  scale_size_continuous(name = "VIP Score") +
  labs(x = soil_var, y = "Abundance", color = "genus") +
  theme_minimal() +
  theme(legend.position = "right") +
  ggtitle("Abundance of predictive organisms vs EC")

