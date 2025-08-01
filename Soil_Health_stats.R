# Soil chemical parameters
library(dplyr)
library(purrr)
library(tidyverse)
library(patchwork) 



chem <- read_delim("~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/soilChem_bothdatasets_matchesOnly.csv")%>%
  select(4:17)%>%
  na.omit

chem_list <- colnames(chem)


summaries <- map_df(chem_list, function(var) {
  chem %>%
    summarise(
      variable = var,
      mean = mean(.data[[var]], na.rm = TRUE),
      median = median(.data[[var]], na.rm = TRUE),
      sd = sd(.data[[var]], na.rm = TRUE),
      min = min(.data[[var]], na.rm = TRUE),
      q25 = quantile(.data[[var]], 0.25, na.rm = TRUE),
      q75 = quantile(.data[[var]], 0.75, na.rm = TRUE),
      max = max(.data[[var]], na.rm = TRUE)
    )
})

print(summaries)

write.csv(summaries,"~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/Soils data/summaries.csv")

## used the min, quartiles, and max to determine low, medium, and high values
data <- read_delim("~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/Soils data/soil_levels.csv")%>%
  filter(selected == "y")%>%
  select(4:17)%>%
  na.omit(.)

#str(data)
plot_list <- colnames(data)


plots <- map(plot_list, function(var) {

  data %>%
    count(level = .data[[var]]) %>%
    mutate(level = factor(level, levels = c("low", "med", "high"))) %>%  # set order
    ggplot(aes(x = level, y = n, fill = level)) +
    geom_bar(stat = "identity") +
    labs(title = var,
         x = "Level",
         y = "Count") +
    theme_minimal() +
    theme(legend.position = "none") +  
    scale_fill_manual(values = c("low" = "#E41A1C", "med" = "#377EB8", "high" = "#4DAF4A"))
})
plots

(plots[[1]] | plots[[2]]) / (plots[[3]] | plots[[4]])
(plots[[5]] | plots[[6]]) / (plots[[7]] | plots[[8]])
(plots[[9]] | plots[[10]]) / (plots[[11]] | plots[[12]])
(plots[[13]] | plots[[14]]) 


## how do the soil chemical parameters differ by key regen practices?

metadata<-read_delim("~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/metadata_server(producer info).csv")%>%
  #filter(grepl("IN-RICHES", sample))%>%
  #filter(seq == "y") #%>%
  filter(selected == "y")

nrow(metadata)

chem <- read_delim("~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/soilChem_bothdatasets_matchesOnly.csv")%>%
  filter(sample %in% metadata$sample)%>%
  na.omit(.)

names(chem)
#examine outliers (I am trying to remove as few as possible as the replication is low)

q1 <- quantile(chem$EC, 0.25)
q3 <- quantile(chem$EC, 0.75)
iqr <- IQR(chem$EC)

# Define outlier thresholds
lower_bound <- q1 - 1.75 * iqr
upper_bound <- q3 + 1.75 * iqr

# Identify outliers using filter()
outliers <- chem %>%
  filter(EC < lower_bound | EC > upper_bound)

print(outliers$sample)


to_remove<-c("IN-RICHES_139_NM_Leyendecker_CTNCNB_p1_r1_0",
"IN-RICHES_143_NM_Leyendecker_CTMCNB_p1_r1_0",
"IN-RICHES_151_NM_Leyendecker_NTNCNB_p1_r1_0",
"CDA_0524_WestOteroCD_23764",                
"CDA_0524_UteMtnUteTribeCD_23590",
"CDA_0524_StaceDavis_23318",                
"CDA_0524_ShavanoCD_23627",                  
"CDA_0524_ShavanoCD_23307",                   
"CDA_0524_ShavanoCD_23089",                  
"CDA_0524_northeastProwersCD_23195",          
"CDA_0524_MoscaHooperCD_23963",                
"CDA_0524_MoscaHooperCD_23297",               
"CDA_0524_MesaCD_23783",                       
"CDA_0524_MesaCD_23524",                      
"CDA_0524_MesaCD_23348",
"CDA_0524_MancosCD_23786",
"CDA_0524_MancosCD_23424",
"CDA_0524_EastOteroCD_23561",                 
"CDA_0524_CACDEE_23420",
"CDA_0524_BentCoCD_23968",                    
"CDA_0524_BentCoCD_23231")

data <- chem%>%
  inner_join(metadata, by ="sample")%>%
 filter(!sample %in% to_remove)

nrow(data)

names(data)

#shapiro.test(data$per_Total_C)
test <- kruskal.test(EC~ TreeShrub, data = data)

pvals <- c(0.000254,0.9472, 0.00337,0.02617,0.3663)
fdr_pvals <- p.adjust(pvals, method = "fdr")

#stats as a group
variables <- c("%_Total_N", "%_Total_C", "Potentially__Mineralizable_N_(mg/kg)", "EC", "MBC(mg/kg)", "%clay")

kruskal_results <- lapply(variables, function(var) {
  # Add backticks around variable name
  formula <- as.formula(paste("Minimize_Soil_Disturbance ~", paste0("`", var, "`")))
  
  test <- kruskal.test(formula, data = data)
  
  data.frame(
    Variable = var,
    Statistic = as.numeric(test$statistic),
    P_Value = as.numeric(test$p.value),
    stringsAsFactors = FALSE
  )
})

# Combine list into dataframe
kruskal_df <- do.call(rbind, kruskal_results)

# Add FDR-corrected p-values
kruskal_df$P_FDR <- p.adjust(kruskal_df$P_Value, method = "fdr")


# View results
print(kruskal_df)
# Get mean Total N for each irrigation type
aggregate(per_Total_N ~ Irrigated, data = data, mean)


#box plots

grob <- grobTree(textGrob("FDR p_adj: 0.044", x=0.025,  y=0.90, hjust=0,
                          gp=gpar(col="black", fontsize=12)))

plot <- ggplot(data, aes(x = Mulching, y =EC,  fill = Mulching)) +
  geom_boxplot() +
  # facet_wrap(~Type, scales="free") +
  theme_bw() +
  annotation_custom(grob)+
  theme(
    axis.text.x = element_text(angle = 0, size = 12),      # tick labels on x-axis
    axis.text.y = element_text(size = 12),                 # tick labels on y-axis
    axis.title.x = element_text(size = 14, face = "bold"), # x-axis title
    axis.title.y = element_text(size = 14, face = "bold"), # y-axis title
    legend.position = "none"
  )
plot



###correlogram
library(corrgram)
library(ggcorrplot)
library(reshape2)

metadata<- read_delim("~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/soilChem_bothdatasets_matchesOnly.csv")%>%
  select(4:17)%>%
  na.omit
#basic
#correlation_matrix <- cor(metadata, method = "spearman", use = "pairwise.complete.obs")
#write.csv(correlation_matrix,"~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/correlations.csv")
#ggcorrplot(correlation_matrix,
#           type = "upper", 
#           outline.color = "black",
#           colors = c("blue","white","red"),
#           lab_size = 2,
#           ggtheme = theme_classic()) + 
#  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
#  theme(axis.text.x = element_text(size = 8),  # Smaller x-axis text
#        axis.text.y = element_text(size = 8))

# advanced
cor.mtest <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], method = "spearman", conf.level = conf.level)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

p.mat <- cor.mtest(metadata)

corr_melt <- melt(correlation_matrix)
colnames(corr_melt) <- c("Var1", "Var2", "correlation")

# Melt p-value matrix
pval_melt <- melt(p.mat)
colnames(pval_melt) <- c("Var1", "Var2", "pvalue")

corr_melt$stars <- ifelse(pval_melt$pvalue < 0.001, "***",
                          ifelse(pval_melt$pvalue < 0.01, "**",
                                 ifelse(pval_melt$pvalue < 0.05, "*", "")))

corr_melt$label <- paste0(round(corr_melt$correlation, 2), corr_melt$stars)

corr_melt <- corr_melt[as.numeric(corr_melt$Var1) < as.numeric(corr_melt$Var2), ]
write.csv(corr_melt,"~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/correlations_pvalues.csv")

ggplot(corr_melt, aes(x = Var2, y = Var1, fill = correlation)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limit = c(-1, 1)) +
  geom_text(aes(label = label), size = 3, color = "black") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8)) +
  labs(fill = "Spearman\nCorrelation")
