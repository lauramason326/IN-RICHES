library(tidyverse)
library(dplyr)
library(ComplexHeatmap)

data <- read.csv("~/Documents/Wrighton_lab/Agribiome/AG ML project/metadata_for_upset.csv", header = TRUE) %>%
  filter(seq =="FALSE")%>%
select(3:33)
names(data)

m <- make_comb_mat(data)

ncol(data)

UpSet(m,
      pt_size = unit(3, 'mm'),
      lwd = 1, 
      set_order = c("program",                        
                   "research",  
                   "Unknown",                       
                   "Cropping",                       
                   "Grazing",   
                    "Alfalfa",                        
                    "Annual.Grasses..Grains..Legumes",
                     "Grazing.Lands",                 
                    "Mixed.Species.Hay",              
                    "Orchard.Vineyard",              
                    "Root.Crops",                    
                    "Specialty.Row.Crops",            
                    "Organic",                        
                    "Irrigated",                      
                    "Conservation_Cover",            
                    "Conservation_Crop_Rotation",     
                    "No_till",                        
                    "Soil_Carbon_Amendment",          
                    "Cover_Crop",                    
                    "Reduced_Till",                   
                    "Mulching",                       
                     "Forage_Harvest_Management",      
                     "Prescribed_Grazing",            
                    "Row_Arrangement",                
                    "Nutrient_Management",            
                    "TreeShrub",                      
                    "Soil_Armor",                    
                    "Livestock_Integration",         
                    "Continuous_Living_Root",         
                    "Plant_Diversity",                
                    "Minimize_Soil_Disturbance"),
      column_names_gp = gpar(fontsize = 10),  # smaller column names
      row_names_gp = gpar(fontsize =10)      # smaller row names
)
