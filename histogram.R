##Histogram of all chemical data
library(readr)
library(dplyr)
library(ggplot2)
library(gridExtra)

# Load and filter your data
data <- read_delim(
  "~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/research_SoilChem.csv"
) %>%
  filter(Microbiome_sample == "Y") %>%
  select(-Texture)  # removes the Texture column
names(data)
# Clean up column names: trim whitespace
colnames(data) <- trimws(colnames(data))

# Define your lists of variables
bioVars <- c("MBC(mg/kg)",  "DOC(unfumigated)(mg/kg)", 
             "Potentially  Mineralizable N (mg/kg)", "BetaG")


chemVars <- c("pH", "EC","% Total C", "% Total N","% Inorganic Carbon",
              "% Organic Carbon")
              

chemVars2 <- c("Ag 328.068\n(mg/L)", "Al 396.153\n(mg/L)", "Ba 233.527\n(mg/L)",
              "Ca 317.933\n(mg/L)", "Cd 228.802\n(mg/L)", "Cr 267.716\n(mg/L)",
              "Cu 327.393\n(mg/L)", "Fe 238.204\n(mg/L)", "K Radial\n(mg/L)",
              "Mg 285.213\n(mg/L)", "Mn 257.610\n(mg/L)", "Mo 202.031\n(mg/L)",
              "Na 589.592\n(mg/L)", "Ni 231.604\n(mg/L)",
              "P 213.617\n(mg/L)", "P 214.914\n(mg/L)", "P 178.221\n(mg/L)",
              "P 177.434\n(mg/L)", "Ti\n(mg/L)", "Zn 206.200\n(mg/L)",
              "Pb 220.353\n(mg/L)", "B 249.677\n(mg/L)", "Be 313.107\n(mg/L)",
              "Sr 407.771\n(mg/L)")

physVars <- c( "%silt", "%sand","%clay", "8mm dry plant mass (g)", "8mm dry rock mass (g)",
               "8mm dry bag mass (g)", "Soil Moisture", "Bulk Density 8mm mass")
 

plot_list <- list()
for (var in physVars) {
  # Check if the column exists in the dataset
  if (!(var %in% names(data))) {
    warning(paste("Column", var, "not found in data"))
    next
  }
  
  # Create the aesthetic mapping
  mapping <- aes(x = !!sym(var))
  
  # Create the plot
  p <- ggplot(data, mapping)
  
  # Determine the number of unique, finite values
  x_data <- data[[var]]
  uni <- length(unique(x_data[is.finite(x_data)]))
  
  # Choose the appropriate geom based on the number of unique values
  if (uni > 20) {
    p <- p + geom_histogram(bins = 30, fill = "steelblue", color = "black")
  } else {
    p <- p + geom_bar(stat = "count", fill = "steelblue", color = "black")
  }
  
  # Add labels and theme
  p <- p + labs(title = var, x = var, y = "Count/Value") + theme_minimal()
  
  # Store the plot in the list
  plot_list[[var]] <- p
}

# Display the plots
grid.arrange(grobs = plot_list, ncol = 4)
