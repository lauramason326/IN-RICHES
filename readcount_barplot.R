##barplot to visualize sequencing results

# Load the data

data <- read_delim("~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/batch5_counts.csv")

library(dplyr)
library(ggplot2)

# 1. Assign fill color for plotting
data <- data %>%
  mutate(
    fill_color = case_when(
      type == "target" ~ "grey",
      type == "actual" & direction == "F" ~ "blue",
      type == "actual" & direction == "R" ~ "red"
    )
  )

# 2. Create a new column that combines Sample and Type
#    This is key to getting two bars per Sample
data <- data %>%
  mutate(
    sample_type = paste(Sample, type, sep = "_")
  )

# 3. Plot with sample_type on x-axis so F/R stack within actual/target
ggplot(data, aes(x = sample_type, y = gbp, fill = fill_color)) +
  geom_bar(stat = "identity") +
  scale_fill_identity(guide = "legend",
                      labels = c("blue" = "Actual - F",
                                 "red" = "Actual - R",
                                 "grey" = "Target (F + R)"),
                      breaks = c("blue", "red", "grey"),
                      name = "Direction") +
  labs(x = "Sample", y = "GBP", title = "Batch 5") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 9))
