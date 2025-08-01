##analyze Kraken outputs
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(stringr)


#braken  - take a look
data <- read_delim("~/Documents/Wrighton_lab/Spring 2025/metaG workshop/data processing/NT30_bracken_out_order.txt")
head(data)
data <- data %>% arrange(desc(new_est_reads))

barplot(data$new_est_reads, names.arg=data$name, las=2, cex.names=0.7, 
        main = "No till 30 years")

# otu table and stats
## upload and data wrangling
### upload "longer" format from concatenating the krona otu files from Bracken
data_longer <- read_delim("~/Documents/Wrighton_lab/Spring 2025/metaG workshop/otu_table_long.txt", 
                   col_names = FALSE, delim = "\t")  

colnames(data_longer) <- c("Sample", "Count", "Taxonomy")

data_longer <- data_longer %>%
  mutate(Taxonomy = str_replace_all(Taxonomy, "\t", ";"))

head(data_longer)


### take a look at the OTU table - you might not actually use this file, but we want to know if all the samples are there and the counts are appropriately arranged
data_wider <- data_longer %>%
  pivot_wider(names_from = Sample, values_from = Count, values_fill = list(Count = 0))


## making a rank abundance curve using "data_longer"
data <- data_longer %>%
  separate(Taxonomy, into = c(
    "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
    sep = ";", fill = "right")%>%
  filter(Sample=="NT10")

# Summarize read counts per order
data_summary <- data %>%
  mutate(Genus = str_trim(Genus))%>%
  group_by(Sample, Genus) %>%
  summarise(total_count = sum(Count), .groups = "drop")


# Convert to relative abundance per sample
data_summary <- data_summary %>%
  group_by(Sample) %>%
  mutate(relative_count = total_count / sum(total_count) * 100)

rank_abundance <- data_summary %>%
  group_by(Genus) %>%
  summarise(total_count = sum(total_count, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_count)) %>%
  mutate(rank = row_number()) %>%
  mutate(normalized_abundance = total_count / sum(total_count))


# Rank Curve 
rank_abundance_plot <- ggplot(rank_abundance, aes(x = rank, y = normalized_abundance, label = Genus)) +
  geom_bar(stat = "identity", aes(fill = Genus), width = 0.7) + 
  geom_text(aes(y = normalized_abundance * 1.05), angle = 45, hjust = 0, size = 3) + 
  # scale_fill_manual(values = random_colors) + 
  labs(title = "Rank Abundance Curve -Josh NT10",
       x = "Rank (Most to Least Abundant Taxa)",
       y = "Relative Abundance",
       fill = "Genus") +
  theme_minimal() + theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

rank_abundance_plot


# rank curve top 20
rank_abundance_top20 <- rank_abundance %>% 
  filter(rank <= 20)

rank_abundance_plot <- ggplot(rank_abundance_top20, aes(x = rank, y = normalized_abundance, label = Genus)) +
  geom_bar(stat = "identity", aes(fill = Genus), width = 0.7) + 
  geom_text(aes(y = normalized_abundance * 1.05), angle = 45, hjust = 0, size = 3) +
  labs(title = "Rank Abundance Curve - Josh NT10",
       x = "Rank (Most to Least Abundant Taxa)",
       y = "Relative Abundance",
       fill = "Genus") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

rank_abundance_plot


## color by phyla
data_summary <- data %>%
  mutate(Genus = str_trim(Genus)) %>%
  group_by(Sample, Genus, Phylum) %>%
  summarise(total_count = sum(Count), .groups = "drop")

rank_abundance <- data_summary %>%
  group_by(Genus, Phylum) %>%
  summarise(total_count = sum(total_count, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_count)) %>%
  mutate(rank = row_number(),
         normalized_abundance = total_count / sum(total_count, na.rm = TRUE))

rank_abundance_plot <- ggplot(rank_abundance, aes(x = rank, y = normalized_abundance, fill = Phylum)) +
  geom_bar(stat = "identity", width = 0.7) + 
  geom_text(aes(label = Genus, y = normalized_abundance * 1.05),
            angle = 45, hjust = 0, size = 3) +
  labs(title = "Rank Abundance Curve - Josh NT10",
       x = "Rank (Most to Least Abundant Taxa)",
       y = "Relative Abundance",
       fill = "Phylum") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


