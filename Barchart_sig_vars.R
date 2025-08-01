##Bar plot to show significant variables from permanova, envfit, mantel etc

data <- read_delim("~/Documents/Wrighton_lab/Agribiome/AG ML project/data processing/kraken/practices_forbar.csv") 

perm <- data %>%
  filter(term != "Residual") %>%
  filter(term!= "Total")

perm$Group <- ifelse(perm$`Pr(>F)` <= 0.05, "Significant", "Not Significant")

# Ensure 'Significant' is ordered before 'Not Significant'
perm$Group <- factor(perm$Group, levels = c("Significant", "Not Significant"))

# Step 2: Order the data first by Group (Significant first) and then by R2 within each group, descending
perm <- perm[with(perm, order(Group, -R2)), ]

# Step 3: Create an ordered factor for plotting
perm$name_Ordered <- factor(perm$term, levels = unique(perm$term))


#With legend
permutation <- ggplot(perm, aes(x = name_Ordered, y = R2, fill = Group)) + 
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste("P =", sprintf("%.3f", `Pr(>F)`))), 
            vjust = -0.3, hjust = 0.5, size = 3) +
  scale_fill_manual(values = c("Significant" = "#b35", "Not Significant" = "#D6D6D6"),
  #scale_fill_manual(values = c("Significant" = "darkblue", "Not Significant" = "#D6D6D6"),
 # scale_fill_manual(values = c("Significant" = "lightblue", "Not Significant" = "#D6D6D6"),
         name = "",
                    labels = c("Significant", "Not Significant")) +
  theme_minimal() +
  labs(title = "",
       x = "",
       y = "R2") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        legend.position = "none")

permutation




