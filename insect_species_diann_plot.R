library(tidyverse)
library(RColorBrewer)
library(patchwork) 
file_path <- "D:/data/DIA-NN_result/cricket"

out_put_path <- "D:/data_analysis/BP_insect_species"

df <- read_tsv(file.path(file_path, "report.pg_matrix.tsv"))

df_simplified <- df %>% 
  mutate(main_protein = gsub(";.*", "", Protein.Ids))

df_clean <- df_simplified %>% 
  select(6:18) %>% 
  mutate_all(~replace(., is.na(.), 0)) %>% 
  t() 

colnames(df_clean) <- df_clean[13, ]
df_clean <- df_clean[1:12, ]

df_clean <- as.data.frame(df_clean)


df_clean <- df_clean %>%
  mutate(across(everything(), as.numeric))

p_type <- c('HD', 'HD','HD', 'HD',
            'GW', 'GW','GW','GW',
            'GD', 'GD','GD','GD') %>% 
  as.data.frame()
colnames(p_type) <- "protein_name"

df_clean <- cbind(p_type, df_clean)

anyDuplicated(colnames(df_clean))

colnames(df_clean) <- make.unique(colnames(df_clean))

write_csv(
  df_clean, file = file.path(out_put_path, "cricket_protein.csv")
)

colnames(df_clean) <- make.names(colnames(df_clean))
# pca analysis ------------------------------------------------------------

df_scaled <- scale(df_clean[, 2:1079]) %>% 
  as.data.frame()

pca_result <- prcomp(df_scaled, center = TRUE, scale. = TRUE)

pca_score <- as.data.frame(pca_result$x)

pca_score <- cbind(p_type, pca_score) %>% 
  as.data.frame()

custom_palette <- colorRampPalette(brewer.pal(12, "Paired"))(16)  # Expand Set3 to 20 colors

var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100

pca_plot<- ggplot(
  pca_score, 
  mapping = aes(x = PC1, y = PC2, colour = protein_name, shape = protein_name)
) + 
  geom_point(size = 4, alpha = 0.7) +
  stat_ellipse(level = 0.95, aes(colour = protein_name), alpha = 0.8)+
  labs(title = "PCA score plot of BSF proteins", 
       x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(var_explained[2], 1), "%)")) +
  theme_minimal()+
  scale_fill_manual(values = custom_palette) +
  theme(
    strip.background = element_rect(fill = "white", colour = "white"),
    plot.background = element_rect(fill = "white", colour = "white"),
    axis.title.x = element_text(size = 14, face = "bold"),  # X-axis title size
    axis.title.y = element_text(size = 14, face = "bold"),  # Y-axis title size
    axis.text.x = element_text(size = 10, angle = 0, face = "bold"),                 # X-axis text size
    axis.text.y = element_text(size = 12, face = "bold"),                 # Y-axis text size
    legend.title = element_text(size = 12),                # Legend title size
    legend.text = element_text(size = 10),                 # Legend text size
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5), 
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5), # Add main ticks
    axis.ticks.length = unit(0.25, "cm")
  )

ggsave(filename = file.path(out_put_path, "cricket_protein_pca_score_plot.png"), 
       plot = pca_plot, width = 10, height = 8, dpi = 300)



# data plotting -----------------------------------------------------------

custom_palette <- colorRampPalette(brewer.pal(12, "Paired"))(50)  # Expand Set3 to 20 colors

genus_data_48_top50 <- df_clean %>% 
  pivot_longer(2:1079)

# Step 1: Perform hierarchical clustering
# Pivot data to a wide format for clustering
#data_wide <- genus_data_48_top50 %>%
  #pivot_wider(names_from = name, values_from = value, values_fill = 0)

data_wide <- df_clean

# Use row-wise clustering (transpose the data for clustering by `Type`)
clustering <- hclust(dist(data_wide[-1]))  # Exclude the Type column from clustering
dendro_data <- as.dendrogram(clustering)

# Step 2: Reorder the x-axis based on the clustering
new_order <- data_wide$protein_name[clustering$order]
genus_data_48_top50 <- genus_data_48_top50 %>%
  mutate(protein_name = factor(protein_name, levels = new_order))

# Step 3: Create the dendrogram plot
dendro_plot <- ggdendrogram(dendro_data, rotate = FALSE) +
  theme_void() +
  theme(
    axis.text.x = element_blank(),  # Hide x-axis text
    axis.ticks.x = element_blank(), # Hide x-axis ticks
    axis.title.x = element_blank(), # Hide x-axis label
    plot.margin = margin(0, 0, 0, 0)  # Remove plot margins
  )+ labs(
    title = "Hierarchical clustering and heatmap for 48h gut microbiome"
  )

# Plotting the stacked bar chart using ggplot2 for the top 20 genera
micro_48h_plot <- ggplot(genus_data_48_top50, 
                         aes(x = protein_name, y = name, fill = value)) + # Mapping fill to 'name' to stack by genus
  geom_tile(aes(fill = value), color = "white") + 
  scale_fill_gradient(low = "white", high = "blue") +
  #geom_bar(stat = "identity", position = "stack", alpha = 0.8) +      # Create stacked bars
  #scale_fill_manual(values = scales::pal_hue()(20)) + # Optional: Custom color palette for 20 genera
  #scale_fill_manual(values = custom_palette) +
  theme_minimal() +                                      # Minimal theme for a cleaner look
  labs(
    x = "sample type",                                       # Label x-axis
    y = "protein",                                              # Label y-axis
    fill = "Relative abundance"                                            # Label the fill legend
  ) +
  theme(strip.background = element_rect(fill = "white", colour = "white"),
        plot.background = element_rect(fill = "white", colour = "white"),
        axis.title.x = element_text(size = 14, face = "bold"),  # X-axis title size
        axis.title.y = element_text(size = 14, face = "bold", angle = 90),  # Y-axis title size
        axis.text.x = element_text(size = 10, angle = 0, face = "bold"),                 # X-axis text size
        axis.text.y = element_text(size = 12, face = "bold"),                 # Y-axis text size
        legend.title = element_text(size = 12),                # Legend title size
        legend.text = element_text(size = 10),                 # Legend text size
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5), 
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", size = 0.5), # Add main ticks
        axis.ticks.length = unit(0.25, "cm"))


# Step 5: Combine the dendrogram and bar chart without gaps
final_48h_plot <- (dendro_plot + micro_48h_plot + plot_layout(heights = c(1, 4), guides = "collect")) & 
  theme(plot.background = element_rect(fill = "white", color = NA))  # Remove plot background


ggsave(filename = file.path(file_path, "microbiome_48h_heatmap_plot.png"), 
       plot = final_48h_plot, width = 12, height = 11, dpi = 300)




