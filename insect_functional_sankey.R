library(tidyverse)
library(networkD3)
library(ggalluvial)
library(RColorBrewer)
library(dplyr)
library(webshot)
library(htmlwidgets)

file_path <- "D:/data_analysis/BP_insect_species"

df_all <- read_csv(file.path(file_path, "BPs_proteininfo.csv"))

df_bsf_functional_protein <- read_csv("D:/data_analysis/BP_insect_species/BSF/BSF_functional_protein_summary_quant.csv")

df_cricket_functional_protein <- read_csv("D:/data_analysis/BP_insect_species/cricket/cricket_functional_protein_summary_quant.csv")

df_mealworm_functional_protein <- read_csv("D:/data_analysis/BP_insect_species/mealworm/mealworm_functional_protein_summary_quant.csv")

df_functional_protein <- rbind(df_bsf_functional_protein, df_cricket_functional_protein)

df_functional_protein <- rbind(df_functional_protein, df_mealworm_functional_protein)

accession_ids <- df_functional_protein[[1]]

df_filter_p <- data.frame()

for (a in accession_ids) {
  
  df_filtered <- df_all %>% 
    filter(str_detect(Accessions, a))
  
  df_quant <- df_functional_protein %>% 
    filter(functional_protein == a)
  
  df_filtered <- cbind(df_filtered, df_quant)
  
  df_filter_p <- rbind(df_filter_p, df_filtered)
  
}


# sankey ------------------------------------------------------------------


df_sankey <- df_filter_p %>% 
  select(Functional, P_sequence, protein.ID, Names, functional_protein, total) %>% 
  filter(total >= 1)


nodes <- data.frame(
  name = unique(c(
    df_sankey$protein.ID,
    df_sankey$functional_protein,
    df_sankey$P_sequence,
    df_sankey$Functional
  ))
)

# Create a dataframe for links (relationships between nodes)
links <- data.frame(
  source = match(df_sankey$protein.ID, nodes$name) - 1,        # Match food_name and get index
  target = match(df_sankey$functional_protein, nodes$name) - 1,       # Match protein_ID and get index
  value = df_sankey$total                        # Use the flow quantity
) %>%
  bind_rows(data.frame(                                       # Add relationships for protein_ID → P_sequence
    source = match(df_sankey$functional_protein, nodes$name) - 1,
    target = match(df_sankey$P_sequence, nodes$name) - 1,
    value = df_sankey$total
  )) %>%
  bind_rows(data.frame(                                       # Add relationships for P_sequence → Functional
    source = match(df_sankey$P_sequence, nodes$name) - 1,
    target = match(df_sankey$Functional, nodes$name) - 1,
    value = df_sankey$total
  ))

# Add a group column to nodes for coloring (example: manually define based on input levels)
nodes$group <- case_when(
  nodes$name %in% unique(df_sankey$protein.ID) ~ "Insect Name",
  nodes$name %in% unique(df_sankey$functional_protein) ~ "Biological functional proteins",
  nodes$name %in% unique(df_sankey$P_sequence) ~ "Bioactive peptides",
  nodes$name %in% unique(df_sankey$Functional) ~ "Functional",
  TRUE ~ "Other"
)


colourScale <- 'd3.scaleOrdinal()
                .domain(["Insect Name", "Biological functional proteins", "Bioactive peptides", "Functional"])
                .range(["#FF5733", "#33FF57", "#3357FF", "#FFD700"])'


# Generate Sankey Diagram
sankey <- sankeyNetwork(
  Links = links,
  Nodes = nodes,
  Source = "source",
  Target = "target",
  Value = "value",
  NodeID = "name",
  fontSize = 12,
  nodeWidth = 15,
  nodePadding = 15,
  colourScale = colourScale,
  NodeGroup = "group"
)

# Print the diagram
sankey

# Save the Sankey diagram as an HTML file temporarily
saveWidget(sankey, "sankey_functional_proteins_1.html")

# (Optional) Save as JPEG
webshot("sankey_functional_proteins_1.html", file = "sankey_functional_proteins_1.jpeg")

webshot("sankey_functional_proteins_1.html", file = "sankey_functional_proteins_1.pdf")


# heatmap -----------------------------------------------------------------

df_heatmap <- df_filter_p %>% 
  select(Functional, protein.ID, total) %>%
  group_by(protein.ID, Functional) %>% 
  summarise(numb = mean(total))

custom_colors <- c("white", "blue", "green","yellow", "red") 

heatmap_plot <- ggplot(df_heatmap, 
                   aes(x = protein.ID, y = Functional, fill = numb)) + # Mapping fill to 'name' to stack by genus
  #geom_tile(aes(fill = pp_number), color = "white") + 
  #scale_fill_gradient(low = "white", high = "blue") +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = custom_colors) +
  theme_minimal() +                                      # Minimal theme for a cleaner look
  labs(
    x = "Samples",                                       # Label x-axis
    y = "Functions",                                              # Label y-axis
    fill = "Protein quantity (%)"                                            # Label the fill legend
  ) + 
  theme(strip.background = element_rect(fill = "white", colour = "white"),
        plot.background = element_rect(fill = "white", colour = "white"),
        axis.title.x = element_text(size = 14, face = "bold"),  # X-axis title size
        axis.title.y = element_text(size = 14, face = "bold", angle = 90),  # Y-axis title size
        axis.text.x = element_text(size = 10, angle = 0, hjust = 1, face = "bold"),                 # X-axis text size
        axis.text.y = element_text(size = 12, face = "bold"),                 # Y-axis text size
        legend.title = element_text(size = 12),                # Legend title size
        legend.text = element_text(size = 10),                 # Legend text size
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5), 
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", size = 0.5), # Add main ticks
        axis.ticks.length = unit(0.25, "cm"))

ggsave(filename = file.path(file_path, "functional_proteins_plot.png"), 
       plot = heatmap_plot, width = 9, height = 8, dpi = 300)

