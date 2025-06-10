library(tidyverse)

file_path <- "D:/data_analysis/BP_insect_species"

df_mealworm <- read.csv("D:/data_analysis/BP_insect_species/mealworm/mealwom_BP.csv") 

df_cricket <- read.csv("D:/data_analysis/BP_insect_species/cricket/cricket_BP.csv") 

df_BSF <- read.csv("D:/data_analysis/BP_insect_species/BSF/BSF_BP.csv") 

df_cricket <- df_cricket[,1:13]

df_all <- rbind(df_mealworm, df_cricket)

df_all <- rbind(df_all, df_BSF)

# total number ------------------------------------------------------------


df_summary <- df_all %>% 
  group_by(protein.ID, gastric.ID) %>% 
  summarise(pp_number = n_distinct(P_sequence)) %>% 
  mutate(Sample = paste0(protein.ID, "_", gastric.ID))

df_summary$gastric.ID <- as.character(df_summary$gastric.ID)

df_summary$Sample <- factor(
  df_summary$Sample, 
  levels = c("BSF_0","BSF_60","BSF_120","BSF_180","BSF_300",
             "cricket_0","cricket_60","cricket_120","cricket_180", 
             "mealworm_0","mealworm_60", "mealworm_120","mealworm_180", "mealworm_300")
)

BPs_plot <- ggplot(
  df_summary, 
  mapping = aes(x = Sample, y = pp_number, fill = gastric.ID)
) + geom_col(alpha = 0.8) +
  geom_text(
    aes(label = pp_number),
    position = position_stack(vjust = 0.5),
    color = "black",
    size = 4
  ) +
  theme_minimal() +
  labs(
    title = "Released peptides during digestion",
    x = "Samples",
    y = "Number of Peptides",
    fill = "Digestive time 
    (minutes)"
  ) +
  theme(
    strip.background = element_rect(fill = "white", colour = "white"),
    plot.background = element_rect(fill = "white", colour = "white"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold", angle = 90),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.25, "cm")
  )
BPs_plot
ggsave(filename = file.path(file_path, "BPs_total_plot.png"), 
       plot = BPs_plot, width = 12, height = 8, dpi = 300)


# in functions ------------------------------------------------------------

df_func_summary <- df_all %>% 
  group_by(protein.ID, gastric.ID, Functional) %>% 
  summarise(pp_number = n_distinct(P_sequence)) %>% 
  mutate(Sample = paste0(protein.ID, "_", gastric.ID))

custom_colors <- c("white", "blue", "green","yellow", "red") 

df_func_summary$Sample <- factor(
  df_func_summary$Sample, 
  levels = c("BSF_0","BSF_60","BSF_120","BSF_180","BSF_300",
             "cricket_0","cricket_60","cricket_120","cricket_180", 
             "mealworm_0","mealworm_60", "mealworm_120","mealworm_180", "mealworm_300")
)

bps_plot <- ggplot(df_func_summary, 
                        aes(x = Sample, y = Functional, fill = pp_number)) + # Mapping fill to 'name' to stack by genus
  #geom_tile(aes(fill = pp_number), color = "white") + 
  #scale_fill_gradient(low = "white", high = "blue") +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = custom_colors) +
  theme_minimal() +                                      # Minimal theme for a cleaner look
  labs(
    x = "Samples",                                       # Label x-axis
    y = "Functions",                                              # Label y-axis
    fill = "Number of peptides"                                            # Label the fill legend
  ) +
  theme(strip.background = element_rect(fill = "white", colour = "white"),
        plot.background = element_rect(fill = "white", colour = "white"),
        axis.title.x = element_text(size = 14, face = "bold"),  # X-axis title size
        axis.title.y = element_text(size = 14, face = "bold", angle = 90),  # Y-axis title size
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1, face = "bold"),                 # X-axis text size
        axis.text.y = element_text(size = 12, face = "bold"),                 # Y-axis text size
        legend.title = element_text(size = 12),                # Legend title size
        legend.text = element_text(size = 10),                 # Legend text size
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5), 
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", size = 0.5), # Add main ticks
        axis.ticks.length = unit(0.25, "cm"))

bps_plot

ggsave(filename = file.path(file_path, "BPs_func_plot.png"), 
       plot = bps_plot, width = 11, height = 8, dpi = 300)


# BPs data reform ---------------------------------------------------------

df_mealworm_p <- read.csv("D:/data_analysis/BP_insect_species/mealworm/mealworm_digested_peptides_fa.csv") %>% 
  select(Sequence, Accessions, Names, time, Header)

df_cricket_p <- read.csv("D:/data_analysis/BP_insect_species/cricket/cricket_digested_peptides_fa.csv") %>% 
  select(Sequence, Accessions, Names, time, Header)

df_BSF_p <- read.csv("D:/data_analysis/BP_insect_species/BSF/BSF_digested_peptides_fa.csv") %>% 
  select(Sequence, Accessions, Names, time, Header)

colnames(df_mealworm_p) <- c('Sequence', 'Accessions', 'Names', 'time', 'Query.ID')
colnames(df_cricket_p) <- c('Sequence', 'Accessions', 'Names', 'time', 'Query.ID')
colnames(df_BSF_p) <- c('Sequence', 'Accessions', 'Names', 'time', 'Query.ID')

df_all_p <- rbind(df_mealworm_p, df_cricket_p)

df_all_p <- rbind(df_all_p, df_BSF_p)

df_all <- left_join(df_all, df_all_p, by = "Query.ID")


df_all %>% 
  write_csv(file.path(file_path, "BPs_proteininfo.csv"))

df_all_p %>% 
  write_csv(file.path(file_path, "proteininfo.csv"))
# sankey plot -------------------------------------------------------------
library(networkD3)
library(ggalluvial)
library(RColorBrewer)
library(dplyr)

df_sankey <- df_all %>% 
  select(Functional, P_sequence, protein.ID, Accessions, Names) %>% 
  mutate(quantity = "1")

df_sankey <- df_sankey %>% 
  mutate(Accessions = gsub(";.*", "", Accessions))


nodes <- data.frame(
  name = unique(c(
    df_sankey$protein.ID,
    df_sankey$Accessions,
    df_sankey$P_sequence,
    df_sankey$Functional
  ))
)

# Create a dataframe for links (relationships between nodes)
links <- data.frame(
  source = match(df_sankey$protein.ID, nodes$name) - 1,        # Match food_name and get index
  target = match(df_sankey$Accessions, nodes$name) - 1,       # Match protein_ID and get index
  value = df_sankey$quantity                        # Use the flow quantity
) %>%
  bind_rows(data.frame(                                       # Add relationships for protein_ID → P_sequence
    source = match(df_sankey$Accessions, nodes$name) - 1,
    target = match(df_sankey$P_sequence, nodes$name) - 1,
    value = df_sankey$quantity
  )) %>%
  bind_rows(data.frame(                                       # Add relationships for P_sequence → Functional
    source = match(df_sankey$P_sequence, nodes$name) - 1,
    target = match(df_sankey$Functional, nodes$name) - 1,
    value = df_sankey$quantity
  ))

# Define Custom Colors Using JS ColourScale

# Add a group column to nodes for coloring (example: manually define based on input levels)
nodes$group <- case_when(
  nodes$name %in% unique(df_sankey$protein.ID) ~ "Insect Name",
  nodes$name %in% unique(df_sankey$Accessions) ~ "Proteins",
  nodes$name %in% unique(df_sankey$P_sequence) ~ "Bioactive peptides",
  nodes$name %in% unique(df_sankey$Functional) ~ "Functional",
  TRUE ~ "Other"
)


colourScale <- 'd3.scaleOrdinal()
                .domain(["Insect Name", "Proteins", "Bioactive peptides", "Functional"])
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

library(webshot)
library(htmlwidgets)

# Save the Sankey diagram as an HTML file temporarily
saveWidget(sankey, "sankey_all1_candidate.html")

# (Optional) Save as JPEG
webshot("sankey_temp_candidate.html", file = "sankey_diagram.jpeg")

webshot("sankey_temp_candidate.html", file = "sankey_temp_candidate.pdf")


# sankey time separation --------------------------------------------------

df_sankey_0 <- df_all %>% 
  select(Functional, P_sequence, protein.ID, Accessions, Names, time) %>% 
  filter(time == "0") %>% 
  mutate(quantity = "5")

df_sankey_0 <- df_sankey_0 %>% 
  mutate(Accessions = gsub(";.*", "", Accessions))


nodes <- data.frame(
  name = unique(c(
    df_sankey_0$protein.ID,
    df_sankey_0$Accessions,
    df_sankey_0$P_sequence,
    df_sankey_0$Functional
  ))
)

# Create a dataframe for links (relationships between nodes)
links <- data.frame(
  source = match(df_sankey_0$protein.ID, nodes$name) - 1,        # Match food_name and get index
  target = match(df_sankey_0$Accessions, nodes$name) - 1,       # Match protein_ID and get index
  value = df_sankey_0$quantity                        # Use the flow quantity
) %>%
  bind_rows(data.frame(                                       # Add relationships for protein_ID → P_sequence
    source = match(df_sankey_0$Accessions, nodes$name) - 1,
    target = match(df_sankey_0$P_sequence, nodes$name) - 1,
    value = df_sankey_0$quantity
  )) %>%
  bind_rows(data.frame(                                       # Add relationships for P_sequence → Functional
    source = match(df_sankey_0$P_sequence, nodes$name) - 1,
    target = match(df_sankey_0$Functional, nodes$name) - 1,
    value = df_sankey_0$quantity
  ))

# Define Custom Colors Using JS ColourScale

# Add a group column to nodes for coloring (example: manually define based on input levels)
nodes$group <- case_when(
  nodes$name %in% unique(df_sankey_0$protein.ID) ~ "Insect Name",
  nodes$name %in% unique(df_sankey_0$Accessions) ~ "Proteins",
  nodes$name %in% unique(df_sankey_0$P_sequence) ~ "Bioactive peptides",
  nodes$name %in% unique(df_sankey_0$Functional) ~ "Functional",
  TRUE ~ "Other"
)


colourScale <- 'd3.scaleOrdinal()
                .domain(["Insect Name", "Proteins", "Bioactive peptides", "Functional", "Other"])
                .range(["#FF5733", "#33FF57", "#3357FF", "#FFD700", "#CCCCCC"])'


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
  nodePadding = 30,
  colourScale = colourScale,
  NodeGroup = "group"
)

# Print the diagram
sankey

# Save the Sankey diagram as an HTML file temporarily
saveWidget(sankey, "sankey__01_candidate.html")


# -------------------------------------------------------------------------

df_sankey_g <- df_all %>% 
  select(Functional, P_sequence, protein.ID, Accessions, Names, time) %>% 
  filter(time %in% c("60", "120")) %>% 
  mutate(quantity = "1")

df_sankey_g <- df_sankey_g %>% 
  mutate(Accessions = gsub(";.*", "", Accessions))


nodes <- data.frame(
  name = unique(c(
    df_sankey_g$protein.ID,
    df_sankey_g$Accessions,
    df_sankey_g$P_sequence,
    df_sankey_g$Functional
  ))
)

# Create a dataframe for links (relationships between nodes)
links <- data.frame(
  source = match(df_sankey_g$protein.ID, nodes$name) - 1,        # Match food_name and get index
  target = match(df_sankey_g$Accessions, nodes$name) - 1,       # Match protein_ID and get index
  value = df_sankey_g$quantity                        # Use the flow quantity
) %>%
  bind_rows(data.frame(                                       # Add relationships for protein_ID → P_sequence
    source = match(df_sankey_g$Accessions, nodes$name) - 1,
    target = match(df_sankey_g$P_sequence, nodes$name) - 1,
    value = df_sankey_g$quantity
  )) %>%
  bind_rows(data.frame(                                       # Add relationships for P_sequence → Functional
    source = match(df_sankey_g$P_sequence, nodes$name) - 1,
    target = match(df_sankey_g$Functional, nodes$name) - 1,
    value = df_sankey_g$quantity
  ))

# Define Custom Colors Using JS ColourScale

# Add a group column to nodes for coloring (example: manually define based on input levels)
nodes$group <- case_when(
  nodes$name %in% unique(df_sankey_g$protein.ID) ~ "Insect Name",
  nodes$name %in% unique(df_sankey_g$Accessions) ~ "Proteins",
  nodes$name %in% unique(df_sankey_g$P_sequence) ~ "Bioactive peptides",
  nodes$name %in% unique(df_sankey_g$Functional) ~ "Functional",
  TRUE ~ "Other"
)


colourScale <- 'd3.scaleOrdinal()
                .domain(["Insect Name", "Proteins", "Bioactive peptides", "Functional", "Other"])
                .range(["#FF5733", "#33FF57", "#3357FF", "#FFD700", "#CCCCCC"])'


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
  nodePadding = 10,
  colourScale = colourScale,
  NodeGroup = "group"
)

# Print the diagram
sankey

# Save the Sankey diagram as an HTML file temporarily
saveWidget(sankey, "sankey__g_candidate.html")


# -------------------------------------------------------------------------

df_sankey_s <- df_all %>% 
  select(Functional, P_sequence, protein.ID, Accessions, Names, time) %>% 
  filter(time %in% c("180", "300")) %>% 
  mutate(quantity = "1")

df_sankey_s <- df_sankey_s %>% 
  mutate(Accessions = gsub(";.*", "", Accessions))


nodes <- data.frame(
  name = unique(c(
    df_sankey_s$protein.ID,
    df_sankey_s$Accessions,
    df_sankey_s$P_sequence,
    df_sankey_s$Functional
  ))
)

# Create a dataframe for links (relationships between nodes)
links <- data.frame(
  source = match(df_sankey_s$protein.ID, nodes$name) - 1,        # Match food_name and get index
  target = match(df_sankey_s$Accessions, nodes$name) - 1,       # Match protein_ID and get index
  value = df_sankey_s$quantity                        # Use the flow quantity
) %>%
  bind_rows(data.frame(                                       # Add relationships for protein_ID → P_sequence
    source = match(df_sankey_s$Accessions, nodes$name) - 1,
    target = match(df_sankey_s$P_sequence, nodes$name) - 1,
    value = df_sankey_s$quantity
  )) %>%
  bind_rows(data.frame(                                       # Add relationships for P_sequence → Functional
    source = match(df_sankey_s$P_sequence, nodes$name) - 1,
    target = match(df_sankey_s$Functional, nodes$name) - 1,
    value = df_sankey_s$quantity
  ))

# Define Custom Colors Using JS ColourScale

# Add a group column to nodes for coloring (example: manually define based on input levels)
nodes$group <- case_when(
  nodes$name %in% unique(df_sankey_s$protein.ID) ~ "Insect Name",
  nodes$name %in% unique(df_sankey_s$Accessions) ~ "Proteins",
  nodes$name %in% unique(df_sankey_s$P_sequence) ~ "Bioactive peptides",
  nodes$name %in% unique(df_sankey_s$Functional) ~ "Functional",
  TRUE ~ "Other"
)


colourScale <- 'd3.scaleOrdinal()
                .domain(["Insect Name", "Proteins", "Bioactive peptides", "Functional", "Other"])
                .range(["#FF5733", "#33FF57", "#8A2BE2", "#FFD700","#CCCCCC" ])'


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
  nodePadding = 10,
  colourScale = colourScale,
  NodeGroup = "group"
)

# Print the diagram
sankey

# Save the Sankey diagram as an HTML file temporarily
saveWidget(sankey, "sankey__s_candidate.html")

webshot("sankey__s_candidate.html", file = "sankey_s_diagram.jpeg")

webshot("sankey__s_candidate.html", file = "sankey_s_candidate.pdf")


# BPs ggven ---------------------------------------------------------------
library(VennDiagram)
library(grid)

df_all <- read_csv(file.path(file_path, "BPs_proteininfo.csv"))

df_ven <- df_all %>% 
  select(protein.ID, Sequence, time, Functional)

df_time_label <- data.frame(
  time = c(0, 60, 120, 180, 300),
  digest = c("0", "gastic", "gastric", "intestine", "intestine")
)

df_ven <- left_join(df_ven, df_time_label, by = "time")

insect_type <- c("mealworm", "cricket", "BSF")

for (a in insect_type) {
  # Filter the data for the specific insect type
  df_ven_plot <- df_ven %>% 
    filter(protein.ID == a)
  
  # Prepare the data for the Venn diagram
  venn_input <- list(
    Set1 = df_ven_plot$Sequence[df_ven_plot$digest == "0"],
    Set2 = df_ven_plot$Sequence[df_ven_plot$digest == "gastric"],
    Set3 = df_ven_plot$Sequence[df_ven_plot$digest == "intestine"]
  )
  
  # Define the colors for the Venn diagram
  color <- c("#00FF00", "#FF0000", "#0000FF")
  
  # Create the Venn diagram
  venn_plot <- venn.diagram(
    x = venn_input,
    category.names = c("non-digest", "gastric", "intestine"),
    filename = NULL,
    fill = color,       # Set the fill colors
    alpha = 0.5,        # Transparency level of the colors
    fontface = "bold",  # Bold font for labels
    fontfamily = "sans",# Font family
    cex = 2.2,          # Font size for counts inside the circles
    cat.fontfamily = "sans",   # Font family for category labels
    cat.cex = 1.8,        # Font size for category labels
    cat.dist = c(0.1, 0.1, 0.1),  # Distance of category labels from the plot
    cat.pos = c(0, 0, 0),  # Position of category labels
    cat.default.pos = "text",  # Set default position of category labels to outer
    cat.col = c("white", "white", "white"),  # Set label text colors
    margin = 0.05,       # Set the margin around the plot
    main = a,  # Title of the diagram
    main.cex = 2,        # Font size for the main title
    main.fontface = "bold",   # Bold font for the main title
    main.fontfamily = "sans"  # Font family for the main title
  )
  
  
  grid.draw(venn_plot)
  
  ggsave(filename = file.path(file_path, paste0(a, "_BPs_venn_plot_001.png")), 
         plot = venn_plot, width = 20, height = 18,
         dpi = 300, units = "cm")
  
}

digestive_time <- c('0', 'gastric', 'intestine')

for (b in digestive_time) {
  # Filter the data for the specific insect type
  df_ven_plot <- df_ven %>% 
    filter(digest == b)
  
  # Prepare the data for the Venn diagram
  venn_input <- list(
    Set1 = df_ven_plot$Sequence[df_ven_plot$protein.ID == "mealworm"],
    Set2 = df_ven_plot$Sequence[df_ven_plot$protein.ID == "cricket"],
    Set3 = df_ven_plot$Sequence[df_ven_plot$protein.ID == "BSF"]
  )
  
  # Define the colors for the Venn diagram
  color <- c("#00FF00", "#FF0000", "#0000FF")
  
  # Create the Venn diagram
  venn_plot <- venn.diagram(
    x = venn_input,
    category.names = c("mealworm", "cricket", "BSF"),
    filename = NULL,
    fill = color,       # Set the fill colors
    alpha = 0.5,        # Transparency level of the colors
    fontface = "bold",  # Bold font for labels
    fontfamily = "sans",# Font family
    cex = 2.2,          # Font size for counts inside the circles
    cat.fontfamily = "sans",   # Font family for category labels
    cat.cex = 1.8,        # Font size for category labels
    cat.dist = c(0.1, 0.1, 0.1),  # Distance of category labels from the plot
    cat.pos = c(0, 0, 0),  # Position of category labels
    cat.default.pos = "text",  # Set default position of category labels to outer
    cat.col = c("white", "white", "white"),  # Set label text colors
    margin = 0.05,       # Set the margin around the plot
    main = b,  # Title of the diagram
    main.cex = 2,        # Font size for the main title
    main.fontface = "bold",   # Bold font for the main title
    main.fontfamily = "sans"  # Font family for the main title
  )
  
  
  grid.draw(venn_plot)
  
  ggsave(filename = file.path(file_path, paste0(b, "_BPs_venn_plot_001.png")), 
         plot = venn_plot, width = 20, height = 18,
         dpi = 300, units = "cm")
  
}


venn_input <- list(
  Set1 = df_ven$Sequence[df_ven$protein.ID == "mealworm"],
  Set2 = df_ven$Sequence[df_ven$protein.ID == "cricket"],
  Set3 = df_ven$Sequence[df_ven$protein.ID == "BSF"]
)

# Define the colors for the Venn diagram
color <- c("#00FF00", "#FF0000", "#0000FF")

# Create the Venn diagram
venn_plot <- venn.diagram(
  x = venn_input,
  category.names = c("mealworm", "cricket", "BSF"),
  filename = NULL,
  fill = color,       # Set the fill colors
  alpha = 0.5,        # Transparency level of the colors
  fontface = "bold",  # Bold font for labels
  fontfamily = "sans",# Font family
  cex = 2.2,          # Font size for counts inside the circles
  cat.fontfamily = "sans",   # Font family for category labels
  cat.cex = 1.8,        # Font size for category labels
  cat.dist = c(0.1, 0.1, 0.1),  # Distance of category labels from the plot
  cat.pos = c(0, 0, 0),  # Position of category labels
  cat.default.pos = "text",  # Set default position of category labels to outer
  cat.col = c("white", "white", "white"),  # Set label text colors
  margin = 0.05,       # Set the margin around the plot
  main = "Bioavtive peptides for insects",  # Title of the diagram
  main.cex = 2,        # Font size for the main title
  main.fontface = "bold",   # Bold font for the main title
  main.fontfamily = "sans"  # Font family for the main title
)


grid.draw(venn_plot)

ggsave(filename = file.path(file_path, paste0("insect_BPs_venn_plot_001.png")), 
       plot = venn_plot, width = 20, height = 18,
       dpi = 300, units = "cm")



