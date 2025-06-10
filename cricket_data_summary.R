library(tidyverse)

file_path <- "D:/data_analysis/BP_insect_species/cricket"


df_raw <- read_csv(file.path(file_path, "cricket_digested_peptide_output.csv"))
colnames(df_raw) <- make.names(colnames(df_raw))

df_Q <- df_raw %>% 
  filter(Identity.Percentage == 100) %>% 
  filter(grepl("^BP", Subject.Description)) %>% 
  mutate(Description = Subject.Description) %>% 
  separate(Description, into = c('BP', 'Number', 'Method', 'Functional', 
                                 'mol_w1', 'mol_w2', 'se'), sep = "\\|")

# get digested peptide sequence -------------------------------------------



bioactive_p_database <- file.path(file_path, "bioactive_peptides_database.fasta")

peptide_file <- file.path(file_path,  "cricket_digested_peptide.fasta")
lines <- readLines(peptide_file)
# Find headers and indices
header_indices <- grep("^>", lines)

# Append an end line index for easy sequence parsing
header_indices <- c(header_indices, length(lines) + 1)

# Initialize lists to store headers and sequences
headers <- character(length(header_indices) - 1)
sequences <- character(length(header_indices) - 1)

# Parse the headers and sequences
for (i in seq_along(header_indices[-length(header_indices)])) {
  headers[i] <- gsub("^>", "", lines[header_indices[i]])
  sequences[i] <- paste0(
    lines[(header_indices[i] + 1):(header_indices[i + 1] - 1)],
    collapse = ""
  )
}

# Create a data frame of headers and sequences
digested_peptides <- data.frame(
  Header = headers,
  Sequence = sequences,
  stringsAsFactors = FALSE
)


lines <- readLines(bioactive_p_database)

headers <- lines[grep("^>", lines)]

sequences <- lines[!grepl("^>", lines)]

bioactive_p_db <- data.frame(
  Header = gsub("^>", "", headers),
  Sequence = sequences
)


# -------------------------------------------------------------------------

colnames(digested_peptides) <- c("Query.ID", "P_sequence")

colnames(bioactive_p_db) <- c("Subject.Description", "db_sequence")

df_sele <- left_join(df_Q, digested_peptides, by = "Query.ID")

df_sele <- left_join(df_sele, bioactive_p_db, by = "Subject.Description")

df_sele <- df_sele %>% 
  mutate(db_sequence_n = str_count(db_sequence, "[a-zA-Z]"))
df_sele <- df_sele %>% 
  mutate(hit_percent = round((Alignment.Length/ db_sequence_n)*100, 2))

df_100 <- df_sele %>% 
  filter(hit_percent >= 80) %>% 
  mutate(protein.ID = Query.ID) %>% 
  separate(protein.ID, into = c('protein.ID', 'gp','gastric.ID', 'sp','small.ID'), 
           sep = '_')

df_selected <- df_100 %>% 
  select(
    "Query.ID", "Alignment.Length", "Hit.Sequence", "Method" ,
    "Functional", "P_sequence" , "db_sequence" , "db_sequence_n" ,
    "hit_percent", "protein.ID", "gp", "gastric.ID", "sp"
  ) 


BP_summary_time <- df_selected %>% 
  group_by(gastric.ID, Functional) %>% 
  summarise(functional_bp = n())

BP_summary_time$gastric.ID <- as.character(BP_summary_time$gastric.ID)

df_summarise$sample_type <- factor(
  df_summarise$sample_type,
  levels = c("0_defatted", "0_whole", "60_defatted",  "60_whole", 
             "120_defatted", "120_whole", "180_defatted", "180_whole", 
             "300_defatted", "300_whole")
)

BP_summary_time$gastric.ID <- factor(
  BP_summary_time$gastric.ID,
  levels = c("0", "60", "120", "180")
)


BP_total_number_summary <- df_selected %>% 
  group_by(gastric.ID) %>% 
  summarise(functional_pp = n_distinct(Hit.Sequence))

BP_total_num_rep_summary <- df_selected %>% 
  group_by(gastric.ID) %>% 
  summarise(functional_pp = n())

df_selected %>% 
  write_csv(file.path(file_path, "cricket_BP.csv"))



# plot --------------------------------------------------------------------

BP_food_plot <- ggplot(
  BP_summary_time, 
  mapping = aes(x = gastric.ID, y = functional_bp, fill = Functional)
) + 
  geom_bar(stat = "identity", position = "stack") +
  geom_text(
    aes(label = functional_bp),
    position = position_stack(vjust = 0.5),
    color = "white",
    size = 3
  ) +
  theme_minimal() +
  labs(
    title = "Bioactive peptides from cricket digestion",
    x = "Digestive time (minutes)",
    y = "Number of Bioactive Peptides",
    fill = "Function"
  ) +
  theme(
    strip.background = element_rect(fill = "white", colour = "white"),
    plot.background = element_rect(fill = "white", colour = "white"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold", angle = 90),
    axis.text.x = element_text(size = 14, angle = 0, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.25, "cm")
  )

BP_food_plot
ggsave(filename = file.path(file_path, "cricket_BP_summary_plot.png"), 
       plot = BP_food_plot, width = 10, height = 8, dpi = 300)



# protein tracing ---------------------------------------------------------

cricket_digested_pp <- read_csv(file.path(file_path, "cricket_digested_peptides_fa.csv"))

colnames(cricket_digested_pp) <- c("Sequence", "Accessions","Names", "id" ,"process", 
                                   "time" ,"rept" ,"x", "Query.ID")

df_selected <- left_join(df_selected, cricket_digested_pp, by = "Query.ID")


df_selected%>% 
  distinct(Accessions)

accessions_list <- df_selected %>% 
  distinct(Accessions) %>% 
  mutate(Accessions = str_split(Accessions, ";\\s*")) %>% 
  unnest(Accessions) %>% 
  distinct()

accessions_list %>% 
  write_csv(file.path(file_path, "cricket_functional_proteins.csv"))

