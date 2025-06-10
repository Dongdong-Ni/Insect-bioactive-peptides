library(tidyverse)
library(RColorBrewer)

file_path <- "D:/data_analysis/BP_insect_species"

df_kinetic <- read_csv(file.path(file_path, "digestive_kinetic_curve.CSV"))

df_1 <- df_kinetic %>% 
  select(1:4) %>% 
  pivot_longer(2:4)

df_2 <- df_kinetic %>% 
  select(1, 5:7) %>% 
  pivot_longer(2:4)

df_3 <- cbind(df_1, df_2[,3]) 

colnames(df_3) <- c('time', 'insect_type', 'digest_mean', 'digest_SD')

Digestive_kinetic_plot <- ggplot(
  df_3,
  mapping = aes(x =time, y = digest_mean, colour = insect_type, shape = factor(insect_type))
) + geom_point(size = 3) + geom_line()+ 
  geom_errorbar(aes(
    ymin = digest_mean - digest_SD, ymax = digest_mean + digest_SD
  ), width = 2) +
  labs(shape = "Insect type", colour = "Insect type", 
       title = "Digestive kinetic curve for insects",
       x = "Digestive time (minutes)",
       y = "Protein digestion percentage (%)")+
  theme_minimal()+
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

ggsave(filename = file.path(file_path, "Digestive_kinetic_plot.png"), 
       plot = Digestive_kinetic_plot, width = 12, height = 8, dpi = 300)


# particle size -----------------------------------------------------------

df_particle <- read_csv(file.path(file_path, "particle_size.CSV"))

df_particle[is.na(df_particle)] <- 0

df_partile_summarise <- df_particle %>% 
  group_by(insect_type, digestion_type, sample_type) %>% 
  summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop")

df_particle_long <- df_partile_summarise %>% 
  pivot_longer(4:104)
df_particle_long$name <- as.numeric(df_particle_long$name)

particle_size_plot <- ggplot(
  df_particle_long, 
  mapping = aes(x = name, y = value, 
                colour = digestion_type, shape = factor(insect_type))
) +geom_point(size = 3) +geom_line() +
  scale_x_log10()+
  theme_minimal() +
  labs(
    title = "Particle size distribution for insects",
    x = "Size (micrometer)",
    y = "Vol.Freq.(%)",
    shape = "Insect type", colour = "Digestive type"
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

ggsave(filename = file.path(file_path, "particle_size_plot.png"), 
       plot = particle_size_plot, width = 12, height = 8, dpi = 300)

# particle compare --------------------------------------------------------

df_par <- read_csv(file.path(file_path, "particle_size_1.CSV"))

df_par_summarise <- df_par %>% 
  group_by(insect_type, digestion_type, sample_type) %>% 
  summarise(
    across(everything(), list(mean = ~mean(.x, na.rm = TRUE), 
                              sd = ~sd(.x, na.rm = TRUE)), 
           .names = "{.col}_{.fn}"),
    .groups = "drop"
  )

colnames(df_par_summarise) <- c('insect_type', "digestion_type", "sample_type",
                                'D_4_3_mean', 'D_4_3_sd', 'uniformity_mean', 
                                'uniformity_sd', 'surface_area_mean', 'surface_area_sd',
                                'D_3_2_mean', 'D_3_2_sd', 'd_0.1_mean', 'd_0.1_sd',
                                'd_0.5_mean', 'd_0.5_sd', 'd_0.9_mean', 'd_0.9_sd')

D90_plot <- ggplot(
  df_par_summarise,
  mapping = aes(x = sample_type, y = d_0.9_mean, colour = insect_type, fill = digestion_type)
) + geom_col(position = position_dodge(0.8), alpha = 0.8)+
  geom_errorbar(aes(
    ymin = d_0.9_mean - d_0.9_sd ,ymax = d_0.9_mean + d_0.9_sd
  ), width = 0.3, position = position_dodge(0.8))+
  labs(fill = "Digestive type", colour = "Insect type", 
       title = "90% particles",
       x = "Sample types",
       y = "D90 particles") +
  theme_minimal()+
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


D43_plot <- ggplot(
  df_par_summarise,
  mapping = aes(x = sample_type, y = D_4_3_mean, colour = insect_type, fill = digestion_type)
) + geom_col(position = position_dodge(0.8), alpha = 0.8)+
  geom_errorbar(aes(
    ymin = D_4_3_mean - D_4_3_sd ,ymax = D_4_3_mean + D_4_3_sd
  ), width = 0.3, position = position_dodge(0.8))+
  labs(fill = "Digestive type", colour = "Insect type", 
       title = "D[4,3] particles",
       x = "Sample types",
       y = "D[4,3] particles") +
  theme_minimal()+
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

ggsave(filename = file.path(file_path, "D90_plot.png"), 
       plot = D90_plot, width = 12, height = 10, dpi = 300)
ggsave(filename = file.path(file_path, "D43_plot.png"), 
       plot = D43_plot, width = 12, height = 10, dpi = 300)

ggplot(
  df_par_summarise,
  mapping = aes(x = sample_type, y = D_3_2_mean, colour = insect_type, fill = digestion_type)
) + geom_col(position = position_dodge(0.8), alpha = 0.8)+
  geom_errorbar(aes(
    ymin = D_3_2_mean - D_3_2_sd ,ymax = D_3_2_mean + D_3_2_sd
  ), width = 0.3, position = position_dodge(0.8))+
  labs(fill = "Digestive type", colour = "Insect type", 
       title = "D[3,2] particles",
       x = "Sample types",
       y = "D[3,2] particles") +
  theme_minimal()+
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

