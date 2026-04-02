## 14_tumor_growth_curves.R — Tumor growth curves from longitudinal measurements
if (!exists("sample_meta")) source("00_load_data.R")
library(patchwork)

tumor_data <- read_excel(
  file.path(repo, "data/metadata.xlsx"),
  sheet = "tumor_measurements"
) %>%
  filter(mouse_id %in% sample_meta$mouse_id) %>%
  left_join(sample_meta %>% select(mouse_id, lib_id, treatment), by = "mouse_id")

# Panel A: Individual growth curves colored by treatment
p_growth <- ggplot(tumor_data,
  aes(x = days_post_implant, y = tumor_volume_mm3,
      color = treatment, group = mouse_id)) +
  geom_line(linewidth = 0.8, alpha = 0.8) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_color_manual(values = c(vehicle = "#4daf4a", mirdametinib = "#e41a1c")) +
  labs(x = "Days post implant", y = expression("Tumor volume (mm"^3*")"),
       color = "Treatment",
       title = "WU-487 PDX tumor growth curves") +
  theme_white +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 11))

# Panel B: Terminal tumor volume bar plot with individual points
terminal <- tumor_data %>%
  group_by(mouse_id) %>%
  filter(days_post_implant == max(days_post_implant)) %>%
  ungroup()

p_terminal <- ggplot(terminal,
  aes(x = reorder(lib_id, tumor_volume_mm3), y = tumor_volume_mm3,
      fill = treatment)) +
  geom_col(alpha = 0.8, width = 0.7) +
  geom_text(aes(label = sprintf("%.0f", tumor_volume_mm3)),
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c(vehicle = "#4daf4a", mirdametinib = "#e41a1c")) +
  labs(x = NULL, y = expression("Terminal tumor volume (mm"^3*")"),
       fill = "Treatment",
       title = "Terminal tumor volumes at blood collection") +
  theme_white +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1))

combined <- p_growth + p_terminal +
  plot_annotation(
    title = "WU-487 MPNST PDX: Tumor Growth and Terminal Volumes",
    subtitle = "6 terminal bleeds (3 vehicle, 3 mirdametinib) | cfDNA collected at endpoint",
    theme = theme(plot.title = element_text(size = 15, face = "bold"),
                  plot.subtitle = element_text(size = 11))
  )

ggsave("plots/tumor_growth_curves.png", combined, width = 14, height = 6, dpi = 200)
cat("Saved plots/tumor_growth_curves.png\n")
