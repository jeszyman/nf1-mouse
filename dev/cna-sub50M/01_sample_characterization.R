## 01_sample_characterization.R — Sample-level covariates and bioanalyzer profiles
if (!exists("sample_meta")) source("00_load_data.R")
# library(bioanalyzeR)  # TODO: install in biotools or use ichor env

# --- Sample covariate summary table ---
sample_summary <- sample_meta %>%
  select(lib_id, mouse_id, treatment, cfdna_conc_pg_ul, tumor_vol_mm3)
write_tsv(sample_summary, "data/sample_summary.tsv")
cat("Sample summary:\n")
print(sample_summary)

# --- cfDNA concentration by treatment ---
p_cfdna <- ggplot(sample_meta, aes(x = lib_id, y = cfdna_conc_pg_ul, fill = treatment)) +
  geom_col(alpha = 0.8) +
  scale_fill_manual(values = c(vehicle = "#4daf4a", mirdametinib = "#e41a1c")) +
  labs(x = NULL, y = "cfDNA concentration (pg/µL)",
       title = "cfDNA input concentration by sample") +
  theme_white +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/cfdna_concentration.png", p_cfdna, width = 8, height = 5, dpi = 150)

# --- Tumor volume by treatment ---
p_tumor <- ggplot(sample_meta, aes(x = lib_id, y = tumor_vol_mm3, fill = treatment)) +
  geom_col(alpha = 0.8) +
  scale_fill_manual(values = c(vehicle = "#4daf4a", mirdametinib = "#e41a1c")) +
  labs(x = NULL, y = "Tumor volume (mm³)",
       title = "Terminal tumor volume by sample") +
  theme_white +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/tumor_volume.png", p_tumor, width = 8, height = 5, dpi = 150)

# --- cfDNA vs tumor volume scatter ---
p_scatter <- ggplot(sample_meta, aes(x = tumor_vol_mm3, y = cfdna_conc_pg_ul,
                                      color = treatment, label = lib_id)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, size = 3) +
  scale_color_manual(values = c(vehicle = "#4daf4a", mirdametinib = "#e41a1c")) +
  labs(x = "Tumor volume (mm³)", y = "cfDNA concentration (pg/µL)",
       title = "cfDNA concentration vs tumor volume") +
  theme_white
ggsave("plots/cfdna_vs_tumor.png", p_scatter, width = 7, height = 5, dpi = 150)

# --- Bioanalyzer profiles ---
# TODO: bioanalyzeR not in biotools env. Add or use separate env.
# WU-487 terminal bleeds collected 2024-08-20 (vehicle) and 2024-08-30 (mirdametinib)
# Bioanalyzer XML data at /mnt/gcs/jeszyman/projects/nf1-mouse/inputs/bioanalyzer/
cat("\nBioanalyzer analysis deferred — bioanalyzeR package not available in current env.\n")

cat("\nScript 01 complete. Plots in plots/\n")
