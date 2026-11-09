# The data used below are private. DO NOT COMMIT generated CSVs.
# FIGURE: DOCUMENTED ACCESS AND CONSUMPTION BY FACULTATIVE SCAVENGERS
# Designed for the muladares / feeding points manuscript.
# Produces figures in RStudio from BASE_DATOS_BUITRES.csv.
# Observed proportions ONLY: not causal inference, and bars do not adjust for
# repeated provisioning events nested within the same hunting estate.

needed <- c("readr", "dplyr", "tidyr", "ggplot2", "scales")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install missing packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2); library(scales)
})

carpeta <- Sys.getenv("SFS_DATA_DIR", unset=".")
file <- file.path(carpeta, "BASE_DATOS_BUITRES.csv")
stopifnot(file.exists(file))
out <- file.path(carpeta, "FIGURAS_MANUSCRITO")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

raw <- readr::read_csv2(file, locale = readr::locale(encoding = "Latin1"),
                        show_col_types = FALSE)

# One record per PROVISIONING EVENT; not one record per camera detection.
base <- raw %>%
  filter(POST_aporte == "SI", !is.na(ID_aporte), ID_aporte != "") %>%
  mutate(opportunist = Clasificación == "Oportunista",
         opportunist = coalesce(opportunist, FALSE),
         # The field notes distinguish consumption of carrion from
         # dog grazing: only carrion consumption counts.
         grass_only = grepl("hierba|grass",
                            coalesce(as.character(Notas), ""),
                            ignore.case=TRUE),
         documented_consumption = opportunist & (behav == "Consumo") & !grass_only,
         documented_consumption = coalesce(documented_consumption, FALSE)) %>%
  group_by(ID_aporte, ID_COTO, Tipo_recinto, Valla_acondicionada) %>%
  summarise(
    opp_recorded = any(opportunist),
    opp_consumed = any(documented_consumption),
    .groups = "drop"
  )

if (anyDuplicated(base$ID_aporte)) stop("An ID_aporte appears under multiple site combinations. Inspect the database.")

long <- base %>%
  pivot_longer(c(opp_recorded, opp_consumed), names_to = "outcome", values_to = "positive") %>%
  mutate(outcome = factor(outcome,
                          levels = c("opp_recorded", "opp_consumed"),
                          labels = c("Opportunist recorded", "Carcass consumption documented")))

by_site <- long %>% transmute(outcome, group = "Feeding-site type",
                              category = recode(Tipo_recinto,
                                                MULADAR = "Muladar", PUNTO = "Feeding point"),
                              positive)
by_fence <- long %>% transmute(outcome, group = "Fence adaptation",
                               category = recode(Valla_acondicionada, NO = "No", SI = "Yes"),
                               positive)
plot_data <- bind_rows(by_site, by_fence) %>%
  group_by(outcome, group, category) %>%
  summarise(n_positive = sum(positive), n_total = n(),
            proportion = n_positive / n_total, .groups = "drop") %>%
  mutate(group = factor(group, levels = c("Feeding-site type", "Fence adaptation")),
         category = factor(category, levels = c("Muladar", "Feeding point", "No", "Yes")),
         count_label = paste0(n_positive, "/", n_total),
         colour_group = ifelse(category %in% c("Muladar", "No"), "A", "B"))

write.csv(base, file.path(out, "APORTE_detected_consumption_flags.csv"), row.names = FALSE)
write.csv(plot_data, file.path(out, "OBSERVED_consumption_frequencies.csv"), row.names = FALSE)

plot <- ggplot(plot_data, aes(x = category, y = proportion, fill = colour_group)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = paste0(scales::percent(proportion, accuracy=1), "\n(", count_label, ")")),
            vjust = -0.36, size = 3.8, lineheight = 0.95) +
  facet_grid(outcome ~ group, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = scales::label_percent(), limits = c(0, 1.13),
                     breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values = c("A" = "#2B5D67", "B" = "#C47742")) +
  labs(title = "Documented facultative-scavenger activity at carcass provisioning sites",
       x = NULL, y = "Proportion of provisioning events",
       caption = "Values are descriptive event-level frequencies (positive / monitored provisioning events).\nRepeated provisioning events were nested within hunting estates; non-detection does not demonstrate absence.") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14),
        strip.background = element_rect(fill = "#EEF1F2", colour = NA),
        strip.text = element_text(face = "bold", size = 11),
        panel.grid.major.y = element_line(colour = "#E4E9EB", linewidth = 0.35),
        axis.text.x = element_text(size = 11),
        legend.position = "none",
        plot.caption = element_text(size = 9, colour = "#53616A", hjust = 0),
        plot.margin = margin(15, 15, 15, 15))

print(plot)
for (ext in c("pdf", "png", "tiff")) {
  target <- file.path(out, paste0("Figure_documented_consumption.", ext))
  if (ext == "pdf") {
    ggsave(target, plot = plot, width = 11, height = 7, units = "in", device = "pdf")
  } else if (ext == "png") {
    ggsave(target, plot = plot, width = 11, height = 7, units = "in", dpi = 600, bg = "white")
  } else {
    ggsave(target, plot = plot, width = 11, height = 7, units = "in", dpi = 600,
           bg = "white", compression = "lzw")
  }
}
cat("Saved Figure_documented_consumption in:", out, "\n")
print(plot_data)
