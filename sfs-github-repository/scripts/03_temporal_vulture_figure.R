################################################################################
# FIGURE 4 | WITHIN-SEASON VARIATION IN VULTURE USE
# Two-panel, journal-style figure (English labels).
#
# A: Time to first recorded vulture (Gamma Bayesian model, in hours)
# B: Number of vultures per recorded detection (negative-binomial Bayesian model)
#
# No model fitting. Uses existing saved posterior predictions and model RDS files.
# Model curves are predictions for one reference management configuration;
# pale points are observed data from ALL configurations and are not adjusted.
#
# Run with Source in RStudio.
################################################################################

needed <- c("ggplot2", "dplyr", "readr", "patchwork", "scales", "lubridate", "brms")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing packages: install.packages(c(",
       paste(sprintf('"%s"', missing), collapse = ", "), "))")
}
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
  library(scales)
  library(lubridate)
  library(brms)
})

# EDIT THIS PATH IF NECESSARY
project_dir <- Sys.getenv("SFS_DATA_DIR", unset=".")
prediction_file <- file.path(project_dir, "FIGURAS_NATURE_V3",
                             "Temporal_predictions.csv")
model_dir <- Sys.getenv("SFS_MODEL_DIR", unset=file.path(project_dir, "models"))
time_model_file <- file.path(model_dir,
                             "modelo_tiempo_buitres_bayes_revisado.rds")
number_model_file <- file.path(model_dir,
                               "modelo_abundancia_bayes_revisado.rds")

required_files <- c(prediction_file, time_model_file, number_model_file)
if (any(!file.exists(required_files))) {
  stop("Cannot find required files:\n",
       paste(required_files[!file.exists(required_files)], collapse = "\n"),
       "\nCheck project_dir and that V3 figures have been generated.")
}

output_dir <- file.path(project_dir, "FIGURAS_NATURE_V3")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 1. SAVED POSTERIOR PREDICTIONS -----------------------------------------------
# The V3 script already generated these values from the corrected brms fits.
pred <- readr::read_csv(prediction_file, show_col_types = FALSE) %>%
  mutate(Date = as.Date(Date)) %>%
  filter(model %in% c("vulture_t", "vulture_n")) %>%
  mutate(
    Response = if_else(model == "vulture_t", "A", "B"),
    Estimate = if_else(model == "vulture_t", Estimate / 60, Estimate),
    Q2.5 = if_else(model == "vulture_t", Q2.5 / 60, Q2.5),
    Q97.5 = if_else(model == "vulture_t", Q97.5 / 60, Q97.5)
  )

if (any(!is.finite(pred$Estimate) | !is.finite(pred$Q2.5) |
        !is.finite(pred$Q97.5) | pred$Q2.5 <= 0)) {
  stop("Non-positive or non-finite posterior predictions; check CSV.")
}
stopifnot(all(c("A","B") %in% pred$Response))

# 2. RECOVER ORIGINAL DATES FROM THE DATABASE -------------------------------
# brms stores fecha_sc, but may NOT retain Fecha_aporte in model$data.
# We reconstruct the original standardized date exactly as in the analysis:
# fecha_sc = scale(as.numeric(Fecha_aporte)) across the entire input dataset.
# Matching uses date values, not a guessed row order.
raw <- readr::read_csv2(
  file.path(project_dir, "BASE_DATOS_BUITRES.csv"),
  locale = readr::locale(encoding = "Latin1"),
  show_col_types = FALSE
) %>%
  filter(!is.na(ID_COTO), trimws(as.character(ID_COTO)) != "") %>%
  mutate(Fecha_aporte = lubridate::dmy(Fecha_aporte))

if (sum(!is.na(raw$Fecha_aporte)) < 2L) {
  stop("Insufficient valid provisioning dates in the CSV.")
}

# A constant date offset does not change z-standardization, so using
# absolute Date numbers produces the same scale as days since first date.
raw$fecha_z <- as.numeric(scale(as.numeric(raw$Fecha_aporte)))
reference_dates <- raw %>%
  filter(!is.na(Fecha_aporte), is.finite(fecha_z)) %>%
  distinct(Fecha_aporte, fecha_z) %>%
  arrange(Fecha_aporte)

recover_dates <- function(fit, response_name) {
  md <- fit$data
  if (!response_name %in% names(md)) {
    stop("Response '", response_name, "' not found in saved fit$data.")
  }

  # If the original date survived in model$data, use it directly.
  if ("Fecha_aporte" %in% names(md)) {
    x <- md$Fecha_aporte
    date <- if (inherits(x, "Date")) x else lubridate::dmy(as.character(x))
  } else {
    if (!"fecha_sc" %in% names(md)) {
      stop("The saved model contains neither Fecha_aporte nor fecha_sc.")
    }
    z <- as.numeric(md$fecha_sc)
    if (any(!is.finite(z))) {
      stop("Non-finite fecha_sc in saved model.")
    }
    closest <- vapply(z, function(a) which.min(abs(reference_dates$fecha_z - a)), integer(1))
    discrepancy <- abs(reference_dates$fecha_z[closest] - z)
    if (any(discrepancy > 1e-4)) {
      stop("Cannot match some model dates to BASE_DATOS_BUITRES.csv. ",
           "Maximum standardized difference: ", signif(max(discrepancy), 4), ". ",
           "Check whether model and CSV originate from the same database.")
    }
    date <- reference_dates$Fecha_aporte[closest]
  }

  y <- suppressWarnings(as.numeric(md[[response_name]]))
  out <- data.frame(Date = as.Date(date), Value = y)
  out <- out[!is.na(out$Date) & is.finite(out$Value) & out$Value > 0, , drop = FALSE]
  out
}

# The observed values below are taken from the exact model data (not all CSV
# rows), preserving the n = 59 / n = 88 used in the models.
time_fit <- readRDS(time_model_file)
number_fit <- readRDS(number_model_file)

obs_t <- recover_dates(time_fit, "tiempo_buitres") %>%
  mutate(Value = Value / 60, Response = "A")
obs_n <- recover_dates(number_fit, "N") %>%
  mutate(Response = "B")

observations <- bind_rows(obs_t, obs_n)
if (nrow(obs_t) == 0L || nrow(obs_n) == 0L) {
  stop("No valid observed data found in the saved fits.")
}
message("Matched model dates: ", nrow(obs_t), " arrival times, ",
        nrow(obs_n), " vulture detections.")

# 3. OPTIONAL MONTHLY SAMPLE SIZES FOR THE FIGURE CAPTION ----------------------
month_counts <- observations %>%
  mutate(Month = format(Date, "%Y-%m")) %>%
  count(Response, Month, name = "n_observations")
write.csv(month_counts,
          file.path(output_dir, "Figure_4_monthly_sample_sizes.csv"),
          row.names = FALSE)

write.csv(pred, file.path(output_dir, "Figure_4_temporal_posterior.csv"),
          row.names = FALSE)
write.csv(observations,
          file.path(output_dir, "Figure_4_observed_data.csv"),
          row.names = FALSE)

# 4. THEME AND PANEL BUILDING --------------------------------------------------
TEAL <- "#315A65"
INK <- "#253640"
LIGHT <- "#DDE5E7"
GREY <- "#88969C"

theme_journal <- function() {
  theme_classic(base_size = 12, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = 12.5,
                                colour = INK, margin = margin(b = 8)),
      axis.title = element_text(colour = INK, size = 11.5),
      axis.text = element_text(colour = INK, size = 10),
      axis.line = element_line(colour = INK, linewidth = 0.35),
      axis.ticks = element_line(colour = INK, linewidth = 0.3),
      panel.grid = element_blank(),
      plot.margin = margin(12, 16, 12, 12)
    )
}

# Explicit English month names, independently of the computer locale.
month_label <- function(x) {
  paste(month.abb[as.integer(format(x, "%m"))], format(x, "%Y"))
}

plot_panel <- function(panel_id, title, ylabel) {
  pd <- pred %>% filter(Response == panel_id)
  od <- observations %>% filter(Response == panel_id)

  ggplot() +
    # Reference configuration prediction: uncertainty ribbon
    geom_ribbon(data = pd,
                aes(x = Date, ymin = Q2.5, ymax = Q97.5),
                fill = TEAL, alpha = .17) +
    # Actual observations over all available site configurations
    geom_point(data = od,
               aes(x = Date, y = Value),
               colour = GREY, alpha = .31, size = 1.55,
               position = position_jitter(width = 0.7, height = 0),
               show.legend = FALSE) +
    geom_line(data = pd, aes(x = Date, y = Estimate),
              colour = TEAL, linewidth = 1.05) +
    # Log scale aids readability across very large observed ranges and CrIs
    scale_y_log10(labels = scales::label_number(accuracy = 0.1),
                  breaks = scales::breaks_log(n = 5)) +
    scale_x_date(date_breaks = "1 month", labels = month_label,
                 expand = expansion(mult = c(.01, .025))) +
    labs(title = title, x = "Carcass provisioning date", y = ylabel) +
    theme_journal() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
}

panel_a <- plot_panel("A",
                      "A  Time to first vulture arrival",
                      "Arrival time (hours, log scale)")
panel_b <- plot_panel("B",
                      "B  Number of vultures detected",
                      "Vultures per detection (log scale)")

combined <- (panel_a | panel_b) +
  patchwork::plot_annotation(
    title = "Within-season variation in vulture use",
    subtitle = paste(
      "Posterior means and 95% credible intervals,",
      "with observed records shown as grey points"
    ),
    caption = paste(
      "Predictions are conditional on a carcass deposition point with an",
      "adapted fence and no livestock; standardized covariates are set to zero,",
      "and estate-level random effects to zero.",
      "Grey points represent unadjusted observations from all site configurations.",
      "Panel B is per recorded vulture detection, not total unique individuals",
      "or an estimate of population abundance."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", colour = INK, size = 15),
      plot.subtitle = element_text(size = 10.5, colour = INK),
      plot.caption = element_text(size = 9, colour = "#52636B", hjust = 0,
                                  margin = margin(t = 12)),
      plot.margin = margin(16, 16, 12, 16)
    )
  )

# 5. EXPORT --------------------------------------------------------------------
base_name <- file.path(output_dir, "Figure_4_temporal_vulture_use")
ggsave(paste0(base_name, ".pdf"), combined,
       width = 12.7, height = 6.2, units = "in", bg = "white")
ggsave(paste0(base_name, ".png"), combined,
       width = 12.7, height = 6.2, units = "in",
       dpi = 600, bg = "white")
ggsave(paste0(base_name, ".tiff"), combined,
       width = 12.7, height = 6.2, units = "in",
       dpi = 600, compression = "lzw", bg = "white")

expected <- paste0(base_name, c(".pdf", ".png", ".tiff"))
stopifnot(all(file.exists(expected)))

print(combined)
message("\nFigure 4 saved successfully:\n", paste(expected, collapse = "\n"))
message("Observed model data: ", nrow(obs_t), " arrival records and ",
        nrow(obs_n), " vulture detections.\n",
        "See Figure_4_monthly_sample_sizes.csv for monthly sample sizes.")
