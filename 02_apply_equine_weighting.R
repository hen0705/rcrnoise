# =============================================================================
# 02_apply_equine_weighting.R
#
# Purpose:
#   Apply the equine H-weighting curve (built by 01_build_equine_model.R) to
#   one or more 24-hr LTSA CSVs, producing per-file:
#     - a broadband "equine-weighted vs. unweighted" time series (CSV + PNG)
#     - a full equine-weighted spectrogram (CSV + PNG heatmap)
#   plus a combined comparison across all files processed in one run.
#
# Usage:
#   Batch mode (processes every *.csv in input_dir):
#       Rscript 02_apply_equine_weighting.R
#
#   Single-file mode (overrides input_dir for one specific file):
#       Rscript 02_apply_equine_weighting.R path/to/one_file.csv
#
# Expected folder layout (folder-based, not filename-pattern-based, so this
# never mistakes the model file or its own outputs for LTSA input):
#   data/ltsa_raw/        <- put ALL your LTSA CSVs here, any file names
#   data/ltsa_processed/  <- this script writes CSV outputs here
#   figures/              <- this script writes PNG outputs here
#   equine_weighting_model.csv   <- from script 01, lives at repo root
#
# Confirmed input LTSA format (from your actual file):
#   CSV, wide format:
#     - column 1: timestamps like "2026-06-29 17:25:02"
#     - remaining columns: frequency in Hz, 10 Hz spacing, 0-24000 Hz
#       (a 0 Hz / DC bin is present and is dropped - see step 2b)
#     - cell values: spectral level in dB (dBFS-derived)
#
# Every output CSV keeps a datetime column - the per-file spectrogram and
# timeseries CSVs have one, and the combined timeseries CSV has both
# datetime and source_file.
#
# Math (matches your notebook's Python weighting step exactly:
# weight_linear = 10**(weight_dB/10), then sum, then 10*log10):
#     P(f)        = 10^(level_dB(f) / 10)
#     H_linear(f) = 10^(H_dB(f)      / 10)
#     Pweighted(f)= P(f) * H_linear(f)
#     Ptotal      = sum_f Pweighted(f)
#     L_equine    = 10 * log10(Ptotal)
#
# NOTE ON dBFS: these levels are relative to full-scale, not an absolute
# SPL reference - "L_equine_dB" is a relative index unless/until the
# recording chain is calibrated to a known SPL reference.
# =============================================================================

library(tidyverse)

# ---- Paths (edit if your folder layout differs) -----------------------------

input_dir  <- "data/ltsa_raw"
output_dir <- "data/ltsa_processed"
figures_dir <- "figures"
model_file <- "equine_weighting_model.csv"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 0. Load the equine model (once, shared across all files) --------------

equine_model <- read_csv(model_file, show_col_types = FALSE)

get_equine_H <- function(freq_hz_query, model = equine_model) {
  approx(
    x    = log10(model$frequency_hz),
    y    = model$H_dB,
    xout = log10(freq_hz_query),
    rule = 2
  )$y
}

# ---- 1. Core per-file processing function -----------------------------------

process_one_ltsa <- function(file_path, equine_model, output_dir, figures_dir) {

  base_name <- tools::file_path_sans_ext(basename(file_path))
  message("Processing: ", basename(file_path))

  ltsa_wide_raw <- read_csv(file_path, show_col_types = FALSE)

  ltsa_wide <- ltsa_wide_raw %>%
    rename(datetime = 1) %>%
    mutate(datetime = ymd_hms(datetime))
    # If this produces NA timestamps for one of your files, that file's
    # date format differs - swap in the matching lubridate parser for it.

  ltsa_long <- ltsa_wide %>%
    pivot_longer(-datetime, names_to = "frequency_raw", values_to = "level_db") %>%
    mutate(frequency_hz = as.numeric(frequency_raw)) %>%
    select(datetime, frequency_hz, level_db) %>%
    filter(frequency_hz > 0)   # drop the 0 Hz / DC bin - see header note

  freq_lookup <- tibble(frequency_hz = unique(ltsa_long$frequency_hz)) %>%
    mutate(H_dB = get_equine_H(frequency_hz))

  ltsa_long <- ltsa_long %>%
    left_join(freq_lookup, by = "frequency_hz")

  # -- Weighted spectrogram (per time x frequency bin) -----------------------

  ltsa_weighted_spectrogram <- ltsa_long %>%
    mutate(level_db_equine_weighted = level_db + H_dB)

  write_csv(
    ltsa_weighted_spectrogram,
    file.path(output_dir, paste0(base_name, "_equine_weighted_spectrogram.csv"))
  )

  # -- Broadband time series: unweighted vs. equine-weighted -----------------

  ltsa_timeseries <- ltsa_long %>%
    mutate(
      P_linear   = 10^(level_db / 10),
      H_linear   = 10^(H_dB / 10),
      P_weighted = P_linear * H_linear
    ) %>%
    group_by(datetime) %>%
    summarise(
      L_unweighted_dB = 10 * log10(sum(P_linear, na.rm = TRUE)),
      L_equine_dB     = 10 * log10(sum(P_weighted, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    arrange(datetime)

  write_csv(
    ltsa_timeseries,
    file.path(output_dir, paste0(base_name, "_equine_weighted_timeseries.csv"))
  )

  # -- Per-file figures, saved to disk -----------------------------------------

  p_heatmap <- ltsa_weighted_spectrogram %>%
    ggplot(aes(datetime, frequency_hz, fill = level_db_equine_weighted)) +
    geom_raster() +
    scale_fill_viridis_c(name = "dB (equine-weighted)") +
    scale_y_continuous(labels = scales::label_number(scale = 1e-3, suffix = "")) +
    labs(title = paste("Equine-Weighted LTSA:", base_name),
         x = "Time", y = "Frequency (kHz)") +
    theme_minimal()

  p_series <- ltsa_timeseries %>%
    pivot_longer(c(L_unweighted_dB, L_equine_dB),
                 names_to = "weighting", values_to = "level_db") %>%
    mutate(weighting = recode(weighting,
      L_unweighted_dB = "Unweighted", L_equine_dB = "Equine-weighted")) %>%
    ggplot(aes(datetime, level_db, color = weighting)) +
    geom_line(linewidth = 0.4) +
    labs(title = paste("Acoustic level comparison:", base_name),
         x = "Time", y = "Relative integrated level (dB)", color = NULL) +
    theme_minimal()

  ggsave(file.path(figures_dir, paste0(base_name, "_heatmap.png")),
         p_heatmap, width = 10, height = 5, dpi = 300)
  ggsave(file.path(figures_dir, paste0(base_name, "_timeseries.png")),
         p_series, width = 10, height = 4, dpi = 300)

  ltsa_timeseries
}

# ---- 2. Decide which file(s) to run - batch, or single-file override ------

cli_args <- commandArgs(trailingOnly = TRUE)

if (length(cli_args) >= 1) {
  files_to_process <- cli_args[1]
} else {
  files_to_process <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
}

if (length(files_to_process) == 0) {
  stop("No CSV files found in ", input_dir,
       " - drop your LTSA files there, or pass a single file path as an argument.")
}

# ---- 3. Run it, collecting all timeseries results for a combined figure ---

all_timeseries <- files_to_process %>%
  set_names(tools::file_path_sans_ext(basename(.))) %>%
  map(process_one_ltsa, equine_model = equine_model,
      output_dir = output_dir, figures_dir = figures_dir) %>%
  list_rbind(names_to = "source_file")

write_csv(all_timeseries, file.path(output_dir, "all_files_equine_weighted_timeseries.csv"))

# ---- 4. Combined figure across all processed files, saved to disk ---------

p_combined <- all_timeseries %>%
  pivot_longer(c(L_unweighted_dB, L_equine_dB),
               names_to = "weighting", values_to = "level_db") %>%
  mutate(weighting = recode(weighting,
    L_unweighted_dB = "Unweighted", L_equine_dB = "Equine-weighted")) %>%
  ggplot(aes(datetime, level_db, color = weighting)) +
  geom_line(linewidth = 0.4) +
  facet_wrap(~source_file, scales = "free_x") +
  labs(title = "Acoustic level comparison - all files",
       x = "Time", y = "Relative integrated level (dB)", color = NULL) +
  theme_minimal()

ggsave(file.path(figures_dir, "all_files_comparison.png"),
       p_combined, width = 12, height = 6, dpi = 300)

message("Done. Processed ", length(files_to_process), " file(s).")
message("CSV outputs in ", output_dir, "/")
message("Figures in ", figures_dir, "/")

# ---- 5. Sanity-check reminder ------------------------------------------------
# Before trusting any of this beyond exploration: pick 2-3 timestamps from
# ONE file, compute L_equine_dB for them independently in Python (same
# model, same raw levels), and confirm the numbers agree to within
# floating-point rounding. Silent unit mismatches (Hz vs kHz, log10 vs
# natural log, power vs amplitude dB) are the most common way this kind of
# pipeline goes quietly wrong, and won't show up as an R error - only as
# numbers that are wrong by a consistent, easy-to-miss factor.
