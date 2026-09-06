# =============================================================================
# 02_apply_equine_weighting.R
#
# Purpose:
#   Apply the equine H-weighting curve (built by 01_build_equine_model.R) to
#   24-hr LTSA recording SESSIONS, producing per-session:
#     - a broadband "equine-weighted vs. unweighted" time series (CSV + PNG)
#     - a full equine-weighted spectrogram (CSV + PNG heatmap)
#   plus a combined comparison across all sessions processed in one run.
#
# IMPORTANT - split-session files:
#   Your LTSA export splits a recording session at midnight, producing TWO
#   CSVs per 24-hr session, sharing a filename prefix but with different
#   trailing dates, e.g.:
#       20260629_184946_20260629.csv   (session start: 6/29 18:49:46,
#                                        rows for 6/29)
#       20260629_184946_20260630.csv   (same session, rows for 6/30)
#   This script detects that shared "<startdate>_<starttime>" prefix,
#   groups files by it, and concatenates them (sorted by datetime, not by
#   filename, so order is robust either way) into ONE continuous session
#   before computing anything. A file whose name doesn't match that pattern
#   is treated as its own single-file session, so nothing gets silently
#   dropped - if your export ever names things differently, check the
#   `extract_session_id()` function below and adjust the regex.
#
# Usage:
#   Batch mode (groups + processes every *.csv in input_dir into sessions):
#       Rscript 02_apply_equine_weighting.R
#
#   Single-file mode (bypasses session grouping - processes exactly the
#   file you name, on its own):
#       Rscript 02_apply_equine_weighting.R path/to/one_file.csv
#
# Expected folder layout:
#   I:/RCR_Acoustics/ltsa_raw/  <- raw LTSA CSVs live here (external drive,
#                                  NOT committed to the GitHub repo - kept
#                                  out per the earlier size/portability
#                                  discussion)
#   data/ltsa_processed/        <- this script writes CSV outputs here
#                                  (relative to the repo - fine to commit,
#                                  these are much smaller)
#   figures/                    <- this script writes PNG outputs here
#   equine_weighting_model.csv  <- from script 01, lives at repo root
#
# NOTE: input_dir below is an absolute, machine-specific path. This script
# will only run correctly on a machine where that drive letter/path exists
# exactly as written - it will NOT work unmodified on a collaborator's
# machine with a different drive mapping. If that becomes a problem,
# consider an environment variable or a config file instead of hardcoding it.
#
# Confirmed input LTSA format:
#   CSV, wide format:
#     - column 1: timestamps like "2026-06-29 17:25:02"
#     - remaining columns: frequency in Hz, 10 Hz spacing, 0-24000 Hz
#       (a 0 Hz / DC bin is present and is dropped - see step 2b)
#     - cell values: spectral level in dB (dBFS-derived)
#
# Every output CSV keeps a datetime column - the per-session spectrogram
# and timeseries CSVs have one, and the combined timeseries CSV has both
# datetime and session_id.
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

input_dir   <- "I:/RCR_Acoustics/ltsa_raw"
output_dir  <- "data/ltsa_processed"
figures_dir <- "figures"
model_file  <- "equine_weighting_model.csv"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 0. Load the equine model (once, shared across all sessions) -----------

equine_model <- read_csv(model_file, show_col_types = FALSE)

get_equine_H <- function(freq_hz_query, model = equine_model) {
  approx(
    x    = log10(model$frequency_hz),
    y    = model$H_dB,
    xout = log10(freq_hz_query),
    rule = 2
  )$y
}

# ---- 1. Group midnight-split files into sessions ---------------------------

extract_session_id <- function(file_path) {
  bn <- basename(file_path)
  m  <- str_extract(bn, "^\\d{8}_\\d{6}")  # e.g. "20260629_184946"
  if (is.na(m)) tools::file_path_sans_ext(bn) else m
  # Files that don't match the expected pattern fall back to using their
  # own full name as the session id, i.e. they're treated as a
  # single-file session rather than being dropped or erroring out.
}

# ---- 2. Core per-session processing function --------------------------------

process_one_session <- function(file_paths, equine_model, output_dir, figures_dir, session_id) {

  if (length(file_paths) > 1) {
    message("Session ", session_id, ": merging ", length(file_paths), " files -> ",
            paste(basename(file_paths), collapse = ", "))
  } else {
    message("Session ", session_id, ": single file -> ", basename(file_paths))
  }

  ltsa_wide <- file_paths %>%
    map(~ read_csv(.x, show_col_types = FALSE) %>% rename(datetime = 1)) %>%
    bind_rows() %>%
    mutate(datetime = ymd_hms(datetime)) %>%
    arrange(datetime)
    # arrange() here means file order/naming doesn't matter - the merge is
    # always sorted correctly by actual timestamp.
    # If datetime comes back all NA, this file's date format differs from
    # "YYYY-MM-DD HH:MM:SS" - swap in the matching lubridate parser.

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
    file.path(output_dir, paste0(session_id, "_equine_weighted_spectrogram.csv"))
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
    file.path(output_dir, paste0(session_id, "_equine_weighted_timeseries.csv"))
  )

  # -- Per-session figures, saved to disk --------------------------------------

  p_heatmap <- ltsa_weighted_spectrogram %>%
    ggplot(aes(datetime, frequency_hz, fill = level_db_equine_weighted)) +
    geom_raster() +
    scale_fill_viridis_c(name = "dB (equine-weighted)") +
    scale_y_continuous(labels = scales::label_number(scale = 1e-3, suffix = "")) +
    labs(title = paste("Equine-Weighted LTSA:", session_id),
         x = "Time", y = "Frequency (kHz)") +
    theme_minimal()

  p_series <- ltsa_timeseries %>%
    pivot_longer(c(L_unweighted_dB, L_equine_dB),
                 names_to = "weighting", values_to = "level_db") %>%
    mutate(weighting = recode(weighting,
      L_unweighted_dB = "Unweighted", L_equine_dB = "Equine-weighted")) %>%
    ggplot(aes(datetime, level_db, color = weighting)) +
    geom_line(linewidth = 0.4) +
    labs(title = paste("Acoustic level comparison:", session_id),
         x = "Time", y = "Relative integrated level (dB)", color = NULL) +
    theme_minimal()

  ggsave(file.path(figures_dir, paste0(session_id, "_heatmap.png")),
         p_heatmap, width = 10, height = 5, dpi = 300)
  ggsave(file.path(figures_dir, paste0(session_id, "_timeseries.png")),
         p_series, width = 10, height = 4, dpi = 300)

  ltsa_timeseries
}

# ---- 3. Decide which file(s)/session(s) to run ------------------------------

cli_args <- commandArgs(trailingOnly = TRUE)

if (length(cli_args) >= 1) {
  # Explicit single-file override: process exactly this file, on its own,
  # bypassing session grouping entirely.
  sessions <- list(cli_args[1])
  names(sessions) <- tools::file_path_sans_ext(basename(cli_args[1]))
} else {
  all_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(all_files) == 0) {
    stop("No CSV files found in ", input_dir,
         " - drop your LTSA files there, or pass a single file path as an argument.")
  }
  session_ids <- map_chr(all_files, extract_session_id)
  sessions <- split(all_files, session_ids)
}

# ---- 4. Run it, collecting all timeseries results for a combined figure ---

all_timeseries <- sessions %>%
  imap(~ process_one_session(.x, equine_model = equine_model,
                              output_dir = output_dir, figures_dir = figures_dir,
                              session_id = .y)) %>%
  list_rbind(names_to = "session_id")

write_csv(all_timeseries, file.path(output_dir, "all_sessions_equine_weighted_timeseries.csv"))

# ---- 5. Combined figure across all processed sessions, saved to disk -------

p_combined <- all_timeseries %>%
  pivot_longer(c(L_unweighted_dB, L_equine_dB),
               names_to = "weighting", values_to = "level_db") %>%
  mutate(weighting = recode(weighting,
    L_unweighted_dB = "Unweighted", L_equine_dB = "Equine-weighted")) %>%
  ggplot(aes(datetime, level_db, color = weighting)) +
  geom_line(linewidth = 0.4) +
  facet_wrap(~session_id, scales = "free_x") +
  labs(title = "Acoustic level comparison - all sessions",
       x = "Time", y = "Relative integrated level (dB)", color = NULL) +
  theme_minimal()

ggsave(file.path(figures_dir, "all_sessions_comparison.png"),
       p_combined, width = 12, height = 6, dpi = 300)

message("Done. Processed ", length(sessions), " session(s).")
message("CSV outputs in ", output_dir, "/")
message("Figures in ", figures_dir, "/")

# ---- 6. Sanity-check reminder ------------------------------------------------
# Before trusting any of this beyond exploration: pick 2-3 timestamps
# straddling the midnight boundary of a merged session, confirm they're in
# the right order and that no rows were dropped or duplicated across the
# two source files. Then independently compute L_equine_dB for a couple of
# timestamps in Python (same model, same raw levels) and confirm the
# numbers agree to within floating-point rounding.
