# =============================================================================
# 01_build_equine_model.R
#
# Purpose:
#   Build an equine auditory perceptual weighting curve ("H-weighting"),
#   analogous in spirit to human A-weighting, from a digitized horse
#   audiogram (threshold vs. frequency).
#
# Input:
#   WebPlotDigitizerHorseAudiogram.csv
#     Two columns, NO header: frequency (kHz), threshold (dB)
#
# Output:
#   equine_weighting_model.csv
#     Columns: frequency_hz, threshold_db, H_dB
#     A smooth lookup table (regular in log10-frequency) that can later be
#     interpolated onto ANY LTSA's frequency bins (see script 02).
#
# Method:
#   1. Read + clean digitized points; convert kHz -> Hz
#   2. Collapse duplicate/near-duplicate frequencies (digitizing jitter)
#   3. Interpolate (linear) in log10(frequency) onto a regular grid
#   4. Light loess smoothing in log10(frequency) to remove residual
#      digitizing "overshoot" (small manual-digitizing wiggles), while
#      preserving the real U-shape of the audiogram
#   5. Compute H(f) = Tmin - T(f)
#        Tmin = best (lowest) threshold across the whole curve
#        T(f) = threshold at frequency f
#      -> H(f) = 0 dB at the most sensitive frequency, increasingly
#         negative (attenuated) elsewhere - same shape convention as the
#         human A-weighting curve (which is 0 dB near 1-4 kHz and negative
#         elsewhere).
# =============================================================================

library(tidyverse)

# ---- Parameters you may want to tune ---------------------------------------

input_file   <- "WebPlotDigitizerHorseAudiogram.csv"
output_file  <- "equine_weighting_model.csv"
n_bins       <- 300   # log-frequency bins used to collapse digitizing jitter
grid_n       <- 1000  # resolution of the final interpolated model
loess_span   <- 0.08  # smoothing span (larger = smoother, more flattening)

# ---- 1. Read raw digitized data --------------------------------------------

raw <- read_csv(
  input_file,
  col_names = c("freq_khz", "threshold_db"),
  show_col_types = FALSE
)

audiogram <- raw %>%
  mutate(freq_hz = freq_khz * 1000) %>%
  filter(freq_hz > 0, is.finite(threshold_db)) %>%
  arrange(freq_hz)

# ---- 2. Collapse duplicate / near-duplicate frequencies --------------------
# WebPlotDigitizer traces frequently produce several points essentially on
# top of one another in x (frequency) with slightly different y (threshold),
# and occasionally a few points that are locally non-monotonic in x. This is
# digitizing jitter, not signal, and is the source of the "overshoot" you
# saw. Binning in log10(f) and taking the median threshold per bin fixes
# both problems (removes jitter, restores monotonicity in x) in one step.

audiogram_binned <- audiogram %>%
  mutate(log10_f = log10(freq_hz)) %>%
  mutate(bin = cut(log10_f, breaks = n_bins)) %>%
  group_by(bin) %>%
  summarise(
    log10_f      = median(log10_f),
    threshold_db = median(threshold_db),
    .groups = "drop"
  ) %>%
  arrange(log10_f)

# ---- 3. Interpolate (linear) on a regular log10(f) grid --------------------

log10_f_grid <- seq(
  min(audiogram_binned$log10_f),
  max(audiogram_binned$log10_f),
  length.out = grid_n
)

interp <- approx(
  x      = audiogram_binned$log10_f,
  y      = audiogram_binned$threshold_db,
  xout   = log10_f_grid,
  method = "linear",
  rule   = 2  # clamp at edges rather than extrapolate
)

model <- tibble(
  log10_f          = interp$x,
  freq_hz          = 10^interp$x,
  threshold_db_raw = interp$y
)

# ---- 4. Smooth to remove residual digitizing overshoot ---------------------

smoothed <- loess(
  threshold_db_raw ~ log10_f,
  data   = model,
  span   = loess_span,
  degree = 2
)

model <- model %>%
  mutate(threshold_db = predict(smoothed, newdata = model))

# ---- 5. Compute H-weighting -------------------------------------------------

Tmin <- min(model$threshold_db)

equine_model <- model %>%
  transmute(
    frequency_hz = freq_hz,
    threshold_db = threshold_db,
    H_dB         = Tmin - threshold_db
  )

write_csv(equine_model, output_file)

# ---- Quick check plot -------------------------------------------------------

p_check <- ggplot() +
  geom_point(
    data = audiogram, aes(freq_hz, threshold_db),
    color = "grey60", alpha = 0.4, size = 0.8
  ) +
  geom_line(
    data = equine_model, aes(frequency_hz, threshold_db),
    color = "darkorange", linewidth = 1
  ) +
  scale_x_log10() +
  labs(
    title = "Raw digitized audiogram (grey) vs. smoothed model (orange)",
    x = "Frequency (Hz, log scale)", y = "Threshold (dB)"
  ) +
  theme_minimal()

p_hweight <- ggplot(equine_model, aes(frequency_hz, H_dB)) +
  geom_line(color = "steelblue", linewidth = 1) +
  scale_x_log10() +
  labs(
    title = "Equine H-Weighting Curve",
    x = "Frequency (Hz, log scale)", y = "H(f) weighting (dB)"
  ) +
  theme_minimal()

print(p_check)
print(p_hweight)
