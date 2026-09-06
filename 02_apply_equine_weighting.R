# =============================================================================
# 02_apply_equine_weighting.R
#
# Purpose:
#   Apply the equine H-weighting curve (built in 01_build_equine_model.R) to
#   a 24-hr LTSA (long-term spectral average) to compute an equine-
#   perceptually-weighted broadband level - directly analogous to how dBA is
#   computed from an unweighted spectrum using the A-weighting curve.
#
# Input:
#   - equine_weighting_model.csv   (frequency_hz, H_dB)   <- from script 01
#   - your LTSA. This script expects LONG format:
#         datetime      : a time-bin identifier (POSIXct or character)
#         frequency_hz  : numeric, Hz - the LTSA's OWN frequency bins
#                         (these do NOT need to match the audiogram's bins;
#                          they're interpolated onto the model below)
#         level_db      : numeric, dB (whatever reference your LTSA uses,
#                          e.g. dB re 1 uPa^2/Hz) - just be consistent
#
#     If your LTSA is in WIDE format (rows = time, one column per frequency
#     bin, column names are the frequency in Hz), pivot it first - see the
#     commented example below.
#
# Output:
#   ltsa_equine_weighted.csv: one integrated, equine-weighted broadband
#   level per time bin (L_equine_dB), plus (optionally) a full per-bin
#   weighted spectrogram if you want to re-plot the LTSA itself.
#
# ----------------------------------------------------------------------------
# A NOTE ON THE MATH (this differs slightly from the original prompt, and
# it's worth flagging why):
#
#   H(f) as built in script 01 is a difference of two dB values
#   (Tmin - T(f)), so it lives on a DECIBEL (logarithmic) scale, just like
#   the human A-weighting curve A(f) does. You cannot multiply a linear
#   power value by a dB-scale number directly - the units don't match.
#
#   The standard way this is resolved for A-weighting (and what's used here
#   for H-weighting) is:
#       P(f)        = 10^(level_db(f) / 10)   # LTSA power -> linear
#       H_linear(f) = 10^(H_dB(f)      / 10)   # H(f) dB   -> linear factor
#       Pweighted(f)= P(f) * H_linear(f)
#       Ptotal      = sum_f Pweighted(f)
#       L_equine    = 10 * log10(Ptotal)
#
#   This is the same 4-step pattern you described (power -> weight -> sum ->
#   log), just with H(f) converted to linear first so the multiplication is
#   dimensionally valid. It reduces to your dBA-analogue exactly.
# =============================================================================

library(tidyverse)

# ---- 0. Load the equine model -----------------------------------------------

equine_model <- read_csv("equine_weighting_model.csv", show_col_types = FALSE)

# ---- 1. Load your LTSA (EDIT THIS SECTION for your actual file/format) -----

# --- Example: LTSA already in long format ------------------------------------
# ltsa_long <- read_csv("my_24hr_ltsa.csv") %>%
#   rename(
#     datetime     = your_time_column,
#     frequency_hz = your_freq_column,
#     level_db     = your_level_column
#   )

# --- Example: LTSA in wide format (time in rows, freq in columns, column ---
#     names are frequency in Hz, e.g. "100", "125", "160", ...) ---------------
# ltsa_wide <- read_csv("my_24hr_ltsa.csv")
# ltsa_long <- ltsa_wide %>%
#   pivot_longer(
#     -datetime,
#     names_to  = "frequency_hz",
#     values_to = "level_db"
#   ) %>%
#   mutate(frequency_hz = as.numeric(frequency_hz))

# ---- 2. Interpolate H(f) onto the LTSA's own frequency bins ----------------
# Interpolation happens in log10(frequency), consistent with how the model
# itself was built. rule = 2 clamps at the edges (holds the nearest
# audiogram-derived value) rather than extrapolating past the frequency
# range the horse audiogram actually covered - edit this if you'd rather
# treat out-of-range bins as inaudible (see commented alternative below).

get_equine_H <- function(freq_hz_query, model = equine_model) {
  approx(
    x    = log10(model$frequency_hz),
    y    = model$H_dB,
    xout = log10(freq_hz_query),
    rule = 2
  )$y

  # Alternative: treat frequencies outside the digitized audiogram's range
  # as effectively inaudible instead of clamping, e.g. by setting H to a
  # large negative number (e.g. -80) for freq_hz_query below/above the
  # model's min/max frequency_hz. Uncomment/adapt if that better matches
  # your assumptions about the animal's true audible range.
}

# ---- 3. Core weighting function ---------------------------------------------

apply_equine_weighting <- function(ltsa_long) {
  ltsa_long %>%
    mutate(H_dB = get_equine_H(frequency_hz)) %>%
    mutate(
      P_linear   = 10^(level_db / 10),  # LTSA power, linear
      H_linear   = 10^(H_dB / 10),      # equine sensitivity factor, linear
      P_weighted = P_linear * H_linear
    ) %>%
    group_by(datetime) %>%
    summarise(
      P_total_weighted = sum(P_weighted, na.rm = TRUE),
      L_equine_dB       = 10 * log10(P_total_weighted),
      .groups = "drop"
    )
}

# ---- 4. (Optional) full weighted spectrogram, per frequency bin -----------
# Useful if you want to re-plot the LTSA itself with equine weighting
# applied, rather than only the single integrated broadband number per time
# bin. Because dB + dB is the same operation as (linear power * 10^(H/10))
# converted back to dB, this per-bin version is just simple addition -
# no need for the linear round-trip unless you're also summing across freq.

weight_spectrogram <- function(ltsa_long) {
  ltsa_long %>%
    mutate(
      H_dB                       = get_equine_H(frequency_hz),
      level_db_equine_weighted   = level_db + H_dB
    )
}

# ---- 5. Run it ---------------------------------------------------------------

# ltsa_weighted <- apply_equine_weighting(ltsa_long)
# write_csv(ltsa_weighted, "ltsa_equine_weighted.csv")

# ltsa_weighted_spectrogram <- weight_spectrogram(ltsa_long)
# write_csv(ltsa_weighted_spectrogram, "ltsa_equine_weighted_spectrogram.csv")

# ---- 6. Quick plot -----------------------------------------------------------

# ggplot(ltsa_weighted, aes(datetime, L_equine_dB)) +
#   geom_line(color = "firebrick") +
#   labs(
#     title = "24-hr Equine-Weighted Broadband Level",
#     x = "Time", y = "H-weighted Level (dB)"
#   ) +
#   theme_minimal()
