# find peaks for each COVID time-series

library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times", base_size = 14))
library(lubridate)
library(data.table)
library(purrr)
library(MMWRweek)
library(segmented)
library(ggpubr)
library(zoo)
library(pracma)

setwd("") # add directory

locs <- read.csv("ED.csv") %>%
  dplyr::filter(!is.na(percent_visits_smoothed_covid)) %>% # covid
  distinct(geography, hsa, hsa_counties, hsa_nci_id)

ed <- read.csv("ED.csv") %>%
  distinct(geography, hsa, hsa_counties, hsa_nci_id, week_end, .keep_all = TRUE) %>%
  dplyr::select(-dplyr::any_of(c("county", "fips"))) %>%
  right_join(locs, by = c("geography", "hsa", "hsa_counties", "hsa_nci_id")) %>%
  mutate(week_end = as.Date(week_end)) %>%
  arrange(hsa_nci_id, week_end)

ed <- ed %>%
  group_by(hsa_nci_id) %>%
  arrange(week_end) %>%
  ungroup()

# Helper: compute true prominence for each peak
prominence <- function(x, peak_idx) {
  n <- length(x)
  prom <- numeric(length(peak_idx))
  for (i in seq_along(peak_idx)) {
    p <- peak_idx[i]
    higher <- peak_idx[x[peak_idx] > x[p]]
    left_bound  <- max(c(1, higher[higher < p]))
    right_bound <- min(c(n, higher[higher > p]))
    prom[i] <- x[p] - max(min(x[left_bound:p]), min(x[p:right_bound]))
  }
  prom
}

find_peaks_prominent <- function(x, min_prom_abs = 0.3, min_dist = 8) {
  candidates <- findpeaks(x, minpeakheight = 0.5, minpeakdistance = min_dist, npeaks = 0)

  # --- sliding window pass to catch broad hills ---
  window_size <- 10
  n <- length(x)
  window_peaks <- c()
  for (i in (window_size+1):(n-window_size)) {
    if ((!is.na(x[i]) && x[i] == max(x[(i-window_size):(i+window_size)], na.rm = TRUE)) && x[i] > 0.5) {
      window_peaks <- c(window_peaks, i)
    }
  }

  # --- enforce min_dist on window peaks ---
  if (length(window_peaks) > 1) {
    keep <- c(TRUE, diff(window_peaks) >= min_dist)
    window_peaks <- window_peaks[keep]
  }

  # --- combine with findpeaks candidates ---
  all_idx <- sort(unique(c(if (!is.null(candidates)) candidates[, 2] else NULL, window_peaks)))
  if (length(all_idx) == 0) return(NULL)

  # --- build matrix and apply prominence filter ---
  mat <- matrix(cbind(x[all_idx], all_idx, all_idx - window_size, all_idx + window_size), ncol = 4)
  prom <- prominence(x, all_idx)
  mat[prom >= min_prom_abs, , drop = FALSE]
}

peaks <- ed %>%
  group_by(hsa_nci_id) %>%
  arrange(week_end, .by_group = TRUE) %>%
  reframe(
    peak_date = {
      p <- find_peaks_prominent(percent_visits_smoothed_covid, min_prom_abs = 0.3, min_dist = 8) # covid
      if (is.null(p) || nrow(p) == 0) {
        as.Date(character(0))
      } else {
        as.Date(week_end[p[, 2]], origin = "1970-01-01")
      }
    },
    peak_size= {
      p <- find_peaks_prominent(percent_visits_smoothed_covid, min_prom_abs = 0.3, min_dist = 8) # covid
      if (is.null(p) || nrow(p) == 0) {
        numeric(0)
      } else {
        p[, 1]
      }
    }
  ) %>%
  mutate(peak_date = format(peak_date, "%Y-%m-%d"))

write.csv(peaks, "peaks_covid.csv", row.names = FALSE)

