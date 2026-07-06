# spatial coherence test using breakpoint dates
library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times", base_size = 14))
library(lubridate)
library(geosphere)
library(purrr)

setwd("") # add directory

start_date<- "2025-08-01"
end_date <- "2026-05-01"

# start_date <- "2023-03-20" # for pd1/covid
# end_date <- "2023-09-01"
# start_date<- "2023-09-01" # for pd2/covid
# end_date <- "2024-03-01"
# start_date <- "2024-03-01" # for pd3/covid
# end_date <- "2024-09-01"
# start_date <- "2024-09-01" # for pd4/covid
# end_date <- "2025-03-20"
# start_date <- "2025-03-20" # for pd5/covid
# end_date <- "2025-09-01"
# start_date <- "2025-09-01" # for pd6/covid
# end_date <- "2026-03-01"


hsa_geo <- read.csv("HSA_GEO.csv") %>%
  mutate(hsa_id = as.numeric(hsa_id)) %>%
  filter(total_population > 5000) %>%
  arrange(hsa_id)

breakpoints <- read.csv("breakpoints_rsv_30excluded.csv") %>%
  filter(hsa_nci_id != "All") %>%
  mutate(
    hsa_nci_id = as.numeric(hsa_nci_id),
    breakpoint = as.Date(breakpoint, format = "%m/%d/%Y")
  ) %>%
  filter(!(hsa_nci_id %in% c(821, 822, 823))) %>%
  filter(hsa_nci_id %in% hsa_geo$hsa_id) %>%
  filter(!is.na(breakpoint)) %>%
  filter(breakpoint >= start_date & breakpoint <= end_date) %>%
  group_by(hsa_nci_id) %>%
  slice_min(breakpoint, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(hsa_nci_id)

hsa_geo <- hsa_geo %>%
  filter(hsa_id %in% breakpoints$hsa_nci_id) %>%
  filter(weighted_lon>-90) %>% # limit to east coast
  arrange(hsa_id)

breakpoints <- breakpoints %>%
  filter(hsa_nci_id %in% hsa_geo$hsa_id) %>%
  left_join(hsa_geo, by = c("hsa_nci_id" = "hsa_id")) %>%
  arrange(hsa_nci_id)

lat <- hsa_geo$weighted_lat
long <- hsa_geo$weighted_lon

n <- nrow(hsa_geo)

pair_df <- combn(n, 2, simplify = FALSE) %>%
  map_dfr(function(idx) {
    i <- idx[1]; j <- idx[2]
    tibble(
      dist_km = distHaversine(c(long[i], lat[i]), c(long[j], lat[j]))/1000,
      onset_diff_days = abs(as.numeric(breakpoints$breakpoint[i] -breakpoints$breakpoint[j]))
    )
  })

pair_df <- pair_df %>%
  mutate(dist_bin = cut(dist_km, breaks = seq(0, max(dist_km) + 50, by = 50)))

bin_summary <- pair_df %>%
  group_by(dist_bin) %>%
  summarise(
    mean_diff = mean(onset_diff_days),
    se        = sd(onset_diff_days) / sqrt(n()),
    ci_lo     = mean_diff - 1.96 * se,
    ci_hi     = mean_diff + 1.96 * se,
    dist_mid  = mean(dist_km),
    .groups   = "drop"
  )

## smoothed line 
# ggplot(pair_df, aes(x = dist_km, y = onset_diff_days)) +
#   geom_point(alpha = 0.15, size = 0.6, color = "blue") +
#   geom_smooth(method = "gam", formula = y ~ s(x),
#               color = "red", linewidth = 0.9, se = TRUE) +
#   labs(
#     x = "Distance (km)",
#     y = "Onset difference (days)",
#     title = "Pairwise RSV onset difference vs. geographic distance"
#   )

# confidence interval ribbon plot
ggplot(bin_summary, aes(x = dist_mid, y = mean_diff)) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), fill = "steelblue", alpha = 0.4) +
  geom_line(color = "steelblue4", linewidth = 0.9) +
  labs(
    x     = "Distance (km)",
    y     = "Onset difference (days)",
    title = "Pairwise COVID onset difference vs. geographic distance"
  )