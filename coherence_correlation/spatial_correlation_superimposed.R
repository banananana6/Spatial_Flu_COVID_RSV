# spatial coherence test using time-series in multiple seasons

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
library(ncf)
library(reshape2)

setwd("") # add directory

## influenza/rsv

# seasons <- list(
#   "Winter 2023-2024" = list(start = "2023-07-01", end = "2024-05-01"),
#   "Winter 2024-2025" = list(start = "2023-07-01", end = "2024-05-01"),
#   "Winter 2025-2026" = list(start = "2023-07-01", end = "2024-05-01")
# )


## covid

# seasons <- list(
#   "Winter 2023-2024" = list(start = "2023-09-01", end = "2024-03-01"),
#   "Winter 2024-2025" = list(start = "2024-09-01", end = "2025-03-20"),
#   "Winter 2025-2026" = list(start = "2025-09-01", end = "2026-03-01")
# )

seasons <- list(
  "Summer 2023" = list(start = "2023-03-20", end = "2023-09-01"),
  "Summer 2024" = list(start = "2024-03-01", end = "2024-09-01"),
  "Summer 2025" = list(start = "2025-03-20", end = "2025-09-01")
)


locs <- read.csv("ED.csv") %>%
  dplyr::filter(!is.na(percent_visits_covid)) %>%
  distinct(geography, hsa, hsa_counties, hsa_nci_id)

ed_full <- read.csv("ED.csv") %>%
  distinct(geography, hsa, hsa_counties, hsa_nci_id, week_end, .keep_all = TRUE) %>%
  dplyr::select(-dplyr::any_of(c("county", "fips"))) %>%
  right_join(locs, by = c("geography", "hsa", "hsa_counties", "hsa_nci_id")) %>%
  mutate(week_end = as.Date(week_end),
         hsa_nci_id = as.numeric(hsa_nci_id)) %>%
  arrange(hsa_nci_id, week_end)

hsa_geo <- read.csv("HSA_GEO.csv") %>%
  filter(hsa_id %in% locs$hsa_nci_id) %>%
  filter(!(hsa_id %in% c(821, 822, 823))) %>%
  filter(total_population > 5000) %>%
  mutate(hsa_id = as.numeric(hsa_id)) %>%
  arrange(hsa_id)

lat <- hsa_geo$weighted_lat
long <- hsa_geo$weighted_lon

ncf_test <- function(long, lat, df) {
  test <- Sncf(long, lat, df, type= "boot", resamp = 1000, latlon = TRUE, na.rm  = TRUE)
  tibble(
    cbar = test$real$cbar,
    x = test$real$predicted$x[1, ],
    y = test$real$predicted$y[1, ]
  )
}

results <- map_dfr(names(seasons), function(s) {
  d <- seasons[[s]]
  
  ed_s <- ed_full %>%
    filter(hsa_nci_id %in% hsa_geo$hsa_id) %>%
    filter(!(hsa_nci_id %in% c(821, 822, 823))) %>%
    filter(week_end >= as.Date(d$start) & week_end <= as.Date(d$end)) %>%
    arrange(hsa_nci_id)
  
  z_mat <- acast(ed_s, hsa_nci_id ~ week_end, value.var = "percent_visits_covid")
  
  keep <- hsa_geo$hsa_id %in% as.numeric(rownames(z_mat))
  
  ncf_test(long[keep], lat[keep], z_mat) %>%
    mutate(season=s)
})

cbar_df <- results %>%
  distinct(season, cbar)

ggplot(results, aes(x = x, y = y, color = season)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(data = cbar_df, aes(yintercept = cbar, color = season),
             linetype = "dotted", linewidth = 0.6) +
  scale_color_brewer(palette = "Dark2", name = "Season") +
  scale_y_continuous(limits = c(-0.1, 1)) +
  labs(
    x = "Distance (km)",
    y = "Correlation",
    title = "Spatial Coherence — COVID"
  )