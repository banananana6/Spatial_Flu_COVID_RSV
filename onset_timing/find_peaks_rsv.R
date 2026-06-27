# find peaks for each RSV time-series

library(tidyr)
library(dplyr)
library(lubridate)

setwd("") # add directory

locs <- read.csv("ED.csv") %>%
  dplyr::filter(!is.na(percent_visits_rsv)) %>%
  distinct(geography, hsa, hsa_counties, hsa_nci_id)

ed <- read.csv("ED.csv") %>%
  distinct(geography, hsa, hsa_counties, hsa_nci_id, week_end, .keep_all = TRUE) %>%
  dplyr::select(-dplyr::any_of(c("county", "fips"))) %>%
  right_join(locs, by = c("geography", "hsa", "hsa_counties", "hsa_nci_id")) %>%
  mutate(week_end = as.Date(week_end)) %>%
  arrange(hsa_nci_id, week_end)

seasons <- list(
  "2023-2024" = list(start = "2023-08-01", end = "2024-05-01"),
  "2024-2025" = list(start = "2024-08-01", end = "2025-05-01"),
  "2025-2026" = list(start = "2025-08-01", end = "2026-05-01")
)

results <- purrr::map_dfr(names(seasons), function(season_name) {
  s <- seasons[[season_name]]
  ed %>%
    filter(week_end >= as.Date(s$start), week_end < as.Date(s$end)) %>%
    group_by(hsa_nci_id, hsa, hsa_counties, geography) %>%
    slice_max(percent_visits_rsv, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(season = season_name, hsa_nci_id, hsa, hsa_counties, geography,
              peak_date = week_end,
              peak_size= percent_visits_rsv)
})

write.csv(results, "peaks_rsv.csv", row.names = FALSE)
