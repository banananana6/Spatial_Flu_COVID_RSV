# n-week rolling mean to smooth COVID time series

library(dplyr)
library(lubridate)
library(zoo)

setwd("") # add directory

locs <- read.csv("ED.csv") %>%
  dplyr::filter(!is.na(percent_visits_covid)) %>% 
  distinct(geography, hsa, hsa_counties, hsa_nci_id)

ed_smooth <- read.csv("ED.csv") %>%
  distinct(
    geography,
    hsa,
    hsa_counties,
    hsa_nci_id,
    week_end,
    .keep_all = TRUE
  ) %>%
  dplyr::select(-county, -fips) %>%
  right_join(
    locs,
    by = c("geography", "hsa", "hsa_counties", "hsa_nci_id")
  ) %>%
  mutate(
    week_end = as.Date(week_end),
    hsa_nci_id = as.numeric(hsa_nci_id)
  ) %>%
  arrange(hsa_nci_id, week_end) %>%
  group_by(hsa_nci_id) %>%
  arrange(week_end, .by_group = TRUE) %>%
  mutate(
    covid_smooth7 = rollapply(
      percent_visits_covid,
      width = 7,
      FUN = mean,
      partial = TRUE,
      fill = NA,
      align = "center",
      na.rm = TRUE
    )
  ) %>%
  ungroup()

write.csv(
  ed_smooth,
  "ED_covid_smoothed.csv",
  row.names = FALSE
)