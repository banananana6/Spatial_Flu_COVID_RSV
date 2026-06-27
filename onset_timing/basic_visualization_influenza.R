# visualize flu time-series, with peaks and onsets displayed

library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times", base_size = 14))
library(lubridate)
library(data.table)
library(purrr)
library(MMWRweek)

setwd("") # add directory

start_date <- "2022-10-14"
end_date <- "2026-03-01"

# start_date <- "2023-07-01"
# end_date <- "2024-05-01"
# start_date <- "2024-08-01"
# end_date <- "2025-05-01"
# start_date <- "2025-08-01"
# end_date <- "2026-05-01"

locs <- read.csv("ED.csv") %>%
  dplyr::filter(!is.na(percent_visits_influenza)) %>%
  distinct(geography, hsa, hsa_counties, hsa_nci_id)

locs.all <- read.csv("ED.csv") %>%
  distinct(geography, hsa, hsa_counties, hsa_nci_id)

ed <- read.csv("ED.csv") %>%
  distinct(geography, hsa, hsa_counties, hsa_nci_id, week_end, .keep_all = TRUE) %>%
  dplyr::select(-dplyr::any_of(c("county", "fips"))) %>%
  right_join(locs, by = c("geography", "hsa", "hsa_counties", "hsa_nci_id")) %>%
  mutate(
    week_end = as.Date(week_end),
    hsa_nci_id = as.numeric(hsa_nci_id),
    panel_label = paste0("(", hsa_nci_id, ") ", hsa)
  ) %>%
  arrange(hsa_nci_id, week_end) %>%
  filter(week_end > ymd(start_date) & week_end < ymd(end_date))

peaks <- read.csv("peaks_influenza.csv", stringsAsFactors = FALSE) %>%
  mutate(
    hsa_nci_id = as.numeric(hsa_nci_id),
    peak_date = as.Date(peak_date, format = "%m/%d/%Y")
  ) %>%
  filter(!is.na(hsa_nci_id), !is.na(peak_date))

breakpoints <- read.csv("breakpoints_influenza.csv") %>%
  filter(!is.na(breakpoint)) %>%
  mutate(
    breakpoint = as.Date(breakpoint, format = "%m/%d/%Y"),
    hsa_nci_id = as.numeric(hsa_nci_id)
  )

colors <- c("COVID" = "red", "Flu" = "blue", "RSV" = "green")

# selected_ids <- c(814, 806, 662, 525, 813, 531, 675, 797, 805, 626) # bottom10
selected_ids <- c(287, 408, 453, 22, 274, 153, 688, 83, 66, 16)     # top10
# selected_ids <- c(378, 573, 294, 773, 343, 362, 364, 389, 46, 92)     # median10

# one shared lookup: exactly one panel_label per hsa_nci_id
hsa_labels <- ed %>%
  distinct(hsa_nci_id, panel_label) %>%
  group_by(hsa_nci_id) %>%
  slice(1) %>%
  ungroup()

vlines_df <- peaks %>%
  filter(hsa_nci_id %in% selected_ids) %>%
  group_by(hsa_nci_id, season) %>%
  slice_max(order_by = peak_size, n = 1) %>%
  ungroup() %>%
  filter(peak_date >= as.Date(start_date) & peak_date <= as.Date(end_date)) %>%
  left_join(hsa_labels, by = "hsa_nci_id")

vlines_bp <- breakpoints %>%
  filter(hsa_nci_id %in% selected_ids) %>%
  filter(breakpoint >= as.Date(start_date) & breakpoint <= as.Date(end_date)) %>%
  left_join(hsa_labels, by = "hsa_nci_id")

g0 <- ggplot(subset(ed, hsa_nci_id %in% selected_ids)) +
  
  geom_line(aes(week_end, percent_visits_influenza, color = "Flu"), lwd = 0.8) +
  
  geom_vline(
    data = vlines_df,
    aes(xintercept = peak_date),
    linetype = "dotted",
    color = "gray",
    linewidth = 0.8
  ) +
  geom_vline(
    data = vlines_bp,
    aes(xintercept = breakpoint),
    linetype = "dashed",
    color = "orange",
    linewidth = 0.8
  ) +
  
  labs(x = "Date", y = "Weekly % visits", color = "Legend") +
  scale_color_manual(values = colors) +
  facet_wrap(~panel_label, ncol = 5, nrow = 2) +
  theme(legend.position = "bottom") +
  scale_x_date(
    date_breaks = "12 months",
    date_labels = "%b %Y"
  )

g0
