# correlation test for onset date vs. peak size

library(ggplot2)
library(dplyr)
library(tidyr)

setwd("C:/Work/research/nih_epi/models/breakpoint_extraction")

# start_date <- "2023-07-01" # for 2023-2024/influenza
# end_date <- "2024-05-01"
# start_date <- "2024-08-01" # for 2024-2025/influenza
# end_date <- "2025-05-01"
# start_date <- "2025-08-01" # for 2025-2026/influenza
# end_date <- "2026-05-01"

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

start_date <- "2023-07-01" # for 2023-2024/rsv
end_date <- "2024-05-01"
# start_date <- "2024-07-01" # for 2024-2025/rsv
# end_date <- "2025-05-01"
# start_date <- "2025-07-01" # for 2025-2026/rsv
# end_date <- "2026-05-01"

breakpoints <- read.csv("breakpoints_rsv.csv") #change
breakpoints <-breakpoints%>%
  filter(!is.na(breakpoint)) %>%
  mutate(breakpoint=as.Date(breakpoint,format="%m/%d/%Y")) %>%
  filter(breakpoint>=as.Date(start_date) & breakpoint<=as.Date(end_date)) %>%
  mutate(hsa_nci_id=as.integer(hsa_nci_id)) %>%
  arrange(hsa_nci_id,breakpoint)

peaks<- read.csv("peaks_rsv.csv") # change
peaks <- peaks %>%
  filter(!is.na(peak_date))%>%
  mutate(peak_date=as.Date(peak_date,format="%m/%d/%Y")) %>%
  filter(peak_date>=as.Date(start_date) & peak_date<=as.Date(end_date)) %>%
  mutate(hsa_nci_id=as.integer(hsa_nci_id))%>%
  arrange(hsa_nci_id,peak_date) 

hsa_geo <- read.csv("HSA_GEO.csv")
hsa_geo <- hsa_geo %>%
  # filter(hsa_id %in% peaks$hsa_nci_id) %>%
  mutate(hsa_id=as.integer(hsa_id)) %>%
  arrange(hsa_id)

# --- peak vs onset ---

find_next_peak <- function(breakpoint, peaks_df,hsa_id) {
  peaks_df %>%
    filter(hsa_nci_id==hsa_id,peak_date>as.Date(breakpoint))%>%
    arrange(peak_date) %>%
    slice(1) %>%
    dplyr::select(peak_date,peak_size)
}

peak_table <- breakpoints %>%
  rowwise() %>%
  mutate(result=list(find_next_peak(breakpoint, peaks,hsa_nci_id))) %>%
  unnest(result) %>%
  ungroup()

# --- weighted lat/lon vs onset ---

peak_table <-breakpoints %>%
  left_join(hsa_geo %>% dplyr::select(hsa_id, weighted_lat, weighted_lon),
            by = c("hsa_nci_id" = "hsa_id"))

# --- population vs onset ---

peak_table <-breakpoints %>%
  left_join(hsa_geo %>% dplyr::select(hsa_id, total_population),
            by = c("hsa_nci_id" = "hsa_id")) %>%
  filter(!is.na(total_population))

# --- peak size vs onset date correlation/plot ---

cor(
  as.numeric(as.Date(peak_table$breakpoint)-as.Date(start_date)),
  peak_table$peak_size,
  use="complete.obs"
)

ggplot(peak_table, aes(x = peak_table$breakpoint, y =peak_table$peak_size)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE) +
  labs(
    x = "Onset date",
    y = "Peak size",
    title = "Peak size vs onset date"
  ) +
  theme_bw()

# --- lat vs onset date correlation/plot ---

cor(
  as.numeric(as.Date(peak_table$breakpoint)-as.Date(start_date)),
  peak_table$weighted_lat,
  use="complete.obs"
)
ggplot(peak_table, aes(x = peak_table$breakpoint, y =peak_table$weighted_lat)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE) +
  labs(
    x = "Onset date",
    y = "Latitude",
    title = "Latitude vs onset date"
  ) +
  theme_bw()

# --- lon vs. onset date correlation/plot ---

cor(
  as.numeric(as.Date(peak_table$breakpoint)-as.Date(start_date)),
  peak_table$weighted_lon,
  use="complete.obs"
)
ggplot(peak_table, aes(x = peak_table$breakpoint, y =peak_table$weighted_lon)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE) +
  labs(
    x = "Onset date",
    y = "Longitude",
    title = "Longitude vs onset date"
  ) +
  theme_bw()

# population vs. onset correlation/plot ---

cor(
  as.numeric(as.Date(peak_table$breakpoint)-as.Date(start_date)),
  log10(as.numeric(peak_table$total_population)),
  use="complete.obs"
)

ggplot(peak_table, aes(x = peak_table$breakpoint, y =log10(as.numeric(peak_table$total_population)))) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE) +
  labs(
    x = "Onset date",
    y = "log(population)",
    title = "log(population) vs onset date"
  ) +
  theme_bw()
