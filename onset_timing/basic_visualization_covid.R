# visualize time-series

library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_family = "Times", base_size = 14))
library(lubridate)
library(data.table)
library(purrr)
library(MMWRweek)
#library(cdcfluview)

###############################
# download ED  dataset from HHS 
###############################

setwd("C:/Work/research/nih_epi/models/breakpoint_extraction")
##download CDC dataset from API, 2020-present. Use the HSA-level time series version

# download.file("https://data.cdc.gov/api/views/rdmq-nq56/rows.csv?accessType=DOWNLOAD","ED.csv")

start_date <- "2022-10-14" # for covid/influenza
end_date <- "2026-03-01"

locs=read.csv("ED.csv") %>% dplyr::filter(!is.na(percent_visits_smoothed_covid))  %>% 
  distinct(geography, hsa, hsa_counties, hsa_nci_id)

locs.all=read.csv("ED.csv") %>% distinct(geography, hsa, hsa_counties, hsa_nci_id)

ed=read.csv("ED.csv") %>% distinct(geography, hsa, hsa_counties, hsa_nci_id, week_end, .keep_all = TRUE) %>%
  dplyr::select(-dplyr::any_of(c("county", "fips"))) %>%
  right_join(locs, by=c("geography", "hsa", "hsa_counties", "hsa_nci_id")) %>%
  mutate(week_end=as.Date(week_end))%>%
  arrange(hsa_nci_id, week_end)%>%
  filter(week_end>ymd(start_date) & week_end<ymd(end_date))

peaks <- read.csv("peaks_covid.csv", stringsAsFactors = FALSE) %>% # covid
  mutate(
    hsa_nci_id = as.numeric(hsa_nci_id),
    peak_date = as.Date(peak_date, format = "%m/%d/%Y")
  ) %>%
  filter(!is.na(hsa_nci_id), !is.na(peak_date))

breakpoints <- read.csv("breakpoints_covid.csv") %>% # covid
  filter(!is.na(breakpoint)) %>%
  mutate(breakpoint = as.Date(breakpoint, format = "%m/%d/%Y"),hsa_nci_id=as.numeric(hsa_nci_id))


#-------------


colors <- c("COVID" = "red","Flu" = "blue", 
            "RSV"="green")

# selected_ids <- c(814, 806, 662, 525,813,531,675,797,805,626) # bottom10 in population
# selected_ids <- c(287, 408, 453, 22,274,153,688,83,66,16) # top10 in population
# selected_ids <- c(378,573,294,773,343,362,364,389,46,92) # median10 in population

# selected_ids <- c(123,55,45,105,149,215,217,223,226,229) # top10 in breakpoint frequency
# selected_ids <- c(943,176,178,290,359,559,580,611,616) # bottom10 in breakpoint frequency
selected_ids <- c(474,477,478,481,489,491,495,498,501,502) # median10 in breakpoint frequency

vlines_df <- peaks %>% 
  mutate(hsa_nci_id = as.numeric(hsa_nci_id)) %>%
  filter(hsa_nci_id %in% selected_ids) %>%
  left_join(
    ed %>% 
      mutate(hsa_nci_id = as.numeric(hsa_nci_id)) %>%
      distinct(hsa_nci_id, hsa),
    by = "hsa_nci_id"
  ) 

vlines_bp <- breakpoints %>%
  filter(hsa_nci_id %in% selected_ids) %>%
  left_join(
    ed %>% 
      mutate(hsa_nci_id = as.numeric(hsa_nci_id)) %>%
      distinct(hsa_nci_id, hsa),
    by = "hsa_nci_id"
  )

g0 <- ggplot(subset(ed, hsa_nci_id %in% selected_ids)) + 

  geom_line(aes(week_end, percent_visits_smoothed_covid, color="COVID"), lwd=0.8) +
  # geom_line(aes(week_end, percent_visits_influenza, color="Flu"), lwd=0.8) +
  # geom_line(aes(week_end, percent_visits_rsv, color="RSV"), lwd=0.8) +

  geom_vline(
    data = vlines_df %>% 
      filter(peak_date >= as.Date(start_date) & 
               peak_date <= as.Date(end_date)),
    aes(xintercept = peak_date),
    linetype = "dotted",
    color = "gray",
    linewidth = 0.8
  )+
  geom_vline(
    data = vlines_bp %>% filter(breakpoint >= as.Date(start_date) &
                                  breakpoint <= as.Date(end_date)),
    aes(xintercept = breakpoint),
    linetype = "dashed",
    color = "orange",
    linewidth = 0.8
  ) +
  
  labs(x = "Date",
       y = "Weekly % visits",
       color = "Legend")+
  scale_color_manual(values = colors) +
  facet_wrap(~paste0("(", hsa_nci_id, ") ",hsa),ncol=5,nrow=2)+
  theme(legend.position = "bottom")+
  scale_x_date(
    date_breaks = "12 months",
    date_labels = "%b %Y"
  )

g0