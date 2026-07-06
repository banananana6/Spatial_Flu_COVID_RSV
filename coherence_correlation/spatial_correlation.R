# spatial coherence test

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

start_date <- "2023-07-01" # for 2023-2024/influenza
end_date <- "2024-05-01"
# start_date <- "2024-08-01" # for 2024-2025/influenza
# end_date <- "2025-05-01"
# start_date <- "2025-08-01" # for 2025-2026/influenza
# end_date <- "2026-05-01"
# 
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
# 
# start_date <- "2023-07-01" # for 2023-2024/rsv
# end_date <- "2024-05-01"
# start_date <- "2024-07-01" # for 2024-2025/rsv
# end_date <- "2025-05-01"
# start_date <- "2025-07-01" # for 2025-2026/rsv
# end_date <- "2026-05-01"

locs <- read.csv("ED.csv") %>%
  dplyr::filter(!is.na(percent_visits_influenza)) %>%
  distinct(geography, hsa, hsa_counties, hsa_nci_id)

ed <- read.csv("ED.csv") %>%
  distinct(geography, hsa, hsa_counties, hsa_nci_id, week_end, .keep_all = TRUE) %>%
  dplyr::select(-dplyr::any_of(c("county", "fips"))) %>%
  right_join(locs, by = c("geography", "hsa", "hsa_counties", "hsa_nci_id")) %>%
  mutate(week_end = as.Date(week_end)) %>%
  mutate(hsa_nci_id=as.numeric(hsa_nci_id))%>%
  arrange(hsa_nci_id, week_end)

hsa_geo <- read.csv("HSA_GEO.csv") %>%
  filter(hsa_id %in% locs$hsa_nci_id) %>%
  filter(!(hsa_id %in% c(821,822,823))) %>% # filters out hawaii
  filter(total_population>5000) %>%
  mutate(hsa_id=as.numeric(hsa_id)) %>%
  arrange(hsa_id)

ed <- ed %>%
  filter(hsa_nci_id %in% hsa_geo$hsa_id) %>%
  filter(!(hsa_nci_id %in% c(821,822,823))) %>% # filters out hawaii
  filter(week_end>=as.Date(start_date) & week_end<=as.Date(end_date)) %>%
  arrange(hsa_nci_id) 

lat <- hsa_geo$weighted_lat
long <- hsa_geo$weighted_lon

z_mat <- acast(ed, hsa_nci_id ~ week_end, value.var = "percent_visits_influenza")

ncf_test <- function(long, lat, df) {
  test <- Sncf(long, lat, df, type= "boot", resamp = 1000, latlon = TRUE, na.rm  = TRUE)
  
  df.test <- tibble(
    cbar = test$real$cbar,
    x = test$real$predicted$x[1, ],
    y = test$real$predicted$y[1, ]
  )
  
  return(df.test)  
}

result<- ncf_test(long,lat,z_mat)

ggplot(result, aes(x=x,y=y))+
  geom_line()+
  geom_hline(yintercept=0,linetype="dashed",color="red")+
  geom_hline(yintercept=unique(result$cbar),
             linetype="dashed",color="blue")+
  labs(
    x="Distance (km)",
    y="Correlation",
    title="Spatial Coherence"
  )
