# aggregate hsa with county population data to compute hsa population, weighted lat/long
# census dataset: https://simplemaps.com/data/us-counties

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
library(readxl)
#library(cdcfluview)

setwd("") # add directory

hsa=read_excel("HSA.xls")
census=read.csv("CENSUS.csv")

hsa <- hsa |> rename(hsa_id=1, hsa_name=2,fips=4)
hsa <- hsa |> mutate(fips=as.integer(fips))
census <- census |> rename(fips=4, population=9, lat=7, lon=8)

hsa_full <- hsa |>
  left_join(census, by = "fips")

hsa_geo <- hsa_full |>
  group_by(hsa_id,hsa_name) |>
  summarise(
    total_population = sum(population, na.rm = TRUE),
    weighted_lat = sum(lat * population, na.rm = TRUE) / sum(population, na.rm = TRUE),
    weighted_lon = sum(lon * population, na.rm = TRUE) / sum(population, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(hsa_geo, "HSA_GEO.csv", row.names = FALSE)