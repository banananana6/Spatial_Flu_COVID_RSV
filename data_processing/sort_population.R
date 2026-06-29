# find lowest, highest, median populations

library(dplyr)

setwd("") # add directory

breakpoints <- read.csv("breakpoints_covid.csv") %>%
  filter(!is.na(breakpoint)) %>%
  mutate(breakpoint = as.Date(breakpoint, format = "%m/%d/%Y"),hsa_nci_id=as.numeric(hsa_nci_id))

pop <- read.csv("HSA_GEO.csv") %>%
  mutate(total_population = as.numeric(gsub(",", "", total_population))) %>%
  filter(hsa_id %in% breakpoints$hsa_nci_id)

bottom_10 <- pop %>% slice_min(total_population, n=10)
top_10 <- pop %>% slice_max(total_population, n=10)
med <- median(pop$total_population, na.rm = TRUE)
median_10 <- pop %>%
  mutate(dist_to_median = abs(total_population-med)) %>%
  slice_min(dist_to_median, n=10, with_ties=FALSE) %>%
  dplyr::select(-dist_to_median)


print(bottom_10)
print(top_10)
print(median_10)