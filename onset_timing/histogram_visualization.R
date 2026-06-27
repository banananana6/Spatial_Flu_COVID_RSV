# visualization of breakpoints with histogram

library(ggplot2)
library(dplyr)

setwd("C:/Work/research/nih_epi/models/breakpoint_extraction")

# start_date <-"2025-09-01"
# end_date <- "2026-03-01"

breakpoints <- read.csv("breakpoints_2024_2025_rsv.csv")
breakpoints <- breakpoints %>%
  filter(!is.na(breakpoint)) %>%
  mutate(breakpoint = as.Date(breakpoint, format = "%m/%d/%Y"))
  # filter(breakpoint >= as.Date(start_date) & breakpoint <= as.Date(end_date))

breakpoints <- breakpoints[order(breakpoints$breakpoint), ]

#--- to find earliest breakpoints ---
# breakpoints <- breakpoints[1:10,]
# breakpoints$hsa_nci_id

ggplot(breakpoints, aes(x = breakpoint)) +
  geom_histogram(binwidth = 7, fill = "steelblue", color = "white") +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  labs(x = "Breakpoint date", y = "Number of HSAs") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

breakpoints$breakpoint <- format(breakpoints$breakpoint, "%m/%d/%Y")
write.csv(breakpoints, "breakpoints_covid.csv", row.names = FALSE)