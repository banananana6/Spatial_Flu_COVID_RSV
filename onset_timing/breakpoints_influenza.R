# find breakpoints through segmented regression, for influenza

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
#library(cdcfluview)

###############################
# download ED  dataset from HHS 
###############################

setwd("C:/Work/research/nih_epi/models/breakpoint_extraction")
##download CDC dataset from API, 2020-present. Use the HSA-level time series version

# download.file("https://data.cdc.gov/api/views/rdmq-nq56/rows.csv?accessType=DOWNLOAD","ED.csv")

locs=read.csv("ED.csv") %>% dplyr::filter(!is.na(percent_visits_influenza))  %>% # influenza
  distinct(geography, hsa, hsa_counties, hsa_nci_id)

locs.all=read.csv("ED.csv") %>% distinct(geography, hsa, hsa_counties, hsa_nci_id)

# start_date <- "2023-07-01" # for 2023-2024/influenza
# end_date <- "2024-05-01"
start_date <- "2024-08-01" # for 2024-2025/influenza
end_date <- "2025-05-01"
# start_date <- "2025-08-01" # for 2025-2026/influenza
# end_date <- "2026-05-01"

ed=read.csv("ED.csv") %>% distinct(geography, hsa, hsa_counties, hsa_nci_id, week_end, .keep_all = TRUE) %>%
  dplyr::select(-county, -fips) %>%
  right_join(locs, by=c("geography", "hsa", "hsa_counties", "hsa_nci_id")) %>%
  mutate(week_end=as.Date(week_end))%>%
  # filter(week_end >= as.Date(start_date) & week_end <= as.Date(end_date)) %>%
  arrange(hsa_nci_id, week_end)

ed <- ed %>%
  group_by(hsa_nci_id) %>%
  arrange(week_end) %>%
  mutate(cases_corrected_smoothed = rollmean(percent_visits_influenza, k = 3, fill = NA, align = "right")) # influenza
  # mutate(cases_corrected_smoothed = rollmean(percent_visits_covid, k = 3, fill = NA, align = "right")) # covid

#----------------------------

#--- function to find peak ---
peak_dates <- ed %>%
  filter(week_end>ymd(start_date) & week_end<ymd(end_date)) %>%
  group_by(hsa_nci_id) %>%
  filter(max(percent_visits_influenza, na.rm = TRUE) > 2) %>%
  slice_max(percent_visits_influenza, n=1, with_ties=FALSE) %>% # influenza
  # slice_max(percent_visits_covid, n=1, with_ties=FALSE) %>% # covid
  dplyr::select(hsa_nci_id, peak_date=week_end)

col2 <- c("black", "red")

#--- function to find breakpoints ---
find.breakpoints <- function(df) {
  df %>% 
    # subset(pop > 50000) %>%
    mutate(
      days.from.start = as.numeric(week_end - ymd(start_date))
    ) %>% 
    group_by(hsa_nci_id) %>% 
    nest() %>% 
    mutate(
      # start.date = map(data, 
      #                  ~ .$days.from.start[15])[[1]],
      start.date = map(data, ~ .$days.from.start[max(1, floor(nrow(.)/2))]),
      segmented.analysis = map2(data, start.date,
                                ~ segmented(lm(percent_visits_influenza ~ days.from.start, data = .), # influenza
                                # ~ segmented(lm(percent_visits_covid ~ days.from.start, data = .), # covid
                                            seg.Z = ~ days.from.start,
                                            alpha = 0.0001,
                                            psi = list(days.from.start = start.date))),
      breakpoint = map(segmented.analysis,
                       ~ ymd(start_date) + days(round(.$psi[2], 0)))[[1]], 
      
      loglikelihood = map(segmented.analysis, 
                          ~ logLik(.)), 
      
      segmented.ci = map(segmented.analysis, 
                         ~ confint.segmented(., 
                                             method = "delta", 
                                             level = 0.95)),
      early.ci = map(segmented.ci,
                     ~ ymd(start_date) + days(round(.[2], 0)))[[1]],
      
      late.ci = map(segmented.ci,
                    ~ ymd(start_date) + days(round(.[3], 0)))[[1]],
      
      initial.slope = map(segmented.analysis, 
                          ~ slope(.)[[1]][1] )[[1]],
      second.slope = map(segmented.analysis, 
                         ~ slope(.)[[1]][2])[[1]],
      
      # correcting the ones where the second slope is negative 
      breakpoint = case_when(as.numeric(second.slope) < 0 ~ NA_Date_,
                             as.numeric(second.slope) < as.numeric(initial.slope) ~ NA_Date_,
                             TRUE ~ breakpoint),
      
      my.fitted = map(segmented.analysis, 
                      ~ fitted(.)),
      my.model = map2(data, my.fitted, 
                      ~ data.frame(days.from.start = .x$days.from.start, 
                                   model.fit = .y)),
      
      plots = map2(data, my.model, 
                   ~ ggplot() + 
                     geom_line(data = .x, 
                                aes(x = days.from.start,
                                    y = percent_visits_influenza, # influenza
                                    # y = percent_visits_covid, # covid
                                    color = "Data"),
                               linewidth = 0.6) + 
                     geom_line(data = .y, 
                                aes(x = days.from.start,
                                    y = model.fit, 
                                    color = "Segmented Model"),
                                linewidth = 0.6) + 
                     labs(x = paste0("Days after ",start_date), y = "Weekly % visits", color = "Legend")+
                     theme_bw() + 
                     scale_color_manual(values = c(col2[1],
                                                   col2[2])) +
                     ggtitle(paste0("(",hsa_nci_id, ") ",.x$hsa)))
    ) #%>%
  # drop_na(breakpoint)
} 

#--- epidemic start function: dataset for start of an epidemic ---

epidemic.start1.function <- function(df, date1, date2){
  
  epidemic.start1 <- df %>% 
    dplyr::select(hsa_nci_id, breakpoint, early.ci, late.ci) %>% 
    mutate(
      time.to.outcome = ymd(breakpoint) - ymd(date1),
      early.tto = ymd(early.ci) - ymd(date1), 
      late.tto = ymd(late.ci) - ymd(date1),
      
      outcome = case_when(
        is.na(breakpoint) ~ 0, 
        !is.na(breakpoint) ~ 1
      ), 
      
      time.to.outcome = case_when(
        is.na(breakpoint) ~ as.numeric(ymd(date2) - ymd(date1)) + 1, 
        !is.na(breakpoint) ~ as.numeric(time.to.outcome))
      
    )
  
    # left_join(geo2, by = "hsa_nci_id") #%>% 
  # subset(pop > 50000)
}

#----------------------

# results <- find.breakpoints(ed)
# time_data <- epidemic.start1.function(results, date1=start_date, date2=end_date)
# str(time_data)

epidemic.start <- 
  find.breakpoints(ed %>% 
      # filter(hsa_nci_id >= 1 & hsa_nci_id <= 900)%>%
        
      # if peak is lower than 2%, then the id is removed
      
      # subset(!(hsa_nci_id %in% c(197,250,329,490,595,804,806,942,950,965,977))) %>%  # for 2023-2024/influenza, 652 has double outbreaks & first onset
      # subset(!(hsa_nci_id %in% c(176,178,197,230,236,250,253,329,595,804,942,950,977,580,682))) %>% # for 2024-2025/influenza, 230 mostly flat
      # subset(!(hsa_nci_id %in% c(176,178,197,230,236,250,329,365,595,804))) %>% # for 2025-2026/influenza, 197 and 250 mostly flat
      
      subset(week_end > ymd(start_date) & week_end < ymd(end_date)) %>%
      left_join(peak_dates, by = "hsa_nci_id") %>%
      group_by(hsa_nci_id) %>%
      filter(week_end <= peak_date) %>% 
      ungroup()
) 

epidemic.start %>% 
  # filter(hsa_nci_id %in% c(652, 821, 890, 823, 936, 180, 538, 427, 929, 557)) %>% # 2023-2024 influenza (earliest 10 breakpoints)
  filter(hsa_nci_id %in% c(538, 505, 419, 530, 792, 743, 528, 131, 141,525)) %>% # 2024-2025 influenza
  # filter(hsa_nci_id %in% c(484, 717, 518, 756, 408, 934, 784, 482,389,490)) %>% # 2025-2026 influenza
  
  # filter(hsa_nci_id %in% c(1,10,11,654)) %>%
  
  ungroup %>% 
  pull(plots) %>% 
  ggarrange(plotlist = ., nrow = 2, ncol = 5, 
            common.legend = TRUE,
            legend = "bottom")

epidemic.start_flat <- epidemic.start %>%
  ungroup() %>%
  mutate(
    breakpoint = as.Date(unlist(breakpoint), origin = "1970-01-01"),
    early.ci   = as.Date(unlist(early.ci),   origin = "1970-01-01"),
    late.ci    = as.Date(unlist(late.ci),    origin = "1970-01-01")
  )

time_data <- epidemic.start1.function(epidemic.start_flat, date1 = start_date, date2 = end_date)

write.csv(time_data, "breakpoints_2025_2026_influenza.csv", row.names = FALSE)


