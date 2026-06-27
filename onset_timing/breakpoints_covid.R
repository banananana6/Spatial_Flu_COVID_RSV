# find breakpoints through segmented regression, for covid

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
#library(cdcfluview)

setwd("") # add directory

##download CDC dataset from API, 2020-present. Use the HSA-level time series version
# download.file("https://data.cdc.gov/api/views/rdmq-nq56/rows.csv?accessType=DOWNLOAD","ED.csv")

locs=read.csv("ED.csv") %>% dplyr::filter(!is.na(percent_visits_covid))  %>% # covid
  distinct(geography, hsa, hsa_counties, hsa_nci_id)

locs.all=read.csv("ED.csv") %>% distinct(geography, hsa, hsa_counties, hsa_nci_id)

ed=read.csv("ED.csv") %>% distinct(geography, hsa, hsa_counties, hsa_nci_id, week_end, .keep_all = TRUE) %>%
  dplyr::select(-county, -fips) %>%
  right_join(locs, by=c("geography", "hsa", "hsa_counties", "hsa_nci_id")) %>%
  mutate(week_end=as.Date(week_end))%>%
  arrange(hsa_nci_id, week_end)

ed <- ed %>%
  group_by(hsa_nci_id) %>%
  arrange(week_end)

start_date <- "2022-10-22"
end_date <- "2026-05-01"

all_peak_dates <- read.csv("peaks_covid.csv", stringsAsFactors = FALSE) %>%
  filter(!hsa_nci_id %in% c("All", ""), !is.na(hsa_nci_id)) %>%
  mutate(
    hsa_nci_id = as.numeric(hsa_nci_id),
    peak_date = as.Date(peak_date, format = "%m/%d/%Y")
  ) %>%
  filter(!is.na(hsa_nci_id), !is.na(peak_date)) %>%
  arrange(hsa_nci_id, peak_date) %>%
  group_by(hsa_nci_id) %>%
  mutate(
    peak_num       = row_number(),
    prev_peak_date = lag(peak_date),
    window_start   = if_else(
      is.na(prev_peak_date),
      as.Date(start_date),
      prev_peak_date
    )
  ) %>%
  ungroup()

col2 <- c("black", "red")

#--- function to find breakpoints ---
find.breakpoints <- function(df,b_start_date,b_end_date) {
  df %>%
    # subset(pop > 50000) %>%
    mutate(
      days.from.start = as.numeric(week_end - ymd(b_start_date))
    ) %>%
    group_by(hsa_nci_id)%>%
    nest() %>%
    mutate(
      # start.date = map(data,
      #                  ~ .$days.from.start[15])[[1]],
      # start.date = map(data, ~ .$days.from.start[max(1, floor(nrow(.)/2))]),
      start.date = map(data, ~ median(.$days.from.start, na.rm = TRUE)),
      segmented.analysis = map2(data, start.date,
                                # ~ segmented(lm(percent_visits_influenza ~ days.from.start, data = .), # influenza
                                ~ segmented(lm(percent_visits_smoothed_covid ~ days.from.start, data = .), # covid
                                            seg.Z = ~ days.from.start,
                                            alpha = 0.0001,
                                            psi = list(days.from.start = start.date))),
      breakpoint = map(segmented.analysis,
                       ~ ymd(b_start_date) + days(round(.$psi[2], 0)))[[1]],

      loglikelihood = map(segmented.analysis,
                          ~ logLik(.)),

      segmented.ci = map(segmented.analysis,
                         ~ confint.segmented(.,
                                             method = "delta",
                                             level = 0.95)),
      early.ci = map(segmented.ci,
                     ~ ymd(b_start_date) + days(round(.[2], 0)))[[1]],

      late.ci = map(segmented.ci,
                    ~ ymd(b_start_date) + days(round(.[3], 0)))[[1]],

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
                                   # y = percent_visits_influenza, # influenza
                                   y = percent_visits_smoothed_covid, # covid
                                   color = "Data"),
                               linewidth = 0.6) +
                     geom_line(data = .y,
                               aes(x = days.from.start,
                                   y = model.fit,
                                   color = "Segmented Model"),
                               linewidth = 0.6) +
                     labs(x = paste0("Days from ",b_start_date), y = "Weekly % visits", color = "Legend")+
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

epidemic.start <- all_peak_dates %>%
  # filter((hsa_nci_id %in% c(395))) %>%
  mutate(
    results = pmap(
      list(hsa_nci_id, peak_num, window_start, peak_date),
      function(id, pk_num, w_start, pk_date) {
        
        segment_df <- ed %>%
          filter(
            hsa_nci_id == id,
            week_end > w_start,
            week_end <= pk_date
          ) %>%
          mutate(
            peak_num     = pk_num,
            window_start = w_start,
            peak_date    = pk_date
          )
        
        if (nrow(segment_df) < 10) return(NULL)
        
        result <- find.breakpoints(segment_df, w_start, pk_date)
        
        # --- Breakpoint correction ---
        bp <- as.Date(unlist(result$breakpoint), origin = "1970-01-01")
        
        if (!is.na(bp)) {
          # Get the valley (minimum) between the two peaks
          valley_df <- ed %>%
            filter(
              hsa_nci_id == id,
              week_end > w_start,
              week_end <= pk_date
            )
          
          if (nrow(valley_df) > 0) {
            min_idx  <- which.min(valley_df$percent_visits_smoothed_covid)
            min_date <- valley_df$week_end[min_idx]
            
            # If trough is after breakpoint, breakpoint is on the declining
            # portion of the previous wave — re-fit from trough to current peak
            if (min_date > bp) {
              
              segment_df2 <- ed %>%
                filter(
                  hsa_nci_id == id,
                  week_end >= min_date,  # start from the trough
                  week_end <= pk_date    # end at the current peak
                ) %>%
                mutate(
                  peak_num = pk_num,
                  window_start = min_date,
                  peak_date = pk_date
                )
              
              if (nrow(segment_df2) >= 10) {
                result2 <- find.breakpoints(segment_df2, min_date, pk_date)
                new_bp<- as.Date(unlist(result2$breakpoint), origin = "1970-01-01")
                new_bp<- new_bp[!is.na(new_bp)][1]
                
                if (!is.na(new_bp)) {
                  result$breakpoint <- list(new_bp)
                  result$early.ci <- list(as.Date(unlist(result2$early.ci), origin = "1970-01-01")[1])
                  result$late.ci <- list(as.Date(unlist(result2$late.ci), origin = "1970-01-01")[1])
                }
              }
            }
          }
        }
        
        result
      }
    )
  ) %>%
  filter(!sapply(results, is.null)) %>%
  # Normalize breakpoint/ci columns to Date before unnesting to avoid type conflicts
  mutate(results = map(results, ~ {
    .x$breakpoint <- list(as.Date(unlist(.x$breakpoint), origin = "1970-01-01"))
    .x$early.ci <- list(as.Date(unlist(.x$early.ci), origin = "1970-01-01"))
    .x$late.ci <- list(as.Date(unlist(.x$late.ci), origin = "1970-01-01"))
    .x
  })) %>%
  dplyr::select(peak_num, results) %>%
  unnest(results)


epidemic.start %>%
  filter(hsa_nci_id %in% c(10)) %>%

  ungroup %>%
  pull(plots) %>%
  ggarrange(plotlist = ., nrow = 2, ncol = 2,
            common.legend = TRUE,
            legend = "bottom")

epidemic.start_flat <- epidemic.start %>%
  ungroup() %>%
  mutate(
    breakpoint = as.Date(unlist(breakpoint), origin = "1970-01-01"),
    early.ci = as.Date(unlist(early.ci), origin = "1970-01-01"),
    late.ci = as.Date(unlist(late.ci), origin = "1970-01-01")
  ) %>%
  arrange(hsa_nci_id, breakpoint)

time_data <- epidemic.start1.function(epidemic.start_flat, date1=start_date, date2=end_date)

write.csv(time_data, "breakpoints_covid.csv", row.names = FALSE)
