# simulation of epidemic with comparison to observed data

library(EpiILM)
library(dplyr)
library(ggplot2)

setwd("") #add directory

sim_start_date <- as.Date("2023-09-15")
# seed_end_date  <- sim_start_date + 6  

breakpoints <- read.csv("breakpoints_2023_2024_influenza.csv") %>%
  mutate(
    hsa_nci_id = as.numeric(hsa_nci_id),
    breakpoint = as.Date(breakpoint, format = "%m/%d/%Y")
  ) %>%
  filter(
    !is.na(breakpoint),
    !(hsa_nci_id %in% c(821, 822, 823))
  ) %>%
  arrange(breakpoint)

hsa_geo <- read.csv("HSA_GEO.csv") %>%
  mutate(
    hsa_nci_id = as.numeric(hsa_nci_id),
    total_population = as.numeric(total_population)
  ) %>%
  filter(
    !(hsa_nci_id %in% c(821, 822, 823)),
    !is.na(total_population),
    !is.na(weighted_lon),
    !is.na(weighted_lat)
  )

n <- nrow(hsa_geo)
coords <- cbind(hsa_geo$weighted_lon, hsa_geo$weighted_lat)

nsim <- 50
seed_time <- 1 # infected immediately at sim_start_date

set.seed(42)

# sim_counts_all <- data.frame()
# 
# seed_breakpoints <- breakpoints %>%
#   filter(breakpoint <= sim_start_date + 6)
# 
# init_inftime <- rep(0, n)
# 
# init_inftime[match(seed_breakpoints$hsa_nci_id, hsa_geo$hsa_nci_id)] <-
#   as.integer((seed_breakpoints$breakpoint - sim_start_date) / 7) + 1

init_inftime <- rep(0, n)
seed_idx <- sample(seq_len(n), 1)
init_inftime[seed_idx] <- 1  # one random seed

sim_counts_all <- data.frame()

# --- model1---

alpha <- 3.309518e-03
beta <- 0.5
spark <- 3.265915e-05

t_end <- 100

for (i in 1:nsim) {
  sim1 <- epidata(type = "SI", n = n, tmax = t_end, sus.par = alpha,
    beta = beta, spark = spark, x = coords[, 1], y = coords[, 2], inftime=init_inftime)
  
  sim_breakpoints <- data.frame(
    hsa_nci_id = hsa_geo$hsa_nci_id,
    inftime = sim1$inftime,
    breakpoint = sim_start_date+7*(sim1$inftime-1)
  ) %>%
    filter(inftime > 0) %>%
    semi_join(breakpoints, by = "hsa_nci_id")
  
  sim_counts <- sim_breakpoints %>%
    mutate(week = as.Date(cut(breakpoint, breaks = "week"))) %>%
    count(week, name = "count") %>%
    mutate(sim = i)
  
  sim_counts_all <- bind_rows(sim_counts_all, sim_counts)
}

obs_counts <- breakpoints %>%
  mutate(week = as.Date(cut(breakpoint, breaks = "week"))) %>%
  count(week, name = "count")

ggplot() +
  geom_line(
    data = sim_counts_all,
    aes(x = week, y = count, group = sim),
    color = "gray80"
  ) +
  geom_line(
    data = obs_counts,
    aes(x = week, y = count),
    color = "red",
    linewidth = 1.2
  ) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  labs(x = "Breakpoint date", y = "Number of HSAs") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# --- model 2 ---
set.seed(42)

alpha <- 1.824705e-02
beta<- 1
spark<- 3.469454e-05
t_end <- 100

log10pop <- log10(hsa_geo$total_population)

for (i in 1:nsim){
  sim2 <- epidata(type= "SI", n= n, tmax= t_end, sus.par=alpha, beta=beta, spark=spark, Sformula=~-1+log10pop,
                  x=coords[, 1], y=coords[, 2])
  
  sim_breakpoints <- data.frame(
    hsa_nci_id = hsa_geo$hsa_nci_id,
    inftime = sim2$inftime,
    breakpoint = as.Date(sim_start_date) + 7 * (sim2$inftime - 1)
  ) %>%
    filter(inftime > 0) %>%
    semi_join(breakpoints, by = "hsa_nci_id")
  
  sim_counts <- sim_breakpoints %>%
    mutate(week = as.Date(cut(breakpoint, breaks = "week"))) %>%
    count(week, name = "count") %>%
    mutate(sim=i)
  
  sim_counts_all <-bind_rows(sim_counts_all,sim_counts)
}

sim_breakpoints <- data.frame(
  hsa_nci_id = hsa_geo$hsa_nci_id,
  inftime = sim2$inftime,
  breakpoint = as.Date(sim_start_date) + 7*(sim2$inftime - 1)
) %>%
  filter(inftime > 0) %>%
  semi_join(breakpoints, by = "hsa_nci_id")

obs_counts <- breakpoints %>%
  mutate(week = as.Date(cut(breakpoint, breaks = "week"))) %>%
  count(week, name = "count")

sim_counts <- sim_breakpoints %>%
  mutate(week = as.Date(cut(breakpoint, breaks = "week"))) %>%
  count(week, name = "count")

ggplot() +
  geom_line(
    data = sim_counts_all,
    aes(x = week, y = count, group = sim),
    color = "gray80"
  ) +
  geom_line(
    data = obs_counts,
    aes(x = week, y = count),
    color = "red",
    linewidth = 1.2
  ) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "1 month") +
  labs(x = "Breakpoint date", y = "Number of HSAs") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

