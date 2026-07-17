# MCMC for the two-parameter model (alpha, spark) with fixed beta

library(doRNG)
library(foreach)
library(doParallel)
library(dplyr)
library(tidyverse)
library(purrr)
library(EpiILM)

start_time <- Sys.time()

setwd("") # add directory

cl <- makeCluster(parallel::detectCores() - 1)
registerDoParallel(cl)

hsa_geo <- read.csv("HSA_GEO.csv")

# --- seasons ---
season_files <- list(
  "2023-2024_influenza" = "breakpoints_2023_2024_influenza.csv"
  # "2024-2025_influenza" = "breakpoints_2024_2025_influenza.csv"
  # "2025-2026_influenza" = "breakpoints_2025_2026_influenza.csv",
  # "2023-2024_rsv" = "breakpoints_2023_2024_rsv.csv",
  # "2024-2025_rsv" = "breakpoints_2024_2025_rsv.csv",
  # "2025-2026_rsv" = "breakpoints_2025_2026_rsv.csv"
)

n.iterations <- 5000
n.burnin <- 2000
n.samples <- 1
fixed_b <- 0.5

# --- makes parameter table ---
table.function <- function(df) {
  df %>%
    as_tibble() %>%
    mutate(
      rownumber = row_number(),
      parameter = case_when(
        rownumber %% 2 == 1 ~ "mean",
        rownumber %% 2 == 0 ~ "sd",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select(-rownumber) %>%
    pivot_longer(cols = alpha.1:spark,
                 names_to  = "variable",
                 values_to = "value") %>%
    group_by(variable, parameter) %>%
    arrange(variable, parameter) %>%
    nest() %>%
    mutate(mean.value = map(data, ~ mean(.$value))[[1]]) %>%
    dplyr::select(-data) %>%
    ungroup() %>%
    pivot_wider(names_from = parameter, values_from = mean.value) %>%
    dplyr::select(variable, mean, sd)
}

for (season_label in names(season_files)) {
  
  cat("\n===== Season:", season_label, "=====\n")
  out_file <- paste0("epimcmc_results_", gsub("-", "_", season_label), ".csv")
  wrote_header <-FALSE
  
  df_raw <- read.csv(season_files[[season_label]]) %>%
    mutate(
      hsa_nci_id = as.numeric(hsa_nci_id),
      time.to.outcome = as.numeric(time.to.outcome),
      early.tto = as.numeric(early.tto),
      late.tto = as.numeric(late.tto)
    ) %>%
    left_join(hsa_geo, by="hsa_nci_id")
  
  df_raw <- df_raw %>%
    filter(!(hsa_nci_id %in% c(821, 822, 823)))
  # filter(weighted_lon>-90) # limits to east coast
  
  for (i in 1:n.samples) {
    df.sample <- df_raw %>%
      rowwise() %>%
      mutate(
        time.to.outcome = case_when(
          !is.na(breakpoint) ~ round(runif(1, early.tto, late.tto), 0),
          TRUE ~ time.to.outcome
        )
      ) %>%
      ungroup() %>%
      filter(!is.na(hsa_nci_id), !is.na(time.to.outcome), !is.na(total_population))
    
    log10pop <- log10(as.numeric(df.sample$total_population))
    
    SI.model <- EpiILM::as.epidata(type = "SI", n = nrow(df.sample), x = as.numeric(df.sample$weighted_lon),
                                   y = as.numeric(df.sample$weighted_lat), inftime = as.numeric(df.sample$time.to.outcome))
    
    t_start <- 15
    t_end <- max(as.numeric(df.sample$time.to.outcome))
    
    # --- Model 1: no susceptibility covariate ---
    model1 <- epimcmc(SI.model,
                      tmin = t_start, tmax = t_end,
                      niter = n.iterations,
                      
                      sus.par.ini  = 0.1, beta.ini = fixed_b, spark.ini = 0.001,
                      pro.sus.var  = 0.01, pro.beta.var = 0, pro.spark.var = 0.00005,
                      prior.sus.dist  = "uniform", prior.sus.par  = c(0, 0.3),
                      prior.beta.dist = "uniform", prior.beta.par = c(0.5, 2.5),
                      prior.spark.dist = "uniform", prior.spark.par = c(0, 0.01),
                      
                      adapt = TRUE, acc.rate = 0.4)
    
    sus.par <- mean(model1$Estimates[(n.burnin+1):n.iterations, "alpha.1"])
    beta.par<- fixed_b
    spark.par <- mean(model1$Estimates[(n.burnin+1):n.iterations, "spark"])
    
    loglike1 <- epilike(SI.model, tmax = t_end,
                        sus.par = sus.par, beta = beta.par, spark = spark.par)
    dic1 <- epidic(burnin = n.burnin, niter = n.iterations,
                   LLchain = model1$Loglikelihood, LLpostmean = loglike1)
    
    parameters1 <- rbind(
      apply(model1$Estimates[n.burnin:n.iterations, ], 2, mean),
      apply(model1$Estimates[n.burnin:n.iterations, ], 2, sd)
    )
    
    # --- Model 2: log10(population) susceptibility covariate ---
    model2 <- epimcmc(SI.model,
                      tmin = t_start, tmax = t_end,
                      Sformula = ~ -1 + log10pop,
                      niter = n.iterations,
                      
                      sus.par.ini  = 0.0003, beta.ini = fixed_b, spark.ini = 0.0001,
                      pro.sus.var  = 0.00001, pro.beta.var = 0, pro.spark.var = 0.00001,
                      prior.sus.dist  = "uniform", prior.sus.par  = c(0, 0.1),
                      prior.beta.dist = "uniform", prior.beta.par = c(0.5, 2.5),
                      prior.spark.dist = "uniform", prior.spark.par = c(0, 0.01),
                      
                      adapt = TRUE, acc.rate = 0.4)
    
    sus.par <- mean(model2$Estimates[(n.burnin+1):n.iterations, "alpha.1"])
    beta.par<- fixed_b
    spark.par <- mean(model2$Estimates[(n.burnin+1):n.iterations, "spark"])
    
    loglike2 <- epilike(SI.model, tmax = t_end, Sformula = ~ -1+log10pop,
                        sus.par = sus.par, beta = beta.par, spark = spark.par)
    dic2 <- epidic(burnin = n.burnin, niter = n.iterations,
                   LLchain = model2$Loglikelihood, LLpostmean = loglike2)
    
    parameters2 <- rbind(
      apply(model2$Estimates[n.burnin:n.iterations, ], 2, mean),
      apply(model2$Estimates[n.burnin:n.iterations, ], 2, sd)
    )
    
    # --- results ---
    sample_results <- table.function(parameters1) %>%
      rbind(c("variable" = "DIC", "mean" = dic1, "sd" = NA)) %>%
      mutate(model = "model1", .before = "variable") %>%
      rbind(
        table.function(parameters2) %>%
          rbind(c("variable" = "DIC", "mean" = dic2, "sd" = NA)) %>%
          mutate(model = "model2", .before = "variable")
      ) %>%
      mutate(
        wave   = season_label,
        sample = i,
        .before = "model"
      )
    
    #--- append to csv ---
    write.table(sample_results, file = out_file, append = wrote_header,
                sep = ",", row.names = FALSE, col.names = !wrote_header)
    wrote_header <- TRUE
    
  }
  cat("  Season", season_label, "complete", out_file, "\n")
}

stopCluster(cl)
end_time <-Sys.time()
print(end_time-start_time)


## --- various plots below ----

# for trace plot
plot(model1, partype = "parameter", start = 2001, density = FALSE)

# ll plot
ll_trace <- model1$Loglikelihood
plot(ll_trace, type = "l", xlab = "iteration", ylab = "log-likelihood")

## alpha profile likelihood, model1
alpha_grid <- seq(0.0001, 0.002, by = 0.00002)

posterior_mean_beta <- fixed_b
posterior_mean_spark <- mean(model1$Estimates[(n.burnin+1):n.iterations, "spark"])

ll_profile_alpha <- sapply(alpha_grid, function(a) {
  epilike(object = SI.model, tmin = t_start, tmax = t_end, 
          sus.par = a, beta = fixed_b, spark = posterior_mean_spark)
})

plot(
  alpha_grid, ll_profile_alpha,
  type = "l",
  xlab = "alpha",
  ylab = "log-likelihood"
)

# --- mark alpha value in grid ---
alpha_mark <- 0.0003250
idx <- which.min(abs(alpha_grid - alpha_mark))

x_mark <- alpha_grid[idx]
y_mark <- ll_profile_alpha[idx]

plot(
  alpha_grid, ll_profile_alpha,
  type = "l",
  xlab = "alpha",
  ylab = "log-likelihood"
)

abline(v = x_mark, lty = 2, col = "red")
abline(h = y_mark, lty = 2, col = "blue")

points(x_mark, y_mark, pch = 19, col = "black", cex = 1.3)

text(
  x_mark, y_mark,
  labels = sprintf("(%.6f, %.2f)", x_mark, y_mark),
  pos = 4,
  cex = 0.8
)

# --------------

## spark profile likelihood, model1
spark_grid <- seq(0.00001, 0.0005, by = 0.000002)

posterior_mean_alpha <- mean(model1$Estimates[(n.burnin+1):n.iterations, "alpha.1"])
posterior_mean_beta <- fixed_b

ll_profile_spark <- sapply(spark_grid, function(a) {
  epilike(
    object = SI.model,
    tmin = t_start,
    tmax = t_end,
    sus.par = posterior_mean_alpha,
    beta = posterior_mean_beta,
    spark=a
    
  )
})

plot(
  spark_grid, ll_profile_spark,
  type = "l",
  xlab = "spark",
  ylab = "log-likelihood"
)

#---
spark_mark <- 7.835e-5

idx <- which.min(abs(spark_grid - spark_mark))

x_mark <- spark_grid[idx]
y_mark <- ll_profile_spark[idx]

abline(v = x_mark, lty = 2, col = "red")
abline(h = y_mark, lty = 2, col = "blue")

points(x_mark, y_mark, pch = 19, col = "black", cex = 1.3)

text(
  x_mark, y_mark,
  labels = sprintf("(%.6f, %.2f)", x_mark, y_mark),
  pos = 4,
  cex = 0.8
)

#=======================

## alpha profile likelihood, model2
alpha_grid_m2 <- seq(0.00001, 0.0002, by = 0.000002)

posterior_mean_beta_m2 <- fixed_b
posterior_mean_spark_m2 <- mean(model2$Estimates[(n.burnin+1):n.iterations, "spark"])

ll_profile_alpha_m2 <- sapply(alpha_grid_m2, function(a) {
  epilike(object = SI.model, tmin = t_start, tmax = t_end, Sformula = ~ -1 + log10pop,
          sus.par = a, beta = posterior_mean_beta_m2, spark = posterior_mean_spark_m2)
})

plot(
  alpha_grid_m2, ll_profile_alpha_m2,
  type = "l",
  xlab = "alpha",
  ylab = "log-likelihood"
)

# --- mark alpha in grid ---
alpha_mark <- 6.534e-05
idx_m2 <- which.min(abs(alpha_grid_m2 - alpha_mark))

x_mark_m2 <- alpha_grid_m2[idx_m2]
y_mark_m2 <- ll_profile_alpha_m2[idx_m2]

abline(v = x_mark_m2, lty = 2, col = "red")
abline(h = y_mark_m2, lty = 2, col = "blue")

points(x_mark_m2, y_mark_m2, pch = 19, col = "black", cex = 1.3)

text(
  x_mark_m2, y_mark_m2,
  labels = sprintf("(%.6f, %.2f)", x_mark_m2, y_mark_m2),
  pos = 4,
  cex = 0.8
)

#---------------------

## spark profile likelihood, model2
spark_grid_m2 <- seq(0.00001, 0.0001, by = 0.000002)

posterior_mean_alpha_m2 <- mean(model2$Estimates[(n.burnin+1):n.iterations, "alpha.1"])
posterior_mean_beta_m2  <- fixed_b

ll_profile_spark_m2 <- sapply(spark_grid_m2, function(a) {
  epilike(object = SI.model, tmin = t_start, tmax = t_end,
    Sformula = ~ -1+log10pop, sus.par = posterior_mean_alpha_m2,
    beta = posterior_mean_beta_m2, spark = a)
})

plot(
  spark_grid_m2, ll_profile_spark_m2,
  type = "l",
  xlab = "spark",
  ylab = "log-likelihood"
)

# -- mark on grid --- 
spark_mark <- 7.8357e-5
idx_m2 <- which.min(abs(spark_grid_m2-spark_mark))

x_mark_m2 <- spark_grid_m2[idx_m2]
y_mark_m2 <- ll_profile_spark_m2[idx_m2]

abline(v = x_mark_m2, lty = 2, col = "red")
abline(h = y_mark_m2, lty = 2, col = "blue")

points(x_mark_m2, y_mark_m2, pch = 19, col = "black", cex = 1.3)

text(
  x_mark_m2, y_mark_m2,
  labels = sprintf("(%.6f, %.2f)", x_mark_m2, y_mark_m2),
  pos = 4,
  cex = 0.8
)

