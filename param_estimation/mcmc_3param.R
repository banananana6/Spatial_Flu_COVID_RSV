# parameter estimation with metropolis-hastings MCMC

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
    # mutate(mean.value = map(data, ~ mean(.$value))[[1]]) %>%
    mutate(mean.value = map_dbl(data, ~ mean(.$value)))
    dplyr::select(-data) %>%
    ungroup() %>%
    pivot_wider(names_from = parameter, values_from = mean.value) %>%
    dplyr::select(variable, mean, sd)
}

for (season_label in names(season_files)) {
  
  cat("\n===== Season:", season_label, "=====\n")
  out_file <- paste0("epimcmc_results_", gsub("-", "_", season_label), ".csv")
  wrote_header <- FALSE
  
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
                      
                      sus.par.ini  = 0.0012, beta.ini = 1.2, spark.ini = 0.001,
                      pro.sus.var  = 0.0006, pro.beta.var = 0.4, pro.spark.var = 0.0002,
                      prior.sus.dist  = "uniform", prior.sus.par  = c(0, 0.3),
                      prior.beta.dist = "uniform", prior.beta.par = c(0.5, 2.5),
                      prior.spark.dist = "uniform", prior.spark.par = c(0, 0.01),

                      adapt = TRUE, acc.rate = 0.4)
    
    postburnin1 <- model1$Estimates[(n.burnin + 1):n.iterations, ]
    
    sus.par <- mean(postburnin1[, 1])
    beta.par<- mean(postburnin1[, 2])
    spark.par <- mean(postburnin1[, 3])
    
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

                      sus.par.ini  = 0.0003, beta.ini = 1.2, spark.ini = 0.001,
                      pro.sus.var  = 0.0001, pro.beta.var = 0.6, pro.spark.var = 0.00002,
                      prior.sus.dist  = "uniform", prior.sus.par  = c(0, 0.1),
                      prior.beta.dist = "uniform", prior.beta.par = c(0.5, 2.5),
                      prior.spark.dist = "uniform", prior.spark.par = c(0, 0.01),

                      adapt = TRUE, acc.rate = 0.4)

    postburnin1 <- model1$Estimates[(n.burnin + 1):n.iterations, ]
    
    sus.par <- mean(postburnin1[, 1])
    beta.par<- mean(postburnin1[, 2])
    spark.par <- mean(postburnin1[, 3])

    loglike2 <- epilike(SI.model, tmax = t_end, Sformula = ~ -1 + log10pop,
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

posterior_mean_beta <- mean(model1$Estimates[n.burnin:n.iterations, 2])
posterior_mean_spark <- mean(model1$Estimates[n.burnin:n.iterations, 3])

ll_profile_alpha <- sapply(alpha_grid, function(a) {
  epilike(object = SI.model, tmin = t_start, tmax = t_end, 
          sus.par = a, beta = posterior_mean_beta, spark = posterior_mean_spark)
})

plot(
  alpha_grid, ll_profile_alpha,
  type = "l",
  xlab = "alpha",
  ylab = "log-likelihood"
)

# --- mark alpha value in grid ---
alpha_mark <- 0.00031817247
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

## beta profile likelihood, model1
beta_grid <- seq(0.1, 3, by = 0.1)

posterior_mean_sus <- mean(model1$Estimates[n.burnin:n.iterations, 1])
posterior_mean_spark<- mean(model1$Estimates[n.burnin:n.iterations, 3])

ll_profile <- sapply(beta_grid, function(b) {
  epilike(object = SI.model, tmin = t_start, tmax = t_end, 
          sus.par = posterior_mean_sus, beta = b, spark = posterior_mean_spark)
})

plot(
  beta_grid, ll_profile,
  type = "l",
  xlab = "beta",
  ylab = "log-likelihood"
)

# --- mark beta value in grid ---
beta_mark <-1.301

idx <- which.min(abs(beta_grid - beta_mark))

x_mark <- beta_grid[idx]
y_mark <- ll_profile[idx]

abline(v = x_mark, lty = 2, col = "red")
abline(h = y_mark, lty = 2, col = "blue")

points(x_mark, y_mark, pch = 19, col = "black", cex = 1.3)

text(
  x_mark, y_mark,
  labels = sprintf("(%.1f, %.2f)", x_mark, y_mark),
  pos = 4,
  cex = 0.8
)

#---------------------

## spark profile likelihood, model1
spark_grid <- seq(0.00001, 0.0001, by = 0.000002)

posterior_mean_alpha  <- mean(model1$Estimates[n.burnin:n.iterations, 1])
posterior_mean_beta <- mean(model1$Estimates[n.burnin:n.iterations, 2])

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
spark_mark <- 7.40766101656709e-05

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

## for model 2 profile likelihoods refer to the mcmc_fixedbeta.R file 

# #--- save params ---
# samples <- model1$Estimates[(n.burnin+1):nrow(model1$Estimates), ]
# 
# for (i in 1:ncol(samples)) {
#   
#   param_name <- colnames(samples)[i]
#   
#   write.csv(
#     data.frame(value = samples[, i]),
#     file = paste0("MCMC_parameters/", param_name, "_chain_postburnin.csv"),
#     row.names = FALSE
#   )
# }
