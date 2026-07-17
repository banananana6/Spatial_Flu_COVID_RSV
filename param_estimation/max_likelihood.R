# MLE for the three-parameter model (alpha, beta, eps) and the fixed-beta model (alpha, eps)

library(dplyr)
library(geosphere)

setwd("") #add directory

hsa_geo <- read.csv("HSA_GEO.csv")

# start_date <- as.Date("2025-07-01")
# end_date   <- as.Date("2026-05-01")

# start_date <- as.Date("2023-03-20") # for pd1/covid
# end_date <- as.Date("2023-09-01")
# start_date<- as.Date("2023-09-01") # for pd2/covid
# end_date <- as.Date("2024-03-01")
# start_date <- as.Date("2024-03-01") # for pd3/covid
# end_date <- as.Date("2024-09-01")
# start_date <- as.Date("2024-09-01") # for pd4/covid
# end_date <- as.Date("2025-03-20")
# start_date <- as.Date("2025-03-20") # for pd5/covid
# end_date <- as.Date("2025-09-01")
start_date <- as.Date("2025-09-01") # for pd6/covid
end_date <- as.Date("2026-03-01")


# ---- data loading ----
breakpoints <- read.csv("breakpoints_covid.csv") %>%
  filter(hsa_nci_id != "All") %>%
  mutate(
    hsa_nci_id = as.numeric(hsa_nci_id),
    breakpoint = as.Date(breakpoint, format = "%m/%d/%Y"),
  ) %>%
  filter(!(hsa_nci_id %in% c(821, 822, 823))) %>%
  filter(hsa_nci_id %in% hsa_geo$hsa_nci_id) %>%
  filter(!is.na(breakpoint)) %>%
  filter(breakpoint >= start_date & breakpoint <= end_date) %>%
  group_by(hsa_nci_id) %>%
  slice_min(breakpoint, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(hsa_nci_id)

hsa_geo <- hsa_geo %>%
  filter(hsa_nci_id %in% breakpoints$hsa_nci_id) %>%
  mutate(
    total_population=as.numeric(total_population),
  )%>%
  # filter(weighted_lon > -90) %>%
  filter(total_population > 5000) %>%
  arrange(hsa_nci_id)

breakpoints <- breakpoints %>%
  filter(hsa_nci_id %in% hsa_geo$hsa_nci_id) %>%
  left_join(hsa_geo, by = c("hsa_nci_id" = "hsa_nci_id")) %>%
  arrange(hsa_nci_id)

locs <- hsa_geo %>%
  left_join(
    breakpoints %>% dplyr::select(hsa_nci_id, breakpoint),
    by = c("hsa_nci_id" = "hsa_nci_id")
  ) %>%
  mutate(
    onset = as.integer(breakpoint-start_date)+1   # day 1 = start_date
  ) %>%
  transmute(
    id = hsa_nci_id,
    pop = total_population,
    onset,
    lon = weighted_lon,
    lat = weighted_lat
  ) %>%
  arrange(id)

build_distance_matrix <- function(locs) {
  coords <- as.matrix(locs[, c("lon", "lat")])
  D <- distm(coords, coords, fun = distHaversine) / 1000  # metres -> km
  diag(D) <- NA
  D
}
D <- build_distance_matrix(locs)

# ---- three parameter MLE ----

# ilm_loglik also used in fixed beta MLE
ilm_loglik <- function(par, locs, D, log_scale = TRUE) {
  
  if (log_scale) {
    alpha <- exp(par[1])
    beta <- exp(par[2])
    eps <- exp(par[3])
  } else {
    alpha <- par[1]; beta <- par[2]; eps <- par[3]
  }
  
  n <- nrow(locs)
  Xi <- log10(locs$pop)
  # Xi <- rep(1,nrow(locs))
  onset <- locs$onset
  Tmax <- max(onset, na.rm = TRUE)
  
  Dbeta <- D^(-beta)
  diag(Dbeta) <- 0
  
  FOI <- rep(0, n)
  infected_already <- rep(FALSE, n)
  ll <- 0
  
  for (t in 1:(Tmax - 1)) {
    
    susceptible <- is.na(onset) | onset > t
    if (!any(susceptible)) break
    
    x <- alpha * Xi[susceptible] * FOI[susceptible] + eps
    
    logP <- log(-expm1(-x))
    log1mP <- -x
    
    Y <- as.numeric(onset[susceptible] == (t + 1))
    Y[is.na(Y)] <- 0
    
    ll <- ll + sum(Y*logP+(1-Y)*log1mP)
    
    newly_infected <- which(!infected_already & !is.na(onset) & onset == (t + 1))
    if (length(newly_infected) > 0) {
      FOI <- FOI + rowSums(Dbeta[, newly_infected, drop = FALSE])
      infected_already[newly_infected] <- TRUE
    }
  }
  ll
}

neg_ll <- function(par, locs, D) {
  val <- -ilm_loglik(par, locs, D)
  val
}

fit_ilm <-function(locs, D,
                   start = log(c(alpha = 4e-3, beta = 0.5, eps = 4e-3))) {
  fit <- optim(
    par     = start,
    fn      = neg_ll,
    locs    = locs,
    D       = D,
    method  = "Nelder-Mead", # alternative BGFS
    control = list(maxit = 5000, reltol = 1e-10),
    hessian = TRUE
  )
  
  est    <- exp(fit$par)
  names(est) <- c("alpha", "beta", "eps")
  
  se_log <- sqrt(diag(solve(fit$hessian)))
  ci_log <- cbind(lower = fit$par-1.96*se_log,
                  upper = fit$par+1.96*se_log)
  ci<- exp(ci_log)
  rownames(ci) <- c("alpha", "beta", "eps")
  
  list(estimate = est,
       ci95 = ci,
       loglik = -fit$value,
       convergence = fit$convergence,
       raw_fit = fit)
}

result_full <- fit_ilm(locs, D)

print(result_full$estimate)
print(result_full$ci95)
print(result_full$loglik)
aic_full <- 2*3-2*result_full$loglik
print(aic_full)

# ---- fixed beta MLE ----
f_beta=1.5

neg_ll_fixed_beta <- function(par, locs, D, beta_fixed=f_beta) {
  full_par <- c(
    alpha = par[1],
    beta = log(beta_fixed),
    eps = par[2]
  )
  
  -ilm_loglik(
    par = full_par,
    locs = locs,
    D = D
  )
}

fit_ilm_fixed_beta <- function(
    locs,
    D,
    beta_fixed = f_beta,
    start = log(c(alpha = 1e-2, eps = 1e-3))
) {
  
  fit <- optim(
    par = start,
    fn = neg_ll_fixed_beta,
    locs = locs,
    D = D,
    beta_fixed = beta_fixed,
    method = "Nelder-Mead",
    control = list(maxit = 5000, reltol = 1e-10),
    hessian = TRUE
  )
  
  est <- c(
    alpha = unname(exp(fit$par[1])),
    beta = beta_fixed,
    eps = unname(exp(fit$par[2]))
  )

  se_log <- sqrt(diag(solve(fit$hessian)))
  
  ci <- matrix(
    NA_real_,
    nrow = 3,
    ncol = 2,
    dimnames = list(
      c("alpha", "beta", "eps"),
      c("lower", "upper")
    )
  )
  
  ci[c("alpha", "eps"), ] <- exp(
    cbind(
      lower = fit$par-1.96*se_log,
      upper = fit$par+1.96*se_log
    )
  )
  
  ci["beta", ] <- beta_fixed
  
  k <- 2  # num params
  aic <- 2*k-2*(-fit$value)
  
  list(
    estimate    = est,
    ci95        = ci,
    loglik      = -fit$value,
    aic         = aic,
    convergence = fit$convergence,
    raw_fit     = fit
  )
}

result_fixed <- fit_ilm_fixed_beta(
  locs = locs,
  D = D,
  beta_fixed = f_beta
)

print(result_fixed$estimate)
print(result_fixed$ci95)
print(result_fixed$loglik)
print(result_fixed$convergence)
print(result_fixed$aic)

# --- save values for multiple beta vals ---

beta_values <- c(0, 0.5, 1, 1.5) # define fixed beta values

results_list <- lapply(beta_values, function(b) {
  message("Fitting beta = ", b)
  fit <- fit_ilm_fixed_beta(locs = locs, D = D, beta_fixed = b)
  
  data.frame(
    beta        = b,
    alpha_est   = unname(fit$estimate["alpha"]),
    alpha_low   = fit$ci95["alpha", "lower"],
    alpha_high  = fit$ci95["alpha", "upper"],
    spark_est   = unname(fit$estimate["eps"]),
    spark_low   = fit$ci95["eps", "lower"],
    spark_high  = fit$ci95["eps", "upper"],
    loglik      = fit$loglik,
    aic         = fit$aic,
    convergence = fit$convergence
  )
})

results_df <- bind_rows(results_list)

print(results_df)

write.csv(results_df, "mle_fixedbeta.csv", row.names = FALSE)

# ---- profile likelihoods, run this first -----
mle_alpha <- result_full$estimate["alpha"]
mle_beta <- result_full$estimate["beta"]
mle_eps <- result_full$estimate["eps"]

mle_alpha <- result_fixed$estimate["alpha"]
mle_beta <- result_fixed$estimate["beta"]
mle_eps <- result_fixed$estimate["eps"]

# ---- alpha profile likelihood ----

alpha_grid <- seq(mle_alpha*0.3, mle_alpha*3, length.out=95)

ll_profile_alpha <- sapply(alpha_grid, function(a) {
  ilm_loglik(
    par = c(a, mle_beta, mle_eps),
    locs = locs,
    D = D,
    log_scale = FALSE
  )
})

plot(
  alpha_grid, ll_profile_alpha,
  type = "l",
  xlab = "alpha",
  ylab = "log-likelihood"
)
#-----
alpha_mark <- mle_alpha
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

# ---- beta profile likelihood ----

beta_grid <- seq(mle_beta * 0.3, mle_beta * 3, length.out = 95)

ll_profile_beta <- sapply(beta_grid, function(b) {
  ilm_loglik(
    par = c(mle_alpha, b, mle_eps),
    locs = locs,
    D = D,
    log_scale = FALSE
  )
})

plot(
  beta_grid, ll_profile_beta,
  type = "l",
  xlab = "beta",
  ylab = "log-likelihood"
)
#-----
beta_mark <- mle_beta
idx <- which.min(abs(beta_grid - beta_mark))
x_mark <- beta_grid[idx]
y_mark <- ll_profile_beta[idx]

plot(
  beta_grid, ll_profile_beta,
  type = "l",
  xlab = "beta",
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

# ---- eps profile likelihood ----

eps_grid <- seq(mle_eps * 0.3, mle_eps * 3, length.out = 95)

ll_profile_eps <- sapply(eps_grid, function(e) {
  ilm_loglik(
    par = c(mle_alpha, mle_beta, e),
    locs = locs,
    D = D,
    log_scale = FALSE
  )
})

plot(
  eps_grid, ll_profile_eps,
  type = "l",
  xlab = "spark",
  ylab = "log-likelihood"
)
#-----
eps_mark <- mle_eps
idx <- which.min(abs(eps_grid - eps_mark))
x_mark <- eps_grid[idx]
y_mark <- ll_profile_eps[idx]

plot(
  eps_grid, ll_profile_eps,
  type = "l",
  xlab = "spark",
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
