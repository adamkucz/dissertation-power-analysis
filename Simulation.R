library(mnormt)
library(dplyr)
library(nlme)
library(doParallel)
library(foreach)

# Set seed for reproducibility
set.seed(10729)

# Data Generator ----------------------------------------------------------

#' @param n_subjects number of subjects enrolled
#' @param n_prompts number of prompts per subject total
#' @return sampled data with N 
sampleData <- function(n_subjects, n_prompts, gamma, G, sigma2){
  
  # Create participant ID numbers
  pids <- rep(seq_len(n_subjects), each =  n_prompts)
  
  # Gamma: fixed-effects
  gamma <- c("intercept" = 2,
             "loneliness_gmc" = 0.2,
             "lononeliness_pmc" = 0.2,
             "lon_pmc_l1" = 0.2)
  
  # G matrix: random effects variance/covariance matrix
  G <- matrix(c(1.00, 0.00, 0.00,
                0.00, 0.05, 0.00,
                0.00, 0.00, 0.05),
              nrow = 3, byrow = TRUE)
  
  # Within-person (level 1) variance
  sigma2 <- 1.60
  
  # Fixed effects design matrix
  X <- tibble("intercept" = rep(1, n_subjects * n_prompts),
              "loneliness" = unlist(replicate(n_subjects, arima.sim(list(order = c(1, 0, 0), ar = 0.5),
                                                                    n = n_prompts) + rnorm(1, 0, 3),
                                              simplify = F)))
  
  # Create between- and within-person centered vars and lagged loneliness
  X <- X %>%
    mutate(pid = pids,
           day = rep(1:14, each = 5, times = n_subjects),
           loneliness_gm = mean(loneliness)) %>%
    group_by(pid) %>%
    arrange(pid, day) %>%
    mutate(loneliness_pmc = loneliness - mean(loneliness),
           loneliness_gmc = mean(loneliness) - loneliness_gm,
           lon_pmc_l1 = dplyr::lag(loneliness_pmc, n = 1)) %>%
    ungroup() %>%
    select(intercept, loneliness_gmc, loneliness_pmc, lon_pmc_l1) %>%
    as.matrix()
    
  # Between-person (level 2) random effects
  # multivariate normal with mu = 0 and sigma = G
  uj <- rmnorm(n = n_subjects,
               mean = rep(0, ncol(G)),
               varcov = G)
  
  # Add column of 0s in random effect of loneliness_gmc
  # (between-subjects factors cannot be random)
  uj <- matrix(c(uj[, 1], rep(0, nrow(uj)), uj[, 2:ncol(uj)]),
               nrow = nrow(uj),
               byrow = TRUE)
  
  # repeated-measures (level 1) autocorrelated errors
  eij <- arima.sim(list(order = c(1, 0, 0), ar = 0.25),
                   n = n_subjects*n_prompts,
                   sd = sqrt(sigma2))
  
  # Matrix of beta js
  # gamma (fixed effects) with random effects added to each
  betaj <- matrix(gamma, nrow = n_subjects, ncol = length(gamma), byrow = TRUE) + uj
  
  # Simulated data
  yij <- rowSums(X * betaj[pids, ]) + eij
  
  dat <- X %>%
    as_tibble() %>%
    mutate(pid = pids,
           depressedmood = yij) %>%
    select(pid, depressedmood, loneliness_gmc, loneliness_pmc, lon_pmc_l1)
  
  return(dat)
}

# Simulation --------------------------------------------------------------

nsims <- 10000
alpha <- .05

ncores <- detectCores(logical = F)
cl <- makeCluster(ncores)
registerDoParallel(cl, cores = ncores)

fits <- foreach(i = 1:nsims,
                .packages = c("mnormt", "dplyr", "nlme"),
                .export = c("sampleData")) %dopar% {
                  
                  # Options
                  n_subj <- 100 # Number of subjects
                  n_days <- 14 # Number of days
                  n_pings <- 5 # Number of pings per day
                  
                  # Sample data
                  d <- sampleData(n_subjects = n_subj,
                                  n_prompts = n_days*n_pings)
                  
                  # Fit model
                  fit <- try({
                    lme(fixed = depressedmood ~ 1 + loneliness_gmc + loneliness_pmc + lon_pmc_l1,
                        random = ~ 1 + loneliness_pmc + lon_pmc_l1 | pid,
                        data = d,
                        correlation = corAR1(),
                        na.action = na.omit,
                        control = lmeControl(opt = "optim"))
                  })
                  
                  # If model does not converge, resample and fit again
                  while(class(fit) == "try-error"){
                    d <- sampleData(n_subjects = n_subj,
                                    n_prompts = n_days*n_pings)
                    
                    fit <- try({
                      lme(fixed = depressedmood ~ 1 + loneliness_gmc + loneliness_pmc + lon_pmc_l1,
                          random = ~ 1 + loneliness_pmc + lon_pmc_l1 | pid,
                          data = d,
                          correlation = corAR1(),
                          na.action = na.omit,
                          control = lmeControl(opt = "optim"))
                      
                    })   
                  }
                  
                  return(fit)
                  
                }

stopImplicitCluster()
stopCluster(cl)

anovas <- lapply(fits, anova)
fixedefs <- lapply(fits, fixed.effects)

pval_lonpmc <- sapply(anovas, function(x) x["loneliness_pmc", "p-value"])
pval_longmc <- sapply(anovas, function(x) x["loneliness_gmc", "p-value"])
pval_lonpmc_l1 <- sapply(anovas, function(x) x["lon_pmc_l1", "p-value"])

power_lonpmc <- sum(pval_lonpmc < alpha) / length(pval_lonpmc)
power_longmc <- sum(pval_longmc < alpha) / length(pval_longmc)
power_lonpmc_l1 <- sum(pval_lonpmc_l1 < alpha) / length(pval_lonpmc_l1)

message("Power Loneliness PMC: ", power_lonpmc)
message("Power Loneliness GMC: ", power_longmc)
message("Power Loneliness LAG1 PMC: ", power_lonpmc_l1)
