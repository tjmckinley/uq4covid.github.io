## load libraries
library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(parallel)

## source simulation model
sourceCpp("discreteStochModel.cpp")

## read in truncated Skellam sampler
source("trSkellam.R")

## source function to run PF and return log-likelihood
source("PF1.R")

## set wave
wave <- 1

## read in simulated data and generate incidence curves
data <- readRDS("outputs/disSims.rds")

## read in parameters, remove guff and reorder
pars <- readRDS(paste0("wave", wave, "/disease.rds")) %>%
    rename(nu = `beta[1]`, nuA = `beta[6]`) %>%
    select(!c(starts_with("beta"), repeats)) %>%
    select(nu, nuA, !output)

## read in contact matrix
contact <- read_csv("inputs/POLYMOD_matrix.csv", col_names = FALSE) %>%
    as.matrix()

## solution to round numbers preserving sum
## adapted from:
## https://stackoverflow.com/questions/32544646/round-vector-of-numerics-to-integer-while-preserving-their-sum
smart_round <- function(x) {
    y <- floor(x)
    indices <- tail(order(x - y), round(sum(x)) - sum(y))
    y[indices] <- y[indices] + 1
    y
}

## set up number of initial individuals in each age-class
N <- 10000
N <- smart_round(read_csv("inputs/age_seeds.csv", col_names = FALSE)$X2 * N)
I0 <- smart_round(read_csv("inputs/age_seeds.csv", col_names = FALSE)$X2 * 1)
S0 <- N - I0

## set initial counts
u <- matrix(0, 12, 8)
u[1, ] <- S0
u[2, ] <- I0

## run PF with some model discrepancy
runs_md <- PF1(pars, C = contact, data = data, u = u, ndays = 30, npart = 100, MD = TRUE, 
    a_dis = 0.05, b_dis = 0.05, a1 = 0.01, a2 = 0.2, b = 0.1, saveAll = NA)

## save outputs
saveRDS(runs_md, paste0("wave", wave, "/runs_md.rds"))
