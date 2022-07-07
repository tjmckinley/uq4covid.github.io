## load libraries
library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(sitmo)
library(abind)
library(parallel)

## check if being run in batch mode
args <- commandArgs(TRUE)
if(length(args) != 0) {
    ## extract command line arguments
    args <- commandArgs(TRUE)
    if(length(args) > 0) {
        stopifnot(length(args) == 2)
        wave <- as.numeric(args[1])
        hash <- as.numeric(args[2])
    } else {
        stop("No arguments")
    }
} else {
    ## set wave number
    wave <- 1
}

## source Rcpp PF code
sourceCpp("APF1.cpp")

## source function to run PF and return log-likelihood
source("APF1.R")

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

## read in initial conditions
#u1 <- readRDS("outputs/u1.rds")
u1_moves <- readRDS("outputs/u1_moves.rds")

## solution to round numbers preserving sum
## adapted from:
## https://stackoverflow.com/questions/32544646/round-vector-of-numerics-to-integer-while-preserving-their-sum
smart_round <- function(x) {
    y <- floor(x)
    indices <- tail(order(x - y), round(sum(x)) - sum(y))
    y[indices] <- y[indices] + 1
    y
}

## add age probabilities
ageProbs <- read_csv("inputs/age_seeds.csv", col_names = FALSE)$X2

## read in commuter data
EW19 <- read_delim("inputs/EW19.dat", delim = " ", col_names = FALSE)

## expand to deal with age-classes
u1 <- apply(EW19, 1, function(x, ageProbs) {
        u <- matrix(0, 12, length(ageProbs))
        u[1, ] <- smart_round(ageProbs * x[3])
        list(u)
    }, ageProbs = ageProbs) %>%
    map(1) %>%
    abind(along = 3)

## run PF with some model discrepancy
if(exists("hash")) {
    runs_md <- APF1(pars[hash, ], C = contact, data = data, u1_moves = u1_moves,
        u1 = u1, ndays = 100, npart = 50, a_dis = 0.05, b_dis = 0.05,
        a1 = 0.01, a2 = 0.2, b = 0.001, saveAll = NA, ncores = 1)
    ## save outputs
    saveRDS(runs_md, paste0("wave", wave, "/runs_md_", hash, ".rds"))
} else {
    runs_md <- APF1(pars, C = contact, data = data, u1_moves = u1_moves,
        u1 = u1, ndays = 100, npart = 50, a_dis = 0.05, b_dis = 0.05,
        a1 = 0.01, a2 = 0.2, b = 0.001, saveAll = NA)
    ## save outputs
    saveRDS(runs_md, paste0("wave", wave, "/runs_md.rds"))
}

