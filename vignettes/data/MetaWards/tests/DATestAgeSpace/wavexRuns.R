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

## set case specific values
tstart <- 0
tstop <- 50
lockdown_day <- 20
npart <- 50
niter <- 10
a1 <- 0
a2 <- 0
b1 <- 0.1
b2 <- 0.1
b_dis <- 0.01
sigma2_lad <- 1
sigma2_age_region <- 1
sigma2_nhsregion <- 1
sigma2_age_nhsregion <- 1
saveAll <- TRUE
snapshot <- TRUE
writeExt <- TRUE

## source Rcpp PF code
sourceCpp("BPF.cpp")

## source function to run PF and return log-likelihood
source("BPF.R")

## read in simulated data
cumDeath_lad <- readRDS("outputs/cumDeath_lad.rds")
cumDeath_age_region <- readRDS("outputs/cumDeath_age_region.rds")
hosp_nhsregion <- readRDS("outputs/hosp_nhsregion.rds")
cumHospAd_age_nhsregion <- readRDS("outputs/cumHospAd_age_nhsregion.rds")

## read in parameters, remove guff and reorder
pars <- readRDS(paste0("wave", wave, "/disease.rds")) %>%
    rename(nu = `beta[1]`, nuA = `beta[6]`) %>%
    select(!c(starts_with("beta["), repeats)) %>%
    select(nu, nuA, !c(beta_scale, p_move, a_dis, output), beta_scale, p_move, a_dis)

## read in contact matrix
contact1 <- read_csv("inputs/POLYMOD_matrix.csv", col_names = FALSE) %>%
    as.matrix()
contact2 <- read_csv("inputs/coMix_matrix.csv", col_names = FALSE) %>%
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

## read in player data
PM19 <- read_delim("inputs/PlayMatrix19.dat", delim = " ", col_names = FALSE)
PlaySize19 <- read_delim("inputs/PlaySize19.dat", delim = " ", col_names = FALSE)
  
## expand to deal with age-classes
u2 <- apply(PlaySize19, 1, function(x, ageProbs) {
        u <- matrix(0, 12, length(ageProbs))
        u[1, ] <- smart_round(ageProbs * x[2])
        list(u)
    }, ageProbs = ageProbs) %>%
    map(1) %>%
    abind(along = 3)
    
## load lookups
lookup <- readRDS("outputs/lookup.rds")
age_lookup <- readRDS("outputs/age_lookup.rds")

## set up inputs
u <- list(u1 = u1, u2 = u2, u1_moves = u1_moves, u2_moves = as.matrix(PM19))

## run PF with some model discrepancy
if(exists("hash")) {
    runs_md <- BPF(pars[hash, ], C1 = contact1, C2 = contact2, lockdown_day = lockdown_day,
        cumDeath_lad = cumDeath_lad, cumDeath_age_region = cumDeath_age_region,
        hosp_nhsregion = hosp_nhsregion, cumHospAd_age_nhsregion = cumHospAd_age_nhsregion,
        lookup = lookup, age_lookup = age_lookup, u = u,
        tstart = tstart, tstop = tstop, npart = npart, niter = niter,
        a1 = a1, a2 = a2, b1 = b1, b2 = b2, b_dis = b_dis,
        sigma2_lad = sigma2_lad, sigma2_age_region = sigma2_age_region, 
        sigma2_nhsregion = sigma2_nhsregion, sigma2_age_nhsregion = sigma2_age_nhsregion,
        saveAll = saveAll, snapshot = snapshot, writeExt = writeExt, 
        outputName = paste0("saveOut_", hash),
        ncores = 1)
    ## save outputs
    saveRDS(runs_md, paste0("wave", wave, "/runs_md_", hash, ".rds"))
    if(writeExt) system(paste0("mv saveOut_", hash, " wave", wave))
} else {
    runs_md <- BPF(pars, C1 = contact1, C2 = contact2, lockdown_day = 20,
        cumDeath_lad = cumDeath_lad, cumDeath_age_region = cumDeath_age_region, 
        hosp_nhsregion = hosp_nhsregion, cumHospAd_age_nhsregion = cumHospAd_age_nhsregion, 
        lookup = lookup, age_lookup = age_lookup, u = u,
        tstart = tstart, tstop = tstop, npart = npart, niter = niter,
        a1 = a1, a2 = a2, b1 = b1, b2 = b2, b_dis = b_dis,
        sigma2_lad = sigma2_lad, sigma2_age_region = sigma2_age_region, 
        sigma2_nhsregion = sigma2_nhsregion, sigma2_age_nhsregion = sigma2_age_nhsregion,
        saveAll = NA)
    ## save outputs
    saveRDS(runs_md, paste0("wave", wave, "/runs_md.rds"))
}

