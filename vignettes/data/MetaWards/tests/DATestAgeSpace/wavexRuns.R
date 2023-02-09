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
        stopifnot(length(args) == 3)
        wave <- args[1]
        hash <- as.numeric(args[2])
        outputs <- args[3]
    } else {
        stop("No arguments")
    }
} else {
    ## set wave number
    wave <- "1"
    hash <- 1
    outputs <- "outputs"
}

## read in fixed input data
fixedInputs <- readLines(paste0("wave", wave, "/fixedInputs.txt"))

## set case specific values
tstart <- as.numeric(fixedInputs[1])
tstop <- as.numeric(fixedInputs[2])
lockdown_day <- as.numeric(fixedInputs[3])
npart <- as.numeric(fixedInputs[4])
niter <- as.numeric(fixedInputs[5])
a1 <- as.numeric(fixedInputs[6])
a2 <- as.numeric(fixedInputs[7])
b1 <- as.numeric(fixedInputs[8])
b2 <- as.numeric(fixedInputs[9])
a_dis <- as.numeric(fixedInputs[10])
b_dis <- as.numeric(fixedInputs[11])
b_dis_age8 <- as.numeric(fixedInputs[12])
sigma2_lad <- as.numeric(fixedInputs[13])
sigma2_age_region <- as.numeric(fixedInputs[14])
sigma2_nhsregion <- as.numeric(fixedInputs[15])
sigma2_age_nhsregion <- as.numeric(fixedInputs[16])
saveAll <- as.logical(as.numeric(fixedInputs[17]))
snapshot <- as.logical(as.numeric(fixedInputs[18]))
writeExt <- as.logical(as.numeric(fixedInputs[19]))

## source Rcpp PF code
sourceCpp("BPF.cpp")

## source function to run PF and return log-likelihood
source("BPF.R")

## read in simulated data
cumDeath_lad <- readRDS(paste0(outputs, "/cumDeath_lad.rds"))
cumDeath_age_region <- readRDS(paste0(outputs, "/cumDeath_age_region.rds"))
hosp_nhsregion <- readRDS(paste0(outputs, "/hosp_nhsregion.rds"))
cumHospAd_age_nhsregion <- readRDS(paste0(outputs, "/cumHospAd_age_nhsregion.rds"))

## read in parameters, remove guff and reorder
pars <- readRDS(paste0("wave", wave, "/disease.rds")) %>%
    select(!output)

## read in contact matrix
contact1 <- read_csv("inputs/POLYMOD_matrix.csv", col_names = FALSE) %>%
    as.matrix()
contact2 <- read_csv("inputs/coMix_matrix.csv", col_names = FALSE) %>%
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

## add age probabilities
ageProbs <- read_csv("inputs/age_seeds.csv", col_names = FALSE)$X2

## read in commuter data
EW19 <- read_delim("inputs/EW19.dat", delim = " ", col_names = FALSE)
u1_moves <- as.matrix(EW19[, 1:2])

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
lookup <- readRDS(paste0(outputs, "/lookup.rds"))
age_lookup <- readRDS(paste0(outputs, "/age_lookup.rds"))

## set up inputs
u <- list(u1 = u1, u2 = u2, u1_moves = u1_moves, u2_moves = as.matrix(PM19))

## create model discrepancy matrices
a_dis <- a_dis * exp(-(pars$MD_time[hash] - tstart:tstop) * ifelse(tstart:tstop < pars$MD_time[hash], 1, 0))
a_dis <- array(rep(a_dis, nrow(age_lookup) * nrow(lookup)), c(tstop - tstart + 1, nrow(age_lookup), nrow(lookup)))
b_dis <- b_dis * exp(-(pars$MD_time[hash] - tstart:tstop) * ifelse(tstart:tstop < pars$MD_time[hash], 1, 0))
b_dis <- array(rep(b_dis, nrow(age_lookup) * nrow(lookup)), c(tstop - tstart + 1, nrow(age_lookup), nrow(lookup)))
b_dis_age8 <- b_dis_age8 * exp(-(pars$MD_time[hash] - tstart:tstop) * ifelse(tstart:tstop < pars$MD_time[hash], 1, 0))
b_dis[, 8, ] <- rep(b_dis_age8, dim(b_dis)[3])

## run PF with some model discrepancy
if(exists("hash")) {
    runs_md <- BPF(pars[hash, ], C1 = contact1, C2 = contact2, lockdown_day = lockdown_day,
        cumDeath_lad = cumDeath_lad, cumDeath_age_region = cumDeath_age_region,
        hosp_nhsregion = hosp_nhsregion, cumHospAd_age_nhsregion = cumHospAd_age_nhsregion,
        lookup = lookup, age_lookup = age_lookup, u = u,
        tstart = tstart, tstop = tstop, npart = npart, niter = niter,
        a1 = a1, a2 = a2, b1 = b1, b2 = b2, a_dis = a_dis, b_dis = b_dis,
        sigma2_lad = sigma2_lad, sigma2_age_region = sigma2_age_region, 
        sigma2_nhsregion = sigma2_nhsregion, sigma2_age_nhsregion = sigma2_age_nhsregion,
        saveAll = saveAll, snapshot = snapshot, writeExt = writeExt, 
        outputName = paste0("saveOut_wave", wave, "_", hash),
        ncores = 1)
    ## save outputs
    saveRDS(runs_md, paste0("wave", wave, "/runs_md_", hash, ".rds"))
    if(writeExt) system(paste0("mv saveOut_wave", wave, "_", hash, " wave", wave))
} else {
    runs_md <- BPF(pars, C1 = contact1, C2 = contact2, lockdown_day = 20,
        cumDeath_lad = cumDeath_lad, cumDeath_age_region = cumDeath_age_region, 
        hosp_nhsregion = hosp_nhsregion, cumHospAd_age_nhsregion = cumHospAd_age_nhsregion, 
        lookup = lookup, age_lookup = age_lookup, u = u,
        tstart = tstart, tstop = tstop, npart = npart, niter = niter,
        a1 = a1, a2 = a2, b1 = b1, b2 = b2, a_dis = a_dis, b_dis = b_dis,
        sigma2_lad = sigma2_lad, sigma2_age_region = sigma2_age_region, 
        sigma2_nhsregion = sigma2_nhsregion, sigma2_age_nhsregion = sigma2_age_nhsregion,
        saveAll = NA)
    ## save outputs
    saveRDS(runs_md, paste0("wave", wave, "/runs_md.rds"))
}

