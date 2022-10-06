## load libraries
library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(sitmo)
library(abind)
library(parallel)
library(data.table)
library(R.utils)

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
    hash <- 1
}

## set case specific values
tstart <- 50
tstop <- 75
npart <- 50
a1 <- 0
a2 <- 0
b <- 0.1
b_dis <- 0.05
sigma2_lad <- 1
sigma2_age_region <- 1
sigma2_nhsregion <- 1
sigma2_age_nhsregion <- 1

## source Rcpp PF code
sourceCpp("BPF.cpp")

## source function to run PF and return log-likelihood
source("BPF.R")

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

## load lookups
lookup <- readRDS("outputs/lookup.rds")
age_lookup <- readRDS("outputs/age_lookup.rds")

## set up output folder
folder <- paste0("wave", wave, "/saveOut_", hash)

## set up inputs
u1_moves <- read_csv(paste0(folder, "/snapshot_u1moves_t", tstart, ".csv.bz2"), col_names = FALSE) %>%
    as.matrix()
colnames(u1_moves) <- NULL
ncohorts1 <- read_csv(paste0(folder, "/snapshot_ncohorts1_t", tstart, ".csv.bz2"), col_names = FALSE)$X1 %>%
    as.vector()
ncohorts2 <- read_csv(paste0(folder, "/snapshot_ncohorts2_t", tstart, ".csv.bz2"), col_names = FALSE)$X1 %>%
    as.vector()
playprobs <- read_csv(paste0(folder, "/snapshot_playprobs_t", tstart, ".csv.bz2"), col_names = FALSE)$X1 %>%
    as.vector()

## read in snapshots
u1 <- list()
u2 <- list()
for(i in 1:npart) {
    u1[[i]] <- fread(paste0(folder, "/snapshot_u1_t", tstart, "_", i - 1, ".csv.bz2"))
    udims <- c(max(u1[[i]]$class) + 1, ncol(u1[[i]]) - 2, max(u1[[i]]$cohort))
    u1[[i]][, class := NULL]
    u1[[i]][, cohort := NULL]
    u1[[i]] <- as.matrix(u1[[i]])
    u1[[i]] <- aperm(array(t(u1[[i]]), c(udims[2], udims[1], udims[3])), c(2, 1, 3))
        
    u2[[i]] <- fread(paste0(folder, "/snapshot_u2_t", tstart, "_", i - 1, ".csv.bz2"))
    udims <- c(max(u2[[i]]$class) + 1, ncol(u2[[i]]) - 2, max(u2[[i]]$lad))
    u2[[i]][, class := NULL]
    u2[[i]][, lad := NULL]
    u2[[i]] <- as.matrix(u2[[i]])
    u2[[i]] <- aperm(array(t(u2[[i]]), c(udims[2], udims[1], udims[3])), c(2, 1, 3))
}

u <- list(u1_moves = u1_moves, u1 = u1, u2 = u2, ncohorts1 = ncohorts1, ncohorts2 = ncohorts2, playprobs = playprobs)

## run forecasts        
runs_md <- BPF(pars[hash, ], C1 = contact1, C2 = contact2, lockdown_day = 20,
    lookup = lookup, age_lookup = age_lookup, u = u,
    tstart = tstart, tstop = tstop, npart = npart,
    a1 = a1, a2 = a2, b = b, b_dis = b_dis,
    sigma2_lad = sigma2_lad, sigma2_age_region = sigma2_age_region, 
    sigma2_nhsregion = sigma2_nhsregion, sigma2_age_nhsregion = sigma2_age_nhsregion,
    saveAll = TRUE, writeExt = TRUE, outputName = paste0("wave", wave, "/saveOut_", hash),
    PF = FALSE, ncores = 1)

## save outputs
saveRDS(runs_md, paste0("wave", wave, "/runs_md_cont_", hash, ".rds"))
