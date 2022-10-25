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
tstart <- 50
tstop <- 75
print("Need to figure out best way to pass forecast stuff in")

lockdown_day <- as.numeric(fixedInputs[3])
npart <- as.numeric(fixedInputs[4])
niter <- as.numeric(fixedInputs[5])
a1 <- as.numeric(fixedInputs[6])
a2 <- as.numeric(fixedInputs[7])
b1 <- as.numeric(fixedInputs[8])
b2 <- as.numeric(fixedInputs[9])
a_dis <- as.numeric(fixedInputs[10])
b_dis <- as.numeric(fixedInputs[11])
b_dis_london <- as.numeric(fixedInputs[12])
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

## read in parameters, remove guff and reorder
pars <- readRDS(paste0("wave", wave, "/disease.rds")) %>%
    select(!output)

## read in contact matrix
contact1 <- read_csv("inputs/POLYMOD_matrix.csv", col_names = FALSE) %>%
    as.matrix()
contact2 <- read_csv("inputs/coMix_matrix.csv", col_names = FALSE) %>%
    as.matrix()

## load lookups
lookup <- readRDS(paste0(outputs, "/lookup.rds"))
age_lookup <- readRDS(paste0(outputs, "/age_lookup.rds"))

## set up output folder
folder <- paste0("wave", wave, "/saveOut_", outputs, "_", hash)

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

## create model discrepancy matrices
region_lookup <- readRDS("data/region_lookup.rds")
london_FID <- lookup$FID[!is.na(lookup$FID_region) & lookup$FID_region == region_lookup$FID[region_lookup$RGN19NM == "London"]]
b_dis <- matrix(rep(b_dis, nrow(age_lookup) * nrow(lookup)), nrow(age_lookup), nrow(lookup))
b_dis[, london_FID] <- b_dis_london

## run forecasts        
runs_md <- BPF(pars[hash, ], C1 = contact1, C2 = contact2, lockdown_day = 20,
    lookup = lookup, age_lookup = age_lookup, u = u,
    tstart = tstart, tstop = tstop, npart = npart,
    a1 = a1, a2 = a2, b1 = b1, b2 = b2, a_dis = a_dis, b_dis = b_dis,
    sigma2_lad = sigma2_lad, sigma2_age_region = sigma2_age_region, 
    sigma2_nhsregion = sigma2_nhsregion, sigma2_age_nhsregion = sigma2_age_nhsregion,
    saveAll = TRUE, writeExt = TRUE, outputName = paste0("wave", wave, "/saveOut_", hash),
    PF = FALSE, ncores = 1)

## save outputs
saveRDS(runs_md, paste0("wave", wave, "/runs_md_cont_", hash, ".rds"))
