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
fixedInputs <- readLines(paste0("wave", wave, "/fixedInputs_cont.txt"))

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
a_dis_ini <- as.numeric(fixedInputs[10])
b_dis_ini <- as.numeric(fixedInputs[11])
b_dis_8_ini <- as.numeric(fixedInputs[12])
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

## load lookups
lookup <- readRDS(paste0(outputs, "/lookup.rds"))
age_lookup <- readRDS(paste0(outputs, "/age_lookup.rds"))

## set up output folder
folder <- paste0("wave", wave, "/saveOut_wave", wave, "_", hash)

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
region_lookup <- readRDS("data/region_lookup.rds") %>%
    mutate(initials = NA) %>%
    mutate(initials = ifelse(RGN19NM == "North East", "NE", initials)) %>%
    mutate(initials = ifelse(RGN19NM == "North West", "NW", initials)) %>%
    mutate(initials = ifelse(RGN19NM == "Yorkshire and The Humber", "YH", initials)) %>%
    mutate(initials = ifelse(RGN19NM == "East Midlands", "EM", initials)) %>%
    mutate(initials = ifelse(RGN19NM == "West Midlands", "WM", initials)) %>%
    mutate(initials = ifelse(RGN19NM == "East of England", "EE", initials)) %>%
    mutate(initials = ifelse(RGN19NM == "London", "L", initials)) %>%
    mutate(initials = ifelse(RGN19NM == "South East", "SE", initials)) %>%
    mutate(initials = ifelse(RGN19NM == "South West", "SW", initials))

## set some initial matrices
a_dis <- a_dis_ini * exp(-pars$MD_scale[hash] * (pars$MD_time_L[hash] - tstart:tstop) * ifelse(tstart:tstop < pars$MD_time_L[hash], 1, 0))
a_dis <- array(rep(a_dis, nrow(age_lookup) * nrow(lookup)), c(tstop - tstart + 1, nrow(age_lookup), nrow(lookup)))
b_dis <- b_dis_ini * exp(-pars$MD_scale[hash] * (pars$MD_time_L[hash] - tstart:tstop) * ifelse(tstart:tstop < pars$MD_time_L[hash], 1, 0))
b_dis <- array(rep(b_dis, nrow(age_lookup) * nrow(lookup)), c(tstop - tstart + 1, nrow(age_lookup), nrow(lookup)))

## use spatially-explicit matrices
for(i in 1:nrow(region_lookup)) {
    FID <- lookup$FID[!is.na(lookup$FID_region) & lookup$FID_region == region_lookup$FID[region_lookup$RGN19NM == region_lookup$RGN19NM[i]]]
    MD_time <- pluck(pars, paste0("MD_time_", region_lookup$initials[i]))
    a_dis_temp <- a_dis_ini * exp(-pars$MD_scale[hash] * (MD_time[hash] - tstart:tstop) * ifelse(tstart:tstop < MD_time[hash], 1, 0))
    a_dis[, , FID] <- rep(a_dis_temp, length(FID))
    b_dis_temp <- b_dis_ini * exp(-pars$MD_scale[hash] * (MD_time[hash] - tstart:tstop) * ifelse(tstart:tstop < MD_time[hash], 1, 0))
    b_dis[, , FID] <- rep(b_dis_temp, length(FID))
    b_dis_8 <- b_dis_8_ini * exp(-pars$MD_scale[hash] * (MD_time[hash] - tstart:tstop) * ifelse(tstart:tstop < MD_time[hash], 1, 0))
    b_dis[, 8, FID] <- rep(b_dis_8, length(FID))
}

## run forecasts
if(exists("hash")) {
    runs_md <- BPF(pars[hash, ], C1 = contact1, C2 = contact2, lockdown_day = lockdown_day,
        cumDeath_lad = cumDeath_lad, cumDeath_age_region = cumDeath_age_region,
        hosp_nhsregion = hosp_nhsregion, cumHospAd_age_nhsregion = cumHospAd_age_nhsregion,
        lookup = lookup, age_lookup = age_lookup, u = u,
        tstart = tstart, tstop = tstop, npart = npart, niter = niter,
        a1 = a1, a2 = a2, b1 = b1, b2 = b2, a_dis = a_dis, b_dis = b_dis,
        sigma2_lad = sigma2_lad, sigma2_age_region = sigma2_age_region, 
        sigma2_nhsregion = sigma2_nhsregion, sigma2_age_nhsregion = sigma2_age_nhsregion,
        saveAll = saveAll, snapshot = snapshot, writeExt = writeExt, PF = FALSE, 
        outputName = paste0("wave", wave, "/saveOut_wave", wave, "_", hash),
        ncores = 1)
} else {
    runs_md <- BPF(pars, C1 = contact1, C2 = contact2, lockdown_day = lockdown_day,
        cumDeath_lad = cumDeath_lad, cumDeath_age_region = cumDeath_age_region, 
        hosp_nhsregion = hosp_nhsregion, cumHospAd_age_nhsregion = cumHospAd_age_nhsregion, 
        lookup = lookup, age_lookup = age_lookup, u = u,
        tstart = tstart, tstop = tstop, npart = npart, niter = niter,
        a1 = a1, a2 = a2, b1 = b1, b2 = b2, a_dis = a_dis, b_dis = b_dis,
        sigma2_lad = sigma2_lad, sigma2_age_region = sigma2_age_region, 
        sigma2_nhsregion = sigma2_nhsregion, sigma2_age_nhsregion = sigma2_age_nhsregion,
        saveAll = NA, PF = FALSE)
    ## save outputs
    saveRDS(runs_md, paste0("wave", wave, "/runs_md_cont.rds"))
}

print("Finished")
