## script to remove a run

## load libraries
library(tidyverse)

## check if being run in batch mode
args <- commandArgs(TRUE)
if(length(args) != 0) {
    ## extract command line arguments
    args <- commandArgs(TRUE)
    if(length(args) > 0) {
        stopifnot(length(args) == 2)
        wave <- args[1]
        run <- as.numeric(args[2])
    } else {
        stop("No arguments")
    }
} else {
    ## set name of directory to search for outcomes
    wave <- "1Feb"
    run <- 1
}

## read in input file
pars <- readRDS(paste0("../wave", wave, "/disease.rds"))

## remove run 
pars <- pars[-run, ]
saveRDS(pars, paste0("../wave", wave, "/disease.rds"))

## read in input file
pars <- readRDS(paste0("../wave", wave, "/inputs.rds"))

## remove run 
pars <- pars[-run, ]
saveRDS(pars, paste0("../wave", wave, "/inputs.rds"))

## remove run from all simulated outputs
inds <- 1:(nrow(pars) + 1)
inds <- inds[inds > run]
if(length(inds) > 0) {
    if(dir.exists(paste0("../wave", wave, "/saveOut_wave", wave, "_", run))) {
        system(paste0("rm -r ../wave", wave, "/saveOut_wave", wave, "_", run))
    }
    for(i in inds) {
        system(paste0("mv ../wave", wave, "/runs_md_", i, ".rds ../wave", wave, "/runs_md_", i - 1, ".rds"))
        system(paste0("mv ../wave", wave, "/wave", wave, "Runs_", i, ".Rout ../wave", wave, "/wave", wave, "Runs_", i - 1, ".Rout"))
        system(paste0("mv ../wave", wave, "/saveOut_wave", wave, "_", i, " ../wave", wave, "/saveOut_wave", wave, "_", i - 1))
    }
} else {
    system(paste0("rm ../wave", wave, "/wave", wave, "Runs_", run, ".Rout"))
    if(file.exists(paste0("../wave", wave, "/runs_md_", run, ".rds"))) {
        system(paste0("rm ../wave", wave, "/runs_md_", run, ".rds"))
    }
    if(dir.exists(paste0("../wave", wave, "/saveOut_wave", wave, "_", run))) {
        system(paste0("rm -r ../wave", wave, "/saveOut_wave", wave, "_", run))
    }
}

cat("Finished\n")


