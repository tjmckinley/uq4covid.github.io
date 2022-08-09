## load libraries
library(tidyverse)

## check if being run in batch mode
args <- commandArgs(TRUE)
if(length(args) != 0) {
    ## extract command line arguments
    args <- commandArgs(TRUE)
    if(length(args) > 0) {
        stopifnot(length(args) >= 2)
        wave <- as.numeric(args[1])
        updateJobLookup <- as.logical(args[2])
        if(length(args) > 2) {
            time <- args[3]
        } else {
            time <- "00:35:00"
        }
    } else {
        stop("No arguments")
        time <- "00:35:00"
    }
} else {
    ## set name of directory to search for outcomes
    wave <- 1
    updateJobLookup <- TRUE
}

## read in input file
pars <- readRDS(paste0("../wave", wave, "/disease.rds"))

## extract log-likelihoods
ll <- map_dbl(1:nrow(pars), function(i, wave) {
    print(i)
    ifelse(file.exists(paste0("../wave", wave, "/runs_md_cont_", i, ".rds")), 1, 0)
}, wave = wave)

## check all runs have completed
if(any(ll == 0)) {
    cat("Missing runs:\n")
    print(which(ll == 0))
    if(updateJobLookup) {
        system("rm job_lookup.txt")
        writeLines(as.character(which(is.na(ll))), "job_lookup.txt")
        code <- readLines("submit_job_template.sbatch")
        code <- gsub("RANGES", paste0("1-", sum(is.na(ll))), code)
        code <- gsub("FILEDIR", wave, code)
        code <- gsub("RUNCODE", "runDesign", code)
        code <- gsub("TIME", time, code)
        writeLines(code, "submit_job.sbatch")
    }
    stop("Stopped")
}

## cleanup
map(1:nrow(pars), function(i, wave) {
    system(paste0("rm ../wave", wave, "/wave", wave, "Forecasts_", i, ".Rout"))
}, wave = wave)
