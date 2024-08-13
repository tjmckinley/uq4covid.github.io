## load libraries
library(tidyverse)

## check if being run in batch mode
args <- commandArgs(TRUE)
if(length(args) != 0) {
    ## extract command line arguments
    args <- commandArgs(TRUE)
    if(length(args) > 0) {
        stopifnot(length(args) >= 3)
        wave <- args[1]
        outputs <- args[2]
        updateJobLookup <- as.logical(args[3])
        if(length(args) > 3) {
            time <- args[4]
        } else {
            time <- "00:35:00"
        }
    } else {
        stop("No arguments")
    }
} else {
    ## set name of directory to search for outcomes
    wave <- "1Feb"
    outputs <- "Feb"
    updateJobLookup <- TRUE
    time <- "00:35:00"
}

## append outputs
outputs <- paste0("outputs", outputs)

## read in input file
pars <- readRDS(paste0("../wave", wave, "/disease.rds"))

## extract log-likelihoods
ll <- map_dbl(1:nrow(pars), function(i, wave) {
    print(i)
    if(file.exists(paste0("../wave", wave, "/runs_md_", i, ".rds"))) {
        run <- readRDS(paste0("../wave", wave, "/runs_md_", i, ".rds"))$ll
    } else {
        run <- NA
    }
    run
}, wave = wave)

## check all runs have completed
if(any(is.na(ll))) {
    cat("Missing runs:\n")
    print(which(is.na(ll)))
    if(updateJobLookup) {
        system(paste0("rm job_lookup_wave", wave, ".txt"))
        writeLines(as.character(which(is.na(ll))), paste0("job_lookup_wave", wave, ".txt"))
        code <- readLines("submit_job_template.sbatch")
        code <- gsub("RANGES", paste0("1-", sum(is.na(ll))), code)
        code <- gsub("FILEDIR", wave, code)
        code <- gsub("RUNCODE", "runDesign_full", code)
        code <- gsub("OUTPUTS", outputs, code)
        code <- gsub("TIME", time, code)
        writeLines(code, paste0("submit_job_wave", wave, ".sbatch"))
    }
    stop("Stopped")
}

## extract runtimes
times <- map_dbl(1:nrow(pars), function(i, wave) {
    time <- readLines(paste0("../wave", wave, "/wave", wave, "Runs_full_", i, ".Rout"))
    time <- time[length(time)]
    time <- strsplit(time, " ")[[1]]
    as.numeric(time[length(time)])
}, wave = wave)
times <- times / 60
times <- times / 60
saveRDS(times, paste0("../wave", wave, "/times.rds"))

## cleanup
#map(1:nrow(pars), function(i, wave) {
#    system(paste0("rm ../wave", wave, "/wave", wave, "Runs_", i, ".Rout"))
#}, wave = wave)

## save
saveRDS(ll, paste0("../wave", wave, "/ll.rds"))
