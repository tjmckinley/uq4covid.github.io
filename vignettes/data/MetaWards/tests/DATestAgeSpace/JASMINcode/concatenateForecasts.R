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
    wave <- 1
    outputs <- "outputs"
    updateJobLookup <- TRUE
    time <- "00:35:00"
}

## read in input file
pars <- readRDS(paste0("../wave", wave, "/disease.rds"))

## check runs
runs <- map_dbl(1:nrow(pars), function(i, wave) {
    print(i)
    file <- paste0("../wave", wave, "/wave", wave, "Forecasts_", i, ".Rout")
    out <- NA
    if(file.exists(file)) {
        file <- readLines(file)
        if(length(grep("Finished", file)) == 1) {
            out <- 1
        }
    }
    out
}, wave = wave)

## check all runs have completed
if(any(is.na(runs))) {
    cat("Missing runs:\n")
    print(which(is.na(runs)))
    if(updateJobLookup) {
        system(paste0("rm job_lookup_wave", wave, ".txt"))
        writeLines(as.character(which(is.na(runs))), paste0("job_lookup_wave", wave, ".txt"))
        code <- readLines("submit_job_template.sbatch")
        code <- gsub("RANGES", paste0("1-", sum(is.na(runs))), code)
        code <- gsub("FILEDIR", wave, code)
        code <- gsub("RUNCODE", "runForecasts", code)
        code <- gsub("OUTPUTS", outputs, code)
        code <- gsub("TIME", time, code)
        writeLines(code, paste0("submit_job_wave", wave, ".sbatch"))
    }
    stop("Stopped")
}

## cleanup
#map(1:nrow(pars), function(i, wave) {
#    system(paste0("rm ../wave", wave, "/wave", wave, "Forecasts_", i, ".Rout"))
#}, wave = wave)

