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
    }
} else {
    ## set name of directory to search for outcomes
    wave <- 1
    updateJobLookup <- TRUE
    time <- "00:35:00"
}

## read in input file
pars <- readRDS(paste0("../wave", wave, "/disease.rds"))

## check that all summaries are present
runs <- map_lgl(1:nrow(pars), function(i, wave) {
    run <- TRUE
    run <- ifelse(file.exists(paste0("../wave", wave, "/plotSum_", i, "_natFull.rds")), run, FALSE)
    run <- ifelse(file.exists(paste0("../wave", wave, "/plotSum_", i, "_natDeaths.rds")), run, FALSE)
    run <- ifelse(file.exists(paste0("../wave", wave, "/plotSum_", i, "_ageRegionDeaths.rds")), run, FALSE)
    run <- ifelse(file.exists(paste0("../wave", wave, "/plotSum_", i, "_nhsregionHosp.rds")), run, FALSE)
    run <- ifelse(file.exists(paste0("../wave", wave, "/plotSum_", i, "_ageNhsregionHosp.rds")), run, FALSE)
    ## produce lad-level plots if required
    if(file.exists("lads.txt")) {
        run <- ifelse(file.exists(paste0("../wave", wave, "/plotSum_", i, "_lads.rds")), run, FALSE)
    }
    run
}, wave = wave)

if(!all(runs)) {
    cat("Missing runs:\n")
    print(which(!runs))
    if(updateJobLookup) {
        system("rm job_lookup.txt")
        writeLines(as.character(which(!runs)), "job_lookup.txt")
        code <- readLines("submit_job_template.sbatch")
        code <- gsub("RANGES", paste0("1-", sum(!runs)), code)
        code <- gsub("FILEDIR", wave, code)
        code <- gsub("RUNCODE", "runPlotSum", code)
        code <- gsub("TIME", time, code)
        writeLines(code, "submit_job.sbatch")
    }
    stop("Stopped")
}

## cleanup
system(paste0("rm ../wave", wave, "/plot", wave, "Sum*.Rout"))

print("Finished")
