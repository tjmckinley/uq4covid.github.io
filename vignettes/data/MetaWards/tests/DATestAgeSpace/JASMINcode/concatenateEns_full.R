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
    wave <- "1"
    outputs <- "outputs"
    updateJobLookup <- TRUE
    time <- "00:35:00"
}

## read in job_lookup
jobs <- as.numeric(readLines(paste0("job_lookup_wave", wave, ".txt")))

## check that all summaries are present
runs <- map_lgl(jobs, function(i, wave) {
    run <- TRUE
    run <- ifelse(file.exists(paste0("../wave", wave, "/plotAgg_T", i, "_natFull.rds")), run, FALSE)
    run <- ifelse(file.exists(paste0("../wave", wave, "/plotAgg_T", i, "_natAgeDeathsHosp.rds")), run, FALSE)
    ## produce lad-level plots if required
    if(file.exists(paste0("lads_", outputs, ".txt"))) {
        run <- ifelse(file.exists(paste0("../wave", wave, "/plotAgg_T", i, "_age_lads.rds")), run, FALSE)
    }
    run
}, wave = wave)

if(!all(runs)) {
    cat("Missing runs:\n")
    print(jobs[which(!runs)])
    if(updateJobLookup) {
        system(paste0("rm job_lookup_wave", wave, ".txt"))
        writeLines(as.character(jobs[which(!runs)]), paste0("job_lookup_wave", wave, ".txt"))
        code <- readLines("submit_job_template.sbatch")
        code <- gsub("RANGES", paste0("1-", sum(!runs)), code)
        code <- gsub("FILEDIR", wave, code)
        code <- gsub("RUNCODE", "runPlotAgg_full", code)
        code <- gsub("OUTPUTS", outputs, code)
        code <- gsub("TIME", time, code)
        writeLines(code, paste0("submit_job_wave", wave, ".sbatch"))
    }
    stop("Stopped")
}

## cleanup
system(paste0("rm ../wave", wave, "/plot", wave, "Agg*.Rout"))

## load runs in and concatenate
runs <- map(jobs, function(i, wave) {
        readRDS(paste0("../wave", wave, "/plotAgg_T", i, "_natFull.rds"))
    }, wave = wave) %>%
    bind_rows()
    
## save output
saveRDS(runs, paste0("../wave", wave, "/sumEns_natFull.rds"))

## load runs in and concatenate
runs <- map(jobs, function(i, wave) {
        readRDS(paste0("../wave", wave, "/plotAgg_T", i, "_natAgeDeathsHosp.rds"))
    }, wave = wave) %>%
    bind_rows()
    
## save output
saveRDS(runs, paste0("../wave", wave, "/sumEns_natAgeDeathsHosp.rds"))

if(file.exists(paste0("../", outputs, "/lads_", outputs, ".txt"))) {
    ## load runs in and concatenate
    runs <- map(jobs, function(i, wave) {
            readRDS(paste0("../wave", wave, "/plotAgg_T", i, "_age_lads.rds"))
        }, wave = wave) %>%
        bind_rows()
        
    ## save output
    saveRDS(runs, paste0("../wave", wave, "/sumEns_age_lads.rds"))
}

## cleanup
system(paste0("rm ../wave", wave, "/plotAgg_T*.rds"))

