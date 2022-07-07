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

## read in job_lookup
jobs <- as.numeric(readLines("job_lookup.txt"))

## check that all summaries are present
runs <- map_lgl(jobs, function(i, wave, lads) {
    if(is.na(lads[1])) {
        if(file.exists(paste0("../wave", wave, "/plotAgg_T", i, ".rds"))) {
            run <- TRUE
        } else {
            run <- FALSE
        }
    } else {
        stop("Not yet implemented for individual LADs")
    }
    run
}, wave = wave, lads = NA)

if(!all(runs)) {
    cat("Missing runs:\n")
    print(jobs[which(!runs)])
    if(updateJobLookup) {
        system("rm job_lookup.txt")
        writeLines(as.character(jobs[which(!runs)]), "job_lookup.txt")
        code <- readLines("submit_job_template.sbatch")
        code <- gsub("RANGES", paste0("1-", sum(!runs)), code)
        code <- gsub("FILEDIR", wave, code)
        code <- gsub("RUNCODE", "runPlotAgg", code)
        code <- gsub("TIME", time, code)
        writeLines(code, "submit_job.sbatch")
    }
    stop("Stopped")
}

## cleanup
system(paste0("rm ../wave", wave, "/plot", wave, "Agg*.Rout"))

## load runs in and concatenate
runs <- map(jobs, function(i, wave, lads) {
        if(is.na(lads[1])) {
            run <- readRDS(paste0("../wave", wave, "/plotAgg_T", i, ".rds"))
        } else {
            stop("Not yet implemented for individual LADs")
        }
        run
    }, wave = wave, lads = NA) %>%
    bind_rows()
    
## save output
saveRDS(runs, paste0("../wave", wave, "/sumEns.rds"))

## cleanup
system(paste0("rm ../wave", wave, "/plotAgg_T*.rds"))

