## load libraries
library(tidyverse)

## check if being run in batch mode
args <- commandArgs(TRUE)
if(length(args) != 0) {
    ## extract command line arguments
    args <- commandArgs(TRUE)
    if(length(args) > 0) {
        stopifnot(length(args) == 2)
        wave <- as.numeric(args[1])
        t <- as.numeric(args[2])
    } else {
        stop("No arguments")
    }
} else {
    ## set name of directory to search for outcomes
    wave <- 1
    t <- 1
}

## read in input file
pars <- readRDS(paste0("wave", wave, "/disease.rds"))

## concatenate runs over ensemble
runs <- map(1:nrow(pars), function(i, t, wave) {
        readRDS(paste0("wave", wave, "/plotSum_", i, ".rds")) %>%
            filter(time == t) %>%
            select(n, age, time, var)
    }, t = t, wave = wave) %>%
    bind_rows() %>%
    group_by(age, time, var) %>%
    summarise(
        LCI = quantile(n, probs = 0.025),
        LQ = quantile(n, probs = 0.25),
        Median = quantile(n, probs = 0.5),
        UQ = quantile(n, probs = 0.75),
        UCI = quantile(n, probs = 0.975),
        .groups = "drop"
    ) %>%
    rename(t = time)
    
## save output
saveRDS(runs, paste0("wave", wave, "/plotAgg_T", t, ".rds"))

