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
runs <- map(1:nrow(pars), function(i, time, wave) {
        readRDS(paste0("wave", wave, "/plotSum_", i, "_natFull.rds")) %>%
            filter(t == time) %>%
            select(n, age, t, var)
    }, time = t, wave = wave) %>%
    bind_rows() %>%
    group_by(age, t, var) %>%
    summarise(
        LCI = quantile(n, probs = 0.025),
        LQ = quantile(n, probs = 0.25),
        Median = quantile(n, probs = 0.5),
        UQ = quantile(n, probs = 0.75),
        UCI = quantile(n, probs = 0.975),
        .groups = "drop"
    )
    
## save output
saveRDS(runs, paste0("wave", wave, "/plotAgg_T", t, "_natFull.rds"))

## concatenate runs over ensemble
runs <- map(1:nrow(pars), function(i, time, wave) {
        readRDS(paste0("wave", wave, "/plotSum_", i, "_natDeaths.rds")) %>%
            filter(t == time) %>%
            select(t, n)
    }, time = t, wave = wave) %>%
    bind_rows() %>%
    group_by(t) %>%
    summarise(
        LCI = quantile(n, probs = 0.025),
        LQ = quantile(n, probs = 0.25),
        Median = quantile(n, probs = 0.5),
        UQ = quantile(n, probs = 0.75),
        UCI = quantile(n, probs = 0.975),
        .groups = "drop"
    )
    
## save output
saveRDS(runs, paste0("wave", wave, "/plotAgg_T", t, "_natDeaths.rds"))

## concatenate runs over ensemble
runs <- map(1:nrow(pars), function(i, time, wave) {
        readRDS(paste0("wave", wave, "/plotSum_", i, "_ageRegionDeaths.rds")) %>%
            filter(t == time) %>%
            select(!particle)
    }, time = t, wave = wave) %>%
    bind_rows() %>%
    group_by(t, age, RGN19NM) %>%
    summarise(
        LCI = quantile(n, probs = 0.025),
        LQ = quantile(n, probs = 0.25),
        Median = quantile(n, probs = 0.5),
        UQ = quantile(n, probs = 0.75),
        UCI = quantile(n, probs = 0.975),
        .groups = "drop"
    )
    
## save output
saveRDS(runs, paste0("wave", wave, "/plotAgg_T", t, "_ageRegionDeaths.rds"))

## concatenate runs over ensemble
runs <- map(1:nrow(pars), function(i, time, wave) {
        readRDS(paste0("wave", wave, "/plotSum_", i, "_nhsregionHosp.rds")) %>%
            filter(t == time) %>%
            select(!particle)
    }, time = t, wave = wave) %>%
    bind_rows() %>%
    group_by(t, areaName) %>%
    summarise(
        LCI = quantile(n, probs = 0.025),
        LQ = quantile(n, probs = 0.25),
        Median = quantile(n, probs = 0.5),
        UQ = quantile(n, probs = 0.75),
        UCI = quantile(n, probs = 0.975),
        .groups = "drop"
    )
    
## save output
saveRDS(runs, paste0("wave", wave, "/plotAgg_T", t, "_nhsregionHosp.rds"))

## concatenate runs over ensemble
runs <- map(1:nrow(pars), function(i, time, wave) {
        readRDS(paste0("wave", wave, "/plotSum_", i, "_ageNhsregionHosp.rds")) %>%
            filter(t == time) %>%
            select(!particle)
    }, time = t, wave = wave) %>%
    bind_rows() %>%
    group_by(t, age, areaName) %>%
    summarise(
        LCI = quantile(n, probs = 0.025),
        LQ = quantile(n, probs = 0.25),
        Median = quantile(n, probs = 0.5),
        UQ = quantile(n, probs = 0.75),
        UCI = quantile(n, probs = 0.975),
        .groups = "drop"
    )
    
## save output
saveRDS(runs, paste0("wave", wave, "/plotAgg_T", t, "_ageNhsregionHosp.rds"))

## produce lad-level plots if required
if(file.exists("JASMINcode/lads.txt")) {
    ## concatenate runs over ensemble
    runs <- map(1:nrow(pars), function(i, time, wave) {
            readRDS(paste0("wave", wave, "/plotSum_", i, "_lads.rds")) %>%
                filter(t == time) %>%
                select(!particle)
        }, time = t, wave = wave) %>%
        bind_rows() %>%
        group_by(t, lad, var, age) %>%
        summarise(
            LCI = quantile(n, probs = 0.025),
            LQ = quantile(n, probs = 0.25),
            Median = quantile(n, probs = 0.5),
            UQ = quantile(n, probs = 0.75),
            UCI = quantile(n, probs = 0.975),
            .groups = "drop"
        )
    
    ## save output
    saveRDS(runs, paste0("wave", wave, "/plotAgg_T", t, "_lads.rds"))
}

