## load libraries
library(tidyverse)

## check if being run in batch mode
args <- commandArgs(TRUE)
if(length(args) != 0) {
    ## extract command line arguments
    args <- commandArgs(TRUE)
    if(length(args) > 0) {
        stopifnot(length(args) == 3)
        wave <- args[1]
	    outputs <- args[2]
        t <- as.numeric(args[3])
    } else {
        stop("No arguments")
    }
} else {
    ## set name of directory to search for outcomes
    wave <- "1"
    outputs <- "outputs"
    t <- 1
}

###############################################
#######        LAD-level truth          #######
###############################################

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

###############################################
#######     LAD-level observations      #######
###############################################

## collapse data for plotting
data <- readRDS(paste0(outputs, "/cumDeath_lad.rds"))
data <- data[apply(select(data, !t), 1, function(x) any(!is.na(x))), ]
mint <- min(data$t)

## concatenate runs over ensemble
if(mint <= t) {
    runs <- map(1:nrow(pars), function(i, time, wave, mint) {
        if(mint == 0) {
            runs <- readRDS(paste0("wave", wave, "/plotSum_", i, "_natDeaths.rds")) %>%
                filter(t == time) %>%
                select(t, n)
        } else {
            runs <- readRDS(paste0("wave", wave, "/plotSum_", i, "_natDeaths.rds")) %>%
                filter(t <= time & t >= mint) %>%
                group_by(particle) %>%
                mutate(n = n - min(n)) %>%
                ungroup() %>%
                filter(t == time) %>%
                select(t, n)
        }
        runs
    }, time = t, wave = wave, mint = mint) %>%
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
} else {
    runs <- tibble(t = t, LCI = NA, LQ = NA, Median = NA, UQ = NA, UCI = NA)
}
        
## save output
saveRDS(runs, paste0("wave", wave, "/plotAgg_T", t, "_natDeaths.rds"))
              
###############################################
#######  age/region-level observations  #######
###############################################

## collapse data for plotting
data <- readRDS(paste0(outputs, "/cumDeath_age_region.rds"))  
data <- data[apply(select(data, !t), 1, function(x) any(!is.na(x))), ]
mint <- min(data$t)

if(mint <= t) {
    ## concatenate runs over ensemble
    runs <- map(1:nrow(pars), function(i, time, wave, mint) {
        if(mint == 0) {
            runs <- readRDS(paste0("wave", wave, "/plotSum_", i, "_ageRegionDeaths.rds")) %>%
                filter(t == time) %>%
                select(!particle)
        } else {
            runs <- readRDS(paste0("wave", wave, "/plotSum_", i, "_ageRegionDeaths.rds")) %>%
                filter(t <= time & t >= mint) %>%
                group_by(particle, age, RGN19NM) %>%
                mutate(n = n - min(n)) %>%
                ungroup() %>%
                filter(t == time) %>%
                select(!particle)
        }
        runs
    }, time = t, wave = wave, mint = mint) %>%
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
} else {
    runs <- readRDS(paste0("wave", wave, "/plotSum_1_ageRegionDeaths.rds")) %>%
        filter(t == 0) %>%
        select(!c(particle, n)) %>%
        distinct() %>%
        arrange(age, RGN19NM) %>%
        mutate(LCI = NA, LQ = NA, Median = NA, UQ = NA, UCI = NA) %>%
        mutate(across(!c(t, age, RGN19NM), as.numeric))
    runs$t <- t
}
        
## save output
saveRDS(runs, paste0("wave", wave, "/plotAgg_T", t, "_ageRegionDeaths.rds"))
        
###############################################
#######     NHS region observations     #######
###############################################

## collapse data for plotting
data <- readRDS(paste0(outputs, "/hosp_nhsregion.rds"))
data <- data[apply(select(data, !t), 1, function(x) any(!is.na(x))), ]
mint <- min(data$t)
 
if(mint <= t) {
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
} else {
    runs <- readRDS(paste0("wave", wave, "/plotSum_1_nhsregionHosp.rds")) %>%
        filter(t == 0) %>%
        select(!c(particle, n)) %>%
        distinct() %>%
        arrange(areaName) %>%
        mutate(LCI = NA, LQ = NA, Median = NA, UQ = NA, UCI = NA) %>%
        mutate(across(!c(t, areaName), as.numeric))
    runs$t <- t
}
        
## save output
saveRDS(runs, paste0("wave", wave, "/plotAgg_T", t, "_nhsregionHosp.rds"))
        
###############################################
#####  NHS age/region-level observations  #####
###############################################

## collapse data for plotting
data <- readRDS(paste0(outputs, "/cumHospAd_age_nhsregion.rds"))
data <- data[apply(select(data, !t), 1, function(x) any(!is.na(x))), ]
mint <- min(data$t)
 
if(mint <= t) {
    ## concatenate runs over ensemble
    runs <- map(1:nrow(pars), function(i, time, wave, mint) {
            if(mint == 0) {
                runs <- readRDS(paste0("wave", wave, "/plotSum_", i, "_ageNhsregionHosp.rds")) %>%
                    filter(t == time) %>%
                    select(!particle)
            } else  {
                runs <- readRDS(paste0("wave", wave, "/plotSum_", i, "_ageNhsregionHosp.rds")) %>%
                    filter(t <= time & t >= mint) %>%
                    group_by(particle, age, areaName) %>%
                    mutate(n = n - min(n)) %>%
                    ungroup() %>%
                    filter(t == time) %>%
                    select(!particle)
            } 
        }, time = t, wave = wave, mint = mint) %>%
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
} else {
    runs <- readRDS(paste0("wave", wave, "/plotSum_1_ageNhsregionHosp.rds")) %>%
        filter(t == 0) %>%
        select(!c(particle, n)) %>%
        distinct() %>%
        arrange(age, areaName) %>%
        mutate(LCI = NA, LQ = NA, Median = NA, UQ = NA, UCI = NA) %>%
        mutate(across(!c(t, areaName), as.numeric))
    runs$t <- t
}
        
## save output
saveRDS(runs, paste0("wave", wave, "/plotAgg_T", t, "_ageNhsregionHosp.rds"))

## produce lad-level plots if required
if(file.exists(paste0(outputs, "/lads_", outputs, ".txt"))) {
    data <- readRDS(paste0(outputs, "/cumDeath_lad.rds"))
    data <- data[apply(select(data, !t), 1, function(x) any(!is.na(x))), ]
    mint <- min(data$t)

    ## concatenate runs over ensemble
    if(mint <= t) {
        runs <- map(1:nrow(pars), function(i, time, wave, mint) {
            if(mint == 0) {
                runs <- readRDS(paste0("wave", wave, "/plotSum_", i, "_lads.rds")) %>%
                    filter(t == time) %>%
                    select(t, n, lad)
            } else {
                runs <- readRDS(paste0("wave", wave, "/plotSum_", i, "_lads.rds")) %>%
                    filter(t <= time & t >= mint) %>%
                    group_by(particle) %>%
                    mutate(n = n - min(n)) %>%
                    ungroup() %>%
                    filter(t == time) %>%
                    select(t, n, lad)
            }
	    runs
        }, time = t, wave = wave, mint = mint) %>%
        bind_rows() %>%
        group_by(t, lad) %>%
        summarise(
	    LCI = quantile(n, probs = 0.025),
            LQ = quantile(n, probs = 0.25),
            Median = quantile(n, probs = 0.5),
            UQ = quantile(n, probs = 0.75),
            UCI = quantile(n, probs = 0.975),
            .groups = "drop"
        )
    } else {
        runs <- tibble(t = t, LCI = NA, LQ = NA, Median = NA, UQ = NA, UCI = NA)
    }

    ## save output
    saveRDS(runs, paste0("wave", wave, "/plotAgg_T", t, "_lads.rds"))
}

