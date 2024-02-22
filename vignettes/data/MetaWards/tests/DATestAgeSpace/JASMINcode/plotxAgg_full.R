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
        t <- as.numeric(args[3])
        if(length(args) > 3) {
            stopifnot(length(args) == 4)
            ind <- args[4]
        } else {
            ind <- NA
        }
    } else {
        stop("No arguments")
    }
} else {
    ## set name of directory to search for outcomes
    wave <- "1"
    outputs <- "outputs"
    t <- 1
    ind <- NA
}

###############################################
#######        LAD-level truth          #######
###############################################

## read in input file
if(is.na(ind)) {
    pars <- readRDS(paste0("wave", wave, "/disease.rds"))
    pars <- 1:nrow(pars)
    pars <- paste0("wave", wave, "/plotSum_", pars)
} else {
    pars <- read_csv(paste0("wave", wave, "/", ind))
    stopifnot(identical(colnames(pars), c("wave", "ind")))
    pars <- paste0("wave", pars$wave, "/plotSum_", pars$ind)
}
print(pars)

## concatenate runs over ensemble
runs <- map(pars, function(par, time) {
        readRDS(paste0(par, "_natFull.rds")) %>%
            filter(t == time) %>%
            select(n, age, t, var)
    }, time = t) %>%
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
data <- readRDS(paste0(outputs, "/cumDI_age_lad.rds"))
data <- data[apply(select(data, !t), 1, function(x) any(!is.na(x))), ]
mint <- min(data$t)

## concatenate runs over ensemble
if(mint <= t) {
    runs <- map(pars, function(par, time, mint) {
        if(mint == 0) {
            runs <- readRDS(paste0(par, "_natAgeDeathsHosp.rds")) %>%
                filter(t == time)
        } else {
            runs <- readRDS(paste0(par, "_natAgeDeathsHosp.rds")) %>%
                filter(t <= time & t >= mint) %>%
                group_by(particle, age) %>%
                mutate(DI = DI - min(DI)) %>%
                mutate(DH = DH - min(DH)) %>%
                mutate(Hcum = Hcum - min(Hcum)) %>%
                ungroup() %>%
                filter(t == time)
        }
        runs
    }, time = t, mint = mint) %>%
    bind_rows() %>%
    group_by(t, age) %>%
    summarise(
        across(c(DI, DH, H, Hcum), list(
            LCI =~quantile(., probs = 0.025),
            LQ = ~quantile(., probs = 0.25),
            Median = ~quantile(., probs = 0.5),
            UQ = ~quantile(., probs = 0.75),
            UCI = ~quantile(., probs = 0.975)
        )),
        .groups = "drop"
    )
} else {
    runs <- matrix(NA, 1, 22)
    colnames(runs) <- c("t", "age", paste0(rep(c("DI", "DH", "H", "Hcum"), each = 5), "_", c("LCI", "LQ", "Median", "UQ", "UCI")))
    runs <- as_tibble(runs)
}
        
## save output
saveRDS(runs, paste0("wave", wave, "/plotAgg_T", t, "_natAgeDeathsHosp.rds"))

## produce lad-level plots if required
if(file.exists(paste0(outputs, "/lads_", outputs, ".txt"))) {
    data <- readRDS(paste0(outputs, "/cumDI_age_lad.rds"))
    data <- data[apply(select(data, !t), 1, function(x) any(!is.na(x))), ]
    mint <- min(data$t)

    if(mint <= t) {
        runs <- map(pars, function(par, time, mint) {
            if(mint == 0) {
                runs <- readRDS(paste0(par, "_age_lads.rds")) %>%
                    filter(t == time)
            } else {
                runs <- readRDS(paste0(par, "_age_lads.rds")) %>%
                    filter(t <= time & t >= mint) %>%
                    group_by(particle, age, lad) %>%
                    mutate(DI = DI - min(DI)) %>%
                    mutate(DH = DH - min(DH)) %>%
                    mutate(Hcum = Hcum - min(Hcum)) %>%
                    ungroup() %>%
                    filter(t == time)
            }
            runs
        }, time = t, mint = mint) %>%
        bind_rows() %>%
        group_by(t, age, lad) %>%
        summarise(
            across(c(DI, DH, H, Hcum), list(
                LCI =~quantile(., probs = 0.025),
                LQ = ~quantile(., probs = 0.25),
                Median = ~quantile(., probs = 0.5),
                UQ = ~quantile(., probs = 0.75),
                UCI = ~quantile(., probs = 0.975)
            )),
            .groups = "drop"
        )
    } else {
        runs <- matrix(NA, 1, 23)
        colnames(runs) <- c("t", "age", "lad", paste0(rep(c("DI", "DH", "H", "Hcum"), each = 5), "_", c("LCI", "LQ", "Median", "UQ", "UCI")))
        runs <- as_tibble(runs)
    }

    ## save output
    saveRDS(runs, paste0("wave", wave, "/plotAgg_T", t, "_age_lads.rds"))
}

