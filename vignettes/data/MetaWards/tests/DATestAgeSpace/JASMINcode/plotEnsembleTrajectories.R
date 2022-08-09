## load libraries
library(tidyverse)
library(patchwork)

## check if being run in batch mode
args <- commandArgs(TRUE)
if(length(args) != 0) {
    ## extract command line arguments
    args <- commandArgs(TRUE)
    if(length(args) > 0) {
        stopifnot(length(args) == 3)
        wave <- as.numeric(args[1])
        tstart <- as.numeric(args[2])
        tstop <- as.numeric(args[3])
    } else {
        stop("No arguments")
    }
} else {
    ## set wave number and hash
    wave <- 1
    tstart <- NA
    tstop <- NA
}
cont <- ifelse(is.na(tstart), FALSE, TRUE)

###############################################
#######        LAD-level truth          #######
###############################################

## load in runs
sims_md <- readRDS(paste0("../wave", wave, "/sumEns_natFull.rds"))
if(is.na(tstop)) tstop <- max(sims_md$t)

## load in data
data <- readRDS("../outputs/disSims.rds") %>%
    filter(t <= tstop) %>%
    pivot_longer(!t, names_to = "var", values_to = "n") %>%
    mutate(age = gsub('^(?:[^_]*_)(.*)', '\\1', var)) %>%
    mutate(LAD = gsub('^(?:[^_]*_)(.*)', '\\1', age)) %>%
    mutate(age = gsub('(.*)_[0-9]*', '\\1', age)) %>%
    mutate(var = gsub('^(.*)_[0-9]*_.*', '\\1', var)) %>%
    mutate(var = gsub("one", "1", var)) %>%
    mutate(var = gsub("two", "2", var))
    
p1 <- list()     
p1[[1]] <- ggplot(sims_md, aes(x = t)) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5) +
    geom_ribbon(aes(ymin = LQ, ymax = UQ), alpha = 0.5) +
    geom_line(aes(y = Median)) +
    geom_line(
        aes(y = n), 
        data = group_by(data, t, var, age) %>%
            summarise(n = sum(n), .groups = "drop"),
        col = "red", linetype = "dashed"
    ) +
    facet_grid(var ~ age, scales = "free") +
    xlab("Days") + 
    ylab("Counts") +
    ggtitle("Truth")
    
if(cont) p1[[1]] <- p1[[1]] + geom_vline(xintercept = tstart, linetype = "dashed", colour = "blue")
        
###############################################
#######     LAD-level observations      #######
###############################################
        
## collapse data for plotting
data <- readRDS("../outputs/cumDeath_lad.rds") %>%
    filter(t <= tstop) %>%
    pivot_longer(!t, names_to = "lad", values_to = "n") %>%
    mutate(lad = as.numeric(gsub("deaths_", "", lad))) %>%
    group_by(t) %>%
    summarise(n = sum(n), .groups = "drop")

## load in runs
sims_md <- readRDS(paste0("../wave", wave, "/sumEns_natDeaths.rds"))
     
p1[[2]] <- ggplot(sims_md, aes(x = t)) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5) +
    geom_ribbon(aes(ymin = LQ, ymax = UQ), alpha = 0.5) +
    geom_line(aes(y = Median)) +
    geom_line(
        aes(y = n), 
        data = data,
        col = "red", linetype = "dashed"
    ) +
    xlab("Days") + 
    ylab("Counts") +
    ggtitle("Observed cumulative deaths (aggregated over LADs)")
    
if(cont) p1[[2]] <- p1[[2]] + geom_vline(xintercept = tstart, linetype = "dashed", colour = "blue")
              
###############################################
#######  age/region-level observations  #######
###############################################

## load lookup for labels
region_lookup <- readRDS("../data/region_lookup.rds")

## collapse data for plotting
data <- readRDS("../outputs/cumDeath_age_region.rds") %>%
    filter(t <= tstop) %>%
    pivot_longer(!t, names_to = "var", values_to = "n") %>%
    mutate(age = gsub('^(?:[^_]*_)(.*)', '\\1', var)) %>%
    mutate(region = gsub('^(?:[^_]*_)(.*)', '\\1', age)) %>%
    mutate(age = gsub('(.*)_[0-9]*', '\\1', age)) %>%
    mutate(across(c(age, region), as.numeric)) %>%
    inner_join(region_lookup, by = c("region" = "FID")) %>%
    select(t, n, age, RGN19NM)
    
## load in runs
sims_md <- readRDS(paste0("../wave", wave, "/sumEns_ageRegionDeaths.rds"))
    
p1[[3]] <- ggplot(sims_md, aes(x = t)) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5) +
    geom_ribbon(aes(ymin = LQ, ymax = UQ), alpha = 0.5) +
    geom_line(aes(y = Median)) +
    geom_line(
        aes(y = n), 
        data = data,
        col = "red", linetype = "dashed"
    ) +
    facet_grid(RGN19NM ~ age, labeller = label_wrap_gen(width = 10)) +
    xlab("Days") + 
    ylab("Counts") +
    ggtitle("Observed cumulative deaths (age / region)")
    
if(cont) p1[[3]] <- p1[[3]] + geom_vline(xintercept = tstart, linetype = "dashed", colour = "blue")
        
###############################################
#######     NHS region observations     #######
###############################################

## load lookup for labels
nhsregion_lookup <- readRDS("../data/nhsregion_lookup.rds")
        
## collapse data for plotting
data <- readRDS("../outputs/hosp_nhsregion.rds") %>%
    filter(t <= tstop) %>%
    pivot_longer(!t, names_to = "region", values_to = "n") %>%
    mutate(region = as.numeric(gsub("hosp_", "", region))) %>%
    inner_join(nhsregion_lookup, by = c("region" = "FID")) %>%
    select(t, n, areaName)
    
## load in runs
sims_md <- readRDS(paste0("../wave", wave, "/sumEns_nhsregionHosp.rds"))
     
p1[[4]] <- ggplot(sims_md, aes(x = t)) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5) +
    geom_ribbon(aes(ymin = LQ, ymax = UQ), alpha = 0.5) +
    geom_line(aes(y = Median)) +
    geom_line(
        aes(y = n), 
        data = data,
        col = "red", linetype = "dashed"
    ) +
    facet_wrap(~ areaName, nrow = 1, labeller = label_wrap_gen(width = 10)) +
    xlab("Days") + 
    ylab("Counts") +
    ggtitle("Observed hospital cases (NHS region)")
    
if(cont) p1[[4]] <- p1[[4]] + geom_vline(xintercept = tstart, linetype = "dashed", colour = "blue")
        
###############################################
#####  NHS age/region-level observations  #####
###############################################

## collapse data for plotting
data <- readRDS("../outputs/cumHospAd_age_nhsregion.rds") %>%
    filter(t <= tstop) %>%
    pivot_longer(!t, names_to = "var", values_to = "n") %>%
    mutate(age = gsub('^(?:[^_]*_)(.*)', '\\1', var)) %>%
    mutate(region = gsub('^(?:[^_]*_)(.*)', '\\1', age)) %>%
    mutate(age = gsub('(.*)_[0-9]*', '\\1', age)) %>%
    mutate(across(c(age, region), as.numeric)) %>%
    inner_join(nhsregion_lookup, by = c("region" = "FID")) %>%
    select(t, n, age, areaName)

## load in runs
sims_md <- readRDS(paste0("../wave", wave, "/sumEns_ageNhsregionHosp.rds"))
     
p1[[5]] <- ggplot(sims_md, aes(x = t)) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5) +
    geom_ribbon(aes(ymin = LQ, ymax = UQ), alpha = 0.5) +
    geom_line(aes(y = Median)) +
    geom_line(
        aes(y = n), 
        data = data,
        col = "red", linetype = "dashed"
    ) +
    facet_grid(areaName ~ age, labeller = label_wrap_gen(width = 10)) +
    xlab("Days") + 
    ylab("Counts") +
    ggtitle("Observed cumulative hospital incidence (NHS age / region)")
    
if(cont) p1[[5]] <- p1[[5]] + geom_vline(xintercept = tstart, linetype = "dashed", colour = "blue")
        
###############################################
#####   combine plots and save outputs    #####
###############################################

## combine plots
p1[[4]] <- p1[[4]] / p1[[2]]
p1 <- p1[-2]
p1 <- wrap_plots(p1, nrow = 2, heights = c(0.8, 0.5))
ggsave(paste0("../wave", wave, "/simsBPFEns.pdf"), p1, width = 15, height = 15)

###############################################
#####          LAD-level plots            #####
###############################################

if(file.exists("lads.txt")) {

    ## read in lads
    lads <- as.numeric(readLines("lads.txt"))
    
    ## load in data
    data <- readRDS("../outputs/disSims.rds") %>%
        filter(t <= tstop) %>%
        pivot_longer(!t, names_to = "var", values_to = "n") %>%
        mutate(age = gsub('^(?:[^_]*_)(.*)', '\\1', var)) %>%
        mutate(lad = gsub('^(?:[^_]*_)(.*)', '\\1', age)) %>%
        mutate(age = gsub('(.*)_[0-9]*', '\\1', age)) %>%
        mutate(var = gsub('^(.*)_[0-9]*_.*', '\\1', var)) %>%
        mutate(var = gsub("one", "1", var)) %>%
        mutate(var = gsub("two", "2", var)) %>%
        mutate(age = as.numeric(age)) %>%
        filter(lad %in% lads) 

    ## load in runs
    sims_md <- readRDS(paste0("../wave", wave, "/sumEns_lads.rds"))
        
    p1 <- ggplot(sims_md, aes(x = t, colour = lad, fill = lad)) +
        geom_ribbon(aes(ymin = LCI, ymax = UCI), colour = NA, alpha = 0.5) +
        geom_ribbon(aes(ymin = LQ, ymax = UQ), colour = NA, alpha = 0.5) +
        geom_line(aes(y = Median)) +
        geom_line(
            aes(y = n), 
            data = data,
            linetype = "dashed"
        ) +
        facet_grid(var ~ age, scales = "free") +
        xlab("Days") + 
        ylab("Counts") +
        ggtitle("Truth")

    if(cont) p1 <- p1 + geom_vline(xintercept = tstart, linetype = "dashed", colour = "blue")

    ggsave(paste0("../wave", wave, "/simsTopLADsBPFEns.pdf"), p1, width = 10, height = 10)
}

