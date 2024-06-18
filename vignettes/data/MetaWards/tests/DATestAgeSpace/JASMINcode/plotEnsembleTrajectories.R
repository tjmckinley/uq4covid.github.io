## load libraries
library(tidyverse)
library(patchwork)
library(sf)
library(gganimate)

## set theme
theme_set(theme_minimal())

## check if being run in batch mode
args <- commandArgs(TRUE)
if(length(args) != 0) {
    ## extract command line arguments
    args <- commandArgs(TRUE)
    if(length(args) > 0) {
        stopifnot(length(args) == 5)
        wave <- args[1]
        outputs <- args[2]
        tstart <- as.numeric(args[3])
        tstop <- as.numeric(args[4])
        simData <- as.logical(args[5])
    } else {
        stop("No arguments")
    }
} else {
    ## set wave number and hash
    wave <- "1"
    outputs <- "outputs"
    tstart <- NA
    tstop <- NA
    simData <- FALSE
}
cont <- ifelse(is.na(tstart), FALSE, TRUE)

## set up facet labeller vector
age_label <- as.character(1:8)
age_label <- c("<5", "5-17", 
    "18-29", "30-39", "40-49", "50-59", 
    "60-69", "70+")
names(age_label) <- as.character(1:8)
age_nhs_label <- c("<5", "5-17",
    "18-59", "60+")
names(age_nhs_label) <- as.character(1:4)

###############################################
#######        LAD-level truth          #######
###############################################

## load in runs
sims_md <- readRDS(paste0("../wave", wave, "/sumEns_natFull.rds"))
if(is.na(tstop)) tstop <- max(sims_md$t)

## load in data
if(simData) {
    data <- readRDS(paste0("../", outputs, "/disSims.rds")) %>%
        filter(t <= tstop) %>%
        pivot_longer(!t, names_to = "var", values_to = "n") %>%
        mutate(age = gsub('^(?:[^_]*_)(.*)', '\\1', var)) %>%
        mutate(LAD = gsub('^(?:[^_]*_)(.*)', '\\1', age)) %>%
        mutate(age = gsub('(.*)_[0-9]*', '\\1', age)) %>%
        mutate(var = gsub('^(.*)_[0-9]*_.*', '\\1', var)) %>%
        mutate(var = gsub("one", "1", var)) %>%
        mutate(var = gsub("two", "2", var))
}
    
p1 <- list()
p1[[1]] <- ggplot(sims_md, aes(x = t)) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5) +
    geom_ribbon(aes(ymin = LQ, ymax = UQ), alpha = 0.5) +
    geom_line(aes(y = Median)) +
    facet_grid(var ~ age, scales = "free", labeller = labeller(age = age_label)) +
    xlab("Days") + 
    ylab("Counts") +
    ggtitle("Hidden states")

if(simData) {
    p1[[1]] <- p1[[1]] +
    geom_line(
        aes(y = n), 
        data = group_by(data, t, var, age) %>%
            summarise(n = sum(n), .groups = "drop"),
        col = "red", linetype = "dashed"
    )
}
    
if(cont) p1[[1]] <- p1[[1]] + geom_vline(xintercept = tstart, linetype = "dashed")
        
###############################################
#######     LAD-level observations      #######
###############################################
        
## collapse data for plotting
data <- readRDS(paste0("../", outputs, "/cumDeath_lad.rds")) %>%
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
        col = "blue", linetype = "dashed"
    ) +
    xlab("Days") + 
    ylab("Counts") +
    ggtitle("Observed cumulative deaths (aggregated over LTLAs)")
    
if(cont) p1[[2]] <- p1[[2]] + geom_vline(xintercept = tstart, linetype = "dashed")
              
###############################################
#######  age/region-level observations  #######
###############################################

## load lookup for labels
region_lookup <- readRDS("../data/region_lookup.rds")

## collapse data for plotting
data <- readRDS(paste0("../", outputs, "/cumDeath_age_region.rds")) %>%
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
        col = "blue", linetype = "dashed"
    ) +
    facet_grid(RGN19NM ~ age,
        labeller = labeller(RGN19NM = label_wrap_gen(width = 10), age = age_label)) +
    xlab("Days") + 
    ylab("Counts") +
    ggtitle("Observed cumulative deaths (age / region)")
    
if(cont) p1[[3]] <- p1[[3]] + geom_vline(xintercept = tstart, linetype = "dashed")
        
###############################################
#######     NHS region observations     #######
###############################################

## load lookup for labels
nhsregion_lookup <- readRDS("../data/nhsregion_lookup.rds")
        
## collapse data for plotting
data <- readRDS(paste0("../", outputs, "/hosp_nhsregion.rds")) %>%
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
        col = "blue", linetype = "dashed"
    ) +
    facet_wrap(~ areaName, nrow = 1, labeller = label_wrap_gen(width = 10)) +
    xlab("Days") + 
    ylab("Counts") +
    ggtitle("Observed hospital cases (NHS region)")
    
if(cont) p1[[4]] <- p1[[4]] + geom_vline(xintercept = tstart, linetype = "dashed")
        
###############################################
#####  NHS age/region-level observations  #####
###############################################

## collapse data for plotting
data <- readRDS(paste0("../", outputs, "/cumHospAd_age_nhsregion.rds")) %>%
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
        col = "blue", linetype = "dashed"
    ) +
    facet_grid(areaName ~ age, 
        labeller = labeller(areaName = label_wrap_gen(width = 10), age = age_nhs_label)) +
    xlab("Days") + 
    ylab("Counts") +
    ggtitle("Observed cumulative hospital incidence (NHS age / region)")
    
if(cont) p1[[5]] <- p1[[5]] + geom_vline(xintercept = tstart, linetype = "dashed")
        
###############################################
#####   combine plots and save outputs    #####
###############################################

## save outputs
saveRDS(p1, paste0("../wave", wave, "/plots.rds"))

## combine plots
p1[[4]] <- p1[[4]] / p1[[2]]
p1 <- p1[-2]
p1 <- wrap_plots(p1, nrow = 2, heights = c(0.8, 0.5))
ggsave(paste0("../wave", wave, "/simsBPFEns.pdf"), p1, width = 18, height = 18)

###############################################
#####           spatial plots             #####
###############################################

## load in shapefile
lad19 <- st_read(paste0("../", outputs, "/Local_Authority_Districts_(December_2019)_Boundaries_UK_BUC.shp"))

## load in lookup
death_lookup <- readRDS(paste0("../", outputs, "/death_lookup.rds"))

## load in data
data <- readRDS(paste0("../", outputs, "/cumDeath_lad.rds")) %>%
    filter(t <= tstop) %>%
    pivot_longer(!t, names_to = "lad", values_to = "n") %>%
    mutate(lad = as.numeric(gsub("deaths_", "", lad))) %>%
    inner_join(death_lookup, by = c("lad" = "FID")) %>%
    rename(Data = n)

## load in runs
sims_md <- readRDS(paste0("../wave", wave, "/sumEns_lads.rds")) %>%
    dplyr::select(t, lad, Median) %>%
    rename(Prediction = Median)

## join runs and data
data <- inner_join(data, sims_md, by = c("lad", "t"))

## join to shapefile
lad19 <- inner_join(lad19, data, by = c("lad19cd" = "areaCode"))

### animation
#p <- pivot_longer(lad19, c(Data, Prediction), names_to = "type", values_to = "Count") %>%
#    ggplot() +
#        geom_sf(aes(fill = Count), colour = NA) +
#        facet_wrap(~ type) +
#        scale_fill_viridis_c() +
#        transition_time(t) +
#        ggtitle("Deaths at t = {frame_time}")
#anim_save(paste0("../wave", wave, "/simsspanimation.gif"), p)

## static plot
p <- filter(lad19, t == max(t)) %>%
    pivot_longer(c(Data, Prediction), names_to = "type", values_to = "Count") %>%
    ggplot() +
        geom_sf(aes(fill = Count), colour = NA) +
        facet_wrap(~ type) +
        scale_fill_viridis_c() +
        ggtitle(paste0("Deaths at t = ", max(lad19$t)))
ggsave(paste0("../wave", wave, "/simsspstatic.pdf"), p)

## save output
p1 <- p

###############################################
#####             LTLA plots              #####
###############################################

## load in data
data <- readRDS(paste0("../", outputs, "/cumDeath_lad.rds")) %>%
    filter(t <= tstop) %>%
    pivot_longer(!t, names_to = "lad", values_to = "n") %>%
    mutate(lad = as.numeric(gsub("deaths_", "", lad))) %>%
    inner_join(death_lookup, by = c("lad" = "FID")) %>%
    rename(Data = n) %>%
    select(t, lad, Data)

## load in runs
sims_md <- readRDS(paste0("../wave", wave, "/sumEns_lads.rds")) %>%
    dplyr::select(t, lad, Median, LCI, UCI) %>%
    rename(Prediction = Median)

## join runs and data
data <- inner_join(data, sims_md, by = c("lad", "t"))

## produce plot
p <- filter(data, t == max(t)) %>%
    arrange(desc(Prediction)) %>%
    mutate(lad = 1:n()) %>%
    mutate(inside = ifelse(Data >= LCI & Data <= UCI, "Inside 95% CI", "Outside 95% CI")) %>%
    pivot_longer(c(Data, Prediction), names_to = "type", values_to = "Count") %>%
    mutate(inside = ifelse(type == "Prediction", "Prediction", inside)) %>%
    mutate(inside = factor(inside, levels = c("Inside 95% CI", "Outside 95% CI", "Prediction"))) %>%
    arrange(desc(inside)) %>%
    ggplot() +
        geom_errorbar(aes(x = lad, ymin = LCI, ymax = UCI), colour = "#52854C") +
        geom_point(aes(x = lad, y = Count, colour = inside)) +
        ggtitle(paste0("Cumulative deaths at t = ", max(data$t))) +
        xlab("LTLA (in decreasing order of deaths)") +
        ylab("Cumulative Deaths") +
        scale_colour_manual(name = "", values = c("#E69F00", "#D55E00", "#52854C")) +
        theme(axis.text.x = element_blank(), legend.position = "bottom")
ggsave(paste0("../wave", wave, "/simsLTLA.pdf"), p, width = 15, height = 5)

## save outputs
p1 <- list(p1, p)
saveRDS(p1, paste0("../wave", wave, "/plots_sp.rds"))

###############################################
#####          LAD-level plots            #####
###############################################

if(file.exists(paste0("../", outputs, "/lads_", outputs, ".txt"))) {

    ## read in lads
    lads <- as.numeric(readLines(paste0("../", outputs, "/lads_", outputs, ".txt")))
    
    ## load in data
    data <- readRDS(paste0("../", outputs, "/cumDeath_lad.rds")) %>%
        filter(t <= tstop) %>%
        pivot_longer(!t, names_to = "lad", values_to = "n") %>%
        mutate(lad = gsub("deaths_", "", lad)) %>%
        filter(lad %in% lads)

    ## load in runs
    sims_md <- readRDS(paste0("../wave", wave, "/sumEns_lads.rds")) %>%
        filter(lad %in% lads)
        
    p1 <- ggplot(sims_md, aes(x = t)) +
        geom_ribbon(aes(ymin = LCI, ymax = UCI), colour = NA, alpha = 0.5) +
        geom_ribbon(aes(ymin = LQ, ymax = UQ), colour = NA, alpha = 0.5) +
        geom_line(aes(y = Median)) +
        geom_line(
            aes(y = n), 
            data = data,
            linetype = "dashed",
            col = "blue"
        ) +
        facet_wrap(~ lad, scales = "free") +
        xlab("Days") + 
        ylab("Counts") +
        ggtitle(paste0("Observed deaths in top ", length(lads), " LTLAs"))

    if(cont) p1 <- p1 + geom_vline(xintercept = tstart, linetype = "dashed")

    ## save plot
    ggsave(paste0("../wave", wave, "/simsTopLADsBPFEns.pdf"), p1, width = 10, height = 10)
    
    ## save outputs
    saveRDS(p1, paste0("../wave", wave, "/plots_lad.rds"))
}

