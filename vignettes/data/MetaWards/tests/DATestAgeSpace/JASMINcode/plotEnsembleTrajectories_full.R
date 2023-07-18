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
        stopifnot(length(args) == 4)
        wave <- args[1]
	    outputs <- args[2]
        tstart <- as.numeric(args[3])
        tstop <- as.numeric(args[4])
    } else {
        stop("No arguments")
    }
} else {
    ## set wave number and hash
    wave <- "1"
    outputs <- "outputs"
    tstart <- NA
    tstop <- NA
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
data <- readRDS(paste0("../", outputs, "/disSims.rds")) %>%
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
    facet_grid(var ~ age, scales = "free", labeller = labeller(age = age_label)) +
    xlab("Days") + 
    ylab("Counts") +
    ggtitle("Hidden states")
    
if(cont) p1[[1]] <- p1[[1]] + geom_vline(xintercept = tstart, linetype = "dashed")
        
###############################################
#######     LAD-level observations      #######
###############################################
        
## collapse data for plotting
DI <- readRDS(paste0("../", outputs, "/cumDI_age_lad.rds")) %>%
    filter(t <= tstop) %>%
    pivot_longer(!t, values_to = "n") %>%
    mutate(age = gsub('^(?:[^_]*_)(.*)', '\\1', name)) %>%
    mutate(lad = gsub('^(?:[^_]*_)(.*)', '\\1', age)) %>%
    mutate(age = gsub('(.*)_[0-9]*', '\\1', age)) %>%
    mutate(across(c(age, lad), as.numeric)) %>%
    select(!name) %>%
    group_by(age, lad) %>%
    mutate(n = cumsum(n)) %>%
    ungroup()
DH <- readRDS(paste0("../", outputs, "/cumDH_age_lad.rds")) %>%
    filter(t <= tstop) %>%
    pivot_longer(!t, values_to = "n") %>%
    mutate(age = gsub('^(?:[^_]*_)(.*)', '\\1', name)) %>%
    mutate(lad = gsub('^(?:[^_]*_)(.*)', '\\1', age)) %>%
    mutate(age = gsub('(.*)_[0-9]*', '\\1', age)) %>%
    mutate(across(c(age, lad), as.numeric)) %>%
    select(!name) %>%
    group_by(age, lad) %>%
    mutate(n = cumsum(n)) %>%
    ungroup()
H <- readRDS(paste0("../", outputs, "/H_age_lad.rds")) %>%
    filter(t <= tstop) %>%
    pivot_longer(!t, values_to = "n") %>%
    mutate(age = gsub('^(?:[^_]*_)(.*)', '\\1', name)) %>%
    mutate(lad = gsub('^(?:[^_]*_)(.*)', '\\1', age)) %>%
    mutate(age = gsub('(.*)_[0-9]*', '\\1', age)) %>%
    mutate(across(c(age, lad), as.numeric)) %>%
    select(!name)

## load in runs
sims_md <- readRDS(paste0("../wave", wave, "/sumEns_natAgeDeathsHosp.rds")) %>%
    pivot_longer(!c(t, age)) %>%
    mutate(var = gsub('^([^_]*)_(.*)', '\\1', name)) %>%
    mutate(name = gsub('^(?:[^_]*_)(.*)', '\\1', name))

## plot deaths
p1[[2]] <- filter(sims_md, var == "DI") %>%
    select(!var) %>%
    pivot_wider(names_from = name, values_from = value) %>%
    ggplot(aes(x = t)) +
        geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5) +
        geom_ribbon(aes(ymin = LQ, ymax = UQ), alpha = 0.5) +
        geom_line(aes(y = Median)) +
        geom_line(
            aes(y = n), 
            data = group_by(deaths, t, age) %>%
                summarise(n = sum(n), .groups = "drop"),
            col = "blue", linetype = "dashed"
        ) +
        facet_wrap(~age) +
        xlab("Days") + 
        ylab("Counts") +
        ggtitle("Observed cumulative community deaths (aggregated over LTLAs)")
    
if(cont) p1[[2]] <- p1[[2]] + geom_vline(xintercept = tstart, linetype = "dashed")

## plot hospitalisations
p1[[3]] <- filter(sims_md, var == "DH") %>%
    select(!var) %>%
    pivot_wider(names_from = name, values_from = value) %>%
    ggplot(aes(x = t)) +
        geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5) +
        geom_ribbon(aes(ymin = LQ, ymax = UQ), alpha = 0.5) +
        geom_line(aes(y = Median)) +
        geom_line(
            aes(y = n), 
            data = group_by(hosp, t, age) %>%
                summarise(n = sum(n), .groups = "drop"),
            col = "blue", linetype = "dashed"
        ) +
        facet_wrap(~age) +
        xlab("Days") + 
        ylab("Counts") +
        ggtitle("Observed cumulative hospital deaths (aggregated over LTLAs)")
    
if(cont) p1[[3]] <- p1[[3]] + geom_vline(xintercept = tstart, linetype = "dashed")

## plot hospital cases
p1[[4]] <- filter(sims_md, var == "H") %>%
    select(!var) %>%
    pivot_wider(names_from = name, values_from = value) %>%
    ggplot(aes(x = t)) +
        geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5) +
        geom_ribbon(aes(ymin = LQ, ymax = UQ), alpha = 0.5) +
        geom_line(aes(y = Median)) +
        geom_line(
            aes(y = n), 
            data = group_by(H, t, age) %>%
                summarise(n = sum(n), .groups = "drop"),
            col = "blue", linetype = "dashed"
        ) +
        facet_wrap(~age) +
        xlab("Days") + 
        ylab("Counts") +
        ggtitle("Observed hospital cases (aggregated over LTLAs)")
    
if(cont) p1[[4]] <- p1[[4]] + geom_vline(xintercept = tstart, linetype = "dashed")
            
###############################################
#####   combine plots and save outputs    #####
###############################################

## save outputs
saveRDS(p1, paste0("../wave", wave, "/plots.rds"))

## combine plots
p1 <- (p1[[1]] + p1[[2]]) / (p1[[3]] + p1[[4]])
ggsave(paste0("../wave", wave, "/simsBPFEns.pdf"), p1, width = 15, height = 15)

###############################################
#####           spatial plots             #####
###############################################

## load in shapefile
lad19 <- st_read(paste0("../", outputs, "/Local_Authority_Districts_(December_2019)_Boundaries_UK_BUC.shp"))

## load in lookup
death_lookup <- readRDS(paste0("../", outputs, "/death_lookup.rds"))

## load in data
data <- group_by(deaths, t, lad) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    rename(Data = n) %>%
    inner_join(death_lookup, by = c("lad" = "FID"))

## load in runs
sims_md <- readRDS(paste0("../wave", wave, "/sumEns_age_lads.rds")) %>%
    dplyr::select(t, lad, age, deaths_Median) %>%
    rename(Prediction = deaths_Median) %>%
    group_by(t, lad) %>%
    summarise(Prediction = sum(Prediction), .groups = "drop")

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

###############################################
#####          LAD-level plots            #####
###############################################

if(file.exists(paste0("../", outputs, "/lads_", outputs, ".txt"))) {

    ## read in lads
    lads <- as.numeric(readLines(paste0("../", outputs, "/lads_", outputs, ".txt")))
    
    ## load in data
    deaths <- filter(deaths, lad %in% lads)
    hosp <- filter(hosp, lad %in% lads)
    H <- filter(H, lad %in% lads)

    ## load in runs
    sims_md <- readRDS(paste0("../wave", wave, "/sumEns_age_lads.rds")) %>%
        filter(lad %in% lads) %>%
        pivot_longer(!c(t, age, lad)) %>%
        mutate(var = gsub('^([^_]*)_(.*)', '\\1', name)) %>%
        mutate(name = gsub('^(?:[^_]*_)(.*)', '\\1', name)) %>%
        pivot_wider(names_from = name, values_from = value)
        
    p1 <- list()
    p1[[1]] <- filter(sims_md, var == "deaths") %>%
        ggplot(aes(x = t)) +
            geom_ribbon(aes(ymin = LCI, ymax = UCI), colour = NA, alpha = 0.5) +
            geom_ribbon(aes(ymin = LQ, ymax = UQ), colour = NA, alpha = 0.5) +
            geom_line(aes(y = Median)) +
            geom_line(
                aes(y = n), 
                data = deaths,
                linetype = "dashed",
                col = "blue"
            ) +
            facet_grid(lad ~ age, scales = "free") +
            xlab("Days") + 
            ylab("Counts") +
            ggtitle(paste0("Observed deaths in top ", length(lads), " LTLAs"))
    p1[[2]] <- filter(sims_md, var == "hosp") %>%
        ggplot(aes(x = t)) +
            geom_ribbon(aes(ymin = LCI, ymax = UCI), colour = NA, alpha = 0.5) +
            geom_ribbon(aes(ymin = LQ, ymax = UQ), colour = NA, alpha = 0.5) +
            geom_line(aes(y = Median)) +
            geom_line(
                aes(y = n), 
                data = hosp,
                linetype = "dashed",
                col = "blue"
            ) +
            facet_grid(lad ~ age, scales = "free") +
            xlab("Days") + 
            ylab("Counts") +
            ggtitle(paste0("Observed hospitalisations in top ", length(lads), " LTLAs"))
    p1[[3]] <- filter(sims_md, var == "H") %>%
        ggplot(aes(x = t)) +
            geom_ribbon(aes(ymin = LCI, ymax = UCI), colour = NA, alpha = 0.5) +
            geom_ribbon(aes(ymin = LQ, ymax = UQ), colour = NA, alpha = 0.5) +
            geom_line(aes(y = Median)) +
            geom_line(
                aes(y = n), 
                data = H,
                linetype = "dashed",
                col = "blue"
            ) +
            facet_grid(lad ~ age, scales = "free") +
            xlab("Days") + 
            ylab("Counts") +
            ggtitle(paste0("Observed hospital cases in top ", length(lads), " LTLAs"))

    if(cont) {
        p1[[1]] <- p1[[1]] + geom_vline(xintercept = tstart, linetype = "dashed")
        p1[[2]] <- p1[[2]] + geom_vline(xintercept = tstart, linetype = "dashed")
        p1[[3]] <- p1[[3]] + geom_vline(xintercept = tstart, linetype = "dashed")
    }
    p1 <- (p1[[1]] + p1[[2]]) / (p1[[3]] + plot_spacer())

    ggsave(paste0("../wave", wave, "/simsTopLADsBPFEns.pdf"), p1, width = 25, height = 25)
}

