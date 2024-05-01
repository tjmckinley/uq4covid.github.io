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
    select(!name)
DH <- readRDS(paste0("../", outputs, "/cumDH_age_lad.rds")) %>%
    filter(t <= tstop) %>%
    pivot_longer(!t, values_to = "n") %>%
    mutate(age = gsub('^(?:[^_]*_)(.*)', '\\1', name)) %>%
    mutate(lad = gsub('^(?:[^_]*_)(.*)', '\\1', age)) %>%
    mutate(age = gsub('(.*)_[0-9]*', '\\1', age)) %>%
    mutate(across(c(age, lad), as.numeric)) %>%
    select(!name)
H <- readRDS(paste0("../", outputs, "/H_age_lad.rds")) %>%
    filter(t <= tstop) %>%
    pivot_longer(!t, values_to = "n") %>%
    mutate(age = gsub('^(?:[^_]*_)(.*)', '\\1', name)) %>%
    mutate(lad = gsub('^(?:[^_]*_)(.*)', '\\1', age)) %>%
    mutate(age = gsub('(.*)_[0-9]*', '\\1', age)) %>%
    mutate(across(c(age, lad), as.numeric)) %>%
    select(!name)
Hcum <- readRDS(paste0("../", outputs, "/cumH_age_lad.rds")) %>%
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

## plot community deaths
p1[[2]] <- filter(sims_md, var == "DI") %>%
    select(!var) %>%
    pivot_wider(names_from = name, values_from = value) %>%
    ggplot(aes(x = t)) +
        geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5) +
        geom_ribbon(aes(ymin = LQ, ymax = UQ), alpha = 0.5) +
        geom_line(aes(y = Median)) +
        geom_line(
            aes(y = n), 
            data = group_by(DI, t, age) %>%
                summarise(n = sum(n), .groups = "drop"),
            col = "blue", linetype = "dashed"
        ) +
        facet_wrap(~age, labeller = labeller(age = age_label)) +
        xlab("Days") + 
        ylab("Counts") +
        ggtitle("Observed cumulative community deaths (aggregated over LTLAs)")
    
if(cont) p1[[2]] <- p1[[2]] + geom_vline(xintercept = tstart, linetype = "dashed")

## plot hospital deaths
p1[[3]] <- filter(sims_md, var == "DH") %>%
    select(!var) %>%
    pivot_wider(names_from = name, values_from = value) %>%
    ggplot(aes(x = t)) +
        geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5) +
        geom_ribbon(aes(ymin = LQ, ymax = UQ), alpha = 0.5) +
        geom_line(aes(y = Median)) +
        geom_line(
            aes(y = n), 
            data = group_by(DH, t, age) %>%
                summarise(n = sum(n), .groups = "drop"),
            col = "blue", linetype = "dashed"
        ) +
        facet_wrap(~age, labeller = labeller(age = age_label)) +
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
        facet_wrap(~age, labeller = labeller(age = age_label)) +
        xlab("Days") + 
        ylab("Counts") +
        ggtitle("Observed hospital cases (aggregated over LTLAs)")
    
if(cont) p1[[4]] <- p1[[4]] + geom_vline(xintercept = tstart, linetype = "dashed")

## plot cumulative hospital cases
p1[[5]] <- filter(sims_md, var == "Hcum") %>%
    select(!var) %>%
    pivot_wider(names_from = name, values_from = value) %>%
    ggplot(aes(x = t)) +
        geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5) +
        geom_ribbon(aes(ymin = LQ, ymax = UQ), alpha = 0.5) +
        geom_line(aes(y = Median)) +
        geom_line(
            aes(y = n), 
            data = group_by(Hcum, t, age) %>%
                summarise(n = sum(n), .groups = "drop"),
            col = "blue", linetype = "dashed"
        ) +
        facet_wrap(~age, labeller = labeller(age = age_label)) +
        xlab("Days") + 
        ylab("Counts") +
        ggtitle("Observed cumulative hospital cases (aggregated over LTLAs)")
    
if(cont) p1[[5]] <- p1[[5]] + geom_vline(xintercept = tstart, linetype = "dashed")
            
###############################################
#####   combine plots and save outputs    #####
###############################################

## save outputs
saveRDS(p1, paste0("../wave", wave, "/plots.rds"))

## combine plots
layout <- "
AAB
CDE
"
p1 <- p1[[1]] + p1[[2]] + p1[[3]] + p1[[4]] + p1[[5]] + plot_layout(design = layout)
ggsave(paste0("../wave", wave, "/simsBPFEns.pdf"), p1, width = 20, height = 15)

###############################################
#####           spatial plots             #####
###############################################

## load in shapefile
lad19 <- st_read(paste0("../", outputs, "/Local_Authority_Districts_(December_2019)_Boundaries_UK_BUC.shp"))

## load in lookup
death_lookup <- readRDS(paste0("../", outputs, "/death_lookup.rds"))

## load in data
DH <- rename(DH, DH = n)
Hcum <- rename(Hcum, Hcum = n)
data <- inner_join(DH, Hcum, by = c("lad", "t", "age")) %>%
    inner_join(death_lookup, by = c("lad" = "FID"))
DH <- rename(DH, n = DH)
Hcum <- rename(Hcum, n = Hcum)

## load in runs
sims_md <- readRDS(paste0("../wave", wave, "/sumEns_age_lads.rds")) %>%
    dplyr::select(t, lad, age, DH_Median, Hcum_Median)

## join runs and data
data <- inner_join(data, sims_md, by = c("lad", "t", "age"))

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
temp <- filter(lad19, t == max(t))
p <- list()
p[[1]] <- rename(temp, Count = DH) %>%
    ggplot() +
        geom_sf(aes(fill = Count), colour = NA) +
        facet_wrap(~age, labeller = labeller(age = age_label)) +
        scale_fill_viridis_c() +
        ggtitle("Data (cumulative hospital deaths)")
p[[2]] <- rename(temp, Count = DH_Median) %>%
    ggplot() +
        geom_sf(aes(fill = Count), colour = NA) +
        facet_wrap(~age, labeller = labeller(age = age_label)) +
        scale_fill_viridis_c() +
        ggtitle("Predictions (cumulative hospital deaths)")
p1 <- p[[1]] + p[[2]] 
#p1 <- p1 & theme(legend.position = "bottom")
p1 <- p1 & scale_fill_viridis_c(limits = range(c(temp$DH, temp$DH_Median)))
p1 <- p1 + plot_layout(guides = "collect")
p1 <- p1 + plot_annotation(title = paste0("Cumulative hospital deaths at t = ", max(lad19$t)))
p2 <- list()
p2[[1]] <- p1
p[[1]] <- rename(temp, Count = Hcum) %>%
    ggplot() +
        geom_sf(aes(fill = Count), colour = NA) +
        facet_wrap(~age, labeller = labeller(age = age_label)) +
        scale_fill_viridis_c() +
        ggtitle("Data (cumulative hospital cases)")
p[[2]] <- rename(temp, Count = Hcum_Median) %>%
    ggplot() +
        geom_sf(aes(fill = Count), colour = NA) +
        facet_wrap(~age, labeller = labeller(age = age_label)) +
        scale_fill_viridis_c() +
        ggtitle("Predictions (cumulative hospital cases)")
p1 <- p[[1]] + p[[2]] 
#p1 <- p1 & theme(legend.position = "bottom")
p1 <- p1 & scale_fill_viridis_c(limits = range(c(temp$Hcum, temp$Hcum_Median)))
p1 <- p1 + plot_layout(guides = "collect")
p1 <- p1 + plot_annotation(title = paste0("Cumulative hospital cases at t = ", max(lad19$t)))
p2[[2]] <- p1
p2 <- p2[[1]] / p2[[2]]
ggsave(paste0("../wave", wave, "/simsspstatic.pdf"), p2, height = 14, width = 10)

###############################################
#####             LTLA plots              #####
###############################################

## load in runs
sims_md <- readRDS(paste0("../wave", wave, "/sumEns_age_lads.rds"))

## join runs and data
data <- inner_join(
    select(data, !c(areaCode, areaName, DH_Median, Hcum_Median)), 
    sims_md, 
    by = c("lad", "t", "age")
)

## produce plots
p <- filter(data, t == max(t)) %>%
    select(age, lad, Data = DH, Prediction = DH_Median, LCI = DH_LCI, UCI = DH_UCI) %>%
    group_by(age) %>%
    arrange(desc(Prediction)) %>%
    mutate(lad = 1:n()) %>%
    ungroup() %>%
    mutate(inside = ifelse(Data >= LCI & Data <= UCI, "Inside 95% CI", "Outside 95% CI")) %>%
    pivot_longer(c(Data, Prediction), names_to = "type", values_to = "Count") %>%
    mutate(inside = ifelse(type == "Prediction", "Prediction", inside)) %>%
    group_by(age) %>%
    arrange(desc(inside)) %>%
    ungroup() %>%
    ggplot() +
        geom_errorbar(aes(x = lad, ymin = LCI, ymax = UCI), colour = "#52854C") +
        geom_point(aes(x = lad, y = Count, colour = inside)) +
        facet_wrap(~ age, ncol = 1, labeller = labeller(age = age_label), scales = "free_x") +
        ggtitle(paste0("Cumulative deaths at t = ", max(data$t))) +
        xlab("LTLA (in decreasing order of deaths)") +
        ylab("Cumulative deaths") +
        scale_colour_manual(name = "", values = c("#E69F00", "#D55E00", "#52854C")) +
        theme(axis.text.x = element_blank(), legend.position = "bottom")
ggsave(paste0("../wave", wave, "/simsLTLA_cDeaths.pdf"), p, width = 15, height = 15)

p <- filter(data, t == max(t)) %>%
    select(age, lad, Data = Hcum, Prediction = Hcum_Median, LCI = Hcum_LCI, UCI = Hcum_UCI) %>%
    group_by(age) %>%
    arrange(desc(Prediction)) %>%
    mutate(lad = 1:n()) %>%
    ungroup() %>%
    mutate(inside = ifelse(Data >= LCI & Data <= UCI, "Inside 95% CI", "Outside 95% CI")) %>%
    pivot_longer(c(Data, Prediction), names_to = "type", values_to = "Count") %>%
    mutate(inside = ifelse(type == "Prediction", "Prediction", inside)) %>%
    group_by(age) %>%
    arrange(desc(inside)) %>%
    ungroup() %>%
    ggplot() +
        geom_errorbar(aes(x = lad, ymin = LCI, ymax = UCI), colour = "#52854C") +
        geom_point(aes(x = lad, y = Count, colour = inside)) +
        facet_wrap(~ age, ncol = 1, labeller = labeller(age = age_label), scales = "free_x") +
        ggtitle(paste0("Cumulative hospital cases at t = ", max(data$t))) +
        xlab("LTLA (in decreasing order of deaths)") +
        ylab("Cumulative hospital cases") +
        scale_colour_manual(name = "", values = c("#E69F00", "#D55E00", "#52854C")) +
        theme(axis.text.x = element_blank(), legend.position = "bottom")
ggsave(paste0("../wave", wave, "/simsLTLA_cHosp.pdf"), p, width = 15, height = 15)

## save outputs
p2 <- list(p2, p)
saveRDS(p2, paste0("../wave", wave, "/plots_sp.rds"))

###############################################
#####          LAD-level plots            #####
###############################################

if(file.exists(paste0("../", outputs, "/lads_", outputs, ".txt"))) {

    ## read in lads
    lads <- as.numeric(readLines(paste0("../", outputs, "/lads_", outputs, ".txt")))
    
    ## load in data
    DI <- filter(DI, lad %in% lads)
    DH <- filter(DH, lad %in% lads)
    H <- filter(H, lad %in% lads)
    Hcum <- filter(Hcum, lad %in% lads)

    ## load in runs
    sims_md <- readRDS(paste0("../wave", wave, "/sumEns_age_lads.rds")) %>%
        filter(lad %in% lads) %>%
        pivot_longer(!c(t, age, lad)) %>%
        mutate(var = gsub('^([^_]*)_(.*)', '\\1', name)) %>%
        mutate(name = gsub('^(?:[^_]*_)(.*)', '\\1', name)) %>%
        pivot_wider(names_from = name, values_from = value)
        
    p1 <- list()
    p1[[1]] <- filter(sims_md, var == "DI") %>%
        ggplot(aes(x = t)) +
            geom_ribbon(aes(ymin = LCI, ymax = UCI), colour = NA, alpha = 0.5) +
            geom_ribbon(aes(ymin = LQ, ymax = UQ), colour = NA, alpha = 0.5) +
            geom_line(aes(y = Median)) +
            geom_line(
                aes(y = n), 
                data = DI,
                linetype = "dashed",
                col = "blue"
            ) +
            facet_grid(lad ~ age, scales = "free", labeller = labeller(age = age_label)) +
            xlab("Days") + 
            ylab("Counts") +
            ggtitle(paste0("Observed community deaths in top ", length(lads), " LTLAs"))
    p1[[2]] <- filter(sims_md, var == "DH") %>%
        ggplot(aes(x = t)) +
            geom_ribbon(aes(ymin = LCI, ymax = UCI), colour = NA, alpha = 0.5) +
            geom_ribbon(aes(ymin = LQ, ymax = UQ), colour = NA, alpha = 0.5) +
            geom_line(aes(y = Median)) +
            geom_line(
                aes(y = n), 
                data = DH,
                linetype = "dashed",
                col = "blue"
            ) +
            facet_grid(lad ~ age, scales = "free", labeller = labeller(age = age_label)) +
            xlab("Days") + 
            ylab("Counts") +
            ggtitle(paste0("Observed hospital deaths in top ", length(lads), " LTLAs"))
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
            facet_grid(lad ~ age, scales = "free", labeller = labeller(age = age_label)) +
            xlab("Days") + 
            ylab("Counts") +
            ggtitle(paste0("Observed hospital cases in top ", length(lads), " LTLAs"))
    p1[[4]] <- filter(sims_md, var == "Hcum") %>%
        ggplot(aes(x = t)) +
            geom_ribbon(aes(ymin = LCI, ymax = UCI), colour = NA, alpha = 0.5) +
            geom_ribbon(aes(ymin = LQ, ymax = UQ), colour = NA, alpha = 0.5) +
            geom_line(aes(y = Median)) +
            geom_line(
                aes(y = n), 
                data = Hcum,
                linetype = "dashed",
                col = "blue"
            ) +
            facet_grid(lad ~ age, scales = "free", labeller = labeller(age = age_label)) +
            xlab("Days") + 
            ylab("Counts") +
            ggtitle(paste0("Observed cumulative hospital cases in top ", length(lads), " LTLAs"))

    if(cont) {
        p1[[1]] <- p1[[1]] + geom_vline(xintercept = tstart, linetype = "dashed")
        p1[[2]] <- p1[[2]] + geom_vline(xintercept = tstart, linetype = "dashed")
        p1[[3]] <- p1[[3]] + geom_vline(xintercept = tstart, linetype = "dashed")
        p1[[4]] <- p1[[4]] + geom_vline(xintercept = tstart, linetype = "dashed")
    }
    p1 <- (p1[[1]] + p1[[2]]) / (p1[[3]] + p1[[4]])
    
    ## save plot
    ggsave(paste0("../wave", wave, "/simsTopLADsBPFEns.pdf"), p1, width = 25, height = 25)
    
    ## save outputs
    saveRDS(p1, paste0("../wave", wave, "/plots_lad.rds"))
}

