## load libraries
library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(parallel)
library(abind)
library(sitmo)
library(sf)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) > 0) {
    if(length(args) != 3) {
        stop("Must be three arguments")
    }
    seed <- as.numeric(args[1])
    inputdir <- paste0("wave1", args[2])
    id <- as.numeric(args[3])
} else {
    seed <- 456
    inputdir <- "wave1Feb"
    id <- 1
}

## set seed
set.seed(seed)

## create output directory
outputdir <- paste0("outputs", seed)
if(dir.exists(outputdir)) {
    stop(outputdir, " directory already exists. Please delete if you wish to overwrite.")
}
dir.create(outputdir)

## source simulation function
sourceCpp("BPF.cpp")
source("BPF.R")

## read in fixed input data
fixedInputs <- readLines(paste0(inputdir, "/fixedInputs.txt"))

## set case specific values
tstart <- as.numeric(fixedInputs[1])
tstop <- as.numeric(fixedInputs[2])
lockdown_day <- as.numeric(fixedInputs[3])
npart <- as.numeric(fixedInputs[4])
niter <- as.numeric(fixedInputs[5])
a1 <- as.numeric(fixedInputs[6])
a2 <- as.numeric(fixedInputs[7])
b1 <- as.numeric(fixedInputs[8])
b2 <- as.numeric(fixedInputs[9])
a_dis_ini <- as.numeric(fixedInputs[10])
b_dis_ini <- as.numeric(fixedInputs[11])
b_dis_8_ini <- as.numeric(fixedInputs[12])
sigma2_lad <- as.numeric(fixedInputs[13])
sigma2_age_region <- as.numeric(fixedInputs[14])
sigma2_nhsregion <- as.numeric(fixedInputs[15])
sigma2_age_nhsregion <- as.numeric(fixedInputs[16])
saveAll <- as.logical(as.numeric(fixedInputs[17]))
snapshot <- as.logical(as.numeric(fixedInputs[18]))
writeExt <- as.logical(as.numeric(fixedInputs[19]))

## read in parameters, remove guff and reorder
pars <- readRDS(paste0(inputdir, "/disease.rds"))

## read in contact matrices
contact1 <- read_csv("inputs/POLYMOD_matrix.csv", col_names = FALSE) %>%
    as.matrix()
contact2 <- read_csv("inputs/coMix_matrix.csv", col_names = FALSE) %>%
    as.matrix()

## extract parameters for simulation   
pars <- select(slice(pars, id), !output)
pars_save <- readRDS(paste0(inputdir, "/inputs.rds")) %>%
    slice(id)
saveRDS(pars_save, paste0(outputdir, "/pars.rds"))

## solution to round numbers preserving sum
## adapted from:
## https://stackoverflow.com/questions/32544646/round-vector-of-numerics-to-integer-while-preserving-their-sum
smart_round <- function(x) {
    y <- floor(x)
    indices <- tail(order(x - y), round(sum(x)) - sum(y))
    y[indices] <- y[indices] + 1
    y
}

## add age probabilities
ageProbs <- read_csv("inputs/age_seeds.csv", col_names = FALSE)$X2

## read in commuter data
EW19 <- read_delim("inputs/EW19.dat", delim = " ", col_names = FALSE)
u1_moves <- as.matrix(EW19[, 1:2])

## expand to deal with age-classes
u1 <- apply(EW19, 1, function(x, ageProbs) {
        u <- matrix(0, 12, length(ageProbs))
        u[1, ] <- smart_round(ageProbs * x[3])
        list(u)
    }, ageProbs = ageProbs) %>%
    map(1) %>%
    abind(along = 3)

## read in player data
PM19 <- read_delim("inputs/PlayMatrix19.dat", delim = " ", col_names = FALSE)
PlaySize19 <- read_delim("inputs/PlaySize19.dat", delim = " ", col_names = FALSE)
  
## expand to deal with age-classes
u2 <- apply(PlaySize19, 1, function(x, ageProbs) {
        u <- matrix(0, 12, length(ageProbs))
        u[1, ] <- smart_round(ageProbs * x[2])
        list(u)
    }, ageProbs = ageProbs) %>%
    map(1) %>%
    abind(along = 3)

## set up stage names
stageNms <- map(c("S", "E", "A", "RA", "P", "Ione", "DI", "Itwo", "RI", "H", "RH", "DH"), ~paste0(., "_", 1:8)) %>%
    map(~map(., ~paste0(., "_", 1:max(EW19[, 1])))) %>%
    reduce(c) %>%
    reduce(c)
    
## load lookups
lookup <- readRDS("data/lookup.rds")
lookup <- select(lookup, starts_with("FID"))
age_lookup <- readRDS("data/age_lookup.rds")
death_lookup <- readRDS("data/death_lookup.rds")
saveRDS(lookup, paste0(outputdir, "/lookup.rds"))
saveRDS(age_lookup, paste0(outputdir, "/age_lookup.rds"))
saveRDS(death_lookup, paste0(outputdir, "/death_lookup.rds"))

## set up input object
u <- list(u1_moves = u1_moves, u1 = u1, u2 = u2, u2_moves = as.matrix(PM19))

## create region lookup
region_lookup <- readRDS("data/region_lookup.rds") %>%
    mutate(initials = NA) %>%
    mutate(initials = ifelse(RGN19NM == "North East", "NE", initials)) %>%
    mutate(initials = ifelse(RGN19NM == "North West", "NW", initials)) %>%
    mutate(initials = ifelse(RGN19NM == "Yorkshire and The Humber", "YH", initials)) %>%
    mutate(initials = ifelse(RGN19NM == "East Midlands", "EM", initials)) %>%
    mutate(initials = ifelse(RGN19NM == "West Midlands", "WM", initials)) %>%
    mutate(initials = ifelse(RGN19NM == "East of England", "EE", initials)) %>%
    mutate(initials = ifelse(RGN19NM == "London", "L", initials)) %>%
    mutate(initials = ifelse(RGN19NM == "South East", "SE", initials)) %>%
    mutate(initials = ifelse(RGN19NM == "South West", "SW", initials))

## set some initial matrices
a_dis <- a_dis_ini * exp(-pars$MD_scale[1] * (pars$MD_time_L[1] - tstart:tstop) * ifelse(tstart:tstop < pars$MD_time_L[1], 1, 0))
a_dis <- array(rep(a_dis, nrow(age_lookup) * nrow(lookup)), c(tstop - tstart + 1, nrow(age_lookup), nrow(lookup)))
b_dis <- b_dis_ini * exp(-pars$MD_scale[1] * (pars$MD_time_L[1] - tstart:tstop) * ifelse(tstart:tstop < pars$MD_time_L[1], 1, 0))
b_dis <- array(rep(b_dis, nrow(age_lookup) * nrow(lookup)), c(tstop - tstart + 1, nrow(age_lookup), nrow(lookup)))

## use spatially-explicit matrices
for(i in 1:nrow(region_lookup)) {
    FID <- lookup$FID[!is.na(lookup$FID_region) & lookup$FID_region == region_lookup$FID[region_lookup$RGN19NM == region_lookup$RGN19NM[i]]]
    MD_time <- pluck(pars, paste0("MD_time_", region_lookup$initials[i]))
    a_dis_temp <- a_dis_ini * exp(-pars$MD_scale[1] * (MD_time[1] - tstart:tstop) * ifelse(tstart:tstop < MD_time[1], 1, 0))
    a_dis[, , FID] <- rep(a_dis_temp, length(FID))
    b_dis_temp <- b_dis_ini * exp(-pars$MD_scale[1] * (MD_time[1] - tstart:tstop) * ifelse(tstart:tstop < MD_time[1], 1, 0))
    b_dis[, , FID] <- rep(b_dis_temp, length(FID))
    b_dis_8 <- b_dis_8_ini * exp(-pars$MD_scale[1] * (MD_time[1] - tstart:tstop) * ifelse(tstart:tstop < MD_time[1], 1, 0))
    b_dis[, 8, FID] <- rep(b_dis_8, length(FID))
}

## simulate discrete-time model
disSims_full <- BPF(pars, C1 = contact1, C2 = contact2, lockdown_day = lockdown_day, 
    lookup = lookup, age_lookup = age_lookup, u = u, a1 = a1, a2 = a2,
    b1 = b1, b2 = b2, a_dis = a_dis, b_dis = b_dis, tstart = tstart, tstop = tstop,
    sigma2_lad = sigma2_lad, sigma2_age_region = sigma2_age_region, 
    sigma2_nhsregion = sigma2_nhsregion, sigma2_age_nhsregion = sigma2_age_nhsregion,
    npart = 8, PF = FALSE)
    
###############################################
#######        LAD-level truth          #######
###############################################
        
## collapse to data frame
disSims <- map(1:length(disSims_full$particles[[1]]$full), function(i, x) {
        ## loop over particles
        map(1:length(x[[i]]), function(i, x) {
            ## collapse to vector in order: stage, age, LAD
            aperm(x[[i]], 3:1) %>%
            as.vector() %>%
            c(i)
        }, x = x[[i]]) %>%
        {do.call("rbind", .)} %>%
        cbind(rep(i, nrow(.)))
    }, x = disSims_full$particles[[1]]$full) %>%
    {do.call("rbind", .)}
colnames(disSims) <- c(stageNms, "rep", "t")
disSims <- as_tibble(disSims) %>%
    mutate(t = t - 1)

## extract simulation closest to median
medRepInd <- pivot_longer(disSims, !c(rep, t), names_to = "var", values_to = "n") %>%
    group_by(t, var) %>%
    summarise(
        median = median(n),
        .groups = "drop"
    ) %>%
    inner_join(
        pivot_longer(disSims, !c(rep, t), names_to = "var", values_to = "n"),
        by = c("t", "var")
    ) %>%
    mutate(diff = (median - n)^2) %>%
    group_by(rep) %>%
    summarise(diff = sum(diff), .groups = "drop") %>%
    arrange(diff) %>%
    slice(1) 
medRep <- inner_join(medRepInd, disSims, by = "rep") %>%
    select(!c(diff, rep))

## plot replicates at national level
p <- pivot_longer(disSims, !c(rep, t), names_to = "var", values_to = "n") %>%
    mutate(var = gsub("_[0-9]*$", "", var)) %>%
    group_by(rep, t, var) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    group_by(t, var) %>%
    summarise(
        LCI = quantile(n, probs = 0.025),
        LQ = quantile(n, probs = 0.25),
        median = median(n),
        UQ = quantile(n, probs = 0.75),
        UCI = quantile(n, probs = 0.975),
        .groups = "drop"
    ) %>%
    mutate(age = gsub("[^0-9]", "", var)) %>%
    mutate(var = gsub("[^a-zA-Z]", "", var)) %>%
    mutate(var = gsub("one", "1", var)) %>%
    mutate(var = gsub("two", "2", var))
    
p1 <- list()     
p1[[1]] <- ggplot(p, aes(x = t)) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5) +
    geom_ribbon(aes(ymin = LQ, ymax = UQ), alpha = 0.5) +
    geom_line(aes(y = median)) +
    geom_line(
        aes(y = n), 
        data = pivot_longer(medRep, !t, names_to = "var", values_to = "n") %>%
            mutate(var = gsub("_[0-9]*$", "", var)) %>%
            group_by(t, var) %>%
            summarise(n = sum(n), .groups = "drop") %>%
            mutate(age = gsub("[^0-9]", "", var)) %>%
            mutate(var = gsub("[^a-zA-Z]", "", var)) %>%
            mutate(var = gsub("one", "1", var)) %>%
            mutate(var = gsub("two", "2", var)),
        col = "red", linetype = "dashed"
    ) +
    facet_grid(var ~ age, scales = "free") +
    xlab("Days") + 
    ylab("Counts") +
    ggtitle("Truth")
        
###############################################
#######     LAD-level observations      #######
###############################################
        
## collapse to data frame
disSims <- map(1:length(disSims_full$particles[[1]]$lads), function(i, x) {
        ## loop over particles
        map(1:length(x[[i]]), function(i, x) {
            as.vector(x[[i]]) %>%
            c(i)
        }, x = x[[i]]) %>%
        {do.call("rbind", .)} %>%
        cbind(rep(i, nrow(.)))
    }, x = disSims_full$particles[[1]]$lads) %>%
    {do.call("rbind", .)}
colnames(disSims) <- c(paste0("deaths_", 1:(ncol(disSims) - 2)), "rep", "t")
disSims <- as_tibble(disSims) %>%
    mutate(t = t - 1)
    
## convert to cumulative count for plotting
disSims <- group_by(disSims, rep) %>%
    mutate(across(starts_with("deaths_"), cumsum)) %>%
    ungroup()

## extract simulation closest to median
medRep <- inner_join(medRepInd, disSims, by = "rep") %>%
    select(!c(diff, rep))

## plot replicates at national level
p <- pivot_longer(disSims, !c(rep, t), names_to = "var", values_to = "n") %>%
    group_by(rep, t) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    group_by(t) %>%
    summarise(
        LCI = quantile(n, probs = 0.025),
        LQ = quantile(n, probs = 0.25),
        median = median(n),
        UQ = quantile(n, probs = 0.75),
        UCI = quantile(n, probs = 0.975),
        .groups = "drop"
    )
     
p1[[2]] <- ggplot(p, aes(x = t)) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5) +
    geom_ribbon(aes(ymin = LQ, ymax = UQ), alpha = 0.5) +
    geom_line(aes(y = median)) +
    geom_line(
        aes(y = n), 
        data = pivot_longer(medRep, !t, names_to = "var", values_to = "n") %>%
            group_by(t) %>%
            summarise(n = sum(n), .groups = "drop"),
        col = "blue", linetype = "dashed"
    ) +
    xlab("Days") + 
    ylab("Counts") +
    ggtitle("Observed cumulative deaths (aggregated over LADs)")
       
###############################################
#######  age/region-level observations  #######
###############################################

## load lookup for labels
region_lookup <- readRDS("data/region_lookup.rds")
        
## collapse to data frame
disSims <- map(1:length(disSims_full$particles[[1]]$age_region), function(i, x) {
        ## loop over particles
        map(1:length(x[[i]]), function(i, x) {
            as.vector(x[[i]]) %>%
            c(i)
        }, x = x[[i]]) %>%
        {do.call("rbind", .)} %>%
        cbind(rep(i, nrow(.)))
    }, x = disSims_full$particles[[1]]$age_region) %>%
    {do.call("rbind", .)}
temp <- disSims_full$particles[[1]]$age_region[[1]][[1]]
tempNames <- rep(1:nrow(temp), times = ncol(temp))
tempNames <- cbind(tempNames, rep(1:ncol(temp), each = nrow(temp)))
tempNames <- paste0("deaths_", apply(tempNames, 1, paste, collapse = "_"))
colnames(disSims) <- c(tempNames, "rep", "t")
disSims <- as_tibble(disSims) %>%
    mutate(t = t - 1)
    
## convert to cumulative count for plotting
disSims <- group_by(disSims, rep) %>%
    mutate(across(starts_with("deaths_"), cumsum)) %>%
    ungroup()

## extract simulation closest to median
medRep <- inner_join(medRepInd, disSims, by = "rep") %>%
    select(!c(diff, rep))

## plot replicates at national level
p <- pivot_longer(disSims, !c(rep, t), names_to = "var", values_to = "n") %>%
    mutate(var = gsub("deaths_", "", var)) %>%
    separate(var, c("age", "region"), sep = "_") %>%
    mutate(across(c(age, region), as.numeric)) %>%
    inner_join(region_lookup, by = c("region" = "FID")) %>%
    group_by(t, age, RGN19NM) %>%
    summarise(
        LCI = quantile(n, probs = 0.025),
        LQ = quantile(n, probs = 0.25),
        median = median(n),
        UQ = quantile(n, probs = 0.75),
        UCI = quantile(n, probs = 0.975),
        .groups = "drop"
    )
     
p1[[3]] <- ggplot(p, aes(x = t)) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5) +
    geom_ribbon(aes(ymin = LQ, ymax = UQ), alpha = 0.5) +
    geom_line(aes(y = median)) +
    geom_line(
        aes(y = n), 
        data = pivot_longer(medRep, !t, names_to = "var", values_to = "n") %>%
            mutate(var = gsub("deaths_", "", var)) %>%
            separate(var, c("age", "region"), sep = "_") %>%
            mutate(across(c(age, region), as.numeric)) %>%
            inner_join(region_lookup, by = c("region" = "FID")),
        col = "blue", linetype = "dashed"
    ) +
    facet_grid(RGN19NM ~ age, labeller = label_wrap_gen(width = 10), scales = "free_y") +
    xlab("Days") + 
    ylab("Counts") +
    ggtitle("Observed cumulative deaths (age / region)")
        
###############################################
#######     NHS region observations     #######
###############################################

## load lookup for labels
nhsregion_lookup <- readRDS("data/nhsregion_lookup.rds")
        
## collapse to data frame
disSims <- map(1:length(disSims_full$particles[[1]]$nhsregion), function(i, x) {
        ## loop over particles
        map(1:length(x[[i]]), function(i, x) {
            as.vector(x[[i]]) %>%
            c(i)
        }, x = x[[i]]) %>%
        {do.call("rbind", .)} %>%
        cbind(rep(i, nrow(.)))
    }, x = disSims_full$particles[[1]]$nhsregion) %>%
    {do.call("rbind", .)}
colnames(disSims) <- c(paste0("hosp_", 1:(ncol(disSims) - 2)), "rep", "t")
disSims <- as_tibble(disSims) %>%
    mutate(t = t - 1)

## extract simulation closest to median
medRep <- inner_join(medRepInd, disSims, by = "rep") %>%
    select(!c(diff, rep))
    
## plot replicates at national level
p <- pivot_longer(disSims, !c(rep, t), names_to = "region", values_to = "n") %>%
    mutate(region = as.numeric(gsub("hosp_", "", region))) %>%
    inner_join(nhsregion_lookup, by = c("region" = "FID")) %>%
    group_by(t, areaName) %>%
    summarise(
        LCI = quantile(n, probs = 0.025),
        LQ = quantile(n, probs = 0.25),
        median = median(n),
        UQ = quantile(n, probs = 0.75),
        UCI = quantile(n, probs = 0.975),
        .groups = "drop"
    )
     
p1[[4]] <- ggplot(p, aes(x = t)) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5) +
    geom_ribbon(aes(ymin = LQ, ymax = UQ), alpha = 0.5) +
    geom_line(aes(y = median)) +
    geom_line(
        aes(y = n), 
        data = pivot_longer(medRep, !t, names_to = "region", values_to = "n") %>%
            mutate(region = as.numeric(gsub("hosp_", "", region))) %>%
            inner_join(nhsregion_lookup, by = c("region" = "FID")),
        col = "blue", linetype = "dashed"
    ) +
    facet_wrap(~ areaName, nrow = 1, labeller = label_wrap_gen(width = 10), scales = "free_y") +
    xlab("Days") + 
    ylab("Counts") +
    ggtitle("Observed hospital cases (NHS region)")
        
###############################################
#####  NHS age/region-level observations  #####
###############################################
        
## collapse to data frame
disSims <- map(1:length(disSims_full$particles[[1]]$age_nhsregion), function(i, x) {
        ## loop over particles
        map(1:length(x[[i]]), function(i, x) {
            as.vector(x[[i]]) %>%
            c(i)
        }, x = x[[i]]) %>%
        {do.call("rbind", .)} %>%
        cbind(rep(i, nrow(.)))
    }, x = disSims_full$particles[[1]]$age_nhsregion) %>%
    {do.call("rbind", .)}
temp <- disSims_full$particles[[1]]$age_nhsregion[[1]][[1]]
tempNames <- rep(1:nrow(temp), times = ncol(temp))
tempNames <- cbind(tempNames, rep(1:ncol(temp), each = nrow(temp)))
tempNames <- paste0("hospInc_", apply(tempNames, 1, paste, collapse = "_"))
colnames(disSims) <- c(tempNames, "rep", "t")
disSims <- as_tibble(disSims) %>%
    mutate(t = t - 1)
    
## convert to cumulative count for plotting
disSims <- group_by(disSims, rep) %>%
    mutate(across(starts_with("hospInc_"), cumsum)) %>%
    ungroup()

## extract simulation closest to median
medRep <- inner_join(medRepInd, disSims, by = "rep") %>%
    select(!c(diff, rep))

## plot replicates at national level
p <- pivot_longer(disSims, !c(rep, t), names_to = "var", values_to = "n") %>%
    mutate(var = gsub("hospInc_", "", var)) %>%
    separate(var, c("age", "region"), sep = "_") %>%
    mutate(across(c(age, region), as.numeric)) %>%
    inner_join(nhsregion_lookup, by = c("region" = "FID")) %>%
    group_by(t, age, areaName) %>%
    summarise(
        LCI = quantile(n, probs = 0.025),
        LQ = quantile(n, probs = 0.25),
        median = median(n),
        UQ = quantile(n, probs = 0.75),
        UCI = quantile(n, probs = 0.975),
        .groups = "drop"
    )
     
p1[[5]] <- ggplot(p, aes(x = t)) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5) +
    geom_ribbon(aes(ymin = LQ, ymax = UQ), alpha = 0.5) +
    geom_line(aes(y = median)) +
    geom_line(
        aes(y = n), 
        data = pivot_longer(medRep, !t, names_to = "var", values_to = "n") %>%
            mutate(var = gsub("hospInc_", "", var)) %>%
            separate(var, c("age", "region"), sep = "_") %>%
            mutate(across(c(age, region), as.numeric)) %>%
            inner_join(nhsregion_lookup, by = c("region" = "FID")),
        col = "blue", linetype = "dashed"
    ) +
    facet_grid(areaName ~ age, labeller = label_wrap_gen(width = 10), scales = "free_y") +
    xlab("Days") + 
    ylab("Counts") +
    ggtitle("Observed cumulative hospital incidence (NHS age / region)")
        
###############################################
#####   combine plots and save outputs    #####
###############################################

## combine plots
p1[[4]] <- p1[[4]] / p1[[2]]
p1 <- p1[-2]
p1 <- wrap_plots(p1, nrow = 2, heights = c(0.8, 0.5))
ggsave(paste0(outputdir, "/simsNational.pdf"), p1, width = 15, height = 15)

## truth
disSims <- map(1:length(disSims_full$particles[[1]]$full), function(i, x) {
        ## loop over particles
        map(1:length(x[[i]]), function(i, x) {
            ## collapse to vector in order: stage, age, LAD
            aperm(x[[i]], 3:1) %>%
            as.vector() %>%
            c(i)
        }, x = x[[i]]) %>%
        {do.call("rbind", .)} %>%
        cbind(rep(i, nrow(.)))
    }, x = disSims_full$particles[[1]]$full) %>%
    {do.call("rbind", .)}
colnames(disSims) <- c(stageNms, "rep", "t")
disSims <- as_tibble(disSims) %>%
    mutate(t = t - 1)

## extract simulation closest to median
medRep <- inner_join(medRepInd, disSims, by = "rep") %>%
    select(!c(diff, rep)) %>%
    select(t, everything())
    
## save outputs
saveRDS(medRep, paste0(outputdir, "/disSims.rds"))

## lad-level observed deaths
disSims <- map(1:length(disSims_full$particles[[1]]$lads), function(i, x) {
        ## loop over particles
        map(1:length(x[[i]]), function(i, x) {
            as.vector(x[[i]]) %>%
            c(i)
        }, x = x[[i]]) %>%
        {do.call("rbind", .)} %>%
        cbind(rep(i, nrow(.)))
    }, x = disSims_full$particles[[1]]$lads) %>%
    {do.call("rbind", .)}
colnames(disSims) <- c(paste0("deaths_", 1:(ncol(disSims) - 2)), "rep", "t")
disSims <- as_tibble(disSims) %>%
    mutate(t = t - 1)
    
## convert to cumulative count for plotting
disSims <- group_by(disSims, rep) %>%
    mutate(across(starts_with("deaths_"), cumsum)) %>%
    ungroup()

## extract simulation closest to median
medRep <- inner_join(medRepInd, disSims, by = "rep") %>%
    select(!c(diff, rep)) %>%
    select(t, everything())
    
## save outputs
saveRDS(medRep, paste0(outputdir, "/cumDeath_lad.rds"))

## age / region observed deaths
disSims <- map(1:length(disSims_full$particles[[1]]$age_region), function(i, x) {
        ## loop over particles
        map(1:length(x[[i]]), function(i, x) {
            as.vector(t(x[[i]])) %>%
            c(i)
        }, x = x[[i]]) %>%
        {do.call("rbind", .)} %>%
        cbind(rep(i, nrow(.)))
    }, x = disSims_full$particles[[1]]$age_region) %>%
    {do.call("rbind", .)}
temp <- t(disSims_full$particles[[1]]$age_region[[1]][[1]])
tempNames <- rep(1:ncol(temp), each = nrow(temp))
tempNames <- cbind(tempNames, rep(1:nrow(temp), times = ncol(temp)))
tempNames <- paste0("deaths_", apply(tempNames, 1, paste, collapse = "_"))
colnames(disSims) <- c(tempNames, "rep", "t")
disSims <- as_tibble(disSims) %>%
    mutate(t = t - 1)
    
## convert to cumulative count for plotting
disSims <- group_by(disSims, rep) %>%
    mutate(across(starts_with("deaths_"), cumsum)) %>%
    ungroup()

## extract simulation closest to median
medRep <- inner_join(medRepInd, disSims, by = "rep") %>%
    select(!c(diff, rep)) %>%
    select(t, everything())
    
## save outputs
saveRDS(medRep, paste0(outputdir, "/cumDeath_age_region.rds"))

## NHS region observed hospital cases
disSims <- map(1:length(disSims_full$particles[[1]]$nhsregion), function(i, x) {
        ## loop over particles
        map(1:length(x[[i]]), function(i, x) {
            as.vector(x[[i]]) %>%
            c(i)
        }, x = x[[i]]) %>%
        {do.call("rbind", .)} %>%
        cbind(rep(i, nrow(.)))
    }, x = disSims_full$particles[[1]]$nhsregion) %>%
    {do.call("rbind", .)}
colnames(disSims) <- c(paste0("hosp_", 1:(ncol(disSims) - 2)), "rep", "t")
disSims <- as_tibble(disSims) %>%
    mutate(t = t - 1)

## extract simulation closest to median
medRep <- inner_join(medRepInd, disSims, by = "rep") %>%
    select(!c(diff, rep)) %>%
    select(t, everything())
    
## save outputs
saveRDS(medRep, paste0(outputdir, "/hosp_nhsregion.rds"))

## NHS region and age observed hospital incidence
disSims <- map(1:length(disSims_full$particles[[1]]$age_nhsregion), function(i, x) {
        ## loop over particles
        map(1:length(x[[i]]), function(i, x) {
            as.vector(t(x[[i]])) %>%
            c(i)
        }, x = x[[i]]) %>%
        {do.call("rbind", .)} %>%
        cbind(rep(i, nrow(.)))
    }, x = disSims_full$particles[[1]]$age_nhsregion) %>%
    {do.call("rbind", .)}
temp <- t(disSims_full$particles[[1]]$age_nhsregion[[1]][[1]])
tempNames <- rep(1:ncol(temp), each = nrow(temp))
tempNames <- cbind(tempNames, rep(1:nrow(temp), times = ncol(temp)))
tempNames <- paste0("hospInc_", apply(tempNames, 1, paste, collapse = "_"))
colnames(disSims) <- c(tempNames, "rep", "t")
disSims <- as_tibble(disSims) %>%
    mutate(t = t - 1)
    
## convert to cumulative count for plotting
disSims <- group_by(disSims, rep) %>%
    mutate(across(starts_with("hospInc_"), cumsum)) %>%
    ungroup()

## extract simulation closest to median
medRep <- inner_join(medRepInd, disSims, by = "rep") %>%
    select(!c(diff, rep)) %>%
    select(t, everything())
    
## save outputs
saveRDS(medRep, paste0(outputdir, "/cumHospAd_age_nhsregion.rds"))

###############################################
#####          LAD-level plots            #####
###############################################

## collapse to data frame
disSims <- map(1:length(disSims_full$particles[[1]]$full), function(i, x) {
        ## loop over particles
        map(1:length(x[[i]]), function(i, x) {
            ## collapse to vector in order: stage, age, LAD
            aperm(x[[i]], 3:1) %>%
            as.vector() %>%
            c(i)
        }, x = x[[i]]) %>%
        {do.call("rbind", .)} %>%
        cbind(rep(i, nrow(.)))
    }, x = disSims_full$particles[[1]]$full) %>%
    {do.call("rbind", .)}
colnames(disSims) <- c(stageNms, "rep", "t")
disSims <- as_tibble(disSims) %>%
    mutate(t = t - 1)

## extract simulation closest to median
medRep <- inner_join(medRepInd, disSims, by = "rep") %>%
    select(!c(diff, rep))

## plot medRep at LAD level for LADs with largest epidemic load
top5 <- pivot_longer(medRep, !t, names_to = "var", values_to = "n") %>%
    mutate(age = gsub('^(?:[^_]*_)(.*)', '\\1', var)) %>%
    mutate(LAD = gsub('^(?:[^_]*_)(.*)', '\\1', age)) %>%
    mutate(age = gsub('(.*)_[0-9]*', '\\1', age)) %>%
    mutate(var = gsub('^(.*)_[0-9]*_.*', '\\1', var)) %>%
    filter(var == "DI") %>%
    filter(t == max(t)) %>%
    group_by(LAD) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    arrange(desc(n)) %>%
    slice(1:5) %>%
    select(-n)
p <- inner_join(top5,
        pivot_longer(medRep, !t, names_to = "var", values_to = "n") %>%
            mutate(age = gsub('^(?:[^_]*_)(.*)', '\\1', var)) %>%
            mutate(LAD = gsub('^(?:[^_]*_)(.*)', '\\1', age)),
        by = "LAD"
    ) %>%
    mutate(age = gsub('(.*)_[0-9]*', '\\1', age)) %>%
    mutate(var = gsub('^(.*)_[0-9]*_.*', '\\1', var)) %>%
    mutate(var = gsub("one", "1", var)) %>%
    mutate(var = gsub("two", "2", var))
p1 <- ggplot(p, aes(x = t, y = n, colour = LAD)) +
    geom_line() +
    facet_grid(var ~ age, scales = "free") +
    xlab("Days") + 
    ylab("Counts") +
    ggtitle("Truth")
ggsave(paste0(outputdir, "/simsTopLADs.pdf"), p1, width = 10, height = 10)

## copy shapefiles in
system(paste0("cp data/Local_Authority_Districts_\\(December_2019\\)_Boundaries_UK_BUC.zip ", outputdir))

## write out top LADs
readRDS(paste0(outputdir, "/cumDeath_lad.rds")) %>%
    filter(t == tstop) %>%
    pivot_longer(!t, names_to = "lad") %>%
    mutate(lad = as.numeric(gsub("deaths_", "", lad))) %>%
    arrange(desc(value)) %>%
    slice(1:10) %>%
    select(lad) %>%
    write_delim(file = paste0(outputdir, "/lads_outputs", seed, ".txt"), col_names = FALSE)

