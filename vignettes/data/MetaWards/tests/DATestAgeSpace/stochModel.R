## load libraries
library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(parallel)
library(abind)
library(sitmo)
library(sf)
library(gganimate)
library(viridis)
library(patchwork)

## set seed
set.seed(666)

## create output directory
dir.create("outputs")

## source simulation function
sourceCpp("BPF.cpp")
source("BPF.R")

## read in parameters, remove guff and reorder
pars <- readRDS("wave1/disease.rds") %>%
    rename(nu = `beta[1]`, nuA = `beta[6]`) %>%
    select(!c(starts_with("beta["), repeats)) %>%
    select(nu, nuA, !c(beta_scale, output), beta_scale, output)

## read in contact matrices
contact1 <- read_csv("inputs/POLYMOD_matrix.csv", col_names = FALSE) %>%
    as.matrix()
contact2 <- read_csv("inputs/coMix_matrix.csv", col_names = FALSE) %>%
    as.matrix()

## extract parameters for simulation   
pars <- select(slice(pars, 100), !output)

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

## add seeding information
EW19 <- mutate(EW19, X4 = 0)
for(i in 1:10) {
    val <- 0
    while(val == 0) {
        seed <- sample(which(EW19[, 1] == 317), 1)
        if((EW19$X4[seed] + 1) <= EW19$X3[seed]) {
            EW19$X4[seed] <- EW19$X4[seed] + 1
            val <- 1
        }
    }
}
## expand to deal with age-classes
u1 <- apply(EW19, 1, function(x, ageProbs) {
        u <- matrix(0, 12, length(ageProbs))
        u[2, ] <- smart_round(ageProbs * x[4])
        u[1, ] <- smart_round(ageProbs * x[3]) - u[2, ]
        list(u)
    }, ageProbs = ageProbs) %>%
    map(1) %>%
    abind(along = 3)
u1_moves <- as.matrix(EW19[, 1:2])

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

## write inputs out
saveRDS(u1, "outputs/u1.rds")
saveRDS(u1_moves, "outputs/u1_moves.rds")
saveRDS(u2, "outputs/u2.rds")
saveRDS(as.matrix(PM19), "outputs/u2_moves.rds")

## set up stage names
stageNms <- map(c("S", "E", "A", "RA", "P", "Ione", "DI", "Itwo", "RI", "H", "RH", "DH"), ~paste0(., "_", 1:8)) %>%
    map(~map(., ~paste0(., "_", 1:max(EW19[, 1])))) %>%
    reduce(c) %>%
    reduce(c)
    
## load lookups
lookup <- readRDS("data/lookup.rds")
lookup <- select(lookup, starts_with("FID"))
age_lookup <- readRDS("data/age_lookup.rds")
saveRDS(lookup, "outputs/lookup.rds")
saveRDS(age_lookup, "outputs/age_lookup.rds")

## simulate discrete-time model
disSims_full <- BPF(pars, C1 = contact1, C2 = contact2, lockdown_day = 20, 
    lookup = lookup, age_lookup = age_lookup, u1_moves = u1_moves,
    u1 = u1, u2_moves = as.matrix(PM19), u2 = u2, ndays = 100, npart = 8, PF = FALSE)
    
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
        col = "red", linetype = "dashed"
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
        col = "red", linetype = "dashed"
    ) +
    facet_grid(RGN19NM ~ age, labeller = label_wrap_gen(width = 10)) +
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
        col = "red", linetype = "dashed"
    ) +
    facet_wrap(~ areaName, nrow = 1, labeller = label_wrap_gen(width = 10)) +
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
        col = "red", linetype = "dashed"
    ) +
    facet_grid(areaName ~ age, labeller = label_wrap_gen(width = 10)) +
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
ggsave("outputs/simsNational.pdf", p1, width = 15, height = 15)

## save parameters
saveRDS(pars, "outputs/pars.rds")

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
saveRDS(medRep, "outputs/disSims.rds")

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

## extract simulation closest to median
medRep <- inner_join(medRepInd, disSims, by = "rep") %>%
    select(!c(diff, rep)) %>%
    select(t, everything())
    
## save outputs
saveRDS(medRep, "outputs/cumDeath_lad.rds")

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

## extract simulation closest to median
medRep <- inner_join(medRepInd, disSims, by = "rep") %>%
    select(!c(diff, rep)) %>%
    select(t, everything())
    
## save outputs
saveRDS(medRep, "outputs/cumDeath_age_region.rds")

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
saveRDS(medRep, "outputs/hosp_nhsregion.rds")

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

## extract simulation closest to median
medRep <- inner_join(medRepInd, disSims, by = "rep") %>%
    select(!c(diff, rep)) %>%
    select(t, everything())
    
## save outputs
saveRDS(medRep, "outputs/cumHospAd_age_nhsregion.rds")

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
ggsave("outputs/simsTopLADs.pdf", p1, width = 10, height = 10)

## spatial animation of simulation

## read in shapefile
lad19 <- st_read("inputs/LAD19_shapefile/LAD19_shapefile.shp")

## extract cases over time in each age-group
p <- select(medRep, t, starts_with("DH_")) %>%
    pivot_longer(!t, values_to = "counts", names_to = "var") %>%
    separate(var, c("var", "age", "lad"), sep = "_") %>%
    select(!var) %>%
    mutate(age = as.numeric(age), lad = as.numeric(lad)) %>%
    group_by(lad) %>%
    mutate(tot = sum(counts)) %>%
    ungroup() %>%
    filter(tot > 0) %>%
    select(!tot) %>%
    group_by(lad, age) %>%
    mutate(tot = cumsum(counts)) %>%
    group_by(lad, t) %>%
    mutate(tot = sum(tot)) %>%
    ungroup() %>%
    mutate(counts = ifelse(tot == 0, NA, counts)) %>%
    select(!tot)
max_p <- max(p$counts, na.rm = TRUE)
p <- inner_join(lad19, p, by = c("objectid" = "lad")) %>%
    ggplot() +
    geom_sf(aes(fill = counts), colour = NA) +
    facet_wrap(~age) +
    scale_fill_viridis_c(limits = c(0, max_p)) +
    theme_bw()

## add transitions
p <- p + transition_time(t) + ggtitle("Day = {frame_time}")
spatial_gif <- animate(p, nframes = 50, fps = 1, renderer = gifski_renderer())
anim_save("outputs/simsSpatialDH.gif", spatial_gif)

## extract cases over time in each age-group
p <- select(medRep, t, starts_with("E")) %>%
    pivot_longer(!t, values_to = "counts", names_to = "var") %>%
    separate(var, c("var", "age", "lad"), sep = "_") %>%
    select(!var) %>%
    mutate(age = as.numeric(age), lad = as.numeric(lad)) %>%
    group_by(lad) %>%
    mutate(tot = sum(counts)) %>%
    ungroup() %>%
    filter(tot > 0) %>%
    select(!tot) %>%
    group_by(lad, age) %>%
    mutate(tot = cumsum(counts)) %>%
    group_by(lad, t) %>%
    mutate(tot = sum(tot)) %>%
    ungroup() %>%
    mutate(counts = ifelse(tot == 0, NA, counts)) %>%
    select(!tot)
max_p <- max(p$counts, na.rm = TRUE)
p <- inner_join(lad19, p, by = c("objectid" = "lad")) %>%
    ggplot() +
    geom_sf(aes(fill = counts), colour = NA) +
    facet_wrap(~age) +
    scale_fill_viridis_c(limits = c(0, max_p)) +
    theme_bw()

## add transitions
p <- p + transition_time(t) + ggtitle("Day = {frame_time}")
spatial_gif <- animate(p, nframes = 50, fps = 1, renderer = gifski_renderer())
anim_save("outputs/simsSpatialE.gif", spatial_gif)

