## load libraries
library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(parallel)
library(abind)
library(sitmo)
library(patchwork)
library(data.table)

## source Rcpp PF code
sourceCpp("../BPF.cpp")

## source function to run PF and return log-likelihood
source("../BPF.R")

## read in simulated data
cumDeath_lad <- readRDS("../outputs/cumDeath_lad.rds")
cumDeath_age_region <- readRDS("../outputs/cumDeath_age_region.rds")
hosp_nhsregion <- readRDS("../outputs/hosp_nhsregion.rds")
cumHospAd_age_nhsregion <- readRDS("../outputs/cumHospAd_age_nhsregion.rds")

## read in parameters, remove guff and reorder
pars <- readRDS("../wave1/disease.rds") %>%
    rename(nu = `beta[1]`, nuA = `beta[6]`) %>%
    select(!c(starts_with("beta["), repeats)) %>%
    select(nu, nuA, !c(beta_scale, p_move, output), beta_scale, p_move) %>%
    as.data.frame()

## read in contact matrices
contact1 <- read_csv("../inputs/POLYMOD_matrix.csv", col_names = FALSE) %>%
    as.matrix()
contact2 <- read_csv("../inputs/coMix_matrix.csv", col_names = FALSE) %>%
    as.matrix()

## read in initial conditions
#u1 <- readRDS("../outputs/u1.rds")
u1_moves <- readRDS("../outputs/u1_moves.rds")

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
ageProbs <- read_csv("../inputs/age_seeds.csv", col_names = FALSE)$X2

## read in commuter data
EW19 <- read_delim("../inputs/EW19.dat", delim = " ", col_names = FALSE)

## expand to deal with age-classes
u1 <- apply(EW19, 1, function(x, ageProbs) {
        u <- matrix(0, 12, length(ageProbs))
        u[1, ] <- smart_round(ageProbs * x[3])
        list(u)
    }, ageProbs = ageProbs) %>%
    map(1) %>%
    abind(along = 3)

## read in player data
PM19 <- read_delim("../inputs/PlayMatrix19.dat", delim = " ", col_names = FALSE)
PlaySize19 <- read_delim("../inputs/PlaySize19.dat", delim = " ", col_names = FALSE)
  
## expand to deal with age-classes
u2 <- apply(PlaySize19, 1, function(x, ageProbs) {
        u <- matrix(0, 12, length(ageProbs))
        u[1, ] <- smart_round(ageProbs * x[2])
        list(u)
    }, ageProbs = ageProbs) %>%
    map(1) %>%
    abind(along = 3)
    
## load lookups
lookup <- readRDS("../outputs/lookup.rds")
age_lookup <- readRDS("../outputs/age_lookup.rds")

## set seed for reproducibility
set.seed(42)

## set number of days to run for
## and number of particles
tstart <- 0
tstop <- 50
npart <- 10

## set save folder
folder <- "saveOut"

## lookup table
class_lookup <- data.frame(var = c("S", "E", "A", "RA", "P", "I1", "DI", "I2", "RI", "H", "RH", "DH")) %>%
    mutate(class = 0:(n() - 1))
    
## set up inputs
u <- list(u1_moves = u1_moves, u1 = u1, u2_moves = as.matrix(PM19), u2 = u2)
    
## run model with model discrepancy
runs_md <- BPF(pars[100, ], C1 = contact1, C2 = contact2, lockdown_day = 20,
    cumDeath_lad = cumDeath_lad, cumDeath_age_region = cumDeath_age_region, 
    hosp_nhsregion = hosp_nhsregion, cumHospAd_age_nhsregion = cumHospAd_age_nhsregion, 
    lookup = lookup, age_lookup = age_lookup, u = u, tstart = tstart, tstop = tstop, npart = npart, 
    a1 = 0, a2 = 0, b = 0.1, a_dis = 0.05, b_dis = 0.05, niter = 10,
    sigma2_lad = 1, sigma2_age_region = 1, sigma2_nhsregion = 1, sigma2_age_nhsregion = 1,
    saveAll = TRUE, writeExt = TRUE, snapshot = TRUE)
    
## run continuation if required
cont <- TRUE
if(cont) {
    ## continuation method if required
    tstart <- 50
    tstop <- 75

    ## set up inputs
    u1_moves <- read_csv(paste0(folder, "/snapshot_u1moves_t", tstart, ".csv"), col_names = FALSE) %>%
        as.matrix()
    colnames(u1_moves) <- NULL
    ncohorts1 <- read_csv(paste0(folder, "/snapshot_ncohorts1_t", tstart, ".csv"), col_names = FALSE)$X1 %>%
        as.vector()
    ncohorts2 <- read_csv(paste0(folder, "/snapshot_ncohorts2_t", tstart, ".csv"), col_names = FALSE)$X1 %>%
        as.vector()
    playprobs <- read_csv(paste0(folder, "/snapshot_playprobs_t", tstart, ".csv"), col_names = FALSE)$X1 %>%
        as.vector()

    ## read in snapshots
    u1 <- list()
    u2 <- list()
    for(i in 1:npart) {
        u1[[i]] <- fread(paste0(folder, "/snapshot_u1_t", tstart, "_", i - 1, ".csv"))
        udims <- c(max(u1[[i]]$class) + 1, ncol(u1[[i]]) - 2, max(u1[[i]]$cohort))
        u1[[i]][, class := NULL]
        u1[[i]][, cohort := NULL]
        u1[[i]] <- as.matrix(u1[[i]])
        u1[[i]] <- aperm(array(t(u1[[i]]), c(udims[2], udims[1], udims[3])), c(2, 1, 3))
            
        u2[[i]] <- fread(paste0(folder, "/snapshot_u2_t", tstart, "_", i - 1, ".csv"))
        udims <- c(max(u2[[i]]$class) + 1, ncol(u2[[i]]) - 2, max(u2[[i]]$lad))
        u2[[i]][, class := NULL]
        u2[[i]][, lad := NULL]
        u2[[i]] <- as.matrix(u2[[i]])
        u2[[i]] <- aperm(array(t(u2[[i]]), c(udims[2], udims[1], udims[3])), c(2, 1, 3))
    }

    u <- list(u1_moves = u1_moves, u1 = u1, u2 = u2, ncohorts1 = ncohorts1, ncohorts2 = ncohorts2, playprobs = playprobs)
        
    ## run model with model discrepancy
    runs_md <- BPF(pars[100, ], C1 = contact1, C2 = contact2, lockdown_day = 20,
        lookup = lookup, age_lookup = age_lookup, u = u, tstart = tstart, tstop = tstop, 
        npart = npart, a1 = 0, a2 = 0, b = 0.1, a_dis = 0.05, b_dis = 0.05, 
        sigma2_lad = 1, sigma2_age_region = 1, sigma2_nhsregion = 1, sigma2_age_nhsregion = 1,
        saveAll = TRUE, writeExt = TRUE, PF = FALSE)   
} 
    
###############################################
#######        LAD-level truth          #######
###############################################

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

## extract file names
files <- list.files(folder)
files <- files[grep("p_[0-9]*.csv", files)]
files <- paste0(folder, "/", files)
if(cont) {
    ## extract file names
    files1 <- list.files(paste0(folder, "_cont"))
    files1 <- files1[grep("p_[0-9]*.csv", files1)]
    files1 <- paste0(folder, "_cont/", files1)
    files <- c(files, files1)
}

sims_md <- map(files, read.csv, header = TRUE)

## load in runs and group at the national level
sims_md <- map(files, function(y) {
        read.csv(y, header = TRUE) %>%
        group_by(time, class) %>%
        summarise(across(!lad, sum), .groups = "drop")
    })
names(sims_md) <- gsub("^.*_([0-9]*).csv", "\\1", files)    
sims_md <- bind_rows(sims_md, .id = "particle") %>%    
    inner_join(class_lookup, by = "class") %>%
    select(!class) %>%
    rename(t = time) %>%
    pivot_longer(!c(particle, t, var), names_to = "age", values_to = "n") %>%
    mutate(age = as.numeric(gsub("age", "", age))) %>%
    group_by(t, var, age) %>%
    summarise(
        LCI = quantile(n, probs = 0.025),
        LQ = quantile(n, probs = 0.25),
        median = median(n),
        UQ = quantile(n, probs = 0.75),
        UCI = quantile(n, probs = 0.975),
        .groups = "drop"
    )
    
p1 <- list()     
p1[[1]] <- ggplot(sims_md, aes(x = t)) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5) +
    geom_ribbon(aes(ymin = LQ, ymax = UQ), alpha = 0.5) +
    geom_line(aes(y = median)) +
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
data <- filter(cumDeath_lad, t <= tstop) %>%
    pivot_longer(!t, names_to = "lad", values_to = "n") %>%
    mutate(lad = as.numeric(gsub("deaths_", "", lad))) %>%
    group_by(t) %>%
    summarise(n = sum(n), .groups = "drop")

## extract file names
files <- list.files(folder)
files <- files[grep("p_lads_[0-9]*.csv", files)]
files <- paste0(folder, "/", files)
if(cont) {
    ## extract file names
    files1 <- list.files(paste0(folder, "_cont"))
    files1 <- files1[grep("p_lads_[0-9]*.csv", files1)]
    files1 <- paste0(folder, "_cont/", files1)
    files <- c(files, files1)
}

## load in runs and group at the national level
sims_md <- map(files, read.csv, header = TRUE)
names(sims_md) <- gsub("^.*_([0-9]*).csv", "\\1", files)
sims_md <- bind_rows(sims_md, .id = "particle") %>%
    arrange(particle, time, lad) %>%
    group_by(particle, lad) %>%
    mutate(deaths = cumsum(deaths)) %>%
    group_by(particle, time) %>%
    summarise(deaths = sum(deaths), .groups = "drop") %>%
    rename(t = time, n = deaths) %>%
    group_by(t) %>%
    summarise(
        LCI = quantile(n, probs = 0.025),
        LQ = quantile(n, probs = 0.25),
        median = median(n),
        UQ = quantile(n, probs = 0.75),
        UCI = quantile(n, probs = 0.975),
        .groups = "drop"
    )
     
p1[[2]] <- ggplot(sims_md, aes(x = t)) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5) +
    geom_ribbon(aes(ymin = LQ, ymax = UQ), alpha = 0.5) +
    geom_line(aes(y = median)) +
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
data <- filter(cumDeath_age_region, t <= tstop) %>%
    pivot_longer(!t, names_to = "var", values_to = "n") %>%
    mutate(age = gsub('^(?:[^_]*_)(.*)', '\\1', var)) %>%
    mutate(region = gsub('^(?:[^_]*_)(.*)', '\\1', age)) %>%
    mutate(age = gsub('(.*)_[0-9]*', '\\1', age)) %>%
    mutate(across(c(age, region), as.numeric)) %>%
    inner_join(region_lookup, by = c("region" = "FID")) %>%
    select(t, n, age, RGN19NM)
    
## extract file names
files <- list.files(folder)
files <- files[grep("p_age_region_[0-9]*.csv", files)]
files <- paste0(folder, "/", files)
if(cont) {
    ## extract file names
    files1 <- list.files(paste0(folder, "_cont"))
    files1 <- files1[grep("p_age_region_[0-9]*.csv", files1)]
    files1 <- paste0(folder, "_cont/", files1)
    files <- c(files, files1)
}

## load in runs and group at the national level
sims_md <- map(files, read.csv, header = TRUE)
names(sims_md) <- gsub("^.*_([0-9]*).csv", "\\1", files)
sims_md <- bind_rows(sims_md, .id = "particle") %>%
    pivot_longer(!c(particle, time, region), names_to = "age", values_to = "n") %>%
    mutate(age = as.numeric(gsub("age", "", age))) %>%
    rename(t = time) %>%
    arrange(particle, region, age, t) %>%
    group_by(particle, region, age) %>%
    mutate(n = cumsum(n)) %>%
    ungroup() %>%
    inner_join(region_lookup, by = c("region" = "FID")) %>%
    select(t, n, age, RGN19NM) %>%
    group_by(t, age, RGN19NM) %>%
    summarise(
        LCI = quantile(n, probs = 0.025),
        LQ = quantile(n, probs = 0.25),
        median = median(n),
        UQ = quantile(n, probs = 0.75),
        UCI = quantile(n, probs = 0.975),
        .groups = "drop"
    )
    
p1[[3]] <- ggplot(sims_md, aes(x = t)) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5) +
    geom_ribbon(aes(ymin = LQ, ymax = UQ), alpha = 0.5) +
    geom_line(aes(y = median)) +
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
data <- filter(hosp_nhsregion, t <= tstop) %>%
    pivot_longer(!t, names_to = "region", values_to = "n") %>%
    mutate(region = as.numeric(gsub("hosp_", "", region))) %>%
    inner_join(nhsregion_lookup, by = c("region" = "FID")) %>%
    select(t, n, areaName)

## extract file names
files <- list.files(folder)
files <- files[grep("p_nhsregion_[0-9]*.csv", files)]
files <- paste0(folder, "/", files)
if(cont) {
    ## extract file names
    files1 <- list.files(paste0(folder, "_cont"))
    files1 <- files1[grep("p_nhsregion_[0-9]*.csv", files1)]
    files1 <- paste0(folder, "_cont/", files1)
    files <- c(files, files1)
}

## load in runs and group at the national level
sims_md <- map(files, read.csv, header = TRUE)
names(sims_md) <- gsub("^.*_([0-9]*).csv", "\\1", files)

## load in runs and group at the national level
sims_md <- bind_rows(sims_md, .id = "particle") %>%
    rename(t = time, n = hosp) %>%
    inner_join(nhsregion_lookup, by = c("nhsregion" = "FID")) %>%
    select(t, n, areaName) %>%
    arrange(areaName, t) %>%
    group_by(t, areaName) %>%
    summarise(
        LCI = quantile(n, probs = 0.025),
        LQ = quantile(n, probs = 0.25),
        median = median(n),
        UQ = quantile(n, probs = 0.75),
        UCI = quantile(n, probs = 0.975),
        .groups = "drop"
    )
     
p1[[4]] <- ggplot(sims_md, aes(x = t)) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5) +
    geom_ribbon(aes(ymin = LQ, ymax = UQ), alpha = 0.5) +
    geom_line(aes(y = median)) +
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
data <- filter(cumHospAd_age_nhsregion, t <= tstop) %>%
    pivot_longer(!t, names_to = "var", values_to = "n") %>%
    mutate(age = gsub('^(?:[^_]*_)(.*)', '\\1', var)) %>%
    mutate(region = gsub('^(?:[^_]*_)(.*)', '\\1', age)) %>%
    mutate(age = gsub('(.*)_[0-9]*', '\\1', age)) %>%
    mutate(across(c(age, region), as.numeric)) %>%
    inner_join(nhsregion_lookup, by = c("region" = "FID")) %>%
    select(t, n, age, areaName)

## extract file names
files <- list.files(folder)
files <- files[grep("p_age_nhsregion_[0-9]*.csv", files)]
files <- paste0(folder, "/", files)
if(cont) {
    ## extract file names
    files1 <- list.files(paste0(folder, "_cont"))
    files1 <- files1[grep("p_age_nhsregion_[0-9]*.csv", files1)]
    files1 <- paste0(folder, "_cont/", files1)
    files <- c(files, files1)
}

## load in runs and group at the national level
sims_md <- map(files, read.csv, header = TRUE)
names(sims_md) <- gsub("^.*_([0-9]*).csv", "\\1", files)

## load in runs and group at the national level
sims_md <- bind_rows(sims_md, .id = "particle") %>%
    pivot_longer(!c(particle, time, nhsregion), names_to = "age", values_to = "n") %>%
    mutate(age = as.numeric(gsub("age", "", age))) %>%
    rename(t = time) %>%
    arrange(particle, nhsregion, age, t) %>%
    group_by(particle, nhsregion, age) %>%
    mutate(n = cumsum(n)) %>%
    ungroup() %>%
    inner_join(nhsregion_lookup, by = c("nhsregion" = "FID")) %>%
    select(particle, t, n, age, areaName) %>%
    group_by(t, age, areaName) %>%
    summarise(
        LCI = quantile(n, probs = 0.025),
        LQ = quantile(n, probs = 0.25),
        median = median(n),
        UQ = quantile(n, probs = 0.75),
        UCI = quantile(n, probs = 0.975),
        .groups = "drop"
    )
     
p1[[5]] <- ggplot(sims_md, aes(x = t)) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5) +
    geom_ribbon(aes(ymin = LQ, ymax = UQ), alpha = 0.5) +
    geom_line(aes(y = median)) +
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
ggsave("simsBPF.pdf", p1, width = 15, height = 15)

###############################################
#####          LAD-level plots            #####
###############################################

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
    mutate(age = as.numeric(age))
    
## extract LADs with largest epidemic load
top5 <- filter(data, var == "DI") %>%
    filter(t == max(t)) %>%
    group_by(lad) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    arrange(desc(n)) %>%
    slice(1:5) %>%
    select(-n)
    
## extract top 5 lads
data <- filter(data, lad %in% top5$lad) 

## extract file names
files <- list.files(folder)
files <- files[grep("p_[0-9]*.csv", files)]
files <- paste0(folder, "/", files)
if(cont) {
    ## extract file names
    files1 <- list.files(paste0(folder, "_cont"))
    files1 <- files1[grep("p_[0-9]*.csv", files1)]
    files1 <- paste0(folder, "_cont/", files1)
    files <- c(files, files1)
}

## load in runs and group at the national level
sims_md <- map(files, function(y, lads) {
        read.csv(y, header = TRUE) %>%
            filter(lad %in% lads)
    }, lads = top5$lad)
names(sims_md) <- gsub("^.*_([0-9]*).csv", "\\1", files)
sims_md <- bind_rows(sims_md, .id = "particle") %>% 
    inner_join(class_lookup, by = "class") %>%
    select(!class) %>%
    rename(t = time) %>%
    pivot_longer(!c(particle, t, var, lad), names_to = "age", values_to = "n") %>%
    mutate(age = as.numeric(gsub("age", "", age))) %>%
    mutate(var = gsub("one", "1", var)) %>%
    mutate(var = gsub("two", "2", var)) %>%
    group_by(t, var, age, lad) %>%
    summarise(
        LCI = quantile(n, probs = 0.025),
        LQ = quantile(n, probs = 0.25),
        median = median(n),
        UQ = quantile(n, probs = 0.75),
        UCI = quantile(n, probs = 0.975),
        .groups = "drop"
    ) %>%
    mutate(lad = as.character(lad))
    
p1 <- ggplot(sims_md, aes(x = t, colour = lad, fill = lad)) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), colour = NA, alpha = 0.5) +
    geom_ribbon(aes(ymin = LQ, ymax = UQ), colour = NA, alpha = 0.5) +
    geom_line(aes(y = median)) +
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

ggsave("simsTopLADsBPF.pdf", p1, width = 10, height = 10)

