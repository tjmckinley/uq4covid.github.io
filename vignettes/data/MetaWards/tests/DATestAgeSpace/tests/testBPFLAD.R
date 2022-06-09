## load libraries
library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(parallel)
library(abind)
library(sitmo)
library(patchwork)
library(BH)
Sys.setenv("PKG_LIBS" = "-lgmp")

## source Rcpp PF code
sourceCpp("../BPF.cpp")

## source function to run PF and return log-likelihood
source("../BPF.R")

## read in simulated data and generate incidence curves
data <- readRDS("../outputs/disSims.rds")

## read in parameters, remove guff and reorder
pars <- readRDS("../wave1/disease.rds") %>%
    rename(nu = `beta[1]`, nuA = `beta[6]`) %>%
    select(!c(starts_with("beta"), repeats)) %>%
    select(nu, nuA, !output) %>%
    as.data.frame()

## read in contact matrix
contact <- read_csv("../inputs/POLYMOD_matrix.csv", col_names = FALSE) %>%
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

## set seed for reproducibility
set.seed(42)

## set up plot data
plot_data <- pivot_longer(filter(data, t <= 100), !t, names_to = "var", values_to = "n") %>%
    mutate(age = gsub('^(?:[^_]*_)(.*)', '\\1', var)) %>%
    mutate(LAD = gsub('^(?:[^_]*_)(.*)', '\\1', age)) %>%
    mutate(age = gsub('(.*)_[0-9]*', '\\1', age)) %>%
    mutate(obs = grepl("obs", var)) %>%
    mutate(var = gsub('^(.*)_[0-9]*_.*', '\\1', var)) %>%
    mutate(var = gsub("one", "1", var)) %>%
    mutate(var = gsub("two", "2", var)) %>%
    mutate(var = gsub("obs", "", var))
    
## run model with model discrepancy
runs_md <- BPF(pars[6, ], C = contact, data = data, u1_moves = u1_moves,
    u1 = u1, ndays = 100, npart = 10, a_dis = 0.05, b_dis = 0.05,
    a1 = 0.01, a2 = 0.2, b = 0.001, saveAll = TRUE, writeExt = TRUE)

## extract file names
folder <- "saveOut"
files <- list.files(folder)

## lookup table
lookup <- data.frame(var = c("S", "E", "A", "RA", "P", "I1", "DI", "I2", "RI", "H", "RH", "DH", "DIobs", "DHobs")) %>%
    mutate(class = 0:(n() - 1))
    
## aggregate data
p <- filter(plot_data, !obs) %>%
    group_by(t, var, age) %>%
    summarise(n = sum(n), .groups = "drop")

## aggregate observed data   
pobs <- filter(plot_data, obs) %>%
    group_by(t, var, age) %>%
    summarise(n = sum(n), .groups = "drop")
    
## plot particle estimates of states at the national level
sims_md <- map(files, function(y, folder, lookup) {
        read.csv(paste0(folder, "/", y), header = TRUE) %>%
            group_by(time, class) %>%
            summarise(across(!lad, sum), .groups = "drop") %>%
            inner_join(lookup, by = "class") %>%
            select(!class) %>%
            rename(t = time) %>%
            pivot_longer(!c(t, var), names_to = "age", values_to = "n") %>%
            mutate(age = as.numeric(gsub("age", "", age)))
    }, folder = folder, lookup = lookup) %>%
    bind_rows(.id = "particle")

## summarise simulations
sims_obs <- filter(sims_md, var %in% paste0(unique(pobs$var), "obs")) %>%
    group_by(var, age, t) %>%
    summarise(
        LCI = quantile(n, probs = 0.025),
        LQ = quantile(n, probs = 0.25),
        Median = quantile(n, probs = 0.5),
        UQ = quantile(n, probs = 0.75),
        UCI = quantile(n, probs = 0.975),
        .groups = "drop"
    ) %>%
    mutate(var = gsub("obs", "", var))
    
## summarise simulations
sims_md <- filter(sims_md, !(var %in% paste0(unique(pobs$var), "obs"))) %>%
    group_by(var, age, t) %>%
    summarise(
        LCI = quantile(n, probs = 0.025),
        LQ = quantile(n, probs = 0.25),
        Median = quantile(n, probs = 0.5),
        UQ = quantile(n, probs = 0.75),
        UCI = quantile(n, probs = 0.975),
        .groups = "drop"
    )

## produce particle trajectories plots
p1 <- list()
p1[[1]] <- ggplot(sims_md, aes(x = t)) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.3, colour = NA) +
    geom_ribbon(aes(ymin = LQ, ymax = UQ), alpha = 0.3, colour = NA) +
    geom_line(aes(y = Median)) +
    geom_line(aes(y = n), data = p, linetype = "dashed") +
    facet_grid(var ~ age, scales = "free") +
    xlab("Days") + 
    ylab("Counts") +
    ggtitle("Truth")
p1[[2]] <- ggplot(sims_obs, aes(x = t)) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.3, colour = NA) +
    geom_ribbon(aes(ymin = LQ, ymax = UQ), alpha = 0.3, colour = NA) +
    geom_line(aes(y = Median)) +
    geom_line(aes(y = n), data = pobs, linetype = "dashed") +
    facet_grid(var ~ age, scales = "free") +
    xlab("Days") + 
    ylab("Counts") +
    ggtitle("Observed")
p1 <- wrap_plots(p1, nrow = 2, heights = c(0.8, 0.2))
ggsave("simsBPF.pdf", p1, width = 10, height = 10)

## extract LADs with largest epidemic load at day 100
p <- select(data, t, !contains("obs") & starts_with("DI")) %>%
    pivot_longer(!t, names_to = "var", values_to = "n") %>%
    mutate(age = gsub('^(?:[^_]*_)(.*)', '\\1', var)) %>%
    mutate(LAD = gsub('^(?:[^_]*_)(.*)', '\\1', age)) %>%
    mutate(age = gsub('(.*)_[0-9]*', '\\1', age)) %>%
    mutate(var = gsub('^(.*)_[0-9]*_.*', '\\1', var)) %>%
    filter(t == 100) %>%
    group_by(LAD) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    arrange(desc(n)) %>%
    slice(1:5) %>%
    select(-n)
    
## extract observed data
pobs <- inner_join(p,
    filter(plot_data, obs) %>%
        mutate(LAD = gsub("obs", "", LAD)),
    by = "LAD"
)

## plot particle estimates of states at the LAD level
sims_md <- map(files, function(y, folder, lookup, toplads) {
        read.csv(paste0(folder, "/", y), header = TRUE) %>%
            inner_join(toplads, by = c("lad" = "LAD")) %>%
            inner_join(lookup, by = "class") %>%
            select(!class) %>%
            rename(t = time) %>%
            pivot_longer(!c(t, var, lad), names_to = "age", values_to = "n") %>%
            mutate(age = as.numeric(gsub("age", "", age))) %>%
            rename(LAD = lad) %>%
            mutate(LAD = as.character(LAD))
    }, folder = folder, lookup = lookup, toplads = mutate(p, LAD = as.numeric(LAD))) %>%
    bind_rows(.id = "particle")

## extract simulations for observed states
sims_obs <- inner_join(p, sims_md, by = "LAD") %>%
    filter(var %in% paste0(unique(pobs$var), "obs")) %>%
    group_by(var, LAD, age, t) %>%
    summarise(
        LCI = quantile(n, probs = 0.025),
        LQ = quantile(n, probs = 0.25),
        Median = quantile(n, probs = 0.5),
        UQ = quantile(n, probs = 0.75),
        UCI = quantile(n, probs = 0.975),
        .groups = "drop"
    ) %>%
    mutate(var = gsub("obs", "", var))
    
## extract simulations for all states
sims_md <- inner_join(p, sims_md, by = "LAD") %>%
    filter(!(var %in% paste0(unique(pobs$var), "obs"))) %>%
    group_by(var, LAD, age, t) %>%
    summarise(
        LCI = quantile(n, probs = 0.025),
        LQ = quantile(n, probs = 0.25),
        Median = quantile(n, probs = 0.5),
        UQ = quantile(n, probs = 0.75),
        UCI = quantile(n, probs = 0.975),
        .groups = "drop"
    )

## data for all states    
p <- inner_join(p, select(plot_data, !obs), by = "LAD")

## plot of all states against simulations
p1 <- list()
p1[[1]] <- ggplot(sims_md, aes(x = t)) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = LAD), alpha = 0.3, colour = NA) +
    geom_ribbon(aes(ymin = LQ, ymax = UQ, fill = LAD), alpha = 0.3, colour = NA) +
    geom_line(aes(y = Median, colour = LAD)) +
    geom_line(aes(y = n, colour = LAD), data = p, linetype = "dashed") +
    facet_grid(var ~ age, scales = "free") +
    xlab("Days") + 
    ylab("Counts") +
    ggtitle("Truth")
p1[[2]] <- ggplot(sims_obs, aes(x = t)) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = LAD), alpha = 0.3, colour = NA) +
    geom_ribbon(aes(ymin = LQ, ymax = UQ, fill = LAD), alpha = 0.3, colour = NA) +
    geom_line(aes(y = Median, colour = LAD)) +
    geom_line(aes(y = n, colour = LAD), data = pobs, linetype = "dashed") +
    facet_grid(var ~ age, scales = "free") +
    xlab("Days") + 
    ylab("Counts") +
    ggtitle("Observed")
p1 <- wrap_plots(p1, nrow = 2, heights = c(0.8, 0.2)) +
    plot_layout(guides = "collect")
ggsave("simsTopLADsBPF.pdf", p1, width = 10, height = 10)


