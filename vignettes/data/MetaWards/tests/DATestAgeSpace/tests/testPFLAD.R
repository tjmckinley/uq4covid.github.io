## load libraries
library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(parallel)
library(abind)
library(sitmo)
library(patchwork)

## source Rcpp PF code
sourceCpp("../PF.cpp")

## source function to run PF and return log-likelihood
source("../PF.R")

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
u1 <- readRDS("../outputs/u1.rds")
u1_moves <- readRDS("../outputs/u1_moves.rds")

## set seed for reproducibility
set.seed(666)

## set up plot data
plot_data <- pivot_longer(filter(data, t <= 100), !t, names_to = "var", values_to = "n") %>%
    mutate(age = gsub('^(?:[^_]*_)(.*)', '\\1', var)) %>%
    mutate(LAD = gsub('^(?:[^_]*_)(.*)', '\\1', age)) %>%
    mutate(age = gsub('(.*)_[0-9]*', '\\1', age)) %>%
    mutate(obs = grepl("obs", var)) %>%
    mutate(var = gsub('^(.*)_[0-9]*_.*', '\\1', var)) %>%
    mutate(var = gsub("one", "1", var)) %>%
    mutate(var = gsub("two", "2", var)) %>%
    mutate(age = gsub("obs", "", age))
    
## run model with model discrepancy
runs_md <- PF(pars[6, ], C = contact, data = data, u1_moves = u1_moves,
    u1 = u1, ndays = 100, npart = 20, MD = TRUE, a_dis = 1, b_dis = 1, 
    a1 = 1, a2 = 1, b = 1, saveAll = TRUE)

## plot particle estimates of states (unweighted)
sims_md <- map(runs_md$particles[[1]], ~{
        map(., function(y) {
            x <- matrix(y, prod(dim(y)[c(1, 3)]), dim(y)[2])
            colnames(x) <- paste0("age", 1:ncol(x))
            as_tibble(x) %>%
            mutate(var = rep(c("S", "E", "A", "RA", "P", "Ione", "DI", "Itwo", "RI", "H", "RH", "DH"), times = dim(y)[3])) %>%
            mutate(LAD = rep(1:dim(y)[3], each = dim(y)[1]))
        }) %>%
        bind_rows(.id = "particle")
    }) %>%
    bind_rows(.id = "t") %>%
    pivot_longer(!c(particle, t, var, LAD), names_to = "age", values_to = "n") %>%
    mutate(age = as.numeric(gsub("age", "", age))) %>%
    mutate(t = as.numeric(t) - 1) %>%
    mutate(LAD = as.character(LAD)) %>%
    mutate(var = gsub("one", "1", var)) %>%
    mutate(var = gsub("two", "2", var))

## extract LADs with largest epidemic load at day 100
p <- select(data, t, !ends_with("obs") & starts_with("DI")) %>%
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

## extract simulations for observed states
sims_obs <- inner_join(p, sims_md, by = "LAD") %>%
    filter(var %in% unique(pobs$var)) %>%
    group_by(var, LAD, age, t) %>%
    summarise(
        LCI = quantile(n, probs = 0.025),
        LQ = quantile(n, probs = 0.25),
        Median = quantile(n, probs = 0.5),
        UQ = quantile(n, probs = 0.75),
        UCI = quantile(n, probs = 0.975),
        .groups = "drop"
    )
    
## extract simulations for all states
sims_md <- inner_join(p, sims_md, by = "LAD") %>%
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
ggsave("simsTopLADs.pdf", p1, width = 10, height = 10)


