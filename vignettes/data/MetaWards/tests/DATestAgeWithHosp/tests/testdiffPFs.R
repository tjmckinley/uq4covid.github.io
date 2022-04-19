## load libraries
library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(parallel)
library(patchwork)

## source simulation model
sourceCpp("../discreteStochModel.cpp")

## read in truncated Skellam sampler
source("../trSkellam.R")

## source function to run PF and return log-likelihood
source("../PF.R")
source("../PF1.R")

## read in simulated data and generate incidence curves
data <- readRDS("../outputs/disSims.rds")

## read in parameters, remove guff and reorder
pars <- readRDS("../wave1/disease.rds") %>%
    rename(nu = `beta[1]`, nuA = `beta[6]`) %>%
    select(!c(starts_with("beta"), repeats)) %>%
    select(nu, nuA, !output)

## read in contact matrix
contact <- read_csv("../inputs/POLYMOD_matrix.csv", col_names = FALSE) %>%
    as.matrix()

## solution to round numbers preserving sum
## adapted from:
## https://stackoverflow.com/questions/32544646/round-vector-of-numerics-to-integer-while-preserving-their-sum
smart_round <- function(x) {
    y <- floor(x)
    indices <- tail(order(x - y), round(sum(x)) - sum(y))
    y[indices] <- y[indices] + 1
    y
}

## set up number of initial individuals in each age-class
N <- 10000
N <- smart_round(read_csv("../inputs/age_seeds.csv", col_names = FALSE)$X2 * N)
I0 <- smart_round(read_csv("../inputs/age_seeds.csv", col_names = FALSE)$X2 * 1)
S0 <- N - I0

## set initial counts
u <- matrix(0, 12, 8)
u[1, ] <- S0
u[2, ] <- I0

## set seed for reproducibility
set.seed(666)

## set number of replicates
nreps <- 20

## extract parameters
pars1 <- slice(pars, 1:9) %>%
    mutate(id = 1:n())

## expand to replicates
pars <- pars1
if(nreps > 1) {
    for(i in 2:nreps) {
        pars <- rbind(pars, pars1)
    }
}
rm(pars1)
pars <- arrange(pars, id)

## set number of particles
npart <- 10

runs <- list()
for(j in 1:length(npart)) {
    ## repeat but adding some model discrepancy
    runs_md <- PF(select(pars, -id), C = contact, data = data, u = u, ndays = 20, 
        npart = npart[j], MD = TRUE, a1 = 0.01, a2 = 0.2, b = 0.1,
        a_dis = 0.05, b_dis = 0.05, saveAll = NA)

    ## repeat but adding some model discrepancy
    runs_mdt100 <- PF(select(pars, -id), C = contact, data = data, u = u, ndays = 20, 
        npart = npart[j] * 100, MD = TRUE, a1 = 0.01, a2 = 0.2, b = 0.1,
        a_dis = 0.05, b_dis = 0.05, saveAll = NA)
    
    ## repeat but adding some model discrepancy using new updating scheme
    runs_md1 <- PF1(select(pars, -id), C = contact, data = data, u = u, ndays = 20, 
        npart = npart[j], MD = TRUE, a1 = 0.01, a2 = 0.2, b = 0.1,
        a_dis = 0.05, b_dis = 0.05, saveAll = NA)
    
    ## repeat but adding some model discrepancy using new updating scheme
    runs_md1t100 <- PF1(select(pars, -id), C = contact, data = data, u = u, ndays = 20, 
        npart = npart[j] * 100, MD = TRUE, a1 = 0.01, a2 = 0.2, b = 0.1,
        a_dis = 0.05, b_dis = 0.05, saveAll = NA)
    
    ## collapse to data frame and plot
    runs[[j]] <- tibble(MD = runs_md, MDt100 = runs_mdt100, newMD = runs_md1, newMDt100 = runs_md1t100) %>%
        cbind(select(pars, id)) %>%
        pivot_longer(!id, names_to = "Type")
}
runs <- bind_rows(runs, .id = "npart1") %>%
    mutate(npart1 = npart[as.numeric(npart1)]) %>%
    rename(npart = npart1) %>%
    mutate(npart = as.character(npart))

## density plot of PF estimates
p <- ggplot(runs) +
    geom_density(aes(x = value, fill = Type), alpha = 0.5) +
    facet_wrap(~id) +
    ggtitle(paste0("Num. of replicates = ", nreps, " Num. of particles = ", npart)) +
    ylab("Density")
ggsave("PF_dens_ll.pdf", p, width = 7, height = 7)
