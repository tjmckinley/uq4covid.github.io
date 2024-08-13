## load libraries
library(lhs)
library(mclust)
library(tidyverse)

## check if being run in batch mode
args <- commandArgs(TRUE)
if(length(args) != 0) {
    ## extract command line arguments
    args <- commandArgs(TRUE)
    if(length(args) > 0) {
        stopifnot(length(args) == 1)
        wave <- args[1]
    } else {
        stop("No arguments")
    }
} else {
    wave <- "1Feb"
}

## create directory to save samples
if(dir.exists(paste0("wave", wave))) {
    stop("wave", wave, " directory already exists. Please delete if you wish to overwrite.")
}
dir.create(paste0("wave", wave))

## set case specific values for simulations
tstart <- 0
tstop <- 37
lockdown_day <- 37
npart <- 50
niter <- 10
a1 <- 0.001
a2 <- 0.001
b1 <- 0.025
b2 <- 0.025
a_dis <- 0.05
b_dis <- 0.01
b_dis_8 <- 0.5
sigma2_lad <- 0.00295858
sigma2_age_region <- 0.01388889
sigma2_nhsregion <- 0.14285714
sigma2_age_nhsregion <- 0.03571429
saveAll <- TRUE
snapshot <- TRUE
writeExt <- TRUE

## write to file
writeLines(as.character(c(tstart, tstop, lockdown_day, npart, niter, a1, a2,
    b1, b2, a_dis, b_dis, b_dis_8, sigma2_lad, sigma2_age_region, sigma2_nhsregion, 
    sigma2_age_nhsregion, saveAll, snapshot, writeExt)), paste0("wave", wave, "/fixedInputs.txt"))

## source dataTools
source("inputs/dataTools.R")

## set up parameter ranges for uniform ranges
parRanges <- data.frame(
    parameter = c("R0", "TE", "TP", "TI1", "TI2", "nuA", "beta_scale", "p_move", "MD_scale", "MD_time_EM", "MD_time_EE", "MD_time_L", "MD_time_NE", "MD_time_NW", "MD_time_SE", "MD_time_SW", "MD_time_WM", "MD_time_YH"),
    lower = c(2, 0.1, 1.2, 2.8, 0.0001, 0, 0, 0, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    upper = c(4.5, 2, 3, 4.5, 0.5, 1, 1, 1, 1, 37, 37, 37, 37, 37, 37, 37, 37, 37),
    stringsAsFactors = FALSE
) 

## read in contact matrix to use for NGM
C <- as.matrix(read.csv("inputs/POLYMOD_matrix.csv", header = FALSE))

## generate LHS design (200 + 50 validation)
ndesign <- 200
design <- randomLHS(ndesign, nrow(parRanges))
nval <- 50
design <- rbind(design, randomLHS(nval, nrow(parRanges)))
colnames(design) <- parRanges$parameter
design <- as_tibble(design)

## convert to input space
inputs <- convertDesignToInput(design, parRanges, "zero_one")

## generate space-filling design for other parameters

## load FMM objects
hospStays <- readRDS("inputs/hospStays.rds")
pathways <- readRDS("inputs/pathways.rds")
hospThresh <- readRDS("inputs/hospThresh.rds")
pathThresh <- readRDS("inputs/pathThresh.rds")

## generate design points for hospital stay lengths
hospStaysInput <- FMMmaximin(
    hospStays, 
    ndesign + 20, 
    matrix(c(-Inf, Inf, 0, Inf), ncol = 2, byrow = TRUE)
) %>%
    as_tibble() %>%
    rename(alphaTH = x1, etaTH = x2)
    
## check against prior density restrictions
hospStaysInput <- hospStaysInput[dens(as.matrix(hospStaysInput), hospStays$modelName, hospStays$parameters, logarithm = TRUE) > hospThresh, ]
if(nrow(hospStaysInput) < ndesign) stop("Can't generate enough valid hospital points")
hospStaysInput <- hospStaysInput[1:ndesign, ]

## add validation points
hospStaysVal <- FMMmaximin(
    hospStays, 
    nval + 20, 
    matrix(c(-Inf, Inf, 0, Inf), ncol = 2, byrow = TRUE)
) %>%
    as_tibble() %>%
    rename(alphaTH = x1, etaTH = x2)

## check against prior density restrictions  
hospStaysVal <- hospStaysVal[dens(as.matrix(hospStaysVal), hospStays$modelName, hospStays$parameters, logarithm = TRUE) > hospThresh, ]
if(nrow(hospStaysVal) < nval) stop("Can't generate enough valid hospital points")
hospStaysVal <- hospStaysVal[1:nval, ]

## bind together
hospStaysInput <- rbind(hospStaysInput, hospStaysVal)

## generate design points for other transition probabilities

## this function checks validity of inputs
pathwaysLimitFn <- function(x, ages) {
    ## check all parameters give valid probabilities
    ## in (0, 1)
    singleProbs <- apply(x, 1, function(x, ages) {
        eta <- x[5]
        alphas <- x[c(1, 4)]
        y <- sapply(alphas, function(a, eta, ages) {
            y <- exp(a + eta * ages)
            all(y >= 0 & y <= 1)
        }, eta = eta, ages = ages)
        y <- all(y)
        eta <- x[5] * x[6]
        a <- x[2]
        y1 <- exp(a + eta * ages)
        y1 <- all(y1 >= 0 & y1 <= 1)
        eta <- x[5] * x[7]
        a <- x[3]
        y2 <- exp(a + eta * ages)
        y2 <- all(y2 >= 0 & y2 <= 1)
        y & y1 & y2
    }, ages = ages)
    ## check multinomial probabilities sum to one
    multiProbs <- apply(x[, -c(1, 3)], 1, function(x, ages) {
        alphaI1D <- x[1]
        alphaI1H <- x[2]
        eta <- x[3]
        eta_scale <- x[4]
        pI1D <- exp(alphaI1D + eta * eta_scale * ages)
        pI1H <- exp(alphaI1H + eta * ages)
        p <- pI1D + pI1H
        all(p >= 0 & p <= 1)
    }, ages = ages)
    multiProbs & singleProbs
}

## produces design points subject to constraints
pathwaysInput <- FMMmaximin(
    pathways, 
    ndesign + 20,
    matrix(c(rep(c(-20, 0), times = 4), 0, 1), ncol = 2, byrow = TRUE),
    pathwaysLimitFn,
    eta_scale = c(0.5, 2),
    ages = c(2.5, 11, 23.5, 34.5, 44.5, 55.5, 65.5, 75.5)
) %>%
    as_tibble() %>%
    rename(alphaEP = x1, alphaI1D = x2, alphaHD = x3, alphaI1H = x4, eta = x5, etaI_scale = x6, etaH_scale = x7)
    
## check against prior density restrictions
pathwaysInput <- pathwaysInput[dens(as.matrix(pathwaysInput)[, -c(6, 7)], pathways$modelName, pathways$parameters, logarithm = TRUE) > pathThresh, ]
pathwaysInput <- pathwaysInput[pathwaysInput[, 6] > 0.5 & pathwaysInput[, 6] < 2, ]
pathwaysInput <- pathwaysInput[pathwaysInput[, 7] > 0.5 & pathwaysInput[, 7] < 2, ]
if(nrow(pathwaysInput) < ndesign) stop("Can't generate enough valid pathways points")
pathwaysInput <- pathwaysInput[1:ndesign, ]
    
## add validation points
pathwaysVal <- FMMmaximin(
    pathways, 
    nval + 20,
    matrix(c(rep(c(-20, 0), times = 4), 0, 1), ncol = 2, byrow = TRUE),
    pathwaysLimitFn,
    eta_scale = c(0.5, 2),
    ages = c(2.5, 11, 23.5, 34.5, 44.5, 55.5, 65.5, 75.5)
) %>%
    as_tibble() %>%
    rename(alphaEP = x1, alphaI1D = x2, alphaHD = x3, alphaI1H = x4, eta = x5, etaI_scale = x6, etaH_scale = x7)
    
## check against prior density restrictions
pathwaysVal <- pathwaysVal[dens(as.matrix(pathwaysVal)[, -c(6, 7)], pathways$modelName, pathways$parameters, logarithm = TRUE) > pathThresh, ]
pathwaysVal <- pathwaysVal[pathwaysVal[, 6] > 0.5 & pathwaysVal[, 6] < 2, ]
pathwaysVal <- pathwaysVal[pathwaysVal[, 7] > 0.5 & pathwaysVal[, 7] < 2, ]
if(nrow(pathwaysVal) < nval) stop("Can't generate enough valid pathways points")
pathwaysVal <- pathwaysVal[1:nval, ]

## bind together
pathwaysInput <- rbind(pathwaysInput, pathwaysVal)

## bind to design
inputs <- cbind(inputs, hospStaysInput, pathwaysInput)

## add unique hash identifier
## (at the moment don't use "a0" type ensembleID, because MetaWards
## parses to dates)
inputs$output <- ensembleIDGen(ensembleID = "Ens1", nrow(inputs))

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
N <- smart_round(read_csv("inputs/age_seeds.csv", col_names = FALSE)$X2 * 56082077)
S0 <- N - smart_round(read_csv("inputs/age_seeds.csv", col_names = FALSE)$X2 * 100)
ages <- c(2.5, 11, 23.5, 34.5, 44.5, 55.5, 65.5, 75.5)

## convert input to disease
disease <- convertInputToDisease(inputs, C, N, S0, ages)
stopifnot(nrow(disease) == ndesign + nval)

## reorder samples
inputs <- arrange(inputs, output)
disease <- arrange(disease, output)

## save samples
saveRDS(inputs, paste0("wave", wave, "/inputs.rds"))
saveRDS(disease, paste0("wave", wave, "/disease.rds"))

## plot inputs
library(GGally)
p <- select(inputs, -output) %>%
    ggpairs(upper = "blank")
ggsave(paste0("wave", wave, "/design.pdf"), p, width = 10, height = 10)

