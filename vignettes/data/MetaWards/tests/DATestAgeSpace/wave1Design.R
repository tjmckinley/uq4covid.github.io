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
    wave <- "1"
}

## create directory to save samples
if(dir.exists(paste0("wave", wave))) {
    stop("Can't overwrite existing directory")
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
b_dis_london <- 0.01
sigma2_lad <- 0.00295858
sigma2_age_region <- 0.01388889
sigma2_nhsregion <- 0.14285714
sigma2_age_nhsregion <- 0.03571429
saveAll <- TRUE
snapshot <- TRUE
writeExt <- TRUE

## write to file
writeLines(as.character(c(tstart, tstop, lockdown_day, npart, niter, a1, a2,
    b1, b2, a_dis, b_dis, b_dis_london, sigma2_lad, sigma2_age_region, sigma2_nhsregion, 
    sigma2_age_nhsregion, saveAll, snapshot, writeExt)), paste0("wave", wave, "/fixedInputs.txt"))

## source dataTools
source("inputs/dataTools.R")

## set up parameter ranges for uniform ranges
parRanges <- data.frame(
    parameter = c("R0", "TE", "TP", "TI1", "TI2", "nuA", "alphaEP", "alphaI1D",
    "alphaI1H", "alphaHD", "etaEP", "etaI1D", "etaI1H", "etaHD", "beta_scale", "p_move", "MD_scale", "MD_time"),
    lower = c(2, 0.1, 1.2, 2.8, 0.0001, 0, -30, -30, -30, -30, 0, 0, 0, 0, 0, 0, 0.1, 0),
    upper = c(4.5, 2, 3, 4.5, 0.5, 1, 0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 1, 1, 1, 37),
    stringsAsFactors = FALSE
)

## read in contact matrix to use for NGM
C <- as.matrix(read.csv("inputs/POLYMOD_matrix.csv", header = FALSE))

## this function checks validity of inputs
pathwaysLimitFn <- function(x, ages) {
    ## check all parameters give valid probabilities
    ## in (0, 1)
    singleProbs <- apply(x, 1, function(x, ages) {
        etas <- x[5:8]
        alphas <- x[1:4]
        y <- map2_lgl(alphas, etas, function(a, e, ages) {
            y <- exp(a + e * ages)
            all(y >= 0 & y <= 1)
        }, ages = ages)
        all(y)
    }, ages = ages)
    ## check multinomial probabilities sum to one
    multiProbs <- apply(x[, c(2, 3, 6, 7)], 1, function(x, ages) {
        alphaI1D <- x[1]
        alphaI1H <- x[2]
        etaI1D <- x[3]
        etaI1H <- x[4]
        pI1D <- exp(alphaI1D + etaI1D * ages)
        pI1H <- exp(alphaI1H + etaI1H * ages)
        p <- pI1D + pI1H
        all(p >= 0 & p <= 1)
    }, ages = ages)
    multiProbs & singleProbs
}

## generate LHS design (200 + 50 validation)
ndesign <- 225
design <- randomLHS(ndesign * 50, nrow(parRanges))
colnames(design) <- parRanges$parameter
design <- as_tibble(design)

## convert to input space
inputs <- convertDesignToInput(design, parRanges, "zero_one")

## check probabilities are valid
inds <- select(inputs, starts_with("alpha"), starts_with("eta")) %>%
    as.matrix %>%
    pathwaysLimitFn(ages = c(2.5, 11, 23.5, 34.5, 44.5, 55.5, 65.5, 75.5))
inputs <- filter(inputs, inds)
stopifnot(nrow(inputs) >= ndesign)
inputs <- slice(inputs, 1:ndesign)

## generate LHS design (200 + 50 validation)
nval <- 75
design <- randomLHS(nval * 50, nrow(parRanges))
colnames(design) <- parRanges$parameter
design <- as_tibble(design)

## convert to input space
design <- convertDesignToInput(design, parRanges, "zero_one")

## check probabilities are valid
inds <- select(design, starts_with("alpha"), starts_with("eta")) %>%
    as.matrix %>%
    pathwaysLimitFn(ages = c(2.5, 11, 23.5, 34.5, 44.5, 55.5, 65.5, 75.5))
design <- filter(design, inds)
stopifnot(nrow(design) >= nval)
design <- slice(design, 1:nval)

## bind validation and training set together
inputs <- rbind(inputs, design)

## generate space-filling design for other parameters

## load FMM objects
hospStays <- readRDS("inputs/hospStays.rds")
hospThresh <- readRDS("inputs/hospThresh.rds")

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

## bind to design
inputs <- cbind(inputs, hospStaysInput)

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
disease <- convertInputToDisease(slice(inputs, 1:ndesign), C, N, S0, ages)
disease_val <- convertInputToDisease(slice(inputs, -c(1:ndesign)), C, N, S0, ages)
ndesign <- ndesign - 25
nval <- nval - 25
stopifnot(nrow(disease) >= ndesign & nrow(disease_val) >= nval)
disease <- slice(disease, 1:ndesign)
disease_val <- slice(disease_val, 1:nval)
disease <- rbind(disease, disease_val)

## match to inputs
inputs <- semi_join(inputs, disease, by = "output")

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

