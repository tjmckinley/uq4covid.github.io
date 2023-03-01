## load libraries
library(mclust)
library(tidyverse)

## check if being run in batch mode
args <- commandArgs(TRUE)
if(length(args) != 0) {
    ## extract command line arguments
    args <- commandArgs(TRUE)
    if(length(args) > 0) {
        stopifnot(length(args) == 3)
        wave <- args[1]
	prevwave <- args[2]
        file <- args[3]
    } else {
        stop("No arguments")
    }
} else {
    wave <- "2"
    prevwave <- "1"
    file <- "inputsWave2.csv"
}
print(wave)
print(prevwave)
print(file)

## copy fixed parameters over
stopifnot(file.copy(paste0("wave", prevwave, "/fixedInputs.txt"), paste0("wave", wave, "/fixedInputs.txt")))

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

##################################################################
####                  read in design                          ####
##################################################################

inputs <- read_csv(paste0("wave", wave, "/", file))

## check inputs are in correct ranges
ndesign <- nrow(inputs)
    
#### AT THIS POINT CHECK THAT inputs ARE IN THE CORRECT RANGES
#### AS GIVEN IN parRanges e.g.
for(par in parRanges$parameter) {
    inputs <- filter(inputs, between(!!sym(par), parRanges$lower[parRanges$parameter == par], parRanges$upper[parRanges$parameter == par]))
}
                  
if(nrow(inputs) != ndesign) stop("Design fails range checks")

## load mixture models and density threshold
hospStays <- readRDS("inputs/hospStays.rds")
hospThresh <- readRDS("inputs/hospThresh.rds")

hosp <- select(inputs, alphaTH, etaTH)

## check ranges
hosp <- hosp[hosp$etaTH > 0, ]

## check against prior density restrictions
hosp <- hosp[dens(as.matrix(hosp), hospStays$modelName, hospStays$parameters, logarithm = TRUE) > hospThresh, ]

if(nrow(hosp) != ndesign) stop("Design fails hosp checks")

## load ages
ages <- c(2.5, 11, 23.5, 34.5, 44.5, 55.5, 65.5, 75.5)

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

## check against prior density restrictions
pathwaysMod <- readRDS("inputs/pathways.rds")
pathThresh <- readRDS("inputs/pathThresh.rds")

pathways <- select(inputs, alphaEP, alphaI1D, alphaHD, alphaI1H, eta, etaI_scale, etaH_scale)

## check ranges
pathways <- pathways[pathways$alphaEP > -20 & pathways$alphaEP < 0, ]
pathways <- pathways[pathways$alphaI1D > -20 & pathways$alphaI1D < 0, ]
pathways <- pathways[pathways$alphaHD > -20 & pathways$alphaHD < 0, ]
pathways <- pathways[pathways$alphaI1H > -20 & pathways$alphaI1H < 0, ]
pathways <- pathways[pathways$eta > 0 & pathways$eta < 1, ]
pathways <- pathways[pathways$etaI_scale > 0.5 & pathways$etaI_scale < 2, ]
pathways <- pathways[pathways$etaH_scale > 0.5 & pathways$etaH_scale < 2, ]

## check probabilities valid
pathways <- pathways[pathwaysLimitFn(pathways, ages), ]

## check against prior density region
pathways <- pathways[dens(as.matrix(pathways)[, -c(6, 7)], pathwaysMod$modelName, pathwaysMod$parameters, logarithm = TRUE) > pathThresh, ]

if(nrow(pathways) != ndesign) stop("Design fails pathways checks")

#######################################################
####   NOW YOU HAVE DESIGN THAT PASSES ALL TESTS   ####
####   PROCEED TO BIND TOGETHER AND CONVERT        ####
#######################################################

## add unique hash identifier
## (at the moment don't use "a0" type ensembleID, because MetaWards
## parses to dates)
inputs$output <- ensembleIDGen(ensembleID = paste0("Ens", wave), nrow(inputs))

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
stopifnot(nrow(disease) == ndesign)

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

