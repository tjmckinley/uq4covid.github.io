## load libraries
library(mclust)
library(tidyverse)

## check if being run in batch mode
args <- commandArgs(TRUE)
if(length(args) != 0) {
    ## extract command line arguments
    args <- commandArgs(TRUE)
    if(length(args) > 0) {
        stopifnot(length(args) == 2)
        wave <- args[1]
        file <- args[2]
    } else {
        stop("No arguments")
    }
} else {
    wave <- "2"
    file <- "inputsWave2.csv"
}

## source dataTools
source("inputs/dataTools.R")

## set up parameter ranges for uniform ranges
parRanges <- data.frame(
    parameter = c("R0", "TE", "TP", "TI1", "TI2", "nuA", "beta_scale", "p_move", "a_dis"),
    lower = c(2, 0.1, 1.2, 2.8, 0.0001, 0, 0, 0, 0.00014),
    upper = c(4.5, 2, 3, 4.5, 0.5, 1, 1, 1, 0.05),
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
inputs <- inputs[inputs$nuA >= parRanges$lower[parRanges$parameter == "nuA"] & 
                  inputs$nuA <= parRanges$upper[parRanges$parameter == "nuA"], ]
inputs <- inputs[inputs$R0 >= parRanges$lower[parRanges$parameter == "R0"] & 
                  inputs$R0 <= parRanges$upper[parRanges$parameter == "R0"], ]
inputs <- inputs[inputs$TE >= parRanges$lower[parRanges$parameter == "TE"] & 
                  inputs$TE <= parRanges$upper[parRanges$parameter == "TE"], ]
inputs <- inputs[inputs$TI1 >= parRanges$lower[parRanges$parameter == "TI1"] & 
                  inputs$TI1 <= parRanges$upper[parRanges$parameter == "TI1"], ]
inputs <- inputs[inputs$TI2 >= parRanges$lower[parRanges$parameter == "TI2"] & 
                  inputs$TI2 <= parRanges$upper[parRanges$parameter == "TI2"], ]
inputs <- inputs[inputs$TP >= parRanges$lower[parRanges$parameter == "TP"] & 
                  inputs$TP <= parRanges$upper[parRanges$parameter == "TP"], ]
inputs <- inputs[inputs$beta_scale >= parRanges$lower[parRanges$parameter == "beta_scale"] & 
                  inputs$beta_scale <= parRanges$upper[parRanges$parameter == "beta_scale"], ]
inputs <- inputs[inputs$p_move >= parRanges$lower[parRanges$parameter == "p_move"] & 
                  inputs$p_move <= parRanges$upper[parRanges$parameter == "p_move"], ]
inputs <- inputs[inputs$a_dis >= parRanges$lower[parRanges$parameter == "a_dis"] & 
                  inputs$a_dis <= parRanges$upper[parRanges$parameter == "a_dis"], ]
                  
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

## this function checks validity of inputs as probabilities
pathwaysLimitFn <- function(x, ages) {
    ## check all parameters give valid probabilities
    ## in (0, 1)
    singleProbs <- apply(x, 1, function(x, ages) {
        eta <- x[5]
        alphas <- x[-5]
        y <- sapply(alphas, function(a, eta, ages) {
            y <- exp(a + eta * ages)
            all(y >= 0 & y <= 1)
        }, eta = eta, ages = ages)
        all(y)
    }, ages = ages)
    ## check multinomial probabilities sum to one
    multiProbs <- apply(x[, -c(1, 3)], 1, function(x, ages) {
        alphaI1D <- x[1]
        alphaI1H <- x[2]
        eta <- x[3]
        pI1D <- exp(alphaI1D + eta * ages)
        pI1H <- exp(alphaI1H + eta * ages)
        p <- pI1D + pI1H
        all(p >= 0 & p <= 1)
    }, ages = ages)
    multiProbs & singleProbs
}

## check against prior density restrictions
pathwaysMod <- readRDS("inputs/pathways.rds")
pathThresh <- readRDS("inputs/pathThresh.rds")

pathways <- select(inputs, alphaEP, alphaI1D, alphaHD, alphaI1H, eta)

## check ranges
pathways <- pathways[pathways$alphaEP > -20 & pathways$alphaEP < 0, ]
pathways <- pathways[pathways$alphaI1D > -20 & pathways$alphaI1D < 0, ]
pathways <- pathways[pathways$alphaHD > -20 & pathways$alphaHD < 0, ]
pathways <- pathways[pathways$alphaI1H > -20 & pathways$alphaI1H < 0, ]
pathways <- pathways[pathways$eta > 0 & pathways$eta < 1, ]

## check probabilities valid
pathways <- pathways[pathwaysLimitFn(pathways, ages), ]

## check against prior density region
pathways <- pathways[dens(as.matrix(pathways), pathwaysMod$modelName, pathwaysMod$parameters, logarithm = TRUE) > pathThresh, ]

if(nrow(pathways) != ndesign) stop("Design fails pathways checks")

#######################################################
####   NOW YOU HAVE DESIGN THAT PASSES ALL TESTS   ####
####   PROCEED TO BIND TOGETHER AND CONVERT        ####
#######################################################

## add unique hash identifier
## (at the moment don't use "a0" type ensembleID, because MetaWards
## parses to dates)
inputs$output <- ensembleIDGen(ensembleID = paste0("Ens", wave), nrow(inputs))
inputs$repeats <- 1

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
p <- select(inputs, -output, -repeats) %>%
    ggpairs(upper = "blank")
ggsave(paste0("wave", wave, "/design.pdf"), p, width = 10, height = 10)

