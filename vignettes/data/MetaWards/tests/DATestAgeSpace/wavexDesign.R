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
        prevwave <- args[2]
    } else {
        stop("No arguments")
    }
} else {
    wave <- "2"
    prevwave <- "1"
}

## create directory to save samples
if(dir.exists(paste0("wave", wave))) {
    stop("Can't overwrite existing directory")
}
dir.create(paste0("wave", wave))

## copy fixed parameters over
file.copy(paste0("wave", prevwave, "/fixedInputs.txt"), paste0("wave", wave, "/fixedInputs.txt"))

## source dataTools
source("inputs/dataTools.R")

## set up parameter ranges for uniform ranges
parRanges <- data.frame(
    parameter = c("R0", "TE", "TP", "TI1", "TI2", "nuA", "beta_scale", "p_move", "MD_scale", "MD_time"),
    lower = c(2, 0.1, 1.2, 2.8, 0.0001, 0, 0, 0, 0.1, 0),
    upper = c(4.5, 2, 3, 4.5, 0.5, 1, 1, 1, 1, 50),
    stringsAsFactors = FALSE
) 

## read in contact matrix to use for NGM
C <- as.matrix(read.csv("inputs/POLYMOD_matrix.csv", header = FALSE))

##################################################################
####   the next section just generates samples randomly       ####
####   to illustrate input checks that need to be undertaken  ####
####   before building your next design - you can put         ####
####   whatever sampling mechanism you like in here           ####
##################################################################

## generate design (200 + 50 validation)
ndesign <- 250
valid <- 0
while(valid == 0) {
    ## generate samples
    design <- matrix(runif(ndesign * nrow(parRanges), 0, 1), nrow = ndesign)
    colnames(design) <- parRanges$parameter
    design <- as_tibble(design)

    ## convert to input space from (0, 1)
    inputs <- convertDesignToInput(design, parRanges, "zero_one")
    
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
    inputs <- inputs[inputs$MD_scale >= parRanges$lower[parRanges$parameter == "MD_scale"] & 
                      inputs$MD_scale <= parRanges$upper[parRanges$parameter == "MD_scale"], ]
    inputs <- inputs[inputs$MD_time >= parRanges$lower[parRanges$parameter == "MD_time"] & 
                      inputs$MD_time <= parRanges$upper[parRanges$parameter == "MD_time"], ]
    if(exists("temp_inputs")) {
        temp_inputs <- rbind(temp_inputs, inputs)
    } else {
        temp_inputs <- inputs
    }
    if(nrow(temp_inputs) >= ndesign) {
        valid <- 1
        inputs <- temp_inputs[1:ndesign, ]
        rm(temp_inputs)
    }
}

## generate design for other parameters: alphaTH, etaTH

## load mixture models and density threshold
hospStays <- readRDS("inputs/hospStays.rds")
hospThresh <- readRDS("inputs/hospThresh.rds")

valid <- 0
while(valid == 0) {
    hosp <- matrix(c(rnorm(ndesign, 0, 0.1), abs(rnorm(ndesign, 0, 0.1))), nrow = ndesign)
    colnames(hosp) <- c("alphaTH", "etaTH")
    hosp <- as_tibble(hosp)
    
    ## check ranges
    hosp <- hosp[hosp$etaTH > 0, ]
    
    ## check against prior density restrictions
    hosp <- hosp[dens(as.matrix(hosp), hospStays$modelName, hospStays$parameters, logarithm = TRUE) > hospThresh, ]
    if(exists("temp_inputs")) {
        temp_inputs <- rbind(temp_inputs, hosp)
    } else {
        temp_inputs <- hosp
    }
    if(nrow(temp_inputs) >= ndesign) {
        valid <- 1
        hosp <- temp_inputs[1:ndesign, ]
        rm(temp_inputs)
    }
}

## generate design for other parameters: alphaHD, alphaI1D, alphaI1H, alphaEP, eta

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

valid <- 0
while(valid == 0) {
    pathways <- matrix(c(runif(ndesign * 4, -20, 0), abs(rnorm(ndesign, 0, 0.1)), runif(ndesign, 0.5, 2), runif(ndesign, 0.5, 2)), nrow = ndesign)
    colnames(pathways) <- c("alphaEP", "alphaI1D", "alphaHD", "alphaI1H", "eta", "etaI_scale", "etaH_scale")
    pathways <- as_tibble(pathways)
    
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
    
    if(exists("temp_inputs")) {
        temp_inputs <- rbind(temp_inputs, pathways)
    } else {
        temp_inputs <- pathways
    }
    if(nrow(temp_inputs) >= ndesign) {
        valid <- 1
        pathways <- temp_inputs[1:ndesign, ]
        rm(temp_inputs)
    }
}

#######################################################
####   NOW YOU HAVE DESIGN THAT PASSES ALL TESTS   ####
####   PROCEED TO BIND TOGETHER AND CONVERT        ####
#######################################################

## bind to design
inputs <- cbind(inputs, hosp, pathways)

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

