## load libraries
library(tidyverse)

## check if being run in batch mode
args <- commandArgs(TRUE)
if(length(args) != 0) {
    ## extract command line arguments
    args <- commandArgs(TRUE)
    if(length(args) > 0) {
        stopifnot(length(args) == 2)
        wave <- as.numeric(args[1])
        hash <- as.numeric(args[2])
    } else {
        stop("No arguments")
    }
} else {
    ## set wave number and hash
    wave <- 1
    hash <- 1
}

## read in LADs from file
## (lad = NA gives national plots, else give vector of lads)
if(file.exists("JASMINcode/lads.txt")) {
    lads <- as.numeric(readLines("JASMINcode/lads.txt"))
} else {
    lads <- NA
}

## set writeExt (eventually pass in as argument)
writeExt <- FALSE

## extract file names
folder <- paste0("wave", wave)

## lookup table for classes
lookup <- data.frame(var = c("S", "E", "A", "RA", "P", "I1", "DI", "I2", "RI", "H", "RH", "DH", "DIobs", "DHobs")) %>%
    mutate(class = 0:(n() - 1))

if(!writeExt) {
    file <- paste0(folder, "/runs_md_", hash, ".rds")
    runs <- readRDS(file)$particles[[1]]
    
    ## aggregate to national level
    runs1 <- map(1:length(runs), function(j, y) {
        y <- y[[j]]
        map(1:length(y), function(i, x) {
            x <- apply(x[[i]], c(1, 2), sum)
            as.vector(x) %>%
            cbind(rep(1:dim(x)[2], each = dim(x)[1])) %>%
            cbind(rep(1:dim(x)[1], times = dim(x)[2])) %>%
            cbind(rep(i, nrow(.)))
        }, x = y) %>%
        {do.call("rbind", .)} %>%
        cbind(rep(j, nrow(.)))
    }, y = runs) %>%
    {do.call("rbind", .)}
    colnames(runs1) <- c("n", "age", "class", "particle", "time")
    runs1 <- as_tibble(runs1) %>%
        mutate(class = class - 1) %>%
        inner_join(lookup, by = "class") %>%
        select(!class) %>%
        mutate(lad = "agg") %>%
        select(n, lad, age, particle, time, var)
        
    ## extract subset of LADs also if required
    if(!is.na(lads[1])) {
        ## extract subset of LADs
        runs <- map(1:length(runs), function(j, y, lads) {
            y <- y[[j]]
            map(1:length(y), function(i, x, lads) {
                x <- x[[i]][, , lads]
                as.vector(x) %>%
                cbind(rep(lads, each = prod(dim(x)[-3]))) %>%
                cbind(rep(rep(1:dim(x)[2], each = dim(x)[1]), times = dim(x)[3])) %>%
                cbind(rep(1:dim(x)[1], times = prod(dim(x)[-1]))) %>%
                cbind(rep(i, nrow(.)))
            }, x = y, lads = lads) %>%
            {do.call("rbind", .)} %>%
            cbind(rep(j, nrow(.)))
        }, y = runs, lads = lads) %>%
        {do.call("rbind", .)}
        colnames(runs) <- c("n", "lad", "age", "class", "particle", "time")
        runs <- as_tibble(runs) %>%
            mutate(class = class - 1) %>%
            inner_join(lookup, by = "class") %>%
            select(!class)
        runs1 <- rbind(runs1, runs)
    }
    saveRDS(runs1, paste0(folder, "/plotSum_", hash, ".rds"))
} else {
    stop("Not yet implemented when writeExt == FALSE")
}

print("Finished")


