## load mixture models and density threshold
hospStays <- readRDS("../inputs/hospStays.rds")
hospThresh <- readRDS("../inputs/hospThresh.rds")
pathwaysMod <- readRDS("../inputs/pathways.rds")
pathwaysThresh <- readRDS("../inputs/pathThresh.rds")

## solution to round numbers preserving sum
## adapted from:
## https://stackoverflow.com/questions/32544646/round-vector-of-numerics-to-integer-while-preserving-their-sum
smart_round <- function(x) {
    y <- floor(x)
    indices <- tail(order(x - y), round(sum(x)) - sum(y))
    y[indices] <- y[indices] + 1
    y
}

## load ages
ages <- c(2.5, 11, 23.5, 34.5, 44.5, 55.5, 65.5, 75.5)

## read in objects for calculating NGM
C <- as.matrix(read.csv("../inputs/POLYMOD_matrix.csv", header = FALSE))
N <- smart_round(read_csv("../inputs/age_seeds.csv", col_names = FALSE)$X2 * 56082077)
S0 <- N - smart_round(read_csv("../inputs/age_seeds.csv", col_names = FALSE)$X2 * 100)

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

## Matern kernel function for Voronoi region
matern <- function(x1, x2, nu, l) {

    ## calculate distance between x1 and x2
    d <- rdist(x1, x2)
    d <- sqrt(2 * nu) * d / l
    
    ## calculate Matern function
    (2^(1 - nu)) * (d^nu) * besselK(d, nu) / gamma(nu)
}

## set up custom implausibility
implausibility <- function(x, z, cutoff = NULL, em, hospStays, hospThresh, pathwaysMod, pathwaysThresh, ages, pathwaysLimitFn, epsilon, range_list, actives, NGM, N, S0, C, target, pcutoff, voronoi) {
    if (!is.numeric(z) && !is.null(z$val)) {
    
        if (!is.null(cutoff)) {
            if(cutoff != 1) {
                print("'cutoff' must be = 1 for pseudo-implausibility to work")
                cutoff <- 1
            }
        }
        
        ## extract active variables
        y <- as.matrix(x)
        y <- y[, actives, drop = FALSE]
        
        ## rescale y to (0, 1)
        for(i in 1:ncol(y)) {
            temp_ranges <- range_list[[colnames(y)[i]]]
            y[, i] <- (y[, i] - temp_ranges[1]) / diff(temp_ranges)
        }
        stopifnot(all(apply(y, 2, function(x) all(x > 0 & x < 1))))
        
        ## set up implausibilities
        imp <- rep(NA, nrow(y))

        ## check against prior criteria
        hosp <- select(x, alphaTH, etaTH)
        
        ## check against prior density restrictions
        imp[dens(as.matrix(hosp), hospStays$modelName, hospStays$parameters, logarithm = TRUE) < hospThresh] <- 2

        ## check probabilities valid
        pathways <- select(x, alphaEP, alphaI1D, alphaHD, alphaI1H, eta, etaI_scale, etaH_scale)
        imp[!pathwaysLimitFn(pathways, ages)] <- 2
        
        ## check against prior density restrictions
        imp[dens(as.matrix(pathways)[, 1:5, drop = FALSE], pathwaysMod$modelName, pathwaysMod$parameters, logarithm = TRUE) < pathwaysThresh] <- 2
        imp[pathways[, 6] < 0.5 | pathways[, 6] > 2] <- 2
        imp[pathways[, 7] < 0.5 | pathways[, 7] > 2] <- 2
        
        if(sum(is.na(imp)) > 0) {
            ## check nu valid
            nu <- filter(x, is.na(imp)) %>%
                mutate(gammaE = 1 / TE) %>%
                mutate(gammaP = 1 / TP) %>%
                mutate(gammaI1 = 1 / TI1) %>%
                mutate(gammaI2 = 1 / TI2) %>%
                mutate(gammaA = 1 / (TP + TI1 + TI2)) %>%
                select(R0, nuA, gammaE, gammaP, gammaI1, gammaI2, gammaA, alphaEP, alphaI1D, alphaI1H, eta, etaI_scale) %>%
                as.matrix() %>%
                apply(1, function(x, ages, N, S0, C) {
                    
                    ## create pathway probabilities
                    pEP <- exp(x[8] + x[11] * ages)
                    pI1D <- exp(x[9] + x[11] * x[12] * ages)
                    pI1H <- exp(x[10] + x[11] * ages)
                    
                    ## expand rates
                    gammaE <- rep(x[3], length(ages))
                    gammaP <- rep(x[4], length(ages))
                    gammaI1 <- rep(x[5], length(ages))
                    gammaI2 <- rep(x[6], length(ages))
                    gammaA <- rep(x[7], length(ages))
                
                    ## calculate nu from NGM
                    NGM(R0 = x[1], nu = NA, C, S0, N, x[2], gammaE, pEP, gammaA, gammaP, gammaI1, 
                        pI1H, pI1D, gammaI2)$nu
                }, ages = ages, N = N, S0 = S0, C = C)
            imp[is.na(imp)][nu < 0 | nu > 1] <- 2
        
            ## check against Voronoi region
            if(!is.null(voronoi) & sum(is.na(imp)) > 0) {
            
                ## extract subset of candidate points
                y1 <- y[is.na(imp), , drop = FALSE]
                
                ## correlations against doubt points
                yXD <- voronoi$kernel(
                    y1,
                    em$data$X[voronoi$XD, , drop = FALSE], 
                    2.5,
                    em$specs$layer1$node1$lengthscales
                )
                ## correlations against normal points
                yXN <- voronoi$kernel(
                    y1,
                    em$data$X[-voronoi$XD, , drop = FALSE], 
                    2.5,
                    em$specs$layer1$node1$lengthscales
                )
                yXD <- apply(cbind(yXD, yXN), 1, function(x, nXD) {
                    x[is.na(x)] <- max(x, na.rm = TRUE) + 1
                    sort.list(x, decreasing = TRUE)[1] %in% 1:nXD
                }, nXD = length(voronoi$XD))
                
                ## set points in Voronoi region to be NROY
                imp[is.na(imp)][yXD] <- 0
            }
        
            ## get variance on current scale
            y <- y[is.na(imp), , drop = FALSE]
            imp_var <- predict(em, y)$results
            imp_pred <- imp_var$mean
            imp_var <- imp_var$var
            
            ## get estimate of probability
            imp1 <- pnorm(log(-log(epsilon) - target$ll$val), mean = imp_pred, sd = sqrt(imp_var))
            
            ## check against cutoff and return pseudo-implausibility
            imp1 <- ifelse((1 - imp1) <= pcutoff, 0, 2)
            
            imp[is.na(imp)] <- imp1
        }
    } else {
        stop("Not sure what to put here")
    }
    if (is.null(cutoff)) return(imp)
    return(imp <= cutoff)
}

## mean function
pred_func <- function(x, em, range_list, actives) {
    x <- x[, actives, drop = FALSE]
    for(i in 1:ncol(x)) {
        temp_ranges <- range_list[[colnames(x)[i]]]
        x[, i] <- (x[, i] - temp_ranges[1]) / diff(temp_ranges)
    }
    stopifnot(all(apply(x, 2, function(x) all(x > 0 & x < 1))))
    predict(em, as.matrix(x))$results$mean
}  

## variance function
pred_var_func <- function(x, em, range_list, actives) {
    x <- x[, actives, drop = FALSE]
    for(i in 1:ncol(x)) {
        temp_ranges <- range_list[[colnames(x)[i]]]
        x[, i] <- (x[, i] - temp_ranges[1]) / diff(temp_ranges)
    }
    stopifnot(all(apply(x, 2, function(x) all(x > 0 & x < 1))))
    predict(em, as.matrix(x))$results$var
}

NGM <- function(R0 = NA, nu = NA, C, S0, N, nuA, gammaE, pEP, gammaA, gammaP, gammaI1, pI1H, pI1D, gammaI2) {
  
    ## check that exactly one of R0 or nu is specified
    stopifnot(length(R0) == 1 & length(nu) == 1 & length(nuA) == 1)
    stopifnot(!is.na(R0) | !is.na(nu))
    stopifnot(is.na(R0) | is.na(nu))
    
    ## check inputs
    stopifnot(!missing(C) & !missing(N) & !missing(S0))
    if(missing(nuA)) nuA <- 0
    if(missing(gammaE)) gammaE <- rep(0, length(N))
    if(missing(pEP)) pEP <- rep(0, length(N))
    if(missing(gammaA)) gammaA <- rep(0, length(N))
    if(missing(gammaP)) gammaP <- rep(0, length(N))
    if(missing(gammaI1)) gammaI1 <- rep(0, length(N))
    if(missing(pI1H)) pI1H <- rep(0, length(N))
    if(missing(pI1D)) pI1D <- rep(0, length(N))
    if(missing(gammaI2)) gammaI2 <- rep(0, length(N))
    stopifnot(all(c(nuA, gammaE, pEP, gammaA, 
                    gammaP, gammaI1, pI1H, pI1D, gammaI2) >= 0))
    
    ## extract lengths
    nage <- length(N)
    stopifnot(all(dim(C) == nage))
    stopifnot(length(S0) == nage & all(S0 <= N) & any(S0 >= 1))
    stopifnot(all(map_lgl(list(gammaE, pEP, gammaA, gammaP, gammaI1, pI1H, pI1D, gammaI2), ~{length(.) == nage})))
    
    ## check for empty pathways
    stopifnot(all(gammaE > 0))
    pathA <- ifelse(pEP == 1, FALSE, TRUE)
    if(any(pathA)) stopifnot(gammaA[pathA] > 0)
    pathI1 <- ifelse(pEP == 0, FALSE, TRUE)
    if(any(pathI1)) stopifnot(gammaP[pathI1] > 0 & gammaI1[pathI1] > 0)
    pathI2 <- ifelse((1 - pI1H - pI1D) == 0 | !pathI1, FALSE, TRUE)
    if(any(pathI2)) stopifnot(gammaI2[pathI2] > 0)
    
    ## set up empty F matrix
    ## (order: E, A, P, I1, I2)
    F <- matrix(0, nage * 5, nage * 5)
    
    ## add A to E
    F[1:nage, (nage + 1):(2 * nage)] <- t(t(nuA * S0 * C) / ifelse(N == 0, 1, N))
    ## add P to E
    F[1:nage, (2 * nage + 1):(3 * nage)] <- t(t(S0 * C) / ifelse(N == 0, 1, N))
    ## add I1 to E
    F[1:nage, (3 * nage + 1):(4 * nage)] <- t(t(S0 * C) / ifelse(N == 0, 1, N))
    ## add I2 to E
    F[1:nage, (4 * nage + 1):(5 * nage)] <- t(t(S0 * C) / ifelse(N == 0, 1, N))
    
    ## set up empty V matrix
    V <- matrix(0, nage * 5, nage * 5)
    
    ## add E to E terms
    diag(V)[1:nage] <- gammaE
    
    ## add E to A terms
    V[(nage + 1):(2 * nage), 1:nage] <- diag(nage) * -(1 - pEP) * gammaE
    ## add A to A terms
    V[(nage + 1):(2 * nage), (nage + 1):(2 * nage)] <- diag(nage) * gammaA
    
    ## add E to P terms
    V[(2 * nage + 1):(3 * nage), 1:nage] <- diag(nage) * -pEP * gammaE
    ## add P to P terms
    V[(2 * nage + 1):(3 * nage), (2 * nage + 1):(3 * nage)] <- diag(nage) * gammaP
    
    ## add P to I1 terms
    V[(3 * nage + 1):(4 * nage), (2 * nage + 1):(3 * nage)] <- diag(nage) * -gammaP
    ## add I1 to I1 terms
    V[(3 * nage + 1):(4 * nage), (3 * nage + 1):(4 * nage)] <- diag(nage) * gammaI1
    
    ## add I1 to I2 terms
    V[(4 * nage + 1):(5 * nage), (3 * nage + 1):(4 * nage)] <- diag(nage) * -(1 - pI1H - pI1D) * gammaI1
    ## add I2 to I2 terms
    V[(4 * nage + 1):(5 * nage), (4 * nage + 1):(5 * nage)] <- diag(nage) * gammaI2
    
    ## remove pathways that don't exist
    rem <- NULL
    if(!all(pathA)) rem <- nage + which(!pathA)
    if(!all(pathI1)) {
        rem <- c(rem, reduce(map(c(2 * nage, 3 * nage, 4 * nage), function(x, pathI1) {
            x + which(!pathI1)
        }, pathI1 = pathI1), c))
    }
    if(!all(pathI2)) rem <- c(rem, 4 * nage + which(!pathI2))
    if(!is.null(rem)) {
        rem <- unique(rem)
        F <- F[-rem, -rem]
        V <- V[-rem, -rem]
    }
    
    if(!is.na(R0)) {
        ## calculate NGM
        K <- F %*% solve(V)
        
        ## extract eigenvalues
        eig <- eigen(K)$values
        if(is.complex(eig)) {
            ## extract real eigenvalues
            eig <- as.numeric(eig[Im(eig) == 0])
            stopifnot(length(eig) >= 1)
        }
        
        ## calculate and return nu
        nu <- R0 / max(eig)
        return(list(nu = nu, K = nu * K))
    } else {        
        ## calculate NGM
        K <- nu * F %*% solve(V)
        
        ## extract eigenvalues
        eig <- eigen(K)$values
        if(is.complex(eig)) {
            ## extract real eigenvalues
            eig <- as.numeric(eig[Im(eig) == 0])
            stopifnot(length(eig) >= 1)
        }
        
        ## calculate and return R0
        R0 <- max(eig)
        return(list(R0 = R0, K = K))
    }
}

## function to generate maximin samples given an arbitrary set of inputs
maximin <- function(candidate_points, nsamp, ranges) {

    ## check nsamp < nrow(candidate_points)
    stopifnot(nsamp < nrow(candidate_points))
    stopifnot(identical(sort(names(ranges)), sort(colnames(candidate_points))))
    orig_points <- candidate_points
    
    ## scale inputs to (0, 1)
    for(i in 1:ncol(candidate_points)) {
        temp_ranges <- ranges[[colnames(candidate_points)[i]]]
        candidate_points[, i] <- (candidate_points[, i] - temp_ranges[1]) / diff(temp_ranges)
    }
    stopifnot(all(apply(candidate_points, 2, function(x) all(x > 0 & x < 1))))
    
    ## sample initial choice at random
    cand_ind <- sample.int(nrow(candidate_points), 1)
    cand_rem <- c(1:nrow(candidate_points))[-cand_ind]
    
    ## choose other points using maximin
    while(length(cand_ind) < nsamp) {
        dists <- rdist(candidate_points[cand_ind, , drop = FALSE], candidate_points[-cand_ind, , drop = FALSE])
        dists <- apply(dists, 2, min)
        dists <- which(dists == max(dists))
        cand_ind <- c(cand_ind, cand_rem[dists[1]])
        cand_rem <- cand_rem[-dists[1]]
    }
    
    ## return indices
    cand_ind
}

