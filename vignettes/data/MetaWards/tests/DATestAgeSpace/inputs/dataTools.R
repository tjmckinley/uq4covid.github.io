## @knitr convertDesignToInput
## function to convert from design to input space
convertDesignToInput <- function(design, parRanges, scale = c("zero_one", "negone_one")) {
  
    require(tidyverse)
  
    ## convert from design space to input space
    input <- mutate(design, ind = 1:n()) %>%
        gather(parameter, value, -ind) %>%
        left_join(parRanges, by = "parameter")
    if(scale[1] == "negone_one") {
        input <- mutate(input, value = value / 2 + 1)
    }
    input <- mutate(input, value = value * (upper - lower) + lower) %>%
        dplyr::select(ind, parameter, value) %>%
        spread(parameter, value) %>%
        arrange(ind) %>%
        dplyr::select(-ind)
    
    ## return inputs
    input
}

## @knitr convertInputToDisease
## function to convert from input to disease space
convertInputToDisease <- function(input, C, N, S0, ages) {
   
    require(tidyverse)
    require(magrittr)
  
    stopifnot(all(c("R0", "nuA", "TE", "TP", "TI1", "TI2", "alphaEP", 
        "alphaI1H", "alphaI1D", "alphaHD", "eta",  
        "alphaTH", "etaTH", "output", "beta_scale", "p_move", "MD_time") %in% colnames(input)))
    
    ## check unique ID
    stopifnot(length(unique(input$output)) == length(input$output))
    
    ## scaling for asymptomatics and lockdown
    disease <- select(input, nuA, beta_scale, p_move, MD_time, output)
  
    ## progressions out of the E class
    for(j in 1:length(ages)) {
        disease <- mutate(disease, ".pE_{j}" := 1 - exp(-1 / input$TE))
    }
    disease <- mutate(disease, 
        temp = map2(input$alphaEP, input$eta, function(alpha, eta, ages) {
            exp(alpha + eta * ages) %>%
                matrix(nrow = 1) %>%
                set_colnames(paste0(".pEP_", 1:length(ages))) %>%
                as_tibble()
        }, ages = ages)) %>%
        unnest(cols = temp)
    
    ## progressions out of the A class
    for(j in 1:length(ages)) {
        disease <- mutate(disease, ".pA_{j}" := 1 - exp(-1 / (input$TP + input$TI1 + input$TI2)))
    }
    
    ## progressions out of the P class
    for(j in 1:length(ages)) {
      disease <- mutate(disease, ".pP_{j}" := 1 - exp(-1 / input$TP))
    }
    
    ## progressions out of the I1 class
    for(j in 1:length(ages)) {
        disease <- mutate(disease, ".pI1_{j}" := 1 - exp(-1 / input$TI1))
    }
    disease <- mutate(disease, 
        temp = map2(input$alphaI1H, input$eta, function(alpha, eta, ages) {
            exp(alpha + eta * ages) %>%
                matrix(nrow = 1) %>%
                set_colnames(paste0(".pI1H_", 1:length(ages))) %>%
                as_tibble()
        }, ages = ages)) %>%
        unnest(cols = temp)
    disease <- mutate(disease, 
        temp = map2(input$alphaI1D, input$eta, function(alpha, eta, ages) {
            exp(alpha + eta * ages) %>%
                matrix(nrow = 1) %>%
                set_colnames(paste0(".pI1D_", 1:length(ages))) %>%
                as_tibble()
        }, ages = ages)) %>%
        unnest(cols = temp)
    ## check for multinomial validity
    temp <- select(disease, starts_with(".pI1")) %>%
        select(-starts_with(".pI1_")) %>%
        mutate(ind = 1:n()) %>%
        gather(par, prob, -ind) %>%
        separate(par, c("par", "age"), sep = "_") %>%
        group_by(age, ind) %>% 
        summarise(prob = sum(prob)) %>%
        pluck("prob")
    if(any(temp < 0) | any(temp > 1)) {
        stop("Some multinomial pI1* probs invalid")
    }
    
    ## progressions out of the I2 class
    for(j in 1:length(ages)) {
        disease <- mutate(disease, ".pI2_{j}" := 1 - exp(-1 / input$TI2))
    }
    
    ## progressions out of the H class
    disease <- mutate(disease, 
        temp = map2(input$alphaTH, input$etaTH, function(alpha, eta, ages) {
          exp(alpha + eta * ages) %>%
              {1 - exp(-1 / .)} %>%
              matrix(nrow = 1) %>%
              set_colnames(paste0(".pH_", 1:length(ages))) %>%
              as_tibble()
        }, ages = ages)) %>%
        unnest(cols = temp)
    disease <- mutate(disease, 
        temp = map2(input$alphaHD, input$eta, function(alpha, eta, ages) {
            exp(alpha + eta * ages) %>%
                matrix(nrow = 1) %>%
                set_colnames(paste0(".pHD_", 1:length(ages))) %>%
                as_tibble()
        }, ages = ages)) %>%
        unnest(cols = temp)
    
    ## remove any invalid inputs
    stopifnot(any(disease$nuA >= 0 | disease$nuA <= 1))
    stopifnot(
        select(disease, starts_with(".p")) %>%
        summarise_all(~{all(. >= 0 & . <= 1)}) %>%
        apply(1, all)
    )

    ## set up data for calculating transmission parameter
    temp <- select(disease, output, starts_with(".pEP"), starts_with(".pI1H"), starts_with(".pI1D")) %>%
        gather(prob, value, -output) %>%
        separate(prob, c("prob", "age"), sep = "_") %>%
        spread(age, value) %>%
        group_by(prob, output) %>%
        nest() %>%
        ungroup() %>%
        spread(prob, data) %>%
        mutate_at(vars(-output), ~purrr::map(., as.data.frame)) %>%
        mutate_at(vars(-output), ~purrr::map(., unlist))
    colnames(temp) <- gsub("\\.", "", colnames(temp))
    
    temp <- select(input, output, R0, nuA, TE, TP, TI1, TI2) %>%
        inner_join(temp, by = "output")

    temp <- mutate(temp, nu = pmap_dbl(as.list(temp)[-1], function(R0, nuA, TE, TP, TI1, TI2, pEP, pI1H, pI1D, C, N, S0) {
            ## transformations
            gammaE <- rep(1 / TE, length(N))
            gammaP <- rep(1 / TP, length(N))
            gammaI1 <- rep(1 / TI1, length(N))
            gammaI2 <- rep(1 / TI2, length(N))
            gammaA <- rep(1 / (TP + TI1 + TI2), length(N))
            
            ## calculate nu from NGM
            NGM(R0 = R0, nu = NA, C, S0, N, nuA, gammaE, pEP, gammaA, gammaP, gammaI1, 
                pI1H, pI1D, gammaI2)$nu
        }, C = C, N = N, S0 = S0))
    disease <- inner_join(disease, select(temp, nu, output), by = "output")
    
    ## checks on nu
    stopifnot(all(disease$nu > 0 & disease$nu < 1))
    
    ## finalise data set
    disease <- mutate(disease, nuA = nuA * nu) %>%
        inner_join(select(input, output), by = "output")
        
    ## reorder
    disease <- select(disease, nu, nuA, !c(nu, nuA, beta_scale, p_move, MD_time, output), beta_scale, p_move, MD_time, output)
    
    print(paste0(nrow(input) - nrow(disease), " invalid inputs removed"))
    print(paste0(nrow(disease), " samples remaining"))
    
    ## return disease file
    disease
}

## @knitr maximin
## function to generate maximin samples given an arbitrary FMM
FMMmaximin <- function(model, nsamp, limits, limitFn = NULL, nseed = 10000, ...) {
    
    ## check inputs and dependencies
    require(mclust)
    require(fields)
  
    stopifnot(!missing(model) & !missing(nsamp) & !missing(limits))
    stopifnot(class(model)[1] == "densityMclust")
    stopifnot(is.numeric(nsamp) & round(nsamp) == nsamp & nsamp > 1)
    stopifnot(is.numeric(nseed) & round(nseed) == nseed & nseed > 1 & nsamp < nseed)
    stopifnot(is.matrix(limits))
    stopifnot(ncol(limits) == 2 & nrow(limits) == nrow(model$parameters$mean))
    stopifnot(all(limits[, 1] < limits[, 2]))
    if(!is.null(limitFn)) {
        stopifnot(class(limitFn) == "function")
        stopifnot(formalArgs(limitFn)[1] == "x")
        cat("When using 'limitFn' ensure that function takes matrix as first argument 'x'
            and returns vector of logicals for inclusion of length 'nrow(x)'\n")
    }
    
    ## produce large number of samples from model and ensure they are
    ## consistent with limits
    sims <- matrix(NA, 1, 1)
    while(nrow(sims) < nseed) {
        
        ## sample from FMM
        sims1 <- sim(model$modelName, model$parameters, nseed)[, -1]
    
        ## check against limits
        for(i in 1:nrow(limits)) {
            sims1 <- sims1[sims1[, i] >= limits[i, 1], ]
            sims1 <- sims1[sims1[, i] <= limits[i, 2], ]
        }
        ## check against limitFn
        if(!is.null(limitFn)) {
            sims1 <- sims1[limitFn(sims1, ...), ]
        }
        if(nrow(sims1) > 0) {
            if(nrow(sims) == 1) {
                sims <- sims1
            } else {
                sims <- rbind(sims, sims1)
            }
        }
    }
    sims <- sims[1:nseed, ]
    
    ## rescale
    simsnorm <- apply(sims, 2, function(x) (x - min(x)) / diff(range(x)))
    
    ## sample initial choice at random
    simind <- sample.int(nrow(sims), 1)
    simrem <- c(1:nrow(sims))[-simind]
    
    ## choose other points using maximin
    while(length(simind) < nsamp) {
        dists <- rdist(simsnorm[simind, , drop = FALSE], simsnorm[-simind, , drop = FALSE])
        dists <- apply(dists, 2, min)
        dists <- which(dists == max(dists))
        simind <- c(simind, simrem[dists[1]])
        simrem <- simrem[-dists[1]]
    }
    
    ## return design
    sims[simind, ]
}

## @knitr ensembleIDGen
## function for generating ensemble hash IDs (from Danny)
ensembleIDGen <- function(ensembleID = "a0", ensembleSize) {
    HexString <- c(as.character(0:9), letters)
    counter1 <- 1
    counter2 <- 1
    counter3 <- 1
    TotalGenerated <- 0
    tIDs <- rep(ensembleID, ensembleSize)
    while(counter1 <= 36 & TotalGenerated < ensembleSize){
      while(counter2 <= 36 & TotalGenerated < ensembleSize){
        while(counter3 <= 36 & TotalGenerated < ensembleSize){
          TotalGenerated <- TotalGenerated + 1
          tIDs[TotalGenerated] <- paste0(tIDs[TotalGenerated], HexString[counter1], HexString[counter2], HexString[counter3])
          counter3 <- counter3 + 1
        }
        counter2 <- counter2 + 1
        counter3 <- 1
      }
      counter1 <- counter1 + 1
      counter2 <- 1
      counter3 <- 1
    }
    stopifnot(!any(duplicated(tIDs)))
    tIDs
}

## @knitr reconstruct
## write Rcpp function to reconstruct counts from incidence
## first element of Einc etc. must be initial states and the following
## must be the incidence at a given time
library(Rcpp)
cppFunction('IntegerMatrix reconstruct(IntegerVector Einc, IntegerVector Pinc, 
    IntegerVector I1inc, IntegerVector I2inc, IntegerVector RIinc, IntegerVector DIinc, 
    IntegerVector Ainc, IntegerVector RAinc, 
    IntegerVector Hinc, IntegerVector RHinc, IntegerVector DHinc) {
    
    // extract sizes
    int n = Einc.size();
    
    // set up output matrix
    IntegerMatrix output(n, 17);
    
    // set up initial counts
    output(0, 0) = 0;
    output(0, 1) = Einc[0];
    output(0, 2) = 0;
    output(0, 3) = Pinc[0];
    output(0, 4) = 0;
    output(0, 5) = I1inc[0];
    output(0, 6) = 0;
    output(0, 7) = I2inc[0];
    output(0, 8) = RIinc[0];
    output(0, 9) = DIinc[0];
    output(0, 10) = 0;
    output(0, 11) = Ainc[0];
    output(0, 12) = RAinc[0];
    output(0, 13) = 0;
    output(0, 14) = Hinc[0];
    output(0, 15) = RHinc[0];
    output(0, 16) = DHinc[0];
    
    int E = Einc[0]; 
    int P = Pinc[0]; 
    int I1 = I1inc[0]; 
    int I2 = I2inc[0];
    int RI = RIinc[0];
    int DI = DIinc[0];
    int A = Ainc[0];
    int RA = RAinc[0];
    int H = Hinc[0];
    int RH = RHinc[0];
    int DH = DHinc[0];
    
    // reconstruct counts
    for(int i = 1; i < n; i++) {
    
        E += Einc[i] - Pinc[i] - Ainc[i];
        output(i, 0) = Einc[i];
        output(i, 1) = E;
    
        P += Pinc[i] - I1inc[i];
        output(i, 2) = Pinc[i];
        output(i, 3) = P;
        
        I1 += I1inc[i] - I2inc[i] - Hinc[i] - DIinc[i];
        output(i, 4) = I1inc[i];
        output(i, 5) = I1;
        
        I2 += I2inc[i] - RIinc[i];
        output(i, 6) = I2inc[i];
        output(i, 7) = I2;
        
        RI += RIinc[i];
        output(i, 8) = RI;
        
        DI += DIinc[i];
        output(i, 9) = DI;
        
        A += Ainc[i] - RAinc[i];
        output(i, 10) = Ainc[i];
        output(i, 11) = A;
        
        RA += RAinc[i];
        output(i, 12) = RA;
        
        H += Hinc[i] - RHinc[i] - DHinc[i];
        output(i, 13) = Hinc[i];
        output(i, 14) = H;
        
        RH += RHinc[i];
        output(i, 15) = RH;
        
        DH += DHinc[i];
        output(i, 16) = DH;
    }
    
    // return counts
    return(output);
}')

## @knitr NGM
## NGM (order: E, A, P, I1, I2)
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
