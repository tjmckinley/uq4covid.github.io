## log-sum-exp function
log_sum_exp <- function(x, mn = FALSE) {
    maxx <- max(x)
    y <- maxx + log(sum(exp(x - maxx)))
    if(mn) y <- y - log(length(x))
    y
}

## function to check counts (mainly useful for error checking)
## u: matrix of counts in order:
##       S, E, A, RA, P, I1, DI, I2, RI, H, RH, DH
## cu: matrix of cumulative counts
## N: population size
checkCounts <- function(u, cu, N) {
    
    ## check states match population size in each age-class
    stopifnot(all((colSums(u) - N) == 0))
    
    ## check states are positive
    stopifnot(all(u >= 0) & all(cu >= 0))
    
    ## check cumulative counts match up
    cu <- t(cu)
    stopifnot(all((cu[, 1] + cu[, 2] - N) == 0))
    stopifnot(all((cu[, 2] - rowSums(cu[, c(3, 5)])) >= 0))
    stopifnot(all((cu[, 3] - cu[, 4]) >= 0))
    stopifnot(all((cu[, 5] - cu[, 6]) >= 0))
    stopifnot(all((cu[, 6] - rowSums(cu[, c(7, 8, 10)])) >= 0))
    stopifnot(all((cu[, 8] - cu[, 9]) >= 0))
    stopifnot(all((cu[, 10] - rowSums(cu[, c(11, 12)])) >= 0))
    
    # ## check infective states and new incidence are valid
    ## KEEPING FOR POSTERITY BUT DON'T NEED DUE TO IMPORTS
    # Einc <- cu[, 2] - cuprev[, 2]
    # print(which(Einc > 0))
    # stopifnot(all(rowSums(uprev[Einc > 0, c(3, 5, 6, 8)]) > 0))
}

## run model for each set of design points in pars
## pars: matrix / data frame of parameters in order: 
##       nu, nuA, pE, pEP, pA, pP, pI1, pI1H, pI1D, pI2, pH, pHD
##       each parameter except nu/nuA have entries for
##       each age-group e.g. nu, nuA, pE1, pE2, pE3 etc.
## C:    contact matrix for mixing between age-classes
## data: data frame of data to fit to, with columns
##       DI1, DI2, ..., DI8, DH1, ..., DH8
## u:    matrix of initial states (rows) in order:
##       S, E, A, RA, P, I1, DI, I2, RI, H, RH, DH
##       with age-groups along the columns
## ndays: the number of days to fit to
## npart: the number of particles
## MD:    turn model discrepancy process on (TRUE) or off (FALSE)
## obsScale: scaling parameter for Poisson observation process (see code)
## a1, a2, b: parameters for Skellam observation process
## saveAll: a logical specifying whether to return all states (if NA then returns just
##          log-likelihood estimate, if TRUE then returns all states, if FALSE then 
##          returns just observed states))
## PF:      a logical, if TRUE then does particle filtering, else just simulates
##          from model

PF1 <- function(pars, C, data, u, ndays, npart = 10, MD = TRUE, a1 = 0.01, a2 = 0.2, b = 0.1, a_dis = 0.05, b_dis = 0.05, saveAll = NA, PF = TRUE, ncores = detectCores()) {
    runs <- mclapply(1:nrow(pars), function(k, pars, C, u, npart, ndays, data, MD, a1, a2, b, a_dis, b_dis, saveAll, PF) {
        
        ## set initial log-likelihood
        ll <- 0
        
        ## set population size vector (used as simple check on processes)
        N <- colSums(u)
        
        ## set up particles
        nages <- ncol(u)
        u <- rep(list(u), npart)
        
        ## capture cumulative counts (used in MD process)
        cu <- u
        
        ## set up proposed particles
        disSims <- u
        
        ## vector of particle weights
        weights <- numeric(npart)
        
        ## set default for saveAll if PF = FALSE
        if(!PF & is.na(saveAll)) saveAll <- TRUE
        
        ## list to store outputs
        if(!is.na(saveAll)) out <- list()
        
        ## extract observations
        if(PF) {
            data <- select(data, t, (starts_with("DI") | starts_with("DH")) & ends_with("obs")) %>%
                {rbind(rep(0, ncol(.)), .)} %>%
                mutate(across(!t, ~. - lag(.))) %>%
                slice(-1) %>%
                as.matrix()
        }
        
        ## loop over time
        for(t in 1:ndays) {
            
            ## extract data
            if(PF) obsInc <- data[t, -1]
            
            ## loop over particles
            for(i in 1:npart) {
                
                ## run model
                temp <- discreteStochModel(unlist(pars[k, ]), t - 1, t, u[[i]], C)[2, -1]
                disSims[[i]] <- matrix(temp, ncol = length(N), byrow = TRUE)
                
                ## cols: c("S", "E", "A", "RA", "P", "I1", "DI", "I2", "RI", "H", "RH", "DH")
                ##          1,   2,   3,   4,    5,   6,    7,    8,    9,    10,  11,   12
                
                ## set weights
                weights[i] <- 0
                
                ## adjust states according to model discrepancy
                if(MD) {
                    if(PF) {
                        ## DH (MD on incidence)
                        DHinc <- disSims[[i]][12, ] - u[[i]][12, ]
                        for(j in 1:length(DHinc)) {
                            temp <- dtskellam(-DHinc[j]:(u[[i]][10, j] - DHinc[j]), a_dis + b_dis * DHinc[j], a_dis + b_dis * DHinc[j], -DHinc[j], u[[i]][10, j] - DHinc[j], log = TRUE)
                            tempcounts <- DHinc[j] + (-DHinc[j]:(u[[i]][10, j] - DHinc[j]))
                            obsDiffs <- obsInc[length(DHinc) + j] - tempcounts
                            temp <- temp + dtskellam(obsDiffs, a1 + b * tempcounts, a2 + b * tempcounts, LB = -tempcounts, UB = obsInc[length(DHinc) + j], log = TRUE)
                            tempnorm <- log_sum_exp(temp)
                            temp <- temp - tempnorm
                            tempind <- which(rmultinom(1, 1, exp(temp)) == 1)
                            DHinc[j] <- tempcounts[tempind]
                            weights[i] <- weights[i] + tempnorm
                        }
                    } else {
                        ## DH (MD on incidence)
                        DHinc <- disSims[[i]][12, ] - u[[i]][12, ]
                        DHinc <- DHinc + rtskellam(ncol(disSims[[i]]), a_dis + b_dis * DHinc, a_dis + b_dis * DHinc, -DHinc, u[[i]][10, ] - DHinc)
                    }
                    disSims[[i]][12, ] <- u[[i]][12, ] + DHinc
                    cu[[i]][12, ] <- cu[[i]][12, ] + DHinc
                    
                    ## RH given DH (MD on incidence)
                    RHinc <- disSims[[i]][11, ] - u[[i]][11, ]
                    RHinc <- RHinc + rtskellam(ncol(disSims[[i]]), a_dis + b_dis * RHinc, a_dis + b_dis * RHinc, -RHinc, u[[i]][10, ] - DHinc - RHinc)
                    disSims[[i]][11, ] <- u[[i]][11, ] + RHinc
                    cu[[i]][11, ] <- cu[[i]][11, ] + RHinc
                    
                    ## H
                    disSims[[i]][10, ] <- disSims[[i]][10, ] + rtskellam(ncol(disSims[[i]]), 
                         a_dis + b_dis * disSims[[i]][10, ], 
                         a_dis + b_dis * disSims[[i]][10, ],
                         -disSims[[i]][10, ] + u[[i]][10, ] - DHinc - RHinc,
                         u[[i]][6, ] - disSims[[i]][10, ] + u[[i]][10, ] - DHinc - RHinc)
                    Hinc <- disSims[[i]][10, ] - u[[i]][10, ] + DHinc + RHinc
                    cu[[i]][10, ] <- cu[[i]][10, ] + Hinc
                    
                    ## DI given H (MD on incidence)
                    if(PF) {
                        DIinc <- disSims[[i]][7, ] - u[[i]][7, ]
                        for(j in 1:length(DIinc)) {
                            temp <- dtskellam(-DIinc[j]:(u[[i]][6, j] - Hinc[j] - DIinc[j]), a_dis + b_dis * DIinc[j], a_dis + b_dis * DIinc[j], -DIinc[j], u[[i]][6, j] - Hinc[j] - DIinc[j], log = TRUE)
                            tempcounts <- DIinc[j] + (-DIinc[j]:(u[[i]][6, j] - Hinc[j] - DIinc[j]))
                            obsDiffs <- obsInc[j] - tempcounts
                            temp <- temp + dtskellam(obsDiffs, a1 + b * tempcounts, a2 + b * tempcounts, LB = -tempcounts, UB = obsInc[j], log = TRUE)
                            tempnorm <- log_sum_exp(temp)
                            temp <- temp - tempnorm
                            tempind <- which(rmultinom(1, 1, exp(temp)) == 1)
                            DIinc[j] <- tempcounts[tempind]
                            weights[i] <- weights[i] + tempnorm
                        }
                    } else {
                        ## DI given H (MD on incidence)
                        DIinc <- disSims[[i]][7, ] - u[[i]][7, ]
                        DIinc <- DIinc + rtskellam(ncol(disSims[[i]]), a_dis + b_dis * DIinc, a_dis + b_dis * DIinc, -DIinc, u[[i]][6, ] - Hinc - DIinc)
                        disSims[[i]][7, ] <- u[[i]][7, ] + DIinc
                        cu[[i]][7, ] <- cu[[i]][7, ] + DIinc
                    }
                    disSims[[i]][7, ] <- u[[i]][7, ] + DIinc
                    cu[[i]][7, ] <- cu[[i]][7, ] + DIinc            
                    
                    ## RI (MD on incidence)
                    RIinc <- disSims[[i]][9, ] - u[[i]][9, ]
                    RIinc <- RIinc + rtskellam(ncol(disSims[[i]]), a_dis + b_dis * RIinc, a_dis + b_dis * RIinc, -RIinc, u[[i]][8, ] - RIinc)
                    disSims[[i]][9, ] <- u[[i]][9, ] + RIinc
                    cu[[i]][9, ] <- cu[[i]][9, ] + RIinc
                    
                    ## I2 given H, RI and DI
                    disSims[[i]][8, ] <- disSims[[i]][8, ] + 
                        rtskellam(ncol(disSims[[i]]), a_dis + b_dis * disSims[[i]][8, ], 
                        a_dis + b_dis * disSims[[i]][8, ], 
                        -disSims[[i]][8, ] + u[[i]][8, ] - RIinc,
                        u[[i]][6, ] - DIinc - Hinc - disSims[[i]][8, ] + u[[i]][8, ] - RIinc)
                    I2inc <- disSims[[i]][8, ] - u[[i]][8, ] + RIinc
                    cu[[i]][8, ] <- cu[[i]][8, ] + I2inc
                    
                    ## I1 given later
                    disSims[[i]][6, ] <- disSims[[i]][6, ] + 
                        rtskellam(ncol(disSims[[i]]), a_dis + b_dis * disSims[[i]][6, ], 
                        a_dis + b_dis * disSims[[i]][6, ], 
                        -disSims[[i]][6, ] + u[[i]][6, ] - I2inc - Hinc - DIinc,
                        u[[i]][5, ] - disSims[[i]][6, ] + u[[i]][6, ] - I2inc - Hinc - DIinc)
                    I1inc <- disSims[[i]][6, ] - u[[i]][6, ] + DIinc + Hinc + I2inc
                    cu[[i]][6, ] <- cu[[i]][6, ] + I1inc
                    
                    ## P given later
                    disSims[[i]][5, ] <- disSims[[i]][5, ] + 
                        rtskellam(ncol(disSims[[i]]), a_dis + b_dis * disSims[[i]][5, ], 
                        a_dis + b_dis * disSims[[i]][5, ], 
                        -disSims[[i]][5, ] + u[[i]][5, ] - I1inc,
                        u[[i]][2, ] - disSims[[i]][5, ] + u[[i]][5, ] - I1inc)
                    Pinc <- disSims[[i]][5, ] - u[[i]][5, ] + I1inc
                    cu[[i]][5, ] <- cu[[i]][5, ] + Pinc
                    
                    ## RA (MD on incidence)
                    RAinc <- disSims[[i]][4, ] - u[[i]][4, ]
                    RAinc <- RAinc + rtskellam(ncol(disSims[[i]]), a_dis + b_dis * RAinc, a_dis + b_dis * RAinc, -RAinc, u[[i]][3, ] - RAinc)
                    disSims[[i]][4, ] <- u[[i]][4, ] + RAinc
                    cu[[i]][4, ] <- cu[[i]][4, ] + RAinc
                    
                    ## A given later
                    disSims[[i]][3, ] <- disSims[[i]][3, ] + 
                        rtskellam(ncol(disSims[[i]]), a_dis + b_dis * disSims[[i]][3, ], 
                        a_dis + b_dis * disSims[[i]][3, ], 
                        -disSims[[i]][3, ] + u[[i]][3, ] - RAinc,
                        u[[i]][2, ] - Pinc - disSims[[i]][3, ] + u[[i]][3, ] - RAinc)
                    Ainc <- disSims[[i]][3, ] - u[[i]][3, ] + RAinc
                    cu[[i]][3, ] <- cu[[i]][3, ] + Ainc
                    
                    ## E
                    disSims[[i]][2, ] <- disSims[[i]][2, ] + 
                        rtskellam(ncol(disSims[[i]]), a_dis + b_dis * disSims[[i]][2, ], a_dis + b_dis * disSims[[i]][2, ], 
                        -disSims[[i]][2, ] + u[[i]][2, ] - Ainc - Pinc,
                        u[[i]][1, ] - disSims[[i]][2, ] + u[[i]][2, ] - Ainc - Pinc)
                    Einc <- disSims[[i]][2, ] - u[[i]][2, ] + Ainc + Pinc
                    cu[[i]][2, ] <- cu[[i]][2, ] + Einc
                    
                    ## S
                    disSims[[i]][1, ] <- u[[i]][1, ] - Einc
                    cu[[i]][1, ] <- cu[[i]][1, ] - Einc
                } else {
                    ## DH
                    DHinc <- disSims[[i]][12, ] - u[[i]][12, ]
                    cu[[i]][12, ] <- cu[[i]][12, ] + DHinc
                    
                    ## RH 
                    RHinc <- disSims[[i]][11, ] - u[[i]][11, ]
                    cu[[i]][11, ] <- cu[[i]][11, ] + RHinc
                    
                    ## H
                    Hinc <- disSims[[i]][10, ] - u[[i]][10, ] + DHinc + RHinc
                    cu[[i]][10, ] <- cu[[i]][10, ] + Hinc
                    
                    ## DI 
                    DIinc <- disSims[[i]][7, ] - u[[i]][7, ]
                    cu[[i]][7, ] <- cu[[i]][7, ] + DIinc
                    
                    ## RI
                    RIinc <- disSims[[i]][9, ] - u[[i]][9, ]
                    cu[[i]][9, ] <- cu[[i]][9, ] + RIinc
                    
                    ## I2 
                    I2inc <- disSims[[i]][8, ] - u[[i]][8, ] + RIinc
                    cu[[i]][8, ] <- cu[[i]][8, ] + I2inc
                    
                    ## I1 
                    I1inc <- disSims[[i]][6, ] - u[[i]][6, ] + DIinc + Hinc + I2inc
                    cu[[i]][6, ] <- cu[[i]][6, ] + I1inc
                    
                    ## P given later
                    Pinc <- disSims[[i]][5, ] - u[[i]][5, ] + I1inc
                    cu[[i]][5, ] <- cu[[i]][5, ] + Pinc
                    
                    ## RA
                    RAinc <- disSims[[i]][4, ] - u[[i]][4, ]
                    cu[[i]][4, ] <- cu[[i]][4, ] + RAinc
                    
                    ## A 
                    Ainc <- disSims[[i]][3, ] - u[[i]][3, ] + RAinc
                    cu[[i]][3, ] <- cu[[i]][3, ] + Ainc
                    
                    ## E
                    Einc <- disSims[[i]][2, ] - u[[i]][2, ] + Ainc + Pinc
                    cu[[i]][2, ] <- cu[[i]][2, ] + Einc
                    
                    ## S
                    disSims[[i]][1, ] <- u[[i]][1, ] - Einc
                    cu[[i]][1, ] <- cu[[i]][1, ] - Einc
                    
                    if(PF) {
                        ## calculate log observation error weights
                        obsDiffs <- obsInc - c(DIinc, DHinc)
                        weights[i] <- sum(dtskellam(obsDiffs, a1 + b * c(DIinc, DHinc), a2 + b * c(DIinc, DHinc), LB = -c(DIinc, DHinc), UB = obsInc, log = TRUE))
                    }
                }
                # ## check counts
                # checkCounts(disSims[[i]], cu[[i]], N)
            }
            if(!is.na(saveAll)) {
                if(saveAll) {
                    out[[t]] <- disSims
                } else {
                    out[[t]] <- map(disSims, ~.[c(7, 12), ])
                }
            }
            
            if(PF) {
                ## calculate log-likelihood contribution
                ll <- ll + log_sum_exp(weights, mn = TRUE)
                
                ## if zero likelihood then return
                if(!is.finite(ll)) {
                    if(!is.na(saveAll)) {
                        return(list(ll = ll, particles = out))
                    } else {
                        return(ll)
                    }
                }
                
                ## normalise weights
                weights <- exp(weights - log_sum_exp(weights))
                
                ## resample
                inds <- apply(rmultinom(npart, 1, weights), 2, function(x) which(x == 1))
            } else {
                inds <- 1:npart
            }
            u <- disSims[inds]
            cu <- cu[inds]
        }
        if(!is.na(saveAll)) {
            if(PF) {
                return(list(ll = ll, particles = out))
            } else {
                return(list(particles = out))
            }
        } else {
            return(ll)
        }
    }, pars = pars, C = C, u = u, npart = npart, ndays = ndays, data = data, MD = MD, 
       a1 = a1, a2 = a2, b = b, a_dis = a_dis, b_dis = b_dis, saveAll = saveAll, PF = PF, mc.cores = ncores)
    if(!is.na(saveAll)) {
        if(PF) {
            ll <- map(runs, "ll")
            runs <- map(runs, "particles")
            ll <- do.call("c", ll)
            return(list(ll = ll, particles = runs))
        } else {
            runs <- map(runs, "particles")
            return(list(particles = runs))
        }
    } else {
        ll <- do.call("c", runs)
        return(ll)
    }
}
