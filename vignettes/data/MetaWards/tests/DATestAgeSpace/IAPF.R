## log-sum-exp function
log_sum_exp <- function(x, mn = FALSE) {
    maxx <- max(x)
    y <- maxx + log(sum(exp(x - maxx)))
    if(mn) y <- y - log(length(x))
    y
}

## run model for each set of design points in pars
## pars: matrix / data frame of parameters in order: 
##       nu, nuA, pE, pEP, pA, pP, pI1, pI1H, pI1D, pI2, pH, pHD
##       each parameter except nu/nuA have entries for
##       each age-group e.g. nu, nuA, pE1, pE2, pE3 etc.
## C:    contact matrix for mixing between age-classes
## data: data frame of data to fit to, with columns
##       DI1, DI2, ..., DI8, DH1, ..., DH8
## u1_moves: matrix of commuter data, with columns:
##            home LAD, work LAD
## u1: 3D array of commuter data, with dimensions:
##            nclasses x nages x nrow(u1_moves)
## ndays: the number of days to fit to
## npart: the number of particles
## kstop: stopping criteria for IAPF
## obsScale: scaling parameter for Poisson observation process (see code)
## a1, a2, b: parameters for Skellam observation process
## saveAll: a logical specifying whether to return all states (if FALSE then returns just observed states))
## PF:      a logical denoting whether to run a particle filter, or just simulate from the model
## ncores:  the number of cores for OpenMP parallelisation (if NA then defaults to all available cores)

IAPF <- function(pars, C, data, u1_moves, u1, ndays, npart = 10, kstop = 3, tau = 1, a1 = 0.01, a2 = 0.2, b = 0.1, 
               a_dis = 0.05, b_dis = 0.5, saveAll = NA, PF = TRUE, ncores = NA) {
               
    ## set default for saveAll if PF = FALSE
    if(!PF & is.na(saveAll)) saveAll <- TRUE
    
    ## convert indicators
    if(is.na(saveAll)) {
        saveAllint <- 0
    } else {
        saveAllint <- ifelse(saveAll, 2, 1)
    }
    PFint <- ifelse(PF, 1, 0)
    
    ## check number of requested cores
    if(is.na(ncores)) {
        ncores <- parallel::detectCores()
    }
    
    ## set up dummy "data" if just simulation required
    if(!PF) {
        data <- matrix(NA, 1, 1)
    }
    
    ## check u1_moves are ordered
    u1_moves <- u1_moves[sort.list(u1_moves[, 1]), ]
    
    ## generate number of cohorts
    ncohorts <- tapply(u1_moves[, 1], u1_moves[, 1], length)
    ncohorts <- c(0, cumsum(ncohorts))
    names(ncohorts) <- NULL
    
    print("Reminder to write code to not hard-code sizes of objects and data")
    
    ## run particle filter for each set of inputs
    runs <- lapply(1:nrow(pars), function(k, pars, C, u1_moves, ncohorts, u1, npart, kstop, tau, ndays, data, a1, a2, b, a_dis, b_dis, saveAll, PF, ncores) {
        
        if(PF == 1) {
            ## extract observations
            data <- select(data, t, (starts_with("DI") | starts_with("DH")) & ends_with("obs")) %>%
                {rbind(rep(0, ncol(.)), .)} %>%
                mutate(across(!t, ~. - lag(.))) %>%
                slice(-1) %>%
                as.matrix()
        
            ## remove time since not used in code
            data <- data[, -1]
        } else {
            ## dummy matrix
            data <- matrix(0, 1, 1)
        }
        
        ## set pars
        pars <- unlist(pars[k, ])
    
        ## do garbage collection (seems to solve allocation issue)
        gc()
        
        ## run particle filter
        particles <- TPF_cpp(pars, C, data, 12L, 8L, 339L, u1_moves, ncohorts, u1, ndays,
            npart, matrix(0, 1, 1), 0, a1, a2, b, a_dis, b_dis, 0, 1, PF, ncores)
            
        ## extract particles
        ll <- particles$ll
        particles <- particles$psi
            
        ## set regularised Poisson functions
        regPois <- function(pars, eta, psi_it) {
            if(any(pars <= 0)) return(NA)
            x1 <- dpois(eta, lambda = pars[1], log = TRUE)
            x2 <- log(pars[2]) + psi_it
            mx <- pmax(x1, x2)
            x <- 2 * mx + log((exp(x1 - mx) - exp(x2 - mx))^2)
            if(max(x) == -Inf) return(-Inf)
            log_sum_exp(x)
        }
        
        ## set update loop
        kcurr <- 2
        valid <- 0
        npart1 <- npart
        while(valid == 0) {
        
            ## now calculate twisting functions
            psi <- list()
            
            ## print progress        
            cat(paste0("Calculating twisting functions: IPF iteration = ", kcurr, "\n"))
            ptm <- proc.time()
            
            ## calculate observation likelihoods
            temp <- particles[(npart * (ndays - 1)  * 4 + 1):(npart * ndays * 4), , ]
            temp <- map(1:npart, function(i, x, data, a1, a2, b) {
                ## extract particle
                x <- x[((i - 1) * 4 + 1):(i * 4), , ]
            
                ## reorder particles to match data
                x <- as.vector(aperm(x[1:2, , ], 3:1))
                
                ## calculate observation densities
                rbind(x, dtskellam(data - x, a1 + b * x, a2 + b * x, -x, data, log = TRUE))
            }, x = temp, data = data[ndays, ], a1 = a1, a2 = a2, b = b)
            temp <- do.call("rbind", temp)
            
            ## optimise twisting functions
            output <- mclapply(1:ncol(temp), function(i, x, fn) {
                x <- x[, i]
                eta <- x[seq(1, length(x) - 1, by = 2)]
                psi_it <- x[seq(2, length(x), by = 2)]
                temp <- optim(c(mean(eta) + 0.1, 1), fn, eta = eta, psi_it = psi_it, control = list(maxit = 5000))
                if(temp$convergence != 0) temp <- optim(temp$par, fn, eta = eta, psi_it = psi_it, control = list(maxit = 5000))
                if(any(is.na(temp$par)) | temp$convergence != 0) browser()
                temp
            }, x = temp, fn = regPois, mc.cores = ncores)
            
            ## check convergence
            conv <- map_int(output, "convergence")
            if(any(conv > 0)) print(table(conv)) #stop(paste0("optim not converged - ", ndays))
            
            ## extract Poisson parameters
            rates <- map(output, "par")
            rates <- map_dbl(rates, 1)
            
            ## save rates
            psi[[ndays]] <- rates
                
            cat(paste0("Day: ", ndays, " t = ", round(as.numeric((proc.time() - ptm)["elapsed"]), 2), " secs\n"))
            ptm <- proc.time()
            
            ## now loop over remaining time points
            for(t in (ndays - 1):1) {
            
                ## generate filter rates
                temp <- particles[(npart * (t - 1)  * 4 + 1):(npart * t * 4), , ]
                temp <- TPF_rates_cpp(pars, data[t, ], 8, 339, temp, npart, rates, a1, a2, b, a_dis, b_dis, ncores)
#                browser()
            
                ## optimise twisting functions
                output <- mclapply(1:ncol(temp), function(i, x, fn) {
                    x <- x[, i]
                    eta <- x[seq(1, length(x) - 1, by = 2)]
                    psi_it <- x[seq(2, length(x), by = 2)]
                    ## guess initial conditions
#                    mn <- mean(eta) + 0.1
#                    ini <- (eta - mn)^2
#                    ini <- which(ini == min(ini))[1]
#                    ini <- c(mn, exp(dpois(eta[ini], lambda = mn, log = TRUE) - psi_it[ini]))
#                    print(ini)
#                    print(fn(ini, eta, psi_it))
#                    if(is.na(fn(ini, eta, psi_it))) browser()
                    temp <- optim(c(mean(eta) + 0.1, 1), fn, eta = eta, psi_it = psi_it, control = list(maxit = 5000))
                    k <- 1
                    while(temp$convergence != 0 & k < 3) {
                        temp <- optim(temp$par, fn, eta = eta, psi_it = psi_it, control = list(maxit = 5000))
                        k <- k + 1
                    }
#                     else {
#                        if(any(is.na(temp$par)) | temp$convergence != 0) browser()
#                    }
                    temp
                }, x = temp, fn = regPois, mc.cores = ncores)
                
                ## check convergence
                conv <- map_int(output, "convergence")
                if(any(conv > 0)) print(table(conv)) #browser()#stop(paste0("optim not converged - ", t))
                
                ## extract Poisson parameters
                rates <- map(output, "par")
                rates <- map_dbl(rates, 1)
                
                ## save rates
                psi[[t]] <- rates
                
                cat(paste0("Day: ", t, " t = ", round(as.numeric((proc.time() - ptm)["elapsed"]), 2), " secs\n"))
                ptm <- proc.time()
            }
            
            ## reduce to rate matrix
            psi <- do.call("rbind", psi)
            
            ## garbage collect just in case
            gc()
            
            ## run particle filter
            npart <- npart1
            cat(paste0("\nRunning model (npart = ", npart, ")\n"))
            particles <- TPF_cpp(pars, C, data, 12L, 8L, 339L, u1_moves, ncohorts, u1, ndays,
                npart, psi, 1, a1, a2, b, a_dis, b_dis, 0, 1, PF, ncores)
            
            ## extract particles
            ll <- c(ll, particles$ll)
            particles <- particles$psi
            print(ll)
            
            ## check
            npart1 <- npart
            if(kcurr >= kstop) {
#                browser()
            
                ## check log-likelihoods
                templl <- ll[(kcurr - kstop + 1):kcurr]
                mnll <- log_sum_exp(templl, mn = TRUE)
                sdll <- max(c(mnll, templl))
                sdll <- 0.5 * (2 * sdll - log(length(templl) - 1) + log(sum((exp(templl - sdll) - exp(mnll - sdll))^2)))
                cv <- round(exp(sdll - mnll), 2)
                cat(paste0("Previous ", kstop, " log-likelihood estimates: ", paste0(round(templl, 2), collapse = ", "), "\n"))
                cat(paste0("CV of likelihood estimate = ", cv, " tau = ", tau, "\n"))
                ## if coefficient of variation < tau, then exit
                if(cv < tau) {
                    valid <- 1
                } else {
                    ## if log-likelihoods not monotonically increasing
                    ## then increase the number of particles
                    if(!all(diff(templl) < 0)) {
                        npart1 <- npart * 2
                    }
                }
            }   
            
            ## update counter
            kcurr <- kcurr + 1
        }
        if(saveAllint != 0) {
            particles <- TPF_cpp(pars, C, data, 12L, 8L, 339L, u1_moves, ncohorts, u1, ndays,
                npart, psi, 1, a1, a2, b, a_dis, b_dis, saveAllint, 1, PF, ncores)
            ll <- particles$ll
            particles <- particles$particles
        }
        list(ll = ll[length(ll)], particles = particles)
    }, pars = pars, C = C, u1_moves = u1_moves, ncohorts = ncohorts, u1 = u1, npart = npart, kstop = kstop, tau = tau, ndays = ndays, data = data, a1 = a1, a2 = a2, b = b, a_dis = a_dis, b_dis = b_dis, saveAll = saveAllint, PF = PFint, ncores = ncores)
    browser()
    
    if(!is.na(saveAll)) {
        if(PF) {
            ll <- map(runs, "ll")
            runs <- map(runs, "particles") %>%
                map(function(runs, ndays, npart) {
                    x <- list()
                    for(i in 1:ndays) {
                        x[[i]] <- runs[(i - 1) * npart + 1:npart]
                    }
                    x
                }, ndays = ndays, npart = npart)
            ll <- do.call("c", ll)
            return(list(ll = ll, particles = runs))
        } else {
            runs <- map(runs, "particles") %>%
                map(function(runs, ndays, npart) {
                    x <- list()
                    for(i in 1:ndays) {
                        x[[i]] <- runs[(i - 1) * npart + 1:npart]
                    }
                    x
                }, ndays = ndays, npart = npart)
            return(list(particles = runs))
        }
    } else {
        ll <- map(runs, "ll")
        ll <- do.call("c", runs)
        return(ll)
    }
}
