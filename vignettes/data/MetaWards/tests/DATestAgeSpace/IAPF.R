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
## obsScale: scaling parameter for Skellam observation process (see code)
## a1, a2, b: parameters for Skellam observation process
## pmix: mixing proportion for twisting functions
## saveAll: a logical specifying whether to return all states (if FALSE then returns just observed states))
## writeExt: a logical denoting whether to save particles externally or not
## PF:      a logical denoting whether to run a particle filter, or just simulate from the model
## ncores:  the number of cores for OpenMP parallelisation (if NA then defaults to all available cores)

IAPF <- function(pars, C, data, u1_moves, u1, ndays, npart = 10, kstop = 3, kmax = 3, tau = 1, a1 = 0.01, a2 = 0.2, b = 0.1, 
               a_dis = 0.05, b_dis = 0.5, pmix = 0.9, saveAll = NA, writeExt = FALSE, PF = TRUE, ncores = NA) {
               
    ## set default for saveAll if PF = FALSE
    if(!PF & is.na(saveAll)) saveAll <- TRUE
    
    ## convert indicators
    if(is.na(saveAll)) {
        saveAllint <- 0
    } else {
        saveAllint <- ifelse(saveAll, 2, 1)
    }
    PFint <- ifelse(PF, 1, 0)
    writeExtint <- ifelse(writeExt, 1, 0)
    
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
    
    ## set up output folder
    if(writeExt) {
        print("Reminder to write code to pass save folder out")
        if(dir.exists("saveOut")) system("rm -rf saveOut")
        dir.create("saveOut")
    }
    
    ## run particle filter for each set of inputs
    runs <- lapply(1:nrow(pars), function(k, pars, C, u1_moves, ncohorts, u1, npart, kstop, tau, ndays, data, a1, a2, b, a_dis, b_dis, pmix, saveAll, writeExt, PF, ncores) {
        
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
        
        ## if just simulating
        if(PF == 0) {
            ## check return outputs
            if(saveAll == 0) stop("Must set 'saveAll' to something if not running a PF")
            
            ## run particle filter
            particles <- TPF_cpp(pars, C, data, 12L, 8L, 339L, u1_moves, ncohorts, u1, ndays,
                npart, matrix(0, 1, 1), matrix(1, 1, 1), pmix, 0, a1, a2, b, a_dis, b_dis, saveAll, writeExt, 0, PF, ncores)
            return(particles)
        }
        
        ## run particle filter
        particles <- TPF_cpp(pars, C, data, 12L, 8L, 339L, u1_moves, ncohorts, u1, ndays,
            npart, matrix(0, 1, 1), matrix(1, 1, 1), pmix, 0, a1, a2, b, a_dis, b_dis, saveAll, writeExt, 1, PF, ncores)
            
        ## extract particles
        ll <- particles$ll
        if(saveAll != 0 & writeExt == 0) saveParticles <- particles$particles
        particles <- particles$psi
            
        ## set regularised Gaussian functions
        regNorm <- function(pars, eta, psi_it) {
            if(any(pars[-1] <= 0)) return(NA)
            x1 <- dnorm(eta, mean = pars[1], sd = pars[2])
            x2 <- pars[3] * exp(psi_it)
            sum((x1 - x2)^2)
        }
        
        ## set update loop
        kcurr <- 1
        valid <- 0
        while(valid == 0) {
            
            ## check
            if(kcurr >= kstop) {
                        
                ## check log-likelihoods
                templl <- ll[(kcurr - kstop + 1):kcurr]
                mnll <- log_sum_exp(templl, mn = TRUE)
                sdll <- max(c(mnll, templl))
                sdll <- 0.5 * (2 * sdll - log(length(templl) - 1) + log(sum((exp(templl - sdll) - exp(mnll - sdll))^2)))
                cv <- round(exp(sdll - mnll), 2)
                cat(paste0("Previous ", kstop, " log-likelihood estimates: ", paste0(round(templl, 2), collapse = ", "), "\n"))
                cat(paste0("CV of likelihood estimate = ", cv, " tau = ", tau, "\n"))
                ## if coefficient of variation < tau, then exit
                if(cv < tau | kcurr >= kmax) {
                    valid <- 1
                } else {
                    ## if log-likelihoods not monotonically increasing
                    ## then increase the number of particles
                    if(!all(diff(templl) > 0) & npart[(kcurr - kstop + 1)] == npart[kcurr]) {
                        npart <- c(npart, npart[kcurr] * 2)
                    } else {
                        npart <- c(npart, npart[kcurr])
                    }
                }
            } else {
                npart <- c(npart, npart[kcurr])
            }
            
            if(valid == 0) {
            
                ## now calculate twisting functions
                psi <- list()
                
                ## print progress        
                cat(paste0("\nCalculating twisting functions: IPF iteration = ", kcurr, "\n"))
                ptm <- proc.time()
                
                ## calculate observation likelihoods
                temp <- particles[(npart[kcurr] * (ndays - 1)  * 4 + 1):(npart[kcurr] * ndays * 4), , ]
                temp <- map(1:npart[kcurr], function(i, x, data, a1, a2, b) {
                    ## extract particle
                    x <- x[((i - 1) * 4 + 1):(i * 4), , ]
                
                    ## reorder particles to match data
                    x <- as.vector(aperm(x[1:2, , ], 3:1))
                    
                    ## calculate observation densities
                    rbind(x, dtskellam(data - x, a1 + b * x, a2 + b * x, -x, log = TRUE))
                }, x = temp, data = data[ndays, ], a1 = a1, a2 = a2, b = b)
                temp <- do.call("rbind", temp)
                
                ## optimise twisting functions
                output <- mclapply(1:ncol(temp), function(i, x, fn) {
                    x <- x[, i]
                    eta <- x[seq(1, length(x) - 1, by = 2)]
                    psi_it <- x[seq(2, length(x), by = 2)]
                    temp <- optim(c(mean(eta) + 0.1, sd(eta) + 0.1, 1), fn, eta = eta, psi_it = psi_it, control = list(maxit = 5000))
                    k <- 1
                    while(temp$convergence != 0 & k < 10) {
                        temp <- optim(temp$par, fn, eta = eta, psi_it = psi_it, control = list(maxit = 5000))
                        k <- k + 1
                    }
                    if(any(is.na(temp$par))) browser()
                    temp
                }, x = temp, fn = regNorm, mc.cores = ncores)
                
                ## check convergence
                conv <- map_int(output, "convergence")
                if(any(conv > 0)) print(table(conv))
                
                ## extract Gaussian parameters
                munorm <- map(output, "par")
                varnorm <- (map_dbl(munorm, 2))^2
                munorm <- map_dbl(munorm, 1)
                
                ## save rates
                psi[[ndays]] <- list(munorm = munorm, varnorm = varnorm)
                    
                cat(paste0("Day: ", ndays, " t = ", round(as.numeric((proc.time() - ptm)["elapsed"]), 2), " secs\n"))
                ptm <- proc.time()
                
                ## now loop over remaining time points
                for(t in (ndays - 1):1) {
                
                    ## generate filter rates
                    temp <- particles[(npart[kcurr] * (t - 1)  * 4 + 1):(npart[kcurr] * t * 4), , ]
                    temp <- TPF_rates_cpp(pars, data[t, ], 8, 339, temp, npart[kcurr], munorm, varnorm, pmix, a1, a2, b, a_dis, b_dis, ncores)
                    
                    if(any(!is.finite(temp))) browser()
                
                    ## optimise twisting functions
                    output <- mclapply(1:ncol(temp), function(i, x, fn) {
                        x <- x[, i]
                        eta <- x[seq(1, length(x) - 1, by = 2)]
                        psi_it <- x[seq(2, length(x), by = 2)]
                        temp <- try(optim(c(mean(eta) + 0.1, sd(eta) + 0.1, 1), fn, eta = eta, psi_it = psi_it, control = list(maxit = 5000)), silent = TRUE)
                        k <- 1
                        while(class(temp) == "try-error" & k < 1000) {
                            temp <- try(optim(c(rnorm(1, 0, 10), rexp(1, 0.1), rexp(1, 100)), fn, eta = eta, psi_it = psi_it, control = list(maxit = 5000)), silent = TRUE)
                            k <- k + 1
                        }
                        k <- 1
                        while(temp$convergence != 0 & k < 10) {
                            temp <- optim(temp$par, fn, eta = eta, psi_it = psi_it, control = list(maxit = 5000))
                            k <- k + 1
                        }
                       if(any(is.na(temp$par))) browser()
                        temp
                    }, x = temp, fn = regNorm, mc.cores = ncores)
                    
                    ## check convergence
                    conv <- map_int(output, "convergence")
                    if(any(conv > 0)) print(table(conv))
                    
                    ## extract Gaussian parameters
                    munorm <- map(output, "par")
                    varnorm <- (map_dbl(munorm, 2))^2
                    munorm <- map_dbl(munorm, 1)
                    
                    ## save rates
                    psi[[t]] <- list(munorm = munorm, varnorm = varnorm)
                    
                    cat(paste0("Day: ", t, " t = ", round(as.numeric((proc.time() - ptm)["elapsed"]), 2), " secs\n"))
                    ptm <- proc.time()
                }
                
                ## reduce to rate matrix
                psimu <- do.call("rbind", map(psi, "munorm"))
                psivar <- do.call("rbind", map(psi, "varnorm"))
                
                ## print summaries to screen
                cat("\nPsi means:\n")
                temp <- apply(psimu, 2, min)
                print(summary(psimu[, which(temp == min(temp))[1]]))
                temp <- apply(psimu, 2, max)
                print(summary(psimu[, which(temp == max(temp))[1]]))
                cat("\nPsi SDs:\n")
                temp <- apply(psivar, 2, min)
                print(summary(sqrt(psivar[, which(temp == min(temp))[1]])))
                temp <- apply(psivar, 2, max)
                print(summary(sqrt(psivar[, which(temp == max(temp))[1]])))
                
                ## garbage collect just in case
                gc()
                
                ## update counter
                kcurr <- kcurr + 1
                
                ## run particle filter
                cat(paste0("\nRunning model (npart = ", npart[kcurr], ")\n"))
                particles <- TPF_cpp(pars, C, data, 12L, 8L, 339L, u1_moves, ncohorts, u1, ndays,
                    npart[kcurr], psimu, psivar, pmix, 1, a1, a2, b, a_dis, b_dis, saveAll, writeExt, 1, PF, ncores)
            
                ## extract particles
                ll <- c(ll, particles$ll)
                if(saveAll != 0 & writeExt == 0) saveParticles <- particles$particles
                if(saveAll != 0 & writeExt == 1) {
                    if(dir.exists(paste0("saveOut_", kcurr))) system(paste0("rm -r saveOut_", kcurr))
                    system(paste0("cp -rf saveOut saveOut_", kcurr))
                }
                particles <- particles$psi
                print(ll)
            }
        }
        if(saveAll != 0 & writeExt == 0) {
            return(list(ll = ll[length(ll)], particles = saveParticles))
        } else {
            return(list(ll = ll[length(ll)]))
        }
    }, pars = pars, C = C, u1_moves = u1_moves, ncohorts = ncohorts, u1 = u1, npart = npart, kstop = kstop, tau = tau, ndays = ndays, data = data, a1 = a1, a2 = a2, b = b, a_dis = a_dis, b_dis = b_dis, pmix = pmix, saveAll = saveAllint, writeExt = writeExtint, PF = PFint, ncores = ncores)
    if(!is.na(saveAll)) {
        if(!writeExt) {
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
            if(PF) {
                ll <- map(runs, "ll")
                ll <- do.call("c", ll)
                return(list(ll = ll))
            } else {
                return(NULL)
            }
        }       
    } else {
        ll <- map(runs, "ll")
        ll <- do.call("c", runs)
        return(ll)
    }
}
