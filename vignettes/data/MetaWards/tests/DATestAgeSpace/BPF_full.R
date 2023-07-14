## log-sum-exp function
log_sum_exp <- function(x, mn = FALSE) {
    maxx <- max(x)
    y <- maxx + log(sum(exp(x - maxx)))
    if(mn) y <- y - log(length(x))
    y
}

## run model for each set of design points in pars
## pars: matrix / data frame of parameters in order: 
##       nu, nuA, pE, pEP, pA, pP, pI1, pI1H, pI1D, pI2, pH, pHD, beta_scale, p_move
##       each parameter except nu/nuA/beta_scale/p_move have entries for
##       each age-group e.g. nu, nuA, pE1, pE2, pE3 etc.
## C1:    contact matrix for mixing between age-classes
## C2:    contact matrix for mixing between age-classes post lockdown
## lockdown_day: the day of first lockdown
## cumDH/DI_age_lad: matrix of form D_a1_l1, D_a2_l1, ..., D_a8_lL along
##               columns with cumulative deaths per age/lad
## H_age_lad: matrix of form H_a1_l1, H_a2_l1, ..., H_a8_lL along
##               columns with hospital cases per age/lad
## lookup: a lookup table mapping lads to different spatial hierarchies
##               of form: LAD, death ID, region ID, NHS region ID
## u: list with elements:
##      u1_moves: matrix of commuter data, with columns:
##              home LAD, work LAD
##      u1: 3D array of commuter data, with dimensions:
##            nclasses x nages x nrow(u1_moves)
##      u2_moves: matrix of commuter data, with columns:
##            home LAD, play LAD, prob of moving
##      u2: 3D array of player data, with dimensions:
##            nclasses x nages x nlads
##  or a list with elements (generated from a previous run with snapshot = TRUE):
##      u1_moves: matrix of commuter data, with columns:
##              home LAD, work LAD
##      u1: list of length npart with each element a
##            3D array of commuter data, with dimensions:
##            nclasses x nages x nrow(u1_moves)
##      u2: list of length npart with each element a
##            3D array of commuter data, with dimensions:
##            nclasses x nages x nlads
##      playprobs: a vector of play probabilities
##      ncohorts1: a vector of movement lookups
##      ncohorts2: a vector of movement lookups
## tstart: when to start the outbreak
## tstop: when to stop the outbreak
## npart: the number of particles
## niter: number of Metropolis-Hastings steps to counter particle impoverishment
## a1, a2, b1, b2: parameters for Skellam observation process
## a_dis: intercept for Skellam MD process (can be nages x nlads matrix)
## b_dis: multiplicative scaling for truncated Gaussian MD process (can be nages x nlads matrix)
## saveAll: a logical specifying whether to return all states (if FALSE then returns just observed states))
## writeExt: a logical denoting whether to save particles externally or not
## snapshot: a logical determining whether to save a snapshot of the whole system at 'tstop' 
##           (for use in forwards simulations)
## PF:      a logical denoting whether to run a particle filter, or just simulate from the model
## ncores:  the number of cores for OpenMP parallelisation (if NA then defaults to all available cores)
## parEnsemble: decides whether to parallelise across or within ensemble

BPF <- function(pars, C1, C2, lockdown_day, cumDI_age_lad, cumDH_age_lad, H_age_lad,
        lookup, u, tstart, tstop, npart = 10, niter = 0, 
        a1 = 0, a2 = 0, b1 = 0.1, b2 = 0.1, a_dis = 0.05, b_dis = 0.01,
        saveAll = NA, writeExt = FALSE, snapshot = FALSE, outputName = "saveOut", PF = TRUE, 
        ncores = NA, parEnsemble = FALSE) {
               
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
    snapshotint <- ifelse(snapshot, 1, 0)
    
    ## check number of requested cores
    if(is.na(ncores)) {
        ncores <- parallel::detectCores()
    }
    
    ## set up dummy "data" if just simulation required
    if(!PF) {
        cumDI_age_lad <- matrix(NA, 1, 1)
        cumDH_age_lad <- matrix(NA, 1, 1)
        H_age_lad <- matrix(NA, 1, 1)
    } else {
        ## extract relevant time points to filter against
        cumDI_age_lad <- filter(cumDI_age_lad, t >= tstart & t <= tstop)
        cumDH_age_lad <- filter(cumDH_age_lad, t >= tstart & t <= tstop)
        H_age_lad <- filter(H_age_lad, t >= tstart & t <= tstop)
        stopifnot(all((cumDI_age_lad$t - tstart:tstop) == 0))
        stopifnot(all((cumDH_age_lad$t - tstart:tstop) == 0))
        stopifnot(all((H_age_lad$t - tstart:tstop) == 0))
    }
    
    if(!is.list(u)) stop("'u' is not a list")
    if(!"u1" %in% names(u)) stop("'u' must have 'u1' element")
    if(!is.list(u$u1)) {
        if(!all(c("u1", "u2", "u1_moves", "u2_moves") %in% names(u))) {
            stop("'u' must have elements: 'u1', 'u2', 'u1_moves', 'u2_moves'")
        }
        u1_moves <- u$u1_moves
        u2_moves <- u$u2_moves
        u1 <- u$u1
        u2 <- u$u2
        if(!is.matrix(u1_moves) & ncol(u1_moves) != 2) stop("'u$u1_moves' must be matrix with two columns")
        if(!is.matrix(u2_moves) & ncol(u2_moves) != 2) stop("'u$u2_moves' must be matrix with two columns")
        if(!is.array(u1) & length(dim(u1)) != 3) stop("'u$u1' must be 3D array")
        if(!is.array(u2) & length(dim(u2)) != 3) stop("'u$u2' must be 3D array")
    
        ## check no missing lads
        if(
            !identical(sort(unique(u1_moves[, 1])), sort(unique(u1_moves[, 1]))) |
            !identical(sort(unique(u1_moves[, 1])), sort(unique(u1_moves[, 2]))) |
            !identical(sort(unique(u1_moves[, 1])), sort(unique(u2_moves[, 1]))) |
            !identical(sort(unique(u1_moves[, 1])), sort(unique(u2_moves[, 2])))
        ) {
            stop("Can't match LADS")
        }
        if(any((1:length(unique(u1_moves[, 1])) - sort(unique(u1_moves[, 1]))) != 0)) stop("Can't match all LADs")
      
        ## check probabilities
        if(!all((tapply(u2_moves[, 3], u2_moves[, 1], sum) - 1) < 1e-15)) stop("u2 probs don't sum to one")
        ## check entries
        if(dim(u2)[3] != length(unique(u2_moves[, 1]))) stop("u2 does not have entries for all LADS")
        if(!identical(dim(u1)[1:2], dim(u2)[1:2])) stop("u1 and u2 have different numbers of ages / classes")
        
        ## create dummy cohorts for players and bind with workers
        temp <- array(0, c(dim(u1)[1], dim(u1)[2], nrow(u2_moves)))
        u1 <- abind(u1, temp, along = 3)
        u1_moves <- rbind(cbind(u1_moves, rep(NA, nrow(u1_moves))), u2_moves)
        inds <- order(u1_moves[, 1], -u1_moves[, 3])
        u1_moves <- u1_moves[inds, ]
        u1 <- u1[, , inds, drop = FALSE]
        
        ## create list to pass to R
        u1_list <- list(NULL)
        u2_list <- list(NULL)
        for(i in 1:npart) {
            u1_list[[i]] <- u1
            u2_list[[i]] <- u2
        }
        
        ## convert play probs to conditional probs
        playprobs <- tapply(u1_moves[, 3], u1_moves[, 1], function(y) {
            x <- y[!is.na(y)]
            tp <- x[1]
            for(i in 2:length(x)) {
                x[i:length(x)] <- x[i:length(x)] / (1 - tp)
                tp <- x[i]
                ## deal with rounding errors
                if(sum(x[i:length(x)]) != 1) {
#                    print(sum(x[i:length(x)]))
                    x[i:length(x)] <- x[i:length(x)] / sum(x[i:length(x)])
                }
            }
            x[length(x)] <- 1
            y[!is.na(y)] <- x
            y
        })
        playprobs <- do.call("c", playprobs)
        
        ## generate number of cohorts
        ncohorts1 <- tapply(u1_moves[, 1], u1_moves[, 1], length)
        ncohorts1 <- c(0, cumsum(ncohorts1))
        names(ncohorts1) <- NULL
        
        ## generate number of play cohorts
        ncohorts2 <- tapply(u1_moves[, 3], u1_moves[, 1], function(x) sum(!is.na(x)))
        names(ncohorts2) <- NULL
        
        ## remove extraneous column
        u1_moves <- u1_moves[, -3]   
    } else {
        ## NOTE TO REALLY ADD MORE CHECKS IN FOR FUTURE DEVELOPMENT
        if(!all(c("u1", "u2", "u1_moves", "playprobs", "ncohorts1", "ncohorts2") %in% names(u))) {
            stop("'u' must have elements: 'u1', 'u2', 'u1_moves', 'playprobs', 'ncohorts1', 'ncohorts2'")
        }
        u1_moves <- u$u1_moves
        u1_list <- u$u1
        u2_list <- u$u2
        playprobs <- u$playprobs
        ncohorts1 <- u$ncohorts1
        ncohorts2 <- u$ncohorts2
        if(length(u1) != length(u2) | length(u1) != npart) stop("'u1'/'u2' must have 'npart' elements")
    }
    
    ## set up output folder
    if(writeExt) {
        if(!dir.exists(outputName) & is.list(u$u1)) stop("How are you doing continuation without output folder existing?")
        if(is.list(u$u1)) outputName <- paste0(outputName, "_cont")
        if(dir.exists(outputName)) {
            system(paste0("rm -rf ", outputName))
        }
        dir.create(outputName, recursive = TRUE)
    }
    
    ## set up auxiliary parameters guiding parallelisation
    if(parEnsemble) {
        ncoresEns <- ncores
        ncores <- 1
    } else {
        ncoresEns <- 1
    }
    
    ## run particle filter for each set of inputs
    runs <- mclapply(1:nrow(pars), function(k, pars, C1, C2, lockdown_day, u1_moves, ncohorts1, u1, u2, playprobs, ncohorts2, npart, niter, tstart, tstop, cumDI_age_lad, cumDH_age_lad, H_age_lad, lookup, a1, a2, b1, b2, a_dis, b_dis, saveAll, writeExt, snapshot, outputName, PF, ncores) {
    
        ## set up lookups and numbers of regions
        lookup <- as.matrix(lookup)
        ndeathlads <- max(lookup[, 2], na.rm = TRUE)
        lookup[is.na(lookup)] <- -1
        
        ## extract number of stages, age classes and lads
        nclasses <- dim(u1[[1]])[1]
        nages <- dim(u1[[1]])[2]
        nlads <- max(u1_moves[, 1])
        
        if(PF == 1) {
            ## reformat observations
            temp <- mutate(cumDI_age_lad, across(!t, ~. - lag(., default = 0))) %>%
                select(!t) %>%
                as.matrix()
            ## set missing values to be negative for Rcpp code
            temp[is.na(as.matrix(select(cumDI_age_lad, !t)))] <- -1
            DIinc_age_lad <- aperm(array(temp, c(nrow(temp), ndeathlads, nages)), c(1, 3, 2))
            
            temp <- mutate(cumDH_age_lad, across(!t, ~. - lag(., default = 0))) %>%
                select(!t) %>%
                as.matrix()
            ## set missing values to be negative for Rcpp code
            temp[is.na(as.matrix(select(cumDH_age_lad, !t)))] <- -1
            DHinc_age_lad <- aperm(array(temp, c(nrow(temp), ndeathlads, nages)), c(1, 3, 2))
            
            temp <- select(H_age_lad, !t) %>%
                as.matrix()
            ## set missing values to be negative for Rcpp code
            temp[is.na(temp)] <- -1
            H_age_lad <- aperm(array(temp, c(nrow(temp), ndeathlads, nages)), c(1, 3, 2))
            
            rm(temp)
        } else {
            ## set dummies if required
            DIinc_age_lad <- array(0, c(1, 1, 1))
            DHinc_age_lad <- array(0, c(1, 1, 1))
            H_age_lad <- array(0, c(1, 1, 1))
        }

        ## check parameters are in correct order
        cnames <- c("nu", "nuA", 
            paste0(".", paste0(
                rep(c("pE", "pEP", "pA", "pP", "pI1", "pI1H", "pI1D", "pI2", "pH", "pHD"), each = nages),
                "_",
                rep(1:nages, times = 10)
            )),
        "beta_scale", "p_move", "MD_scale", "MD_time_EE", "MD_time_EM", "MD_time_L", "MD_time_NE", "MD_time_NW", "MD_time_SE", "MD_time_SW", "MD_time_WM", "MD_time_YH")
        stopifnot(identical(colnames(pars), cnames))
        pars <- select(pars, !c(MD_scale, starts_with("MD_time")))
        
        ## set pars
        pars <- unlist(pars[k, ])
        
        ## set up a_dis and b_dis arrays
        if(!is.array(a_dis)) {
            stopifnot(is.numeric(a_dis) & length(a_dis) == 1)
            a_dis <- array(rep(a_dis, nages * nlads * (tstop - tstart + 1)), c(tstop - tstart + 1, nages, nlads))
        } else {
            stopifnot(is.numeric(a_dis) & dim(a_dis)[1] == (tstop - tstart + 1) & dim(a_dis)[2] == nages & dim(a_dis)[3] == nlads)
        }
        if(!is.array(b_dis)) {
            stopifnot(is.numeric(b_dis) & length(b_dis) == 1)
            b_dis <- array(rep(b_dis, nages * nlads* (tstop - tstart + 1)), c(tstop - tstart + 1, nages, nlads))
        } else {
            stopifnot(is.numeric(b_dis) & dim(b_dis)[1] == (tstop - tstart + 1) & dim(b_dis)[2] == nages & dim(b_dis)[3] == nlads)
        }
    
        ## do garbage collection (seems to solve allocation issue)
        gc()
        
        ## if just simulating
        if(PF == 0) {
            ## check return outputs
            if(saveAll == 0) stop("Must set 'saveAll' to something if not running a PF")
            
            ## run particle filter
            particles <- BPF_cpp(pars, C1, C2, lockdown_day, DIinc_age_lad, DHinc_age_lad, 
                H_age_lad, lookup, nclasses, nages, nlads, ndeathlads, 
                u1_moves, ncohorts1, u1, u2, playprobs, 
                ncohorts2, tstart, tstop, npart, niter, a1, a2, b1, b2, a_dis, b_dis,
                saveAll, writeExt, snapshot, outputName, PF, ncores)
            return(particles)
        }
        
        ## run particle filter
        particles <- BPF_cpp(pars, C1, C2, lockdown_day, DIinc_age_lad, DHinc_age_lad, 
            H_age_lad, lookup, nclasses, nages, nlads, ndeathlads, 
            u1_moves, ncohorts1, u1, u2, playprobs, 
            ncohorts2, tstart, tstop, npart, niter, a1, a2, b1, b2, a_dis, b_dis,
            saveAll, writeExt, snapshot, outputName, PF, ncores)
        
        if(saveAll != 0 & writeExt == 0) {
            return(list(ll = particles$ll, particles = particles$particles))
        } else {
            return(list(ll = particles$ll))
        }
    }, pars = pars, C1 = C1, C2 = C2, lockdown_day = lockdown_day, u1_moves = u1_moves, ncohorts1 = ncohorts1, u1 = u1_list, u2 = u2_list, playprobs = playprobs, ncohorts2 = ncohorts2, npart = npart, niter = niter, tstart = tstart, tstop = tstop, cumDI_age_lad = cumDI_age_lad, cumDH_age_lad = cumDH_age_lad, H_age_l-ad = H_age_lad, lookup = lookup, a1 = a1, a2 = a2, b1 = b1, b2 = b2, a_dis = a_dis, b_dis = b_dis, saveAll = saveAllint, writeExt = writeExtint, snapshot = snapshotint, outputName = outputName, PF = PFint, ncores = ncores, mc.cores = ncoresEns)
    if(!is.na(saveAll)) {
        ndays <- tstop - tstart
        if(!writeExt) {
            if(PF) {
                ll <- map(runs, "ll")
                runs <- map(function(runs, ndays, npart) {
                        age_lads <- map(runs, "age_lads")
                        xage_lads <- list()
                        for(i in 1:(ndays + 1)) {
                            xage_lads[[i]] <- age_lads[(i - 1) * npart + 1:npart]
                        }
                        list(age_lads = xage_lads)
                    }, ndays = ndays, npart = npart)
                ll <- do.call("c", ll)
                return(list(ll = ll, particles = runs))
            } else {
                runs <- map(runs, "particles") %>%
                    map(function(runs, ndays, npart, saveAll) {
                        age_lads <- map(runs, "age_lads")
                        xage_lads <- list()
                        if(saveAll) {
                            full <- map(runs, "full")
                            xfull <- list()
                        }
                        for(i in 1:(ndays + 1)) {
                            xage_lads[[i]] <- age_lads[(i - 1) * npart + 1:npart]
                            if(saveAll) {
                                xfull[[i]] <- full[(i - 1) * npart + 1:npart]
                            }
                        }
                        if(saveAll) {
                            return(list(full = xfull, age_lads = xage_lads))
                        } else {
                            return(list(lads = xage_lads))
                        }
                    }, ndays = ndays, npart = npart, saveAll = saveAll)
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
        ll <- do.call("c", ll)
        return(list(ll = ll))
    }
}
