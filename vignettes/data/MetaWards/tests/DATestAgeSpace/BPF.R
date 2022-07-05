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
## u2_moves: matrix of commuter data, with columns:
##            home LAD, play LAD, prob of moving
## u2: 3D array of player data, with dimensions:
##            nclasses x nages x nlads
## ndays: the number of days to fit to
## npart: the number of particles
## niter: number of Metropolis-Hastings steps to counter particle impoverishment
## obsScale: scaling parameter for Skellam observation process (see code)
## a1, a2, b: parameters for Skellam observation process
## saveAll: a logical specifying whether to return all states (if FALSE then returns just observed states))
## writeExt: a logical denoting whether to save particles externally or not
## PF:      a logical denoting whether to run a particle filter, or just simulate from the model
## ncores:  the number of cores for OpenMP parallelisation (if NA then defaults to all available cores)
## parEnsemble: decides whether to parallelise across or within ensemble

BPF <- function(pars, C, data, u1_moves, u1, u2_moves, u2, ndays, npart = 10, niter = 0, a1 = 0.01, a2 = 0.2, b = 0.1, 
        a_dis = 0.05, b_dis = 0.05, saveAll = NA, writeExt = FALSE, PF = TRUE, ncores = NA, parEnsemble = FALSE) {
               
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
    if(dim(u2)[3] != length(unique(u2_moves[, 1]))) stop("u2 does not have entries for all LADS")
    if(!identical(dim(u1)[1:2], dim(u2)[1:2])) stop("u1 and u2 have different numbers of ages / classes")
    
    ## check probabilities
    if(!all((tapply(u2_moves[, 3], u2_moves[, 1], sum) - 1) < 1e-15)) stop("u2 probs don't sum to one")
    
    ## create dummy cohorts for players and bind with workers
    temp <- array(0, c(dim(u1)[1], dim(u1)[2], nrow(u2_moves)))
    u1 <- abind(u1, temp, along = 3)
    u1_moves <- rbind(cbind(u1_moves, rep(NA, nrow(u1_moves))), u2_moves)
    inds <- order(u1_moves[, 1], -u1_moves[, 3])
    u1_moves <- u1_moves[inds, ]
    u1 <- u1[, , inds, drop = FALSE]
    
    ## convert play probs to conditional probs
    playprobs <- tapply(u1_moves[, 3], u1_moves[, 1], function(y) {
        x <- y[!is.na(y)]
        tp <- x[1]
        for(i in 2:length(x)) {
            x[i:length(x)] <- x[i:length(x)] / (1 - tp)
            tp <- x[i]
            ## deal with rounding errors
            if(sum(x[i:length(x)]) != 1) {
#                print(sum(x[i:length(x)]))
                x[i:length(x)] <- x[i:length(x)] / sum(x[i:length(x)])
            }
        }
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
    
    ## set up output folder
    if(writeExt) {
        print("Reminder to write code to pass save folder out")
        if(dir.exists("saveOut")) system("rm -rf saveOut")
        dir.create("saveOut")
    }
    
    ## set up auxiliary parameters guiding parallelisation
    if(parEnsemble) {
        ncoresEns <- ncores
        ncores <- 1
    } else {
        ncoresEns <- 1
    }
    
    ## run particle filter for each set of inputs
    runs <- mclapply(1:nrow(pars), function(k, pars, C, u1_moves, ncohorts1, u1, u2, playprobs, ncohorts2, npart, niter, ndays, data, a1, a2, b, a_dis, b_dis, saveAll, writeExt, PF, ncores) {
        
        if(PF == 1) {
            ## extract observations
            data <- select(data, t, (starts_with("DI") | starts_with("DH")) & contains("obs")) %>%
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
        
        ## extract number of stages, age classes and lads
        nclasses <- dim(u1)[1]
        nages <- dim(u1)[2]
        nlads <- max(u1_moves[, 1])
    
        ## do garbage collection (seems to solve allocation issue)
        gc()
        
        ## if just simulating
        if(PF == 0) {
            ## check return outputs
            if(saveAll == 0) stop("Must set 'saveAll' to something if not running a PF")
            
            ## run particle filter
            particles <- BPF_cpp(pars, C, data, nclasses, nages, nlads, u1_moves, ncohorts1, u1,
                u2, playprobs, ncohorts2, ndays, npart, niter, a1, a2, b, a_dis, b_dis, saveAll, 
                writeExt, PF, ncores)
            return(particles)
        }
        
        ## run particle filter
        particles <- BPF_cpp(pars, C, data, nclasses, nages, nlads, u1_moves, ncohorts1, u1, 
            u2, playprobs, ncohorts2, ndays, npart, niter, a1, a2, b, a_dis, b_dis, saveAll, 
            writeExt, PF, ncores)
        
        if(saveAll != 0 & writeExt == 0) {
            return(list(ll = particles$ll, particles = particles$particles))
        } else {
            return(list(ll = particles$ll))
        }
    }, pars = pars, C = C, u1_moves = u1_moves, ncohorts1 = ncohorts1, u1 = u1, u2 = u2, playprobs = playprobs, ncohorts2 = ncohorts2, npart = npart, niter = niter, ndays = ndays, data = data, a1 = a1, a2 = a2, b = b, a_dis = a_dis, b_dis = b_dis, saveAll = saveAllint, writeExt = writeExtint, PF = PFint, ncores = ncores, mc.cores = ncoresEns)
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
        ll <- do.call("c", ll)
        return(ll)
    }
}
