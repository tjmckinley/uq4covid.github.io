## load libraries
library(skellam)

## truncated Skellam density
dtskellam <- function(x, lambda1, lambda2, LB = -Inf, UB = Inf, log = FALSE) {
    
    ## checks on inputs
    if(missing(x) | missing(lambda1) | missing(lambda2)) stop("'x', 'lambda1' and 'lambda2' must be specified")
    if(any(!is.numeric(c(x, lambda1, lambda2, LB, UB)))) stop("Inputs must be numeric")
    if(!is.logical(log)) stop("'log' must be logical")
    if(length(log) > 1) {
        cat("Only first element of 'log' used")
        log <- log[1]
    }
    if(length(lambda1) != length(x) & length(lambda1) > 1) stop("'lambda1' must be of length 1 or length(x)")
    if(length(lambda2) != length(x) & length(lambda2) > 1) stop("'lambda2' must be of length 1 or length(x)")
    if(length(LB) != length(x) & length(LB) > 1) stop("'LB' must be of length 1 or length(x)")
    if(length(UB) != length(x) & length(UB) > 1) stop("'UB' must be of length 1 or length(x)")
    
    ## expand inputs
    if(length(lambda1) != length(x)) lambda1 <- rep(lambda1, length(x))
    if(length(lambda2) != length(x)) lambda2 <- rep(lambda2, length(x))
    if(length(LB) != length(x)) LB <- rep(LB, length(x))
    if(length(UB) != length(x)) UB <- rep(UB, length(x))
    
    ## more checks
    if(any(LB > UB)) stop("'LB' can't be > 'UB'")
    if(any(round(LB) != LB & is.finite(LB)) | any(round(UB) != UB & is.finite(UB))) stop("'LB' and 'UB' must be integers")
    if(any(apply(cbind(LB, UB), 1, function(x) identical(x[1], x[2])) & !is.finite(LB))) stop("'LB' and 'UB' misspecified")
    if(any(c(lambda1, lambda2) <= 0)) stop("'lambda1' and 'lambda2' must be > 0")
    
    ## generate non-truncated densities
    ldens <- dskellam(x, lambda1, lambda2, log = TRUE)
    ldens[ldens > 0] <- -Inf
    
    ## adjust components for truncation as required
    UBind <- which(!is.finite(LB) & is.finite(UB))
    LBind <- which(is.finite(LB) & !is.finite(UB))
    bothind <- which(is.finite(LB) & is.finite(UB))
    if(length(UBind) > 0) {
        ldens[UBind] <- ldens[UBind] - pskellam(UB[UBind], lambda1[UBind], lambda2[UBind], lower.tail = TRUE, log.p = TRUE)
        ldens[UBind] <- ifelse(x[UBind] > UB[UBind], -Inf, ldens[UBind])
    }
    if(length(LBind) > 0) {
        ldens[LBind] <- ldens[LBind] - pskellam(LB[LBind] - 1, lambda1[LBind], lambda2[LBind], lower.tail = FALSE, log.p = TRUE)
        ldens[LBind] <- ifelse(x[LBind] < LB[LBind], -Inf, ldens[LBind])
    }
    if(length(bothind) > 0) {
        LBnorm <- pskellam(LB[bothind] - 1, lambda1[bothind], lambda2[bothind], lower.tail = TRUE, log.p = TRUE)
        UBnorm <- pskellam(UB[bothind], lambda1[bothind], lambda2[bothind], lower.tail = TRUE, log.p = TRUE)
        if(any(LBnorm > 0 | UBnorm > 0)) browser()
        norm <- UBnorm + log(1 - exp(LBnorm - UBnorm))
        ldens[bothind] <- ldens[bothind] - norm
        ldens[bothind] <- ifelse(x[bothind] < LB[bothind] | x[bothind] > UB[bothind], -Inf, ldens[bothind])
    }
    if(log) {
        return(ldens)
    } else {
        return(exp(ldens))
    }
}

## truncated Skellam sampling
rtskellam <- function(n, lambda1, lambda2, LB = -Inf, UB = Inf, ntries = 1000) {
    
    ## checks on inputs
    if(missing(n) | missing(lambda1) | missing(lambda2)) stop("'n', 'lambda1' and 'lambda2' must be specified")
    if(any(!is.numeric(c(n, lambda1, lambda2, LB, UB, ntries)))) stop("Inputs must be numeric")
    if(length(n) != 1) stop("'n' must be length 1")
    if(round(n) != n | n < 1) stop("'n' must be a positive integer")
    if(length(lambda1) != n & length(lambda1) > 1) stop("'lambda1' must be of length 1 or n")
    if(length(lambda2) != n & length(lambda2) > 1) stop("'lambda2' must be of length 1 or n")
    if(length(LB) != n & length(LB) > 1) stop("'LB' must be of length 1 or n")
    if(length(UB) != n & length(UB) > 1) stop("'UB' must be of length 1 or n")
    
    ## expand inputs
    if(length(lambda1) != n) lambda1 <- rep(lambda1, n)
    if(length(lambda2) != n) lambda2 <- rep(lambda2, n)
    if(length(LB) != n) LB <- rep(LB, n)
    if(length(UB) != n) UB <- rep(UB, n)
    
    ## more checks
    if(any(LB > UB)) stop("'LB' can't be > 'UB'")
    if(any(round(LB) != LB & is.finite(LB)) | any(round(UB) != UB & is.finite(UB))) stop("'LB' and 'UB' must be integers")
    if(any(apply(cbind(LB, UB), 1, function(x) identical(x[1], x[2])) & !is.finite(LB))) stop("'LB' and 'UB' misspecified")
    if(any(c(lambda1, lambda2) <= 0)) stop("'lambda1' and 'lambda2' must be > 0")
    
    if(all(LB == -Inf) & all(UB == Inf)) {
        return(rskellam(n, lambda1, lambda2))
    }
    x <- rskellam(n, lambda1, lambda2)
    missind <- which(x < LB | x > UB)
    k <- 1
    while(length(missind) > 0 & k < ntries) {
        x[missind] <- rskellam(length(missind), lambda1[missind], lambda2[missind])
        missind <- which(x < LB | x > UB)
        k <- k + 1
    }
    if(length(missind) > 0) {
        ## if brute force doesn't work, then try inverse transform sampling
        u <- runif(n, 0, 1)
        if(any(is.finite(LB) & !is.finite(UB))) {
            stop("Need to check this")
            browser()
            # ## normalising constant
            # norm <- pskellam(LB, lambda1, lambda2, lower.tail = FALSE, log.p = TRUE)
            # UB <- 1000
            # x <- cumsum(exp(dskellam(LB:UB, lambda1, lambda2, log = TRUE) - norm))
            # ind <- which(x >= u)
            # while(length(ind) == 0) {
            #     LB <- UB
            #     UB <- UB + 1000
            #     x <- x[length(x)] + cumsum(exp(dskellam(LB:UB, lambda1, lambda2, log = TRUE) - norm))
            #     ind <- which(x >= u)
            # }
            # x <- c(LB:UB)[ind[1]]
        } else {
            if(any(!is.finite(LB) & is.finite(UB))) {
                stop("Need to check this")
                browser()
                # ## normalising constant
                # norm <- pskellam(UB, lambda1, lambda2, lower.tail = TRUE, log.p = TRUE)
                # LB <- -1000
                # x <- cumsum(rev(exp(dskellam(LB:UB, lambda1, lambda2, log=  TRUE) - norm)))
                # ind <- which(x >= u)
                # while(length(ind) == 0) {
                #     UB <- LB
                #     LB <- LB - 1000
                #     x <- x[length(x)] + cumsum(rev(exp(dskellam(LB:UB, lambda1, lambda2, log.p = TRUE) - norm)))
                #     ind <- which(x >= u)
                # }
                # x <- rev(LB:UB)[ind[1]]
            } else {
                ## normalising constant
                x[missind] <- sapply(missind, function(i, u, UB, LB, lambda1, lambda2) {
                    x <- dskellam(LB[i]:UB[i], lambda1[i], lambda2[i], log = TRUE)
                    maxn <- max(x)
                    norm <- maxn + log(sum(exp(x - maxn)))
                    x <- cumsum(exp(x - norm))
                    ind <- which(x >= u[i])
                    c(LB[i]:UB[i])[ind[1]]
                }, u = u, LB = LB, UB = UB, lambda1 = lambda1, lambda2 = lambda2)
            }
        }
#        print(paste("lambda1 =", lambda1, "lambda2 =", lambda2, "LB =", LB, "UB =", UB))
#        stop(paste0("Can't sample from truncated Skellam after ", ntries, " tries"))
    }
    x
}
