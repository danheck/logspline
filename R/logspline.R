# Copyright (c) 1993 by Charles Kooperberg.                                    *
# All rights reserved. This function is part of the Logspline Density          *
# Estimation Program. It was written by Charles Kooperberg at the University   *
# of Washington and the University of California at Berkeley between 1988 and  *
# 1991. It is described in the paper: Kooperberg, Charles and Stone, Charles J.*
# `Logspline density estimation for censored data', Journal of Computational   *
# and Graphical Statistics, 1 (1992) 301-328.                                  *
# You are free to use this program, for non-commercial purposes only,          *
# under two conditions:                                                        *
# (1) This note is not to be removed                                           *
# (2) Publications using logspline computations should reference the           *
#  publication  mentioned above.                                               *
# For questions, please email clk@stat.washington.edu                          *
#                       Charles Kooperberg, April 21, 1993                     *
plogspline <- function(q, fit)
{
    sq <- rank(q)
    q <- sort(q)
    z <- .C("pqlsd",
            as.double(fit$coef),
            as.double(fit$knots),
            as.double(fit$bound),
            as.integer(1),
            pp = as.double(q),
            as.double(q),
            as.integer(length(fit$knots)),
            as.integer(length(q)),
            PACKAGE="logspline"
            )
    zz <- z$pp[sq]
    if(fit$bound[1] > 0)
        zz[q<fit$bound[2]] <- 0
    if(fit$bound[3] > 0)
        zz[q>fit$bound[4]] <- 1
    zz
}

qlogspline <- function(p, fit)
{
    sp <- rank(p)
    p <- sort(p)
    z <- .C("pqlsd",
            as.double(fit$coef),
            as.double(fit$knots),
            as.double(fit$bound),
            as.integer(0),
            as.double(p),
            qq = as.double(p),
            as.integer(length(fit$knots)),
            as.integer(length(p)),
            PACKAGE="logspline")
    zz <- z$qq[sp]
    zz[p<0] <- NA
    zz[p>1] <- NA
    zz
}

rlogspline <- function(n, fit)
{
    pp <- runif(n)
    qlogspline(pp, fit)
}

dlogspline <- function(x, fit)
{
    y <- fit$coef[1] + x * fit$coef[2]
    for(i in 1:length(fit$knots)) {
        if(fit$coef[i+2] != 0)
            y <- y + fit$coef[i+2] * ((abs(x - fit$knots[i]) +
                                x - fit$knots[i])/2)^3
    }
    y <- exp(y)
    if(fit$bound[1] > 0)
            y[x < fit$bound[2]] <- 0
    if(fit$bound[3] > 0)
            y[x > fit$bound[4]] <- 0
    y
}

logspline.plot <- function(fit, n = 100, what = "d", xlim, xlab, ylab, type, ...)
{
        if(missing(xlim)) {
                u1 <- qlogspline(0.01, fit)
                u2 <- qlogspline(0.99, fit)
                u3 <- 1.1 * u1 - 0.1 * u2
                u4 <- 1.1 * u2 - 0.1 * u1
        }
        else {
                u3 <- xlim[1]
                u4 <- xlim[2]
        }
        xx <- (0:(n - 1))/(n - 1) * (u4 - u3) + u3
        if(what == "d" || what == "D")
                yy <- dlogspline(xx, fit)
        if(what == "f" || what == "F" || what == "p" || what == "P")
                yy <- plogspline(xx, fit)
        if(what == "s" || what == "S")
                yy <- 1 - plogspline(xx, fit)
        if(what == "h" || what == "H")
                yy <- dlogspline(xx, fit)/(1 - plogspline(xx, fit))
        if(missing(xlab))
                xlab <- ""
        if(missing(ylab))
                ylab <- ""
        if(missing(type))
                type <- "l"
        plot(xx, yy, xlab = xlab, ylab = ylab, type = type, ...)
}

logspline.summary <- function(fit)
{
    if(fit$delete==FALSE)
        stop(paste("logspline.summary can only provide",
       "information if delete in logspline.fit is T"))
    ul <- fit$penalty
    um <- fit$sample
    ll <- fit$logl
    kk <- (1:length(ll))
    kk <- kk[ll != 0] + 2
    ll <- ll[ll != 0]
    error<-FALSE
    rr <- ll[1:(length(ll)-1)]-ll[2:length(ll)]
    if(max(rr)>0)error<-TRUE
    bb <- -2 * ll + ul * kk
    cc1 <- bb
    cc2 <- bb
    cc2[1] <- 5/0
    cc1[length(bb)] <- 0
    if(length(bb) > 1) {
        for(i in 1:(length(bb) - 1)) {
            cc1[i] <- max((ll[(i + 1):(length(bb))] - ll[i])/(
                    kk[(i + 1):(length(bb))] - kk[i]))
            cc2[i + 1] <- min((ll[1:i] - ll[i + 1])/(kk[1:i] - kk[i + 1]))
        }
    }
    c3 <- cc2 - cc1
    cc1[c3 < 0] <- NA
    cc2[c3 < 0] <- NA
    uu <- cbind(kk, ll, bb, 2 * cc1, 2 * cc2)
    ww <- rep("", length(bb))
    if(error){
    cat("Warning - imprecision in loglikelihood (possibly due to heavy tails)\n")
    cat("the output of logspline.summary might not be correct\n")
    }
    dimnames(uu) <- list(ww, c("knots", "loglik", "AIC", "minimum penalty",
        "maximum penalty"))
    print(round(uu, 2))
    cat(paste("the present optimal number of knots is ", kk[bb== min(bb)],"\n"))
    if(ul == log(um))
        cat(paste("penalty(AIC) was the default: BIC=log(samplesize): log(",
                um, ")=", round(ul, 2),"\n"))
    else
        cat(paste("penalty(AIC) was ", round(ul, 2),", the default (BIC) ",
                "would have been", round(log(um), 2),"\n"))
    if(min(kk) > 3 && fit$delete==TRUE){
        cat(paste( "models with fewer than", kk[1],"knots ",
                  "can be fitted, but they are not optimal for\n"))
        cat(paste("the present choice of penalty - choose penalty in",
                  "logspline.fit larger\nto see these fits\n"))
    }
    if(min(kk) > 3 && fit$delete==3)
        cat(paste("models with fewer than", kk[1],"knots ",
                    "were not fitted because of convergence problems\n"))

    invisible()
}

logspline.fit <- function(uncensored, right, left, interval, lbound, ubound,
        nknots, knots, penalty, delete = TRUE)
{
    nsample <- rep(0, 6)
    # interval is the nterval censored data - a matrix with two columns
    if(!missing(interval)) {
        if(length(interval[1,  ]) != 2)
            stop("interval must have two columns")
        if(min(abs(interval[, 1] - interval[, 2])) < 0) stop(
                   "not all lower bounds smaller than upper bounds")
        nsample[3] <- length(interval)/2
        nsample[1] <- length(interval)/2
        # grouping boundaries can not be beyond the boundaries of the density
        if(!missing(lbound))
            interval[interval[, 1] < lbound, 1] <- lbound
        if(!missing(ubound))
            interval[interval[, 2] > ubound, 2] <- ubound
        sample <- as.vector(t(interval))
        ror <- order(interval[,1],interval[,2])
        if(nsample[3]>1){
      ro1 <- interval[ror[(1:(nsample[3]-1))],1]==interval[ror[2:nsample[3]],1]
      ro2 <- interval[ror[(1:(nsample[3]-1))],2]==interval[ror[2:nsample[3]],2]
            nsample[6] <- nsample[3]-sum(ro1+ro2==2)
        }
        else nsample[6] <- 1
    }
# uncensored is the uncensored data
    if(!missing(uncensored)) {
        uncensored2 <- uncensored[!is.na(uncensored)]
        u2 <- length(uncensored) - length(uncensored2)
        if(u2 > 0)
            print(paste("***", u2, " NAs ignored in uncensored"))
        uncensored <- uncensored2
        if(nsample[1] > 0)
            sample <- c(uncensored, sample)
        if(nsample[1] == 0)
            sample <- uncensored
        nsample[1] <- length(uncensored) + nsample[1]
        nsample[2] <- length(uncensored)
        uncensored <- sort(uncensored)
        if(nsample[2]>1)
            nsample[6] <- sum(uncensored[2:nsample[2]] !=
                uncensored[1:(nsample[2]-1)]) + 1 + nsample[6]
        else
            nsample[6] <- nsample[6]+1
    }
# we can not run on only right or left censored data
        if(nsample[1] == 0) stop("you either need uncensored or interval censored data")
        # right is the right censored data
        if(!missing(right)) {
                if(nsample[1] > 0)
                        sample <- c(sample, right)
                if(nsample[1] == 0)
                        sample <- right
                nsample[1] <- length(right) + nsample[1]
                nsample[4] <- length(right)
                right <- sort(right)
                if(nsample[4]>1){
                nsample[6] <- sum(right[2:nsample[4]]!=right[1:(nsample[4]-1)])+
                          1 + nsample[6]
                }
                else nsample[6] <- nsample[6]+1
        }
# left is the left censored data
        if(!missing(left)) {
                if(nsample[1] > 0)
                        sample <- c(sample, left)
                if(nsample[1] == 0)
                        sample <- left
                nsample[1] <- length(left) + nsample[1]
                nsample[5] <- length(left)
                left <- sort(left)
                if(nsample[5]>1){
                nsample[6] <- sum(left[2:nsample[5]]!=left[1:(nsample[5]-1)])+
                          1 + nsample[6]
                }
                else nsample[6] <- nsample[6]+1
        }
# the default for penalty is bic: log(length(sample))
        if(missing(penalty)) penalty <- log(nsample[1])
        n1 <- 4 * nsample[1]^0.2 + 1
        if(!missing(nknots))
                n1 <- nknots + 1
        if(!missing(knots)) n1 <- length(knots) + 1      # user provides knots
        if(!missing(knots)) {
                nknots <- length(knots)
                knots <- sort(knots)
                iautoknot <- 0
                if(knots[1] > min(sample))
                        stop("first knot must be smaller than smallest sample")

                if(knots[nknots] < max(sample))
                        stop("last knot should be larger than largest sample")

        }
        else {
                if(missing(nknots))
                        nknots <- 0
                knots <- vector(mode = "double", length = max(nknots, 50))
                iautoknot <- 1
        }
        xbound <- c(1, 0, 0, 0, 0)
        if(!missing(lbound)) {
                xbound[2] <- 1
                xbound[3] <- lbound
                if(lbound > min(sample))
                        stop("lbound should be smaller than smallest sample")
        }
        if(!missing(ubound)) {
                xbound[4] <- 1
                xbound[5] <- ubound
                if(ubound < max(sample))
                        stop("ubound should be larger than largest sample")
        }
# SorC will carry the error messages - in code form
        SorC <- vector(mode = "integer", length = 35)
        SorC[1] <- 1    # the actual function call
        nsample[6] <- nsample[6]-1
        z <- .C("logcensor",
                as.integer(delete),
                as.integer(iautoknot),
                as.double(sample),
                as.integer(nsample),
                bd = as.double(xbound),
                SorC = as.integer(SorC),
                nk = as.integer(nknots),
                kt = as.double(knots),
                cf = as.double(c(knots, 0, 0)),
                as.double(penalty),
                as.double(sample),
                as.double(sample),
                logl = as.double(rep(0, n1 + 1)),
                PACKAGE="logspline")
        bound <- c(z$bd[2], z$bd[3], z$bd[4], z$bd[5])
        SorC <- z$SorC  # error messages
        if(abs(SorC[1]) > 2) {
                for(i in 3:abs(SorC[1]))
                        cat(paste("===> warning: knot ", SorC[i - 1],
                                " removed - double knot\n"))
                if(SorC[1] < 0)
                        SorC[1] <- -1
                if(SorC[1] == 23)
                        SorC[1] <- -3
        }
        if(abs(SorC[1]) > 3) {
                cat("* several double knots suggests that your data is *\n")
                cat("* strongly rounded, attention might be required   *\n")
                SorC[1] <- 1
        }
        if(SorC[1] == -3)
                stop("* too many double knots")
        if(SorC[1] == -1 && SorC[28] == 0)
                stop("* no convergence")
        if(SorC[28] > 0)
                cat(paste("* convergence problems, smallest number of knots",
                        " tried is ", SorC[28] + 1," *\n"))
        if(SorC[1] == 2)
                stop("* sample is too small")
        if(SorC[1] == -2)
                stop(paste("* too many knots, at most ", SorC[2],
                        "knots possible"))
        if(SorC[22] == 1) {
                cat("possible discontinuity at lower end\n")
                cat(paste("consider rerunning with lbound=", z$kt[1],
			"\n"))

        }
        if(SorC[22] == 3) {
                cat("possible infinite density at lower end\n")
                cat("running program with fewer knots\n")
        }
        if(SorC[21] == 1)
                cat("running with maximum degrees of freedom\n")
        if(SorC[25] >0)
               cat("* problems are possibly due to a very heavy right tail *\n")
        if(SorC[24] >0)
                cat("* problems are possibly due to a very heavy left tail *\n")
        if(SorC[23] == 3) {
                cat("possible infinite density at upper end\n")
                cat("running program with fewer knots\n")
        }
        if(SorC[23] == 1) {
                cat("possible discontinuity at upper end\n")
                cat(paste("consider rerunning with ubound=", z$kt[z$nk],
			"\n"))

        }
        if(delete && SorC[28]>0)delete<-3
        coef <- z$cf[1:(z$nk + 2)]
        uu <- 3:z$nk
        if(delete == FALSE) uu <- 1
        list(coef = coef, knots = z$kt[1:z$nk], bound = bound, logl = z$logl[
                uu], penalty = penalty, sample = nsample[1], delete = delete)
}

.First.lib <- function(lib, pkg)
{
   library.dynam("logspline", pkg, lib)
}
