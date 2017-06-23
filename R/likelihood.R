#' Limited information likelihood
#' @param d object of class \code{RDData}
IVregLI.fit <- function(d) {
    be <- (d$T[1, 2]- d$ei[1] * d$S[1, 2]) / (d$T[2, 2]-d$ei[1]*d$S[2, 2])
    Om <- (d$nu*d$S + d$n*(d$T-d$ei[2]*(c(be, 1) %o% c(be, 1)) /
                           aoa(be, d$S))) / (d$n-d$l)
    lam <- (d$n-d$l) * d$ei[2] / d$nu
    ## Information matrix at ML estimates
    se <- sqrt(aoa(be, d$S)*bob(be, d$S)/ (d$n*lam))
    list(beta=be, se=se, lam=lam, Om=Om)
}

#' Random-effects likelihood
#' @param d object of class \code{RDData}
IVregRE.fit <- function(d) {
    ## Parameter estimates
    be <- (d$T[1, 2]- d$ei[1] * d$S[1, 2]) / (d$T[2, 2]-d$ei[1]*d$S[2, 2])
    lam <- max(d$ei[2] - d$k/d$n, 0)
    Om <- (d$nu*d$S + d$n*(d$T-lam*c(be, 1) %o% c(be, 1) /
                           aoa(be, d$S))) / (d$n-d$l)

    ## SE based on RE likelihood
    Qs <- bob(be, d$T) / bob(be, Om)
    c.re <- lam * Qs / ((d$k/d$n+lam) * (1-d$l/d$n))
    H.re <- bob(be, Om) * (lam+d$k/d$n) /
        (d$n*lam * (Qs*Om[2, 2] - d$T[2, 2] +
                       c.re*Qs/((1-c.re)*aoa(be, Om))))

    list(beta=be, se=sqrt(-H.re), lam=lam, Om=Om)
}


#' Invariant likelihood
#' @param d object of class \code{RDData}
IVregIL.fit <- function(d) {

    ## log(G_k(t))
    lG <- function(t, k)
        log(besselI(t, k/2-1, expon.scaled=TRUE))+t-(k/2-1)*log(t/2)

    ## minus log-likelihood
    logl <- function(be, lam, Om) {
        Qt <- drop(crossprod(solve(Om, c(be, 1)), d$T %*% solve(Om, c(be, 1))) /
                   aoa(be, Om))
        ((d$n-d$l)*log(det(Om)) + sum(diag(solve(Om, d$nu*d$S+d$n*d$T)))
            + d$n*lam - 2*lG(d$n*sqrt(lam*Qt), d$k)) / 2
    }
    ## Reformat inputs
    ff <- function(be, t)
        logl(be, exp(t[1]), cbind(c(exp(t[2]), t[3]), c(t[3],
        (exp(t[4])+t[3]^2)/exp(t[2]))))
    ## ff <- function(be, t)
    ##     logl(be, exp(t[1]), cbind(c(t[2], t[3]), c(t[3], t[4])))

    r1 <- IVregRE.fit(d)
    ## start <- c(log(r1$lam), log(r1$Om[1, 1]), r1$Om[1, 2], log(det(r1$Om)))
    start <- c(log(r1$lam), log(d$S[1, 1]), d$S[1, 2], log(det(d$S)))

    il <- stats::optim(start, function(t) ff(r1$beta, t))
    se.il <- sqrt(solve(numDeriv::hessian(function(t) ff(t[1], t[2:5]),
                                          c(r1$beta, il$par)))[1, 1])
    Om <- cbind(c(exp(il$par[2]), il$par[3]),
                c(il$par[3], (exp(il$par[4])+il$par[3]^2)/exp(il$par[2])))

    list(beta=r1$beta, se=se.il, lam=exp(il$par[1]), Om=Om,
         lik=-logl(r1$beta, exp(il$par[1]), Om))
}


#' Invariant likelihood
#' @param d object of class \code{RDData}
IVregIL2.fit <- function(d) {

    ## log(G_k(t))
    lG <- function(t, k)
        log(besselI(t, k/2-1, expon.scaled=TRUE))+t-(k/2-1)*log(t/2)

    ## minus log-likelihood
    logl <- function(be, lam, Om) {
        Qt <- drop(crossprod(solve(Om, c(be, 1)), d$T %*% solve(Om, c(be, 1))) /
                   aoa(be, Om))
        ((d$n-d$l)*log(det(Om)) + sum(diag(solve(Om, d$nu*d$S+d$n*d$T)))
            + d$n*lam - 2*lG(d$n*sqrt(lam*Qt), d$k)) / 2
    }

    Q.lam <- function(la)
        d$ei[2]* (d$n-d$l+d$n*la) / (d$nu+d$n*d$ei[2])

    ff <- function(lam)
        besselI(d$n*sqrt(lam*Q.lam(lam)), d$k/2, expon.scaled=TRUE) /
            besselI(d$n*sqrt(lam*Q.lam(lam)), d$k/2-1, expon.scaled=TRUE) -
            sqrt(lam/Q.lam(lam))
    lam0 <- max(d$ei[2]-d$k/d$n)
    lam <- stats::uniroot(ff, c(lam0/3+1e-6, 3*lam0+1), tol=1e-10)$root
    be <- (d$T[1, 2]- d$ei[1] * d$S[1, 2]) / (d$T[2, 2]-d$ei[1]*d$S[2, 2])
    Om <- (d$nu*d$S+d$n*d$T - d$n*lam*d$ei[2]/Q.lam(lam) *
           (c(be, 1) %o% c(be, 1)) / (aoa(be, d$S))) / (d$n-d$l)

    ## Hessian
    ff <- function(t) logl(t[1], t[2], cbind(c(t[3], t[4]), c(t[4], t[5])))
    se <- sqrt(solve(numDeriv::hessian(ff, c(be, lam, Om[1, 1],
                                             Om[1, 2], Om[2, 2])))[1, 1])

    list(beta=be, se=se, lam=lam, Om=Om, lik=-logl(be, lam, Om))
}
