#' Limited information likelihood
#' @param d object of class \code{RDData}
#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
IVregIL.fit <- function(d) {

    be <- (d$T[1, 2]- d$ei[1] * d$S[1, 2]) / (d$T[2, 2]-d$ei[1]*d$S[2, 2])

    ## log(G_k(t))
    lG <- function(t, k)
        log(besselI(t, k/2-1, expon.scaled=TRUE))+t-(k/2-1)*log(t/2)

    ## minus log-likelihood
    logl <- function(be, lam, Om) {
        Qt <- drop(crossprod(solve(Om, c(be, 1)), d$T %*% solve(Om, c(be, 1))) /
                   aoa(be, Om))
        ((d$n-d$l)*log(det(Om)) + sum(diag(solve(Om, d$nu*d$S+d$n*d$T)))
            + d$n*lam - 2*lG(d$n*sqrt(max(lam*Qt, 1e-10)), d$k)) / 2
    }

    Q.lam <- function(la)
        d$ei[2]* (d$n-d$l+d$n*la) / (d$nu+d$n*d$ei[2])

    ff <- function(lam)
        besselI(d$n*sqrt(lam*Q.lam(lam)), d$k/2, expon.scaled=TRUE) /
            besselI(d$n*sqrt(lam*Q.lam(lam)), d$k/2-1, expon.scaled=TRUE) -
            sqrt(lam/Q.lam(lam))

    if (ff(1e-8)<0 | is.nan(ff(1e-8)))
        return(list(beta=be, se=Inf, lam=1e-6, Om=d$S,
                    lik=-logl(be, 1e-6, d$S)))

    lam <- stats::uniroot(ff, c(1e-8, 3*d$ei[2]+1), tol=1e-10)$root
    Om <- (d$nu*d$S+d$n*d$T - d$n*lam*d$ei[2]/Q.lam(lam) *
           (c(be, 1) %o% c(be, 1)) / (aoa(be, d$S))) / (d$n-d$l)

    ## Hessian
    fh <- function(t) logl(t[1], t[2], cbind(c(t[3], t[4]), c(t[4], t[5])))
    Hes <- numDeriv::hessian(fh, c(be, lam, Om[1, 1],
                                             Om[1, 2], Om[2, 2]))
    se <- if (sum(is.nan(Hes))>0) {
              Inf
          } else if (Matrix::rankMatrix(Hes)<5) {
              Inf
          } else {
              sqrt(solve(Hes)[1, 1])
          }

    list(beta=be, se=se, lam=lam, Om=Om, lik=-logl(be, lam, Om))
}
