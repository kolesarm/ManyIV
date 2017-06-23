library('numDeriv')

## Get features of the design matrices
MDDesignFeatures <- function(Z, W, Mfeat=FALSE) {
    d <- list()
    d$k <- ncol(Z)                        # number of instruments
    d$l <- ncol(W)                        # number of controls
    d$n <- nrow(Z)                        # number of observations
    d$nu <- d$n-d$k-d$l                   # degrees of freedom

    d$qr.W <- qr(W)
    d$Zp <- qr.Q(qr(qr.resid(d$qr.W, as.matrix(Z)))) # as.matrix in case Z is sparse

    ## Compute M matrix
    if (Mfeat) {
        PZp <- tcrossprod(d$Zp)
        d$M <- diag(d$n)-W%*%solve(crossprod(W))%*% t(W)-PZp
        d$mii <- sum(diag(d$M)^2)
        d$m3 <- sum(d$M^3)
        d$m4 <- sum(d$M^4)
    }

    d
}

MDivreg <- function(Y, Z, W, theta=NULL, desFeat=NULL) {
    ## Y is an [n x 2] matrix of LHS variables, Z is an [n x k] matrix of
    ## instruments, and W an [n x l] matrix of exogenous regressors.

    d <- if(!is.null(desFeat)) desFeat else MDDesignFeatures(Z, W)

    Yt <- qr.resid(d$qr.W, Y)            # \tilde{Y}
    Rhat <- d$Zp %*% crossprod(d$Zp, Yt)
    T <- crossprod(Rhat) / d$n
    S <- crossprod(Yt-Rhat) / d$nu
    ei <- sort(eigen(solve(S, T))$values) # [m_min, m_max]

    nl <- d$n - d$l
    ak <- d$k/d$n

    Qs <- function(be, Om) bob(be, T) / bob(be, Om)

    Gamma <- function(be) matrix(c(1, -be, 0, 1), nrow=2)
    Sig <- function(Om, be) crossprod(Gamma(be), Om %*% Gamma(be))

    ## 4. Standard errors based on RE likelihood, or uncorrelated RE likelihood
    ## (for mbtsls, with estimate of Lambda_11 set to zero)
    H.RE <- function(be, Om, lam) {
        q <- Qs(be, Om)
        ch <- lam*q / ((lam+ak)*(nl/d$n))
        bob(be, Om) * (lam+ak) /
            (d$n*lam * (q*Om[2,2]-T[2,2]+ch*q/((1-ch)*aoa(be, Om))) )
    }

    Vmbtsls <- function(Lam, Sig)
        (Sig[1, 1]/Lam + (ak*nl/d$nu)*(det(Sig)+2*Sig[1, 2]^2) / Lam^2)/d$n

    Lam <- if(ei[1] > ak) be.den(mkap[4]) else lam.re/aoa(be['liml'], Om.re)
    Sig.re <- Sig(if (ei[1] > ak) S else Om.re, be['mbtsls'])
    se$re <- sqrt(c(NA, NA, -H.RE(be['liml'], Om.re, lam.re), Vmbtsls(Lam, Sig.re)))

    ## 5. plug-in Standard errors robust to many instruments, assuming Gaussian
    ## errors
    Vliml.N <- function(be, Om, lam)    # Many iv and Normality
        (1+d$k*(nl/d$n)/(lam*d$nu)) * V.info(be, Om, lam)

    Sig.umd <- Sig(S, be['mbtsls'])
    se$manyivN <- sqrt(c(NA, NA, Vliml.N(be['liml'], Om.re, lam.re),
                         Vmbtsls(be.den(mkap[4]), Sig.umd)))

    ## 7. Standard errors based on MD objective function, no rank restriction,
    ## with and without Gaussian errors
    Lam11 <- max(0, c(1,-be['mbtsls']) %*% (T-ak*S) %*% c(1,-be['mbtsls']))
    Lam <- be.den(mkap[4])

    se$invalidivN <- sqrt(c(NA, NA, NA, Vmbtsls(Lam, Sig.umd) +
                                (Lam11*Sig.umd[2, 2] + Lam11*Lam/ak)/(d$n*Lam^2)))

    ################## OTHER OUTPUT
    F <- T[2, 2]/(ak*S[2, 2])              # first-stage F

    ## Sargan and Cragg-Donald
    overid <- if (d$k==1) c(NA, NA) else c(d$n*ei[1]/(d$nu/d$n+ei[1]), d$n*ei[1])
    p.value <- c(1 - pchisq(overid[1], d$k-1),
                 1-pnorm(sqrt(d$nu/nl)*qnorm(pchisq(overid[2],d$k-1))))
    names(overid) <- names(p.value) <- c('Sargan', 'Modified-CD')

    th.re <- list(be=be['liml'],lam=lam.re, Om=Om.re)


    return(list(beta=be, se=se, F=F, overid.test=overid, overid.pvalue=p.value,
                re=th.re, S=S, T=T, df=c(d$n,d$k,d$l)))
}


MDDelta <- function(Lam, be, Om, n, k, l, Psi3=matrix(rep(0, 8), ncol=2),
                      Psi4=matrix(rep(0, 16), ncol=4), delta=0, mu=0) {
    tau <- (k/n) * (n-l)/(n-k-l)

    a <- c(be, 1)
    aa <- a %o% a
    Delta1 <- 2*N2%*% (Lam*kronecker(aa, Om) + Lam*kronecker(Om, aa)+
                          tau*kronecker(Om, Om))
    Delta2 <- (k/n)*delta*(Psi4 - drop(vec(Om)) %o% drop(vec(Om)) -
                               2*N2*kronecker(Om, Om))
    Delta3 <- 2*N2*sqrt(k/n)*mu* kronecker(t(Psi3), a)

    L2 %*% (Delta1+Delta2+Delta3+t(Delta3)) %*% t(L2)
}

MDPsi <- function(Y, d) {
    V <- d$M %*% Y
    Rhat <- d$Zp %*% crossprod(d$Zp, Y)

    vi <- function(i) kronecker(V[i,]%o%V[i,], V[i,])
    Psi3 <- matrix(rowSums(sapply(1:d$n, vi)), ncol=2)/d$m3
    Om <- crossprod(V)/d$nu

    vi <- function(i) kronecker(V[i,]%o%V[i,], V[i,]%o% V[i,])
    Psi4 <- (d$mii-d$m4)*(2*N2*kronecker(Om, Om)+tcrossprod(vec(Om)))
    Psi4 <- (matrix(rowSums(sapply(1:d$n, vi)),ncol=4) - Psi4)/ d$m4

    dH <- rowSums(d$Zp^2)-diag(d$M)*d$k/d$nu
    delta <- drop(crossprod(dH))/d$k
    mu <- drop(crossprod(dH, Rhat[, 2])) / sqrt(d$n*d$k)

    list(Psi3=Psi3, Psi4=Psi4, delta=delta, mu=mu)
}

## Standard errors based on minimum distance objective function
MDse <- function(Wm, Lam, be, Om, n, k, l, Psi3=matrix(rep(0, 8), ncol=2),
                      Psi4=matrix(rep(0, 16), ncol=4), delta=0, mu=0) {

    a <- c(be, 1)
    Del <- MDDelta(Lam, be, Om, n, k, l,
                     Psi3=Psi3, Psi4=Psi4, delta=delta, mu=mu)

    G <- L2 %*% cbind(Lam*(kronecker(a, c(1,0))+kronecker(c(1,0), a)),
                     kronecker(a, a))
    iGWG <- solve(crossprod(G, Wm%*% G))

    sqrt((iGWG %*% crossprod(G, Wm%*%Del%*%Wm %*% G) %*% iGWG)[1,1]/n)
}

UMDse <- function(Xi, be, Om){
    G <- matrix(c(1, Xi[2,2], 0, 1, 0, 0, 0, be, 1),ncol=3)
}

Q.md <- function(be, Lam, Wm, S, T, n, k) {
    M <- vech(T- k/n*S - Lam*c(be, 1)%o%c(be, 1))
    drop(crossprod(M, Wm %*% M))
}

## All possible estimators
MDivregall <- function(Y, Z, W, theta=NULL, desFeat=NULL) {
    d <- if(!is.null(desFeat)) desFeat else MDDesignFeatures(Z, W, Mfeat=TRUE)
    r <- MDivreg(Y, Z, W, theta=theta, desFeat=d)
    n <- r$df[1]
    k <- r$df[2]
    l <- r$df[3]

    ##  MD Normal standard errors

    ## weight matrix
    Wm <- t(D2) %*% kronecker(solve(r$re$Om), solve(r$re$Om)) %*% D2
    Lam <- r$re$lam / aoa(r$re$be, r$re$Om)

    se.md <- MDse(Wm, Lam, r$re$be, r$re$Om, n, k, l)
    r$se$mdN <- c(NA, NA, se.md, NA)

    ## 6. plug-in Standard errors based on MD objective function, not assuming
    ## Gaussian errors

    ## MD non-normal standard errors
    psi <- MDPsi(Y, d)
    se.md <- MDse(Wm, Lam, r$re$be, r$re$Om, n, k, l,
                  Psi3=psi$Psi3, Psi4=psi$Psi4, delta=psi$delta, mu=psi$mu)

    ## TODO: mbtsls: normal, not normal, not valid

    r$se$md <- c(NA, NA, se.md, NA)

    ## RE as MD estimator
    ## W.re <- t(D2) %*% kronecker(solve(r$S), solve(r$S)) %*% D2
    ## qq <- function(t) Q.md(t[1], exp(t[2]), W.re, r$S, r$T, n, k)
    ## il <- optim(c(r$re$be,log(Lam)), qq, method="BFGS")


    r
}

## TODO: Efficient MD, UMD, Testing
