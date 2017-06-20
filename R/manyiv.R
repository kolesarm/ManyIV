## load("/home/kolesarm/research/manyiv-md/code/matchStata/star1.RData")
## dd <- list(Z=Z, W=W, Y=Y)
## rm(Z, W, Y)
## Y <- Matrix::Matrix(dd$Y)
## Z <- Matrix::Matrix(dd$Z)
## W <- Matrix::Matrix(dd$W)
## d <- IVData(Y,Z,W)

bob <- function(be, Om) drop(c(1, -be) %*% Om %*% c(1, -be))
aoa <- function(be, Om) drop(crossprod(c(be, 1), solve(Om, c(be, 1))))
## Duplication, elimination, and commutation matrices, same as
## matrixcalc::duplication.matrix(2),elimination.matrix(2),N.matrix(2)
D2 <- cbind(c(1, 0, 0, 0), c(0, 1, 1, 0), c(0, 0, 0, 1))
L2 <- rbind(c(1, 0, 0, 0), c(0, 1, 0, 0), c(0, 0, 0, 1))
N2 <- rbind(c(1, 0, 0, 0), c(0, 1/2, 1/2, 0 ), c(0, 1/2, 1/2, 0),
            c(0, 0, 0, 1))


## 1. Inference base on invariant and RE likelihood

#' Class constructor for IVData
#'
#' Convert data to standardized format for use with low-level functions. Uses
#' \code{Matrix} package, which speeds up calculations.
#'
#' @param Y [n x 2] class \code{Matrix}
#' @param Z [n x k] Matrix of instruments
#' @param W [n x ell] Matrix of covariates
#' @param moments if \code{TRUE}, compute estimates of third and fourth moments
#'     of the reduced-form errors based on least squares residuals
#' @export
IVData <- function(Y, Z, W, moments=TRUE) {
    ols <- function(X, Y)
        Matrix::solve(Matrix::crossprod(X), Matrix::crossprod(X, Y))

    Zt <- Z-W %*% ols(W, Z)
    d <- list(l=W@Dim[2], k=Zt@Dim[2], n=Zt@Dim[1])
    d$nu <- d$n-d$k-d$l             # degrees of freedom
    X <- Matrix::cBind(Z, W)
    hatPi <- ols(X, Y)
    d$S <- as.matrix(Matrix::crossprod(Y-X%*%hatPi) / d$nu)
    d$T <- as.matrix(Matrix::crossprod(Zt %*% hatPi[1:d$k, ]) / d$n)
    d$ei <- sort(eigen(solve(d$S, d$T))$values) # [m_min, m_max]

    if(moments) {
        d$M <- Matrix::Diagonal(d$n) -
            X %*% Matrix::solve(Matrix::crossprod(X), Matrix::t(X))
        V <- Y-X %*% hatPi
        vi <- function(i) kronecker(V[i, ]%o%V[i, ], V[i, ])
        d$Psi3 <- matrix(rowSums(sapply(1:d$n, vi)), ncol=2)/sum(d$M^3)

        vi <- function(i) kronecker(V[i, ]%o%V[i, ], V[i, ]%o% V[i, ])
        m2 <- sum(Matrix::diag(d$M)^2)
        m4 <- sum(d$M^4)

        d$Psi4 <- (m2-m4)*(2*N2*kronecker(d$S, d$S)+tcrossprod(as.vector(d$S)))
        d$Psi4 <- (matrix(rowSums(sapply(1:d$n, vi)), ncol=4) - d$Psi4)/ m4

        ## Diagonal of H
        h <- Matrix::colSums(Matrix::t(X) *
                             Matrix::solve(Matrix::crossprod(X), Matrix::t(X)))
        h <- h-Matrix::diag(d$M)*d$k/d$nu
        d$delta <- sum(h^2)/d$k
        d$mu <- sum(drop(Zt %*% hatPi[1:d$k, 2]) * h) / sqrt(d$n*d$k)
    }

    structure(d, class="IVData")
}





IVreg.fit <- function(d) {
    ## Estimates of beta
    mkap <- c(-d$nu/d$n, 0, d$ei[1], d$k/d$n)              # m(kappa)
    be.den <- function(m) d$T[2, 2]-m*d$S[2, 2] # denominator
    be <- sapply(mkap, function(m) (d$T[1, 2]- m * d$S[1, 2])/be.den(m))
    names(be) <- c("ols", "tsls", "liml", "mbtsls")

}


MDDelta <- function(be, Om, Xi, d, Gaussian=FALSE) {
    tau <- (d$k/d$n) * (d$n-d$l)/(d$n-d$k-d$l)

    aa <- c(be, 1) %o% c(be, 1)
    Delta1 <- 2*N2%*% (Xi*kronecker(aa, Om) + Xi*kronecker(Om, aa)+
                          tau*kronecker(Om, Om))

    if (Gaussian==TRUE) {
        Delta2 <- Delta3 <- matrix(0, nrow=4, ncol=4)
    } else {
        Delta2 <- (d$k/d$n)*d$delta*(d$Psi4 - as.vector(Om) %o% as.vector(Om) -
                                     2*N2*kronecker(Om, Om))
        Delta3 <- 2*N2*sqrt(d$k/d$n)*d$mu* kronecker(t(d$Psi3), c(be, 1))
    }

    L2 %*% (Delta1+Delta2+Delta3+t(Delta3)) %*% t(L2)
}


#' Minimum distance
IVregMD.fit <- function(d, weight="LIML") {

    ## LIML standard error
    r1 <- IVregRE.fit(d)
    Xi.re <- r1$lam.re / aoa(r1$beta["liml"], r1$Om.re)
    Del <- MDDelta(r1$beta["liml"], r1$Om.re, Xi.re, d)
    DelN <- MDDelta(r1$beta["liml"], r1$Om.re, Xi.re, d, Gaussian=TRUE)
    a.re <- c(r1$beta["liml"], 1)

    G <- L2 %*% cbind(Xi.re*(kronecker(a.re, c(1, 0)) +
                             kronecker(c(1, 0), a.re)), kronecker(a.re, a.re))
    Wm <- t(D2) %*% kronecker(solve(r1$Om.re), solve(r1$Om.re)) %*% D2
    iGWG <- solve(crossprod(G, Wm%*% G))

    ## Without and with assuming normal errors
    se1 <- function(Del) (iGWG %*%
                          crossprod(G, Wm%*%Del%*%Wm %*% G) %*% iGWG)[1, 1]

    liml.se <- sqrt(c(se1(Del), se1(DelN))/d$n)

    obj <- function(be, Xi) {
        r1 <- drop(L2%*%as.vector(d$T-d$k/d$n*d$S-Xi*c(be, 1) %o% c(be, 1)))
        sum(r1*solve(Del, r1))
    }
    start <- c(r1$beta["liml"], Xi.re)
    r2 <- stats::optim(start, function(t) obj(t[1], t[2]))
    se2 <- function(Del) solve(crossprod(G, solve(Del, G)))[1, 1]
    umd.se <- sqrt(c(se2(Del), se2(DelN))/d$n)

    list(beta=r2$par[1], liml.se=liml.se, umd.se=umd.se)
}
