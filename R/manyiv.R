bob <- function(be, Om) drop(c(1, -be) %*% Om %*% c(1, -be))
aoa <- function(be, Om) drop(crossprod(c(be, 1), solve(Om, c(be, 1))))
## Duplication, elimination, and commutation matrices, same as
## matrixcalc::duplication.matrix(2),elimination.matrix(2),N.matrix(2)
D2 <- cbind(c(1, 0, 0, 0), c(0, 1, 1, 0), c(0, 0, 0, 1))
L2 <- rbind(c(1, 0, 0, 0), c(0, 1, 0, 0), c(0, 0, 0, 1))
N2 <- rbind(c(1, 0, 0, 0), c(0, 1/2, 1/2, 0 ), c(0, 1/2, 1/2, 0),
            c(0, 0, 0, 1))


## 1. Inference based on invariant and RE likelihood

#' Class constructor for IVData
#'
#' Convert data to standardized format for use with low-level functions. Uses
#' \code{Matrix} package, which speeds up calculations.
#'
#' @param Y n-vector
#' @param X n-vector
#' @param Z [n x k] Matrix of instruments, class \code{Matrix}
#' @param W [n x ell] Matrix of covariates, class \code{Matrix}
#' @param moments if \code{TRUE}, compute estimates of third and fourth moments
#'     of the reduced-form errors based on least squares residuals
#' @param approx if \code{TRUE}, then estimates of third and fourth moments use
#'     an appriximation to speed up the calculations.
#' @import Matrix
#' @export
IVData <- function(Y, X, Z, W, moments=TRUE, approx=TRUE) {
    ols <- function(X, Y)
        Matrix::solve(Matrix::crossprod(X), Matrix::crossprod(X, Y))
    d <- list(l=W@Dim[2], k=Z@Dim[2], n=Z@Dim[1], Z=Z, W=W,
              Y=cbind(Y, X))
    Zt <- if (d$l==0) Z else (Z-W %*% ols(W, Z))
    d$Yt <- if (d$l==0) d$Y else (d$Y-W %*% ols(W, d$Y))

    d$nu <- d$n-d$k-d$l             # degrees of freedom
    X <- cbind(Z, W)
    hatPi <- ols(X, d$Y)
    d$S <- as.matrix(Matrix::crossprod(d$Y-X%*%hatPi) / d$nu)
    d$Yhatp <- Zt %*% hatPi[1:d$k, ]
    d$T <- as.matrix(Matrix::crossprod(d$Yhatp) / d$n)
    d$ei <- sort(eigen(solve(d$S, d$T))$values) # [m_min, m_max]

    d$F <- d$T[2, 2]/(d$k/d$n*d$S[2, 2])              # first-stage F

    if(moments) {
        diagP <- function(x)
            Matrix::colSums(Matrix::t(x) *
                            Matrix::solve(Matrix::crossprod(x), Matrix::t(x)))
        diagPX <- diagP(X)
        d$m2 <- sum((1-diagPX)^2)

        XX <- Matrix::solve(Matrix::crossprod(X), Matrix::t(X))
        Hj <- function(j, p) sum((X[j, ] %*% XX)^p)

        ## Split computation into s parts
        s <- max((d$n*(d$k+d$l)) %/% 1e5, 1)
        ix <- cbind((d$n%/%s)*(0:(s-1))+1, (d$n%/%s)*(1:s))
        if ((d$n %% s) > 0)
            ix <- rbind(ix, c((d$n%/%s)*s, d$n))

        if (approx) {
            d$m3 <- d$n-3*(d$k+d$l)
            d$m4 <- d$n-4*(d$k+d$l)
        } else {
            m3 <- sum(sapply(1:s, function(j) Hj(ix[j, 1]:ix[j, 2], 3)))
            d$m3 <- sum((1-diagPX)^3)+sum(diagPX^3)-m3
            m4 <- sum(sapply(1:s, function(j) Hj(ix[j, 1]:ix[j, 2], 4)))
            d$m4 <- sum((1-diagPX)^4)+m4-sum(diagPX^4)
        }

        V <- d$Y-X %*% hatPi
        d$Psi3 <- c(sum(V[, 1]^3), sum(V[, 1]^2*V[, 2]), sum(V[, 1]*V[, 2]^2),
                    sum(V[, 2]^3))
        d$Psi3 <- cbind(c(d$Psi3[1], d$Psi3[2], d$Psi3[2], d$Psi3[3]),
                        c(d$Psi3[2], d$Psi3[3], d$Psi3[3], d$Psi3[4])) / d$m3

        d$Psi4 <- c(sum(V[, 1]^4), sum(V[, 1]^3*V[, 2]), sum(V[, 1]^2*V[, 2]^2),
                    sum(V[, 1]*V[, 2]^3), sum(V[, 2]^4))
        d$Psi4 <- cbind(c(d$Psi4[1], d$Psi4[2], d$Psi4[2], d$Psi4[3]),
                        c(d$Psi4[2], d$Psi4[3], d$Psi4[3], d$Psi4[4]),
                        c(d$Psi4[2], d$Psi4[3], d$Psi4[3], d$Psi4[4]),
                        c(d$Psi4[3], d$Psi4[4], d$Psi4[4], d$Psi4[5]))
        d$Psi4 <- (d$Psi4 - (d$m2-d$m4) * (2*N2 %*% kronecker(d$S, d$S)+
                                           tcrossprod(as.vector(d$S)))) / d$m4

        ## Diagonal of H
        h <- diagP(Zt)-d$k/d$nu * (1-diagPX)
        d$delta <- sum(h^2)/d$k
        d$mu <- sum(drop(Zt %*% hatPi[1:d$k, 2]) * h) / sqrt(d$n*d$k)
        ## For invalid IV
        d$mu1 <- sum(drop(Zt %*% hatPi[1:d$k, 1]) * h) / sqrt(d$n*d$k)
    }

    structure(d, class="IVData")
}

#' Fit instrumental-variable regression
#' @param formula specification of the regression relationship and the
#'     instruments of the form \code{y ~ x + w1 + w2 | z1 + z2 + z3}, where
#'     \code{y} is the outcome variable, \code{x} is a scalar endogenous
#'     variable, \code{w1,w2} are exogenous regressors, and \code{z1,z2,z3} are
#'     excluded instruments.
#' @param data optional data frame, list or environment (or object coercible by
#'     \code{as.data.frame} to a data frame) containing the outcome and running
#'     variables in the model. If not found in \code{data}, the variables are
#'     taken from \code{environment(formula)}, typically the environment from
#'     which the function is called.
#' @param subset optional vector specifying a subset of observations to be used
#'     in the fitting process.
#' @param na.action function which indicates what should happen when the data
#'     contain \code{NA}s. The default is set by the \code{na.action} setting of
#'     \code{options} (usually \code{na.omit}).
#' @param approx if \code{TRUE}, then estimates of third and fourth moments used
#'     in inference based on the minimum distance objective function
#'     (\code{inference="md"}) use an appriximation to speed up the
#'     calculations.
#' @inheritParams IVreg.fit
#' @examples
#' ## Specification as in Table V, columns (1) and (2) in Angrist and Krueger
#' IVreg(lwage~education+as.factor(yob)|as.factor(qob)*as.factor(yob),
#'            data=ak80, inference=c("standard", "re", "il", "lil"))
#' ## Only quarter of birth as instrument, add married, black and smsa as exogenous
#' #regressors
#' IVreg(lwage~education+as.factor(yob)+black+smsa+married|as.factor(qob),
#'            data=ak80, inference=c("standard", "re", "il", "lil"))
#' @export
IVreg <- function (formula, data, subset, na.action, inference="standard",
                   approx=TRUE) {
    formula <- Formula::as.Formula(formula)
    cl <- mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$formula <- formula
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())

    Y <- stats::model.response(mf, "numeric")
    mtX <- stats::terms(formula, data = data, rhs = 1)
    Xname <- attr(mtX, "term.labels")[1]
    W <- stats::model.matrix(mtX, mf)
    X <- W[, Xname]
    W <- Matrix::Matrix(W[, !(colnames(W) %in% Xname)])
    mtZ <- stats::delete.response(stats::terms(formula, data = data, rhs = 2))
    Z <- stats::model.matrix(mtZ, mf)
    ## Remove intercept now, using attr(mtZ, "intercept") <- 0 won't help if Z
    ## consists of indicators
    Z <- Matrix::Matrix(Z[, !colnames(Z) %in% colnames(W)])

    d <- IVData(Y, X, Z, W, moments=("md" %in% inference), approx=approx)

    ret <- structure(list(IVData=d, call=cl, formula=formula(formula),
                                            na.action=attr(mf, "na.action")),
                     class="IVResults")

    ret$estimate <- IVreg.fit(d, inference=inference)

    ret
}

#' Low-level interfact to \code{IVreg}
#' @param d Object of class \code{"IVData"}
#' @param inference Vector specifying inference method(s). The elements of
#'     the vector can consist of the following methods:
#' \describe{
#'
#'   \item{"standard"}{Report inference based on tsls, liml, and mbtsls, along
#'                     with homoscedastic heteroscedasticity-robust standard
#'                     errors valid under standard asymptotic sequence}
#'
#'   \item{"re"}{Inference based on Hessian of random effects likelihood}
#'
#'   \item{"il"}{Inference based on Hessian of invariant likelihood, evaluated
#'               numerically}
#'
#'   \item{"lil"}{Inference based on information matrix of limited
#'               information likelihood}
#'
#'   \item{"md"}{Inference based on the minimum distance objective function} }
#' @export
IVreg.fit <- function(d, inference) {
    est <- data.frame(row.names=c("ols", "tsls", "liml", "mbtsls", "emd"))

    if ("standard" %in% inference) {
        r <- IVregSI.fit(d)
        est[-5, "beta"] <- r$beta
        est[-5, "se"] <- r$se
        est[-5, "ser"] <- r$ser
    }

    if ("lil" %in% inference) {
        r <- IVregLI.fit(d)
        est["liml", "beta"] <- r$beta
        est["liml", "lil"] <- r$se
    }

    if ("re" %in% inference) {
        r <- IVregRE.fit(d)
        est["liml", "beta"] <- r$beta
        est["liml", "re"] <- r$se
    }

    if ("il" %in% inference) {
        r <- IVregIL.fit(d)
        est["liml", "beta"] <- r$beta
        est["liml", "il"] <- r$se
    }

    if ("md" %in% inference) {
        r <- IVregMD.fit(d, weight="LIML")
        est["liml", "beta"] <- r$beta
        est["liml", "md"] <- r$se
        r <- IVregMD.fit(d, weight="Optimal")
        est["emd", "beta"] <- r$beta
        est["emd", "md"] <- r$se
        r <- IVregUMD.fit(d, invalid=FALSE)
        est["mbtsls", "beta"] <- r$beta
        est["mbtsls", "md"] <- r$se
        est["mbtsls", "umd"] <- IVregUMD.fit(d, invalid=TRUE)$se
    }

    est
}

#' Test of overidentifying restrictions
#'
#' Report test statistic and p-value for testing of overidentifying
#' restrictions. The Sargan test is valid under few instruments. The Modified
#' Cragg-Donald test (Modified-CD) corresponds to a test due to Cragg and Donald
#' (1993), with a critical value modified to make it robust to many instruments
#' and many exogenous regressors.
#' @param r object of class \code{RDResults}
#' @export
IVoverid <- function(r) IVoverid.fit(r$IVData)


IVoverid.fit <- function(d) {
    if (!("Psi4" %in% names(d)))
        d <- IVData(d$Y[, 1], d$Y[, 2], d$Z, d$W, moments=TRUE)

    ## Sargan and Cragg-Donald
    overid <- if (d$k==1) {
                  c(NA, NA)
              } else {
                  c(d$n*d$ei[1]/(d$nu/d$n+d$ei[1]), d$n*d$ei[1])
              }

    r1 <- IVregRE.fit(d)
    a <- c(r1$beta, 1)
    kap <- drop(crossprod(kronecker(a, a), d$Psi4 %*% kronecker(a, a))) /
        bob(r1$beta, r1$Om)^2 - 3
    rr <- (d$n-d$l)/d$nu+d$delta*kap/2
    p.value <- c(1 - stats::pchisq(overid[1], d$k-1),
                 1-stats::pnorm(stats::qnorm(stats::pchisq(overid[2],
                                                           d$k-1))/sqrt(rr)))
    names(overid) <- names(p.value) <- c("Sargan", "Modified-CD")
    data.frame(statistic=overid, p.value=p.value)
}


IVregSI.fit <- function(d) {
    ## 1. Estimates of beta
    mkap <- c(-d$nu/d$n, 0, d$ei[1], d$k/d$n)   # m(kappa)
    be.den <- function(m) d$T[2, 2]-m*d$S[2, 2] # denominator
    be <- sapply(mkap, function(m) (d$T[1, 2]- m * d$S[1, 2])/be.den(m))

    ## 2. Stata standard errors
    he <- function(be) d$Yt[, 1]- d$Yt[, 2]*be
    hat.sig <- function(be) drop(crossprod(he(be))) / d$n
    se <- sqrt(sapply(be, hat.sig) / (d$n*pmax(be.den(mkap), 0)))
    sm <- d$n/(d$n - d$l-1)
    se[1] <- se[1]*sqrt(sm)

    ## 3. Stata robust standard errors
    num <- sapply(2:4, function(j) drop(crossprod(d$Yhatp[, 2]*he(be[j]))))
    ser <- sqrt(c(drop(crossprod(d$Yt[, 2]*he(be[1])))*sm, num)) /
        pmax(d$n*pmax(be.den(mkap), 0))

    names(be) <- names(se) <- names(ser) <- c("ols", "tsls", "liml", "mbtsls")
    list(beta=be, se=se, ser=ser, lam=NA, Om=d$S)
}


MDDelta <- function(be, Om, Xi22, d, Gaussian=FALSE, invalid=FALSE) {
    tau <- (d$k/d$n) * (d$n-d$l)/(d$n-d$k-d$l)

    Xim <- if (invalid) (d$T-(d$k/d$n)*d$S) else (Xi22*c(be, 1) %o% c(be, 1))
    D1 <- 2*N2 %*%  (kronecker(Xim, Om) + kronecker(Om, Xim) +
                     tau*kronecker(Om, Om))

    if (Gaussian) {
        D2 <- D3 <- matrix(0, nrow=4, ncol=4)
    } else {
        D2 <- (d$k/d$n)*d$delta*(d$Psi4 - as.vector(Om) %o% as.vector(Om) -
                                     2*N2%*%kronecker(Om, Om))
        mu1 <- if (invalid) d$mu1 else be*d$mu
        D3 <- 2*sqrt(d$k/d$n) * N2 %*% kronecker(t(d$Psi3), c(mu1, d$mu))
    }

    L2 %*% (D1+D2+D3+t(D3)) %*% t(L2)
}


#' Minimum distance
#' @keywords internal
IVregMD.fit <- function(d, weight="Optimal") {

    ## LIML standard error
    r1 <- IVregRE.fit(d)
    Xi.re <- r1$lam / aoa(r1$beta, r1$Om)
    Del <-  MDDelta(r1$beta, r1$Om, Xi.re, d)
    a.re <- c(r1$beta, 1)

    G <- L2 %*% cbind(Xi.re*(kronecker(a.re, c(1, 0)) +
                             kronecker(c(1, 0), a.re)), kronecker(a.re, a.re))
    Wm <- crossprod(D2, kronecker(solve(r1$Om), solve(r1$Om)) %*% D2)
    iGWG <- solve(crossprod(G, Wm %*% G))

    if (weight=="LIML") {
        se <- (iGWG %*% crossprod(G, Wm%*%Del%*%Wm %*% G) %*% iGWG)[1, 1]
        return(list(beta=r1$beta, se=sqrt(se/d$n), lam=r1$lam, Om=r1$Om))
    }

    obj <- function(be, Xi) {
        r <- drop(L2%*%as.vector(d$T-d$k/d$n*d$S-Xi*c(be, 1) %o% c(be, 1)))
        sum(r*solve(Del, r))
    }
    start <- c(r1$beta, Xi.re)
    r2 <- stats::optim(start, function(t) obj(t[1], t[2]), method="BFGS")
    se2 <- solve(crossprod(G, solve(Del, G)))[1, 1]

    list(beta=r2$par[1], se=sqrt(se2/d$n), lam=r2$par[2]*aoa(r2$par[1], r1$Om),
         Om=r1$Om)
}

#' Unrestriced minimum distance
#' @keywords internal
IVregUMD.fit <- function(d, invalid=TRUE) {
    Xi <- d$T - (d$k/d$n)*d$S
    be <- Xi[1, 2] / Xi[2, 2]
    Om <- d$S

    G <- c(0, 1, -be) / Xi[2, 2]
    se <- drop(crossprod(G, MDDelta(be, Om, Xi[2, 2], d,
                                    invalid=invalid) %*% G))

    list(beta=be, se=sqrt(se/d$n), lam=Xi[2, 2]*aoa(be, Om), Om=Om)
}


#' @export
print.IVResults <- function(x, digits = getOption("digits"), ...) {
    if (!is.null(x$call))
        cat("Call:\n", deparse(x$call), sep = "", fill=TRUE)

    r <- x$estimate[!is.na(x$estimate$beta), ]
    colnames(r) <- c("Estimate", colnames(r)[-1])

    if ("se" %in% colnames(r))
        colnames(r) <- c(c("Estimate", "Conventional", "Conv. (robust)"),
                         colnames(r[, -(1:3)]))

    cat("\nFirst-stage F: ", x$IVData$F, "\n\n")
    cat("Estimates and standard errors:\n")

    print.data.frame(r, digits = digits, ...)


    invisible(x)
}
