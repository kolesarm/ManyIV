bob <- function(be, Om) drop(c(1, -be) %*% Om %*% c(1, -be))
aoa <- function(be, Om) drop(crossprod(c(be, 1), solve(Om, c(be, 1))))
## Duplication, elimination, and commutation matrices, same as
## matrixcalc::duplication.matrix(2),elimination.matrix(2),N.matrix(2)
D2 <- cbind(c(1, 0, 0, 0), c(0, 1, 1, 0), c(0, 0, 0, 1))
L2 <- rbind(c(1, 0, 0, 0), c(0, 1, 0, 0), c(0, 0, 0, 1))
N2 <- rbind(c(1, 0, 0, 0), c(0, 1/2, 1/2, 0), c(0, 1/2, 1/2, 0),
            c(0, 0, 0, 1))


## Class constructor for \code{"IVData"}
## Convert data to standardized format for use with low-level functions. Uses
## \code{Matrix} package, which speeds up calculations.
## @param Y n-vector
## @param X n-vector
## @param Z [n x k] Matrix of instruments, class \code{Matrix}
## @param W [n x ell] Matrix of covariates, class \code{Matrix}
## @param moments if \code{TRUE}, compute estimates of third and fourth moments
##     of the reduced-form errors based on least squares residuals
## @param approx if \code{TRUE}, then estimates of third and fourth moments use
##     an approximation to speed up the calculations.
IVData <- function(Y, X, Z, W, moments=TRUE, approx=TRUE) {
    ols <- function(X, Y) {
        Matrix::solve(Matrix::crossprod(X), Matrix::crossprod(X, Y))
    }
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

    d$F <- d$T[2, 2] / (d$k/d$n*d$S[2, 2])              # first-stage F

    if (moments) {
        diagP <- function(x) {
            Matrix::colSums(Matrix::t(x) * Matrix::solve(Matrix::crossprod(x),
                                                         Matrix::t(x)))
        }
        diagPX <- diagP(X)
        d$m2 <- sum((1-diagPX)^2)

        XX <- Matrix::solve(Matrix::crossprod(X), Matrix::t(X))
        Hj <- function(j, p) sum((X[j, ] %*% XX)^p)

        ## Split computation into s parts
        s <- max((d$n * (d$k+d$l)) %/% 1e5, 1)
        ix <- cbind((d$n%/%s) * (0:(s-1))+1, (d$n%/%s) * (1:s))
        if ((d$n %% s) > 0)
            ix <- rbind(ix, c((d$n%/%s)*s, d$n))

        if (approx) {
            d$m3 <- d$n - 3 * (d$k+d$l)
            d$m4 <- d$n - 4 * (d$k+d$l)
        } else {
            m3 <- sum(vapply(seq_len(s),
                             function(j) Hj(ix[j, 1]:ix[j, 2], 3), numeric(1)))
            d$m3 <- sum((1-diagPX)^3)+sum(diagPX^3)-m3
            m4 <- sum(vapply(seq_len(s),
                             function(j) Hj(ix[j, 1]:ix[j, 2], 4), numeric(1)))
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
        d$Psi4 <- (d$Psi4 - (d$m2-d$m4) *
                       (2*N2 %*% kronecker(d$S, d$S)+
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

#' Instrumental Variables Regression
#'
#' Fit instrumental variables regression by a number of methods, and compute
#' associated standard errors, as specified by \code{inference}
#' @param formula Specification of the regression relationship and the
#'     instruments of the form \code{y ~ x + w1 + w2 | z1 + z2 + z3}, where
#'     \code{y} is the outcome variable, \code{x} is a scalar endogenous
#'     variable, \code{w1}, \code{w2} are exogenous regressors, and \code{z1},
#'     \code{z2}, and \code{z3} are excluded instruments.
#' @param data An optional data frame, list or environment (or object coercible
#'     by \code{as.data.frame} to a data frame) containing the variables in the
#'     model. If not found in \code{data}, the variables are taken from
#'     \code{environment(formula)}, typically the environment from which the
#'     function is called.
#' @param subset An optional vector specifying a subset of observations to be
#'     used in the fitting process.
#' @param na.action A function indicating what should happen when the data
#'     contain \code{NA}s. The default is set by the \code{na.action} setting of
#'     \code{options} (usually \code{na.omit}).
#' @param approx if \code{TRUE}, then estimates of third and fourth moments used
#'     in inference based on the minimum distance objective function
#'     (\code{inference="md"}) are based on an approximation to speed up the
#'     calculations.
#' @param inference Vector specifying inference method(s). The elements of
#'     the vector can consist of the following methods:
#' \describe{
#'
#'   \item{\code{"standard"}}{Report inference based on TSLS, LIML, and MBTSLS,
#'                along with homoskedastic and heteroskedasticity-robust
#'                standard errors, standard errors that are valid under
#'                heteroskedasticity and treatment effect heterogeneity. All
#'                three standard errors are valid under standard asymptotics
#'                only.}
#'
#'   \item{\code{"re"}}{Report standard errors for LIML based on Hessian of
#'                random effects likelihood}
#'
#'   \item{\code{"il"}}{Report standard errors for LIML based on Hessian of
#'                invariant likelihood, evaluated numerically}
#'
#'   \item{\code{"lil"}}{Report standard errors for LIML based on the
#'                information matrix of limited information likelihood}
#'
#'   \item{\code{"md"}}{Compute the EMD, LIML, and MBTSLS estimators, and report
#'                standard errors for LIML, MBTSLS, and EMD based on the minimum
#'                distance objective function proposed in Kolesár (2018)}
#'
#' }
#' See the vignette \code{vignette("ManyIV", package = "ManyIV")} for a detailed
#' description of these methods.
#' @references {
#'
#' \cite{Kolesár, Michal. Minimum Distance Approach to Inference with Many
#' Instruments.” Journal of Econometrics 204 (1): 86–100.}
#'
#' }
#' @return An object of class \code{"IVResults"}, which is a list with the
#'     following components:
#'
#' \describe{
#'
#' \item{IVData}{An object of class \code{"IVData"}, which is a list with at
#'    least the following components:
#'
#'    \describe{
#'    \item{Z}{Matrix of instruments}
#'    \item{Y}{Matrix with two columns collecting the endogenous variables}
#'    \item{W}{Matrix of exogenous regressors}
#'
#'    \item{n}{Number of observations used, the number of rows of \code{Z},
#'         \code{W}, or \code{Yp}}
#'
#'    \item{l}{Dimension of the exogenous regressors, the number of columns of
#'          \code{W}}
#'
#'    \item{k}{Dimension of the instruments, the number of columns of \code{Z}}
#'
#'    \item{F}{First-stage \eqn{F} statistic}
#'
#'     }
#'
#' }
#'
#' \item{call}{The matched call.}
#'
#' \item{estimate}{A data frame containing the estimation results.}
#' }
#'
#' The \code{print} function can be used to print a summary of the results.
#' @examples
#' ## Use quarter of birth as an instrument for education, controlling for
#' ## marriage and black indicators
#' IVreg(lwage~education+black+married | as.factor(qob),
#'            data=ak80, inference=c("standard", "re", "il", "lil"))
#' @export
IVreg <- function(formula, data, subset, na.action, inference="standard",
                  approx=TRUE) {
    cl <- mf <- match.call(expand.dots = FALSE)
    if (missing(data))
        data <- environment(formula)
    m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE

    formula <- Formula::as.Formula(formula)
    mf$formula <- formula
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())

    Y <- stats::model.response(mf, "numeric")
    mtX <- stats::terms(formula, data = data, rhs = 1)
    Xname <- attr(mtX, "term.labels")[1]
    W <- stats::model.matrix(mtX, mf)
    X <- W[, Xname]
    W <- Matrix::Matrix(W[, !(colnames(W) %in% Xname), drop=FALSE])
    mtZ <- stats::delete.response(stats::terms(formula, data = data, rhs = 2))
    Z <- stats::model.matrix(mtZ, mf)
    ## Remove intercept
    Z <- Matrix::Matrix(Z[, !colnames(Z) %in% colnames(W), drop=FALSE])

    d <- IVData(Y, X, Z, W, moments = ("md" %in% inference), approx=approx)

    structure(list(IVData=d, call=cl,
                   estimate=ivreg.fit(d, inference=inference)),
              class="IVResults")
}

## Low-level computing engine called by \code{IVreg}
## @param d Object of class \code{"IVData"}
## @return A data frame containing the estimation results.
ivreg.fit <- function(d, inference) {
    est <- data.frame(row.names=c("ols", "tsls", "liml", "mbtsls", "emd"))

    if ("standard" %in% inference) {
        r <- si.fit(d)
        est[-5, "estimate"] <- r$estimate
        est[-5, "se"] <- r$se
        est[-5, "ser"] <- r$ser
        est[-5, "seh"] <- r$seh
    }

    if ("lil" %in% inference) {
        r <- li.fit(d)
        est["liml", "estimate"] <- r$estimate
        est["liml", "lil"] <- r$se
    }

    if ("re" %in% inference) {
        r <- re.fit(d)
        est["liml", "estimate"] <- r$estimate
        est["liml", "re"] <- r$se
    }

    if ("il" %in% inference) {
        r <- il.fit(d)
        est["liml", "estimate"] <- r$estimate
        est["liml", "il"] <- r$se
    }

    if ("md" %in% inference) {
        r <- md.fit(d, weight="LIML")
        est["liml", "estimate"] <- r$estimate
        est["liml", "md"] <- r$se
        r <- md.fit(d, weight="Optimal")
        est["emd", "estimate"] <- r$estimate
        est["emd", "md"] <- r$se
        r <- umd.fit(d, invalid=FALSE)
        est["mbtsls", "estimate"] <- r$estimate
        est["mbtsls", "md"] <- r$se
        est["mbtsls", "umd"] <- umd.fit(d, invalid=TRUE)$se
    }

    est
}

#' Test of overidentifying restrictions
#'
#' Report the Sargan and modified Cragg-Donald test statistics and
#' \eqn{p}-values for testing of overidentifying restrictions, assuming
#' homoskedasticity of the reduced form. The Sargan test is valid under few
#' instruments. The Modified Cragg-Donald test (Modified-CD) corresponds to a
#' test due to Cragg and Donald (1993), with a modified critical value. The
#' modification was suggested in Kolesár (2018) to make it robust to many
#' instruments and many exogenous regressors.
#' @param r An object of class \code{RDResults}
#' @examples
#' r1 <- IVreg(lwage~education+black+married | as.factor(qob), data=ak80,
#'             inference="standard")
#' IVoverid(r1)
#' @references {
#'
#' \cite{Kolesár, Michal. Minimum Distance Approach to Inference with Many
#' Instruments.” Journal of Econometrics 204 (1): 86–100.
#' \doi{10.1016/j.jeconom.2018.01.004}.}
#'
#' \cite{Cragg, John G., and Stephen G. Donald. 1993. "Testing Identifiability
#' and Specification in Instrumental Variable Models." Econometric Theory 9 (2):
#' 222–40. \doi{10.1017/S0266466600007519}.}
#'
#' \cite{Sargan, John Denis. 1958. "The Estimation of Economic Relationships
#' Using Instrumental Variables." Econometrica 26 (3): 393–415.
#' \doi{10.2307/1907619}.}
#'
#' }
#' @export
IVoverid <- function(r) ivoverid.fit(r$IVData)

ivoverid.fit <- function(d) {
    if (!("Psi4" %in% names(d)))
        d <- IVData(d$Y[, 1], d$Y[, 2], d$Z, d$W, moments=TRUE)

    ## Sargan and Cragg-Donald
    if (d$k==1) {
        overid <- c(NA, NA)
    } else {
        overid <- c(d$n*d$ei[1] / (d$nu/d$n+d$ei[1]), d$n*d$ei[1])
    }

    r1 <- re.fit(d)
    a <- c(r1$estimate, 1)
    kap <- drop(crossprod(kronecker(a, a), d$Psi4 %*% kronecker(a, a))) /
        bob(r1$estimate, r1$Om)^2 - 3
    rr <- (d$n-d$l)/d$nu+d$delta*kap/2
    p.value <- c(1 - stats::pchisq(overid[1], d$k-1),
                 1-stats::pnorm(stats::qnorm(stats::pchisq(overid[2],
                                                           d$k-1))/sqrt(rr)))
    names(overid) <- names(p.value) <- c("Sargan", "Modified-CD")
    data.frame(statistic=overid, p.value=p.value)
}


si.fit <- function(d) {
    ## 1. Estimates of beta
    mkap <- c(-d$nu/d$n, 0, d$ei[1], d$k/d$n)   # m(kappa)
    be.den <- function(m) d$T[2, 2]-m*d$S[2, 2] # denominator
    be <- vapply(mkap, function(m) (d$T[1, 2]- m * d$S[1, 2])/be.den(m),
                 numeric(1))

    ## 2. Stata standard errors
    he <- function(be) d$Yt[, 1]- d$Yt[, 2]*be
    hat.sig <- function(be) drop(crossprod(he(be))) / d$n
    se <- sqrt(vapply(be, hat.sig, numeric(1)) /
                   (d$n*pmax(be.den(mkap), 0)))
    sm <- d$n / (d$n - d$l-1) # Small-sample adjustment for OLS only
    se[1] <- se[1]*sqrt(sm)

    ## 3. Stata robust standard errors
    num <- vapply(2:4, function(j) sum((d$Yhatp[, 2]*he(be[j]))^2),
                  numeric(1))
    ser <- sqrt(c(drop(crossprod(d$Yt[, 2]*he(be[1])))*sm, num)) /
        pmax(d$n*pmax(be.den(mkap), 0))
    ## 4. Standard errors under TE heteroegeneity
    num <- function(j) {
        sum((d$Yhatp[, 2]*he(be[j])+ (d$Yhatp[, 1]-d$Yhatp[, 2]*be[j]) *
                 (d$Yt[, 2]-d$Yhatp[, 2]))^2)
    }
    seh <- sqrt(c(NA, num(2), NA, num(4)))/pmax(d$n*pmax(be.den(mkap), 0))

    names(be) <- names(se) <- names(ser) <- names(seh) <-
        c("ols", "tsls", "liml", "mbtsls")
    list(estimate=be, se=se, ser=ser, seh=seh, lam=NA, Om=d$S)
}


MDDelta <- function(be, Om, Xi22, d, Gaussian=FALSE, invalid=FALSE) {
    tau <- (d$k/d$n) * (d$n-d$l) / (d$n-d$k-d$l)

    Xim <- if (invalid) (d$T - (d$k/d$n)*d$S) else (Xi22*c(be, 1) %o% c(be, 1))
    D1 <- 2*N2 %*%  (kronecker(Xim, Om) + kronecker(Om, Xim) +
                         tau*kronecker(Om, Om))

    if (Gaussian) {
        D2 <- D3 <- matrix(0, nrow=4, ncol=4)
    } else {
        D2 <- (d$k/d$n) * d$delta * (d$Psi4 - as.vector(Om) %o% as.vector(Om) -
                                         2*N2%*%kronecker(Om, Om))
        mu1 <- if (invalid) d$mu1 else be*d$mu
        D3 <- 2*sqrt(d$k/d$n) * N2 %*% kronecker(t(d$Psi3), c(mu1, d$mu))
    }

    L2 %*% (D1+D2+D3+t(D3)) %*% t(L2)
}


## Minimum distance
md.fit <- function(d, weight="Optimal") {

    ## LIML standard error
    r1 <- re.fit(d)
    Xi <- r1$lam / aoa(r1$estimate, r1$Om)
    Del <-  MDDelta(r1$estimate, r1$Om, Xi, d)
    a.re <- c(r1$estimate, 1)

    G <- L2 %*% cbind(Xi * (kronecker(a.re, c(1, 0)) +
                                kronecker(c(1, 0), a.re)),
                      kronecker(a.re, a.re))
    Wm <- crossprod(D2, kronecker(solve(r1$Om), solve(r1$Om)) %*% D2)
    iGWG <- solve(crossprod(G, Wm %*% G))

    if (weight=="LIML") {
        se <- (iGWG %*% crossprod(G, Wm%*%Del%*%Wm %*% G) %*% iGWG)[1, 1]
        return(list(estimate=r1$estimate, se=sqrt(se/d$n), lam=r1$lam,
                    Om=r1$Om))
    }

    obj <- function(be, Xi) {
        r <- drop(L2%*%as.vector(d$T-d$k/d$n*d$S-Xi*c(be, 1) %o% c(be, 1)))
        sum(r*solve(Del, r))
    }
    start <- c(r1$estimate, Xi)
    r2 <- stats::optim(start, function(t) obj(t[1], t[2]), method="BFGS")
    se2 <- solve(crossprod(G, solve(Del, G)))[1, 1]

    list(estimate=r2$par[1], se=sqrt(se2/d$n),
         lam=r2$par[2]*aoa(r2$par[1], r1$Om),
         Om=r1$Om)
}

## Unrestricted minimum distance
umd.fit <- function(d, invalid=TRUE) {
    Xi <- d$T - (d$k/d$n)*d$S
    be <- Xi[1, 2] / Xi[2, 2]
    Om <- d$S

    G <- c(0, 1, -be) / Xi[2, 2]
    se <- drop(crossprod(G, MDDelta(be, Om, Xi[2, 2], d,
                                    invalid=invalid) %*% G))

    list(estimate=be, se=sqrt(se/d$n), lam=Xi[2, 2]*aoa(be, Om), Om=Om)
}
