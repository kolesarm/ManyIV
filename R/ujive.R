#' Jackknife Instrumental variables regression
#'
#' TODO: full description
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
#' @param tol Numerical tolerance for determining rank of instruments and
#'     covariates
#' @return An object of class \code{"IVResults"}.
#' @export
ujive <- function(formula, data, subset, na.action, tol=1e-8) {
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
    W <- Matrix::Matrix(Matrix::sparse.model.matrix(mtX, mf,
                                                    drop.unused.levels = TRUE))
    Didx <- 1+attr(mtX, "intercept")
    D <- W[, Didx]
    W <- standardize(W[, -Didx, drop=FALSE])
    mtZ <- stats::delete.response(stats::terms(formula, data = data, rhs = 2))
    Z <- Matrix::Matrix(Matrix::sparse.model.matrix(mtZ, mf,
                                                    drop.unused.levels = TRUE))
    ## Remove intercept
    Z <- standardize(Matrix::Matrix(Z[, !colnames(Z) %in% colnames(W),
                                      drop=FALSE]))

    ## qrX and qrW will typically be sparse, so we will work with those
    r <- remove_collinear(W, Z, tol)

   if (length(r$drp) > 0) {
        Y <- Y[-r$drp]
        D <- D[-r$drp]
        Z <- Z[-r$drp, ]
        W <- W[-r$drp, ]
    }
    if (NCOL(r$qrW)<NCOL(W))
        message("Dropped ", NCOL(W)-NCOL(r$qrW), " collinear controls.")
    if (NCOL(r$qrX)-NCOL(r$qrW)<NCOL(Z))
        message("Dropped ", NCOL(Z)-NCOL(r$qrX)+NCOL(r$qrW), " collinear IVs.")

    ret <- ujive.fit(Y, D, r$qrX, r$qrW, r$dX, tol)
    structure(list(IVData=list(Y=Y, D=D, Z=Z, W=W, F=ret$F, k=ret$k, n=ret$n,
                               l=ret$l), call=cl, drop_obs = r$drp,
                   estimate=ret$r),
              class="IVResults")
}

## Standardize non-binary columns
standardize <- function(X) {
    idx <- which(Matrix::colSums(X!=0 & X!=1)>0)
    if (length(idx)>0)
        X[, idx] <- scale(X[, idx])
    X
}

remove_collinear <- function(W, Z, tol) {
    X <- cbind(W, Z)
    rownames(X) <- seq_len(NROW(X))

    ## 1. First remove leverage one observations
    ## 1a. Remove singletons, this is fast
    rs <- drp_singleton(X, tol)
    drp_obs <- rs$drp

    ## 1b. Remove obs based on leverage, first dropping collinear cols
    qrX <- rr_qr(rs$X, tol)
    dX <- unname(Matrix::rowSums(Matrix::qr.Q(qrX)^2))
    drp_lev <- which(dX>1-tol)
    if (length(drp_lev)>0) {
        drp_obs <- c(drp_obs, rownames(rs$X[drp_lev, ]))
        message("Dropping ", length(drp_lev), " obs with leverage 1.")
        ## Update qrX and projection diagonal
        qrX <- rr_qr(drp_singleton(rs$X[-drp_lev, ], tol)$X, tol)
        dX <- unname(Matrix::rowSums(Matrix::qr.Q(qrX)^2))
    }
    ## 2. qr decomposition of W
    if (length(drp_obs)>0) {
        W <- W[-as.numeric(drp_obs), , drop=FALSE]
    }
    qrW <- rr_qr(drp_singleton(W, tol)$X, tol)

    list(qrX=qrX, qrW=qrW, dX=dX, drp=as.numeric(drp_obs))
}

## qr decomposition, allowing for reduced rank matrix. Drops collinear columns
rr_qr <- function(A, tol) {
    qrA <- Matrix::qr(A)
    dr <- drp(qrA, tol)
    if (length(dr)>0) {
        A <- A[, -dr, drop=FALSE]
        qrA <- Matrix::qr(A)
    }
    qrA
}


## Drop singletons
drp_singleton <- function(X, tol) {
    idx1 <- which(Matrix::colSums(X!=0)==1)
    idx_obs <- vector()
    while (length(idx1) > 0) {
        drp <- which(Matrix::rowSums(abs(X[, idx1, drop=FALSE]))>0)
        idx_obs <- c(idx_obs, drp)
        X <- X[-drp, ]
        idx1 <- which(Matrix::colSums(X!=0)==1)
    }
    if (length(idx_obs)>0)
        message("Recursively dropping ", length(unique(names(idx_obs))),
                " obs with singleton covariates or IVs")
    idx0 <- which(Matrix::colSums(X!=0)==0)
    if (length(idx0) > 0) X <- X[, -idx0, drop=FALSE]

    list(X=X, drp=unique(names(idx_obs)), col_idx=names(idx0))
}



drp <- function(qrA, tol) {
    isqrA <- is.qr(qrA)
    ## Diagonal of the R matrix
    dR <- abs(if (isqrA) diag(qrA$qr) else abs(Matrix::diag(qrA@R)))
    d <- max(dim(if (isqrA) qrA$qr else qrA))
    dr <- which(dR < d * tol * max(dR))
    ## For some reason, the permutation vector q starts at 0
    ## ea <- Matrix::expand2(qrA, complete = FALSE)
    ## ea$P1. %*% ea$Q1 %*% ea$R1 == A[, qrA@q+1]
    ret <- if (isqrA) qrA$pivot[dr] else (qrA@q+1)[dr]
    ret
}






ujive.fit <- function(Y, D, qrX, qrW, dX, tol) {
    ## Diagonals of projection matrices
    dW <- Matrix::rowSums(Matrix::qr.Q(qrW)^2)
    stopifnot(max(dW)<1-tol)

    dM <- 1-dX
    MX <- function(A) Matrix::qr.resid(qrX, A)

    ## Treatment residuals
    HD <- Matrix::qr.fitted(qrX, D)-Matrix::qr.fitted(qrW, D)
    DtW <- Matrix::qr.resid(qrW, D)
    MD <- DtW-HD
    ## Zero out cases with leverage one
    if (sum(dX > 1-tol)>0)
        warning(paste("There are observations with leverage one",
                      "cannot compute UJIVE."))

    Dujive <- HD - MD * (dX-dW)/dM
    Djive1 <- DtW-Matrix::qr.resid(qrW, MD/dM)
    Dojive <- DtW / (1-dW) - MD/dM
    Dijive <- DtW-Matrix::qr.resid(qrW, MD / (dM+dW))

    ## TSLS, UJIVE, old UJIVE, IJIVE, JIVE1
    Dlist <- list(DtW, HD, Dujive, Dojive, Dijive, Djive1)
    den <- vapply(Dlist, function(x) sum(x*D), numeric(1))
    num <- vapply(Dlist, function(x) sum(x*Y), numeric(1))
    est <- num/den
    names(est) <- c("ols", "tsls", "ujive", "old ujive", "ijive1", "jive1")

    ## Textbook robust
    YtW <- Matrix::qr.resid(qrW, Y)
    epD <- (YtW - DtW %o% est)*do.call(cbind, Dlist) # epsilon*Dhat
    text <- Matrix::colSums(epD^2)
    ## HTE robust
    HY <- Matrix::qr.fitted(qrX, Y)-Matrix::qr.fitted(qrW, Y)
    GYtsls <- HY-HD*est[2]
    GYujive <- HY-HD*est[3] - MX((Y-D*est[3]) * (dX-dW) / dM)
    GYjive1 <- YtW-DtW*est[6] - MX((YtW-DtW*est[6]) / dM)
    GYojive <- Matrix::qr.resid(qrW, (Y-D*est[4]) / (1-dW))-
        MX((Y-D*est[4]) / dM)
    GYijive <- YtW-DtW*est[5]-MX((YtW-DtW*est[5]) / (dM+dW))

    GY <- cbind(0, GYtsls, GYujive, GYojive, GYijive, GYjive1)
    hte <- Matrix::colSums((GY*MD+epD)^2)
    r <- cbind(estimate=est, se_text=sqrt(text/den^2), se_hte=sqrt(hte/den^2))
    n <- length(Y)
    kl <- if (is.qr(qrX)) qrX$rank else dim(qrX)[2]
    k <- kl - if (is.qr(qrW)) qrW$rank else dim(qrW)[2]
    Fhom <- sum(HD^2)/sum(MD^2) * (n-kl) / k

    list(r=as.data.frame(r), F=Fhom, n=n, k=k, l=kl-k)
}
