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
    Dname <- attr(mtX, "term.labels")[1]
    W <- Matrix::Matrix(Matrix::sparse.model.matrix(mtX, mf))
    D <- W[, Dname]
    W <- W[, !(colnames(W) %in% Dname), drop=FALSE]
    mtZ <- stats::delete.response(stats::terms(formula, data = data, rhs = 2))
    Z <- Matrix::Matrix(Matrix::sparse.model.matrix(mtZ, mf))
    ## Remove intercept
    Z <- Matrix::Matrix(Z[, !colnames(Z) %in% colnames(W), drop=FALSE])

    ## qrX and qrW will typically be sparse, so we will work with those
    r <- remove_collinear(W, Z, tol)
    if (length(r$drp) > 0) {
        Y <- Y[-r$drp]
        D <- D[-r$drp]
        Z <- Z[-r$drp, ]
        W <- W[-r$drp, ]
    }
    W <- W[, !(colnames(W) %in% r$col_idx), drop=FALSE]
    Z <- Z[, !(colnames(Z) %in% r$col_idx), drop=FALSE]

    ret <- ujive.fit(Y, D, r$qrX, r$qrW, r$dX, tol)
    structure(list(IVData=list(Y=Y, D=D, Z=Z, W=W, F=ret$F, k=ret$k, n=ret$n,
                               l=ret$l), call=cl,
                   drop_obs = r$drp, drop_col=r$col_idx,
                   estimate=ret$r),
              class="IVResults")
}

remove_collinear <- function(W, Z, tol) {
    W1 <- W
    Z1 <- Z
    rownames(W1) <- rownames(Z1) <- seq_len(NROW(W1))
    drp <- col_idx <- vector()
    Xdrp <- 0L
    while (length(Xdrp)>0) {
        rs <- drp_singleton(W1, Z1, tol)
        drp <- c(drp, rs$drp)
        col_idx <- c(col_idx, rs$col_idx)
        rc <- drp_collinear(rs$W, rs$Z, tol)
        drp <- c(drp, rc$drp)
        col_idx <- c(col_idx, rc$col_idx)
        ## Leverages
        dX <- unname(Matrix::rowSums(Matrix::qr.Q(rc$qrX)^2))
        Xdrp <- which(dX>1-1e-8)
        if (length(Xdrp)>0) {
            drp <- c(drp, rownames(rc$W[Xdrp, ]))
            W1 <- rc$W[-Xdrp, , drop=FALSE]
            ## Z matrix before dropping collinear cols
            Z1 <- rc$Z[-Xdrp, , drop=FALSE]
            message("Dropping ", length(Xdrp), " obs with leverage 1.")
        } else {
            ## Commit the collinear controls
            if (length(rc$coll_idx)>0) {
                message("Dropping ", length(rc$coll_idx), " collinear IVs.")
                col_idx <- c(col_idx, rc$coll_idx)
            }
        }
    }

    list(qrX=rc$qrX, qrW=rc$qrW, dX=dX, drp=as.numeric(drp), col_idx=col_idx)
}

## Drop singletons
drp_singleton <- function(W, Z, tol) {
    idxW1 <- which(Matrix::colSums(W!=0)==1)
    idxZ1 <- which(Matrix::colSums(Z!=0)==1)
    drpidx <- vector()
    while (length(c(idxW1, idxZ1)) > 0) {
        W1drp <- which(Matrix::rowSums(W[, idxW1, drop=FALSE])>0)
        Z1drp <- which(Matrix::rowSums(Z[, idxZ1, drop=FALSE])>0)
        drp <- c(Z1drp, W1drp)
        stopifnot(length(idxW1)==length(W1drp) && length(idxZ1)==length(Z1drp))
        W <- W[-drp, ]
        Z <- Z[-drp, ]
        idxW1 <- which(Matrix::colSums(W!=0)==1)
        idxZ1 <- which(Matrix::colSums(Z!=0)==1)
        drpidx <- c(drpidx, drp)
    }
    if (length(drpidx)>0)
        message("Recursively dropping ", length(unique(names(drpidx))),
                " obs with singleton covariates or IVs")
    idxW0 <- which(Matrix::colSums(W!=0)==0)
    idxZ0 <- which(Matrix::colSums(Z!=0)==0)
    if (length(idxW0) > 0) W <- W[, -idxW0, drop=FALSE]
    if (length(idxZ0) > 0) Z <- Z[, -idxZ0, drop=FALSE]

    list(W=W, Z=Z, drp=unique(names(drpidx)), col_idx=names(c(idxW0, idxZ0)))
}

## Drop collinear columns
drp_collinear <- function(W, Z, tol) {
    ## Suppress warning matrix rank deficient
    qrW <- suppressWarnings(Matrix::qr(W))
    drpidx <- drp_col <- coll_idx <- vector()
    if (length(idxW <- drp(qrW, tol))>0) {
        message("Dropping ", length(idxW), " collinear controls.")
        coll_idx <- colnames(W)[idxW]
        W <- W[, -idxW, drop=FALSE]
        qrW <- Matrix::qr(W) # TODO: use qrdelete analog
    }

    ## Suppress warning matrix rank deficient
    suppressWarnings(qrX <- Matrix::qr(cbind(W, Z)))
    if (length(drp(qrX, tol))>0) {
        ## If not sparse QR, convert to standard matrix; suppress warnings about
        ## coercion to dense
        ZZ <- if (is.qr(qrW)) as.matrix(Z) else Z
        suppressWarnings(Zt <- Matrix::qr.resid(qrW, ZZ))
        idx0 <- which(Matrix::colSums(abs(Zt)) <= tol)
        if (length(idx0>0)) {
            Zdrop <- which(Matrix::rowSums(Z[, idx0])>0)
            message("Dropping ", length(Zdrop), " obs with ",
                    "no variation in IVs conditional on controls.")
            drp_col <- c(drp_col, colnames(Z)[idx0])
            Z <- Z[-Zdrop, -idx0, drop=FALSE]
            Zt <- Zt[-Zdrop, -idx0, drop=FALSE]
            W <- W[-Zdrop, , drop=FALSE]
            idxW <- which(Matrix::colSums(abs(W))==0)
            drp_col <- c(drp_col, colnames(W)[idxW])
            W <- W[, -idxW]
            qrW <- Matrix::qr(W)
            drpidx <- c(drpidx, names(Zdrop))
        }
        qrZZ <- Matrix::qr(Matrix::crossprod(Zt))
        idxZ <- drp(qrZZ, tol)
        ## Don't update Z, only X
        if (length(idxZ)>0) coll_idx <- colnames(Z)[idxZ]
        ## TODO: use qrdelete analog
        qrX <- Matrix::qr(cbind(W, Z[, -idxZ, drop=FALSE]))
        stopifnot(length(drp(qrX, tol))==0)
    }
    list(qrW=qrW, qrX=qrX, W=W, drp=drpidx, col_idx=drp_col,
         coll_idx=coll_idx, Z=Z)
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
