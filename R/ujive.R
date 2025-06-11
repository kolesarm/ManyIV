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
    r <- remove_collinear(Y, D, W, Z, tol)
    ret <- ujive.fit(r$Y, r$D, r$qrX, r$qrW, tol)
    structure(list(IVData=list(Z=r$Z, D=r$D, W=r$W, Y=r$Y, F=ret$F, k=ret$k,
                               n=ret$n, l=ret$l), call=cl,
                   drop_idx = r$ret_idx, estimate=ret$r),
              class="IVResults")
}

remove_collinear <- function(Y, D, W, Z, tol) {
    ## Drop zero columns

   Wsum <- Matrix::colSums(W!=0)
    idxW1 <- which(Wsum==1)
    if (length(idxW1) > 0) {
        W1drop <- which(Matrix::rowSums(W[, idxW1])>0)
        message("Dropping ", length(W1drop), " obs with singleton",
                " covariates")
        stopifnot(length(idxW1)==length(W1drop))
        W <- W[-W1drop, -idxW1, drop=FALSE]
        Z <- Z[-W1drop, , drop=FALSE]
        Y <- Y[-W1drop]
        D <- D[-W1drop]
    }
    idxW0 <- which(Wsum==0)
    if (length(idxW0) > 0) W <- W[, -idxW0]

    ## Suppress warning matrix rank deficient
    qrW <- suppressWarnings(Matrix::qr(W))
    ret_idx <- vector()
    if (length(idxW <- drp(qrW, tol))>0) {
        message("Dropping the following collinear controls: ",
                paste(colnames(W)[idxW], collapse="\n, "))
        ret_idx <- c(ret_idx, colnames(W)[idxW])
        W <- W[, -idxW]
        qrW <- Matrix::qr(W) # TODO: use qrdelete analog
        stopifnot(length(drp(qrW, tol))==0)
    }

    Zsum <- Matrix::colSums(Z!=0)
    idxZ0 <- which(Zsum==0)
    ret_idx <- names(c(idxW0, idxW1, idxZ0))
    if (length(idxZ0) > 0) Z <- Z[, -idxZ0]
    ## Suppress warning matrix rank deficient
    suppressWarnings(qrX <- Matrix::qr(cbind(W, Z)))
    if (length(drp(qrX, tol))>0) {
        ## If not sparse QR, convert to standard matrix
        Zt <- Matrix::qr.resid(qrW, if (is.qr(qrW)) as.matrix(Z) else Z)
        idx0 <- which(Matrix::colSums(abs(Zt)) <= tol)
        if (length(idx0>0)) {
            Zdrop <- which(Matrix::rowSums(Z[, idx0])>0)
            message("Dropping ", length(Zdrop), " obs with ",
                    "no variation in IV conditional on controls")
            ret_idx <- c(ret_idx, colnames(Z)[idx0])
            Z <- Z[-Zdrop, -idx0]
            Zt <- Zt[-Zdrop, -idx0]
            W <- W[-Zdrop, ]
            W <- W[, Matrix::colSums(abs(W))>0]
            qrW <- Matrix::qr(W)
            Y <- Y[-Zdrop]
            D <- D[-Zdrop]
        }
        qrZZ <- Matrix::qr(Matrix::crossprod(Zt))
        idxZ <- drp(qrZZ, tol)
        if (length(idxZ)>0) {
            message("Dropping ", length(idxZ), " collinear instruments.")
            ret_idx <- c(ret_idx, colnames(Z)[idxZ])
            Z <- Z[, -idxZ]
        }
        qrX <- Matrix::qr(cbind(W, Z)) # TODO: use qrdelete analog
        stopifnot(length(drp(qrX, tol))==0)
    }
    list(qrX=qrX, qrW=qrW, Y=Y, D=D, W=W, Z=Z, ret_idx=ret_idx)
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



ujive.fit <- function(Y, D, qrX, qrW, tol) {
    ## Diagonals of projection matrices
    dW <- Matrix::rowSums(Matrix::qr.Q(qrW)^2)
    stopifnot(max(dW)<1-tol)

    dX <- Matrix::rowSums(Matrix::qr.Q(qrX)^2)
    dM <- 1-dX
    MX <- function(A) Matrix::qr.resid(qrX, A)

    ## Treatment residuals
    HD <- Matrix::qr.fitted(qrX, D)-Matrix::qr.fitted(qrW, D)
    DtW <- Matrix::qr.resid(qrW, D)
    MD <- DtW-HD
    ## Zero out cases with leverage one
    idx0 <- dX > 1-tol
    ddM <- function(A) ifelse(idx0, 0L, A/dM)

    Dujive <- HD - MD * (dX-dW)/dM
    Djive1 <- DtW-Matrix::qr.resid(qrW, ddM(MD))
    Dojive <- DtW / (1-dW) - MD/dM
    Dijive <- DtW-Matrix::qr.resid(qrW, MD / (dM+dW))
    Dujive[idx0] <- Dojive[idx0] <- Dijive[idx0] <- 0L

    ## TSLS, UJIVE, old UJIVE, IJIVE, JIVE1
    Dlist <- list(DtW, HD, Dujive, Dojive, Dijive, Djive1)
    den <- vapply(Dlist, function(x) sum(x*D), numeric(1))
    num <- vapply(Dlist, function(x) sum(x*Y), numeric(1))
    est <- num/den
    names(est) <- c("ols", "tsls", "ujive", "old ujive", "ijive1", "jive1")
    if (sum(idx0)>0)
        message("Zeroing out ", sum(idx0), " observations with leverage one.")

    ## Textbook robust
    YtW <- Matrix::qr.resid(qrW, Y)
    epD <- (YtW - DtW %o% est)*do.call(cbind, Dlist) # epsilon*Dhat
    text <- Matrix::colSums(epD^2)
    ## HTE robust
    HY <- Matrix::qr.fitted(qrX, Y)-Matrix::qr.fitted(qrW, Y)
    GYtsls <- HY-HD*est[2]
    GYujive <- HY-HD*est[3] - MX((Y-D*est[3]) * ddM(dX-dW))
    GYjive1 <- YtW-DtW*est[6] - MX(ddM(YtW-DtW*est[6]))
    GYojive <- Matrix::qr.resid(qrW, (Y-D*est[4]) / (1-dW))-
        MX(ddM(Y-D*est[4]))
    GYijive <- YtW-DtW*est[5]-MX((YtW-DtW*est[5]) / (dM+dW))
    GYujive[idx0] <- GYojive[idx0] <- GYijive[idx0] <- 0L

    GY <- cbind(0, GYtsls, GYujive, GYojive, GYijive, GYjive1)
    hte <- Matrix::colSums((GY*MD+epD)^2)
    r <- cbind(estimate=est, se_text=sqrt(text/den^2), se_hte=sqrt(hte/den^2))
    n <- length(Y)
    kl <- if (is.qr(qrX)) qrX$rank else dim(qrX)[2]
    k <- kl - if (is.qr(qrW)) qrW$rank else dim(qrW)[2]
    Fhom <- sum(HD^2)/sum(MD^2) * (n-kl) / k

    list(r=as.data.frame(r), F=Fhom, n=n, k=k, l=kl-k)
}
