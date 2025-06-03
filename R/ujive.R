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
#' @return An object of class \code{"IVResults"}.
#' @export
ujive <- function(formula, data, subset, na.action) {
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

    d <- IVData(Y, X, Z, W, moments=FALSE)

    structure(list(IVData=d, call=cl, estimate=as.data.frame(ujive.fit(d))),
              class="IVResults")
}


ujive.fit <- function(d) {
    Z <- d$Z
    W <- d$W
    Y <- d$Y[, 1]
    D <- d$Y[, 2]

    X <- cbind(Z, W)
    X0 <- Matrix::solve(Matrix::crossprod(X), Matrix::t(X)) # (X'X)^{-1}X'
    W0 <- Matrix::solve(Matrix::crossprod(W), Matrix::t(W))
    diagH <- function(x, x0) Matrix::colSums(Matrix::t(x) * x0)
    proj <- function(x, x0, y) Matrix::crossprod(x0, Matrix::crossprod(x, y))
    dX <- diagH(X, X0)
    dW <- diagH(W, W0)
    DtX <- D-proj(X, X0, D)
    DtW <- D-proj(W, W0, D)

    Dtsls <- DtW-DtX
    Dujive <- Dtsls - (dX-dW) / (1-dX)*DtX
    Zjive1 <- D-1 / (1-dX)*DtX
    Djive1 <- Zjive1 - proj(W, W0, Zjive1)
    Dojive <- Zjive1 - D + 1 / (1-dW)*DtW
    Zijive <- DtW - 1 / (1-dX+dW)*DtX
    Dijive <- Zijive - proj(W, W0, Zijive)

    ## TSLS, UJIVE, old UJIVE, IJIVE, JIVE1
    Dlist <- list(DtW, Dtsls, Dujive, Dojive, Dijive, Djive1)
    den <- vapply(Dlist, function(x) sum(x*D), numeric(1))
    num <- vapply(Dlist, function(x) sum(x*Y), numeric(1))
    est <- num/den
    names(est) <- c("ols", "tsls", "ujive", "old ujive", "ijive1", "jive1")

    ## Textbook robust
    YtW <- Y-proj(W, W0, Y)
    epD <- (YtW - DtW %o% est)*do.call(cbind, Dlist) # epsilon*Dhat
    text <- Matrix::colSums(epD^2)
    ## HTE robust
    YtX <- Y-proj(X, X0, Y)
    GYtsls <- YtW-YtX-Dtsls*est[2]
    DYujive <- (dX-dW) / (1-dX) * (Y-D*est[3])
    GYujive <- YtW-YtX-Dtsls*est[3] - DYujive + proj(X, X0, DYujive)
    DYjive1 <- 1 / (1-dX) * (YtW-DtW*est[6])
    GYjive1 <- (YtW-DtW*est[6]) - DYjive1 + proj(X, X0, DYjive1)
    DYujivW <- 1 / (1-dW) * (Y-D*est[4])
    DYujivX <- 1 / (1-dX) * (Y-D*est[4])
    GYojive <- DYujivW-DYujivX-proj(W, W0, DYujivW)+proj(X, X0, DYujivX)
    DYijiv1 <- 1 / (1-dX+dW) * (Y-D*est[5])
    DYijiv2 <- (Y-D*est[5])-DYijiv1+proj(X, X0, DYijiv1)-proj(W, W0, DYijiv1)
    GYijive <- DYijiv2 - proj(W, W0, DYijiv2)

    GY <- cbind(0, GYtsls, GYujive, GYojive, GYijive, GYjive1)
    hte <- Matrix::colSums((GY*DtX+epD)^2)
    r <- cbind(estimate=est,
               se_text=sqrt(text/den^2),
               se_hte=sqrt(hte/den^2))
    r
}


jive.fit <- function(Y, D, Z, W) {
    X <- cbind(Z, W)
    diagH <- function(x) {
        Matrix::colSums(Matrix::t(x) *
                            Matrix::solve(Matrix::crossprod(x), Matrix::t(x)))
    }
    ols <- function(X, Y) {
        Matrix::solve(Matrix::crossprod(X), Matrix::crossprod(X, Y))
    }
    proj <- function(x, y) Matrix::drop(x %*% ols(x, y))

    ## Reduced form coeffs
    DX <- diagH(X)
    DW <- diagH(W)
    Zhat <- (proj(X, D)-DX*D) / (1-DX)
    Dhat <- Zhat - (proj(W, D)-DW*D) / (1-DW)

    Yt <- Y-proj(W, Y)
    Dt <- D-proj(W, D)
    Zt <- Z-proj(W, Z)
    Rt <- Zt %*% ols(X, cbind(Y, D))[seq_len(ncol(Z)), ]

    ## (TSLS, JIVE1, UJIVE)
    be  <- c(sum(Rt[, 1]*Rt[, 2]) / sum(Rt[, 2]^2),
             sum(Zhat*Yt)/sum(Zhat*Dt),
             sum(Dhat*Y)/sum(Dhat*D))
    ## Standard errors
    epsilon <- function(beta) Yt-Dt*beta
    U <- cbind(Y, D)-X %*% ols(X, cbind(Y, D))
    num1 <- vapply(be, function(b) mean(epsilon(b)^2), numeric(1))
    den1 <- c(sum(D*proj(Zt, D)), sum(Zhat*Dt), sum(Dhat*D))
    num2 <- vapply(be, function(b) sum(epsilon(b)^2*Rt[, 2]^2), numeric(1))
    ff <- function(b) sum((epsilon(b)*Rt[, 2] + (Rt[, 1]-Rt[, 2]*b)*U[, 2])^2)
    num3 <- vapply(be, ff, numeric(1))
    ## Many weak term
    A <- (U[, 1]-U[, 2]*be[3])*Zt
    B <- U[, 2]*Zt
    AA <- Matrix::solve(Matrix::crossprod(Zt), Matrix::crossprod(A))
    BB <- Matrix::solve(Matrix::crossprod(Zt), Matrix::crossprod(B))
    CC <- Matrix::solve(Matrix::crossprod(Zt), Matrix::crossprod(A, B))
    w1 <- sum(diag(as.matrix(AA) %*% as.matrix(BB)))
    w2 <- sum(diag(as.matrix(CC) %*% as.matrix(CC)))
    mi <- c(NA, w1+w2, w1+w2)
    estimate <- cbind(be, sqrt(num1/ifelse(den1<0, NA, den1)),
                      sqrt(num2/den1^2), sqrt(num3/den1^2),
                      sqrt(mi+num3)/abs(den1))
    rownames(estimate) <- c("TSLS", "JIVE1", "UJIVEold")
    colnames(estimate) <- c("Estimate", "Homosks", "Het. Robust", "HTE Robust",
                            "MW Robust")
    list(estimate=estimate, r=den1)
}
