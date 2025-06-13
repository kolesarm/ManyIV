context("UJIVE vs old UJIVE code")

ujive2.fit <- function(Y, D, Z, W) {
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


test_that("ujive", {
    t1 <- ujive(lwage~education | qob, data=ak80)
    s1 <- jive.fit(t1$IVData$Y, t1$IVData$D, t1$IVData$Z, t1$IVData$W)
    u1 <- ujive2.fit(t1$IVData$Y, t1$IVData$D, t1$IVData$Z, t1$IVData$W)

    expect_lt(max(abs(u1-t1$estimate)), 1e-7)
    expect_lt(max(abs(t1$estimate[c(2, 6, 4), 1]-s1$estimate[, 1])), 1e-7)
    ## Standard error in the numerator uses TSLS
    expect_lt(max(abs(s1$estimate[1, c(1, 3:4)]-t1$estimate[2, ])), 1e-7)

    t2 <- ujive(lwage~education+as.factor(yob) | qob*as.factor(yob), data=ak80)
    s2 <- jive.fit(t2$IVData$Y, t2$IVData$D, t2$IVData$Z, t2$IVData$W)
    u2 <- ujive2.fit(t2$IVData$Y, t2$IVData$D, t2$IVData$Z, t2$IVData$W)

    expect_lt(max(abs(u2-t2$estimate)), 1e-8)
    expect_lt(max(abs(t2$estimate[c(2, 6, 4), 1]-s2$estimate[, 1])), 1e-9)
    ## Standard error in the numerator uses TSLS
    expect_lt(max(abs(s2$estimate[1, c(1, 3:4)]-t2$estimate[2, ])), 3e-10)

    ## Drop collinear
    expect_message(t2a <- ujive(lwage~education+as.factor(yob) |
                                    0+qob*as.factor(yob), data=ak80))
    expect_lt(max(abs(t2$estimate-t2a$estimate)), 1e-9)
    MM <- model.matrix(~0+qob:as.factor(yob), data=ak80)
    colnames(MM) <- NULL
    expect_message(t2b <- ujive(lwage~education+as.factor(yob) | MM, data=ak80))
    expect_lt(max(abs(t2$estimate-t2b$estimate)), 1e-9)

    expect_message(t4 <- ujive(lwage~education+as.factor(yob) |
                                   qob*as.factor(yob),
                               data=ak80, subset=ak80$sob=="AK"))
    expect_equal(NCOL(t4$IVData$W)+NCOL(t4$IVData$Z)-t4$IVData$k-t4$IVData$l,
                 18L)
    expect_equal(as.numeric(t4$estimate[3, ]),
                 c(0.090491169, 0.211386425, 0.378517855))

    ## Single instrument
    ts <- ujive(lwage~education+married | I(qob=="Q1"), data=ak80)
    tm <- IVreg(lwage~education+married | I(qob=="Q1"), data=ak80)
    expect_equal(tm$IVData$F, ts$IVData$F)
    expect_lt(max(abs(ts$estimate[1:2, 1:2]-tm$estimate[1:2, c(1, 3)])), 1e-6)

    ## FHL data
    fm <- droplevels(fhl[1:1000, ])
    t0 <- ujive(ln_total_patents_appl~dallowed+ind_year | examiner, data=fm)
    expect_lt(max(abs(t0$estimate[3, ]-c(2.452386784, 2.151379876,
                                         3.170662565))), 1e-8)
    t1 <- ujive(ln_total_patents_appl~dallowed+ind_year | examiner,
                data=fm[-t0$drop_obs, ])
    expect_lt(max(abs(t1$estimate - t0$estimate)), 1e-10)

})
