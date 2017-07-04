context("Minimum distance calculations")

test_that("Psi3 and Psi4 calculations", {
    r2 <- IVreg(lwage~education+as.factor(yob)|as.factor(qob)*as.factor(yob),
                data=ak80[1:3000, ])
    Y <- r2$IVData$Y[, 1]
    X <- r2$IVData$Y[, 2]
    Z <- r2$IVData$Z
    W <- r2$IVData$W
    d <- IVData(Y, X, Z, W, moments=TRUE)


    ols <- function(X, Y)
        Matrix::solve(Matrix::crossprod(X), Matrix::crossprod(X, Y))
    XX <- cbind(Z, W)
    hatPi <- ols(XX, d$Y)
    V <- d$Y-XX %*% hatPi
    M <- Matrix::diag(d$n) -
        XX %*% Matrix::solve(Matrix::crossprod(XX), Matrix::t(XX))
    m3 <- sum(M^3)
    m4 <- sum(M^4)
    m2 <- sum(Matrix::diag(M)^2)

    vi <- function(i) kronecker(V[i, ]%o%V[i, ], V[i, ])
    Psi3 <- matrix(rowSums(sapply(1:d$n, vi)), ncol=2) / m3

    expect_equal(Psi3, d$Psi3)

    vi <- function(i) kronecker(V[i, ]%o%V[i, ], V[i, ]%o% V[i, ])
    Psi4 <- (m2-m4)*(2*N2 %*% kronecker(d$S, d$S)+tcrossprod(as.vector(d$S)))
    Psi4 <- (matrix(rowSums(sapply(1:d$n, vi)), ncol=4) - Psi4)/ m4

    expect_equal(Psi4, d$Psi4)
    expect_equal(m3, d$m3)
    expect_equal(m4, d$m4)
    expect_equal(m2, d$m2)

})

test_that("LIML, EMD, UMD standard errors are correct", {
    r2 <- IVreg(lwage~education+as.factor(yob)|as.factor(qob)*as.factor(yob),
                data=ak80[1:10000, ], inference="md")
    d <- r2$IVData
    XX <- cbind(d$Z, d$W)
    V <- d$Y-XX %*% Matrix::solve(Matrix::crossprod(XX),
                                  Matrix::crossprod(XX, d$Y))
    r1 <- IVregRE.fit(d)
    tau <- d$k/d$n * (1-d$l/d$n) / (1-d$k/d$n-d$l/d$n)

    ## Try big delta and mu so it matters
    ## LIML
    VV <- function(be, Om, lam, what="liml") {
        VlimlN <- bob(be, Om) * aoa(be, Om) / lam *
            (1+tau/lam)
        ep <- V[, 1]-V[, 2]*be
        gam <- (Om[1, 2]-Om[2, 2]*be)/bob(be, Om)
        v3 <- sum(ep^2*(V[, 2]-gam*ep))/d$m3
        Xi22 <- lam / aoa(be, Om)
        bb <- c(1, -be)
        gg <- c(-gam, 1+gam*be)
        expect_equal(drop(kronecker(bb,  bb) %*% d$Psi3 %*% gg), v3)
        v4 <- drop(kronecker(bb,  bb) %*% d$Psi4 %*% kronecker(gg, gg))

        Vliml <- VlimlN+2*sqrt(d$k/d$n)* d$mu * v3 /Xi22^2 +
            (d$k/d$n) * d$delta* (v4-det(Om))/Xi22^2
        if (what=="liml") return(Vliml)

        ## EMD
        rr <- sqrt(d$k/d$n)*d$mu*sum(ep^3)/d$m3 + (d$k/d$n)*d$delta*
            drop(kronecker(bb,  bb) %*% d$Psi4 %*% kronecker(bb, gg))
        kap <- drop(crossprod(kronecker(bb, bb),
                              d$Psi4 %*% kronecker(bb, bb))) / bob(be, Om)^2 - 3
        ss <- 2*tau+(d$k/d$n)*d$delta*kap
        Vemd <- Vliml - rr^2/ss / (Xi22^2*bob(be, Om)^2)
        if (what=="emd") return(Vemd)

        Vumd <- Vemd + (ss*gam*bob(be, Om)^2+rr)^2/
              (bob(be, Om)^2*ss*Xi22^2)
        if (what=="umd") return(Vumd)

        ## Invalid mbtsls
        Vumd + det(d$T- d$k/d$n * d$S) * Om[2, 2] / Xi22^3 +
            2*sqrt(d$k/d$n)* (d$mu1-d$mu*be) * sum(ep*V[, 2]^2)/d$m3 / Xi22^2
    }

    expect_equal(VV(r1$beta, r1$Om, r1$lam, "liml"),
                 IVregMD.fit(d, weight="LIML")$se^2*d$n)
    ## EMD: standard errors are based on RE starting values
    expect_equal(VV(r1$beta, r1$Om, r1$lam, "emd"),
                 IVregMD.fit(d)$se^2*d$n)
    ## UMD, valid and invalid
    umd <- IVregUMD.fit(d, invalid=FALSE)
    expect_equal(VV(umd$beta, umd$Om, umd$lam, "umd"), umd$se^2*d$n)
    expect_equal(VV(umd$beta, umd$Om, umd$lam, "umdi"),
                 IVregUMD.fit(d, invalid=TRUE)$se^2*d$n)

    ## Try large values so there is a big difference
    d$delta <- 0.1
    d$mu <- -0.2

    expect_equal(VV(r1$beta, r1$Om, r1$lam, "liml"),
                 IVregMD.fit(d, weight="LIML")$se^2*d$n)
    expect_equal(VV(r1$beta, r1$Om, r1$lam, "emd"),
                 IVregMD.fit(d)$se^2*d$n)
    umd <- IVregUMD.fit(d, invalid=FALSE)
    expect_equal(VV(umd$beta, umd$Om, umd$lam, "umd"), umd$se^2*d$n)
    expect_equal(VV(umd$beta, umd$Om, umd$lam, "umdi"),
                 IVregUMD.fit(d, invalid=TRUE)$se^2*d$n)
})

test_that("RE weight matrix is optimal under normality", {

    d <- IVreg(lwage~education+as.factor(yob)|as.factor(qob)*as.factor(yob),
               data=ak80[1:6000, ], inference="md")$IVData

    r1 <- IVregRE.fit(d)
    Xi22 <- r1$lam / aoa(r1$be, r1$Om)
    Del <- MDDelta(r1$beta, r1$Om, Xi22, d, Gaussian=TRUE)
    tau <- d$k/d$n * (1-d$l/d$n) / (1-d$k/d$n-d$l/d$n)

    tt <- rnorm(1)
    m <- c(r1$be, 1) * Xi22^(1/2)
    mm <- outer(m, m)
    Phi <- kronecker(r1$Om, r1$Om) + tt*kronecker(r1$Om, mm) +
        tt*kronecker(mm, r1$Om)
    G <- -L2 %*%  cbind(Xi22^(1/2)* (kronecker(m, c(1, 0)) +
                                      kronecker(c(1, 0), m)),
                         kronecker(m, m)/Xi22)
    moe <- drop(crossprod(m, solve(r1$Om, c(1, 0))))
    Ct <- 2*tau*cbind(c((1+r1$lam/tau)/(1+tt*r1$lam),
                        (2*Xi22^(3/2)*moe) / (1+tt*r1$lam) *
                         (1/tau-tt*(1+2*r1$lam/tau)/(1+2*tt*r1$lam))),
                      c(0, (1+2*r1$lam/tau)/(1+2*tt*r1$lam)))

    expect_equal(solve(Del, G) %*% Ct, crossprod(D2, solve(Phi, D2)) %*% G)
})
