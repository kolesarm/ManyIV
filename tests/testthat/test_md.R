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
