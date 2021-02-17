context("Group instruments calculations")

test_that("H matrix with group instruments", {

    deltan <- function(ns, pis=rep(0, length(ns)-1)) {
        groups <- rep(seq_along(ns), ns)
        ## Z^*
        Zst <- matrix(0, nrow=length(groups), ncol=length(unique(groups)))
        Zst[cbind(seq_along(groups), groups)] <- 1
        ## Drop first column
        Zst <- Zst[, -1]

        W <- matrix(rep_len(1, length.out=nrow(Zst)), ncol=1)
        Ztil <- lm.fit(y=Zst, x=W)$residuals # \tilde{Z}

        ## Diagnonal of H_W is 1/n
        Hw <- diag(W%*% solve(crossprod(W), t(W)))
        ## Diagonal of H_Z is
        Hz <- diag(Ztil%*% solve(crossprod(Ztil), t(Ztil)))
        ## diag(H) from paper
        h <- Hz-ncol(Ztil)/(nrow(Ztil)-1-ncol(Ztil))*(1-Hw-Hz)
        delta <- sum(h^2)/ncol(Ztil)
        Xi <- drop(crossprod(pis, crossprod(Ztil) %*% pis)) / nrow(Ztil)
        mu <- drop(crossprod(drop(Ztil %*% pis), h)/sqrt(ncol(Ztil)*nrow(Ztil)))

        list(delta=delta, Xi=Xi, mu=mu)
    }

    ## vector indicating group membership
    ns <- c(4, 5, 6, 7)
    ## Expression for delta in supplement
    n <- sum(ns)
    k <- length(ns)-1
    a.delta <- (n-1)^2/(n-1-k)^2*(sum(1/ns)-(k+1)^2/n)/k

    expect_equal(deltan(ns)$delta, a.delta)

    ## Example for R2 (page 8 of reply)
    ns <- 3*c(1, 1, 1, 1, 2, 2)
    k <- length(ns)-1
    n <- sum(ns)
    r <- deltan(ns, c(0, 0, 1, 1, 1))
    b.delta <- 1/9*length(ns)^2/(sum(ns)*(length(ns)-1) *
                                 (1-(length(ns)-1)/(sum(ns)-1))^2)
    b.Xi <- 15/64
    expect_equal(r$delta, b.delta)
    expect_equal(r$Xi, b.Xi)
    expect_equal(r$mu, -sqrt(r$delta*r$Xi)*sqrt(3/5))
    r2 <- deltan(c(rep(3, 8), rep(6, 4)),
                 c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1))
    expect_equal(r2$mu^2/(r2$delta*r2$Xi), 3/5)
})
