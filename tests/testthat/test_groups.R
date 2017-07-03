context("Group instruments calculations")

test_that("H matrix with group instruments", {

    ## vector indicating group membership
    ns <- c(4, 5, 6, 7)
    groups <- rep(1:4, ns)
    ## Z^*
    Zst <- matrix(0, nrow=length(groups), ncol=length(unique(groups)))
    Zst[cbind(seq_along(groups), groups)] <- 1
    ## Drop first column
    Zst <- Zst[, -1]
    W <- matrix(rep_len(1, length.out=nrow(Zst)), ncol=1)
    Ztil <- lm.fit(y=Zst, x=W)$residuals

    ## Diagnonal of H_W is 1/n
    Hw <- diag(W%*% solve(crossprod(W), t(W)))
    ## Diagonal of H_Z is
    Hz <- diag(Ztil%*% solve(crossprod(Ztil), t(Ztil)))
    ## diag(H) from paper
    h <- Hz-ncol(Ztil)/(nrow(Ztil)-1-ncol(Ztil))*(1-Hw-Hz)
    delta <- sum(h^2)/ncol(Ztil)

    ## Closed-form solution
    a.Hz <- diag(rep(1/ns, ns)) %*% outer(groups, groups, `==`) - 1/nrow(Ztil)
    Hmean <- (ncol(Ztil)+1)/(sum(1/ns))
    a.delta <- (nrow(Ztil)-1)^2*(ncol(Ztil)+1)/(nrow(Ztil)-1-ncol(Ztil))^2 *
        (1/Hmean-(ncol(Ztil)+1)/nrow(Ztil))/ncol(Ztil)

    expect_equal(delta, a.delta)
})
