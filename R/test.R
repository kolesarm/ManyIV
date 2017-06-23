
## ## vector indicating group membership
## ns <- c(4, 5, 6, 7)
## groups <- rep(1:4, ns)
## n <- length(groups)
## k <- length(unique(groups))-1
## ## turn into a matrix of group indicators
## Zst <- matrix(0, nrow=n, ncol=k+1)
## Zst[cbind(seq_along(groups), groups)] <- 1
## ## Drop first column
## Zst <- Zst[, -1]

## W <- matrix(rep_len(1, length.out=nrow(Zst)), ncol=1)
## Z <- lm.fit(y=Zst, x=W)$residuals

## ## Diagnonal of H_W is 1/n
## Hw <- diag(W%*% solve(crossprod(W), t(W)))
## ## Diagonal of H_Z is
## Hz <- diag(Z%*% solve(crossprod(Z), t(Z)))
## ## diag(H) from paper
## h <- Hz-k/(n-1-k) *(1-Hw-Hz)
## delta <- sum(h^2)/k

## ## Closed-form solution
## a.Hz <- diag(rep(1/ns, ns)) %*% outer(groups, groups, `==`) - 1/n
## Hmean <- (k+1)/(sum(1/ns))
## a.delta <- (n-1)^2*(k+1)/(n-1-k)^2*(1/Hmean-(k+1)/n)/k
