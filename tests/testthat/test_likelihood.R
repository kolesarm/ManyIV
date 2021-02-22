context("Likelihood-based inference")

test_that("Test random effects likelihood", {
    r1 <- IVreg(lwage~education+as.factor(yob) | qob*as.factor(yob),
                data=ak80, inference=c("re", "lil", "il"))

    expect_equal(unname(unlist(r1$estimate["liml", ])),
                 c(0.09287641, 0.01615828, 0.01986003, 0.01978592))

})

test_that("Test Invariant likelihood", {
    ## Old cover for invariant likelihood
    r1 <- IVreg(lwage~education+as.factor(yob)|as.factor(qob)*as.factor(yob),
                data=ak80)
    r3 <- IVreg(lwage~education+as.factor(yob)+as.factor(sob)|
                    as.factor(qob)*as.factor(sob)+
                    as.factor(qob)*as.factor(yob), data=ak80)
    IVregILold.fit <- function(d) {
        lG <- function(t, k)
        log(besselI(t, k/2-1, expon.scaled=TRUE))+t-(k/2-1)*log(t/2)

        ## minus log-likelihood

        logl <- function(be, lam, Om) {
            Qt <- drop(crossprod(solve(Om, c(be, 1)),
                                 d$T %*% solve(Om, c(be, 1))) / aoa(be, Om))
            ((d$n-d$l)*log(det(Om)) + sum(diag(solve(Om, d$nu*d$S+d$n*d$T)))
                + d$n*lam - 2*lG(d$n*sqrt(lam*Qt), d$k)) / 2
        }

        ff <- function(be, t)
            logl(be, exp(t[1]), cbind(c(exp(t[2]), t[3]), c(t[3],
            (exp(t[4])+t[3]^2)/exp(t[2]))))
        r1 <- IVregRE.fit(d)
        start <- c(log(r1$lam), log(d$S[1, 1]), d$S[1, 2], log(det(d$S)))
        il <- stats::optim(start, function(t) ff(r1$beta, t))
        se.il <- sqrt(solve(numDeriv::hessian(function(t) ff(t[1], t[2:5]),
                                              c(r1$beta, il$par)))[1, 1])
        Om <- cbind(c(exp(il$par[2]), il$par[3]),
                    c(il$par[3], (exp(il$par[4])+il$par[3]^2)/exp(il$par[2])))
        list(beta=r1$beta, se=se.il, lam=exp(il$par[1]), Om=Om,
             lik=-logl(r1$beta, exp(il$par[1]), Om))
    }

    expect_true(IVregIL.fit(r1$IVData)$lik-IVregILold.fit(r1$IVData)$lik>0)
    expect_true(IVregIL.fit(r3$IVData)$lik-IVregILold.fit(r3$IVData)$lik>0)

})

test_that("Overid test", {
    r1 <- IVreg(lwage~education | qob=="Q1", data=ak80, inference="il")
    expect_equal(IVoverid(r1),
                 data.frame(statistic=c("Sargan"=NA, "Modified-CD"=NA),
                            p.value=as.numeric(c(NA, NA))))
    r2 <- IVreg(lwage~education | qob, data=ak80, inference="il")
    expect_equal(IVoverid(r2)$p.value, c(0.24119662, 0.241200997))
})
