context("Match Stata output")

test_that("Match stata on AK data", {
    ## This takes 9.8s on XPS13, while Stata (ak.do) takes several mins

    ## Table V in AK
    r1 <- IVreg(lwage~education+as.factor(yob)|as.factor(qob)*as.factor(yob),
                data=ak80)
    expect_equal(as.numeric(r1$estimate["tsls", 1:3]),
                 c(0.0891154613, 0.0161098202, 0.0162120317))
    expect_equal(as.numeric(r1$estimate["liml", 1:3]),
                 c(0.0928764165, 0.0177441446, 0.0196323640))
    expect_equal(as.numeric(r1$estimate["ols", 1:3]),
                 c(0.0710810458, 0.0003390067, 0.0003814625))
    expect_equal(as.numeric(r1$estimate["mbtsls", 1:3]),
                 c(0.0937333665, 0.0180984698, 0.0204147326))

    r2 <- IVreg(lwage~education+as.factor(yob)|as.factor(qob),
                data=ak80[1:100, ])
    expect_equal(as.numeric(r2$estimate["tsls", 1:3]),
                 c(-0.1609028242, 0.3277804536, 0.4120184496))
    expect_equal(as.numeric(r2$estimate["liml", 1:3]),
                 c(-0.2054464334, 0.3788818685, 0.49813860131))
    expect_equal(as.numeric(r2$estimate["ols", 1:3]),
                 c(0.1004438463, 0.0346482285, 0.0233448020))
    ## For mbtsls, covariance matrix is not pd, so we report Inf, unlike stata
    expect_equal(r2$estimate["mbtsls", 1], c(0.3505927902))
    ## Small-sample adjustmenst are not done by default in Stata in IV. One
    ## needs option "small". Then we scale by n/(n-L-1)
    expect_equal(unname(unlist(r2$estimate[3, 2:3]))*sqrt(100/89),
                 c(0.4016139774, 0.5280258614))

    ## Table VII
    r3 <- IVreg(lwage~education+as.factor(yob)+as.factor(sob)|
                    as.factor(qob)*as.factor(sob)+
                    as.factor(qob)*as.factor(yob), data=ak80)
    r0 <- capture.output(print(r3, digits=4))


    expect_equal(as.numeric(r3$estimate["tsls", 1:3]),
                 c(0.0928180625, 0.0093013345, 0.0096641481))
    expect_equal(as.numeric(r3$estimate["liml", 1:3]),
                 c(0.1063979834, 0.0116383738, 0.0149803714))
    expect_equal(as.numeric(r3$estimate["ols", 1:3]),
                 c(0.0673389705, 0.0003464258, 0.0003883466))
    expect_equal(as.numeric(r3$estimate["mbtsls", 1:3]),
                 c(0.1089429483, 0.0120411994, 0.0159979308))
    expect_equal(r0[12],
                 "liml    0.10640    0.0116384      0.0149804         NA")

    ## Table II
    ak80$Q1 <- (ak80$qob==1)+0
    r4 <- IVreg(lwage~education | 0+Q1, data=ak80)
    ## LIML vs TSLS
    expect_equal(unlist(r4$estimate[2, 1:3]),
                 unlist(r4$estimate[3, 1:3]))
    ## TSLS with THE
    expect_equal(r4$estimate[2, 3],
                 r4$estimate[2, 4])

})
