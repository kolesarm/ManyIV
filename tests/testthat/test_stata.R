context("Match Stata output")

test_that("Match stata on AK data", {
    ## This takes 9.8s on XPS13, while Stata (ak.do) takes several mins

    ## Table II
    r0 <- IVreg(lwage~education | qob=="Q1", data=ak80, inference="standard",
                subset=division=="Pacific")
    r1 <- IVreg(ak80$lwage~ak80$education | ak80$qob=="Q1",
                inference="standard", subset=ak80$division=="Pacific")
    expect_identical(r1$estimate, r0$estimate)
    ## LIML vs TSLS
    expect_equal(unlist(r0$estimate[2, 1:3]),
                 unlist(r0$estimate[3, 1:3]))
    ## TSLS with THE
    expect_equal(r0$estimate[2, 3],
                 r0$estimate[2, 4])
    ## Stata
    expect_equal(unname(unlist(r0$estimate[2, 1:2])),
                 c(0.1289102761, 0.0720969950))

    ## Table V in AK
    r1 <- IVreg(lwage~education+as.factor(yob) | qob*as.factor(yob),
                data=ak80)
    expect_equal(as.numeric(r1$estimate["ols", 1:3]),
                 c(0.0710810458, 0.0003390067, 0.0003814625))
    expect_equal(as.numeric(r1$estimate["tsls", 1:3]),
                 c(0.0891154613, 0.0161098202, 0.0162120317))
    expect_equal(as.numeric(r1$estimate["mbtsls", 1:3]),
                 c(0.0937333665, 0.0180984698, 0.0204147326))
    expect_equal(as.numeric(r1$estimate["liml", 1:3]),
                 c(0.0928764165, 0.0177441446, 0.0196323640))

    r2 <- IVreg(lwage~education+as.factor(yob) | qob,
                data=ak80[1:100, ])
    expect_equal(as.numeric(r2$estimate["ols", 1:3]),
                 c(0.1004438463, 0.0346482285, 0.0233448020))
    expect_equal(as.numeric(r2$estimate["tsls", 1:3]),
                 c(-0.1609028242, 0.3277804536, 0.4120184496))
    ## For mbtsls, covariance matrix is not pd, so we report Inf, unlike stata
    expect_equal(r2$estimate["mbtsls", 1], c(0.3505927902))
    expect_equal(as.numeric(r2$estimate["liml", 1:3]),
                 c(-0.2054464334, 0.3788818685, 0.49813860131))

    ## Small-sample adjustmenst are not done by default in Stata in IV. One
    ## needs option "small". Then we scale by n/(n-L-1)
    expect_equal(as.numeric(r2$estimate["tsls", 2:3])*sqrt(100/89),
                 c(0.3474465859, 0.4367386830))
    expect_equal(unname(unlist(r2$estimate[3, 2:3]))*sqrt(100/89),
                 c(0.4016139774, 0.5280258614))

    ## Table VII
    r3 <- IVreg(lwage~education+as.factor(yob)+sob |
                    qob*sob + qob*as.factor(yob), data=ak80)
    r0 <- capture.output(print(r3, digits=4))

    expect_equal(as.numeric(r3$estimate["tsls", 1:3]),
                 c(0.0928180625, 0.0093013345, 0.0096641481))
    expect_equal(as.numeric(r3$estimate["liml", 1:3]),
                 c(0.1063979834, 0.0116383738, 0.0149803714))
    expect_equal(as.numeric(r3$estimate["ols", 1:3]),
                 c(0.0673389705, 0.0003464258, 0.0003883466))
    expect_equal(as.numeric(r3$estimate["mbtsls", 1:3]),
                 c(0.1089429483, 0.0120411994, 0.0159979308))
    expect_equal(r0[length(r0)-1],
                 "liml    0.10640    0.0116384      0.0149804         NA")

})

## Stata file:

## Store as Stata to compare results
## foreign::write.dta(ak80, "ak80.dta")

## Then run:

## timer on 1
## use ak80.dta, replace

## generate Q1 = (qob==1)

## qui: eststo r0: ivregress 2sls lwage (education = Q1) if division==9
## esttab r0, b(%11.10f) se(%11.10f) drop(*_cons*)
## /* 0.1289102761 (0.0720969950) */


## /* Table V Columns (1) and (2) */
## qui: eststo r1: reg  lwage education i.yob
## qui: eststo r2: reg  lwage education i.yob, robust
## qui: eststo r3: reg  lwage education i.yob in 1/100
## qui: eststo r4: reg  lwage education i.yob in 1/100, robust

## esttab r1 r2 r3 r4, b(%11.10f) se(%11.10f) drop(*yob* *cons*)
## /* 0.0710810458 (0.0003390067) [0.0003814625] */
## /* 0.1004438463 (0.0346482285) [0.0233448020] */

## qui: eststo i1: ivregress 2sls lwage (educ=i.yob##i.qob) i.yob
## qui: eststo i2: ivregress 2sls lwage (educ=i.yob##i.qob) i.yob, robust
## qui: eststo i3: ivregress 2sls lwage (educ=i.qob) i.yob in 1/100
## qui: eststo i4: ivregress 2sls lwage (educ=i.qob) i.yob in 1/100, robust
## qui: eststo i5: ivregress 2sls lwage (educ=i.qob) i.yob in 1/100, small
## qui: eststo i6: ivregress 2sls lwage (educ=i.qob) i.yob in 1/100, r small

## esttab i1 i2 i3 i4 i5 i6, b(%11.10f) se(%11.10f) drop(*yob* *cons*)
## /* 0.0891154613 (0.0161098202) [0.0162120317] */
## /*-0.1609028242 (0.3277804536) [0.4120184496]
## /* small: (0.3474465859) [0.4367386830] */

## /* MBTSLS */
## /* first get number of exog and endog variagbles */
## qui: ivreg2 lwage (education=i.qob) i.yob in 1/100
## local k1 = 1+e(exexog_ct)/(e(N)-e(exexog_ct)-1-e(inexog_ct))
## qui: eststo i1: ivreg2 lwage (education=i.qob) i.yob in 1/100, kclass(`k1')
## qui: eststo i2: ivreg2 lwage (education=i.qob) i.yob in 1/100, kclass(`k1') r
## qui: ivreg2 lwage (education=i.yob##i.qob) i.yob
## local k2 = 1+e(exexog_ct)/(e(N)-e(exexog_ct)-1-e(inexog_ct))
## qui: eststo i3: ivreg2 lwage (education=i.yob##i.qob) i.yob, kclass(`k2')
## qui: eststo i4: ivreg2 lwage (education=i.yob##i.qob) i.yob, kclass(`k2') r

## esttab i1 i2 i3 i4, b(%11.10f) se(%11.10f) drop(*yob* *cons*)
## /* 0.3505927902 (0.0504391651) [0.3894877272] */
## /* 0.0937333665 (0.0180984698) [0.0204147326] */

## /* LIML*/
## qui: eststo i1: ivregress liml lwage  (educ=i.yob##i.qob) i.yob
## qui: eststo i2: ivregress liml lwage  (educ=i.yob##i.qob) i.yob, robust
## qui: eststo i3: ivregress liml lwage  (educ=i.qob) i.yob in 1/100
## qui: eststo i4: ivregress liml lwage  (educ=i.qob) i.yob in 1/100, r
## qui: eststo i5: ivregress liml lwage  (educ=i.qob) i.yob in 1/100, small
## qui: eststo i6: ivregress liml lwage  (educ=i.qob) i.yob in 1/100, r small

## esttab i1 i2 i3 i4 i5 i6, b(%11.10f) se(%11.10f) drop(*yob* *cons*)
## /* 0.0928764165 (0.0177441446) [0.0196323640] */
## /*-0.2054464334 (0.3788818685) [0.4981386013] */
## /* small: (0.4016139774) [0.5280258614]*/


## /* Table VII */

## qui: eststo i1: reg  lwage education i.yob i.sob
## qui: eststo i2: reg  lwage education i.yob i.sob, robust
## qui: eststo i3: ivregress 2sls lwage ///
## (education=i.yob##i.qob i.qob##i.sob) i.sob i.yob
## qui: eststo i4: ivregress 2sls lwage (educa=i.yob##i.qob i.qob##i.sob) ///
## i.sob i.yob, robust
## qui: eststo i5: ivregress liml lwage (educa=i.yob##i.qob i.qob##i.sob) ///
## i.sob i.yob
## qui: eststo i6: ivregress liml lwage (educa=i.yob##i.qob i.qob##i.sob) ///
## i.sob i.yob, robust
## esttab i1 i2 i3 i4 i5 i6, b(%11.10f) se(%11.10f) drop(*yob* *sob*)
## /* 0.0673389705 (0.0003464258) [0.0003883466] */
## /* 0.0928180625 (0.0093013345) [0.0096641481] */
## /* 0.1063979834 (0.0116383738) [0.0149803714] */

## /* MBTSLS */
## /* first get number of exog and endog variagbles */
## qui: ivreg2 lwage (education=i.yob##i.qob i.qob##i.sob) i.sob i.yob
## local k3 =  1+e(exexog_ct)/(e(N)-e(exexog_ct)-1-e(inexog_ct))
## qui: eststo i1: ivreg2 lwage (education=i.yob##i.qob i.qob##i.sob) ///
## i.sob i.yob, kclass(`k3')
## qui: eststo i2: ivreg2 lwage (education=i.yob##i.qob i.qob##i.sob) ///
## i.sob i.yob, kclass(`k3') robust
## esttab i1 i2, b(%11.10f) se(%11.10f) drop(*yob* *sob*)
## /* 0.1089429483 (0.0120411994) [0.0159979308] */

## timer off 1
## timer list
