context("m3 and m4 approximations")

if (FALSE) {
inference <- c("standard", "re", "il", "lil", "md")
r2 <- IVreg(lwage~education+as.factor(yob)|as.factor(qob)*as.factor(yob),
            data=ak80[1:10000, ], inference=inference, approx=FALSE)
r2a <- IVreg(lwage~education+as.factor(yob)|as.factor(qob)*as.factor(yob),
            data=ak80[1:10000, ], inference=inference, approx=TRUE)
## 9880 vs 9880.32
## 9840 vs 9840.97

r3a <- IVreg(lwage~education+as.factor(yob)+as.factor(sob)|
                as.factor(qob)*as.factor(sob)+as.factor(qob)*as.factor(yob),
            data=ak80, inference=inference, approx=TRUE)
## m3: 328789 vs 328790
## m4: 328549 vs 328553

## r1 <- IVreg(lwage~education+as.factor(yob)|as.factor(qob)*as.factor(yob),
##             data=ak80, inference=inference, approx=FALSE)
r1a <- IVreg(lwage~education+as.factor(yob)|as.factor(qob)*as.factor(yob),
            data=ak80, inference=inference, approx=TRUE)
## m3: 329389 vs 329389.009
## m4: 329349 vs 329349.029
}
