## ----include=FALSE, cache=FALSE-----------------------------------------------
library("knitr")
knitr::opts_knit$set(self.contained = FALSE)
knitr::opts_chunk$set(tidy = TRUE, collapse=TRUE, comment = "#>",
                      tidy.opts=list(blank=FALSE, width.cutoff=55))

## -----------------------------------------------------------------------------
library("ManyIV")
## Specification as in Table V, columns (1) and (2) in Angrist and Krueger
IVreg(lwage~education+as.factor(yob)|as.factor(qob)*as.factor(yob),
      data=ak80, inference=c("standard", "re", "il", "lil"))

## -----------------------------------------------------------------------------
r1 <- IVreg(lwage~education+as.factor(yob)|as.factor(qob)*as.factor(yob),
            data=ak80, inference="md", approx=TRUE)
print(r1, digits=4)

## -----------------------------------------------------------------------------
IVoverid(r1)

