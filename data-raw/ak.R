## Download NEW7080.rar from
## http://economics.mit.edu/faculty/angrist/data1/data/angkru1991
## Then unrar to get .dta file

## Note v17: sob (FIPS code)
ak <- readstata13::read.dta13("~/teaching/Datasets/AngristKrueger1991/NEW7080.dta")
ak2 <- ak


names(ak) <- c("v1", "age", "v3", "education", "enocent", "esocent", "v7", "v8",
               "lwage", "married", "midatl", "mt", "neweng", "v14", "v15",
               "census", "sob", "qob", "black", "smsa", "soatl", "v22", "v23",
               "wnocent", "wsocent", "v26", "yob")
## Clean up divisions
ak$division <- factor(ak$neweng + 2*ak$midatl + 3*ak$enocent + 4*ak$wnocent +
                      5*ak$soatl+ 6*ak$esocent + 7*ak$wsocent + 8*ak$mt,
                      levels=as.character(0:8),
                      labels=c("Pacific", "New England", "Mid Atlantic",
                          "E N Central", "W N Central", "S Atlantic",
                          "E S Central", "W S Central", "Mountain"))
ak$neweng <- NULL
ak$enocent <- NULL
ak$esocent <- NULL
ak$wnocent <- NULL
ak$wsocent <- NULL
ak$midatl <- NULL
ak$soatl <- NULL
ak$mt <- NULL

ak$age[ak$census==80] <- ak$age[ak$census==80]-1900
ak$yob[ak$yob<1900] <- ak$yob[ak$yob<1900]+1900
ak$v1 <- NULL                           # age. We only use age in quarters
ak$black <- ak$black==1
ak$married <- ak$married==1
ak$smsa <- ak$smsa==1
ak$yob <- as.integer(ak$yob)
## Age: age at quarterly precision

## Keep 1980 census, 1930-1939 cohort, 329,509 observations
ak80 <- data.frame(ak[ak$census==80 & ak$yob<=1939 & ak$yob>=1930, ])
ak80$census <- NULL
ak80$v26 <- NULL                        # uniformly zero
ak80 <- ak80[, c("age", "education", "lwage", "married", "qob",
                 "black", "smsa", "yob", "division")]
devtools::use_data(ak80, overwrite=TRUE, internal=FALSE)
