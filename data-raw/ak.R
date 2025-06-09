## Download NEW7080.rar from
## https://economics.mit.edu/people/faculty/josh-angrist/angrist-data-archive
## Then unrar to get .dta file

## Note v17: sob (FIPS code)
dir <- "~/teaching/Datasets/AngristKrueger1991/NEW7080.dta"
ak <- readstata13::read.dta13(dir)

names(ak) <- c("v1", "age", "v3", "education", "enocent", "esocent", "v7", "v8",
               "lwage", "married", "midatl", "mt", "neweng", "v14", "pacific",
               "census", "sob", "qob", "black", "smsa", "soatl", "v22", "v23",
               "wnocent", "wsocent", "v26", "yob")
## Clean up divisions
ak$division <- factor(ak$neweng + 2*ak$midatl + 3*ak$enocent + 4*ak$wnocent +
                          5*ak$soatl+ 6*ak$esocent + 7*ak$wsocent +
                          8*ak$mt + 9*ak$pacific,
                      levels=as.character(1:9),
                      labels=c("New England", "Mid Atlantic", "E N Central",
                               "W N Central", "S Atlantic", "E S Central",
                               "W S Central", "Mountain", "Pacific"))

ak$age[ak$census==80] <- ak$age[ak$census==80]-1900
ak$yob[ak$yob<1900] <- ak$yob[ak$yob<1900]+1900
ak$v1 <- NULL # Age in years. We only use age in quarters, v2
ak$black <- ak$black==1
ak$married <- ak$married==1
ak$smsa <- ak$smsa==1
ak$yob <- as.integer(ak$yob)
ak$qob <- factor(ak$qob, levels=as.character(1:4),
                 labels=c("Q1", "Q2", "Q3", "Q4"))
ak80 <- data.frame(ak[ak$census==80 & ak$yob<=1939 & ak$yob>=1930, ])

## 03 "American Samoa"
## 07 "Canal Zone"
## 14 "Guam"
## 43 "PR"
## 52 "Virgin Islands"
## 57 "Pacific coast"
## 99
ak80$sob <- factor(ak80$sob, levels=as.character(sort(unique(ak80$sob))),
                   labels=c("AL", "AK", "AZ", "AR", "CA", "CO", "CT", "DE",
                            "DC", "FL", "GA", "HI", "ID", "IL", "IN", "IA",
                            "KS", "KY", "LA", "ME", "MD", "MA", "MI", "MN",
                            "MS", "MO", "MT", "NE", "NV", "NH", "NJ", "NM",
                            "NY", "NC", "ND", "OH", "OK", "OR", "PA", "RI",
                            "SC", "SD", "TN", "TX", "UT", "VT", "VA", "WA",
                            "WV", "WI", "WY"))
## Keep 1980 census, 1930-1939 cohort, 329,509 observations
ak80 <- ak80[, c("age", "education", "lwage", "married", "qob", "sob",
                 "black", "smsa", "yob", "division")]
usethis::use_data(ak80, overwrite=TRUE, internal=FALSE)
