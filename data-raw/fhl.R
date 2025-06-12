## Download Farre-Mensa_Hegde_Ljungqvist_JF2019_MainData_NETSscrambled.dta from
## https://onlinelibrary.wiley.com/doi/full/10.1111/jofi.12867 -> supporting
## information, then unzipjofi12867-sup-0002-Replication-code.zip

dir <- paste0("~/research/jep_invitational/code/Farre-Mensa/",
              "Farre-Mensa_Hegde_Ljungqvist_JF2019_MainData_NETSscrambled.dta")
fm <- readstata13::read.dta13(dir)
fm <- fm[, c("gau", "examiner", "lenience", "pat_appl_year", "dallowed",
             "ln_total_patents_appl", "state_num", "ln_total_patents_approved",
             "ln_all_cites5", "ln_rounds_before_patent_APP")]
attr(fm, "expansion.fields") <- NULL
fm$gau <- as.factor(fm$gau)
fm$examiner <- as.factor(fm$examiner)
fm$ind_year <- factor(dplyr::group_indices(dplyr::group_by(fm, pat_appl_year,
                                                           gau)))
fhl <- fm[with(fm, order(ind_year, examiner)), ]
usethis::use_data(fhl, overwrite=TRUE, internal=FALSE)
