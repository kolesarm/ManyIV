#' Angrist and Krueger (1991) Census data
#'
#' Sample of males born in the US in 1930-39 from 5 percent sample of the 1980
#' US Census
#' @format A data frame with 329,509 observations on 10 variables:
#'
#' \describe{
#'
#' \item{age}{Age, measured at quarterly precision}
#'
#' \item{education}{Years of education}
#'
#' \item{lwage}{Log of weekly wage}
#'
#' \item{married}{Indicator for being married}
#'
#' \item{qob}{Quarter of birth}
#'
#' \item{sob}{State of birth}
#'
#' \item{black}{Indicator for being black}
#'
#' \item{smsa}{SMSA indicator}
#'
#' \item{yob}{Year of birth}
#'
#' \item{division}{Census division}
#'
#' }
#' @source Josh Angrist's website,
#'     \url{http://economics.mit.edu/faculty/angrist/data1/data/angkru1991}
#' @references{
#'
#' \cite{Angrist, Joshua D., and Alan B. Krueger. 1991. “Does Compulsory
#' School Attendance Affect Schooling and Earnings?” The Quarterly Journal
#' of Economics 106 (4): 979–1014. \doi{10.2307/2937954}.}
#'
#' }
"ak80"

#' Farre-Mensa, Hegde, and Ljungqvist (2020) data
#'
#' TODO: describe
#' @format A data frame with 34,435 observations on 11 variables:
#'
#' \describe{
#'
#' \item{gau}{Group Art Unit ID of the worker assigned to application}
#'
#' \item{examiner}{ID of worker assigned to application}
#'
#' \item{lenience}{Leniency measure: number of patents allowed / number of
#'                 patents examined to date.}
#'
#' \item{pat_appl_year}{Patent application year}
#'
#' \item{dallowed}{First patent application approved}
#'
#' \item{ln_total_patents_appl}{log(1 + subsequent patent applications)}
#'
#' \item{state_num}{State code}
#'
#' \item{ln_total_patents_approved}{log(1 + subsequent patent approvals)}
#'
#' \item{ln_all_cites5}{log(1 + total citations to subsequent patent
#'                      applications)}
#'
#' \item{ln_rounds_before_patent_APP}{log(1 + number of prior VC rounds)}
#'
#' \item{ind_year}{Art-unit-by-application-year fixed effects.}
#'
#' }
#' @source Journal of Finance website, \doi{10.1111/jofi.12867}
#' @references{
#'
#' \cite{Farre-Mensa, Hegde, and Ljungqvist. 2020. What Is a Patent Worth?
#' Evidence from the U.S. Patent "Lottery". The Journal of Finance 72
#' (2):639-682. \doi{10.1111/jofi.12867}.}
#'
#' }
"fhl"
