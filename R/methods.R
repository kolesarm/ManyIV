#' @method print IVResults
#' @export
print.IVResults <- function(x, digits = getOption("digits"), ...) {
    if (!is.null(x$call))
        cat("Call:\n", deparse(x$call), sep = "", fill=TRUE)

    r <- x$estimate[!is.na(x$estimate$beta), ]
    colnames(r) <- c("Estimate", colnames(r)[-1])

    if ("se" %in% colnames(r))
        colnames(r) <- c(c("Estimate", "Conventional", "Conv. (robust)"),
                         colnames(r[, -(1:3)]))

    cat("\nFirst-stage F: ", x$IVData$F, "\n\n")
    cat("Estimates and standard errors:\n")

    print.data.frame(r, digits = digits, ...)


    invisible(x)
}
