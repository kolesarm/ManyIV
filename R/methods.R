#' @export
print.IVResults <- function(x, digits = getOption("digits"), ...) {
    if (!is.null(x$call))
        cat("Call:\n", deparse(x$call), sep = "", fill=TRUE)

    r <- x$estimate[!is.na(x$estimate$estimate), ]
    colnames(r) <- c("Estimate", colnames(r)[-1])

    if ("se" %in% colnames(r))
        colnames(r) <- c(c("Estimate", "Conventional", "Conv. (robust)",
                           "HTE robust"),
                         colnames(r)[-(1:4)])
    if ("se_text" %in% colnames(r))
        colnames(r) <- c(c("Estimate", "Conv. robust SE",
                           "HTE robust SE"))

    cat("\nFirst-stage F: ", x$IVData$F, "\n\n")
    cat("Estimates and standard errors:\n")

    print.data.frame(r, digits = digits, ...)

    invisible(x)
}
