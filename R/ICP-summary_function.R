
#' @export

# summary.ICP ==================================================================
print.ICP <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")

  cat("Model and Method Information:")
  print(x$method)
  cat("\n")

  if (!is.null(x$accepted.model)) {
    print.model_analysis(x, call = FALSE)
  }

  cat("\n\nModel Analysis:\n")
  if (is.logical(x$model.analysis$pval)) {
    pval_models <- data.frame(x$model.analysis$pval,
                              row.names = x$model.analysis$model)
    names(pval_models) <- paste0("Invariant at ", x$level)
  } else {
    pval_models <- sapply(x$model.analysis$pval, function(m) {
      if (m >= 0.5) {
        c("III")
      } else if (m >= 0.1) {
        c("II ")
      } else if (m >= 0.05) {
        c("I  ")
      } else if (m >= 0.01) {
        c("i  ")
      } else if (m < 0.01) {
        c("   ")
      }
    })
    nice_ps <- format.pval(x$model.analysis$pval,
                           digits = 3,
                           eps = x$smallest_possible_pval)
    pval_models <- data.frame("p value" = nice_ps,
                              pval_models,
                              row.names = x$model.analysis$model,
                              check.names = FALSE,
                              fix.empty.names = FALSE)
  }

  print(pval_models)
  if (is.logical(x$model.analysis$pval)) {
    cat("\n")
  } else {
    cat("---\nInvariance codes:  1 ","III"," 0.5 ","II"," 0.1 ","I"," 0.05 ",
        "i"," 0.01 ", " ", " 0\n", sep = "'")
  }

  if (!is.null(x$variable.analysis)) {
    cat("\n")
    print.variable_analysis(x, call = FALSE)
  }
  invisible(x)
}

