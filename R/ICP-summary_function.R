

#' Printing an ICP output
#'
#' A simple print function for the \code{\link{ICP}} function.
#'
#' @param x an object of class "ICP", usually, a result of a call to
#'   \code{\link{ICP}}
#' @param ... TODO
#'
#'
#' @seealso The Invariant Causal Prediction function \code{\link{ICP}}.
#'
#' @export

# summary.ICP ==================================================================
print.ICP <- function(x, ...) {
  cat("\n======   Call :   =================================================\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), sep = "")

  cat("\n\n======   Model and Method Informaiton :   =========================")
  print(x$method)

  if (!is.null(x$level)) {
    cat("\n======   Accepted Model At Level ", x$level,
        " :   =========================\n", sep = "")
    if (x$accepted.model == "NO SUBSET OF PREDICTORS WAS ACCEPTED") {
      cat(x$accepted.model, "\n")
    } else if (x$accepted.model == "Empty") {
      cat(deparse(x$call$Y), "~ 1\n")
      if (x$empty.accepted) {
        cat("(Note that the empty set was accepted)\n")
      }
    } else {
      cat(deparse(x$call$Y), "~", x$accepted.model, "\n")
    }
    #print(x$accepted.model)
    #if (x$empty.accepted) {
    #  cat("Note: The empty set was accepted")
    #}
  }

  cat("\n======   Models Analysis :   ======================================\n")

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
  if (is.logical(x$model.analysis$pval)) {
    nice_ps <- x$model.analysis$pval
  } else {
    nice_ps <- format.pval(x$model.analysis$pval,
                           digits = 3,
                           eps = x$smallest_possible_pval)
  }
  pval_models <- data.frame("p value" = nice_ps,
                            pval_models,
                            row.names = x$model.analysis$model,
                            fix.empty.names = F)
  print(pval_models)
  cat("---\nInvariance codes:  1 ","III"," 0.5 ","II"," 0.1 ","I"," 0.05 ",
      "i"," 0.01 ", " ", " 0\n", sep = "'")


  if (!is.null(x$variable.analysis)) {
    cat("\n\n======   Variables Analysis :   ===================================\n")
    pval_var <- sapply(x$variable.analysis$pval, function(v) {
      if (v <= 0.001) {
        c("***")
      } else if (v <= 0.01) {
        c("** ")
      } else if (v <= 0.05) {
        c("*  ")
      } else if (v <= 0.1) {
        c(".  ")
      } else {
        c("   ")
      }
    })
    pval_var <- data.frame("p value" = format.pval(x$variable.analysis$pval,
                                                   digits = 3,
                                                   eps = x$smallest_possible_pval),
                           pval_var,
                           row.names = x$variable.analysis$variables,
                           fix.empty.names = F)
    print(pval_var)
    cat("---\nSignif. codes:  0 ","***"," 0.001 ","**"," 0.01 ","*"," 0.05 ",
        "."," 0.1 ", " ", " 1\n\n", sep = "'")
  }
}

