

#' Printing an ICP output
#'
#' A simple print function for the \code{\link{ICP}} function.
#'
#' @param x an object of class "ICP", usually, a result of a call to
#'   \code{\link{ICP}}
#'
#'
#' @seealso The Invariant Causal Prediction function \code{\link{ICP}}.
#'
#' @examples
#' n <- 100
#' E <- rbinom(n, 3, 0.3)
#' X <- rnorm(n, E, 1)
#' Y <- rpois(n, X)
#' method <- structure(list(model = "glm", family = "poisson"),
#'                     class = "EnvirIrrel")
#' ICP(Y, X, E, method, level = 0.05)
#'
#' @export

# summary.ICP ==================================================================
print.ICP <- function(x) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), sep = "")

  cat("\n")
  print(x$method)


  if (!is.null(x$level)) {
    cat("\n\nAccepted model at level ", x$level, ":\n", sep = "")
    print(x$accepted.model)
    if (x$empty.accepted) {
      cat("Note: The empty set was accepted")
    }
  }

  cat("\n\nModels:\n", sep = "")
  pval_models <- sapply(x$model.analysis$pval, function(m) {
    if (m >= 0.5) {
      c("invariant at 0.5 ")
    } else if (m >= 0.1) {
      c("invariant at 0.1 ")
    } else if (m >= 0.05) {
      c("invariant at 0.05")
    } else if (m < 0.05) {
      c("                 ")
    }
  })
  pval_models <- data.frame("p value" = format.pval(x$model.analysis$pval,
                                                    digits = 3),
                            pval_models,
                            row.names = x$model.analysis$model,
                            fix.empty.names = F)
  print(pval_models)

  cat("\n")

  if (!is.null(x$variable.analysis)) {
    cat("\nVariables:\n", sep = "")
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
                                                   digits = 3),
                           pval_var,
                           row.names = x$variable.analysis$variables,
                           fix.empty.names = F)
    print(pval_var)
    cat("---\nSignif. codes:  0 ","***"," 0.001 ","**"," 0.01 ","*"," 0.05 ",
        "."," 0.1 ", " ", " 1\n", sep = "'")

    cat("\n")
  }
}

