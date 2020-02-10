#' Find Accepted Model From ICP Object
#'
#' The \code{model_analysis} function takes an \code{\link{ICP}} object and a
#' significance level and outputs the accepted model at this level.
#'
#' The function \code{model_analysis} takes an \code{\link{ICP}} object and a
#' significance level and outputs the accepted model at this level:
#' \ifelse{html}{\out{<center>Accepted X columns = &#8745<sub>{S: H<sub>0,S</sub> accepted}</sub> S.</center>}}{\deqn{\textrm{Accepted } X \textrm{columns} = \cap_{\{S:\ H_{0,S} \textrm{ accepted}\}}\ S}}
#' See \code{\link{ICP}} for more details on the hypotheses
#' \ifelse{html}{\out{(H<sub>0,S</sub>)}}{\eqn{(H_{0,S})}}.
#'
#' The \code{gof} parameter protects against making statements when the model is
#' obviously not suitable for the data. If no model reaches the threshold
#' \code{gof} significance level, i.e. the p-values for
#' \ifelse{html}{\out{(H<sub>0,S</sub>)}}{\eqn{(H_{0,S})}} are all smaller then
#' \code{gof}, we report that there is no evidence for an invariant set.
#'
#' This function is also used internally in the \code{\link{ICP}} function
#' itself if \code{\link{ICP}} is called with \code{level} specified.
#'
#' @param x An ICP object
#' @param level The significance level required for a model to be considered
#'   invariant. The higher the significance \code{level} the stronger the
#'   requirements for godness of fit to attain 'invariance' status.
#' @param gof If no set of variables (including the empty set) leads to a
#'   p-value larger than the goodness-of-fit cutoff \code{gof}, the whole model
#'   will be rejected. If the model is correct, this will happen with a
#'   probability of gof and this option protects again making statements when
#'   the model is obviously not suitable for the data.
#'
#' @return \code{model_analysis} returns an object of \code{\link[base]{class}}
#' "\code{model_analysis}" containing the following components:
#' \item{gof }{ goodness-of-fit cufoff.}
#' \item{level }{ significance level of hypothesis tests.}
#' \item{accepted.model }{ the accepted model.}
#' \item{empty.message }{ if the empty set is returned as \code{accepted.model}
#'   then \code{empty.message} will give detailes.}
#' \item{call }{ the matched call.}
#'
#' @seealso
#'   The \code{model_analysis} function is also used internally in the
#'   \code{\link{ICP}} function itself if \code{\link{ICP}} is called with
#'   \code{level} specified.
#'
#'   \code{variable_analysis} is another function for summarizing
#'   \code{\link{ICP}} objects.
#'
#' @examples
#' n <- 100
#' E <- sample(5L, n, replace = TRUE)
#' X <- data.frame(X1 = rnorm(n, E, 1), X2 = rnorm(n, 3, 1))
#' Y <- rnorm(n, X$X1, 1)
#'
#' obj <- ICP(Y, X, E, fullAnalysis = TRUE)
#' model_analysis(obj, level = 0.05, gof = 0.15)
#'
#' @export
model_analysis <- function(x, level = 0.05, gof = 0.01) {
  if (class(x) != "ICP") {
    stop("'x' must be an object of class 'ICP'")
  }
  if (!is.numeric(level)) {
    stop("'level' must be a number between 0 and 1.")
  }
  if (level > 1 | level < 0) {
    stop("'level' must be a number between 0 and 1.")
  }
  if (!is.numeric(gof)) {
    stop("'gof' must be a number between 0 and 1.")
  }
  if (gof > 1 | gof < 0) {
    stop("'gof' must be a number between 0 and 1.")
  }
  if (!is.null(x$accepted.model) & !is.null(x$empty.message)) {
    if (x$level == level & x$gof >= gof) {
      return(structure(list(
        "gof" = gof,
        "level" = level,
        "accepted.model" = x$accepted.model,
        "empty.message" = x$empty.message,
        "call" = match.call()),
        class = "model_analysis"))
    }
  }
  if (!is.null(x$not.tested.list)) {
    if (x$level < level) {
      if (any(x$model.analysis$pval < level & x$model.analysis$pval >= x$level)) {
        stop("The inputted ICP object 'x' might be missing the results of some",
             " hypothesis tests: Rerun the ICP command that produced 'x' with",
             " the options 'fullAnalysis = TRUE' or 'level = ", level, "'.")
      }
    }
  }
  empty_message <- " "
  if (all(x$model.analysis$pval < gof)) {
    accepted_model <- "Empty"
    empty_message <- "Necessary goodness-of-fit not reached."
  } else {
    plausible <- which(x$model.analysis$pval >= level)
    if (length(plausible) == 0) {
      accepted_model <- "Empty"
      empty_message <- "No subset of predictors was accepted !"
    } else if (1 %in% plausible) {
      accepted_model <- paste0(deparse(x$call$Y), " ~ 1")
    } else {
      acc <- Reduce(intersect, x$tested.list[plausible])
      if (length(acc) == 0) {
        accepted_model <- "Empty"
        empty_message <- paste0(
          "Empty set model:",
          "\n    Intersection of plausible predictors is empty,",
          "\n    but the empty set does not result in invariance"
        )
      } else {
        accepted_model <- paste0(deparse(x$call$Y), " ~ ",
                                 paste0(x$covariate.list.code[-1,1][acc], collapse = " + "))
      }
    }
  }
  return(structure(list(
    "gof" = gof,
    "level" = level,
    "accepted.model" = accepted_model,
    "empty.message" = empty_message,
    "call" = match.call()),
    class = "model_analysis"))
}

#' @export
print.model_analysis <- function(x, call = TRUE, ...) {
  if (call) {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n", sep = "")
  }

  cat("\nAccepted Model At Significance Level ", x$level, ":\n", sep = "")
  if (x$accepted.model == "Empty") {
    cat("  ", x$empty.message, "\n")
  } else {
    cat("  ", x$accepted.model, "\n")
  }
}




#' Determine Significance Level of Individual Variables
#'
#' The \code{variable_analysis} function takes an \code{\link{ICP}} object and
#' outputs evidence for each individual variable being a direct cause.
#'
#' The function \code{variable_analysis} takes an \code{\link{ICP}} object and a
#' minimum goodness-of-fit cutoff \code{gof} and outputs the evidence for each
#' individual variable being a direct cause. To conduct a variable analysis the
#' inputted \code{\link{ICP}} object should be the result of a call to the
#' \code{\link{ICP}} function with option \code{fullAnalysis = TRUE}. This
#' ensures that \code{variable_analysis} has access to the p-values for all the
#' hypothses \ifelse{html}{\out{(H<sub>0,S</sub>)}}{\eqn{(H_{0,S})}} (see
#' \code{\link{ICP}} for an explanation of the hypotheses).
#'
#' If the needed \code{gof} level is reached the significance is calculated for
#' each variable \ifelse{html}{\out{X<sup>i</sup>}}{\eqn{X^i}} by
#' \ifelse{html}{\out{<center>p<sub>i</sub> = max{ p-value for H<sub>0,S</sub> &#8739 X<sup>i</sup> is not part of X<sup>S</sup> }.</center>}}{\deqn{p_i = \max\{\textrm{p-values for }H_{0,S}\ |\ X^i \textrm{ is not part of } X^S\}.}}
#'
#' The \code{gof} parameter protects against making statements when the model is
#' obviously not suitable for the data. If no model reaches the threshold
#' \code{gof} significance level, i.e. the p-values for
#' \ifelse{html}{\out{(H<sub>0,S</sub>)}}{\eqn{(H_{0,S})}} are all smaller then
#' \code{gof}, we report that there is no evidence for individual variables, as
#' there is no evidence for an invariant set.
#'
#' The \code{variable_analysis} function is also used internally in the
#' \code{\link{ICP}} function itself if \code{\link{ICP}} is called with
#' \code{level} specified and \code{fullAnalysis = TRUE}.
#'
#' @param x An ICP object
#' @param gof If no set of variables (including the empty set) leads to a
#'   p-value larger than the goodness-of-fit cutoff \code{gof}, the whole model
#'   will be rejected. If the model is correct, this will happen with a
#'   probability of gof and this option protects again making statements when
#'   the model is obviously not suitable for the data.
#'
#' @return
#'   \code{variable_analysis} returns a table with \code{\link[base]{class}}
#'   "\code{variable_analysis}" where each row reports the significance of the
#'   corresponding variable.
#'
#' @seealso
#'   The \code{variable_analysis} function is also used internally in the
#'   \code{\link{ICP}} function itself if \code{\link{ICP}} is called with
#'   \code{fullAnalysis = TRUE}.
#'
#'   \code{model_analysis} is another function for summarizing
#'   \code{\link{ICP}} objects.
#'
#' @examples
#' n <- 100
#' E <- sample(5L, n, replace = TRUE)
#' X <- data.frame(X1 = rnorm(n, E, 1), X2 = rnorm(n, 3, 1))
#' Y <- rnorm(n, X$X1, 1)
#'
#' obj <- ICP(Y, X, E, level = 0.05, fullAnalysis = TRUE)
#' variable_analysis(obj, gof = 0.1)
#'
#' # If we make extreem requirements for gof
#' # all "p-values" will be taken to be 1,
#' # i.e. no variable is a significant causal predictor
#' variable_analysis(obj, gof = 0.9)
#'
#' @export
variable_analysis <- function(x, gof = 0.1) {
  if (class(x) != "ICP") {
    stop("'x' must be an object of class 'ICP'")
  }
  if (! is.null(x$not.tested.list)) {
    stop("'ICP' object 'x' must be the result of a full analysis.\n")
  }
  if (!is.numeric(gof)) {
    stop("'gof' must be a number between 0 and 1.")
  }
  if (gof > 1 | gof < 0) {
    stop("'gof' must be a number between 0 and 1.")
  }
  if (!is.null(x$variable.analysis)) {
    if (x$gof == gof) {
      return(structure(list(
        "gof" = gof,
        "variable.analysis" = x$variable.analysis,
        "smallest_possible_pval" = x$smallest_possible_pval,
        "call" = match.call()),
        class = "variable_analysis"))
    }
  }
  ps <- data.frame(x$covariate.list.code, "confidence" = 1)
  ps <- ps[-1, ]
  if (any(x$model.analysis$pval >= gof)) {
    ps$confidence <- sapply(ps$code, function(i) {
      contained <- sapply(x$tested.list, function(j) {i %in% j})
      return(max(x$model.analysis$pval[!contained]))
    })
  }
  ps <- ps[,-2]
  return(structure(list(
    "gof" = gof,
    "variable.analysis" = ps,
    "smallest_possible_pval" = x$smallest_possible_pval,
    call = match.call()),
    class = "variable_analysis"))
}

#' @export
print.variable_analysis <- function(x, call = TRUE, ...) {
  if (call) {
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n", sep = "")
  }
  cat("\nVariable Analysis:\n")
  pval_var <- sapply(x$variable.analysis$confidence, function(v) {
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
  pval_var <- data.frame("p value" = format.pval(x$variable.analysis$confidence,
                                                 digits = 3,
                                                 eps = x$smallest_possible_pval),
                         pval_var,
                         row.names = x$variable.analysis$variable,
                         fix.empty.names = F)
  print(pval_var)
  cat("---\nSignif. codes:  0 ","***"," 0.001 ","**"," 0.01 ","*"," 0.05 ",
        "."," 0.1 ", " ", " 1\n", sep = "'")
  cat("Required Goodness-Of-Fit is", x$gof, "\n\n")
  invisible(x)
}
