#' Internal regression functions
#'
#' The \code{fit_model} function is a generic function meant for internal use in
#' the \code{ICPSurv} package, and as such it is not exported to the namespace.
#'
#' The \code{fit_model.X} and \code{fit_nonparam_model.X} functions are internal
#' fitter functions. They are usually all quite simple as their main
#' functionality apart from fitting the model is ensureing that the output
#' is compatible with the \code{\link{plausible_predictor_test}}s.
#'
#' All \code{fit_model} methods must return the following:
#' \itemize{
#'   \item{\code{coefficients} }{ The estimated regression coefficients.}
#'   \item{\code{deviance} }{ The deviance of the fitted model, i.e. \eqn{-2 log(likelihood)}}
#'   \item{\code{df} }{ Degrees of freedom in fitted model.}
#'   \item{\code{scale} }{ Scale.}
#' }
#'
#' All \code{fit_nonparam_model} methods must return the following:
#' \itemize{
#'   \item{\code{cum} }{ Cumulative regression coefficient.}
#'   \item{\code{cum.var} }{ Cumulative variance of regression effects.}
#' }
#' Morover if \code{n.sim} in the \code{method} is non-zero then the
#' \code{fit_nonparam_model} also returnes the following
#' \itemize{
#'   \item{\code{sup} }{ Kolmogorov–Smirnov test.}
#'   \item{\code{int} }{ Cramér–von Mises test.}
#' }
#'
#'
#' @param method a \strong{method object} created by the
#'   \code{\link{method_obj}} function.
#' @param Y a vector or \code{\link[survival]{Surv}} object describing the
#'   target variable.
#' @param X a matrix, vector or data frame describing the covariates.
#' @param subset an optional vector specifying a subset of observations to be
#'   used in the fitting process.
#' @param id for timevarying covariates the variable must associate each record
#'   with the id of a subject.
#' @param ... additional arguments to be passed to lower level functions.
#'
#' @return Both the \code{fit_model} and \code{fit_nonparam_model} methods return a list.
#'
#' \code{fit_model} methods must return the following:
#'   \item{\code{coefficients} }{ The estimated regression coefficients.}
#'   \item{\code{deviance} }{ The deviance of the fitted model, i.e. \eqn{-2 log(likelihood)}}
#'   \item{\code{df} }{ Degrees of freedom in fitted model.}
#'   \item{\code{scale} }{ Scale.}
#'
#'
#' \code{fit_nonparam_model} methods return the following:
#'   \item{\code{cum} }{ Cumulative regression coefficient.}
#'   \item{\code{cum.var} }{ Cumulative variance of regression effects.}
#'   \item{\code{sup} }{ Kolmogorov–Smirnov test (only returned if \code{n.sim} not zero).}
#'   \item{\code{int} }{ Cramér–von Mises test (only returned if \code{n.sim} not zero).}
#'
#'
#'
#' @seealso \code{\link{ICP}}, \code{\link{plausible_predictor_test}} for the
#'   full wrapper functions.
#' @keywords internal

fit_model <- function(method, Y, X, subset, ...) {
  UseMethod("fit_model", method)
}

#' @rdname fit_model
fit_model.default <- function(method, Y, X, subset = NULL, ...) {
  warning(method$model,
          "is not a recocnised statistical model in the ICPSurv framework",
          " - the recognised models are linear models 'lm', generalized",
          " linear models 'glm' (distribution family must be specified),",
          " proportional hazard models 'ph' and",
          " additive hazard models 'ah'.")
  return(1)
}

#' @rdname fit_model
fit_model.lm <- function(method, Y, X, subset = NULL, ...) {
  mf <- stats::model.frame(Y ~ ., data = data.frame(Y, X), subset = subset)
  fit <- stats::lm(mf)
  coef <- fit$coefficients
  var <- stats::vcov(fit)
  if (any(is.na(coef))) {
    index <- which(is.na(coef))
    coef[index] <- 0
    var[,index] <- var[index,] <- 0
    diag(var)[index] <- Inf
  }
  return(list("coefficients" = coef,
              "covariance" = var,
              "deviance" = stats::deviance(fit),
              "df" = fit$df.residual,
              "scale" = stats::deviance(fit) / fit$df.residual))
}

#' @rdname fit_model
fit_model.glm <- function(method, Y, X, subset = NULL, ...) {
  if (is.null(method$family)) {
    stop("For model = 'glm' there must be specified a family in the method object")
  }
  mf <- stats::model.frame(Y ~ ., data = data.frame(Y, X), subset = subset)
  fit <- suppressWarnings(stats::glm(mf, family = method$family))
  coef <- fit$coefficients
  var <- stats::vcov(fit)
  if (any(is.na(coef))) {
    index <- which(is.na(coef))
    coef[index] <- 0
    var[,index] <- var[index,] <- 0
    diag(var)[index] <- Inf
  }
  return(list("coefficients" = coef,
              "covariance" = var,
              "deviance" = fit$deviance,
              "df" = fit$df.residual,
              "scale" = summary(fit, dispersion = method$dispersion)$dispersion))
}

#' @rdname  fit_model
fit_model.ph <- function(method, Y, X, subset = NULL, ...) {
  if (! survival::is.Surv(Y)) {
    Y <- survival::Surv(Y)
  }
  if (length(X) == 0) {
    fit <- survival::survreg(Y ~ 1, subset = subset, dist = "exponential")
  } else {
    fit <- survival::survreg(Y ~ ., data = data.frame(X), subset = subset,
                             dist = "exponential")
  }
  if (any(diag(fit$var) == 0)) {
    ind <- which(diag(fit$var) == 0)
    diag(fit$var)[ind] <- Inf
  }
  return(list(
    "coefficients" = fit$coefficients,
    "covariance" = fit$var,
    "deviance" = - 2 * fit$loglik[2],
    "df" = fit$df.residual,
    "scale" = 1
  ))
}

#' @rdname fit_model
fit_model.ah <- function(method, Y, X, subset = NULL, ...) {
  if (! survival::is.Surv(Y)) {
    Y <- survival::Surv(Y)
  }
  dist <- survival::survreg.distributions$exponential
  dist$trans <- dist$itrans <- function(y) y
  dist$dtrans <- function(y) rep(1,length(y))
  if (length(X) == 0) {
    fit <- survival::survreg(Y ~ 1, subset = subset, dist = dist)
  } else {
    fit <- survival::survreg(Y ~ ., data = data.frame(X), subset = subset,
                             dist = dist)
  }
  if (any(diag(fit$var) == 0)) {
    ind <- which(diag(fit$var) == 0)
    diag(fit$var)[ind] <- Inf
  }
  return(list(
    "coefficients" = fit$coefficients,
    "covariance" = fit$var,
    "deviance" = - 2 * fit$loglik[2],
    "df" = fit$df.residual,
    "scale" = 1))
}

#' @rdname fit_model
fit_model.hazard <- function(method, Y, X, subset = NULL, ...) {
  if (! survival::is.Surv(Y)) {
    Y <- survival::Surv(Y)
  }
  if (length(X) == 0) {
    fit <- survival::survreg(Y ~ 1, subset = subset, dist = method$dist)
  } else {
    fit <- survival::survreg(Y ~ ., data = data.frame(X), subset = subset,
                             dist = method$dist)
  }
  return(list(
    "coefficients" = fit$coefficients,
    "covariance" = fit$var,
    "deviance" = - 2 * fit$loglik[2],
    "df" = fit$df.residual,
    "scale" = 1
  ))
}

#' @rdname fit_model
fit_nonparam_model <- function(method, Y, X, subset, ...) {
  UseMethod("fit_nonparam_model", method)
}

#' @rdname fit_model
fit_nonparam_model.default <- function(method, Y, X, subset = NULL, ...) {
  stop("A time varying fitter function is not defined for this model class:",
       class(method))
}

#' @rdname fit_model
fit_nonparam_model.ph <- function(method, Y, X, subset = NULL, id = NULL,  ...) {
  if (! survival::is.Surv(Y)) {
    Y <- survival::Surv(Y)
  }
  if (is.null(method$n.sim)) {
    n.sim <- 50
  } else {
    n.sim <- method$n.sim
  }
  if (length(X) == 0) {
    fit <- timereg::aalen(Y ~ 1, data = data.frame(X), id = id, n.sim = n.sim)
  } else {
    fit <- quiet(timereg::timecox(Y ~ ., data = data.frame(X),
                                  id = id, n.sim = n.sim))
  }
  if (sum(is.na(fit$cum[2,])) > 0) {
    stop("The model can not be fitted due to bugs in the dependencies.\n  ",
         "The error has been reported to the maintainers of the 'timereg' package.")
  }
  if (n.sim != 0) {
    res <- list(
      "cum" = fit$cum,
      "cum.var" = fit$var.cum,
      "sup" = unname(fit$pval.testBeqC),
      "int" = unname(fit$pval.testBeqC.is))
  } else {
    res <- list(
      "cum" = fit$cum,
      "cum.var" = fit$var.cum)
  }
  return(res)
}

#' @rdname fit_model
fit_nonparam_model.ah <- function(method, Y, X, subset = NULL, id = NULL, ...) {
  if (! survival::is.Surv(Y)) {
    Y <- survival::Surv(Y)
  }
  if (is.null(method$n.sim)) {
    n.sim <- 50
  } else {
    n.sim <- method$n.sim
  }
  robust <- ifelse(n.sim == 0, 0, 1)
  if (length(X) == 0) {
    fit <- timereg::aalen(Y ~ 1, data = data.frame(X), id = id,
                          n.sim = n.sim, robust = robust)
  } else {
    fit <- timereg::aalen(Y ~ ., data = data.frame(X), id = id,
                          n.sim = n.sim, robust = robust)
  }
  if (sum(is.na(fit$cum[2,])) > 0) {
    stop("The model can not be fitted due to bugs in the dependencies.\n  ",
         "The error has been reported to the maintainers of the 'timereg' package.")
  }
  if (n.sim != 0) {
    res <- list(
      "cum" = fit$cum,
      "cum.var" = fit$var.cum,
      "sup" = unname(fit$pval.testBeqC),
      "int" = unname(fit$pval.testBeqC.is))
  } else {
    res <- list(
      "cum" = fit$cum,
      "cum.var" = fit$var.cum)
  }
  return(res)
}

# @rdname fit_model
#fit_nonparam_model.hazard <- function(method, Y, X, id = NULL, ...) {
#  if (! survival::is.Surv(Y)) {
#    Y <- survival::Surv(Y)
#  }
#  if (is.null(method$n.sim)) {
#    n.sim <- 50
#  } else {
#    n.sim <- method$n.sim
#  }
#  robust <- ifelse(n.sim == 0, 0, 1)
#  if (length(X) == 0) {
#    fit <- timereg::aalen(Y ~ 1, id = id, n.sim = n.sim, robust = robust)
#  } else if (method$link %in% c("proportional", "log")) {
#    fit <- quiet(timereg::timecox(Y ~ ., data = data.frame(X),
#                                  id = id, n.sim = n.sim, robust = robust))
#  } else if (method$link %in% c("additive", "identity")) {
#    fit <- timereg::aalen(Y ~ ., data = data.frame(X), id = id,
#                          n.sim = n.sim, robust = robust)
#  }
#  if (sum(is.na(fit$cum[2,])) > 0) {
#    stop("The model can not be fitted due to bugs in the dependencies.\n  ",
#         "The error has been reported to the maintainers of the 'timereg' package.")
#  }
#  if (n.sim != 0) {
#    res <- list(
#      "cum" = fit$cum,
#      "cum.var" = fit$var.cum,
#      "sup" = unname(fit$pval.testBeqC),
#      "int" = unname(fit$pval.testBeqC.is))
#  } else {
#    res <- list(
#      "cum" = fit$cum,
#      "cum.var" = fit$var.cum)
#  }
#  return(res)
#
#  stop("The 'nonparam' method is only implementer for the cox and aalen hazard",
#       "models so far.")
#}


quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

