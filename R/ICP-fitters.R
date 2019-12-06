#' Internal regression functions
#'
#' The \code{fit_model} function is a generic function meant for internal use in
#' the \code{ICPSurv} package, and as such it is not exported to the namespace.
#'
#' MORE DESCRIPTION
#'
#' MENTION SPECIAL CASE COX AND AALEN NOT CONST.
#'
#' @param method a \strong{method object} created by the
#'   \code{\link{method_obj}} function.
#' @param Y a vector or \code{\link[survival]{Surv}} object describing the
#'   target variable.
#' @param X a matrix, vector or data frame describing the covariates. WRITE
#'   ABOUT \code{X} WHEN TESTING THE EMPTY SET !!!!
#' @param id TODO
#' @param ... additional arguments. At the moment only the additional argument
#'   \code{const} is allowed.
#'
#' @return Depending on \code{method} the \code{fit_model} function returns
#'   different lists. If \code{method$method} \eqn{!=} "\code{ConstTime}" the output list must
#'   contain the following
#'   \describe{
#'     \item{coefficients}{the estimated regression coefficients.}
#'     \item{covariance}{the estimated covariance matrix.}
#'     \item{deviance}{the deviance, i.e. \eqn{- 2 * LogLik}, where \eqn{LogLik}
#'       is the log likelihood evaluated in the reported estimators.}
#'     \item{df}{degrees of freedom}
#'   }
#'
#'   If \code{method} has class "\code{ConstTime}" then the output list contains
#'   ANNA: DO WE ALWAYS RETURN ALL 4 ??? If we don't need pvals then we may set
#'   n.sim = 0
#'
#'   \describe{
#'     \item{cum}{cumulative timevarying regression coefficients estimated within
#'       the timeinterval with observations.}
#'     \item{var.cum}{rhe martingale based pointwise variance estimated for
#'       cumulatives.}
#'     \item{pval.testBeqC}{Kolmogorov-Smirnov test p-values based on resampling.}
#'     \item{pval.testBeqC.is}{Carmer von Mises test p-values based on resampling.}
#'    }
#'
#' @seealso \code{\link{ICP}}, \code{\link{plausible_predictor_test}} for the
#'   full wrapper functions.

fit_model <- function(method, Y, X, ...) {
  UseMethod("fit_model", method)
}

#' @rdname fit_model
fit_model.default <- function(method, Y, X, ...) {
  warning(method$model,
          "is not a recocnised statistical model in the ICPSurv framework",
          " - the recognised models are linear models 'lm', generalized",
          " linear models 'glm' (distribution family must be specified),",
          " proportional hazard models 'ph' and",
          " additive hazard models 'ah'.")
  return(1)
}

#' @rdname fit_model
fit_model.lm <- function(method, Y, X, ...) {
  fit <- stats::lm(Y ~ ., data = data.frame(Y,X))
  return(list("coefficients" = fit$coefficients,
              "covariance" = stats::vcov(fit),
              "deviance" = stats::deviance(fit),
              "df" = fit$df.residual,
              "scale" = stats::deviance(fit) / fit$df.residual))
}

#' @rdname fit_model
fit_model.glm <- function(method, Y, X, ...) {
  if (is.null(method$family)) {
    stop("For model = 'glm' there must be specified a family in the method object")
  }
  fit <- suppressWarnings(stats::glm(Y ~ .,
                                     data = data.frame(Y,X),
                                     family = method$family))
  #coef <- fit$coefficients
  #cov <- stats::vcov(fit)
  #if (any(is.na(coef))) {
  #  index <- which(in.na(coef))
  #  coef[index] <- rep(0, length(index))
  #  cov[,index] <- cov[index,] <- rep(0, length(index))
  #}
  return(list("coefficients" = fit$coefficients,
              "covariance" = stats::vcov(fit),
              "deviance" = fit$deviance,
              "df" = fit$df.residual,
              "scale" = summary(fit, dispersion = method$dispersion)$dispersion))
}

#' @rdname  fit_model
fit_model.ph <- function(method, Y, X, ...) {
  if (! survival::is.Surv(Y)) {
    Y <- survival::Surv(Y)
  }
  if (length(X) == 0) {
    fit <- survival::survreg(Y ~ 1, dist = "exponential")
  } else {
    fit <- survival::survreg(Y ~ ., data = data.frame(X), dist = "exponential")
  }
  return(list(
    "coefficients" = fit$coefficients,
    "covariance" = fit$var,
    "deviance" = - 2 * fit$loglik[2],
    "df" = fit$df.residual,
    "scale" = 1
  ))
}



# library(cmprsk)
# library(survival)
# n <- 5000
# E <- rbinom(n, 1, 0.2)
# X <- rnorm(n, E, 1)
# Y <- rexp(n, exp(0.5 * X))
# R <- rexp(n, exp(0.5 * X))
# C <- rexp(n, exp(0.5 * X))
# E <- as.factor(E)
# time <- pmin(Y, R, C)
# summary(time)
# status <- ifelse(time == Y, 1, ifelse(time == R, 2, 0))
# table(status)
# str(E)
# mm <- model.matrix(~X+E-1)
# head(mm)
# str(mm)
# fit <- crr(time, status, X, failcode = 2, cencode = 0)
# list("coefficient" = fit$coef,
#      "covariance" = fit$var,
#      "deviance" = - 2 * fit$loglik,
#      "df" = length(fit$coef),
#      "scale" = 1)


#' @rdname fit_model
fit_model.ah <- function(method, Y, X, ...) {
  if (! survival::is.Surv(Y)) {
    Y <- survival::Surv(Y)
  }
  dist <- survival::survreg.distributions$exponential
  dist$trans <- dist$itrans <- function(y) y
  dist$dtrans <- function(y) rep(1,length(y))
  if (length(X) == 0) {
    fit <- survival::survreg(Y ~ 1, dist = dist)
  } else {
    fit <- survival::survreg(Y ~ ., data = data.frame(X), dist = dist)
  }
  return(list(
    "coefficients" = fit$coefficients,
    "covariance" = fit$var,
    "deviance" = - 2 * fit$loglik[2],
    "df" = fit$df.residual,
    "scale" = 1))
}

#' @rdname fit_model
fit_model.hazard <- function(method, Y, X, ...) {
  if (! survival::is.Surv(Y)) {
    Y <- survival::Surv(Y)
  }
  if (length(X) == 0) {
    fit <- survival::survreg(Y ~ 1, dist = method$dist)
  } else {
    fit <- survival::survreg(Y ~ ., data = data.frame(X), dist = method$dist)
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
fit_nonparam_model <- function(method, Y, X, ...) {
  UseMethod("fit_nonparam_model", method)
}

#' @rdname fit_model
fit_nonparam_model.default <- function(method, Y, X, ...) {
  stop("A time varying fitter function is not defined for this model class:",
       class(method))
}

#' @rdname fit_model
fit_nonparam_model.ph <- function(method, Y, X, id = NULL,  ...) {
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
fit_nonparam_model.ah <- function(method, Y, X, id = NULL, ...) {
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

#' @rdname fit_model
fit_nonparam_model.hazard <- function(method, Y, X, id = NULL, ...) {
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
    fit <- timereg::aalen(Y ~ 1, id = id, n.sim = n.sim, robust = robust)
  } else if (method$link %in% c("proportional", "log")) {
    fit <- quiet(timereg::timecox(Y ~ ., data = data.frame(X),
                                  id = id, n.sim = n.sim, robust = robust))
  } else if (method$link %in% c("additive", "identity")) {
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

  stop("The 'nonparam' method is only implementer for the cox and aalen hazard",
       "models so far.")
}


quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

