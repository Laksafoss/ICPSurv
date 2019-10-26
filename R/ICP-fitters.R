
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
#' @param const can be used in the \code{coxph} and
#'   \code{aalen} methods and must be logical. If \code{const = FALSE} then the
#'   survival model fitted will have time varying regression effects.
#'   Furthermore the return object is changed.
#' @param ... additional arguments. At the moment only the additional argument
#'   \code{const} is allowed.
#'
#' @return Depending on \code{method} the \code{fit_model} function returns
#'   different lists. If \code{method$method != 'ConstTime'} the output list must
#'   contain the following
#'
#'   \item{coefficients}{the estimated regression coefficients.}
#'   \item{covariance}{the estimated covariance matrix.}
#'   \item{deviance}{the deviance, i.e. \eqn{- 2 * LogLik}, where \eqn{LogLik}
#'     is the log likelihood evaluated in the reported estimators.}
#'   \item{null.deviance}{the deviance of the null hypothesis, i.e.
#'     \eqn{- 2 * LogLik}, where \eqn{LogLik} is the log likelihood of the null
#'     model  - the model with predictor effect zero.}
#'   \item{df}{degrees of freedom}
#'
#'   If \code{method$method == 'ConstTime'} then the output list contains
#'   ANNA: DO WE ALWAYS RETURN ALL 4 ??? If we don't need pvals then we may set n.sim = 0
#'
#'   \item{cum} cumulative timevarying regression coefficients estimated within
#'     the timeinterval with observations.
#'   \item{var.cum} the martingale based pointwise variance estimated for
#'     cumulatives.
#'   \item{pval.testBeqC} Kolmogorov-Smirnov test p-values based on resampling.
#'   \item{pval.testBeqC.is} Carmer von Mises test p-values based on resampling.
#'
#'
#' @seealso \code{\link{ICP}}, \code{\link{plausible_predictor_test}} for the
#'   full wrapper functions.
#'
#' @examples
#' # create some data, here we simulate from a proportional hazard model
#' X <- rnorm(100)
#' Y <- (-runif(100)/(0.0001 * exp(X)))^(1/2.7)
#' Y <- survival::Surv(Y)
#'
#' # specify method object, we can fit proportional hazard models with coxph
#' method <- method_obj(method = "EnvirIrrel", model = "coxph")
#' # one can note that the 'method' part of the method object
#' # is not acually used in the fitter function;
#' # it is however used by the wrapper function 'plausible_predictor_test'
#'
#' # fit
#' fit_model(method, Y, X)
#'
#'
#'

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
  fit <- stats::lm(Y ~ ., data = data.frame(Y,X), family = method$family)
  return(list("coefficients" = fit$coefficients,
              "covariance" = vcov(fit)))
}


#' @rdname fit_model

fit_model.glm <- function(method, Y, X, ...) {
  if (is.null(method$family)) {
    stop("For model = 'glm' there must be specified a family in the method object")
  }
  fit <- stats::glm(Y ~ ., data = data.frame(Y,X), family = method$family)
  return(list("coefficients" = fit$coefficients,
              "covariance" = vcov(fit)))
}


#' @rdname  fit_model

fit_model.ph <- function(method, Y, X, id,  ...) {
  if (! survival::is.Surv(Y)) {
    Y <- survival::Surv(Y)
  }
  if (method$method == "ConstTime") {
    # TODO: IS THIS CORRECT ?!?
    fit <- timereg::timecox(Y ~ ., data = data.frame(X), id = id, n.sim = 0)
    return(list(
      "cum" = fit$cum,
      "cum.var" = fit$var.cum
    ))
  } else {
    fit <- survival::survreg(Y ~ ., data = data.frame(X), dist = "exponential")
    return(list(
      "coefficients" = fit$coefficients,
      "covariance" = vcov(fit)
    ))
  }

}


#' @rdname fit_model
# @importFrom timereg const

fit_model.ah <- function(method, Y, X, id, ...) {
  if (! survival::is.Surv(Y)) {
    Y <- survival::Surv(Y)
  }
  if (method$method == "ConstTime") {
    # TODO: IS THIS CORRECT ?!?!
    fit <- timereg::aalen(Y ~ ., data = data.frame(X), id = id, n.sim = 0)
    return(list(
      "cum" = fit$cum,
      "cum.var" = fit$var.cum
    ))
  } else {
    # TODO: HOW TO MAKE THE CORRECT
    dist <- survival::survreg.distributions$exponential
    dist$trans <- function(y) y # this does not work in small tests.....
                                # it seems to work now -- but whats with the log(scale)
    fit <- survival::survreg(Y ~ ., data = data.frame(X))
    return(list(
      "coefficients" = fit$coefficients,
      "covariance" = vcov(fit)
    ))
  }
}






# OLD : : : : @rdname fit_model
#
#fit_model.glm <- function(method, Y, X, ...) {
#  if (is.null(method$family)) {
#    stop("For model = 'glm' there must be specified a family in the method object")
#  }
#  fit <- stats::glm(Y ~ ., data = data.frame(Y,X), family = method$family)
#  return(list("coefficients" = fit$coefficients,
#              "covariance" = vcov(fit),
#              "deviance" = fit$deviance,
#              "null.deviance" = fit$null.deviance,
#              "df" = sum(!is.na(fit$coefficients))))
#}


#fit_model.survreg <- function(method, Y, X, ...) {
#  if (survival::is.Surv(Y) == FALSE) {
#    Y <- survival::Surv(Y)
#  }
#  if (is.null(method$dist)) {
#    stop("When using the 'survreg' fitting method a distribution must be",
#         "specified in the method object with the name 'dist'.")
#  }
#  scale <- ifelse(is.null(method$scale), 0, method$scale)
#  fit <- survival::survreg(Y ~ ., data = data.frame(X),
#                           dist = method$dist, scale = scale)
#
#  if (length(fit$coefficients) == fit$df) {
#    coef <- fit$coefficients
#    var <- fit$var
#  } else {
#    coef <- c(fit$coefficients, fit$scale)
#    var <- fit$var
#  }
#
#  return(list("coefficients" = coef,
#              "covariance" = var,
#              "deviance" = (- 2 * fit$loglik[2]),
#              "null.deviance" = ( - 2 * fit$loglik[1]),
#              "df" = fit$df))
#}


# @rdname fit_model
#fit_model.coxph <- function(method, Y, X, const = TRUE, ...) {
#  if (survival::is.Surv(Y) == FALSE) {      # this should be moved out to ICP
#    Y <- survival::Surv(Y)                  # it is too expensive to have in here
#  }
#  if (const) {
#    if (length(X) == 0) {
#      fit <- survival::coxph(Y ~ 1)
#    } else {
#      fit <- survival::coxph(Y ~ ., data = data.frame(X))
#    }
#    res <- list("coefficients" = fit$coefficients,
#                "covariance" = vcov(fit),
#                "deviance" = (- 2 * fit$loglik[2]),
#                "null.deviance" = (- 2 * fit$loglik[1]),
#                "df" = sum(!is.na(fit$coefficients)))
#  } else {
#    if (length(X) == 0) {
#      fit <- timereg::timecox(Y ~ 1, n.sim = 100) # TODO can we set n-sim = 0?
#    } else {
#      fit <- timereg::timecox(Y ~ ., data = data.frame(X), n.sim = 100) # TODO can we set n-sim = 0?
#    }
#    res <- list(pvals = c(fit$pval.testBeqC, fit$pval.testBeqC.is))
#  }
#  return(res)
#}







# fit_model.aalen <- function(method, Y, X, const = TRUE, ...) {
#   if (survival::is.Surv(Y) == FALSE) {  # This should be moved out to ICP
#     Y <- survival::Surv(Y)              # it is too expencive to have in here
#   }
#   if (!is.null(method$n.sim)) {
#     n.sim <- method$n.sim
#   } else {
#     n.sim <- 100
#   }
#   if (const) {
#     if (length(X) == 0) {
#       fit <- survival::aareg(Y ~ 1)
#     } else {
#       fit <- survival::aareg(Y ~ ., data = data.frame(X))
#     }
#     res <- list("coefficients" = fit$coefficients,
#                 "covariance" = vcov(fit),
#                 "deviance" = (- 2 * fit$loglik[2]),
#                 "null.deviance" = (- 2 * fit$loglik[1]),
#                 "df" = sum(!is.na(fit$coefficients)))
#   }
#   if (const) {
#     if (length(X) == 0) {
#       fit <- timereg::aalen(Y ~ 1, n.sim = 0) # is this right ??
#     } else {
#       X <- data.frame(X)
#       fit <- eval(parse(text = paste0(
#         "timereg::aalen(Y ~ ",
#         paste0("const(", colnames(X), ")", collapse = " + "),
#         ", data = X, n.sim = 0)"
#       )))
#     }
#     coef <- fit$gamma
#     var <- fit$robvar.gamma
#     res <- list("coefficients" = coef,
#                 "covariance" = var,
#                 "deviance" = NULL,
#                 "null.deviance" = NULL,
#                 "df" = NULL)
#   } else {
#     if (length(X) == 0) {
#       fit <- timereg::aalen(Y ~ 1, n.sim = n.sim)
#       res <- list(pvals = c(fit$pval.testBeqC, fit$pval.testBeqC.is))
#     } else {
#       fit <- timereg::aalen(Y ~ ., data = data.frame(X), n.sim = n.sim)
#       res <- list(pvals = c(fit$pval.testBeqC[-1], fit$pval.testBeqC.is[-1]))
#     }
#   }
#   return(res)
# }
