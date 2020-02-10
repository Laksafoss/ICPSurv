#' Internal Method for Extracting Data for IMR
#'
#' An internal generic method for extracting fitted response and standerdize
#' residuals.
#'
#' This function finds the relevant model data and returns it in a suitable
#' manner to the \code{\link{IMR}} function.
#'
#' @param m a model object
#'
#' @return a list with the following three entries
#'   \describe{
#'     \item{X}{the fitted values of the model. In hazard models this is the
#'       time points corresponding to observed events.}
#'     \item{Y}{a list with the raw and standerdized residuals. In hazard models
#'       this is the standerdized residuals found when fitting a line through the
#'       estimated cumulative hazard.}
#'     \item{name}{a ("short") string describing the model. This value will be
#'       used as the model name in the final output.}
#'   }
#' @keywords internal
#' @export
find_data <- function(m) {
  UseMethod("find_data", m)
}

#' @rdname find_data
find_data.default <- function(m) {
  stop("The model class '", class(m),
       "' is not recognised as a 'find_data' method.")
}

#' @rdname find_data
find_data.lm <- function(m) {
  name <- paste0("lm(", format(stats::terms(m)), ")")
  return(list(time = m$model[,1], # what to put here
              fitted = m$fitted.values,
              residuals = stats::rstandard(m),
              name = name))
}

#' @rdname find_data
find_data.glm <- function(m) {
  name <- paste0("glm(", format(stats::terms(m)),
                 ", family = ",format(m$call$family),")")
  return(list(time = m$model[,1], # what to put here??
              fitted = m$fitted.values,
              residuals = stats::rstandard(m),
              name = name))
}

#' @rdname find_data
find_data.coxph <- function(m) {
  data <- survival::basehaz(m, centered = FALSE)
  form <- gsub("survival::", "", format(stats::terms(m)))
  if (!is.null(m$call$subset)) {
    form <- paste0(form, ", subset = ", m$call$subset)
  }
  name <- paste0("coxph(", form, ")")
  return(list(time = data$time,
              fitted = data$hazard,
              name = name))
}

#' @rdname find_data
find_data.aalen <- function(m) {
  form <- gsub("survival::", "", format(m$call$formula))
  if (!is.null(m$call$id)) {
    form <- paste0(form, paste0(", id = ", m$call$id))
  }
  if (!is.null(m$call$clusters)) {
    form <- paste0(form, paste0(", clusters = ", m$call$clusters))
  }
  name <- paste0("aalen(", form, ")")
  return(list(time = m$cum[,1],
              fitted = m$cum[,-1],
              name = name))
}

#' @rdname find_data
find_data.timecox <- function(m) {
  form <- gsub("survival::", "", format(m$call$formula))
  if (!is.null(m$call$id)) {
    form <- paste0(form, paste0(", id = ", m$call$id))
  }
  if (!is.null(m$call$clusters)) {
    form <- paste0(form, paste0(", clusters = ", m$call$clusters))
  }
  name <- paste0("timecox(", form, ")")
  if (sum(is.na(m$cum[2,])) > 0) {
    stop("The model", name, "has missing values in esimated cumulative effects")
  }
  return(list(time = m$cum[,1],
              fitted = m$cum[,-1],
              name = name))
}

#' @rdname find_data
find_data.cox.aalen <- function(m) {
  form <- gsub("survival::", "", format(m$call$formula))
  if (!is.null(m$call$id)) {
    form <- paste0(form, paste0(", id = ", m$call$id))
  }
  if (!is.null(m$call$clusters)) {
    form <- paste0(form, paste0(", clusters = ", m$call$clusters))
  }
  name <- paste0("cox.aalen(", form, ")")
  return(list(time = m$cum[,1],
              fitted = m$cum[,-1],
              name = name))
}

#' @rdname find_data
find_data.aareg <- function(m) {
  Y <- apply(m$coefficient, 2, function(x) {cumsum(x)})
  form <- gsub("survival::", "", format(m$call$formula))
  if (!is.null(m$call$subset)) {
    form <- paste0(form, ", subset = ", m$call$subset)
  }
  name <- paste0("aareg(", form, ")")
  return(list(time = m$times,
              fitted = Y,
              name = name))
}
