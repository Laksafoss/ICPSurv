#' Invariant Model Ranking
#'
#' A method for ranking models based on invariance.
#'
#' LONG DESCRIPTION
#'
#' @param M a list of models.
#' @param method the method used for calculating invariance scores. Note that
#'   the different methods of calculations are not necessarily comparable.
#'
#'   The default method is \code{BS}, which fits a b-spline and reports the
#'   curvature as an invariance score. For more detailes see
#'   \code{\link{construct_analysis.BS}}.
#' @param ... additional arguments to be passed to the
#'   \code{\link{construct_analysis}} function
#'
#' @return \code{IMR} returns a \code{\link[base]{data.frame}} of
#'   \code{\link[base]{class}}"\code{IMR}". The \code{\link[base]{data.frame}}
#'   has two columns:
#'   \describe{
#'     \item{model}{The model.}
#'     \item{roughness}{invariance score}
#'   }
#'
#' @export
IMR <- function(M, method = "BS", ...) {
  if (!is.list(M)) {
    stop("'M' must be a list of fitted models.")
  }
  if (length(M) < 2) {
    stop("'M' must contain at least two fitted models.")
  }
  if (!is.character(method)) {
    stop("'method' must be a character recognised as a method for the ",
         "'construct_analysis' function.")
  }
  M_names <- names(M)
  if (is.null(M_names)) {
    names(M) <- paste0("fit", seq_along(M))
  } else if (any(M_names == "")) {
    missing <- which(M_names == "")
    names(M)[missing] <- paste0("fit", missing)
  }
  mc <- sapply(M, function(m) {class(m)})
  class(mc) <- method

  # Find First
  data <- test_find_data(find_data(M[[1]]))
  constructed <- test_constructor(construct_analysis(mc, data = data, ...))
  analysis <- constructed$analysis
  res_first <- c(list("Model" = data$name), analysis(data))

  # Find the rest
  res_rest <- lapply(M[-1], function(m) {
    tmp <- test_find_data(find_data(m))
    if (! identical(tmp$time, data$time)) {
      err <- character()
      if (mc[[1]] == "aareg") {
        err <- paste0("Possible troubleshooting:\n Look up the nmin option in ",
                      "the survival::aareg function documentation.")
      }
      err <- paste(paste0("The models are not comparable - ",
                          "the respons variables are not identical."),
                   err,
                   sep = "\n")
      stop(err)
    }
    res <- c(list("Model" = tmp$name), analysis(tmp))
    return(res)
  })

  res_full <- c(list(res_first), res_rest)
  names(res_full) <- names(M)
  model_sum <- data.frame(t(sapply(res_full, function(r) {
    return(c(Model = r[["Model"]], Score = r[["Score"]]))
  })), stringsAsFactors = FALSE)
  model_sum[,2] <- as.numeric(model_sum[,2])

  tmp <- structure(list(ranking = model_sum,
                        results = res_full,
                        method = constructed$method,
                        call = match.call()),
                   class = c("IMR"))

  return(tmp)
}

test_constructor <- function(constructed) {
  err <- paste0("A 'construct_analysis' method must output a list with three ",
                "enties: 'analysis', 'method' and 'critical'.")
  if (!is.list(constructed)) {
    stop(err)
  }
  if (is.null(constructed$analysis)) {
    stop(paste0(err, "\n An 'analysis' function is not returned."))
  }
  if (!is.function(constructed$analysis)) {
    stop(paste0(err, "\n The constructed 'analysis' is not a function."))
  }
  return(constructed)
}
test_find_data <- function(data) {
  err <- "'find_data' must output a list."
  if (!is.list(data)) { stop(err) }
  if (is.null(data$name)) { stop(paste0(err, "\n'name' is missing.")) }
  if (is.null(data$time)) { stop(paste0(err, "\n'time' is missing.")) }
  return(data)
}

#' @export
print.IMR <- function(x, ...) {
  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")
  cat("\nSummary of Ranking Method:")
  print(x$method$method_expanded, right = FALSE)

  cat("\nInvariance Ranking:\n")
  x$ranking[,2] <- signif(x$ranking[,2], digits = 4)
  print(x$ranking)
  cat("---\n", x$method$critical)
  invisible(x)
}

#' Internal Method for Extracting Data for IMR
#'
#' An internal generic method for extracting fitted response and standerdize
#' residuals.
#'
#' LONG DESCRIPTION
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
#'
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



#' Plot IMR Analysis
#'
#' THIS FUNCTION IS DEFINITELY NOT DONE !!
#'
#' LONG DESCRIPTION OF PLOTTER FUNCTION
#'
#' @param x TODO
#' @param y TODO
#' @param which TODO
#' @param ... TODO
#'
#' @return DESCRIPTION OF RETURNED
#'
#' @examples
#' # Make super good example
#'
#' @export
plot.IMR <- function(x, y = NULL, which = seq_along(x$results), ...) {
  tt <- seq(min(x$method$time), max(x$method$time), length.out = 100)
  X_base <- fda::eval.basis(tt, x$method$basis)
  if (!is.null(x$method$knots)) {
    knot_base <- fda::eval.basis(x$method$knots, x$method$basis)
  }
  graphics::par(mfrow=c(2,2))
  lapply(which, function(i) {
    fit <- x$results[[i]]
    curve <- X_base %*% fit$coef
    graphics::plot(0, type = "n",
         xlim = c(0, max(x$method$time)),
         ylim = c(min(0, min(curve), min(fit$normhaz)),
                  max(1, max(curve), max(fit$normhaz))),
         xlab = "time", ylab = "fitted vals. norm.",
         main = names(x$results)[i])
    graphics::lines(tt, curve, col = 2)
    if (x$method$method == "BS") {
        y_knots <- knot_base %*% fit$coef
        graphics::points(x$method$time, fit$normhaz)
        graphics::points(x$method$knots, y_knots, pch = ".", cex = 6, col = 2)
    } else if (x$method$method == "rss") {
      graphics::segments(x0 = x$method$time,
                 y0 = fit$normhaz,
                 y1 = x$method$time / max(x$method$time), col = 2)
      graphics::points(x$method$time, fit$normhaz)
    } else {
      n <- length(x$method$time)
      graphics::polygon(c(0,x$method$time[rep(seq_len(n), each = 2)]),
              c(0,0,fit$normhaz[rep(seq_len(n-1), each = 2)],1), col = 2)
      graphics::points(x$method$time, fit$normhaz, type = "s")
    }
  })
  graphics::par(mfrow=c(1,1))
}
