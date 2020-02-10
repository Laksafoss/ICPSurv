#' Invariant Model Ranking
#'
#' A method for ranking models based on an invariance score.
#'
#' The \code{IMR} function takes a list of fitted models as input and computes
#' an invariance ranking. The invariance scoring method is specified via
#' \code{method}. It is the users responsibility to input comparable models in
#' \code{M} and choose a suiteble ranking scheem via \code{method}; only
#' rudimentery sanity check of the inputted models in \code{M} is conducted.
#'
#' At the moment the \code{IMR} function only takes hazard models as input, and
#' only two ranking methods have been implementd, namely the \code{BS} and the
#' \code{rss} methods:
#'
#' \strong{The \code{BS} method:}
#' For each model a basis spline is fitted to the cumulative coefficients
#' and the integrated curvature is reported as a measure of invariance. Before
#' the spline is fitted the cumulative coefficients are normalized. A smaller
#' score is evidence of a more invariant model.
#'
#' \strong{The \code{rss} method:}
#' For each model the corresponding invariant model is found and the residual
#' square error is found and reported as a measure of invariance. Before the
#' invariant model and residual square error is found the cumulative
#' coefficients are normalized. A smaller score is evidence of a more invariant
#' model.
#'
#' @param M A list of models. The models may be the result of a call to
#'   \code{\link[survival]{coxph}}, \code{\link[survival]{aareg}},
#'   \code{\link[timereg]{timecox}}, \code{\link[timereg]{aalen}} and
#'   \code{\link[timereg]{cox.aalen}}.
#' @param method The method used for calculating invariance scores. At the
#'   moment only two ranking methods have been implemented, namely the \code{BS}
#'   and \code{rss} method. The default method is \code{BS}, which fits a
#'   b-spline and reports the curvature as an invariance score.
#' @param ... additional arguments to be passed to the internal
#'   \code{construct_analysis} functions.
#'
#' @return \code{IMR} returns an object of \code{\link[base]{class}} "\code{IMR}"
#' with the following components:
#' \item{ranking }{ A data.frame summarizing the findings. For detailes on how
#' to interpret this see the method specific details above.}
#' \item{results }{ A list with the detailed results of the analysis.}
#' \item{method }{ A list with the \code{method}, a data.frame
#' \code{method_expanded} with a more detailed summary of the called method,
#' the character \code{critical} giving a one line interpretation guide to the
#' \code{ranking} data.frame. This list may also contain other method specific
#' components.}
#' \item{call }{ The matched call.}
#'
#'
#' @examples
#' ## Generate Data
#' n <- 1000
#' E <- sample(4L, n, replace = TRUE)
#' X <- data.frame(X1 = rnorm(n, 2 * E, 1),
#'                 X2 = rnorm(n, 2 * (E == 1) + 5 * (E == 4)  , 1),
#'                 X3 = rbinom(n, 1, 0.7))
#' Y <- rexp(n, exp(4 - 0.6 * X$X1 + 0.4 * X$X2))
#' C <- rexp(n, exp(0.5))
#' time <- pmin(Y, C)
#' status <- time == Y
#' # Note that X1 and X2 are the true causal predictors
#'
#' ## Fit Models
#' m1 <- survival::coxph(survival::Surv(time, status) ~ X1, data = X)
#' m2 <- survival::coxph(survival::Surv(time, status) ~ X2, data = X)
#' m3 <- survival::coxph(survival::Surv(time, status) ~ X3, data = X)
#' m12 <- survival::coxph(survival::Surv(time, status) ~ X1 + X2, data = X)
#' m13 <- survival::coxph(survival::Surv(time, status) ~ X1 + X3, data = X)
#'
#'
#' # Ranking via Basis Splines
#' IMR(list(m1, m2, m3, m12, m13), method = "BS")
#'
#' # Ranking via Residual Square Error
#' IMR(list(m1, m2, m3, m12, m13), method = "rss")
#'
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
  cat("---\n", x$method$critical, "\n")
  invisible(x)
}




#' Plot IMR Analysis
#'
#' A plotting method for \code{\link{IMR}} objects.
#'
#' Multiple plots (selectable by \code{which}) are created. The type of plot
#' depend on method of the inputted \code{\link{IMR}} object. At the moment the
#' following plotting methods are avaible:
#'
#' \strong{The \code{BS} method plotted:} The normalized cumulative coefficients
#' are plotted along with the fitted basis spline.
#'
#' \strong{The \code{rss} method plotted:} The normalized cumulative
#' coefficients are plotted along with the invariant model and residuals.
#'
#' @param x An \code{\link{IMR}} object.
#' @param y \code{NULL}. Not used at the moment.
#' @param which Which of the analyzed models to display.
#' @param ncol Number of columns is plotting grid.
#' @param ... additional arguments to be passed to lower level functions.
#'
#' @seealso See \code{\link{IMR}} for the construction of an \code{IMR} object.
#'
#' @examples
#' ## Generate Data
#' n <- 1000
#' E <- sample(4L, n, replace = TRUE)
#' X <- data.frame(X1 = rnorm(n, 2 * E, 1),
#'                 X2 = rnorm(n, 2 * (E == 1) + 5 * (E == 4), 1),
#'                 X3 = rbinom(n, 1, 0.5))
#' Y <- rexp(n, exp(4 - 0.6 * X$X1 + 0.4 * X$X2))
#' C <- rexp(n, exp(0.5))
#' time <- pmin(Y, C)
#' status <- time == Y
#' # Note that X1 and X2 are true causal predictors
#'
#' ## Fit Models
#' m1 <- survival::coxph(survival::Surv(time, status) ~ X1, data = X)
#' m2 <- survival::coxph(survival::Surv(time, status) ~ X2, data = X)
#' m3 <- survival::coxph(survival::Surv(time, status) ~ X3, data = X)
#' m12 <- survival::coxph(survival::Surv(time, status) ~ X1 + X2, data = X)
#' m13 <- survival::coxph(survival::Surv(time, status) ~ X1 + X3, data = X)
#'
#' ## Ranking via Basis Splines
#' BS <- IMR(list(m1, m2, m3, m12, m13), method = "BS")
#' plot(BS, ncol = 3)
#' plot(BS, which = c(1,2,3), ncol = 3) # the 1 variable models
#'
#' ## Ranking via Residual Square Error
#' RSS <- IMR(list(m1, m2, m3, m12, m13), method = "rss")
#' plot(RSS, ncol = 3)
#' @export
plot.IMR <- function(x, y = NULL, which = seq_along(x$results), ncol = 2, ...) {
  tt <- seq(min(x$method$time), max(x$method$time), length.out = 100)
  X_base <- fda::eval.basis(tt, x$method$basis)
  if (!is.null(x$method$knots)) {
    knot_base <- fda::eval.basis(x$method$knots, x$method$basis)
  }
  m <- length(which)
  nr <- ceiling(m / ncol)
  tp <- seq_len(nr * ncol)
  graphics::layout(matrix(tp, nr, ncol, byrow = TRUE))
  res <- lapply(tp, function(w) {
    if (w <= m) {
      i <- which[[w]]
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
    } else {
      graphics::plot(0, type='n', axes=FALSE)
    }
  })
}
