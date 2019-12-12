
#' Internal Method for Constructing IMR Analysis
#'
#' An internal method used in the \code{\link{IMR}} for construction of the
#' model analysis.
#'
#' LONG DESCRIPTION
#'
#' @param mc TODO
#' @param data TODO
#' @param norder TODO
#' @param knots TODO
#' @param L TODO
#' @param ... TODO
#'
#' @seealso \code{\link{IMR}}
#'
#' @return a function that has input \code{Y}.
construct_analysis <- function(mc, data, ...) {
  UseMethod("construct_analysis", mc)
}

#' @rdname construct_analysis
construct_analysis.default <- function(mc, data, ...) {
  stop("'", class(mc),
       "' is not a recognised method for Invariant Model Ranking.")
}


construct_analysis.rss <- function(mc, data, ...) {
  method <- data.frame("Residual Sum of Squares",
                       row.names = "Method:",
                       fix.empty.names = FALSE)
  maxx <- max(data$time)
  foo <- function(data) {
    coef <- solve(crossprod(data$time), crossprod(data$time, data$fitted))
    Y <- data$fitted / (maxx * c(coef))
    rss <- sum((Y - data$time / maxx)^2)
    return(list("Score" = rss, "normhaz" = Y, "coef" = 1/maxx))
  }
  return(list(analysis = foo,
              method = list(
                method = "rss",
                method_expanded = method,
                critical = "Large values indicate a lack of invariance.",
                time = data$time,
                basis = fda::create.monomial.basis(range(data$time), exponents = 1)
              )))
}

#' @rdname construct_analysis
construct_analysis.BS <- function(mc, data, norder, knots, L, ...) {
  if (is.null(data$fitted)) {
    stop("The 'BS' method requires that 'find_data' outputs a 'fitted' vector.")
  }
  if (missing(norder)) { norder <- 4 }
  if (missing(L)) { L <- 2 }
  if (missing(knots)) {
    knots <- max(min(round(length(data$time) / 100), 30), 3)
  }
  if (is.numeric(knots)) {
    if (length(knots) == 1) { qq <- seq(0, 1, length.out = knots)}
  } else { stop("'knots' must be a numeric of length 1.") }
  method <- data.frame(c("Basis-Splines", norder, length(qq)),
                       row.names = c("Method:", "norder:", "Knots:"),
                       fix.empty.names = FALSE)
  breaks <- stats::quantile(data$time, qq)
  basis <- fda::create.bspline.basis(norder = norder, breaks = breaks)
  penmat <- fda::bsplinepen(basis, Lfdobj = L)
  Bt <- fda::eval.basis(data$time, basis)

  foo <- function(data) {
    if (is.null(data$fitted)) {
      stop("The 'BS' method requires that 'find_data' outputs a 'fitted' vector.")
    }
    Y <- data$fitted /
      (c(solve(crossprod(data$time),
               crossprod(data$time, data$fitted))) *
         max(data$time))
    coef <- solve(crossprod(Bt), crossprod(Bt, Y))
    ans <- max(diag(t(coef) %*% penmat %*% coef))
    return(list("Score" = ans, "normhaz" = Y, "coef" = coef))
  }
  return(list(analysis = foo,

              method = list(
                method = "BS",
                method_expanded = method,
                critical = "Large values indicate a lack of invariance.",
                time = data$time,
                basis = basis,
                knots = breaks
              )))
}

