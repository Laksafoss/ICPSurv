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
#'   \code{\link{construct_anaysis.BS}}.
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
  mc <- sapply(M, function(m) {class(m)})
  class(mc) <- method

  # Find First
  data <- test_find_data(find_data(M[[1]]))
  constructed <- test_constructor(construct_analysis(mc, data = data, ...))
  analysis <- constructed$analysis
  res_first <- data.frame("Model" = data$name,
                          "Score" = analysis(data))

  # Find the rest
  res_rest <- lapply(M[-1], function(m) {
    tmp <- test_find_data(find_data(m))
    if (! identical(tmp$target, data$target)) {
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
    res <- data.frame("Model" = tmp$name,
                      "Score" = analysis(tmp))
    return(res)
  })

  res <- rbind(res_first, Reduce(rbind, res_rest))
  if (!is.null(names(M))) {
    rownames(res) <- names(M)
  }
  tmp <- structure(list(ranking = res,
                        method = constructed$method,
                        critical = constructed$critical,
                        call = match.call()),
                   class = "IMR")

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
  if (is.null(data$target)) { stop(paste0(err, "\n'target' is missing.")) }
  return(data)
}

#' @export
print.IMR <- function(x, ...) {
  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")
  cat("\nSummary of Ranking Method:")
  print(x$method, right = FALSE)

  cat("\nInvariance Ranking:\n")
  x$ranking[,2] <- signif(x$ranking[,2], digits = 4)
  print(x$ranking)
  cat("---\n", x$critical)
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
  return(list(target = m$model[,1],
              X = m$fitted.values,
              Y = list(raw = m$residuals, std = stats::rstandard(m)),
              name = name))
}

#' @rdname find_data
find_data.glm <- function(m) {
  name <- paste0("glm(", format(stats::terms(m)),
                 ", family = ",format(m$call$family),")")
  return(list(target = m$model[,1],
              X = m$fitted.values,
              Y = list(raw = m$residuals, std = stats::rstandard(m)),
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
  return(list(target = data$time,
              X = data$time,
              Y = list(raw = data$hazard, std = NULL), # TODO
              beta = m$coefficients,
              var.beta = m$var,
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
  return(list(target = m$cum[,1],
              X = m$cum[,1],
              Y = list(raw = m$cum[,-1], std = NULL), # TODO
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
  return(list(target = m$cum[,1],
              X = m$cum[,1],
              Y = list(raw = m$cum[,-1], std = NULL), #TODO
              name = name))
}

#' @rdname find_data
find_data.cox.aalen <- function(m) {
  X <- m$cum[,1]
  Y <- m$cum[,-1]

  form <- gsub("survival::", "", format(m$call$formula))
  if (!is.null(m$call$id)) {
    form <- paste0(form, paste0(", id = ", m$call$id))
  }
  if (!is.null(m$call$clusters)) {
    form <- paste0(form, paste0(", clusters = ", m$call$clusters))
  }
  name <- paste0("cox.aalen(", form, ")")
  return(list(target = m$cum[,1],
              X = m$cum[,1],
              Y = list(raw = m$cum[,-1], std = NULL), # TODO
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
  return(list(target = m$times,
              X = m$times,
              Y = list(raw = Y, std = NULL), # TODO
              name = name))
}


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
#' @param W TODO
#' @param ... TODO
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


construct_analysis.CR <- function(mc, data, ...) {
  # TODO
}

#' @rdname construct_analysis
construct_analysis.BS <- function(mc, data, norder, knots, L, W, ...) {
  if (missing(norder)) { norder <- 4 }
  if (missing(L)) { L <- 2 }
  if (missing(W)) { W <- max }
  if (missing(knots)) {
    knots <- max(min(round(length(data$target) / 100), 30), 3)
  }
  if (is.numeric(knots)) {
    if (length(knots) == 1) {
      qq <- seq(0, 1, length.out = knots)
    }
  } else {
    stop("'knots' must be a numeric of length 1.")
  }
  nk <- length(qq)
  if (nk <= 3) {
    knots_str <- paste0(nk, " (", paste0(round(qq, 2), collapse = ","), ")")
  } else {
    knots_str <- paste0(nk,
                        " (",
                        paste0(round(qq[c(1,2)], 2), collapse = ","),
                        ",...,",
                        paste0(round(qq[nk],2)),
                        ")")
  }
  method <- data.frame(c("Basis-Splines",
                         norder,
                         knots_str),
                       row.names = c("Method:", "norder:", "Knots:"),
                       fix.empty.names = FALSE)

  foo <- function(data) {
    if (is.null(data$X)) {
      stop("The 'BS' method requires that 'find_data' outputs an 'X'.")
    }
    if (is.null(data$Y$raw)) {
      stop("The 'BS' method requres that 'find_data' outputs a 'Y'.")
    }
    breaks <- stats::quantile(data$X, qq)
    basis <- fda::create.bspline.basis(norder = norder, breaks = breaks)
    penmat <- fda::bsplinepen(basis, Lfdobj = L)
    Y_norm <- normalize(data$Y$raw)
    Bt <- fda::eval.basis(data$X, basis)
    coef <- solve(crossprod(Bt), crossprod(Bt, Y_norm))
    return(W(diag(t(coef) %*% penmat %*% coef)))
  }
  return(list(analysis = foo,
              method = method,
              critical = "Large values indicate a lack of invariance."))
}

normalize <- function(x) {
  if (length(dim(x)) == 2) {
    return(apply(x, 2, function(m) {(m - min(m)) / diff(range(m))}))
  } else {
    return((x - min(x)) / diff(range(x)))
  }
}

#' @rdname construct_analysis
construct_analysis.KS <- function(mc, X, W, ...) {
  if (missing(W)) {
    W <- max
  }
  if (all(mc %in% c("coxph", "aalen", "timecox", "cox.aalen", "aareg"))) {
    foo <- function(Y) {
      B <- solve(crossprod(X), crossprod(X, Y$val))
      resid <- Y$val - X %*% B
      n <- dim(Y$val)[1]
      sd <- sqrt(n * Y$var[n,])
      resid <- t(t(resid) / sd) # standerdize
      return(W(apply(resid, 2, function(r) {max(r)})))
    }
  } else {
    return(1) # TODO
  }
}

#' @rdname construct_analysis
construct_analysis.CM <- function(mc, X, ...) {
  return(1) # TODO
}
