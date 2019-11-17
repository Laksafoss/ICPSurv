#' Invariant Model Ranking
#'
#' A method for ranking models based on invariance.
#'
#' LONG DESCRIPTION
#'
#' @param M a list of models.
#' @param method bla bla bla
#' @param ... TODO
#'
#' @return DESCRIPTION OF THE RETUNED
#'
#' @export
IMR <- function(M, method = "BS", ...) {

  #mf <- as.list(match.call(expand.dots = T))[-1]  # DO THIS DIFFERENTLY ??
  mc <- sapply(M, function(m) {class(m)})
  class(mc) <- method                             #

  # Find First
  data <- find_data(M[[1]])
  analysis <- construct_analysis(mc, X = data$X, ...)
  res_first <- data.frame(model = data$name,
                          roughness = analysis(data$Y))

  # Find the rest
  res_rest <- lapply(M[-1], function(m) {
    tmp <- find_data(m)
    if (! identical(tmp$X, data$X)) {
      stop("The models are not comparable -- ",
           "the respons variables are not identical")
    }
    res <- data.frame(model = tmp$name,
                      roughness = analysis(tmp$Y))
    return(res)
  })

  res <- rbind(res_first, Reduce(rbind, res_rest))
  class(res) <- c("data.frame","IMR")
  return(res)
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
  if (is.null(dim(m))) {
    stop("When the default find_data method is used each model must be an ",
         "object with exactly two dimensions -- such as a matrix or data.frame")
  }
  if(length(dim(m)) != 2) {
    stop("When the default find_data method is used each model must be an ",
         "object with exactly two dimensions -- such as a matrix or data.frame")
  }
  if (dim(m)[2] > 2) {
    warning("Only first and second column is used from the data")
  }
  name <- paste0("raw data :", names(m)[1], " ~ ", names(m)[2])
  return(list(X = m[,1],
              Y = list(raw = m[,2], std = NULL),
              name = name))
}

#' @rdname find_data
find_data.lm <- function(m) {
  name <- paste0("lm(", format(stats::terms(m)), ")")
  return(list(X = m$fitted.values,
              Y = list(raw = m$residuals, std = stats::rstandard(m)),
              name = name))
}

#' @rdname find_data
find_data.glm <- function(m) {
  name <- paste0("glm(", format(stats::terms(m)),
                 ", family = ",format(m$call$family),")")
  return(list(X = m$fitted.values,
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
  return(list(X = data$time,
              Y = list(raw = data$hazard, std = NULL), # TODO
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
  return(list(X = m$cum[,1],
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
  return(list(X = m$cum[,1],
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
  return(list(X = m$cum[,1],
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
  return(list(X = m$times,
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
#' @param X TODO
#' @param norder TODO
#' @param knots TODO
#' @param L TODO
#' @param W TODO
#' @param ... TODO
#'
#' @return a function that has input \code{Y}.
construct_analysis <- function(mc, X, ...) {
  UseMethod("construct_analysis", mc)
}

#' @rdname construct_analysis
construct_analysis.default <- function(mc, X, ...) {
  stop(class(mc), " is not a recognised method for Invariant Model Ranking")
}

#' @rdname construct_analysis
construct_analysis.BS <- function(mc, X, norder, knots, L, W, ...) {
  if (missing(norder)) {
    norder <- 4
  }
  if (missing(knots)) {
    knots <- stats::quantile(X)
  }
  if (missing(L)) {
    L <- 2
  }
  if (missing(W)) {
    W <- max
  }

  basis <- fda::create.bspline.basis(norder = norder, breaks = knots)
  penmat <- fda::bsplinepen(basis, Lfdobj = L)
  X_norm <- apply(X, 2, function(x) {(x - min(x))/ diff(range(x))})

  foo <- function(Y) {
    Y_norm <- apply(Y$raw, 2, function(x) {(x - min(X)) / diff(range(x))}) # OR Y$std ??
    Bt <- fda::eval.basis(X_norm, basis)
    coef <- solve(crossprod(Bt), crossprod(Bt, Y_norm))
    return(W(diag(t(coef) %*% penmat %*% coef)))
  }
  return(foo)
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
