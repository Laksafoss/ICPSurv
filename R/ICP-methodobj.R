
#' Creating a Method Object for ICP
#'
#' Creates a \strong{method object} within the \code{\link{ICP}} function.
#'
#' @param model the regression model used in \code{\link{ICP}}. The model class
#'   specified here must be interpretable by the generic function
#'   \code{\link{fit_model}} or \code{\link{fit_nonparam_model}}. If
#'   \code{model} is not specified it will be set equal to \code{glm} and
#'   \code{family = "gaussian"} will be added to the \strong{method object}.
#' @param method the testing metod used in \code{\link{ICP}}. The method
#'   specified here must be interpretable by the generic function
#'   \code{\link{plausible_predictor_test}}. If nothing is specified
#'   \code{method} is set to \code{EnvirRel}.
#' @param ... further arguments to be passed to the
#'   \code{\link{plausible_predictor_test}}, \code{\link{fit_model}} or
#'   \code{\link{fit_nonparam_model}} functions.
#'
#'   If \code{method} is set to "\code{\link[=plausible_predictor_test]{CR}}"
#'   then a \code{solver} can be specified. The standard solver is
#'   \code{\link[=plausible_predictor_test]{QC}}. However, if \code{method} is
#'   "\code{\link[=plausible_predictor_test]{TimeVar}}" then \code{n.sim} can
#'   be set to any integer larger then 50.
#'
#'   If \code{model} is \code{glm} then a \code{family} should be specified
#'   (see \code{\link[stats]{family}}). If a \code{family} is not specified it
#'   is set to "\code{gaussian}".
#'
#' @param x a \code{method_obj} for printing.
#' @param exclude variables from the \code{method_obj} \code{x} which is not to
#'   be included in the printout.
#'
#' @return The \code{method_obj} function returns a list with all the input
#'   parameters and their values. This list has class "method_obj",
#'   \code{method} and \code{model}. If \code{method} is \code{CR} then the
#'   output list also has class equal to \code{solver}. If, however, the
#'   \code{method} is \code{TimeVar} then the output also has class
#'   \code{nonparamtest}.
#'
#' @seealso \code{\link{ICP}}, \code{\link{plausible_predictor_test}},
#'   \code{\link{fit_model}} or \code{\link{fit_nonparam_model}}  for use of the
#'   method object.
#'
#' @examples
#' # The standard method object is a gaussian glm with method EnvirIrrel
#' ICPSurv:::method_obj()
#'
#' # A method object for Intersecting Confidence Region analysis of linear models
#' ICPSurv:::method_obj(model = "lm", method = "CR", solver = "QC", fullAnalysis = TRUE)
#' # Here 'fullAnalysis' has been set to TRUE which will ensure that QCLP method
#' # does not use rectangle approximations as a first step
#'
#' @keywords internal

method_obj <- function(model = "glm", method = "EnvirRel", ...) {

  mf <- c(as.list(environment()), list(...))

  if (is.null(mf$method)) {
    mf <- c(method = "EnvirRel", mf)
  }

  if (is.null(mf$model)) {
    if (mf$method == "TimeVar") {
      mf <- c(model = "hazard", mf)
    } else {
      mf <- c(model = "glm", mf)
    }
  }

  if (mf$model == "glm") {
    if(is.null(mf$family)) {
      mf$family <- "gaussian"
    }
    if (is.call(mf$family)) {
      mf$family <- eval(mf$family)
    }
  } else if (mf$model == "hazard") {
    if (is.null(mf$dist)) {
      if (mf$dist$dist != "extreme"){
        stop("The basic parent distribution - 'dist$dist' - must be ",
             "'extreme' for 'hazard' models")
      } else {
        stop("'dist' must be specified for 'hazard' models")
      }
    }
  }

  if (mf$method == "TimeVar") {
    if (mf$model == "glm") {
      stop("A 'TimeVar' for glm is not yet avalible in this package")
    }
    if (is.null(mf$nonparamtest)) {
      mf$nonparamtest <- "both"
    }
    if (is.null(mf$n.sim)) {
      mf$n.sim <- 50
    } else if (mf$n.sim < 50) {
      stop("A 'nonparamtest' with n.sim less then 50 not yet implemented")
    }
    if ( ! (mf$nonparamtest %in% c("sup", "int", "both"))) {
      stop("'nonparamtest' must be 'sup', 'int' or 'both'")
    }
    class(mf) <- c("method_obj", method, model, mf$nonparamtest)
  } else if (mf$method == "CR") {
    if (is.null(mf$solver)) {
      mf$solver <- "QC"
    }
    if (is.null(mf$splits)) {
      mf$splits <- "all"
    }
    if ( ! (mf$splits %in% c("LOO", "all"))) {
      stop("The pairwise confidence region solve must be either 'LOO' or 'all'")
    }
    class(mf) <- c("method_obj", method, model, mf$solver)
  } else {
    class(mf) <- c("method_obj", method, model)
  }
  mf$call <- match.call()
  return(mf)
}


#' @rdname method_obj
#' @export
print.method_obj <- function(x, exclude = c("id", "tol", "call", "dist"), ...) {

  # Make all variable names nice and put in correct order
  nam <- names(x)
  nam <- nam[! nam %in% exclude]
  order <- c("model", "family", "dist",
             "method", "solver", "nonparamtest", "splits", "n.sim",
             "Bonferroni", "FullAnalysis", "level")
  sort <- match(order, nam, nomatch = 0)
  notsorted <- seq_along(nam)[-sort]
  nam <- nam[c(sort, notsorted)]
  nam <- unname(sapply(nam, function(x) {
    if (x == "n.sim") {
      return("n.sim: ")
    } else {
      substr(x, 1, 1) <- toupper(substr(x, 1, 1))
      return(paste0(x, ": "))
    }
  }))
  w <- x[c(sort, notsorted)]

  # Make all variable values nice and printable
  w$model <- switch(
    w$model,
    ph = "Proportional Hazard Model",
    ah = "Additive hazard Model",
    glm = "Generalized Linear Model",
    lm = "Linear Model",
    hazard = "Hazard Model",
    w$model)
  w$method <- switch(
    w$method,
    "EnvirRel" = "Environment Irrelevance Test",
    "CR" = "Intersecting Confidence Regions Test",
    "TimeVar" = "Non-parametric Time Variance Test",
    w$method)
  if (!is.null(w$solver)) {
    w$solver <- switch(
      w$solver,
      "QC" = "Quadratically Constrained Program",
      "pairwise" = "Pairwise Intersections",
      "marginal" = "Marginal Overlaps",
      w$solver
    )
  }
  if (!is.null(w$splits)) {
    w$splits <- switch(
      w$splits,
      "LOO" = "Leave One Out CR Constructor",
      "all" = "CR For All Environments",
      w$splits)
  }
  if (!is.null(w$nonparamtest)) {
    w$nonparamtest <- switch(
      w$nonparamtest,
      "sup" = "Kolmogorov-Smirnov Test",
      "int" = "Cramer-von Mises Test",
      "both" = "Both Kolmogorov-Smirnov and Cramer-von Mises Tests",
      w$nonparamtest)
  }

  # Ensure everything is character
  if(!is.null(w$dist)) {
    w$dist <- ifelse(is.character(w$dist), w$dist, w$dist$name)
  }

  w <- lapply(seq_along(w), function(i) {
    if (is.character(w[[i]])) {
      return(w[[i]])
    } else {
      index <- which(names(x$call) == names(w)[i])
      if (length(index) == 1) {
        return(deparse(x$call[[index]]))
      } else {
        return(deparse(w[[i]])[1])
      }
    }
  })

  # print
  dd <- data.frame(unlist(w), row.names = nam, fix.empty.names = FALSE)
  print(dd, right = FALSE)

  invisible(x)
}
