

#' Creating a Method Object for ICP
#'
#' Creates a \strong{method object} within the \code{\link{ICP}} function.
#'
#' @param method the testing metod used in \code{\link{ICP}}. The method
#'   specified here must be interpretable by the generic function
#'   \code{\link{plausible_predictor_test}}. If nothing is specified
#'   \code{method} is set to \code{EnvirIrrel}
#' @param model the regression model used in \code{\link{ICP}}. The model class
#'   specified here must be interpretable by the generic function
#'   \code{\link{fit_model}}. If \code{model} is not specified it will be set
#'   equal to \code{glm} and \code{family = "gaussian"} will be added to the
#'   \strong{method object}.
#' @param ... further arguments to be passed to the
#'   \code{\link{plausible_predictor_test}} or \code{\link{fit_model}}
#'   functions.
#'
#' @return a list with all the input parameters and their values. This list has
#'   class "Method_obj", \code{method} and \code{model}.
#'
#' @seealso \code{\link{ICP}}, \code{\link{plausible_predictor_test}},
#'   \code{\link{fit_model}} for use of the method object.
#'
#' @examples
#' # the standard method object is of method "EnvirIrrel"
#' ICPSurv:::method_obj()
#'
#' # create a method object suitable for "confidence region" testing
#' # with a poisson glm regression model
#' ICPSurv:::method_obj(method = "CRellipsoid", model = "glm", family = "poisson")
#'
#' # create a method object suitable for "environment irrelevance" testing
#' # with a proportional hazard regression model
#' ICPSurv:::method_obj(method = "EnvirIrrel", model = "ph")
#' # if the survival package is attached the 'survival::' part can be dropped

# AS OF NOW THIS FUNCTION DOES NOT ALLOW FOR NEW
# fit_model and plausible_predictor_test methods !!!!!
method_obj <- function(model = "glm", method = "EnvirIrrel", ...) {
  if ( ! (method %in% c("CR", "EnvirIrrel", "ConstTime"))) {
    stop("'method' must be 'CR', 'EnvirIrrel' or 'ConstTime'") # What is the name of time test ?!?!?
  }
  if ( ! (model %in% c("lm", "glm", "ph", "ah"))) {
    stop("'model' must be 'lm', 'glm', 'ph' or 'ah'")
  }

  mf <- as.list(match.call(expand.dots = T))[-1]
  if (missing(method)) {
    mf <- c(method = "EnvirIrrel", mf)
  }
  if (missing(model)) {
    if (mf$method == "ConstTime") {
      mf <- c(model = "ah", mf)
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
  }

  if (mf$method == "ConstTime") {
    if (mf$model == "glm") {
      stop("A time test for glm is not avalible in this package")
    }
    if (is.null(mf$TimeTest)) {
      mf$TimeTest <- "both"
    }
    if ( ! (mf$TimeTest %in% c("sup", "int", "both"))) {
      stop("'TimeTest' must be 'sup', 'int' or 'both'")
    }
  }

  if (mf$method == "CR") {
    if (is.null(mf$solver)) {
      mf$solver <- "QCLP"
    }
    if ( ! (mf$solver %in% c("QCLP", "pairwise", "marginal"))) {
      stop("The confidence region solver 'solver' must be 'QCLP', 'Pairwise' or 'Marginal'")
    }
    if (mf$solver %in% c("pairwise")) {
      if (is.null(mf$splits)) {
        mf$splits <- "LOO"
      }
      if ( ! (mf$splits %in% c("LOO", "all"))) {
        stop("The pairwise confidence region solve must be either 'LOO' or 'all'")
      }
    }
    class(mf) <- c("method_obj", method, model, mf$solver)
  } else {
    class(mf) <- c("method_obj", method, model)
  }
  return(mf)
}



print.method_obj <- function(x) {
  cat("\nMethod Object:\n")

  cat("\n---\n")
  cat("Statistical Model: ")
  if (x$model == "ph") {
    cat("proportional hazard model")
  } else if (x$model == "ad") {
    cat("additive hazard model")
  } else if (x$model == "glm") {
    cat("generalized linear model")
  } else if (x$model == "lm") {
    cat("linear model")
  } else {cat(x$model)}
  if(!is.null(x$family)) {
    cat("\nFamily:             ")
    if (!is.character(x$family)) {
      cat(x$family$family, "(link = '", x$family$link, "')")
    } else {
      cat(x$family)
    }
  }
  cat("\n---\n")
  cat("Testing Method:     ")
  if (x$method == "CR") {
    cat("Intersecting Confidence Regions (CR)")
    cat("\nOverlap Solver:    ", x$solver)
    if(!is.null(x$split)) {
      cat("\nSplits:             ")
      if (x$split == "LOO") {
        cat("Leave One Out")
      } else {
        cat("all splits")
      }
    }
  } else if (x$method == "EnvirIrrel") {
    cat("Environment Irrelevance (EnvirIrrel)")
  } else if (x$method == "ConstTime") {
    cat("Constant Regression Effect Over Time (ConstTime)")
    cat("\nTime Test:         ", x$TimeTest)
    if (x$TimeTest == "both") {
      cat(" (sup and int)")
    }
  }
  cat("\n---")
}


#' Invariant Causal Prediction
#'
#' A method for finding causal predictors of a target variable.
#'
#' The \code{ICP} function implements tooles for the \emph{Invariant Causal
#' Predictor} methodology described when the target variable is either described
#' by a glm, additive hazard model or multiplicative hazard model. As such the
#' target variable \code{Y} is allowed to be both a vector and a survival object.
#'
#' The \code{ICP} function is essentially a wrapper function for the
#' \code{\link{plausible_predictor_test}} which tests the central null
#' hypothesis of the \emph{Invariance Causal Predictor} methodology, namely
#' \deqn{H_{0,S}: S is an invariant set w.r.t. (X,Y)}
#' where an invariant set is defined as a set of indicies \eqn{S} such that
#' \deqn{(Y^e | X^e_S) = (Y^f | X^f_S)}
#' in distribution for all environments \eqn{e,f}. If the data is time dependent
#' we formulate an analog invariance statement for all time points \eqn{t} and
#' \eqn{s}.
#'
#' In the \code{ICPSurv} package 3 standard methods for the generic function
#' \code{\link{plausible_predictor_test}} has been implemented for testing the
#' null hypothesis above. For further discussion on how to implement new testing
#' methods for the above null hypothesis see
#' \code{\link{plausible_predictor_test}}.
#'
#' To make the \code{ICP} function work correctly with the
#' \code{\link{plausible_predictor_test}} function a \strong{method object}
#' \code{method} must be specfied using the \code{\link{method_obj}} function.
#' The class of the \strong{method object} must correspond to an implemeted
#' method of the generic function \code{\link{plausible_predictor_test}}.
#'
#'
#' @param Y an object describing the response variable. \code{Y} will be passed
#'   to the \code{\link{plausible_predictor_test}} for analyzing.
#' @param X a matrix, data.frame or vector describing the covariates
#' @param E a vector, matrix or data frame describing the environments.
#' @param method a \strong{method object} greated by the
#'   \code{\link{method_obj}} function. A method object is a list describing the
#'   regression model, and the list has class equal to the name of the method
#'   that should be invoked by the generic function
#'   \code{\link{plausible_predictor_test}}. The \strong{method object} list
#'   should always contain an entry named \code{model}, which equals a \code{R}
#'   function (e.g. "glm", "coxph").
#'
#' @param level numerical value between 0 and 1 denoting the significance level
#'   used when testing. If not specified the algirithm will not return an
#'   estimated set of identified causal predictors at level \code{level}.
#'
#' @param fullAnalysis If \code{FALSE} those
#'   \code{\link{plausible_predictor_test}} that find p-values based on
#'   iterative test only test the hypothesis at the specified \code{level}.
#'   Hence \code{level} must be specified if \code{fullAnalysis} is set to
#'   \code{FALSE}.
#'
#'   The inbuilt \code{\link{plausible_predictor_test}} for the
#'   \code{CRrectangle} and \code{CRellipsoid} method are both example of tests,
#'   that find the p-values using iterative tests. So if one of these methods is
#'   used, setting \code{fullAnalysis} to \code{FALSE} will save computational
#'   time. This does however also mean that it is not possible to estimate the
#'   p-values of the individual variables.
#'
#' @param maxNoVariables The maximal number of variables in the tested subsets
#'   of \code{X}. A smaller number saves computational time.
#'
#' @param stopIfEmpty If \code{TRUE} the procedure will stop if the null
#'   hypothesis for the empty set has been accepted. Setting to \code{TRUE} will
#'   save computational time in these cases, but means that an analysis of other
#'   subsets of \code{X} is lost.
#'
#' @param ... additional arguments carried to the
#'   \code{\link{plausibel_predictor_test}}.
#'
#' @return Returns an object of class \code{ICP}. Such an object will at least
#'   contain the following
#'
#'   \item{model.analysis }{a data frame listing the different models tested in
#'   the first column and the found p-values in the second column.}
#'   \item{call}{the matched call.}
#'   \item{level}{the significance level applied to all tests. If not specified
#'   this is simply \code{NULL}.}
#'   \item{method}{the method object used for the model fitting and hypothesis
#'   testing.}
#'
#'   If a \code{level} has been specified then the following will also be part
#'   of the \code{ICP} object
#'
#'   \item{empty.accepted  }{is 0 if the empty model was rejected and 1 if the
#'   empty model was accepted at level \code{level}}
#'   \item{accepted.model}{the estimated causal predictors}
#'
#'   Lastly if \code{fullAnalysis} is set to \code{TRUE} then the \code{ICP}
#'   object will also contain
#'
#'   \item{variable.analysis  }{a data frame with the predictor variables in the
#'   first column and their respective p-values found based on the model
#'   analysis in the second column.}
#'
#'
#'
#' @seealso
#'   \code{\link{print.ICP}} for summaries.
#'   \code{\link{plausible_predictor_test}} for customizing the \code{ICP}
#'   function.
#'
#' @references
#'   Peters, Jonas, Peter Bühlmann, and Nicolai Meinshausen. \emph{Causal
#'   inference by using invariant prediction: identification and confidence
#'   intervals.} Journal of the Royal Statistical Society: Series B (Statistical
#'   Methodology) 78.5 (2016): 947-1012.
#'
#' @examples
#' # ===========================================================================
#' # A simple example with normal distributions and no time dependent variables
#' # First we simulate data.
#' G <- matrix(c(0,1,1,1,0,0,0,0,0,0.4,0,0,0,0,0.4,0,0,0,0,0,0,0,0,0,0),
#'             ncol = 5, byrow = T,
#'             dimnames = list(c("E","X1","X2","X3","Y"), c("E","X1","X2","X3","Y")))
#' sim <- c("sample(c(0,5,10),N,replace=T)", rep("rnorm(N,BETA,1)",4))
#' out <- list(Y = 5, X = 2:4, E = 1)
#' SIM <- sim_from_adj(G, 100, sim, out)
#' # note that Y has parents X1 and X2
#'
#' # create list describing regression model
#' method <- list(model = "glm", family = "gaussian")
#'
#' # analyze using the "Environment Irrelevance test"
#' class(method) <- "EnvirIrrel"
#' ICP(SIM$Y, SIM$X, SIM$E, method = method, level = 0.05)
#'
#' # analyze using the old marginal confidence interval test
#' class(method) <- "CRrectangle"
#' ICP(SIM$Y, SIM$X, SIM$E, method = method, level = 0.05)
#'
#' # analyze using the new "Confidence Region" test with correct ellipsoid CR's
#' class(method) <- "CRellipsoid"
#' ICP(SIM$Y, SIM$X, SIM$E, method = method, level = 0.05)
#' # ---------------------------------------------------------------------------
#'
#'
#'
#' # ===========================================================================
#' # A simple example with poisson distributions and no time dependent variables
#' # First we simulate data.
#' G <- matrix(c(0,1,1,1,0,0,0,0,0,0.4,0,0,0,0,0.4,0,0,0,0,0,0,0,0,0,0),
#'             ncol = 5, byrow = T,
#'             dimnames = list(c("E","X1","X2","X3","Y"), c("E","X1","X2","X3","Y")))
#' sim <- c("sample(c(0,5,10),N,replace=T)", rep("rnorm(N,BETA,1)",3), "rpois(N,exp(BETA))")
#' out <- list(Y = 5, X = 2:4, E = 1)
#' SIM <- sim_from_adj(G, 100, sim, out)
#' # note that Y has parents X1 and X2
#'
#' # create list describing regression model
#' method <- list(model = "glm", family = "poisson")
#'
#' # analyze using the "Environment Irrelevance test"
#' class(method) <- "EnvirIrrel"
#' ICP(SIM$Y, SIM$X, SIM$E, method = method, level = 0.05)
#'
#' # analyze using the old marginal confidence interval test
#' class(model) <- "CRrectangle"
#' ICP(SIM$Y, SIM$X, SIM$E, method = method, level = 0.05)
#'
#' # analyze using the new "Confidence Region" test with correct ellipsoid CR's
#' class(method) <- "CRellipsoid"
#' ICP(SIM$Y, SIM$X, SIM$E, method = method, level = 0.05)
#' # ---------------------------------------------------------------------------

#'
#'
#'
#' @export

# ICP ==========================================================================

ICP <- function(Y, X, E = NULL,
                method = structure(list(model = "glm", family = "gaussian"),
                                  class = "EnvirIrrel"),
                level = NULL, fullAnalysis = T,
                maxNoVariables = 8, stopIfEmpty = FALSE, ...) {
#ICP <- function(Y, X, E = NULL, model = "lm", method = "EnvirIrrel",
#                level = 0.05, fullAnalysis = T,
#                maxNoVariables = 8, ...) {
#
#  method <- method_obj(method = method, model = model, ...)

  if (!is.numeric(maxNoVariables)) {
    stop("'maxNoVariables' must be an integer >= 1")
  }
  if (length(maxNoVariables) != 1){
    stop("'maxNoVariables' must be an integer >= 1")
  }
  if (maxNoVariables <= 0) {
    stop("'maxNoVariables' must be an integer >= 1")
  }
  if (!is.null(level)) {
    if (!is.numeric(level)) {
      stop("'level' must either be NULL or a number strictly between 0 and 1")
    }
    if (length(level) > 1) {
      stop("'level' must either be NULL or a number strictly between 0 and 1'")
    }
    if (level >= 1 | level <= 0) {
      stop("'level' must either be NULL or a number strictly between 0 and 1'")
    }
  }
  if (is.vector(X)) {
    X <- data.frame(X)
  }
  if (is.matrix(X)) {
    X <- data.frame(X)
  }
  if (!is.data.frame(X)) {
    stop("'X' must be a vector, matrix or data frame")
  }
  if (!is.null(E)) {
    E <- interaction(E, drop = T, sep = ":")
    if (length(E) != nrow(X)) {
      stop("'E' and 'X' must have same length / number of rows")
    }
    if (length(levels(E)) <= 1) {
      stop(paste0("there is just one environment ('E'=", levels(E),
                  " for all observations) and the method needs at",
                  " least two environments"))
    }
    if (min(table(E)) <= 2) {
      stop("All environment must have at least three (ideally dozens) ",
           "observations in each environment. The given environment ",
           "vector 'E' has environments with only 1 or 2 observations.")
    }
  }

  # CREATE THE TEST_LIST -------------------------------------------------------
  p <- ncol(X)
  max_no_variables <- min(p, maxNoVariables)
  test_list <- unlist(lapply(seq_len(max_no_variables),
                             function(y) {
                               combn(seq_len(p),
                                     y,
                                     simplify = FALSE)}),
                      recursive = FALSE)
  test_list <- unlist(list(0, test_list), recursive = FALSE)

  # LOOP THROUGH TEST_LIST -----------------------------------------------------

  # initialize counter and result matrix
  res_model <- data.frame(models = character(), pval = numeric())
  s <- 1

  # loop through sets
  while (s < length(test_list) + 1) {
    # Find current subset and analyze
    X_sub <- X[ , test_list[[s]], drop = F]
    pval <- plausible_predictor_test(method = method, Y = Y, X = X_sub, E = E,
                                     level = level, fullAnalysis = fullAnalysis,
                                     ...)
    nam <- if (s == 1) {"Empty"} else {paste0(colnames(X_sub), collapse = " + ")}
    res_model <- rbind(res_model, data.frame(model = nam, pval = pval))

    # Clean up test_list if posible
    if (pval >= ifelse(is.null(level), 1, level)) {
      if (nrow(res_model) == 1 && stopIfEmpty) {
        test_list <- test_list[[1]]
      } else {
        remove <- sapply(seq_along(test_list), function(i) {
          if (i <= s) {
            F
          } else {
            all(test_list[[s]] %in% test_list[[i]])
          }
        })
        test_list <- test_list[!remove]
      }
    }
    s <- s + 1
  }



  # SUMMARIZE RESULTS ----------------------------------------------------------
  call <- match.call()
  output <- structure(
    list(
      call = call,
      level = level,
      model.analysis = res_model,
      method = method
      ),
    class = "ICP"
  )


  # FIND ACCEPTED MODEL & VARIABLE 'P-VALUES' ------------------------------------
  if (!is.null(level)) {
    # find accepted set
    output$empty.accepted <- 0
    plausible <- which(res_model$pval >= level)
    if (length(plausible) == 0) {
      res_acc <- "NO SUBSET OF PREDICTORS WAS ACCEPTED"
    } else if (1 %in% plausible) {
      res_acc <- "Empty"
      output$empty.accepted <- 1
    } else {
      acc <- Reduce(intersect, test_list[plausible])
      if (length(acc) == 0) {
        res_acc <- "Empty"
      } else {
        res_acc <- paste0(colnames(X)[acc], collapse = " + ")
      }
    }
    output$accepted.model <- res_acc

    # finde variable 'p-values'
    if (fullAnalysis) {
      if (res_acc == "Empty") { # TODO or is it !any(res_model$pval[-1] >= level)
        pval_var <- rep(1, length.out = ncol(X))
      } else {
        # TODO test correctness of this code !!!!!
        pval_var <- sapply(seq_len(ncol(X)), function(i) {
          contained <- sapply(test_list, function(j) {i %in% j})
          max(res_model$pval[!contained])
        })
      }
      res_var <- data.frame(variables = colnames(X), pval = pval_var)
      output$variable.analysis <- res_var
    }
  }


  return(output)
}


