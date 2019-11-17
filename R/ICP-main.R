

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
#'   \code{method} is set to \code{EnvirIrrel}.
#' @param ... further arguments to be passed to the
#'   \code{\link{plausible_predictor_test}}, \code{\link{fit_model}} or
#'   \code{\link{fit_nonparam_model}} functions.
#'
#'   If \code{method} is set to "\code{\link[=plausible_predictor_test]{CR}}"
#'   then a \code{solver} can be specified. The standard solver is
#'   \code{\link[=plausible_predictor_test]{QCLP}}. However, if \code{method} is
#'   "\code{\link[=plausible_predictor_test]{ConstTime}}" then \code{n.sim} can
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
#'   \code{method} is \code{nonparam} then the output also has class
#'   \code{nonparamtest}.
#'
#' @seealso \code{\link{ICP}}, \code{\link{plausible_predictor_test}},
#'   \code{\link{fit_model}} or \code{\link{fit_nonparam_model}}  for use of the
#'   method object.
#'
#' @examples
#' # The standard method object is a gaussian glm with method EnvirIrrel
#' method_obj()
#'
#' # A method object for Intersecting Confidence Region analysis of linear models
#' method_obj(model = "lm", method = "CR", solver = "QCLP", fullAnalysis = TRUE)
#' # Here 'fullAnalysis' has been set to TRUE which will ensure that QCLP method
#' # does not use rectangle approximations as a first step
#'
#' @export

method_obj <- function(model = "glm", method = "EnvirIrrel", ...) {

  mf <- c(as.list(environment()), list(...))

  if (is.null(mf$method)) {
    mf <- c(method = "EnvirIrrel", mf)
  }

  if (is.null(mf$model)) {
    if (mf$method == "nonparam") {
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
    if (is.null(mf$link)) {
      if (!is.null(mf$dist$name)) {
        mf$link <- mf$dist$name
      } else {
        mf$link <- "proportional"
      }
    }
    if (is.null(mf$dist)) {
      if (mf$link %in% c("proportional", "log")) {
        mf$dist <- survival::survreg.distributions$exponential
      } else if (mf$link %in% c("additive", "identity")) {
        mf$dist <- survival::survreg.distributions$exponential
        mf$dist$trans <- mf$dist$itrans <- function(y) y
        mf$dist$dtrans <- function(y) rep(1, length(y))
      } else if (mf$dist$dist != "extreme"){
        stop("The basic parent distribution - 'dist$dist' - must be ",
             "'extreme' for 'hazard' models")
      } else {
        stop("Either 'link' or 'dist' must be specified for 'hazard' models")
      }
    }
  }

  if (mf$method == "nonparam") {
    if (mf$model == "glm") {
      stop("A 'nonparamtest' for glm is not yet avalible in this package")
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
      mf$solver <- "QCLP"
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
  order <- c("model", "family", "link", "dist",
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
    "EnvirIrrel" = "Environment Irrelevance Test",
    "CR" = "Intersecting Confidence Regions Test",
    "nonparam" = "Non-parametric Test",
    w$method)
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


#' Invariant Causal Prediction
#'
#' A method for finding causal predictors of a target variable.
#'
#' The \code{ICP} function implements tooles for the \emph{Invariant Causal
#' Predictor} methodology when the target variable may be described by a lm,
#' glm or hazard model. To allow for survival type models the target variable
#' \code{Y} is allowed to be both a vector and a survival object.
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
#' In the \code{ICPSurv} package three standard methods for the generic function
#' \code{\link{plausible_predictor_test}} has been implemented for testing the
#' null hypothesis above. These are the
#' \code{\link[=plausible_predictor_test]{EnvirIrrel}},
#' \code{\link[=plausible_predictor_test]{CR}} and
#' \code{\link[=plausible_predictor_test]{nonparam}} method. For further
#' discussion on how to implement new testing methods for the above null
#' hypothesis see \code{\link{plausible_predictor_test}}.
#'
#'
#' @param Y an object describing the response variable. \code{Y} will be passed
#'   to the \code{\link{plausible_predictor_test}} for analyzing.
#' @param X a matrix, data.frame or vector describing the covariates
#' @param E a vector, matrix or data frame describing the environments.
#' @param model a character which must correspont to a \code{\link{fit_model}}
#'   or \code{\link{fit_nonparam_model}} method. In this package the model types
#'   implemented are \code{\link[=fit_model.lm]{lm}},
#'   \code{\link[=fit_model.glm]{glm}}, \code{\link[=fit_model.ph]{ph}},
#'   \code{\link[=fit_model.ah]{ah} and \code{\link[=fit_model.hazard]{hazard}},
#'   but the user may specify new methods for the generic functions
#'   \code{\link{fit_model}} or \code{\link{fit_nonparam_model}} for more
#'   \code{model} options.
#'
#'   If \code{model} is equal to \code{glm} then a \code{family} must be
#'   specified in \code{...}. For more details on valid family inputs see
#'   \code{\ink[stats]{family}}. If \code{family} is not specified, then it is
#'   set to "gaussian". Note that (\code{model = glm, \code{family = "gaussian"}})
#'   is synonymous with (\code{model = lm}).
#'
#'   Is \code{model} is equal to \code{hazard} a \code{link} or \code{dist} must
#'   be specified in \code{...}. A user specified \code{dist} allowes for the
#'   most flexibility. The format of a \code{dist} list is descibed in
#'   \code{\link[survival]{survreg.distributions}}. Possible \code{link} values
#'   include \code{log}, \code{identity}, \code{proportional} and \code{additive}.
#'   Howbeit, \code{link} equal to \code{log} or \code{proportional} is
#'   synonymous with \code{model} equal to \code{ph}, and \code{link} equal to
#'   \code{identity} or \code{additive} is synonymous with \code{model} equal
#'   to \code{ah}.
#' @param method a character which must correspond to a
#'   \code{\link{plausible_predictor_test}} method. In this package the methods
#'   \code{\link[=plausible_predictor_test.EnvirIrrel]{EnvirIrrel}},
#'   \code{\link[=plausible_predictor_test.CR]{CR}} and
#'   \code{\link[=plausible_predictor_test.nonparam]{nonparam}} are already
#'   defined. If the user defines a new S3 method for the
#'   \code{\link{plausible_predictor_test}} then this new method becomes a valid
#'   input for \code{method} here. For more details on how to implement your own
#'   method for use in the \code{ICP} function see
#'   \code{\link{plausible_predictor_test}}.
#' @param level numerical value between 0 and 1 denoting the significance level
#'   used when testing. If not specified the algirithm will only calculate the
#'   p-values of the null hypothesis \eqn{H_{0,S}} and draw no conclusions based
#'   on these values.
#' @param maxNoVariables The maximal number of variables in the tested subsets
#'   of \code{X}. A smaller number saves computational time, but also
#'   corresponds to assuming that the target variable \code{Y} has at most
#'   \code{maxNoVariables} identifiable predictors.
#' @param fullAnalysis if \code{FALSE} different computations time saving
#'   options will be enabled and as a result many methods will not be able to
#'   return proper p-values of the null hypothesis \eqn{H_{0,S}}. In turn this
#'   means that the method will no longer be able to estime the significanse of
#'   a single variable. \code{fullAnalysis = FALSE} will result in different
#'   behavior for different \code{methods}, so for a more detailed discussion of
#'   concrete effects see the help page for the chooosen \code{method}.
#' @param ... additional arguments carried to the
#'   \code{\link{plausible_predictor_test}}.
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
#'   Jonas Peters, Peter Bühlmann, and Nicolai Meinshausen. \emph{Causal
#'   inference by using invariant prediction: identification and confidence
#'   intervals.} Journal of the Royal Statistical Society: Series B (Statistical
#'   Methodology) 78.5 (2016): 947-1012.
#'
#' @examples
#' # ===========================================================================
#' # An example with normal distributions
#' n <- 500
#' E <- sample(5L, n, replace = TRUE)
#' X <- rnorm(n, E, 1)
#' Y <- rnorm(n, X, 1)
#' ICP(Y, X, E) # Environment Irrelevance Test
#' ICP(Y, X, E, method = "CR") # Confidence Region Test
#'
#' # ===========================================================================
#' # An example with a poisson distribution
#' Y <- rpois(n, exp(X))
#' # Environment Irrelevance Test
#' ICP(Y, X, E, model = "glm", family = "poisson")
#' # Intersecting Confidence Region Test
#' ICP(Y, X, E, model = "glm", family = "poisson", method = "CR")
#'
#' # ===========================================================================
#' # An example with right censored survival times
#' Y <- rexp(n, exp(- 1.5 *X))
#' C <- rexp(n, exp(- 2.0 * X))
#' time <- pmin(Y, C)
#' status <- time == Y
#' # Environment Irrelevance Test
#' ICP(survival::Surv(time, status), X, E, model = "ph")
#' # Intersecting Confidence Regions Test
#' ICP(survival::Surv(time, status), X, E, model = "ph", method = "CR")
#' # Non-parametric Constant Effect Over Time Test
#' ICP(survival::Surv(time, status), X, E, model = "ph", method = "nonparam")
#'
#' @export

# ICP ==========================================================================
ICP <- function(Y, X, E = NULL, model = "lm", method = "EnvirIrrel",
                level = 0.05, maxNoVariables = 8, fullAnalysis = TRUE, ...) {


  method <- method_obj(model = model, method = method, ...)

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
    X <- tryCatch(data.frame(X), error = function(x) {0})
    if (identical(X, 0)) {
      stop("When 'X' is a vector it must be coercible to data.frame.")
    }
  }
  if (is.matrix(X)) {
    X <- data.frame(X)
  }
  if (!is.data.frame(X)) {
    stop("'X' must be a vector, matrix or data frame")
  }
  if (!is.null(E)) {

    if (is.matrix(E)) {
      E <- data.frame(E)
    }
    if (is.data.frame(E)) {
      E <- interaction(E, drop = T, sep = ":")
    }
    if (is.list(E)) {
      stop("'E' must be vector, matrix, data frame or factor")
    }
    if (is.vector(E)) {
      E <- as.factor(E)
    }
    if (! is.factor(E)) {
      stop("'E' must be vector, matrix, data frame or factor")
    }
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
                               utils::combn(seq_len(p),
                                            y,
                                            simplify = FALSE)}),
                      recursive = FALSE)
  test_list <- unlist(list(0, test_list), recursive = FALSE)

  # LOOP THROUGH TEST_LIST -----------------------------------------------------

  # initialize
  smallest_possible_pvalue <- .Machine$double.eps
  res_model <- data.frame(models = character(), pval = numeric())
  s <- 1L

  # loop through sets
  while (s <= length(test_list)) {
    # Find current subset and analyze
    X_sub <- X[ , test_list[[s]], drop = F]
    pval <- plausible_predictor_test(method = method,
                                     Y = Y, X = X_sub, E = E,
                                     level = level,
                                     fullAnalysis = fullAnalysis, ...)
    nam <- ifelse(s == 1L, "Empty", paste0(colnames(X_sub), collapse = " + "))
    res_model <- rbind(res_model, data.frame(model = nam, pval = pval))

    # Clean up test_list if posible
    if (! fullAnalysis) {
      if (pval >= ifelse(is.null(level), 1, level)) {
        if (nrow(res_model) == 1) {
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
    }
    s <- s + 1L
  }

  # SUMMARIZE RESULTS ----------------------------------------------------------
  call <- match.call()
  output <- structure(
    list(
      call = call,
      level = level,
      smallest_possible_pvalue = smallest_possible_pvalue,
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
      if (!any(res_model$pval >= level)) { # TODO or is it res_acc == "Empty"
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


