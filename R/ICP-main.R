#' Invariant Causal Prediction
#'
#' A method for finding causal predictors of a target variable described by
#' either a linear, generalized linear or hazard model. The methodology uses
#' heterogeneous data to make causal inference.
#'
#' The \code{ICP} function implements different concrete methods within the
#' methodology of \emph{invariant Causal Predictions} which was first desriced
#' in Peters et al. (2016) (see references below). This implementation of
#' \emph{invariant Causal Predictions} is well suited when the distribution of
#' the target variable may be described by a linear model, generalized linear
#' model or hazard model. There are three different methods for testing
#' invariance implemented in \code{ICP} - \code{EnvirRel}, \code{CR} and
#' \code{TimeVar} - and they are each given a description below under "The
#' Invarince Test Methods".
#'
#' As input the \code{ICP} function takes a target variable \code{Y} which is
#' either a numeric vector or a \code{\link[survival:Surv]{Survival}} object, a
#' matrix or data.frame of covariates \code{X} and possibly - depending on the
#' method - a vector of environments \code{E}. The \code{ICP} function computes
#' a p-value of the following family of null hypotheses:
#' \ifelse{html}{\out{<center>H<sub>0,S</sub> : (Y<sub>i</sub> &#8739 X<sub>i</sub><sup>S</sup> = x) = (Y<sub>j</sub> &#8739 X<sub>j</sub><sup>S</sup> = x) in distribution for all indices i, j and x.</center>}}{\deqn{H_{0,S} :\ (Y_i | X_i^S = x) = (Y_j | X_j^S = s)\ \textrm{ in  distribution for all indices } i,j \textrm{ and } x.}}
#' for every \ifelse{html}{\out{S&#8838{1,...,p}}}{\eqn{S\subseteq \{1,...,p\}}} (where
#' we have assumed that \code{X} encodes \code{p} covariates). The results of
#' these hypothesis tests may be found in \code{model.analysis}.
#'
#' If \code{level} is specified (a subset of) the causal predictors is estimated
#' using the formula (see Peters et al. (2016) for details):
#' \ifelse{html}{\out{<center>A = &#8745<sub>{S: H<sub>0,S</sub> accepted}</sub> S.</center>}}{\deqn{A=\cap_{\{S:\ H_{0,S} \textrm{ accepted}\}} \ S.}}
#' The set \code{A} is outputted under the name \code{accepted.model}. This
#' computation is done by the function \code{\link{model_analysis}}, which is
#' also a function in its own right.
#'
#' Moreover, if both \code{level} is specified and \code{fullAnalysis = TRUE}
#' then the function \code{\link{variable_analysis}} will calculate the
#' significance of each individual variable in \code{X}. This significance table
#' is returned under the name \code{variable.analysis}.
#'
#' The \code{gof} parameter protects against making statements when the model is
#' obviously not suitable for the data. If no model reaches the threshold
#' \code{gof} significance level, i.e. the p-values for
#' \ifelse{html}{\out{(H<sub>0,S</sub>)}}{\eqn{(H_{0,S})}} are all smaller then
#' \code{gof}, we report that there is no evidence for individual variables, as
#' there is no evidence for an invariant set.
#'
#' \strong{The Invarince Test Methods}
#'
#' Three different invariance test methods have been implemented:
#'
#' \code{method = "EnvirRel"} : The invariance test method of \emph{Environment
#' Relevance} is the standard method and can be applied data from to all
#' \code{model} types (lm, glm & hazard). This method requires environments
#' \code{E} as input.
#'
#' \code{method = "CR"} : The invariance test method of \emph{Intersecting
#' Confidence Regions} can be applied to data from to all \code{model} types
#' (lm, glm & hazard). This method requires environments \code{E} as input.
#' Moreover, a solution within the \code{CR} method framework may be found in
#' tree different ways: The standard is \code{solver = "QC"}, which is ususally
#' also the slowest solver. If computational time is an issue the user may need
#' to use the approximate solvers \code{solver = "pairwise"} or
#' \code{solver = "marginal"}.
#'
#' \code{method = "TimeVar"} : The invariance test method of \emph{Time
#' Variability} can only be applied to data from "\code{ph}" or "\code{ah}" type
#' models. This method \emph{does not} require environment information, as it
#' uses time as environment. The "\code{TimeVar}" method has three different
#' concrete \code{nonparamtest}s: a Kolmogorov–Smirnov test type test denotes
#' "\code{sup}", a Cramér–von Mises criterion type test denoted "\code{int}",
#' or simply both tests denoted "\code{test}".
#'
#'
#'
#'
#' @param Y The response or target variable of interest. Either a numeric vector
#'   or \code{\link[survival:Surv]{survival}} object.
#' @param X A matrix (or data frame) with the predictor variables.
#' @param E Indicator of the experiment or the intervention type an observation
#'   belongs to. Can be a vector of the same length as \code{Y} with at least
#'   two unique values.
#' @param model A character indicating how to model the ditribution of the
#'   target variable given covariates. Possible choices are
#'   \itemize{
#'     \item \code{lm} : Linear Model.
#'     \item \code{glm} : Generalized Linear Model. When using
#'       \code{model = "glm"} a \code{\link[stats]{family}} must also be
#'       specified in function options.
#'     \item \code{ph} : Proportional Hazard Model.
#'     \item \code{ah} : Additive Hazard Model.
#'     \item \code{hazard} : Hazad Model. When using \code{model = "hazard"}
#'       a \code{\link[survival:survreg.distributions]{dist}} must be
#'       specified.
#'   }
#' @param method A character indicating which method to use. Possible values are
#'   \itemize{
#'     \item \code{EnvirRel} : Environment Relevance Test.
#'     \item \code{CR} : Intersecting Confidence Regions Test. Using this method
#'       the user may also specify the \code{solver} (see detailes about "The
#'       Invariance Test Methods").
#'     \item \code{TimeVar} : Time Variations Test. Using this method the user
#'       may also specify \code{nonparamtest} (see detailes about "The
#'       Invariance Test Methods").
#'   }
#'   See detailes for more guidence on methods.
#' @param level Numerical value between 0 and 1 denoting the significance level
#'   used when testing. If not specified the algorithm will only calculate the
#'   p-values of the null hypotheses
#'   \ifelse{html}{\out{(H<sub>0,S</sub>>)}}{\eqn{(H_{0,S})}} and draw no
#'   conclusions based on these values.
#' @param gof If no set of variables (including the empty set) leads to a
#'   p-value larger than the goodness-of-fit cutoff \code{gof}, the whole model
#'   will be rejected. If the model is correct, this will happen with a
#'   probability of gof. This option protects again making statements when the
#'   model is obviously not suitable for the data.
#' @param maxNoVariables The maximal number of variables to pre-select (choosing
#'   smaller values saves computational resources but increases approximation
#'   error).
#' @param fullAnalysis If \code{TRUE} p-values for all null hypotheses will be
#'   found. If \code{FALSE} it will often be possible to save computation time:
#'   this depends on the method.
#' @param progress If \code{TRUE} a progress bar will be printed.
#' @param ... Additional arguments carried to the lower level functions.
#'
#'
#' @return The \code{ICP} function returns an object of \code{\link[base]{class}}
#'   \code{ICP}. Such an object will contain the following
#'
#'   \item{model.analysis }{ A data frame listing the different models tested in
#'   the first column and the found p-values in the second column.}
#'   \item{call }{ The matched call.}
#'   \item{level }{ The significance level. If not specified this is \code{NULL}.}
#'   \item{method }{ The method object used for the model fitting and hypothesis
#'     testing.}
#'   \item{accepted.model }{ The estimated causal predictors. Only returned if
#'     \code{level} is specified in the input.}
#'   \item{empty.message }{ If the empty set is returned as \code{accepted.model}
#'     then \code{empty.message} will give detailes. Only returned if
#'     \code{level} is specified in the input.}
#'   \item{variable.analysis  }{ A data.frame with each predictor variables
#'     significance as causal predictors. Will only be returned if
#'     \code{fullAnalysis = TRUE} in the input options.}
#'
#'
#'
#' @seealso
#'   \code{\link{model_analysis}} calculates the accepted model.
#'
#'   \code{\link{variable_analysis}} calculates the individual variables
#'   significance.
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
#' X <- data.frame(X1 = rnorm(n, E, 1), X2 = rnorm(n, 3 * (E %in% c(1,5)), 1))
#' Y <- rnorm(n, X$X1, 1) # X1 is the true parent
#'
#' # Environment Relevance Test:
#' ICP(Y, X, E)
#'
#' # Intersecting Confidence Region Test, Quadratically Constrained Solver:
#' ICP(Y, X, E, method = "CR")
#'
#' # Intersecting Confidence Region Test, Pairwise Solver:
#' ICP(Y, X, E, method = "CR", solver = "pairwise")
#'
#' # Intersecting Confidence Region Test, Marginal Solver:
#' ICP(Y, X, E, method = "CR", solver = "marginal")
#'
#'
#' # ===========================================================================
#' # An example with a poisson distribution
#' Y <- rpois(n, exp(X$X1)) # true causal is X1
#'
#' # Environment Relevance Test
#' ICP(Y, X, E, model = "glm", family = "poisson")
#'
#' # Intersecting Confidence Region Test, Quadratically Constrained Solver:
#' ICP(Y, X, E, model = "glm", family = "poisson", method = "CR")
#'
#' # Intersecting Confidence Region Test, Pairwise Solver:
#' ICP(Y, X, E, model = "glm", family = "poisson",
#'     method = "CR", solver = "pairwise")
#'
#' # Intersecting Confidence Region Test, Marginal Solver:
#' ICP(Y, X, E, model = "glm", family = "poisson",
#'     method = "CR", solver = "marginal")
#'
#'
#' # ===========================================================================
#' # An example with right censored survival times
#' Y <- rexp(n, exp(- 0.5 * X$X1))
#' C <- rexp(n, exp(- 1.5))
#' time <- pmin(Y, C) # trues causal is X1
#' status <- time == Y
#'
#' # Environment Relevance Test
#' ICP(survival::Surv(time, status), X, E, model = "ph")
#'
#' # The user may also define their own link functions, see
#' # ?survival::survreg.distributions
#' my_dist <- survival::survreg.distributions$exponential
#' my_dist$trans <- function(y) log(y / 365)
#' my_dist$dtrans <- function(y) 1 / y
#' my_dist$itrans <- function(y) 365 * exp(y)
#' ICP(survival::Surv(time, status), X, E, model = "hazard", dist = my_dist)
#' # this example is simply a reparametrization and therefore
#' # gives the same results as above.
#'
#' # Intersecting Confidence Regions Test, Quadratically Constrained Solver:
#' ICP(survival::Surv(time, status), X, E, model = "ph", method = "CR")
#'
#' # Intersecting Confidence Regions Test, Pairwise Solver:
#' ICP(survival::Surv(time, status), X, E, model = "ph",
#'     method = "CR", solver = "pairwise")
#'
#' # Intersecting Confidence Regions Test, Marginal Solver:
#' ICP(survival::Surv(time, status), X, E, model = "ph",
#'     method = "CR", solver = "marginal")
#'
#' # Non-parametric Tests of Time Varying Effect
#' ICP(survival::Surv(time, status), X, E, model = "ph", method = "TimeVar")
#'
#' # Non-parametric Tests of Time Varying Effect with n.sim = 1000
#' ICP(survival::Surv(time, status), X, E, model = "ph",
#'     method = "TimeVar", n.sim = 1000)
#'
#' @export
ICP <- function(Y, X, E = NULL, model = "lm", method = "EnvirRel", level = 0.05,
                gof = max(0.01, level), maxNoVariables = 8, fullAnalysis = FALSE,
                progress = FALSE, ...) {

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
    if (is.null(gof)) {
      stop("When 'level' is specified then 'gof' must also be specified.")
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
  no_test_list <- NULL


  # LOOP THROUGH TEST_LIST -----------------------------------------------------

  # initialize
  smallest_possible_pvalue <- .Machine$double.eps
  res_model <- data.frame(models = character(), pval = numeric())
  s <- 1L
  if (progress) {
    pb <- utils::txtProgressBar(min = 0, max = 100,
                                width = 0.8 * getOption("width"),
                                style = 3)
  }

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
      if (! is.null(level)) {
        if (pval >= level) {
          if (nrow(res_model) == 1) { # empty set accepted
            test_list <- test_list[[1]]
            no_test_list <- test_list[-1]
          } else {
            remove <- sapply(seq_along(test_list), function(i) {
              if (i <= s) {FALSE} else {all(test_list[[s]] %in% test_list[[i]])}
            })
            no_test_list <- c(no_test_list, test_list[remove])
            test_list <- test_list[!remove]
          }
        }
      }
    }

    if (progress) {
      utils::setTxtProgressBar(pb, value = (s / length(test_list)) * 100)
    }
    s <- s + 1L
  }

  # SUMMARIZE RESULTS ----------------------------------------------------------
  call <- match.call()
  output <- structure(
    list(
      call = call,
      level = level,
      gof = gof,
      smallest_possible_pvalue = smallest_possible_pvalue,
      covariate.list.codes = data.frame("variable" = c("Empty", names(X)),
                                        "code" = c(0, seq_along(names(X)))),
      model.analysis = res_model,
      tested.list = test_list,
      not.tested.list = no_test_list,
      method = method
      ),
    class = "ICP"
  )

  if (!is.null(level)) {
    tmp <- model_analysis(output, level = level, gof = gof)
    output$accepted.model <- tmp$accepted.model
    output$empty.message <- tmp$empty.message
    if (fullAnalysis) {
      tmp <- variable_analysis(output, gof = gof)
      output$variable.analysis <- tmp$variable.analysis
    }
  }
  if (progress) {
    close(pb)
  }
  return(output)
}


