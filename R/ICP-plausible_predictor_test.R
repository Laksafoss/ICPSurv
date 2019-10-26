#' Test for plausible predictors
#'
#' \code{plausible_predictor_test} is a generic function used to test whether a
#' set of predictors might be plausible causal predictors.
#'
#' The \code{plausible_predictor_test} is meant to be used in the wrapper
#' function \code{\link{ICP}}, but can also be used on its own as a simple
#' hypothesis testing function. The \strong{method object} \code{method}, which
#' is created by the \code{\link{method_obj}} function dictates the
#' \code{plausible_predictor_test} method. As this function is generic which
#' means that new testing methods can easily be added for new classes.
#'
#' For the \code{plausible_predictor_test} to work correctly in the
#' \code{\link{ICP}} wrapper it must be a test of the null hypothesis
#' \deqn{H_{0,S}: S is an invariant set w.r.t. (X,Y)}
#' where an invariant set is defined as a set of indicies \eqn{S} such that
#' \deqn{(Y^e | X^e_S) = (Y^f | X^f_S)}
#' in distribution for all environments \eqn{e,f}. If the data is time dependent
#' we formulate an analog invariance statement for all time points \eqn{t} and
#' \eqn{s}. For more detailes and discussion of this hypothesis and invariance
#' see the references.
#'
#' As the \code{plausible_predictor_test} is a hypothesis test it must return a
#' p value, that is a number between 0 and 1.
#'
#' The \code{plausible_predictor_test}
#'
#' PERHAPS A DESCRIPTION OF THE TACHNICALITIES OF \code{method} AND \code{Y} ?!
#'
#' @param method a method object describing the model class to use in the
#'   fitting procedure and the testing method. When used in the
#'   \code{\link{ICP}} function the method object is inherited from the
#'   \code{\link{ICP}} function.
#' @param Y a vector or \code{\link[survival]{Surv}} object describing the target variable. The \code{Y}
#'   will be passed to the \code{\link{model_fit}} function for fitting, and it
#'   is therefor important that the class of \code{Y} is understood by the
#'   regression method used in \code{\link{model_fit}}. So if \code{method}
#'   specifies that a cox regression is to be used, then \code{Y} must be a
#'   \code{\link[survival]{Surv}} object.
#' @param X a matrix, vector or data frame describing the covariates. WRITE
#'   ABOUT \code{X} WHEN TESTING THE EMPTY SET !!!!
#' @param E a vector describing the environmants
#' @param level the alpha level of testing, but it is only relevant if the
#'   method is iterative \emph{and} \code{fullAnalysis} is \code{FALSE}.
#' @param fullAnalysis only has an effect if the method used is iterative. If
#'   \code{TRUE} the method will go throught the iterative steps to determin a
#'   p-vlaue, however if \code{FALSE} the method will only be tested at the
#'   specified alpha level \code{level}.
#'
#'
#' @return Returns a p-value, i.e. a number between 0 and 1, to be used in the
#'   \code{\link{ICP}} function.
#'
#' @seealso \code{\link{ICP}} for the full wrapper function.
#'
#' @examples
#' # create some data
#' n <- 100
#' E <- rbinom(n, 3, 0.3)
#' X <- rnorm(n, E, 1)
#' Y <- rpois(n, exp(X))
#' E <- as.factor(E)
#'
#' # create the method object
#' method <- structure(list(model = "glm", family = "poisson"),
#'                     class = "EnvirIrrel")
#'
#' # test hypothesis using 'EnvirIrrel' method
#' plausible_predictor_test(method, Y, X, E)
#'
#'
#' @export

plausible_predictor_test <- function(method, Y, X, ...) {
  UseMethod("plausible_predictor_test", method)
}


#' @export

plausible_predictor_test.default <- function(method, Y, X, ...) {
  warning(method$method,
          "is not recocnised as a testing method by the",
          "'plausible_predictor_test'")
  return(0)
}






#' @rdname plausible_predictor_test
#' @export


## p-val of t-test for E = 2 * pnorm(-abs(tval)) with tval = coef[E] / sqrt(covar[E])
plausible_predictor_test.EnvirIrrel <- function(method, Y, X, E, ...) {
  if (is.null(E)) {
    stop("The 'EnvirIrrel' method needs environments")
  }
  if (length(X) == 0) {
    fitE <- fit_model(method, Y, E)
    pval <- 1- pchisq(fitE$null.deviance - fitE$deviance, df = fitE$df)
  } else {
    fitE <- fit_model(method, Y, cbind(X,E))
    fit0 <- fit_model(method, Y, X)
    pval <- 1 - pchisq(fit0$deviance - fitE$deviance, df = fitE$df - fit0$df)
  }
  return(pval)
}




#' @rdname plausible_predictor_test

plausible_predictor_test.CR <- function(method, Y, X, E,
                                        level, fullAnalysis, ...) {

  specialnull <- # TODO this is for coxph and aalen semiparametric

  if (length(X) == 0 & specialnull) {
    return(dHSIC::dhsic.test(as.numeric(E), Y, method = "gamma")$p.value)
  }

  if (method$split == "LOO") {
    models <- lapply(unique(E), function(e) {
      list("E" = fit_model(method, Y[E == e], X[E == e, ]),
           "O" = fit_model(method, Y[E != e], X[E != e, ]))
    })
    # TODO define BonCorr here.....
  } else if (mthod$split == "pairwise") {
    models <- lapply(unique(E), function(e) {
      fit_model(method, Y[E == e], X[E == e, ])
    })
    BonCorr <- length(models[[1]]$coefficients)
  }

  # TODO :
  # Should we simple write intersect_CR(method, models, level, fullAnalysys)
  # or should we have some if statements for fullAnalysis = F, alpha grid and so on..

}

intersect_CR <- function(method, models, q) { # TODO
  UseMethod("intersect_CR", method)
}

intersect_CR.marginal <- function(method, models, q) {
  # THE RECTANGLE CASE
}

intersect_CR.pairwise <- function(method, models, q) {
  # RUNES
}

intersect_CR.QCLP <- function(method, models, q) {
  # Write the QCLP as a SOCP
}




#' @rdname plausible_predictor_test

plausible_predictor_test.CRrectangle <- function(method, Y, X, E,
                                                 level, fullAnalysis, ...) {
  if (length(X) == 0 & length(grep("glm", method$model)) == 0){
    pval <- dHSIC::dhsic.test(as.numeric(E), Y, method = "gamma")$p.value
  } else {
    models <- lapply(unique(E), function(e) {
      fit_model(method, Y[E == e], X[E == e, ])
    })
    BonCorr <- length(models[[1]]$coefficients)

    if (!fullAnalysis & !is.null(level)) {
      q <- qnorm(1-level/(BonCorr*2))
      test <- intersect_CI(models, q)
      pval <- if (test) { 1 } else { 0 }
    } else {
      max_true <- 0
      min_false <- 1
      while ((min_false - max_true) >  2e-10){
        alpha <- (2 * min_false + max_true) / 3  # TODO is this how we want to do it ??
        q <- qnorm(1-alpha/(BonCorr*2))
        test <- intersect_CI(models, q)
        if (test) { max_true <- alpha } else { min_false <- alpha }
      }
      pval <- max_true
    }
  }
  return(pval)
}


intersect_CI <- function(models, q) {
  endpoints <- sapply(models, function(m) {
    L <- m$coefficients - q * sqrt(diag(m$covariance))
    U <- m$coefficients + q * sqrt(diag(m$covariance))
    return(rbind(L,U))
  })
  endpoint_test <- sapply(seq_len(nrow(endpoints)/2), function(u) {
    max(endpoints[u * 2 - 1, ]) <= min(endpoints[u * 2, ])
  })
  return(all(endpoint_test))
}



#' @rdname plausible_predictor_test
#' @export
plausible_predictor_test.CRellipsoid <- function(method, Y, X, E,
                                                 level, fullAnalysis, ...) {
  if (length(X) == 0 & length(grep("glm", method$model)) == 0){
    pval <- dHSIC::dhsic.test(as.numeric(E), Y, method = "gamma")$p.value
  } else {
    models <- lapply(unique(E), function(e) {
      fit_model(method, Y[E == e], X[E == e, ])
    })
    BonCorr <- length(models[[1]]$coefficients)

    if (!fullAnalysis & !is.null(level)) {
      alpha_BonCorr <- level / BonCorr
      q <- qnorm(1 - level / (BonCorr * 2))
      test <- intersect_ellipsoid(models, q)
      pval <- if (test) { 1 } else { 0 }
    } else {
      max_true <- 0
      min_false <- 1
      while ((min_false - max_true) >  2e-10){
        alpha <- (2 * min_false + max_true) / 3  # TODO is this how we want to do it ??
        alpha_BonCorr <- alpha / BonCorr
        test <- intersect_ellipsoid(models, alpha_BonCorr)
        if (test) { max_true <- alpha } else { min_false <- alpha }
      }
      pval <- max_true
    }
  }
  return(pval)
}



intersect_ellipsoid <- function(models, alpha) {
  pairs <- combn(length(models),2)
  tests <- sapply(seq_len(ncol(pairs)), function(k) {
    intersect_ellipsoid_pair(models[[pairs[1,k]]], models[[pairs[2,k]]], alpha)
  })
  return(all(tests))
}


# FROM RUNE :     ===========================================================  #
# checks if two ellipsoid confidence regions with centers ci,
# covariance matrices Ci and coverage 1-alpha intersect
intersect_ellipsoid_pair <- function(M1, M2, alpha){

  c1 <- M1$coefficients
  C1 <- M1$covariance
  c2 <- M2$coefficients
  C2 <- M2$covariance

  if (alpha == 1) return(FALSE)
  if (alpha == 0) return(TRUE)

  m <- length(c1)
  q <- sqrt(qchisq(1-alpha, df = m))

  if (m == 1) {
    return(abs(c1-c2) < (sqrt(C1)+sqrt(C2))*q)
  } else {

    e1 <- eigen(C1)
    d1 <- e1$values
    d1 <- makePositive(d1)
    d1 <- sqrt(d1)
    U1 <- e1$vectors

    c <- 1/q * diag(1/d1) %*% t(U1) %*%(c2-c1)
    C <- diag(1/d1) %*% t(U1) %*% C2 %*% U1 %*% diag(1/d1)

    e <- eigen(C)
    U <- e$vectors
    d <- e$values
    d <- makePositive(d,silent=TRUE)
    d <- sqrt(d)

    y <- -t(U)%*%c
    y <- abs(y) # 0 expressed in coordinate system of l and rotated to the first quadrant

    if(sum((y/d)^2) <= 1){ # y inside the ellipse
      return(TRUE)
    } else { # Newton-Rhapson iterations

      # f goes monotone, quadratically to -1, so sure and fast convergence
      f <- function(t) sum((d*y/(t+d^2))^2)-1
      df <- function(t) -2*sum((y*d)^2/(t+d^2)^3)

      t0 <- 0
      ft0 <- f(t0)

      while(ft0>1e-4){
        t0 <- t0 - f(t0)/df(t0)
        ft0 <- f(t0)
      }

      x0 <- y*d^2/(d^2+t0) # projection of y onto (c,C)
      dist <- sqrt(sum((y-x0)^2))
      return(dist < 1)
    }
  }
}

# (very) adhoc way of dealing with negative eigenvalues
makePositive <- function (v, silent = TRUE){
  w <- which(v < 10^(-14))
  if (length(w) > 0 & !silent)
    warning("Some eigenvalues are below 10^(-14) and will therefore be regularized")
  for (ww in w) {
    if (v[ww] < 0) {
      v[ww] <- v[ww] - 2 * v[ww] + 10^(-14)
    }
    else {
      v[ww] <- v[ww] + 10^(-14)
    }
  }
  return(v)
}
# end FROM RUNE ============================================================== #


#' @rdname plausible_predictor_test
#' @export
plausible_predictor_test.ConstRegParam <- function(method, Y, X, level, ...) {
  fit <- fit_model(method, Y, X, const = F)
  pvals <- fit$pvals
  max_pval <- max(pvals)
  #test <- all(pvals >= level/length(pvals))
  #pval <- if (test) { 1 } else { 0 }
  return(max_pval/length(pvals))
}
