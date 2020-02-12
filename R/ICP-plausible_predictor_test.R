#' Internal Test for plausible predictors
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
#'
#' @param method a method object describing the model class to use in the
#'   fitting procedure and the testing method. When used in the
#'   \code{\link{ICP}} function the method object is inherited from the
#'   \code{\link{ICP}} function.
#' @param Y a vector or \code{\link[survival]{Surv}} object describing the
#'   target variable. The \code{Y} will be passed to the \code{\link{fit_model}}
#'   function for fitting, and it is therefor important that the class of
#'   \code{Y} is understood by the regression method used in
#'   \code{\link{fit_model}}. So if \code{method} specifies that a cox
#'   regression is to be used, then \code{Y} must be a
#'   \code{\link[survival]{Surv}} object.
#' @param X a matrix, vector or data frame describing the covariates.
#' @param E a vector describing the environmants
#' @param level the alpha level of testing, but it is only relevant if the
#'   method is iterative \emph{and} \code{fullAnalysis} is \code{FALSE}.
#' @param fullAnalysis if \code{TRUE} a p-value must be retured. If \code{FALSE}
#'   the method is 'allowed' to only return 0/1 for hypothesis rejected/accepted.
#' @param Bonferroni if \code{TRUE} Bonferroni correcting will be used if relevant.
#' @param ... additional arguments to be passed to lower level functions.
#'
#'
#' @return Returns a p-value, i.e. a number between 0 and 1, to be used in the
#'   \code{\link{ICP}} function.
#'
#' @seealso \code{\link{ICP}} for the full wrapper function.
#'
#' @keywords internal
#' @export

plausible_predictor_test <- function(method, Y, X, ...) {
  UseMethod("plausible_predictor_test", method)
}

#' @rdname plausible_predictor_test
#' @export
plausible_predictor_test.default <- function(method, Y, X, ...) {
  stop(method$method,
          "is not recocnised as a testing method by the",
          "'plausible_predictor_test'")
  return(0)
}


#' @rdname plausible_predictor_test
#' @export
plausible_predictor_test.TimeVar <- function(method, Y, X,
                                              level, Bonferroni = TRUE, ...) {
  fit <- fit_nonparam_model(method, Y, X, ...)
  if (method$n.sim != 0) {
    assign("smallest_possible_pvalue", 1 / method$n.sim, envir = parent.frame())
    if(method$nonparamtest == "sup") {
      pvals <- fit$sup
    } else if (method$nonparamtest == "int") {
      pvals <- fit$int
    } else {
      pvals <- c(fit$sup, fit$int)
    }
  } else {
    # TODO remember to use bonferroni correction
    stop("Khmaladze transformation not yet implemented")
  }
  if (sum(is.na(pvals)) > 0) {
    stop("Faild to compute the necessary p-vlaues:\n",
         "The error is in the 'timereg' package and has been reported.")
  }
  BonCorr <- ifelse(Bonferroni, length(pvals[!is.na(pvals)]), 1)
  res <- min(1, min(pvals) * BonCorr, na.rm = TRUE) # min(pvals) or max(pvals)
  return(res)
}



#' @rdname plausible_predictor_test
#' @export
plausible_predictor_test.EnvirRel <- function(method, Y, X, E, ...) {
  if (is.null(E)) {
    stop("The 'EnvirRel' method needs environments")
  }
  fitE <- fit_model(method, Y, cbind(X,E))
  fit0 <- fit_model(method, Y, X)
  pval <- stats::pchisq((fit0$deviance - fitE$deviance)/fitE$scale,
                        df = fit0$df - fitE$df,
                        lower.tail = FALSE)
  return(pval)
}






check_for_degeneracy <- function(Y, X, E, ...) {
  UseMethod("check_for_degeneracy", Y)
}
check_for_degeneracy.default <- function(Y, X, E, ...) {
  text <- paste0("Target variable of class ", class(Y), " is not recognised by ",
                 "the 'check_for_degeneracy' test for the Intersecting ",
                 "Confidence Regions test (CR).")
  stop(text)
}
check_for_degeneracy.numeric <- function(Y, X, E, ...) {
  index <- lapply(unique(E), function(e) {
    which(E == e)
  })
  Y_check <- X_check <- TRUE
  for (i in seq_along(index)) {
    Y_check <- ifelse(length(unique(Y[index[[i]]])) > 1, TRUE, FALSE)
    if (!Y_check) {break}
    for (j in seq_len(ncol(X))) {
      X_check <- ifelse(length(unique(X[index[[i]],j])) > 1, TRUE, FALSE)
      if (!X_check) {break}
    }
  }
  if (Y_check & X_check) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}
check_for_degeneracy.Surv <- function(Y, X, E, ...) {
  index <- lapply(unique(E), function(e) {
    which(E == e)
  })
  X_check <- TRUE
  for (i in seq_along(index)) {
    for (j in seq_len(ncol(X))) {
      X_check <- ifelse(length(unique(X[index[[i]],j])) > 1, TRUE, FALSE)
      if (!X_check) {break}
    }
  }
  if (X_check) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}




#' @rdname  plausible_predictor_test
#' @export
plausible_predictor_test.CR <- function(method, Y, X, E,
                                        level, fullAnalysis,
                                        Bonferroni = TRUE, ...) {
  if (!check_for_degeneracy(Y, X, E)) {
    stop("One or more variables become degenerated when subsetting data by ",
         "environment 'E'")
  }
  if (is.null(method$splits)) {
    method$splits <- "all"
  }
  if (method$splits == "LOO") {
    models <- lapply(unique(E), function(e) {
      list(fit_model(method, Y, X, subset = (E == e)),
           fit_model(method, Y, X, subset = (E != e)))
    })
    BonCorr <- ifelse(Bonferroni, length(models) * 2, 1)
    alpha <- ifelse(is.null(level), NULL, level / BonCorr)
    res <- sapply(models, function(m) {
      intersect_CR(method, m, alpha, fullAnalysis, ...)
    })
    res <- ifelse(is.logical(res),
                  all(res),
                  min(1, min(res) * BonCorr))
  } else {
    models <- lapply(unique(E), function(e) {
      fit_model(method, Y, X, subset = (E == e))
    })
    BonCorr <- ifelse(Bonferroni, length(models), 1)
    alpha <- ifelse(is.null(level), NULL, level / BonCorr)
    res <- intersect_CR(method, models, alpha, fullAnalysis, ...)
    res <- ifelse(is.logical(res),
                  res,
                  min(1, res * BonCorr))
  }
  return(res)
}

intersect_CR <- function(method, models, alpha, ...) { # TODO
  UseMethod("intersect_CR", method)
}

intersect_CR.marginal <- function(method, models, alpha,
                                  fullAnalysis = FALSE, tol = 2e-10, ...) {
  assign("smallest_possible_pvalue", tol, envir = parent.frame(n = 2))
  if (is.null(alpha) | fullAnalysis) {
    max_true <- 0
    min_false <- 1
    while ((min_false - max_true) > tol) {
      alpha <- (min_false + max_true) / 2
      test <- intersect_marginal_rectangle(models, alpha)
      if (test) {max_true <- alpha} else {min_false <- alpha}
    }
    return(max_true) # alpha value
  } else {
    return(intersect_marginal_rectangle(models, alpha)) # boolean
  }
}

intersect_CR.pairwise <- function(method, models, alpha,
                                  fullAnalysis = FALSE, tol = 2e-10, ...) {
  assign("smallest_possible_pvalue", tol, envir = parent.frame(n = 2))
  pairs <- utils::combn(seq_along(models), 2)
  if (is.null(alpha) | fullAnalysis) {
    max_true <- 0
    min_false <- 1
    while ((min_false - max_true) > tol) {
      alpha <- (min_false + max_true) / 2
      test <- apply(pairs, 2, function(c) {
        intersect_ellipsoid_pair(models[[c[1]]], models[[c[2]]], alpha)
      })
      if (all(test)) {max_true <- alpha} else {min_false <- alpha}
    }
    return(max_true) # alpha value
  } else {
    res <- apply(pairs, 2, function(c) {
      intersect_ellipsoid_pair(models[[c[1]]], models[[c[2]]], alpha)
    })
    return(all(res)) # boolean
  }
}

intersect_CR.QC <- function(method, models, alpha, fullAnalysis = TRUE, ...) {
  m <- length(models[[1]]$coefficients)
  if (fullAnalysis) {
    res <- intersect_ellipsoid(models)
    return(stats::pchisq(res[1], df = m, lower.tail = FALSE)) # alpha value
  } else {
    test <- intersect_marginal_rectangle(models, alpha)
    if (test) {
      res <- intersect_ellipsoid(models)
      return(stats::pchisq(res[1], df = m, lower.tail = FALSE)) # alpha value
    } else {
      return(FALSE)
    }
  }
}


intersect_marginal_rectangle <- function(models, alpha) {
  all_names <- Reduce(union,
                      lapply(models, function(m) {names(m$coefficients)}))
  q <- stats::qchisq(1 - alpha, df = length(all_names))
  endpoints <- sapply(models, function(m) {
    L <- m$coefficients - sqrt(q) * sqrt(diag(m$covariance))
    U <- m$coefficients + sqrt(q) * sqrt(diag(m$covariance))
    #if (! identical(names(m$coefficients), all_names)) {
    #  mis <- which(! all_names %in% names(m$coefficients))
    #  for (i in seq_along(mis)) {
    #    L <- append(L, -Inf, after = mis[i] - 1)
    #    U <- append(U, Inf, after = mis[i] - 1)
    #  }
    #}
    return(rbind(L,U))
  })
  endpoint_test <- sapply(seq_len(nrow(endpoints) / 2), function(u) {
    max(endpoints[u * 2 - 1, ], na.rm = TRUE) <= min(endpoints[u * 2, ], na.rm = TRUE)
  })
  return(all(endpoint_test))
}


intersect_ellipsoid_pair <- function(M1, M2, alpha){

  c1 <- M1$coefficients
  C1 <- M1$covariance
  c2 <- M2$coefficients
  C2 <- M2$covariance

  if (any(!is.finite(C1))) {
    ind <- which(!is.finite(C1))
    C1[ind] <- max(C1[-ind], 1) * 100000
  }
  if (any(!is.finite(C2))) {
    ind <- which(!is.finite(C2))
    C2[ind] <- max(C2[-ind], 1) * 100000
  }

  m <- length(c1)
  q <- sqrt(stats::qchisq(1 - alpha, df = m))

  if (m == 1) {
    return(abs(c1 - c2) < (sqrt(C1) + sqrt(C2)) * q)
  } else {

    e1 <- eigen(C1)
    d1 <- e1$values
    d1 <- make_positive(d1)
    d1 <- sqrt(d1)
    U1 <- e1$vectors

    c <- (1 / q) * diag(1 / d1) %*% t(U1) %*% (c2 - c1)
    C <- diag(1 / d1) %*% t(U1) %*% C2 %*% U1 %*% diag(1 / d1)

    e <- eigen(C)
    U <- e$vectors
    d <- e$values
    d <- make_positive(d, silent = TRUE)
    d <- sqrt(d)

    y <- -t(U) %*% c
    y <- abs(y) # 0 expressed in coordinate system of l and rotated to the first quadrant

    if(sum((y / d) ^ 2) <= 1){ # y inside the ellipse
      return(TRUE)
    } else { # Newton-Rhapson iterations

      # f goes monotone, quadratically to -1, so sure and fast convergence
      f <- function(t) sum(( d * y / (t + d^2))^2) - 1
      df <- function(t) - 2 * sum((y * d)^2 / (t + d^2)^3)

      t0 <- 0
      ft0 <- f(t0)

      while(ft0 > 1e-4){
        t0 <- t0 - f(t0) / df(t0)
        ft0 <- f(t0)
      }

      x0 <- y * d^2 / (d^2 + t0) # projection of y onto (c,C)
      dist <- sqrt(sum((y - x0)^2))
      return(dist < 1)
    }
  }
}

# (very) adhoc way of dealing with negative eigenvalues
make_positive <- function (v, silent = TRUE){
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



intersect_ellipsoid <- function(models) { # should return number
  # -- tuning parameters -------------------------------------------------------
  barrier_iter <- 5 # find better stopping criteria ?
  tol_newton <- 0.1
  small <- 1
  mu <- 20
  linesearch <- TRUE # page 464 Boyd & Vandenbergh
  alpha <- 0.3 # must be between 0 and 0.5
  beta <- 0.8 # must be between 0 and 1

  # -- Find relevant sizes -----------------------------------------------------
  ellipsoids <- lapply(models, function(m) {
    var <- m$covariance
    if (any(!is.finite(var))) {
      ind <- which(!is.finite(var))
      var[ind] <- max(var[-ind], 1) * 1000
    }
    rr <- rcond(var)
    if (rr < 2e-15) {
      add <- diag(diag(var) * 0.01, ncol(var))
      P <- solve(var + add)
    } else {
      P <- solve(var)
    }
    b <- matrix(m$coefficients, ncol = 1)
    q <- - 2 * (P %*% b)
    return(list(P = P, b = b, q = q))
  })

  # -- create relevant functions -----------------------------------------------
  p <- length(ellipsoids)
  f <- function(x) {
    sapply(ellipsoids, function(e) {
      diff <- x - e$b
      t(diff) %*% e$P %*% diff
    })
  }
  gf <- function(x) {
    lapply(ellipsoids, function(e) {
      2 * e$P %*% x + e$q
    })
  }
  hf <- function(x) {
    lapply(ellipsoids, function(e) {
      2 * e$P
    })
  }
  g <- function(s, tune, f) {tune * s - sum(log(s - f)) - log(s)}
  linesearch_test <- function(s, x, f, decrement, tune, delta, step) {
    s_delta <- s + step * delta[1, ]
    if (s_delta <= 0) {return(TRUE)}
    f_delta <- f(x + step * delta[-1, ])
    if (max(f_delta) >= s_delta){return(TRUE)}
    if(g(s_delta, tune, f_delta) >
       g(s, tune, f) + alpha * step * decrement) {
      return(TRUE)
    }
    return(FALSE)
  }
  gradient<- function(s, tune, f, gf) {
    top <- tune - sum(1 / (s - f))
    v <- lapply(1:p, function(i) {
      (1 / (s - f[[i]])) * gf[[i]]
    })
    bottom <- Reduce("+", v)
    return(matrix(c(top, bottom), ncol = 1))
  }
  hessian <- function(s, f, gf, hf) {
    dsds <- sum(1 / (s - f) ^ 2)
    v <- lapply(1:p, function(i) {
      (1 / (s - f[[i]]) ^ 2) * gf[[i]]
    })
    dsdx <- - Reduce("+", v)
    m <- lapply(1:p, function(i) {
      (1 / (s - f[[i]]) ^ 2) * gf[[i]] %*% t(gf[[i]]) +
        (1 / (s - f[[i]])) * hf[[i]]
    })
    dxdx <- Reduce("+", m)
    return(cbind(rbind(dsds, dsdx), rbind(t(dsdx), dxdx)))
  }


  # -- logarithmic barrier method ----------------------------------------------
  #x <- matrix(rowMeans(sapply(ellipsoids, function(e){e$b})), ncol = 1)
  x <- ellipsoids[[1]]$b
  s <- max(f(x)) + small
  tune <- (sum(1 / (s - f(x))) * s + 1) / s

  for (i in seq_len(barrier_iter)) {

    # -- centering :  Newtons method
    lambda <- tol_newton * 5
    while (lambda / 2 > tol_newton) {
      # .. calculate direction and decrement
      f_val <- f(x)
      gf_val <- gf(x)
      hf_val <- hf(x)
      gx <- gradient(s, tune, f_val, gf_val)
      #invhx <- solve(hessian(s, f_val, gf_val, hf_val))
      hg <- hessian(s, f_val, gf_val, hf_val)
      rr <- rcond(hg)
      if (rr < 1e-10) {
        if (rr < 1e-16) {
          break
        } else{
          add <- diag(diag(hg) * 0.01, ncol(hg))
          invhx <- solve(hg + add)
        }
      } else {
        invhx <- solve(hg)
      }
      delta <- - invhx %*% gx
      lambda <- - t(gx) %*% delta

      # .. linesearch
      step <- 1

      while (linesearch_test(s, x, f_val, lambda, tune, delta, step)) {
        step <- beta * step
      }

      # .. update (s,x)
      s <- s + step * delta[1, ]
      x <- x + step * delta[-1, ]
    }
    tune <- tune * mu
  }
  return(c(s,x))
}
