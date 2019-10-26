
#' Simulate data from graph
#'
#' Simulate data from a graph.
#'
#' This function is designed a bit ad hoc to make simulations a bit easier. It
#' might not survive to a publication of the code.
#'
#'
#' @param G A weighted adjecency matrix of a graph. This must be an upper
#'   triangular matrix, whith zeros in the diagonal. The value of \code{G[i,j]}
#'   is multiplied with variable \code{i} to form variable \code{j}. Hence a
#'   non-zero values in \code{G[i,j]} indicates an arrow from variable \code{i}
#'   to variable \code{j} in the graph and the value indicates weight.
#'
#' @param n number of individuals.
#'
#' @param m number of observations for each individual. If not specified taken
#'   to be 1.
#'
#' @param sim a list or vector of strings of length shorter then or equal to the
#'   number of columns in \code{G}. If \code{sim} has length shorter then number
#'   of \code{G} columns then the last entry of \code{sim} will be repeated
#'   until \code{sim} has the same length as \code{G} has columns.
#'
#'   Each entry of \code{sim} may contain a \code{BETA} and an \code{N}, where
#'   \code{BETA} is the linear combination of the parent variables found via
#'   \code{G} and \code{N} is length the \code{sim} function should be.
#'
#' @param OutGroup If \code{NULL} the function will retun a single data frame
#'   with all the simulated data. Otherwise \code{OutGroup} must be a list,
#'   where each list entry is vector of indicies. The indicies must correspond
#'   to either the names og column number in the graph matrix \code{G}. If this
#'   is the case, then the output of the function is a list of data frames of
#'   the simulated data split based on \code{OutGroup}.
#'
#' @return A data frame  or list of data frames with a column for each column in
#'   the adjecency matrix \code{G}. If \code{m} is larger then 1 then there will
#'   also be columns for time and id.
#'
#' @examples
#' ### A simple example ========================================================
#'
#' # make a graph
#' G <- matrix(c(0,0,3,0), ncol = 2, dimnames = list(NULL, c("X", "Y")))
#'
#' # specify distributions
#' sim <- c("rbinom(N, 3, 0.4)", "rnorm(N, BETA, 1)")
#'
#' # simualte
#' SIM <- sim_from_adj(G, n = 100, sim)
#'
#' # ---------------------------------------------------------------------------
#'
#'
#' ### A more advanced example =================================================
#'
#' # make graph
#' G <- matrix(c(0,0,1,3,0,0,
#'               0,0,0,0,1,0,
#'               0,0,0,0,0,1,
#'               0,0,0,0,0,3,
#'               0,0,0,0,0,0,
#'               0,0,0,0,0,0), ncol = 6, byrow = T,
#'               dimnames = list(NULL, c("E1","E2","X1","X2","X3","Y")))
#'
#' # specify grouping in output
#' out <- list(E = 1:2, X = 3:5, Y = 6)
#'
#' # specify distributions
#' sim <- c("rbinom(N, 3, 0.3)", "rbinom(N,1, 0.2)",
#'          "rnorm(N, BETA, 1)", "rpois(N, exp(BETA))", "rnorm(N, BETA, 0.3)",
#'          "rpois(N, exp(BETA))")
#'
#' # simulate
#' SIM <- sim_from_adj(G, n = 100, sim, out)
#'
#' # ---------------------------------------------------------------------------
#'
#' @export

sim_from_adj <- function(G, n, m = 1, sim, OutGroup = NULL) {
  if (missing(G) | missing(n)) {
    stop("Both 'G', 'n' and 'sim' must be supplied")
  }
  if (missing(sim)) {
    if (!missing(m) & !is.numeric(m)) {
      sim <- m
      m <- 1
    } else {
      stop("Both 'G', 'n' and 'sim' must be supplied")
    }
  } else {
    if (missing(OutGroup) & is.list(sim) & !is.numeric(m)) {
      OutGroup <- sim
      sim <- m
      m <- 1
    }
  }
  if (!is.numeric(n) | !is.numeric(m)) {
    stop("'n' and 'm' must be numbers larger then zero")
  }
  if (!(n > 0) |  !(m > 0)) {
    stop("'n' and 'm' must be numbers larger then zero")
  }
  if (ncol(G) != nrow(G)) {
    stop("'G' must be a square matrix")
  }
  if (sum(G[lower.tri(G, diag = T)]) != 0) {
    stop("'G' must be an upper triangular matrix")
  }
  if (length(sim) > ncol(G)) {
    stop("The length of the 'sim' list must not be larger",
         " then the number of columns in G")
  } else if (length(sim) == 0) {
    stop("there must be at least one entry in 'sim'")
  } else {
    dist <- sim[c(seq_along(sim),rep(length(sim), ncol(G) - length(sim)))]
  }

  if (m == 1) {
    N <- n
    data <- matrix(0, ncol = ncol(G), nrow = n, dimnames = list(NULL, colnames(G)))
    for (j in seq_len(ncol(G))) {
      BETA <- data %*% G[ ,j]
      tmp <- eval(parse(text = dist[[j]]))
      data[ ,j] <- tmp
    }
    DATA <- data.frame(data)
  } else if (m > 1) {
    N <- m
    DATA <- lapply(seq_len(n), function(i) {
      data <- matrix(0, ncol = ncol(G), nrow = m, dimnames = list(NULL, colnames(G)))
      for (j in seq_len(ncol(G))) {
        BETA <- data %*% G[ ,j]
        tmp <- eval(parse(text = dist[[j]]))
        data[ ,j] <- tmp
      }
      return(data.frame(data, "id" = i, "time" = 1:m))
    })
    DATA <- Reduce(rbind, DATA)
    DATA <- data.frame(DATA)
  }


  if (!is.null(OutGroup)) {
    tmp <- DATA
    DATA <- lapply(OutGroup, function(u) {
      tmp[ , u]
    })
    if (all(c("id","time") %in% colnames(tmp))) {
      DATA$ID <- tmp[ , c("id", "time")]
    }
  }

  class(DATA) <- c("ICP_sim", class(DATA))
  return(DATA)
}


plot.ICP_sim <- function(x, y = NULL) {
  if (ncol(x) == 3) {
    if (colnames(x)[1] == "time") {
      plot(0, type = "n", xlim = c(min(x[ ,1]),max(x[ ,1])),
           ylim = c(min(x[ ,2]), max(x[ ,2])),
           xlab = colnames(x)[1], ylab = colnames(x)[2])
      for (i in unique(x[ ,3])) {
        lines(x[x[ ,3] == i, -3], col = i)
      }
    } else {
      plot.default(x[ ,1:2], col = x[ ,3])
    }
  } else {
    plot.default(x)
  }
}



