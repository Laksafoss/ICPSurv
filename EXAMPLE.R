


library(ICPSurv)

#  ===  E X A M P L E    1  ====================================================
# DISTRIBUTION OF Y :  Gaussian
# TIME DEPENDENT    :  NO
# TRUE PARENTS      :  X1 + X2
# NUMBER OF OBS.    :  50

# ENVIRONMENT DETAILES:
#   - number of environments :  3
#   - effect size on mean    :  0, 5, 10

# RESULT OF ICP ANALYSIS:
#   - EnvirIrrel  :  X1 + X2 (correct)
#   - CRrectangle :  X1 + X2 (correct)
#   - CRellipsoid :  X1 + X2 (correct)

set.seed(123)
n <- 500  # number of data points

# First we simulate data.
G <- matrix(c(0,1,1,1,0,0,0,0,0,2,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0),
            ncol = 5, byrow = T,
            dimnames = list(c("E","X1","X2","X3","Y"), c("E","X1","X2","X3","Y")))
out <- list(X = 2:4, Y = 5, E = 1, full = 1:5)
sim <- c("sample(c(0,5,10),N,replace=T)", # specify environment
         rep("rnorm(N,BETA,1)",4))        # specify X's and Y
SIM <- sim_from_adj(G, n, sim, out)

# If one wants to check if Y is simulated correctly
cbind("TRUE" = G[c("X1", "X2"),"Y"],
      stats::confint(stats::glm(Y ~ X1 + X2 - 1, family = "gaussian", data = SIM$full)))

# analyze using the "Environment Irrelevance test"
method <- method_obj(method = "EnvirIrrel", model = "glm", family = "gaussian")
ICP(SIM$Y, SIM$X, SIM$E, method, level = 0.05)

# analyze using the old marginal confidence interval test
method <- method_obj(method = "CRrectangle", model = "glm", family = "gaussian")
ICP(SIM$Y, SIM$X, SIM$E, method, level = 0.05)

# analyze using the new "Confidence Region" test with correct ellipsoid CR's
method <- method_obj(method = "CRellipsoid", model = "glm", family = "gaussian")
ICP(SIM$Y, SIM$X, SIM$E, method, level = 0.05)
# ---------------------------------------------------------------------------





#  ===  E X A M P L E    2  ====================================================
# DISTRIBUTION OF Y :  poisson
# TIME DEPENDENT    :  NO
# TRUE PARENTS      :  X1 and X2
# NUMBER OF OBS.    :  100

# ENVIRONMENT DETAILES:
#   - number of environments :  3
#   - effect size on mean    :  0, 2, 5

# RESULT OF ICP ANALYSIS:
#   - EnvirIrrel  :  X1 + X2 (correct)
#   - CRrectangle :  X1 + X2 (correct)
#   - CRellipsoid :  X1 + X2 (correct)

set.seed(123)

n <- 100  # number of data points

# First we simulate data.
G <- matrix(c(0,1,1,1,0,0,0,0,0,0.4,0,0,0,0,0.4,0,0,0,0,0,0,0,0,0,0),
             ncol = 5, byrow = T,
             dimnames = list(c("E","X1","X2","X3","Y"), c("E","X1","X2","X3","Y")))
sim <- c("sample(c(0,2,5),N,replace=T)", # specify Environment
         rep("rnorm(N,BETA,1)",3),       # specify X's
         "rpois(N,exp(BETA))")           # specify Y
out <- list(X = 2:4, Y = 5, E = 1, full = 1:5)
SIM <- sim_from_adj(G, n, sim, out)

# If one wants to check if Y is simulated correctly
cbind("TRUE" = G[c("X1", "X2"),"Y"],
      stats::confint(stats::glm(Y ~ X1 + X2 - 1, family = "poisson", data = SIM$full)))



# analyze using the "Environment Irrelevance test"
method <- method_obj(method = "EnvirIrrel", model = "glm", family = "poisson")
ICP(SIM$Y, SIM$X, SIM$E, method, level = 0.05)

# analyze using the old marginal confidence interval test
method <- method_obj(method = "CRrectangle", model = "glm", family = "poisson")
ICP(SIM$Y, SIM$X, SIM$E, method, level = 0.05)

# analyze using the new "Confidence Region" test with correct ellipsoid CR's
method <- method_obj(method = "CRellipsoid", model = "glm", family = "poisson")
ICP(SIM$Y, SIM$X, SIM$E, method, level = 0.05)
# ---------------------------------------------------------------------------







#  ===  E X A M P L E    3  ====================================================
# DISTRIBUTION OF Y :  hazard given by proportional weibull distribution (we use coxph)
# TIME DEPENDENT    :  NO
# TRUE PARENTS      :  X1 and X2
# NUMBER OF OBS.    :  1000

# ENVIRONMENT DETAILES:
#   - number of environments :  3
#   - effect size on mean    :  0, 2, 5

# RESULT OF ICP ANALYSIS:
#   - EnvirIrrel  :  X1 + X2 (correct)
#   - CRrectangle :  EMPTY   (wrong)
#   - CRellipsoid :  EMPTY   (wrong)

library(survival)
set.seed(123)

n <- 1000  # number of data points

# First we simulate data.
G <- matrix(c(0, 0, 0.1, 0.1,   0, 0,
              0, 0,   0,   0, 0.1, 0,
              0, 0,   0,   0,   0, 3,
              0, 0,   0,   0,   0, 3,
              0, 0,   0,   0,   0, 0,
              0, 0,   0,   0,   0, 0),
            ncol = 6, byrow = T,
            dimnames = list(NULL, c("E1","E2","X1","X2","X3","Y")))
sim <- c(rep("sample(c(0,5,10),N,replace=T)",2),          # specify Environment
         rep("rnorm(N,BETA,0.1)",3),                      # specify X's
         "(-log(runif(N))/(0.0001*exp(BETA)))^(1/2.7)")   # specify Y
out <- list(E = 1:2, X = 3:5, Y = 6)
SIM <- sim_from_adj(G, n, sim, out)
summary(SIM$Y)

# analyze using the "Environment Irrelevance test"
method <- method_obj(method = "EnvirIrrel", model = "survreg", dist = "weibull")
ICP(SIM$Y, SIM$X, SIM$E, method, level = 0.05)

# analyze using the "Environment Irrelevance test"
method <- method_obj(method = "EnvirIrrel",
                     model = "survreg", dist = "weibull", scale = 2.7)
ICP(SIM$Y, SIM$X, SIM$E, method, level = 0.05)



# analyze using the old marginal confidence interval test
method <- method_obj(method = "CRrectangle", model = "survreg", dist = "weibull")
ICP(SIM$Y, SIM$X, SIM$E, method, level = 0.05)

# analyze using the old marginal confidence interval test
method <- method_obj(method = "CRrectangle",
                     model = "survreg", dist = "weibull", scale = 2.7)
ICP(SIM$Y, SIM$X, SIM$E, method, level = 0.05)


# analyze using the new "Confidence Region" test with correct ellipsoid CR's
method <- method_obj(method = "CRellipsoid", model = "survreg", dist = "weibull")
ICP(SIM$Y, SIM$X, SIM$E, method, level = 0.05)


# analyze using the new "Confidence Region" test with correct ellipsoid CR's
method <- method_obj(method = "CRellipsoid",
                     model = "survreg", dist = "weibull", scale = 2.7)
ICP(SIM$Y, SIM$X, SIM$E, method, level = 0.05)
# NOTICE THAT IT ACCEPTS S = X1 and S = X2 - - WHY ??
# ---------------------------------------------------------------------------







#  ===  E X A M P L E    3  ====================================================
# DISTRIBUTION OF Y :  hazard given by proportional weibull distribution (we use coxph)
# TIME DEPENDENT    :  NO
# TRUE PARENTS      :  X1 and X2
# NUMBER OF OBS.    :  1000

# ENVIRONMENT DETAILES:
#   - number of environments :  3
#   - effect size on mean    :  0, 2, 5

# RESULT OF ICP ANALYSIS:
#   - EnvirIrrel  :  X1 + X2 (correct)
#   - CRrectangle :  EMPTY   (wrong)
#   - CRellipsoid :  EMPTY   (wrong)

library(survival)
library(timereg)
set.seed(123)

n <- 1000  # number of data points

# First we simulate data.
G <- matrix(c(0, 0, 0.1, 0.1,   0, 0,
              0, 0,   0,   0, 0.1, 0,
              0, 0,   0,   0,   0, 3,
              0, 0,   0,   0,   0, 3,
              0, 0,   0,   0,   0, 0,
              0, 0,   0,   0,   0, 0),
            ncol = 6, byrow = T,
            dimnames = list(NULL, c("E1","E2","X1","X2","X3","Y")))
sim <- c(rep("sample(c(0,5,10),N,replace=T)",2),          # specify Environment
         rep("rnorm(N,BETA,0.1)",3),                      # specify X's
         "(-log(runif(N))/(0.0001*exp(BETA)))^(1/2.7)")   # specify Y
out <- list(E = 1:2, X = 3:5, Y = 6)
SIM <- sim_from_adj(G, n, sim, out)
summary(SIM$Y)

# analyze using the "Environment Irrelevance test"
method <- method_obj(method = "EnvirIrrel", model = "coxph")
ICP(SIM$Y, SIM$X, SIM$E, method, level = 0.05)

# analyze using the old marginal confidence interval test
method <- method_obj(method = "CRrectangle", model = "coxph")
ICP(SIM$Y, SIM$X, SIM$E, method, level = 0.05)


# analyze using the new "Confidence Region" test with correct ellipsoid CR's
method <- method_obj(method = "CRellipsoid", model = "coxph")
ICP(SIM$Y, SIM$X, SIM$E, method, level = 0.05)

# analyze using the new "Confidence Region" test with correct ellipsoid CR's
method <- method_obj(method = "ConstRegParam", model = "coxph") # add n.sim ?
ICP(SIM$Y, SIM$X, E = NULL, method, level = 0.05)
# ---------------------------------------------------------------------------



#  ===  E X A M P L E    4  ====================================================
# DISTRIBUTION OF Y :  additive hazard given by weibull distribution (we use aalen)
# TIME DEPENDENT    :  NO
# TRUE PARENTS      :  X1 and X2
# NUMBER OF OBS.    :  1000

# ENVIRONMENT DETAILES:
#   - number of environments :  3
#   - effect size on mean    :  0, 2, 5

# RESULT OF ICP ANALYSIS:
#   - EnvirIrrel  :  X1 + X2 (correct)
#   - CRrectangle :  EMPTY   (wrong)
#   - CRellipsoid :  EMPTY   (wrong)

library(survival)
library(timereg)
set.seed(123)

n <- 1000  # number of data points

# First we simulate data.
G <- matrix(c(0,1,1,1,0,0,0,0,0,0.01,0,0,0,0,0.01,0,0,0,0,0,0,0,0,0,0),
            ncol = 5, byrow = T,
            dimnames = list(c("E","X1","X2","X3","Y"), c("E","X1","X2","X3","Y")))
sim <- c("sample(c(1,2,3),N,replace=T)", # specify Environment
         rep("rnorm(N,BETA,0.2)",2),       # specify X's
         "rexp(N, ifelse(BETA==2,1/10,1))",
         "(-log(runif(N))/(0.001 + BETA))^(1/1.2)")           # specify Y
out <- list(E = 1, X = 2:4, Y = 5)
SIM <- sim_from_adj(G, n, sim, out)
summary(SIM$Y)
summary(SIM$X)


# analyze using the "Environment Irrelevance test"
method <- method_obj(method = "EnvirIrrel", model = "aalen")
ICP(SIM$Y, SIM$X, SIM$E, method, level = 0.05)

# analyze using the old marginal confidence interval test
method <- method_obj(method = "CRrectangle", model = "aalen")
ICP(SIM$Y, SIM$X, SIM$E, method, level = 0.05)

# analyze using the new "Confidence Region" test with correct ellipsoid CR's
method <- method_obj(method = "CRellipsoid", model = "aalen")
ICP(SIM$Y, SIM$X, SIM$E, method, level = 0.05)

#
method <- method_obj(method = "ConstRegParam", model = "aalen", n.sim = 100)
ICP(SIM$Y, SIM$X, SIM$E, method, level = 0.05)
# ---------------------------------------------------------------------------



timeSurv <- function(BETA, time, type = "cox", scale = 1/100, shape = 2) {
  ntime <- c(0,time)
  if (length(BETA) > 1) {
    if (length(alpha) == 1) {
      alpha <- rep(alpha, length(BETA))
    }
    if (length(gamma) == 1) {
      gamma <- rep(gamma, length(BETA))
    }
  }
  res <- sapply(seq_along(time), function(i) {
    u <- - log(runif(1))
    if (type == "cox") {
      val <- ntime[i] + (u/(scale[i] * exp(BETA[i])))^(1/shape)
    } else if (type == "aalen") {
      val <- ntime[i] + (u/(scale[i] + BETA[i]))^(1/shape)
    }
    res <- ifelse(val < time[i], 1, 0)
  })

}



#  ===  E X A M P L E    4  ====================================================
# DISTRIBUTION OF Y :  hazard given by proportional weibull distribution (we use coxph)
# TIME DEPENDENT    :  YES
# TRUE PARENTS      :  X1 and X2
# NUMBER OF OBS.    :  500

# ENVIRONMENT DETAILES:
#   - number of environments :  3
#   - effect size on mean    :  0, 2, 5

# RESULT OF ICP ANALYSIS:
#   - EnvirIrrel  :  X1 + X2 (correct)
#   - CRrectangle :  EMPTY   (wrong)
#   - CRellipsoid :  EMPTY   (wrong)

library(survival)
set.seed(123)

n <- 1000  # number of data points

# First we simulate data.
G <- matrix(c(0, 0, 0.1, 0.1,   0, 0,
              0, 0,   0,   0, 0.1, 0,
              0, 0,   0,   0,   0, 3,
              0, 0,   0,   0,   0, 3,
              0, 0,   0,   0,   0, 0,
              0, 0,   0,   0,   0, 0),
            ncol = 6, byrow = T,
            dimnames = list(NULL, c("E1","E2","X1","X2","X3","Y")))
sim <- c(rep("sample(c(0,5,10),N,replace=T)",2),          # specify Environment
         rep("cumsum(rnorm(N,BETA,0.1))",3),              # specify X's (random walks)
         "tcoxph(N, BETA)")   # specify Y
out <- list(E = 1:2, X = 3:5, Y = 6)
SIM <- sim_from_adj(G, n, sim, out)
summary(SIM$Y)

# analyze using the "Environment Irrelevance test"
method <- method_obj(method = "EnvirIrrel", model = "coxph")
ICP(SIM$Y, SIM$X, SIM$E, method, level = 0.05)

# analyze using the old marginal confidence interval test
method <- method_obj(method = "CRrectangle", model = "coxph")
ICP(SIM$Y, SIM$X, SIM$E, method, level = 0.05)


# analyze using the new "Confidence Region" test with correct ellipsoid CR's
method <- method_obj(method = "CRellipsoid", model = "coxph")
ICP(SIM$Y, SIM$X, SIM$E, method, level = 0.05)
# ---------------------------------------------------------------------------



