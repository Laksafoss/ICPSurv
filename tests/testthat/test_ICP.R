context("ICP Results")

n <- 100
E <- rbinom(n, 2, 0.5)
X <- rnorm(n, E, 1)

test_that("EnvirIrrel lm", {
  Y <- rnorm(n, X, 1)
  p1 <- stats::anova(stats::lm(Y ~ 1),
                     stats::lm(Y ~ factor(E)),
                     test = "LRT")$`Pr(>Chi)`[2]
  p2 <- stats::anova(stats::lm(Y ~ X),
                     stats::lm(Y ~ factor(E) + X),
                     test = "LRT")$`Pr(>Chi)`[2]
  analysis <- ICP(Y, X, E, model = "lm")
  expect_is(analysis, "ICP")
  expect_is(analysis$method, c("method_obj", "EnvirIrrel", "lm"))
  expect_is(analysis$model.analysis, "data.frame")
  expect_equal(analysis$model.analysis$pval, c(p1, p2))
})

test_that("EnvirIrrel glm", {
  Y <- rnorm(n, X, 1)
  p1 <- stats::anova(stats::glm(Y ~ 1, family = "gaussian"),
                     stats::glm(Y ~ factor(E), family = "gaussian"),
                     test = "LRT")$`Pr(>Chi)`[2]
  p2 <- stats::anova(stats::glm(Y ~ X, family = "gaussian"),
                     stats::glm(Y ~ factor(E) + X, family = "gaussian"),
                     test = "LRT")$`Pr(>Chi)`[2]
  analysis <- ICP(Y, X, E)
  expect_is(analysis, "ICP")
  expect_is(analysis$method, c("method_obj", "EnvirIrrel", "glm"))
  expect_is(analysis$model.analysis, "data.frame")
  expect_equal(analysis$model.analysis$pval, c(p1, p2))


  Y <- rpois(n, exp(X))
  p1 <- stats::anova(stats::glm(Y ~ 1, family = "poisson"),
                     stats::glm(Y ~ factor(E), family = "poisson"),
                     test = "LRT")$`Pr(>Chi)`[2]
  p2 <- stats::anova(stats::glm(Y ~ X, family = "poisson"),
                     stats::glm(Y ~ factor(E) + X, family = "poisson"),
                     test = "LRT")$`Pr(>Chi)`[2]
  analysis <- ICP(Y, X, E, model = "glm", family = "poisson")
  expect_is(analysis, "ICP")
  expect_is(analysis$method, c("method_obj", "EnvirIrrel", "glm"))
  expect_is(analysis$model.analysis, "data.frame")
  expect_equal(analysis$model.analysis$pval, c(p1, p2))

  analysis <- ICP(Y, X, E, model = "glm", family = poisson(link = "log"))
  expect_is(analysis, "ICP")
  expect_is(analysis$method, c("method_obj", "EnvirIrrel", "glm"))
  expect_is(analysis$model.analysis, "data.frame")
  expect_equal(analysis$model.analysis$pval, c(p1, p2))
})

test_that("EnvirIrrel ph and hazard log link",  {
  Y <- rexp(n, exp(X))
  p1 <- anova(survival::survreg(survival::Surv(Y) ~ 1, dist = "exponential"),
              survival::survreg(survival::Surv(Y) ~ factor(E), dist = "exponential"),
              test = "Chisq")$`Pr(>Chi)`[2]
  p2 <- anova(survival::survreg(survival::Surv(Y) ~ X, dist = "exponential"),
              survival::survreg(survival::Surv(Y) ~ factor(E) + X, dist = "exponential"),
              test = "Chisq")$`Pr(>Chi)`[2]
  analysis <- ICP(Y, X, E, model = "hazard")
  expect_is(analysis, "ICP")
  expect_is(analysis$model.analysis, "data.frame")
  expect_equal(analysis$model.analysis$pval, c(p1, p2))

  analysis <- ICP(Y, X, E, model = "hazard", link = "log")
  expect_is(analysis, "ICP")
  expect_is(analysis$model.analysis, "data.frame")
  expect_equal(analysis$model.analysis$pval, c(p1, p2))


  analysis <- ICP(Y, X, E, model = "hazard", link = "proportional")
  expect_is(analysis, "ICP")
  expect_is(analysis$model.analysis, "data.frame")
  expect_equal(analysis$model.analysis$pval, c(p1, p2))

  analysis <- ICP(Y, X, E, model = "hazard",
                  dist = survival::survreg.distributions$exponential)
  expect_is(analysis, "ICP")
  expect_is(analysis$model.analysis, "data.frame")
  expect_equal(analysis$model.analysis$pval, c(p1, p2))

  analysis <- ICP(Y, X, E, model = "ph")
  expect_is(analysis, "ICP")
  expect_is(analysis$model.analysis, "data.frame")
  expect_equal(analysis$model.analysis$pval, c(p1, p2))

  Y <- survival::Surv(Y)
  analysis <- ICP(Y, X, E, model = "hazard")
  expect_is(analysis, "ICP")
  expect_is(analysis$model.analysis, "data.frame")
  expect_equal(analysis$model.analysis$pval, c(p1, p2))
})

test_that("CR lm", {
  Y <- rnorm(n, X, 1)
  expect_is(ICP(Y, X, E,
                model = "lm",
                method = "CR", solver = "marginal", splits = "all",
                level = 0.05),
            "ICP")
  expect_is(ICP(Y, X, E,
                model = "lm",
                method = "CR", solver = "marginal", splits = "LOO",
                level = 0.05),
            "ICP")
  expect_is(ICP(Y, X, E,
                model = "lm",
                method = "CR", solver = "pairwise", splits = "all",
                level = 0.05),
            "ICP")
  expect_is(ICP(Y, X, E,
                model = "lm",
                method = "CR", solver = "pairwise", splits = "LOO",
                level = 0.05),
            "ICP")
  expect_is(ICP(Y, X, E,
                model = "lm",
                method = "CR", solver = "QCLP", splits = "all",
                level = 0.05),
            "ICP")
  expect_is(ICP(Y, X, E,
                model = "lm",
                method = "CR", solver = "QCLP", splits = "LOO",
                level = 0.05),
            "ICP")
})

test_that("CR unobserved factors in envirs", {
  X2 <- sapply(E, function(e) {sample(seq_len(e + 2), 1)})
  Y <- rnorm(n, X, 1)
  expect_is(ICP(Y, X = data.frame(X1 = X, X2 = factor(X2)), E,
                model = "lm",
                method = "CR", solver = "marginal", splits = "all",
                level = 0.05),
            "ICP")
})

test_that("CR glm", {
  Y <- rpois(n, exp(X))
  expect_is(ICP(Y, X, E,
                model = "glm", family = "poisson",
                method = "CR", solver = "marginal", splits = "all",
                level = 0.05),
            "ICP")
  expect_is(ICP(Y, X, E,
                model = "glm", family = "poisson",
                method = "CR", solver = "marginal", splits = "LOO",
                level = 0.05),
            "ICP")
  expect_is(ICP(Y, X, E,
                model = "glm", family = "poisson",
                method = "CR", solver = "pairwise", splits = "all",
                level = 0.05),
            "ICP")
  expect_is(ICP(Y, X, E,
                model = "glm", family = "poisson",
                method = "CR", solver = "pairwise", splits = "LOO",
                level = 0.05),
            "ICP")
  expect_is(ICP(Y, X, E,
                model = "glm", family = "poisson",
                method = "CR", solver = "QCLP", splits = "all",
                level = 0.05),
            "ICP")
  expect_is(ICP(Y, X, E,
                model = "glm", family = "poisson",
                method = "CR", solver = "QCLP", splits = "LOO",
                level = 0.05),
            "ICP")
})

test_that("CR ph", {
  Y <- rexp(n, exp(X))
  expect_is(ICP(Y, X, E,
                model = "ph",
                method = "CR", solver = "marginal", splits = "all",
                level = 0.05),
            "ICP")
  expect_is(ICP(Y, X, E,
                model = "ph",
                method = "CR", solver = "marginal", splits = "LOO",
                level = 0.05),
            "ICP")
  expect_is(ICP(Y, X, E,
                model = "ph",
                method = "CR", solver = "pairwise", splits = "all",
                level = 0.05),
            "ICP")
  expect_is(ICP(Y, X, E,
                model = "ph",
                method = "CR", solver = "pairwise", splits = "LOO",
                level = 0.05),
            "ICP")
  expect_is(ICP(Y, X, E,
                model = "ph",
                method = "CR", solver = "QCLP", splits = "all",
                level = 0.05),
            "ICP")
  expect_is(ICP(Y, X, E,
                model = "ph",
                method = "CR", solver = "QCLP", splits = "LOO",
                level = 0.05),
            "ICP")
})

test_that("CR hazard", {
  Y <- rexp(n, 10 * exp(X))
  dist <- survival::survreg.distributions$exponential
  dist$trans <- function(y) { log(y / 10) }
  dist$itrans <- function(y) { 10 * exp(y) }

  expect_is(ICP(Y, X, E,
                model = "hazard", dist = dist,
                method = "CR", solver = "marginal", splits = "all",
                level = 0.05),
            "ICP")
  expect_is(ICP(Y, X, E,
                model = "hazard", link = "log",
                method = "CR", solver = "marginal", splits = "LOO",
                level = 0.05),
            "ICP")
  expect_is(ICP(Y, X, E,
                model = "hazard", link = "log",
                method = "CR", solver = "pairwise", splits = "all",
                level = 0.05),
            "ICP")
  expect_is(ICP(Y, X, E,
                model = "hazard", link = "log",
                method = "CR", solver = "pairwise", splits = "LOO",
                level = 0.05),
            "ICP")
  expect_is(ICP(Y, X, E,
                model = "hazard", link = "log",
                method = "CR", solver = "QCLP", splits = "all",
                level = 0.05),
            "ICP")
  expect_is(ICP(Y, X, E,
                model = "hazard", link = "log",
                method = "CR", solver = "QCLP", splits = "LOO",
                level = 0.05),
            "ICP")
})


n <- 1000
E <- rbinom(n, 1, 0.5)
X <- rnorm(n, 5 + E, 1)

test_that("Nonparam", {

  Y <- rexp(n, 1 /(5 * X))
  expect_is(ICP(Y, X, model = "ah",
                method = "nonparam", n.sim = 100, level = 0.05),
            "ICP")
  expect_is(ICP(Y, X, model = "hazard", link = "identity",
                method = "nonparam", n.sim = 50, level = 0.05),
            "ICP")
  expect_is(ICP(Y, X, model = "hazard", link = "additive",
                method = "nonparam", n.sim = 50, level = 0.05),
            "ICP")

  Y <- rexp(n, exp(- 0.4 *X))
  expect_is(ICP(Y, X, model = "ph",
                method = "nonparam", n.sim = 50, level = 0.05),
            "ICP")
  expect_is(ICP(Y, X, model = "hazard", link = "log",
                method = "nonparam", n.sim = 50, level = 0.05),
            "ICP")
  expect_is(ICP(Y, X, model = "hazard", link = "proportional",
                method = "nonparam", n.sim = 50, level = 0.05),
            "ICP")
})


## TODO :::: ADDITIVE HAZARD

# n <- 1000
# X <- rnorm(n, 7, 2)
# Y <- rexp(n, 1 / (1.2 * X))
# my <- survival::survreg.distributions$exponential
# my$trans <- my$itrans <- function(y) y
# my$dtrans <- function(y) rep(1, length(y))
# glm(Y ~ X - 1, family = Gamma(link = "identity"))
# survival::survreg(survival::Surv(Y) ~ X - 1, dist = my)
#
# test_that("EnvirIrrel ah and hazard identity link",  {
#   # TODO
# })
#
# test_that("CR ah", {
#   # TODO !!!!
# })

