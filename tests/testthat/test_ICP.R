context("ICP Results")

n <- 500
E <- rbinom(n, 2, 0.5)
X <- data.frame(X1 = rnorm(n, E, 0.5),
                X2 = rbinom(n, 4, (3 * E + 1) / 8),
                X3 = rbinom(n, E + 1, 0.5))
eta <- X$X1 + 2 * X$X2
X$X3 <- factor(X$X3)


## Linear Models
Y <- rnorm(n, eta, 0.1)

test_that("EnvirRel lm", {
  p1 <- stats::anova(stats::lm(Y ~ 1),
                     stats::lm(Y ~ factor(E)),
                     test = "LRT")$`Pr(>Chi)`[2]
  p2 <- stats::anova(stats::lm(Y ~ X$X1),
                     stats::lm(Y ~ factor(E) + X$X1),
                     test = "LRT")$`Pr(>Chi)`[2]
  p3 <- stats::anova(stats::lm(Y ~ X$X2),
                     stats::lm(Y ~ factor(E) + X$X2),
                     test = "LRT")$`Pr(>Chi)`[2]
  analysis <- ICP(Y, X, E, model = "lm", fullAnalysis = TRUE)
  expect_is(analysis, "ICP")
  expect_is(analysis$method, c("method_obj", "EnvirIrrel", "lm"))
  expect_is(analysis$model.analysis, "data.frame")
  expect_equal(analysis$model.analysis$pval[1:3], c(p1, p2, p3))
})

test_that("CR lm", {
  expect_is(ICP(Y, X, E, model = "lm", method = "CR", solver = "marginal"),
            "ICP")
  expect_is(ICP(Y, X, E, model = "lm", method = "CR", solver = "marginal",
                splits = "LOO"), "ICP")
  expect_is(ICP(Y, X, E, model = "lm", method = "CR", solver = "marginal",
                fullAnalysis = TRUE), "ICP")
  expect_is(ICP(Y, X, E, model = "lm", method = "CR", solver = "pairwise"),
            "ICP")
  expect_is(ICP(Y, X, E, model = "lm", method = "CR", solver = "pairwise",
                splits = "LOO"), "ICP")
  expect_is(ICP(Y, X, E, model = "lm", method = "CR", solver = "pairwise",
                fullAnalysis = TRUE), "ICP")
  expect_is(ICP(Y, X, E, model = "lm", method = "CR", solver = "QC"), "ICP")
  expect_is(ICP(Y, X, E, model = "lm", method = "CR", solver = "QC",
                splits = "LOO"),"ICP")
  expect_is(ICP(Y, X, E, model = "lm", method = "CR", solver = "QC",
                fullAnalysis = TRUE),"ICP")
})

## Generalized Linear Models -- Poisson, log link
Y <- rpois(n, exp(eta))

test_that("EnvirRel glm", {
  p1 <- stats::anova(stats::glm(Y ~ 1, family = "poisson"),
                     stats::glm(Y ~ factor(E), family = "poisson"),
                     test = "LRT")$`Pr(>Chi)`[2]
  p2 <- stats::anova(stats::glm(Y ~ X$X1, family = "poisson"),
                     stats::glm(Y ~ factor(E) + X$X1, family = "poisson"),
                     test = "LRT")$`Pr(>Chi)`[2]
  p3 <- stats::anova(stats::glm(Y ~ X$X2, family = "poisson"),
                     stats::glm(Y ~ factor(E) + X$X2, family = "poisson"),
                     test = "LRT")$`Pr(>Chi)`[2]
  analysis <- ICP(Y, X, E, model = "glm", family = "poisson", fullAnalysis = TRUE)
  expect_is(analysis, "ICP")
  expect_is(analysis$method, c("method_obj", "EnvirIrrel", "glm"))
  expect_is(analysis$model.analysis, "data.frame")
  expect_equal(analysis$model.analysis$pval[1:3], c(p1, p2, p3))
})

test_that("CR glm", {
  expect_is(ICP(Y, X, E, model = "glm", family = "poisson", method = "CR",
                solver = "marginal"), "ICP")
  expect_is(ICP(Y, X, E, model = "glm", family = "poisson", method = "CR",
                solver = "marginal", splits = "LOO"), "ICP")
  expect_is(ICP(Y, X, E, model = "glm", family = "poisson", method = "CR",
                solver = "marginal", fullAnalysis = TRUE), "ICP")
  expect_is(ICP(Y, X, E, model = "glm", family = "poisson", method = "CR",
                solver = "pairwise"), "ICP")
  expect_is(ICP(Y, X, E, model = "glm", family = "poisson", method = "CR",
                solver = "pairwise", splits = "LOO"), "ICP")
  expect_is(ICP(Y, X, E, model = "glm", family = "poisson", method = "CR",
                solver = "pairwise", fullAnalysis = TRUE), "ICP")
  expect_is(ICP(Y, X, E, model = "glm", family = "poisson", method = "CR",
                solver = "QC"), "ICP")
  expect_is(ICP(Y, X, E, model = "glm", family = "poisson", method = "CR",
                solver = "QC", splits = "LOO"),"ICP")
  expect_is(ICP(Y, X, E, model = "glm", family = "poisson", method = "CR",
                solver = "QC", fullAnalysis = TRUE),"ICP")
})



## Hazard Model -- Proportional
Y <- rexp(n, exp(2 - 0.5 * eta))
C <- rexp(n, exp(-3))
Y <- survival::Surv(pmin(Y,C), pmin(Y,C) == Y)

test_that("EnvirRel ph", {
  p1 <- stats::anova(survival::survreg(Y ~ 1, dist = "exponential"),
                     survival::survreg(Y ~ factor(E), dist = "exponential"),
                     test = "Chisq")$`Pr(>Chi)`[2]
  p2 <- stats::anova(survival::survreg(Y ~ X$X1, dist = "exponential"),
                     survival::survreg(Y ~ factor(E) + X$X1, dist = "exponential"),
                     test = "Chisq")$`Pr(>Chi)`[2]
  p3 <- stats::anova(survival::survreg(Y ~ X$X2, dist = "exponential"),
                     survival::survreg(Y ~ factor(E) + X$X2, dist = "exponential"),
                     test = "Chisq")$`Pr(>Chi)`[2]
  analysis <- ICP(Y, X, E, model = "ph", fullAnalysis = TRUE)
  expect_is(analysis, "ICP")
  expect_is(analysis$method, c("method_obj", "EnvirIrrel", "ph"))
  expect_is(analysis$model.analysis, "data.frame")
  expect_equal(analysis$model.analysis$pval[1:3], c(p1, p2, p3))
})

test_that("CR ph", {
  expect_is(ICP(Y, X, E, model = "ph", method = "CR", solver = "marginal"),
            "ICP")
  expect_is(ICP(Y, X, E, model = "ph", method = "CR", solver = "marginal",
                splits = "LOO"), "ICP")
  expect_is(ICP(Y, X, E, model = "ph", method = "CR", solver = "marginal",
                fullAnalysis = TRUE), "ICP")
  expect_is(ICP(Y, X, E, model = "ph", method = "CR", solver = "pairwise"),
            "ICP")
  expect_is(ICP(Y, X, E, model = "ph", method = "CR", solver = "pairwise",
                splits = "LOO"), "ICP")
  expect_is(ICP(Y, X, E, model = "ph", method = "CR", solver = "pairwise",
                fullAnalysis = TRUE), "ICP")
  expect_is(ICP(Y, X, E, model = "ph", method = "CR", solver = "QC"), "ICP")
  expect_is(ICP(Y, X, E, model = "ph", method = "CR", solver = "QC",
                splits = "LOO"),"ICP")
  expect_is(ICP(Y, X, E, model = "ph", method = "CR", solver = "QC",
                fullAnalysis = TRUE),"ICP")
})


## Hazard Model -- Additive
#Y <- rexp(n, 0.1 - min(eta) + eta)
#C <- rexp(n, 0.2)
#Y <- survival::Surv(pmin(Y,C), pmin(Y,C) == Y)
#dist <- survival::survreg.distributions$exponential
#dist$trans <- dist$itrans <- function(y) y
#dist$dtrans <- function(y) rep(1,length(y))

#test_that("EnvirRel ah", {
#  p1 <- stats::anova(survival::survreg(Y ~ 1, dist = dist),
#                     survival::survreg(Y ~ factor(E), dist = dist),
#                     test = "Chisq")$`Pr(>Chi)`[2]
#  p2 <- stats::anova(survival::survreg(Y ~ X$X1, dist = dist),
#                     survival::survreg(Y ~ factor(E) + X$X1, dist = dist),
#                     test = "Chisq")$`Pr(>Chi)`[2]
#  p3 <- stats::anova(survival::survreg(Y ~ X$X2, dist = dist),
#                     survival::survreg(Y ~ factor(E) + X$X2, dist = dist),
#                     test = "Chisq")$`Pr(>Chi)`[2]
#  analysis <- ICP(Y, X, E, model = "ah", fullAnalysis = TRUE)
#  expect_is(analysis, "ICP")
#  expect_is(analysis$method, c("method_obj", "EnvirIrrel", "ah"))
#  expect_is(analysis$model.analysis, "data.frame")
#  expect_equal(analysis$model.analysis$pval[1:3], c(p1, p2, p3))
#})


#test_that("CR ah", {
#  expect_is(ICP(Y, X, E, model = "ah", method = "CR", solver = "marginal"),
#            "ICP")
#  expect_is(ICP(Y, X, E, model = "ah", method = "CR", solver = "marginal",
#                splits = "LOO"), "ICP")
#  expect_is(ICP(Y, X, E, model = "ph", method = "CR", solver = "marginal",
#                fullAnalysis = TRUE), "ICP")
#  expect_is(ICP(Y, X, E, model = "ah", method = "CR", solver = "pairwise"),
#            "ICP")
#  expect_is(ICP(Y, X, E, model = "ah", method = "CR", solver = "pairwise",
#                splits = "LOO"), "ICP")
#  expect_is(ICP(Y, X, E, model = "ah", method = "CR", solver = "pairwise",
#                fullAnalysis = TRUE), "ICP")
#  expect_is(ICP(Y, X, E, model = "ah", method = "CR", solver = "QC"), "ICP")
#  expect_is(ICP(Y, X, E, model = "ah", method = "CR", solver = "QC",
#                splits = "LOO"),"ICP")
#  expect_is(ICP(Y, X, E, model = "ah", method = "CR", solver = "QC",
#                fullAnalysis = TRUE),"ICP")
#})
