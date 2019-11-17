context("Plausible predictor test")

test_that("fit_model", {

  n <- 500
  X <- data.frame(X1 = rnorm(n, 5,1), X2 = rnorm(n, 0, 1))
  Y <- X$X1 - 2*X$X2 + rnorm(n)
  E <- rep(1:4, length.out = n)

  model <- method_obj()
  expect_is(fit_model(model, Y, cbind(X,E)), "list")

  model <- method_obj(family = gaussian(link = "identity"))
  expect_is(fit_model(model, Y, cbind(X,E)), "list")

  Y <- - (-log(runif(n)) - (5 + X$X1 + 2*X$X2))
  model <- method_obj(model = "ah")
  expect_is(fit_model(model, Y, cbind(X,E)), "list")
  expect_is(fit_model(model, survival::Surv(Y), cbind(X,E)), "list")

  model <- method_obj(method = "ConstTime", model = "ah")
  expect_is(fit_model(model, Y, X, id = NULL), "list")


  Y <- -log(runif(n))/(exp(X$X1 - 2*X$X2))
  model <- method_obj(model = "ph")
  expect_is(fit_model(model, Y, cbind(X,E)), "list")
  expect_is(fit_model(model, survival::Surv(Y), cbind(X,E)), "list")

  model <- method_obj(method = "ConstTime", model = "ph")
  expect_is(fit_model(model, Y, X, id = NULL), "list")
})


# == lm ========================================================================

test_that("EnvirIrrel lm", {
  n <- 100
  X <- rnorm(n, 1, 0.3)
  Y <- rnorm(n, X, 0.1)
  E <- factor(rep(1:4, length.out = n))
  method <- method_obj(method = "EnvirIrrel", model = "lm")

  # -- non-empty model ---------------------------------------------------------
  val <- anova(lm(Y ~ X, data = data.frame(Y = Y, X = X)),
               lm(Y ~ X + E, data = data.frame(Y = Y, X = X, E = E)),
               test = "LRT")$`Pr(>Chi)`[2]
  expect_equal(plausible_predictor_test(method, Y, data.frame(X), E), val)

  # -- empty model -------------------------------------------------------------
  val <- anova(lm(Y ~ 1, data = data.frame(Y = Y, X = X)),
               lm(Y ~ E, data = data.frame(Y = Y, E = E)),
               test = "LRT")$`Pr(>Chi)`[2]
  expect_equal(plausible_predictor_test(method, Y, data.frame(X)[,0], E), val)
})


# == glm =======================================================================

test_that("EnvirIrrel glm", {
  n <- 100
  X <- rnorm(n, 1, 0.3)
  E <- factor(rep(1:4, length.out = n))
  Y <- rpois(n, exp(X))
  method <- method_obj(method = "EnvirIrrel", model = "glm", family = "poisson")

  # -- Non-empty model ---------------------------------------------------------
  val <- anova(glm(Y ~ X,
                   data = data.frame(Y = Y, X = X),
                   family = "poisson"),
               glm(Y ~ X + E,
                   data = data.frame(Y = Y, X = X, E = E),
                   family = "poisson"),
               test = "LRT")$`Pr(>Chi)`[2]
  expect_equal(plausible_predictor_test(method, Y, data.frame(X), E), val)

  # -- empty model -------------------------------------------------------------
  val <- anova(glm(Y ~ 1, data = data.frame(Y = Y), family = "poisson"),
               glm(Y ~ E, data = data.frame(Y = Y, E = E), family = "poisson"),
               test = "LRT")$`Pr(>Chi)`[2]
  expect_equal(plausible_predictor_test(method, Y, data.frame(X)[,0], E), val)
})



# == ph ========================================================================

test_that("EnvirIrrel ph", {
  n <- 100
  X <- rnorm(n, 1, 0.3)
  E <- factor(rep(1:4, length.out = n))
  Y <- survival::Surv(rexp(n, exp(-2 * X)))
  method <- method_obj(method = "EnvirIrrel", model = "ph")

  # -- non-empty model ---------------------------------------------------------
  val <- anova(survival::survreg(Y ~ X,
                                 data = data.frame(X = X),
                                 dist = "exponential"),
               survival::survreg(Y ~ X + E,
                                 data = data.frame(X = X, E = E),
                                 dist = "exponential"),
               test = "Chisq")$`Pr(>Chi)`[2]
  expect_equal(plausible_predictor_test(method, Y, data.frame(X), E), val)

  # -- empty model ---------------------------------------------------------
  val <- anova(survival::survreg(Y ~ 1, dist = "exponential"),
               survival::survreg(Y ~ 1 + E,
                                 data = data.frame(E = E),
                                 dist = "exponential"),
               test = "Chisq")$`Pr(>Chi)`[2]
  expect_equal(plausible_predictor_test(method, Y, data.frame(X)[,0], E), val)
})



# == ah ========================================================================
# TODO : DOES NOT WORK YET !!!!!
# test_that("EnvirIrrel ah", {
#   n <- 100
#   X <- rnorm(n, 5, 1)
#   E <- factor(rep(1:4, length.out = n))
#   Y <- survival::Surv(Y)
#   method <- method_obj(method = "EnvirIrrel", model = "ph")
#   str(method)
#   method <- structure(
#     list(method = "EnvirIrrel",
#          model = "hazard",
#          link = "proportional",
#          dist = "exponential"),
#     class = c("method_obj", "EnvirIrrel", "hazard"))
#   plausible_predictor_test(method, Y, data.frame(X), E)
#
#
#   dist <- survival::survreg.distributions$exponential
#   dist$trans <- dist$itrans <- function(y) y
#   dist$dtrans <- function(y) rep(1,length(y))
#
#   # -- non-empty model ---------------------------------------------------------
#   val <- anova(survival::survreg(Y ~ X,
#                                  data = data.frame(X = X),
#                                  dist = dist),
#                survival::survreg(Y ~ X + E,
#                                  data = data.frame(X = X, E = E),
#                                  dist = dist),
#                test = "Chisq")$`Pr(>Chi)`[2]
#
#   expect_equal(plausible_predictor_test(method, Y, data.frame(X), E), val)
#
#   # -- empty model ---------------------------------------------------------
#   val <- anova(survival::survreg(Y ~ 1, dist = "exponential"),
#                survival::survreg(Y ~ 1 + E,
#                                  data = data.frame(E = E),
#                                  dist = "exponential"),
#                test = "Chisq")$`Pr(>Chi)`[2]
#   expect_equal(plausible_predictor_test(method, Y, data.frame(X)[,0], E), val)
# })



