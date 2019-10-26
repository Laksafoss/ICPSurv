context("Plausible predictor test")

test_that("fit_model", {

  n <- 100
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


