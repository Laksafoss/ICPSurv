context("IMR Results")

n <- 100
E <- rbinom(n, 2, 0.5)
X <- data.frame(X1 = rnorm(n, E, 1), X2 = rnorm(n, 5, 1))

test_that("B-Spline coxph", {
  Y <- rexp(n, exp(- 1.5 * X$X1))
  fit1 <- survival::coxph(survival::Surv(Y) ~ X1, data = X)
  fit2 <- survival::coxph(survival::Surv(Y) ~ X2, data = X)
  fit3 <- survival::coxph(survival::Surv(Y) ~ X1 + X2, data = X)
  analysis <- IMR(list(fit1 = fit1, fit2, fit3 = fit3))
  expect_is(analysis, "IMR")
  expect_is(analysis$ranking, "data.frame")
  expect_is(analysis$method, "data.frame")
  expect_is(analysis$call, "call")
})

test_that("B-Spline timecox", {
  Y <- rexp(n, exp(- 1.5 * X$X1))
  fit1 <- ICPSurv:::quiet(timereg::timecox(survival::Surv(Y) ~ X1, data = X))
  fit2 <- ICPSurv:::quiet(timereg::timecox(survival::Surv(Y) ~ X2, data = X))
  fit3 <- ICPSurv:::quiet(timereg::timecox(survival::Surv(Y) ~ X1 + X2, data = X))
  analysis <- IMR(list(fit1 = fit1, fit2, fit3 = fit3))
  expect_is(analysis, "IMR")
  expect_is(analysis$ranking, "data.frame")
  expect_is(analysis$method, "data.frame")
  expect_is(analysis$call, "call")
})

test_that("B-Spline aalen", {
  Y <- rexp(n, 1 + 0.2 * X$X2)
  fit1 <- timereg::aalen(survival::Surv(Y) ~ X1, data = X)
  fit2 <- timereg::aalen(survival::Surv(Y) ~ X2, data = X)
  fit3 <- timereg::aalen(survival::Surv(Y) ~ X1 + X2, data = X)
  analysis <- IMR(list(fit1 = fit1, fit2, fit3 = fit3))
  expect_is(analysis, "IMR")
  expect_is(analysis$ranking, "data.frame")
  expect_is(analysis$method, "data.frame")
  expect_is(analysis$call, "call")
  fit1 <- survival::aareg(survival::Surv(Y) ~ X1, data = X, nmin = 3)
  fit2 <- survival::aareg(survival::Surv(Y) ~ X2, data = X, nmin = 3)
  fit3 <- survival::aareg(survival::Surv(Y) ~ X1 + X2, data = X, nmin = 1)
  analysis <- IMR(list(fit1 = fit1, fit2, fit3 = fit3))
  expect_is(analysis, "IMR")
  expect_is(analysis$ranking, "data.frame")
  expect_is(analysis$method, "data.frame")
  expect_is(analysis$call, "call")
})


test_that("B-Spline Cox-Aalen", {
  suppressWarnings(suppressMessages(library(timereg)))
  Y <- rexp(n, (1 / (5 + 1.5 * X$X2)) * exp(- 0.5 * X$X1))
  fit1 <- ICPSurv:::quiet(timereg::cox.aalen(
    survival::Surv(Y) ~ prop(X1), data = X))
  fit2 <- ICPSurv:::quiet(timereg::cox.aalen(
    survival::Surv(Y) ~ prop(X2), data = X))
  fit3 <- ICPSurv:::quiet(timereg::cox.aalen(
    survival::Surv(Y) ~ X1 + prop(X2), data = X))
  fit4 <- ICPSurv:::quiet(timereg::cox.aalen(
    survival::Surv(Y) ~ prop(X1) + X2, data = X))
  fit5 <- ICPSurv:::quiet(timereg::cox.aalen(
    survival::Surv(Y) ~ prop(X1) + prop(X2), data = X))
  analysis <- IMR(list(fit1 = fit1, fit2, fit3 = fit3))
  expect_is(analysis, "IMR")
  expect_is(analysis$ranking, "data.frame")
  expect_is(analysis$method, "data.frame")
  expect_is(analysis$call, "call")
  detach("package:timereg", unload=TRUE)
})


