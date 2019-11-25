context("IMR Input")

n <- 100
E <- rbinom(n, 2, 0.5)
X <- data.frame(X1 = rnorm(n, E, 1), X2 = rnorm(n, 5, 1))
Y <- rexp(n, exp(- 1.5 * X$X1))

test_that("Models not comparable", {
  fit1 <- survival::coxph(survival::Surv(Y) ~ X1, data = X)
  fit2 <- survival::coxph(survival::Surv(Y) ~ X2, data = X, subset = seq_len(n - 2))
  expect_error(IMR(list(fit1 = fit1, fit2)),
               paste0("The models are not comparable",
                      " - the respons variables are not identical."))
})


test_that("Inputting nonsens", {
  expect_error(IMR(M = list(c(1))),
               "'M' must contain at least two fitted models.")
  expect_error(
    IMR(M = list(c(1), c(1))),
    paste0("The model class 'numeric' is not recognised as a 'find_data' method."))

  expect_error(IMR(M = list(lm(rnorm(10) ~ 1), lm(rnorm(10) ~ 1)), method = "not"),
               "'not' is not a recognised method for Invariant Model Ranking.")
})
