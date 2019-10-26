context("ICP Input")


test_that("method_obj", {
  expect_is(method_obj(), c("method_obj","EnvirIrrel","glm"))
  expect_is(method_obj(method = "ConstTime"), c("method_obj","ConstTime","ah"))
  expect_is(method_obj(model = "glm", family = gaussian())$family, "family")
  expect_is(method_obj(model = "glm")$family, "character")

  expect_error(method_obj(method = "TimeConst", model = "glm"),
               "'method' must be 'CR', 'EnvirIrrel' or 'ConstTime'")

  # when method_obj can allow for new methods these tests should be changed
  expect_error(method_obj(method = "dummy"),
               "'method' must be 'CR', 'EnvirIrrel' or 'ConstTime'")
  expect_error(method_obj(model = "dummy"),
               "'model' must be 'lm', 'glm', 'ph' or 'ah'")

})


test_that("maxNoVariables test", {

  n <- 100
  E <- rbinom(n, 3, 0.4)
  X <- data.frame(X1 = rnorm(n, E,1), X2 = rnorm(n, E==0, 1))
  Y <- X$X1 - 2*X$X2 + rnorm(n)

  expect_error(ICP(Y, X, E, maxNoVariables = 0),
               "'maxNoVariables' must be an integer >= 1")
  expect_error(ICP(Y, X, E, maxNoVariables = c(1,0)),
               "'maxNoVariables' must be an integer >= 1")
  expect_error(ICP(Y, X, E, maxNoVariables = -1),
               "'maxNoVariables' must be an integer >= 1")
  expect_error(ICP(Y, X, E, maxNoVariables = "a"),
               "'maxNoVariables' must be an integer >= 1")
  expect_is(ICP(Y, X, E, maxNoVariables = 2),"ICP")
})


test_that("level is between 0 and 1", {

  n <- 100
  E <- rbinom(n, 3, 0.4)
  X <- data.frame(X1 = rnorm(n, E,1), X2 = rnorm(n, E==0, 1))
  Y <- X$X1 - 2*X$X2 + rnorm(n)

  expect_is(ICP(Y, X, E), "ICP")
  expect_is(ICP(Y, X, E, level = 0.05), "ICP")
  expect_error(ICP(Y, X, E, level = "a"),
               "'level' must either be NULL or a number strictly between 0 and 1")
  expect_error(ICP(Y, X, E, level = c(1,1)),
               "'level' must either be NULL or a number strictly between 0 and 1")
  expect_error(ICP(Y, X, E, level = 0),
               "'level' must either be NULL or a number strictly between 0 and 1")
  expect_error(ICP(Y, X, E, level = 1),
               "'level' must either be NULL or a number strictly between 0 and 1")
  expect_error(ICP(Y, X, E, level = -17),
               "'level' must either be NULL or a number strictly between 0 and 1")
})




test_that("X is vector, matrix or data.frame", {

  n <- 100
  E <- rbinom(n, 3, 0.4)
  X <- data.frame(X1 = rnorm(n, E, 1), X2 = rnorm(n, E == 0, 1))
  Y <- X$X1 - 2*X$X2 + rnorm(n)

  expect_is(ICP(Y, X, E), "ICP")

  X <- matrix(c(rnorm(n, E, 1),rnorm(n, E == 0, 1)), ncol = 2)
  Y <- X[ ,1] - 2*X[,2] + rnorm(n)
  expect_is(ICP(Y, X, E), "ICP")

  X <- rnorm(n, E, 1)
  Y <-  2*X + rnorm(n)
  expect_is(ICP(Y, X, E), "ICP")

  X <- rep(c("a", "b", "c"), length.out = n)
  Y <- rnorm(n)
  Y[X=="a"] <- Y[X=="a"]
  Y[X=="c"] <- Y[X=="c"]
  expect_is(ICP(Y, X, E), "ICP")

  # SHOULD WE ALLOW X TO BE LIST LIKE BELOW
  #X <- list(X1 = rnorm(n, 5, 1), X2 = rnorm(n, 0, 1))
  #Y <- X["X1"] - 2*X["X2"] + rnorm(n)
  #expect_is(ICP(Y, X, E), "ICP")

  X <- matrix(1:80, ncol = 5)
  expect_error(ICP(Y, X, E),
               "'E' and 'X' must have same length / number of rows")

  X <- 1:80
  expect_error(ICP(Y, X, E),
               "'E' and 'X' must have same length / number of rows")

  X <- NULL
  expect_error(ICP(Y, X, E),
               "'X' must be a vector, matrix or data frame")

})



test_that("E is a propper enviroment vector", {

  n <- 200
  E <- rbinom(n, 2, 0.4)
  X <- data.frame(X1 = rnorm(n, E, 1), X2 = rnorm(n, E == 0, 1))
  Y <- X$X1 - 2*X$X2 + rnorm(n)

  expect_is(ICP(Y, X, E), "ICP")

  E <- data.frame("E1" = E,
                  "E2" = sample(c("a", "b"), n, replace = T))
  expect_is(ICP(Y, X, E), "ICP")

  E <- list("E1" = E$E1, "E2" = E$E2)
  expect_error(ICP(Y, X, E), "'E' must be a vector")

  E <- rep(1,20)
  expect_error(ICP(Y, X, E), "'E' and 'X' must have same length / number of rows")

  E <- c(1, rep(2, n - 1))
  expect_error(ICP(Y, X, E), NULL)

})


