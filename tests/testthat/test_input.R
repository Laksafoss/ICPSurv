context("ICP Input")


test_that("method_obj", {
  # standard
  expect_is(method_obj(), c("method_obj","EnvirIrrel","glm"))

  # EnvirIrrel
  mf <- method_obj(model = "lm", method = "EnvirIrrel")
  expect_is(mf, c("method_obj","EnvirIrrel","lm"))

  mf <- method_obj(model = "glm", method = "EnvirIrrel")
  expect_is(mf, c("method_obj","EnvirIrrel","glm"))
  expect_equal(mf$family, "gaussian")

  mf <- method_obj(model = "glm", method = "EnvirIrrel", family = "poisson")
  expect_is(mf, c("method_obj","EnvirIrrel","glm"))
  expect_equal(mf$family, "poisson")

  mf <- method_obj(model = "glm",
                   method = "EnvirIrrel",
                   family = poisson(link = "log"))
  expect_is(mf, c("method_obj","EnvirIrrel","glm"))
  expect_is(mf$family, "family")

  mf <- method_obj(model = "hazard", method = "EnvirIrrel")
  expect_is(mf, c("method_obj", "EnvirIrrel", "hazard"))
  expect_is(mf$dist, "list")
  expect_is(mf$link, "character")

  # CR
  mf <- method_obj(method = "CR")
  expect_is(mf, c("method_obj", "CR", "glm", "QCLP"))
  expect_is(mf$solver, "character")
  expect_is(mf$solver, "character")

  # nonparam
  mf <- method_obj(model = "hazard", method = "nonparam")
  expect_is(mf, c("method_obj","nonparam","hazard"))
  expect_is(mf$link, "character")
  expect_is(mf$dist, "list")
  expect_equal(mf$dist$dist, "extreme")
  expect_equal(mf$nonparamtest, "both")
  expect_equal(mf$n.sim, 50)

  expect_error(method_obj(method = "nonparam"),
               "A 'nonparamtest' for glm is not yet avalible in this package")
})

n <- 20
E <- rep(c(0,1), each = n/2)
X <- data.frame(X1 = rnorm(n, E,1), X2 = rnorm(n, E==0, 1))
Y <- X$X1 - 2*X$X2 + rnorm(n)


test_that("X is vector, matrix or data.frame", {
  expect_is(ICP(Y, X, E), "ICP")

  Xn <- as.matrix(X)
  expect_is(ICP(Y, Xn, E), "ICP")

  Xn <- X[,1]
  expect_is(ICP(Y, Xn, E), "ICP")

  Xn[X[,1] <= -2.9] <- "b"
  Xn[X[,1] >=  1.0] <- "c"
  expect_is(ICP(Y, Xn, E), "ICP")

  expect_is(ICP(Y, X = list(X1 = X[,1], X2 = X[,2]), E), "ICP")
  expect_error(ICP(Y, X = list(X1 = X[,1], X2 = X[1:(n-2),2]), E),
               "When 'X' is a vector it must be coercible to data.frame.")
  expect_error(ICP(Y, X = matrix((n + 5) * 5, ncol = 5), E),
               "'E' and 'X' must have same length / number of rows")
  expect_error(ICP(Y, X = 1:(n + 5), E),
               "'E' and 'X' must have same length / number of rows")
  expect_error(ICP(Y, X = NULL, E), "'X' must be a vector, matrix or data frame")

})

test_that("maxNoVariables test", {
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

test_that("E is a propper enviroment vector", {
  En <- data.frame("E1" = E, "E2" = rep(c(1,2,3,4), each = n/4))
  expect_is(ICP(Y, X, En), "ICP")
  expect_is(ICP(Y, X, E = as.matrix(En)), "ICP")
  expect_is(ICP(Y, X, E = as.factor(E)), "ICP")
  expect_error(ICP(Y, X, E = list("E1" = En$E1, "E2" = En$E2)),
               "'E' must be vector, matrix, data frame or factor")
  expect_error(ICP(Y, X, E = rep(1, n + 5)),
               "'E' and 'X' must have same length / number of rows")
  expect_error(ICP(Y, X, E = c(1, rep(2, n - 1))), NULL)
})
