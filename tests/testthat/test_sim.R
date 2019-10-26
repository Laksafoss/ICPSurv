contex("Simulation function")

test_that("simulation ", {
  n <- 100
  G <- matrix(c(0,1,1,0,
                0,0,0,1,
                0,0,0,0,
                0,0,0,0),
              ncol = 4, byrow = T,
              dimnames = list(c("E","X1","X2","Y"),c("E","X1","X2","Y")))

  sim <- c("sample(c(0,1), N, replace = T)",
           rep("rnorm(N, BETA, 1)",2),
           "rnorm(N, BETA, 0.5)")

  expect_is(sim_from_adj(G, n, sim), "data.frame")
  expect_error(sim_from_adj(G, n), "Both 'G', 'n' and 'sim' must be supplied")
  expect_error(sim_from_adj(G, n, m = 2), "Both 'G', 'n' and 'sim' must be supplied")
  expect_is(sim_from_adj(G, n, m = 2, sim), "data.frame")


  out <- list(Y = 4, X = 2:3, E = 1)
  expect_is(sim_from_adj(G, n, sim, out), "list")
  expect_is(sim_from_adj(G, n, m = 2, sim, out), "list")

  sim <- c("sample(c(0,1), N, replace = T)",
           "rnorm(N, BETA, 1)")
  expect_is(sim_from_adj(G, n, sim), "data.frame")

  sim <- c("sample(c(0,1), 30, replace = T)",
           "rnorm(N, BETA, 1)")
  expect_error(sim_from_adj(G, n, sim), NULL)

  sim <- c("sample(c(0,1), N, replace = T",
           "rnorm(N, BETA, 1)",
           "rnorm(N, BETA, 0.5)")
  expect_error(sim_from_adj(G, n, sim), NULL)

})
