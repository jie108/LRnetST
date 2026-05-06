make_sc_data <- function(n = 80L, p = 3L) {
  set.seed(42)
  Z <- matrix(rnorm(n * p), n, p)
  U <- matrix(rbinom(n * p, 1L, plogis(Z)), n, p)
  list(Z = Z, U = U)
}

is_dag_matrix <- function(adj) {
  pp <- nrow(adj)
  indeg <- colSums(adj != 0)
  queue <- which(indeg == 0)
  seen  <- 0L
  while (length(queue) > 0L) {
    node  <- queue[1L]; queue <- queue[-1L]; seen <- seen + 1L
    for (child in which(adj[node, ] != 0)) {
      indeg[child] <- indeg[child] - 1L
      if (indeg[child] == 0L) queue <- c(queue, child)
    }
  }
  seen == pp
}

# -- SC_prepare ----------------------------------------------------------------

test_that("SC_prepare returns correct structure", {
  dat <- make_sc_data()
  prep <- SC_prepare(dat$Z)
  expect_equal(ncol(prep$Y), 2L * ncol(dat$Z))
  expect_true(all(prep$Y[, (ncol(dat$Z) + 1L):(2L * ncol(dat$Z))] %in% c(0, 1)))
  expect_equal(sum(prep$whiteList), ncol(dat$Z))
})

# -- hcSC basic ----------------------------------------------------------------

test_that("hcSC returns a valid DAG on all-continuous data", {
  dat <- make_sc_data()
  fit <- hcSC(dat$Z, nodeType = rep("c", 3L), maxStep = 20L, seed = 1L)
  expect_true(is_dag_matrix(fit$adjacency))
})

test_that("hcSC with SC_prepare returns a valid DAG", {
  dat <- make_sc_data()
  prep <- SC_prepare(dat$Z)
  fit <- hcSC(prep$Y, nodeType = prep$nodeType,
              whiteList = prep$whiteList, blackList = prep$blackList,
              maxStep = 20L, seed = 1L)
  expect_true(is_dag_matrix(fit$adjacency))
})

# -- hcSC_boot sequential ------------------------------------------------------

test_that("hcSC_boot sequential produces valid bootstrap DAGs", {
  dat <- make_sc_data()
  prep <- SC_prepare(dat$Z)
  arr <- hcSC_boot(prep$Y, n.boot = 3L,
                   nodeType = prep$nodeType,
                   whiteList = prep$whiteList, blackList = prep$blackList,
                   maxStep = 10L, seed = 7L, output_type = "array")
  expect_equal(dim(arr), c(6L, 6L, 3L))
  for (b in seq_len(3L)) expect_true(is_dag_matrix(arr[, , b]))
})

test_that("hcSC_boot 'array' and 'freq' modes agree on frequencies", {
  dat <- make_sc_data()
  prep <- SC_prepare(dat$Z)
  both <- hcSC_boot(prep$Y, n.boot = 4L,
                    nodeType = prep$nodeType,
                    whiteList = prep$whiteList, blackList = prep$blackList,
                    maxStep = 10L, seed = 5L, output_type = "both")
  freq_from_array <- apply(both$adjacency, c(1, 2), mean)
  expect_equal(both$freq, freq_from_array)
})

# -- sequential vs future identity --------------------------------------------

test_that("hcSC_boot sequential and future backends produce identical results", {
  skip_if_not_installed("future")
  dat <- make_sc_data()
  prep <- SC_prepare(dat$Z)
  seq_out <- hcSC_boot(prep$Y, n.boot = 4L,
                       nodeType = prep$nodeType,
                       whiteList = prep$whiteList, blackList = prep$blackList,
                       maxStep = 10L, seed = 11L, backend = "sequential",
                       output_type = "array")
  fut_out <- hcSC_boot(prep$Y, n.boot = 4L,
                       nodeType = prep$nodeType,
                       whiteList = prep$whiteList, blackList = prep$blackList,
                       maxStep = 10L, seed = 11L, backend = "future", workers = 2L,
                       output_type = "array")
  expect_identical(seq_out, fut_out)
})
