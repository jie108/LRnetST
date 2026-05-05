cal_order <- function(adj_matrix) {
  ## Topological order of nodes in a DAG.
  ## adj_matrix[i,j] = 1 means edge i->j.
  ## Returns a length-p permutation vector.

  p <- nrow(adj_matrix)
  if (nrow(adj_matrix) != ncol(adj_matrix))
    stop("adjacency matrix is not square")

  node_order <- NULL
  ind_p      <- numeric(p)
  csums      <- colSums(adj_matrix)

  npar_index <- which(csums == 0)
  hpar_index <- which(csums != 0)
  ind_p[npar_index] <- 1
  node_order <- c(node_order, npar_index)
  current    <- hpar_index

  while (length(current) > 0) {
    j          <- current[1]
    par_index  <- which(adj_matrix[, j] == 1)
    if (all(ind_p[par_index] == 1)) {
      node_order <- c(node_order, j)
      current    <- current[-1]
      ind_p[j]   <- 1
    } else {
      current <- c(current[-1], j)
    }
  }
  node_order
}


moral_graph <- function(adj.matrix) {
  ## Moral graph of a DAG: connect every pair of unmarried parents.
  ## adj.matrix[i,j] = 1 means edge i->j.

  p                <- nrow(adj.matrix)
  moral.adj.matrix <- adj.matrix
  n.parents        <- colSums(adj.matrix)

  for (i in seq_len(p)) {
    if (n.parents[i] >= 2) {
      pars <- which(adj.matrix[, i] > 0)
      for (j in seq_len(n.parents[i] - 1)) {
        for (k in (j + 1):n.parents[i]) {
          moral.adj.matrix[pars[j], pars[k]] <- 1
        }
      }
    }
  }

  moral.adj.matrix <- ((moral.adj.matrix + t(moral.adj.matrix)) > 0)
  moral.adj.matrix
}


vstructures <- function(adj.matrix) {
  ## Find all v-structures (unshielded colliders) in a DAG.
  ## Returns a 3-column matrix (par1, child, par2) with par1 < par2,
  ## or NULL if there are none.

  p   <- ncol(adj.matrix)
  res <- NULL

  for (i in seq_len(p)) {
    parent.node <- which(adj.matrix[, i] != 0)
    n.par       <- length(parent.node)

    if (n.par > 1) {
      for (j in seq_len(n.par - 1)) {
        for (k in (j + 1):n.par) {
          if (adj.matrix[parent.node[j], parent.node[k]] == 0 &&
              adj.matrix[parent.node[k], parent.node[j]] == 0)
            res <- rbind(res, c(parent.node[j], i, parent.node[k]))
        }
      }
    }
  }

  if (!is.null(res)) colnames(res) <- c("par1", "child", "par2")
  res
}


skeleton <- function(adj.matrix) {
  ## Skeleton (undirected) graph of a DAG.
  (adj.matrix + t(adj.matrix)) > 0
}


compare.vstructures <- function(target.vstructures, true.vstructures) {
  ## Subset of true.vstructures that also appear in target.vstructures.

  corr.v <- NULL
  if (!is.null(target.vstructures)) {
    target.vstructures <- matrix(target.vstructures, ncol = 3)
    for (i in seq_len(nrow(target.vstructures))) {
      res <- apply(true.vstructures, 1,
                   function(x) all(x == target.vstructures[i, ]))
      if (any(res))
        corr.v <- rbind(corr.v, true.vstructures[which(res), ])
    }
  }
  corr.v
}
