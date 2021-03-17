#' Impute multilevel dataset
#'
#' This function returns a dataframe with a complete dataset, where the missing values are imputed
#' using a bi-objective kNN method. It assumes that the class variable name is known, and the rest
#' of the variables are numerical.
#'
#' @param df_miss A dataframe with missing values
#' @param className name of the variable that contains the classes
#' @param nIter number of iterations, default = 10
#' @param alpha weight of the kNN values in the objective function, default = 0.5
#' @param k number of nearest neighbours, default = 10
#' @param distance distance function used to get the k-nearest neighbors
#' @return A dataframe with the imputed data
#' @export
biokNN.impute <- function(df_miss, className, nIter = 10, alpha = 0.5, k = 10, distance = "gower"){

  if(check.data(df_miss, className)){
    df_miss <- dplyr::select(df_miss, className, everything())
    vars <- ncol(df_miss)
    rows <- nrow(df_miss)
    neighbors <- matrix(nrow = rows, ncol = k)
    index_miss <- get.index.miss(df_miss)

    M <- mice::complete(mice::mice(df_miss, meth = "sample", m = 1, print = FALSE))
    clusters <- stats::setNames(as.matrix(M[, 1]), 1:rows)

    iter <-1
    while(iter < nIter){
      dist <- as.matrix(cluster::daisy(M, metric = distance))
      neighbors <- get.neighbors(as.matrix(dist), k)
      clusters <- stats::setNames(as.matrix(M[, 1]), 1:rows)

      for(j in 2:vars){
        for(i in index_miss[, j]){
          if(i != 0){
            clust_val <- clusters[i]
              sum_kNN <- sum(M[neighbors[i,], j])
              has_neigh <- has.x.as.neighbor(i, neighbors)
              sum_neigh <- sum(M[has_neigh, j])
              m_knn <- (sum_kNN + sum_neigh)/(k + length(has_neigh))
              m_cluster <- mean(M[as.numeric(names(clusters[clusters==clust_val])), j])
              M[i, j] <- alpha*m_knn + (1-alpha)*m_cluster
          }
        }
      }
      iter <- iter + 1
    }
    M
  } else {
    print("The dataframe is not in the correct format.")
  }
}


#' @importFrom magrittr %>%
check.data <- function(df, className){
  check <- TRUE
  aux <- df  %>% dplyr::select(-className) %>% dplyr::select_if(~!is.numeric(.x))
  if(ncol(aux) > 1){
    check <- FALSE
  }
  check
}



get.neighbors <- function(dist.matrix, kNN){
  n <- nrow(dist.matrix)
  neigh <- matrix(nrow = n, ncol = kNN)
  for(i in 1:n){
    d <- sort(dist.matrix[i,])
    neigh[i,] <- as.numeric(names(d[2:(kNN+1)]))
  }
  neigh
}


has.x.as.neighbor <- function(i, neighbors){
  vec <- c()
  n <- nrow(neighbors)
  for(j in 1:n){
    if(i %in% neighbors[j,]){
      vec <- c(vec, j)
    }
  }
  vec
}


get.index.miss <- function(W){
  index_miss_vars_W <- matrix(NA, ncol = ncol(W), nrow = nrow(W))
  for(i in 1:ncol(W)){
    index_miss_vars_W[,i] <- 1:nrow(W)*is.na(W[,i])
  }
  index_miss_vars_W
}
