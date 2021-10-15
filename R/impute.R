#' Impute multilevel dataset
#'
#' This function returns a dataframe with a complete dataset, where the missing values are imputed
#' using a bi-objective kNN method. It assumes that the class variable is complete and its name is known, and the rest
#' of the variables are numerical.
#'
#' @param data A dataframe with missing values
#' @param className name of the variable that contains the classes
#' @param nIter number of iterations, default = 10
#' @param alpha weight of the kNN values in the objective function, default = 0.5
#' @param k number of nearest neighbours, default = 10
#' @param distance distance function used to get the k-nearest neighbors
#' @return A dataframe with the imputed data
#' @export
#'
#' @examples
#' data(data_example)
#' complete_data <- biokNN_impute(data.example,
#'                className = "class",
#'                nIter = 10,
#'                alpha = 0.9,
#'                k = 15,
#'                distance = "gower")
biokNN_impute <- function(data, className, nIter = 10, alpha = 0.5, k = 10, distance = "gower"){
  df_miss <- data
  if(check.data(df_miss, className)){
    df_miss <- dplyr::select(df_miss, className, dplyr::everything())
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


#' Multiple imputation for a multilevel dataset
#'
#' This function returns a list of m complete datasets, where the missing values are imputed
#' using a bi-objective kNN method. It assumes that the class variable name is known, and the rest
#' of the variables are numerical.
#'
#' @param data A dataframe with missing values
#' @param className name of the variable that contains the classes
#' @param m number of imputations
#' @param nIter number of iterations, default = 10
#' @param alpha weight of the kNN values in the objective function, default = 0.5
#' @param k number of nearest neighbours, default = 10
#' @param distance distance function used to get the k-nearest neighbors
#' @return A dataframe with the imputed data
#' @export
#'
#' @examples
#' data(data_example)
#' complete_data_mi <- biokNN_impute_mi(data.example,
#'                className = "class",
#'                m = 3,
#'                nIter = 10,
#'                alpha = 0.9,
#'                k = 15,
#'                distance = "gower")
#' # View completed data sets
#' str(complete_data_mi)
biokNN_impute_mi <- function(data, className, m =5, nIter = 10, alpha = 0.5, k = 10, distance = "gower"){
  df_miss <- data
  if(check.data(df_miss, className)){

    MI <- vector(mode = "list", length = m)
    for(t in 1:m){
      df_miss <- dplyr::select(df_miss, className, dplyr::everything())
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
      MI[[t]] <- M
    }
    MI
  } else {
    print("The dataframe is not in the correct format.")
  }
}


#' Calibrate parameters
#'
#' This function returns a vector with the two parameters requiered by the biokNN method
#' where the first value is the weighting parameter and the second the number of neighbors
#'
#' @param data A dataframe with missing values
#' @param className name of the variable that contains the classes
#' @param distance distance function used to get the k-nearest neighbors
#' @param nIter number of iterations, default = 10
#' @param prop_valid proportion of missing values
#' @param alpha_space vector with the calibration values to test for the weight parameter
#' @param k_space vector with the calibration values to test for the number of neighbors
#' @param print option to print  the RMSE values of the parameters used for calibration (print = TRUE).
#' @return A dataframe with the imputed data
#' @export
#'
#' @examples
#' data(data_example)
#' calibrate(data_example,
#'           "class",
#'           prop_valid = 0.3,
#'           alpha_space = c(0.5, 0.7, 0.9),
#'           k_space = c(10, 15),
#'           print = TRUE)
calibrate <- function(data, className, prop_valid = 0.1, nIter = 10, distance = "gower", alpha_space = NULL, k_space = NULL, print = FALSE){
  df_miss <- data
  if(check.data(df_miss, className)){
    df_miss <- dplyr::select(df_miss, className, dplyr::everything())}

  if(is.null(alpha_space)){
    alpha <- seq(0, 1, by = 0.1)
  } else {
    alpha <- alpha_space
  }

  if(is.null(k_space)){
    k <- seq(5, 30, by =5)
  } else {
    k <- k_space
  }

  orig_pattern <- get.index.miss(df_miss)
  df_miss_valid <- make.missing(df_miss, make.pattern(df_miss, prop_valid))
  valid_pattern <- get.valid.pattern(orig_pattern, df_miss_valid)
  metrics <- matrix(ncol = length(k), nrow = length(alpha))

  j <- 1
  for(knn in k){
    i <- 1
    for(a in alpha){
      imp_val <- impute.multilevel.num.calibrate(df_miss_valid, distance, valid_pattern, nIter, a, knn)
      RMSE <- get.RMSE(df_miss, imp_val, valid_pattern)
      metrics[i, j] <- RMSE
      i <- i + 1
      if(print){
        print(paste0(a, ",", knn, ",", RMSE))
      }
    }
    j <- j + 1
  }
  best_ind <- which(metrics == min(metrics), arr.ind = TRUE)
  alpha_k <- c(alpha[best_ind[1]], k[best_ind[2]])
}


get.RMSE <- function(x, y, indx){
  RMSE_num <- 0
  vars_num <- 0
  nvar <- ncol(x)

  for(d in 2:(nvar)){
    for(i in indx[,d]){
      if(i != 0){
        RMSE_num <- RMSE_num + (x[i,d] - y[i,d])^2
        vars_num <- vars_num + 1
      }
    }
  }

  RMSE <- sqrt(RMSE_num/vars_num) #+ RMSE_cat/vars_cat)
  RMSE
}



#' Normalize dataset
#'
#' This function returns a dataset with normalized values for numerical variables
#'
#' @param data A dataframe with missing values
#' @return A dataframe with normalized values for the numerical variables
#' @export
#'
#' @examples
#' data(data_example)
#' normalize(data_example)
normalize <- function(df_miss){
  col <- ncol(df_miss)
  for(i in 2:col){
    if(is.numeric(df_miss[,i])){
      df_miss[,i] <- (df_miss[ ,i] - mean(df_miss[ ,i], na.rm = T))/sd(df_miss[ ,i], na.rm = T)
    }
  }
  df_miss
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


make.missing <- function(data, pattern){
  vars <- ncol(data)
  for(i in 1:vars){
    data[pattern[ , i] == 1, i] <- NA
  }
  data
}


#' @importFrom stats rbinom
make.pattern <- function(data, p = 0.1){
  vars <- ncol(data) - 1
  prop <- p/vars
  pattern <- rep(0, nrow(data))
  for(i in 1:vars){
    pattern <- cbind(pattern, rbinom(nrow(data), 1, prop))
  }
  pattern
}


get.valid.pattern <- function(orig_pattern, df_miss_valid){
  n <- nrow(df_miss_valid)
  p <- ncol(df_miss_valid)
  index_miss_valid <- matrix(NA, ncol = p, nrow = n)
  index_miss_joint <- get.index.miss(df_miss_valid)

  for(j in 1:p){
    for(i in 1:n){
      orig <- orig_pattern[i,j]
      joint <- index_miss_joint[i,j]
      if(orig == 0 & joint !=0){
        index_miss_valid[i,j] <- joint
      } else {
        index_miss_valid[i,j] <- 0
      }
    }
  }
  index_miss_valid
}

impute.multilevel.num.calibrate <- function(df_miss, distance, pattern_val, nIter, alpha, k){

  vars <- ncol(df_miss)
  rows <- nrow(df_miss)

  neighbors <- matrix(nrow = rows, ncol = k)
  index_miss <- pattern_val

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
}



normalize <- function(df_miss){
  col <- ncol(df_miss)
  for(i in 2:col){
    if(is.numeric(df_miss[,i])){
      df_miss[,i] <- (df_miss[ ,i] - mean(df_miss[ ,i], na.rm = T))/sd(df_miss[ ,i], na.rm = T)
    }
  }
  df_miss
}

