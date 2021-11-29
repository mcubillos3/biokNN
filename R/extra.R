#' Plot number of missing values by class
#'
#' This function returns a dataframe with a multilevel structure. It generates a dataframe using a varying
#' intercepts/varying slopes linear regression with a single target variable y.
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot aes geom_bar coord_flip facet_wrap labs theme_bw
#' @param df dataframe with missing values
#' @param className name of the variable containing classes
#' @param varNames vector with the name of the variables to be imputed
#' @return A barplot with the number of missing values by class, by variable
#' @export
#'
#' @examples
#' data(data_example)
#' missing_plot(data_example, "class", c("y"))
missing_plot <- function(df, className, varNames){

  variables <- c(className, varNames)
  colnames(df)[which(names(df) == className)] <- "class"

  missClust <- df %>%
    dplyr::select(one_of(variables)) %>%
    dplyr::group_by(class) %>%
    dplyr::summarise_if(is.numeric, ~sum(is.na(.x))) %>%
    tidyr::gather(var, value, -class)

  g <- ggplot(missClust) +
    geom_bar(aes(x = class, y = value), stat="identity") +
    facet_wrap(~var)+
    labs(x = "Class", y ="Number of missing values") +
    theme_bw()
  g
}

#' Plot pattern of missing values by class
#'
#' This function returns a dataframe with a multilevel structure. It generates a dataframe using a varying
#' intercepts/varying slopes linear regression with a single target variable y.
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot aes coord_flip facet_wrap labs theme_bw scale_fill_continuous geom_tile
#' @param df dataframe with missing values
#' @param className name of the variable containing classes
#' @param varNames vector with the name of the variables to be imputed
#' @return A plot with the patter of missing values by class, by variable
#' @export
#'
#' @examples
#' data(data_example)
#' pattern_plot(data_example, "class")
pattern_plot <- function(df, className, varNames){

  variables <- c(className, varNames)
  colnames(df)[which(names(df) == className)] <- "class"

  df_na <- df %>%
    dplyr::select(one_of(variables)) %>%
    dplyr::group_by(class) %>%
    dplyr::mutate(obs = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    tidyr::gather(var, value, -c(class, obs))

  g <- ggplot(df_na, aes(class, obs, fill= value)) +
    facet_wrap(~var, ncol = 1)+
    geom_tile(colour = "black") +
    scale_fill_continuous(high = "gray", na.value = 'white') +
    theme_bw() +
    labs( y ="")
  g

}

#' Plot pattern of missing values by class
#'
#' This function returns a dataframe with a multilevel structure. It generates a dataframe using a varying
#' intercepts/varying slopes linear regression with a single target variable y.
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot aes coord_flip facet_wrap labs theme_bw scale_fill_continuous geom_tile
#' @param df dataframe with missing values
#' @param y target variable
#' @param class name of the variable containing classes
#' @return A boxplot for each class of the target variable
#' @importFrom ggplot2 ggplot aes geom_boxplot theme_bw

#' @export
#'
#' @examples
#' data(data_example)
#' target_boxplot(data_example, y, "class")
target_boxplot <- function(df, y, class){

  ggplot(df %>% tidyr::drop_na(), aes(class, y)) +
    geom_boxplot()+
    theme_bw()
}


#' Generate multilevel dataset
#'
#' This function returns a dataframe with a multilevel structure. It generates a dataframe using a varying
#' intercepts/varying slopes linear regression with a single target variable y.
#'
#' @param nClass number of classes
#' @param nVars number of independent variables (X)
#' @param classMean average number of observations per class
#' @param classSD standard deviation of the number of observations per class
#' @param beta0 intercept parameter
#' @param tau0 variance of the parameter between classes
#' @param beta vector with the slope parameters, one for each independent variable
#' @param tau vector with the variance of the slope parameters, one for each independent variable
#' @param sigma2 error variance
#' @importFrom stats rnorm
#' @importFrom magrittr %>%
#' @return A dataframe with the multilevel dataset
#' @export
#'
#' @examples
#' df <- create_multilevel(nClass = 20,
#'                           nVars = 1,
#'                           classMean = 10,
#'                           classSD = 2)
create_multilevel <- function(nClass = 10, nVars = 1, classMean = 10, classSD = 0, beta0 = 0, tau0 = 1, beta = c(1), tau = c(1), sigma2 = 1) {

  sizeCluster <- as.integer(rnorm(nClass, mean = classMean, sd = classSD))
  sizeCluster[sizeCluster<=5] <- 5
  class <- c()
  beta0_j <- c()
  n <- sum(sizeCluster)
  X <- matrix(nrow = n, ncol = nVars)


  i <- 1
  for(k in sizeCluster){
    class <- c(class, rep(i, k))
    beta0_j <- c(beta0_j, rep(rnorm(1, beta0, tau0), k))
    i <- i + 1
  }

  error <- rnorm(n, mean = 0, sd = sigma2)
  y <- beta0_j + error

  for(j in 1:nVars){
    beta_j <- c()
    for(k in sizeCluster){
      beta_j <- c(beta_j, rep(rnorm(1, beta[j], tau[j]), k))
    }
    x <- rnorm(n)
    X[ ,j] <- x
    y <- y + beta_j*X

  }

  df <- data.frame(class, y, X) %>% dplyr::mutate(class = as.factor(class))
  df

}







