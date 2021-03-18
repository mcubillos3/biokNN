#' Plot number of missing values by class
#'
#' This function returns a dataframe with a multilevel structure. It generates a dataframe using a varying
#' intercepts/varying slopes linear regression with a single target variable y.
#' @importFrom magrittr %>%
#' @param df dataframe with missing values
#' @param className name of the variable containing classes
#' @return A barplot with the number of missing values by class, by variable
#' @export
plot.missing <- function(df, class){

  missClust <- df %>% dplyr::group_by(class) %>%
    dplyr::summarise_if(is.numeric, ~sum(is.na(.x))) %>%
    tidyr::gather(var, value, -class)

  g <- ggplot2::ggplot(missClust) +
    geom_bar(aes(x = class, y = value), stat="identity") +
    facet_wrap(~var)+
    labs(x = "Class", y ="Number of missing values") +
    theme_bw()
  g
}
