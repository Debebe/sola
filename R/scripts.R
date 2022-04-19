# Square runs square root of a function
#
# This is an example function named 'square'
# which prints 'Hello, world!'.
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'
#'  Square a number

#' @param x A number
#'
#' @return The square of \code{x}.
#' @export
#' @examples
#' square(4)
square <- function(x) {
  out  <- x^2
  out
}

#' Logits a number
#'
#' @param p proportion
#'
#' @return The logit of \code{p}.
#' @export
#' @examples
#' logit(0.4)
#'
logit <- function(p) {
  log(p / (1 - p))}


#'  8 times table
#'
#' @param x An input to the function
#'
#' @return The eight times  \code{x}.
#' @export eight_times_table
#'
eight_times_table <- function(x){
out <-8*x
out
}



#' Converts 95% confidence interval to lognormal distribution
#'
#' @param lowerCI  A lower bound of 95% confidence interval
#' @param upperCI  Upper bound of 95% confidence interval
#' @param mu       Mean of confidence interval
#'
#' @return Mean and sd of lognormal distribution
#' @export
#'
#' @examples
#' logNormal(0.02, 0.013, 0.035)

logNormal <- function (mu, lowerCI, upperCI) {
  sd <- (upperCI-lowerCI)/3.92
  # logNorm

  mean_lognorm  <- log(mu^2 / sqrt(sd^2 + mu^2))
  sd_lognorm    <- sqrt(log(1 + (sd^2 / mu^2)))
  return(params = list(mean_lognorm=mean_lognorm ,sd_lognorm=sd_lognorm))
}


#' Generates beta prior from 95% Confidence interval
#'
#' @param mu      Mean of confidence interval
#' @param lower  Upper bound of 95% confidence interval
#' @param upper  Upper bound of 95% confidence interval
#'
#' @return Renders parameters of beta prior distribution
#' @export
#'
#' @examples
#' BetaParams(mu=0.5, lower=0.35, upper=0.72)
BetaParams <- function(mu, lower, upper) {
  sd<- (upper-lower)/3.92
  alpha <- ((1 - mu) / sd^2 - 1 / mu) * mu ^ 2
  beta  <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}


#' Function to conduct logistic population interpolation at time points with missing data.
#'
#' @param population_age_sex data frame containing space, age, sex, year stratified data
#'
#' @return  data frame with imputed (interpolated data)
#' @export
#'
#' @examples
#' log_linear_pop_agesex_interpolation(naomi::demo_population_agesex)

log_linear_pop_agesex_interpolation <- function(population_age_sex){

    if(!"population" %in% colnames(population_age_sex)){
      stop(paste0("population is missing… consider changing variable names in the data"))
    }
    pop_imputed <- population_age_sex %>%
      arrange(area_id, area_name, area_level,sex, age_group, year) %>%
      group_by(area_id, area_name, area_level,sex, age_group) %>%
      # logistic interpolation
      mutate(art_current_interpolated = exp(zoo::na.approx(log(population), x = year,
                                                           rule = 2, na.rm = FALSE)))%>%
      # linear interpolation
      # mutate(art_current_interpolated = zoo::na.approx(population, x = year,
      #                                                      rule = 2, na.rm = FALSE)))%>%
      ungroup()

     return(pop_imputed)

  }

  #pop_interpolated <- interpolate_population_agesex(tmp)

#' Convert list into a data.frame
#'
#' @param list A list
#'
#' @return A dataframe.
#' @export
#'
#' @examples
list_to_df <- function(list){
  data.frame(dplyr::bind_rows(list, .id = "replicate"))
}

#' Title
#'
#' @param data vector that records data inputs from all areal units
#' @param sf_object shape file in sf class
#'
#' @return ggplot of clusters
#' @export
#'
#' @examples
moran_clusters <- function(data, sf_object){

  neighbours_sf <- spdep::poly2nb(sf_object)   # adjacency matrix from sf
  #listw <- nb2listw(neighbours_sf)      # weights matrix

  local <- spdep::localmoran(x = data, listw = nb2listw(neighbours_sf, style = "W"))
  clusters <-as.data.frame(attr(local,"quadr"))

  moran.map <- cbind(sf_object, local,clusters )%>%
    dplyr::rename(cluster=mean)

  gg <-ggplot2::ggplot(data = moran.map) +
    sf::geom_sf(aes(fill = cluster)) +
    scale_fill_viridis_d(option = "plasma")+theme_minimal()+
    theme(axis.ticks = element_blank(),
          axis.text = element_blank())

  gg

}

#moran_clusters(data$chNotifRate, NorthShowa_sf)

#' Title
#'
#' @param data a numeric vector the same length as the neighbours list in sf object
#' @param sf_object spatial data as sf class
#'
#' @return Global moran value with p-value
#' @export
#'
#' @examples
global_moranI <- function(data, sf_object){
    globalMoran <- spdep::moran.test(data, listw)
    global_moran_est <- globalMoran[["estimate"]][["Moran I statistic"]]
    global_moran_pval <- globalMoran[["p.value"]]

return(list(global_moran_estimate=global_moran_est,
            global_moran_pval= global_moran_pval))
}

#global_moranI(data$chNotifRate, NorthShowa)
