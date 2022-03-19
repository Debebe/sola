# Square runs square root of a function
#
# This is an example function named 'square'
# which prints 'Hello, world!'.
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'
#' Title - This function squares a number , any number
#' @param x inputs to the function
#' @return out output from the model
#' @export
square <- function(x) {
  out  <- x^2
  out
}

#' Title
#' @param p proportion
#' @return logit of a proportion
#' @export
logit <- function(p) {
  log(p / (1 - p))}


#' Title gives 8 times table
#'
#' @param x input to the function
#' @return the product of 8 and another number
#' @export eight_times_table
eight_times_table <- function(x){
out <-8*x
out
}
