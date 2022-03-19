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
