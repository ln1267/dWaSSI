# numberOfDays!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' A numberOfDays Function
#'
#' This function allows you to get the number of days for a specific month.
#' @param date A date object.
#' @keywords cats
#' @export
#' @examples
#' date<-as.Date("2001-01-01")
#' numberOfDays(date)



# get the number of days for each month
numberOfDays <- function(date) {
  m <- format(date, format="%m")

  while (format(date, format="%m") == m) {
    date <- date + 1
  }

  return(as.integer(format(date - 1, format="%d")))
}


