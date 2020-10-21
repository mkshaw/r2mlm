#' Helper functions.

#' Return the number in the dataframe associated with a given variable name
#'
#' @param variable_list List of character names of variables.
#' @param data Dataset containing those variable names.
#'

get_covs <- function(variable_list, data) {

  cov_list <- c() # create empty list to fill

  i = 1
  for (variable in variable_list) {
    tmp <- match(variable, names(data)) # get column number
    cov_list[i] <- tmp # add column number to list
    i = i + 1
  }

  return(cov_list)

}

