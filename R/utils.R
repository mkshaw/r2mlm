#' Helper functions.

#' get_covs: Return the number in the dataframe associated with a given variable
#' name.
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

#' get_cwc: Determine whether l1_vars are centered within cluster.
#'
#' @param l1_vars List of level 1 variables.
#' @param cluster_variable Clustering variable.
#' @param data Dataset.

get_cwc <- function(l1_vars, cluster_variable, data) {

  # Group data
  # see "Indirection" here for explanation of this group_by formatting: https://dplyr.tidyverse.org/articles/programming.html
  data_grouped <- data %>%
    dplyr::group_by(.data[[cluster_variable]])

  for (variable in l1_vars) {

    # for each group for the given variable, sum all values
    t <- data_grouped %>%
      dplyr::select(cluster_variable, variable) %>%
      dplyr::group_map(~ sum(.))

    # establish temporary tracker
    temp_tracker <- 0

    # sum all of the sums
    for (i in t) {
      temp_tracker <- temp_tracker + i
    }

    # if the biggie sum is essentially zero (not exactly zero, because floating point), then the variable is CWC
    if (temp_tracker < 0.0000001) {
      centeredwithincluster <- TRUE
    } else {
      centeredwithincluster <- FALSE
      break # break if even one variable is not CWC, because the r2mlm_manual function will need to center everything anyways
    }
  }

  return(centeredwithincluster)

}
