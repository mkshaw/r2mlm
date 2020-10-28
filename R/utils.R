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

#' get_interaction_vars: Return list of interaction variables
#'
#' @param model A model generated using \code{\link[lme4]{lmer}} or
#'   \code{\link[nlme]{nlme}}, passed from the calling function.

get_interaction_vars <- function(model) {

  interaction_vars <- c()
  x <- 1
  for (term in attr(terms(model), "term.labels")) {

    if (grepl(":", term) == TRUE) {
      interaction_vars[x] <- term
      x <- x + 1
    }

  }

  return(interaction_vars)

}

#' sort_variables: Sort predictors into level 1 (variance is non-zero) and level
#' 2 (variance of 0) variables.
#'
#' @param

sort_variables <- function(data, predictors, cluster_variable) {

  l1_vars <- c()
  l2_vars <- c()

  # * Step 4c) Sort variables into l1_vars and l2_vars

  # set counters
  l1_counter <- 1
  l2_counter <- 1

  # * Step 4d) loop through all variables in grouped dataset

  for (variable in predictors) {

    # calculate variance for each cluster
    t <- data %>%
      dplyr::select(tidyselect::all_of(cluster_variable), variable) %>%
      dplyr::group_map(~ var(.))

    # var returns NA if group only 1 row. Replace with 0
    counter = 1

    while (counter < length(t)) {

      if (is.na(t[[counter]])) {
        t[[counter]] <- 0
      }

      counter <- counter + 1

    }

    # variable to track variance
    variance_tracker <- 0

    # add up the variance from each cluster
    for (i in t) {
      variance_tracker <- variance_tracker + i
    }

    # if each cluster has 0 variance, the sum of variance is 0, so it's an L2 variable
    if (variance_tracker == 0) {
      l2_vars[l2_counter] <- variable
      l2_counter <- l2_counter + 1
    } else {
      l1_vars[l1_counter] <- variable
      l1_counter <- l1_counter + 1
    }

  }

  return(list("l1_vars" = l1_vars, "l2_vars" = l2_vars))

}
