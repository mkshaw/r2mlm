#' Pull and prepare data.
#'
#' First, pull dataframe associated with model using broom::augment. Then add
#' interaction terms and other user-created variables to the dataframe.
#'
#' @param model A model generated using \code{\link[lme4]{lmer}} or
#'   \code{\link[nlme]{nlme}}, passed from the calling function.
#' @param calling_function Whether the helper function is r2mlm_lme4 or
#'   r2mlm_nlme.
#' @param cluster_variable Clustering variable in dataframe.
#'
#' @importFrom broom augment

prepare_data <- function(model, calling_function, cluster_variable) {

  # Step 1a: pull dataframe associated with model
  data <- broomExtra::augment(model)

  # Step 1b: remove the dot variables (.fitted, .resid, .fixed, etc.) from data.
  # lme4 returns way more dot variables
  if (calling_function == "nlme") {
    data <- data[1:(length(data) - 3)]
  } else if (calling_function == "lme4") {
    data <- data[1:(length(data) - 11)]
  }


  # Step 2: add interaction terms to the dataframe

  # * Step 2a) pull interaction terms into list
  interaction_vars <- get_interaction_vars(model)

  # * Step 2b) split interaction terms into halves, multiply halves to create new columns in dataframe

  for (whole in interaction_vars) {

    half1 <- str_split_fixed(whole, ":", 2)[1]
    half2 <- str_split_fixed(whole, ":", 2)[2]

    newcol <- dplyr::pull(data[half1] * data[half2])

    data <- data %>%
      dplyr::mutate(!!whole := newcol)

  }

  # Step 3: group data by cluster variable
  data <- data %>%
    dplyr::group_by(.data[[cluster_variable]])

  return(data)

}

