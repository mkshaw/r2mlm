prepare_data <- function(model, calling_function, cluster_variable, second_model = NULL) {

  # Step 1a: pull dataframe associated with model
  if (calling_function == "lme4") {
    data <- model@frame
  } else {
    data <- model[["data"]]
  }

  # Step 2: add interaction terms to the dataframe

  # * Step 2a) pull interaction terms into list
  interaction_vars <- get_interaction_vars(model)

  if(!is.null(second_model)) {
    interaction_vars_2 <- get_interaction_vars(second_model)
    interaction_vars <-  unique(append(interaction_vars, interaction_vars_2))
  }

  # * Step 2b) split interaction terms into halves, multiply halves to create new columns in dataframe

  data <- add_interaction_vars_to_data(data, interaction_vars)

  # Step 3: group data by cluster variable

  data <- group_data(data, cluster_variable)

  return(data)

}

add_interaction_vars_to_data <- function(data, interaction_vars) {

  for (whole in interaction_vars) {

    half1 <- str_split_fixed(whole, ":", 2)[1]
    half2 <- str_split_fixed(whole, ":", 2)[2]

    newcol <- dplyr::pull(data[half1] * data[half2])

    data <- data %>%
      dplyr::mutate(!!whole := newcol)

  }

  return(data)

}

group_data <- function(data, cluster_variable) {

  data <- data %>%
    dplyr::group_by(.data[[cluster_variable]])

  return(data)

}
