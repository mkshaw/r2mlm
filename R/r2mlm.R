#' Compute R-squared values for multilevel models.
#'
#' \code{r2mlm} reads in a multilevel model (MLM) generated using
#' \code{\link[lme4]{lmer}} or \code{\link[nlme]{nlme}}, and outputs all
#' relevant R-square measures and barchart decompositions. Any number of level-1
#' and/or level-2 predictors is supported. Any of the level-1 predictors can
#' have random slopes.
#'
#' \code{r2mlm} first determines whether a given model was generated using
#' \code{\link[lme4]{lmer}} or \code{\link[nlme]{nlme}}, then passes the model
#' to helper functions that pull the raw data and parameter estimates from the
#' model, and pass that information to \code{\link{r2mlm_manual}}.
#'
#' @param model A model generated using \code{\link[lme4]{lmer}} or
#'   \code{\link[nlme]{nlme}}.
#'
#' @return If input is a valid model, then the output will be a list and
#'   associated graphical representation of R-squared decompositions. If model
#'   is not valid, it will return an error prompting the user to input a valid
#'   model.
#'
#' @examples
#' \dontrun{
#' # Using lme4 for your model
#'
#' model_lme4 <- lmer(satisfaction ~ 1 + salary_c + control_c + salary_m +
#' control_m + s_t_ratio + (1 | schoolID), data = teachsat, REML =
#' TRUE)
#'
#' r2mlm(model_lme4)
#'
#' # Using nlme for your model
#'
#' model_nlme <- lme(satisfaction ~ 1 + salary_c + control_c + salary_m +
#'                   control_m + s_t_ratio,
#'                   random = ~ 1 | schoolID,
#'                   data = teachsat,
#'                   method = "REML")
#'
#' r2mlm(model_nlme)
#' }
#'
#' @seealso \href{https://doi.org/10.1037/met0000184}{Rights, J. D., & Sterba,
#'   S. K. (2019). Quantifying explained variance in multilevel models: An
#'   integrative framework for defining R-squared measures. Psychological
#'   Methods, 24(3), 309–338.}
#'
#' @family r2mlm single model functions
#'
#' @importFrom lme4 fortify.merMod ranef fixef VarCorr getME
#' @importFrom nlme asOneFormula
#' @importFrom magrittr %>%
#' @importFrom stats terms formula model.frame
#' @importFrom stringr str_split_fixed
#' @importFrom rlang := .data
#'
#'
#' @export

# 1 r2mlm Master Function -------------------------------------------------

r2mlm <- function(model) {

  if (typeof(model) == "list") {
    r2mlm_nlme(model)
  } else if (typeof(model) == "S4") {
    r2mlm_lmer(model)
  } else {
    stop("You must input a model generated using either the lme4 or nlme package.")
  }

}

# 2 r2mlm_lmer helper function --------------------------------------------

r2mlm_lmer <- function(model) {

  # Step 0: throw error if model contains higher-order terms
  temp_formula <- formula(model)
  grepl_array <- grepl("I(", temp_formula, fixed = TRUE)

  for (bool in grepl_array) {
    if (bool == TRUE) {
      stop("Higher-order terms created with I() syntax are not currently accepted. To include a higher-order term, you must compute it in advance and include it as a column in your dataset.")
    }
  }

  # Step 1: pull data

  data <- fortify.merMod(model)

  # Step 2: has_intercept

  if ((terms(model) %>% attr("intercept")) == 1) {
    has_intercept = TRUE
  } else {
    has_intercept = FALSE
  }

  # Step 3: Pull l1 var names into a variable called l1_vars, and vice versa for l2 var names

  # i) Pull all variable names from the formula
  all_vars <- all.vars(formula(model))
  formula_length <- length(all_vars) # this returns the number of elements in the all_vars list
  cluster_variable <- all_vars[formula_length] # pull cluster, we'll need it later

  # ii) determine whether data is appropriate format. Only the cluster variable can be a factor, for now
  # a) Pull all variables except for cluster

  all_vars_except_cluster <- all_vars[1:(length(all_vars) - 1)]

  # b) If any of those variables is non-numeric, then throw an error

  for (var in all_vars_except_cluster) {

    if (!(class(data[[var]]) == "integer") && !(class(data[[var]]) == "numeric")) {
      stop("Your data must be numeric. Only the cluster variable can be a factor.")
    }

  }

  # iii) Create vectors to fill
  l1_vars <- c()
  l2_vars <- c()

  # iv) Sort variables into l1_vars and l2_vars

  # (a) Pull your temp dataset to work with

  temp_data <- model.frame(model) # this pulls a df of all variables and values used in the model

  # (b) group dataset by clustering variable
  temp_data_grouped <- temp_data %>%
    dplyr::group_by_at(formula_length) #group by the clustering variable, which is the last variable in the df (by virtue of how formula(model) works, where it pulls out the formula, and the last variable is on the other side of the |, i.e., the clustering variable)

  # (c) define variables you'll need

  # all variables to sort into L1 and L2
  all_vars_except_cluster_and_outcome <- all_vars_except_cluster[-1]

  # set counters
  l1_counter <- 1
  l2_counter <- 1

  # (d) loop through all variables in grouped dataset

  for (variable in all_vars_except_cluster_and_outcome) {

    # calculate variance for each cluster
    t <- temp_data_grouped %>%
      dplyr::select(tidyselect::all_of(cluster_variable), variable) %>%
      dplyr::group_map(~ var(.))

    # variable to track variance
    variance_tracker <- 0

    # add up the variance from each cluster
    for (i in t) {
      variance_tracker <- variance_tracker + i
    }

    # if the sum of variance is 0, then each cluster has 0 variance, so it's an L2 variable
    if (variance_tracker == 0) {
      l2_vars[l2_counter] <- variable
      l2_counter <- l2_counter + 1
    } else {
      l1_vars[l1_counter] <- variable
      l1_counter <- l1_counter + 1
    }

  }

  # Step 4: pull variable names for L1 predictors with random slopes into a variable called random_slope_vars

  temp_cov_list <- ranef(model)[[1]]

  # determine where to start pulling from the temp_cov_list, depending on whether you need to bypass the `(Intercept)` column
  if (has_intercept == 1) {
    running_count <- 2
  } else {
    running_count <- 1
  }

  random_slope_vars <- c()
  x <- 1 # counter for indexing in random_slope_vars

  # running count less than or equal to list length, so it traverses the entire list (doesn't leave last element off)
  while (running_count <= length(temp_cov_list)) {
    random_slope_vars[x] <- names(temp_cov_list[running_count])
    x <- x + 1
    running_count <- running_count + 1
  }

  interaction_vars <- c()
  x <- 1
  for (term in attr(terms(model), "term.labels")) {

    if (grepl(":", term) == TRUE) {
      interaction_vars[x] <- term
      x <- x + 1
    }

  }

  data <- data %>% dplyr::ungroup() # need to ungroup in order to create new columns

  for (whole in interaction_vars) {

    half1 <- str_split_fixed(whole, ":", 2)[1]
    half2 <- str_split_fixed(whole, ":", 2)[2]

    if (!is.na(match(half1, l2_vars)) && !is.na(match(half2, l2_vars))) {
      l2_vars[l2_counter] <- whole
      l2_counter <- l2_counter + 1
    } else {
      l1_vars[l1_counter] <- whole
      l1_counter <- l1_counter + 1
    }

    newcol <- dplyr::pull(data[half1] * data[half2])

    data <- data %>%
      dplyr::mutate(!!whole := newcol)

  }

  # Step 5: determine value of centeredwithincluster

  # (a) group data

  data_grouped <- data %>%
    dplyr::group_by(.data[[cluster_variable]]) # see "Indirection" here for explanation of this group_by formatting: https://dplyr.tidyverse.org/articles/programming.html

  if (is.null(l1_vars)) {
    centeredwithincluster <- TRUE
  } else {
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
  }

  # Step 6: pull column numbers for _covs variables
  # 6a) within_covs (l1 variables)
  # for (each value in l1_vars list) {match(value, names(data))}
  within <- c()
  i = 0
  for (var in l1_vars) {
    i = i + 1
    tmp <- match(var, names(data))
    within[i] <- tmp
  }

  # 6b) pull column numbers for between_covs (l2 variables)
  between <- c()
  i = 1
  for (var in l2_vars) {
    tmp <- match(var, names(data))
    between[i] <- tmp
    i = i + 1
  }

  # 6c) pull column numbers for random_covs (l1 variables with random slopes)
  random <- c()
  i = 1
  for (var in random_slope_vars) {
    tmp <- match(var, names(data))
    random[i] <- tmp
    i = i + 1
  }

  # Step 7: pull gamma values (fixed slopes)
  # 7a) gamma_w, fixed slopes for L1 variables (from l1_vars list)
  gammaw <- c()
  i = 1
  for (var in l1_vars) {
    gammaw[i] <- fixef(model)[var]
    i = i + 1
  }

  # 7b) gamma_b, intercept value if hasintercept = TRUE, and fixed slopes for L2 variables (from between list)
  gammab <- c()
  if (has_intercept == TRUE) {
    gammab[1] <- fixef(model)[1]
    i = 2
  } else {
    i = 1
  }
  for (var in l2_vars) {
    gammab[i] <- fixef(model)[var]
    i = i + 1
  }

  # Step 8: Tau matrix, results from VarCorr(model)

  vcov <- VarCorr(model)
  tau <- as.matrix(Matrix::bdiag(vcov))

  # Step 9: sigma^2 value, Rij

  sigma2 <- getME(model, "sigma")^2

  # Step 10: input everything into r2MLM

  r2mlm_manual(as.data.frame(data), within_covs = within, between_covs = between, random_covs = random, gamma_w = gammaw, gamma_b = gammab, Tau = tau, sigma2 = sigma2, has_intercept = has_intercept, clustermeancentered = centeredwithincluster)

}

# 3 r2mlm_nlme helper function --------------------------------------------

r2mlm_nlme <- function(model) {

  # Step 0: throw error if model contains higher-order terms
  temp_formula <- formula(model)
  grepl_array <- grepl("I(", temp_formula, fixed = TRUE)

  for (bool in grepl_array) {
    if (bool == TRUE) {
      stop("Higher-order terms created with I() syntax are not currently accepted. To include a higher-order term, you must compute it in advance and include it as a column in your dataset.")
    }
  }

  # Step 1: pull data

  data <- model$data

  # Step 2: has_intercept

  if ((terms(model) %>% attr("intercept")) == 1) {
    has_intercept = TRUE
  } else {
    has_intercept = FALSE
  }

  # Step 3: Pull l1 var names into a variable called l1_vars, and vice versa for l2 var names

  # i) Pull all variable names from the formula
  all_vars <- all.vars(formula(model))

  # in nlme, formula(model) doesn't return grouping var, but we'll need that later on, so let's grab it here
  full_formula <- all.vars(asOneFormula(model))
  cluster_variable <- full_formula[length(full_formula)]

  # Then add the grouping var to list of all variables, and calculate formula length (for later use, to iterate)
  all_vars[length(all_vars) + 1] <- cluster_variable
  formula_length <- length(all_vars) # this returns the number of elements in the all_vars list

  # ii) determine whether data is appropriate format. Only the cluster variable can be a factor, for now

  # b) If any of those variables is non-numeric, then throw an error

  all_vars_except_cluster <- all.vars(formula(model))

  for (var in all_vars_except_cluster) {

    if (!(class(data[[var]]) == "integer") && !(class(data[[var]]) == "numeric")) {
      stop("Your data must be numeric. Only the cluster variable can be a factor.")
    }

  }

  # iii) Create vectors to fill
  l1_vars <- c()
  l2_vars <- c()

  # iv) Sort variables into l1_vars and l2_vars

  # (a) Pull your temp dataset to work with

  temp_data <- data %>%
    dplyr::select(tidyselect::all_of(all_vars)) # this pulls a df of all variables and values used in the model

  # (b) group dataset by clustering variable
  temp_data_grouped <- temp_data %>%
    dplyr::group_by_at(formula_length) #group by the clustering variable, which is the last variable in the df (by virtue of how formula(model) works, where it pulls out the formula, and the last variable is on the other side of the |, i.e., the clustering variable)

  # (c) define variables you'll need

  # all variables to sort into L1 and L2
  all_vars_except_cluster_and_outcome <- all_vars_except_cluster[-1]

  # set counters
  l1_counter <- 1
  l2_counter <- 1

  # (d) loop through all variables in grouped dataset

  for (variable in all_vars_except_cluster_and_outcome) {

    # calculate variance for each cluster
    t <- temp_data_grouped %>%
      dplyr::select(tidyselect::all_of(cluster_variable), variable) %>%
      dplyr::group_map(~ var(.))

    # variable to track variance
    variance_tracker <- 0

    # add up the variance from each cluster
    for (i in t) {
      variance_tracker <- variance_tracker + i
    }

    # if the sum of variance is 0, then each cluster has 0 variance, so it's an L2 variable
    if (variance_tracker == 0) {
      l2_vars[l2_counter] <- variable
      l2_counter <- l2_counter + 1
    } else {
      l1_vars[l1_counter] <- variable
      l1_counter <- l1_counter + 1
    }

  }

  # Step 4: pull variable names for L1 predictors with random slopes into a variable called random_slope_vars

  temp_cov_list <- nlme::ranef(model)

  # determine where to start pulling from the temp_cov_list, depending on whether you need to bypass the `(Intercept)` column
  if (has_intercept == 1) {
    running_count <- 2
  } else {
    running_count <- 1
  }

  random_slope_vars <- c()
  x <- 1 # counter for indexing in random_slope_vars

  # running count less than or equal to list length, so it traverses the entire list (doesn't leave last element off)
  while (running_count <= length(temp_cov_list)) {
    random_slope_vars[x] <- names(temp_cov_list[running_count])
    x <- x + 1
    running_count <- running_count + 1
  }

  interaction_vars <- c()
  x <- 1
  for (term in attr(terms(model), "term.labels")) {

    if (grepl(":", term) == TRUE) {
      interaction_vars[x] <- term
      x <- x + 1
    }

  }

  data <- data %>% dplyr::ungroup() # need to ungroup in order to create new columns

  for (whole in interaction_vars) {

    half1 <- str_split_fixed(whole, ":", 2)[1]
    half2 <- str_split_fixed(whole, ":", 2)[2]

    if (!is.na(match(half1, l2_vars)) && !is.na(match(half2, l2_vars))) {
      l2_vars[l2_counter] <- whole
      l2_counter <- l2_counter + 1
    } else {
      l1_vars[l1_counter] <- whole
      l1_counter <- l1_counter + 1
    }

    newcol <- dplyr::pull(data[half1] * data[half2])

    data <- data %>%
      dplyr::mutate(!!whole := newcol)

  }

  # Step 5: determine value of centeredwithincluster

  # (a) group data

  data_grouped <- data %>%
    dplyr::group_by(.data[[cluster_variable]]) # annoyingly written, because group_by(!!cluster_variable)) doesn't work

  if (is.null(l1_vars)) {
    centeredwithincluster <- TRUE
  } else {
    for (variable in l1_vars) {

      # for each group for the given variable, sum all values
      t <- data_grouped %>%
        dplyr::select(cluster_variable, variable) %>% # select cluster_variable and variable (the former to prevent "Adding missing grouping variables" printout)
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
  }

  # Step 6: pull column numbers for _covs variables
  # 6a) within_covs (l1 variables)
  # for (each value in l1_vars list) {match(value, names(data))}
  within <- c()
  i = 0
  for (var in l1_vars) {
    i = i + 1
    tmp <- match(var, names(data))
    within[i] <- tmp
  }

  # 6b) pull column numbers for between_covs (l2 variables)
  between <- c()
  i = 1
  for (var in l2_vars) {
    tmp <- match(var, names(data))
    between[i] <- tmp
    i = i + 1
  }

  # 6c) pull column numbers for random_covs (l1 variables with random slopes)
  random <- c()
  i = 1
  for (var in random_slope_vars) {
    tmp <- match(var, names(data))
    random[i] <- tmp
    i = i + 1
  }

  # Step 7: pull gamma values (fixed slopes)
  # 7a) gamma_w, fixed slopes for L1 variables (from l1_vars list)
  gammaw <- c()
  i = 1
  for (var in l1_vars) {
    gammaw[i] <- nlme::fixef(model)[var]
    i = i + 1
  }

  # 7b) gamma_b, intercept value if hasintercept = TRUE, and fixed slopes for L2 variables (from between list)
  gammab <- c()
  if (has_intercept == TRUE) {
    gammab[1] <- nlme::fixef(model)[1]
    i = 2
  } else {
    i = 1
  }
  for (var in l2_vars) {
    gammab[i] <- nlme::fixef(model)[var]
    i = i + 1
  }

  # Step 8: Tau matrix

  tau <- nlme::getVarCov(model)

  # Step 9: sigma^2 value, Rij

  sigma2 <- model$sigma^2

  # Step 10: input everything into r2mlm

  r2mlm_manual(as.data.frame(data), within_covs = within, between_covs = between, random_covs = random, gamma_w = gammaw, gamma_b = gammab, Tau = tau, sigma2 = sigma2, has_intercept = has_intercept, clustermeancentered = centeredwithincluster)

}


