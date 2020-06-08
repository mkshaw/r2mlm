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

  # ii) Create vectors to fill
  l1_vars <- c()
  l2_vars <- c()

  # iii) Sort variables into l1_vars and l2_vars

  # (a) Pull your temp dataset to work with

  temp_data <- model.frame(model) # this pulls a df of all variables and values used in the model

  # (b) Isolate the largest group from the dataframe, which you'll use to test variances of variables to sort into l1 and l2 lists
  number <- temp_data %>%
    dplyr::group_by_at(formula_length) %>% #group by the clustering variable, which is the last variable in the df (by virtue of how formula(model) works, where it pulls out the formula, and the last variable is on the other side of the |, i.e., the clustering variable)
    dplyr::count() %>%
    dplyr::ungroup() %>% #have to ungroup because otherwise top_n will return n rows from each group, rather than n groups
    dplyr::top_n(1) %>% # returns the ID and N of the largest group
    dplyr::pull(1) # returns the number of the group that you'll use for your variance check

  # (c) Filter temp_data by the number you extracted in 3b
  temp_data_number <- temp_data %>%
    dplyr::filter(temp_data[formula_length] == number) #temp_data[formula_length] is the column that holds the clustering variable

  # (d) Iterate through temp_data_number, calculating the variance for each variable in all_vars, and then sorting by whether variance is 0 (l2) or non-zero (l1)
  x <- 2 # setting a counter overall, starting at 2 to skip the outcome variable (which is otherwise var1 in l1_vars)
  l1_counter <- 1 # setting a counter for adding to l1_vars list
  l2_counter <- 1 # setting a counter for adding to l2_vars list
  while (x < formula_length) {
    if (var(temp_data_number[x]) == 0) {
      l2_vars[l2_counter] <- names(temp_data_number[x])
      l2_counter <- l2_counter + 1
    } else {
      l1_vars[l1_counter] <- names(temp_data_number[x])
      l1_counter <- l1_counter + 1
    }
    x <- x + 1 # iterate the counter
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

  # Step 5: determine value of centeredwithincluster

  for(var in l1_vars) {

    # Sum the l1 column at hand (var in l1_vars)
    temp_sum <- temp_data_number %>%
      dplyr::summarize(
        sum = sum(temp_data_number[var])
      ) %>%
      dplyr::select(sum)

    # If that sum is approximately equal to zero (i.e., less than a very small number, to account for floating point issues),
    #   then the column is centered within cluster
    if (temp_sum < 0.0000001) {
      centeredwithincluster <- TRUE
    } else {
      centeredwithincluster <- FALSE # If the sum is non-zero, then the column is not CWC
      break # so break out of the for loop because if at least one L1 var is not CWC, then the variables will need to be centered by the r2MLM function
    }
  }

  # Step 5: pull column numbers for _covs variables
  # 5a) within_covs (l1 variables)
  # for (each value in l1_vars list) {match(value, names(data))}
  within <- c()
  i = 0
  for (var in l1_vars) {
    i = i + 1
    tmp <- match(var, names(data))
    within[[i]] <- tmp
  }

  # 5b) pull column numbers for between_covs (l2 variables)
  between <- c()
  i = 1
  for (var in l2_vars) {
    tmp <- match(var, names(data))
    between[[i]] <- tmp
    i = i + 1
  }

  # 5c) pull column numbers for random_covs (l1 variables with random slopes)
  random <- c()
  i = 1
  for (var in random_slope_vars) {
    tmp <- match(var, names(data))
    random[[i]] <- tmp
    i = i + 1
  }

  # Step 6: pull gamma values (fixed slopes)
  # 6a) gamma_w, fixed slopes for L1 variables (from l1_vars list)
  gammaw <- c()
  i = 1
  for (var in l1_vars) {
    gammaw[[i]] <- fixef(model)[var]
    i = i + 1
  }

  # 6b) gamma_b, intercept value if hasintercept = TRUE, and fixed slopes for L2 variables (from between list)
  gammab <- c()
  if (has_intercept == TRUE) {
    gammab[[1]] <- fixef(model)[1]
    i = 2
  } else {
    i = 1
  }
  for (var in l2_vars) {
    gammab[[i]] <- fixef(model)[var]
    i = i + 1
  }

  # Step 7: Tau matrix, results from VarCorr(model)

  vcov <- VarCorr(model)
  tau <- as.matrix(Matrix::bdiag(vcov))

  # Step 8: sigma^2 value, Rij

  sigma2 <- getME(model, "sigma")^2

  # Step 9: input everything into r2MLM

  r2mlm_manual(as.data.frame(data), within_covs = within, between_covs = between, random_covs = random, gamma_w = gammaw, gamma_b = gammab, Tau = tau, sigma2 = sigma2, has_intercept = has_intercept, clustermeancentered = centeredwithincluster)

}

# 3 r2mlm_nlme helper function --------------------------------------------

r2mlm_nlme <- function(model) {

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
  grouping_var <- full_formula[length(full_formula)]

  # Then add the grouping var to list of all variables, and calculate formula length (for later use, to iterate)
  all_vars[length(all_vars) + 1] <- grouping_var
  formula_length <- length(all_vars) # this returns the number of elements in the all_vars list

  # ii) Create vectors to fill
  l1_vars <- c()
  l2_vars <- c()

  # iii) Sort variables into l1_vars and l2_vars

  # (a) Pull your temp dataset to work with

  temp_data <- data %>%
    dplyr::select(tidyselect::all_of(all_vars)) # this pulls a df of all variables and values used in the model

  # (b) Isolate the largest group from the dataframe, which you'll use to test variances of variables to sort into l1 and l2 lists
  number <- temp_data %>%
    dplyr::group_by_at(formula_length) %>% #group by the clustering variable, which is the last variable in the df (by virtue of how formula(model) works, where it pulls out the formula, and the last variable is on the other side of the |, i.e., the clustering variable)
    dplyr::count() %>%
    dplyr::ungroup() %>% #have to ungroup because otherwise top_n will return n rows from each group, rather than n groups
    dplyr::top_n(1) %>% # returns the ID and N of the largest group
    dplyr::pull(1) # returns the number of the group that you'll use for your variance check

  # (c) Filter temp_data by the number you extracted in 3b
  temp_data_number <- temp_data %>%
    dplyr::ungroup() %>%
    dplyr::filter(temp_data[formula_length] == number) #temp_data[formula_length] is the column that holds the clustering variable

  # (d) Iterate through temp_data_number, calculating the variance for each variable in all_vars, and then sorting by whether variance is 0 (l2) or non-zero (l1)
  x <- 2 # setting a counter overall, starting at 2 to skip the outcome variable (which is otherwise var1 in l1_vars)
  l1_counter <- 1 # setting a counter for adding to l1_vars list
  l2_counter <- 1 # setting a counter for adding to l2_vars list
  while (x < formula_length) {
    if (var(temp_data_number[x]) == 0) {
      l2_vars[l2_counter] <- names(temp_data_number[x])
      l2_counter <- l2_counter + 1
    } else {
      l1_vars[l1_counter] <- names(temp_data_number[x])
      l1_counter <- l1_counter + 1
    }
    x <- x + 1 # iterate the counter
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

  # Step 5: determine value of centeredwithincluster

  for(var in l1_vars) {

    # Sum the l1 column at hand (var in l1_vars)
    temp_sum <- temp_data_number %>%
      dplyr::summarize(
        sum = sum(temp_data_number[var])
      ) %>%
      dplyr::select(sum)

    # If that sum is approximately equal to zero (i.e., less than a very small number, to account for floating point issues),
    #   then the column is centered within cluster
    if (temp_sum < 0.0000001) {
      centeredwithincluster <- TRUE
    } else {
      centeredwithincluster <- FALSE # If the sum is non-zero, then the column is not CWC
      break # so break out of the for loop because if at least one L1 var is not CWC, then the variables will need to be centered by the r2mlm function
    }
  }

  # Step 5: pull column numbers for _covs variables
  # 5a) within_covs (l1 variables)
  # for (each value in l1_vars list) {match(value, names(data))}
  within <- c()
  i = 0
  for (var in l1_vars) {
    i = i + 1
    tmp <- match(var, names(data))
    within[[i]] <- tmp
  }

  # 5b) pull column numbers for between_covs (l2 variables)
  between <- c()
  i = 1
  for (var in l2_vars) {
    tmp <- match(var, names(data))
    between[[i]] <- tmp
    i = i + 1
  }

  # 5c) pull column numbers for random_covs (l1 variables with random slopes)
  random <- c()
  i = 1
  for (var in random_slope_vars) {
    tmp <- match(var, names(data))
    random[[i]] <- tmp
    i = i + 1
  }

  # Step 6: pull gamma values (fixed slopes)
  # 6a) gamma_w, fixed slopes for L1 variables (from l1_vars list)
  gammaw <- c()
  i = 1
  for (var in l1_vars) {
    gammaw[[i]] <- nlme::fixef(model)[var]
    i = i + 1
  }

  # 6b) gamma_b, intercept value if hasintercept = TRUE, and fixed slopes for L2 variables (from between list)
  gammab <- c()
  if (has_intercept == TRUE) {
    gammab[[1]] <- nlme::fixef(model)[1]
    i = 2
  } else {
    i = 1
  }
  for (var in l2_vars) {
    gammab[[i]] <- nlme::fixef(model)[var]
    i = i + 1
  }

  # Step 7: Tau matrix

  tau <- nlme::getVarCov(model)

  # Step 8: sigma^2 value, Rij

  sigma2 <- model$sigma

  # Step 9: input everything into r2mlm

  r2mlm_manual(as.data.frame(data), within_covs = within, between_covs = between, random_covs = random, gamma_w = gammaw, gamma_b = gammab, Tau = tau, sigma2 = sigma2, has_intercept = has_intercept, clustermeancentered = centeredwithincluster)

}


