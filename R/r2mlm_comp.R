#' Compute R-squared values for two hierarchical multilevel models.
#'
#' \code{r2mlm_comp} reads in two multilevel models (MLM) generated using
#' \code{\link[lme4]{lmer}} or \code{\link[nlme]{nlme}}, and outputs all
#' R-squared measures for both models. It also produces side-by-side graphical
#' comparisons of the R-squared measures for Model A vs. Model B, that can be
#' used to visualize changes in each measure across models.
#'
#' #' Details: assumes cwc, hierarchy/nesting, both models lmer or both nlme, no mixing and matching
#'
#' @param modelA,modelB Models generated using \code{\link[lme4]{lmer}} or
#'   \code{\link[nlme]{nlme}}.
#'
#' @return If inputs are valid models, then the output will be a list and
#'   associated graphical representation of R-squared decompositions. If models
#'   are not valid, it will return an error prompting the user to input valid
#'   models.
#'
#' @examples
#' # Using lme4 for your model
#'
#' modelA_lme4 <- lmer(popular ~ 1 + extravCWC + texp + (extravCWC|class), data
#' = popularity, REML = TRUE)
#'
#' modelB_lme4 <- lmer(popular ~ 1 + extravCWC + sexCWC + texp +
#' (extravCWC|class), data = popularity, REML = TRUE)
#'
#' r2mlm_comp(modelA_lme4, modelB_lme4)
#'
#' # Using nlme for your model
#'
#' modelA_nlme <- lme(popular ~ 1 + extravCWC + texp,
#'                    random = ~ 1 + extravCWC|class,
#'                    data = popularity,
#'                    method = "REML")
#'
#' modelB_nlme <- lme(popular ~ 1 + extravCWC + sexCWC + texp,
#'                    random = ~ 1 + extravCWC|class,
#'                    data = popularity,
#'                    method = "REML")
#'
#' r2mlm_comp(modelA_nlme, modelB_nlme)
#'
#'
#' @family r2mlm model comparison functions
#'
#' @importFrom lme4 fortify.merMod ranef fixef VarCorr getME
#' @importFrom nlme asOneFormula
#' @importFrom magrittr %>%
#' @importFrom stats terms formula model.frame
#'
#' @export

# 1 r2mlm_comp_wrapper -----------------------------------------------------

r2mlm_comp <- function(modelA, modelB) {

  if (typeof(modelA) != typeof(modelB)) {
    stop("Your models must be of the same type, either both lme4 or both nlme models.")
  }

  if (typeof(modelA) == "list") {
    r2mlm_comp_nlme(modelA, modelB)
  } else if (typeof(modelA) == "S4") {
    r2mlm_comp_lmer(modelA, modelB)
  } else {
    stop("Error. You must input models generated using either the lme4 or nlme package.")
  }

}

# 2 r2mlm_comp_lmer --------------------------------------------------------

r2mlm_comp_lmer <- function(modelA, modelB) {

  # r2mlm_ wrapper sub-functions, but for modelA

  # Step 1: pull data

  data <- lme4::fortify.merMod(modelA)

  # Step 2: has_intercept

  if ((terms(modelA) %>% attr("intercept")) == 1) {
    has_intercept = TRUE
  } else {
    has_intercept = FALSE
  }

  # Step 3: Pull l1 var names into a variable called l1_vars_A, and vice versa for l2 var names

  # i) Pull all variable names from the formula
  all_vars <- all.vars(formula(modelA))
  formula_length <- length(all_vars) # this returns the number of elements in the all_vars list

  # ii) Create vectors to fill
  l1_vars_A <- c()
  l2_vars_A <- c()

  # iii) Sort variables into l1_vars_A and l2_vars_A

  # (a) Pull your temp dataset to work with

  temp_data <- model.frame(modelA) # this pulls a df of all variables and values used in the modelA

  # (b) Isolate the largest group from the dataframe, which you'll use to test variances of variables to sort into l1 and l2 lists
  number <- temp_data %>%
    dplyr::group_by_at(formula_length) %>% #group by the clustering variable, which is the last variable in the df (by virtue of how formula(modelA) works, where it pulls out the formula, and the last variable is on the other side of the |, i.e., the clustering variable)
    dplyr::count() %>%
    dplyr::ungroup() %>% #have to ungroup because otherwise top_n will return n rows from each group, rather than n groups
    dplyr::top_n(1) %>% # returns the ID and N of the largest group
    dplyr::pull(1) # returns the number of the group that you'll use for your variance check

  # (c) Filter temp_data by the number you extracted in 3b
  temp_data_number <- temp_data %>%
    dplyr::filter(temp_data[formula_length] == number) #temp_data[formula_length] is the column that holds the clustering variable

  # (d) Iterate through temp_data_number, calculating the variance for each variable in all_vars, and then sorting by whether variance is 0 (l2) or non-zero (l1)
  x <- 2 # setting a counter overall, starting at 2 to skip the outcome variable (which is otherwise var1 in l1_vars_A)
  l1_counter <- 1 # setting a counter for adding to l1_vars_A list
  l2_counter <- 1 # setting a counter for adding to l2_vars_A list
  while (x < formula_length) {
    if (var(temp_data_number[x]) == 0) {
      l2_vars_A[l2_counter] <- names(temp_data_number[x])
      l2_counter <- l2_counter + 1
    } else {
      l1_vars_A[l1_counter] <- names(temp_data_number[x])
      l1_counter <- l1_counter + 1
    }
    x <- x + 1 # iterate the counter
  }

  # Step 4: pull variable names for L1 predictors with random slopes into a variable called random_slope_vars_A

  temp_cov_list <- lme4::ranef(modelA)[[1]]

  # determine where to start pulling from the temp_cov_list, depending on whether you need to bypass the `(Intercept)` column
  if (has_intercept == 1) {
    running_count <- 2
  } else {
    running_count <- 1
  }

  random_slope_vars_A <- c()
  x <- 1 # counter for indexing in random_slope_vars_A

  # running count less than or equal to list length, so it traverses the entire list (doesn't leave last element off)
  while (running_count <= length(temp_cov_list)) {
    random_slope_vars_A[x] <- names(temp_cov_list[running_count])
    x <- x + 1
    running_count <- running_count + 1
  }

  # Step 5: determine value of centeredwithincluster

  for(var in l1_vars_A) {

    # Sum the l1 column at hand (var in l1_vars_A)
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
      break # so break out of the for loop because if at least one L1 var is not CWC, then the variables will need to be centered by the r2mlm_ function
    }
  }

  # Step 5: pull column numbers for _covs variables
  # 5a) within_covs (l1 variables)
  # for (each value in l1_vars_A list) {match(value, names(data))}
  within_A <- c()
  i = 0
  for (var in l1_vars_A) {
    i = i + 1
    tmp <- match(var, names(data))
    within_A[[i]] <- tmp
  }

  # 5b) pull column numbers for between_covs (l2 variables)
  between_A <- c()
  i = 1
  for (var in l2_vars_A) {
    tmp <- match(var, names(data))
    between_A[[i]] <- tmp
    i = i + 1
  }

  # 5c) pull column numbers for random_covs (l1 variables with random slopes)
  random_A<- c()
  i = 1
  for (var in random_slope_vars_A) {
    tmp <- match(var, names(data))
    random_A[[i]] <- tmp
    i = i + 1
  }

  # Step 6: pull gamma values (fixed slopes)
  # 6a) gamma_w, fixed slopes for L1 variables (from l1_vars_A list)
  gammaw_A <- c()
  i = 1
  for (var in l1_vars_A) {
    gammaw_A[[i]] <- lme4::fixef(modelA)[var]
    i = i + 1
  }

  # 6b) gamma_b, intercept value if hasintercept = TRUE, and fixed slopes for L2 variables (from between list)
  gammab_A <- c()
  if (has_intercept == TRUE) {
    gammab_A[[1]] <- lme4::fixef(modelA)[1]
    i = 2
  } else {
    i = 1
  }
  for (var in l2_vars_A) {
    gammab_A[[i]] <- lme4::fixef(modelA)[var]
    i = i + 1
  }

  # Step 7: Tau matrix, results from VarCorr(modelA)

  vcov <- lme4::VarCorr(modelA)
  tau_A <- as.matrix(Matrix::bdiag(vcov))

  # Step 8: sigma^2 value, Rij

  sigma2_A <- lme4::getME(modelA, "sigma")^2





  # r2mlm_ sub-functions, for model B

  # Step 2: has_intercept

  if ((terms(modelB) %>% attr("intercept")) == 1) {
    has_intercept = TRUE
  } else {
    has_intercept = FALSE
  }

  # Step 3: Pull l1 var names into a variable called l1_vars_B, and vice versa for l2 var names

  # i) Pull all variable names from the formula
  all_vars <- all.vars(formula(modelB))
  formula_length <- length(all_vars) # this returns the number of elements in the all_vars list

  # ii) Create vectors to fill
  l1_vars_B <- c()
  l2_vars_B <- c()

  # iii) Sort variables into l1_vars_B and l2_vars_B

  # (a) Pull your temp dataset to work with

  temp_data <- model.frame(modelB) # this pulls a df of all variables and values used in the modelB

  # (b) Isolate the largest group from the dataframe, which you'll use to test variances of variables to sort into l1 and l2 lists
  number <- temp_data %>%
    dplyr::group_by_at(formula_length) %>% #group by the clustering variable, which is the last variable in the df (by virtue of how formula(modelB) works, where it pulls out the formula, and the last variable is on the other side of the |, i.e., the clustering variable)
    dplyr::count() %>%
    dplyr::ungroup() %>% #have to ungroup because otherwise top_n will return n rows from each group, rather than n groups
    dplyr::top_n(1) %>% # returns the ID and N of the largest group
    dplyr::pull(1) # returns the number of the group that you'll use for your variance check

  # (c) Filter temp_data by the number you extracted in 3b
  temp_data_number <- temp_data %>%
    dplyr::filter(temp_data[formula_length] == number) #temp_data[formula_length] is the column that holds the clustering variable

  # (d) Iterate through temp_data_number, calculating the variance for each variable in all_vars, and then sorting by whether variance is 0 (l2) or non-zero (l1)
  x <- 2 # setting a counter overall, starting at 2 to skip the outcome variable (which is otherwise var1 in l1_vars_B)
  l1_counter <- 1 # setting a counter for adding to l1_vars_B list
  l2_counter <- 1 # setting a counter for adding to l2_vars_B list
  while (x < formula_length) {
    if (var(temp_data_number[x]) == 0) {
      l2_vars_B[l2_counter] <- names(temp_data_number[x])
      l2_counter <- l2_counter + 1
    } else {
      l1_vars_B[l1_counter] <- names(temp_data_number[x])
      l1_counter <- l1_counter + 1
    }
    x <- x + 1 # iterate the counter
  }

  # Step 4: pull variable names for L1 predictors with random slopes into a variable called random_slope_vars_B

  temp_cov_list <- lme4::ranef(modelB)[[1]]

  # determine where to start pulling from the temp_cov_list, depending on whether you need to bypass the `(Intercept)` column
  if (has_intercept == 1) {
    running_count <- 2
  } else {
    running_count <- 1
  }

  random_slope_vars_B <- c()
  x <- 1 # counter for indexing in random_slope_vars_B

  # running count less than or equal to list length, so it traverses the entire list (doesn't leave last element off)
  while (running_count <= length(temp_cov_list)) {
    random_slope_vars_B[x] <- names(temp_cov_list[running_count])
    x <- x + 1
    running_count <- running_count + 1
  }

  # Step 5: determine value of centeredwithincluster

  for(var in l1_vars_B) {

    # Sum the l1 column at hand (var in l1_vars_B)
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
      break # so break out of the for loop because if at least one L1 var is not CWC, then the variables will need to be centered by the r2mlm_ function
    }
  }

  # Step 5: pull column numbers for _covs variables
  # 5a) within_covs (l1 variables)
  # for (each value in l1_vars_B list) {match(value, names(data))}
  within_B <- c()
  i = 0
  for (var in l1_vars_B) {
    i = i + 1
    tmp <- match(var, names(data))
    within_B[[i]] <- tmp
  }

  # 5b) pull column numbers for between_covs (l2 variables)
  between_B <- c()
  i = 1
  for (var in l2_vars_B) {
    tmp <- match(var, names(data))
    between_B[[i]] <- tmp
    i = i + 1
  }

  # 5c) pull column numbers for random_covs (l1 variables with random slopes)
  random_B<- c()
  i = 1
  for (var in random_slope_vars_B) {
    tmp <- match(var, names(data))
    random_B[[i]] <- tmp
    i = i + 1
  }

  # Step 6: pull gamma values (fixed slopes)
  # 6a) gamma_w, fixed slopes for L1 variables (from l1_vars_B list)
  gammaw_B <- c()
  i = 1
  for (var in l1_vars_B) {
    gammaw_B[[i]] <- lme4::fixef(modelB)[var]
    i = i + 1
  }

  # 6b) gamma_b, intercept value if hasintercept = TRUE, and fixed slopes for L2 variables (from between list)
  gammab_B <- c()
  if (has_intercept == TRUE) {
    gammab_B[[1]] <- lme4::fixef(modelB)[1]
    i = 2
  } else {
    i = 1
  }
  for (var in l2_vars_B) {
    gammab_B[[i]] <- lme4::fixef(modelB)[var]
    i = i + 1
  }

  # Step 7: Tau matrix, results from VarCorr(modelB)

  vcov <- lme4::VarCorr(modelB)
  tau_B <- as.matrix(Matrix::bdiag(vcov))

  # Step 8: sigma^2 value, Rij

  sigma2_B <- lme4::getME(modelB, "sigma")^2


  r2mlm_comp_manual(as.data.frame(data), within_covs_modA = within_A, between_covs_modA = between_A, random_covs_modA = random_A, gamma_w_modA = gammaw_A, gamma_b_modA = gammab_A, Tau_modA = tau_A, sigma2_modA = sigma2_A, within_covs_modB = within_B, between_covs_modB = between_B, random_covs_modB = random_B, gamma_w_modB = gammaw_B, gamma_b_modB = gammab_B, Tau_modB = tau_B, sigma2_modB = sigma2_B)

}

# 3 r2mlm_comp_nlme --------------------------------------------------------

r2mlm_comp_nlme <- function(modelA, modelB) {

  # EXTRACT FOR MODEL A

  # Step 1: pull data

  data <- modelA$data

  # Step 2: has_intercept

  if ((terms(modelA) %>% attr("intercept")) == 1) {
    has_intercept = TRUE
  } else {
    has_intercept = FALSE
  }

  # Step 3: Pull l1 var names into a variable called l1_vars_A, and vice versa for l2 var names

  # i) Pull all variable names from the formula
  all_vars <- all.vars(formula(modelA))

  # in nlme, formula(modelA) doesn't return grouping var, but we'll need that later on, so let's grab it here
  full_formula <- all.vars(asOneFormula(modelA))
  grouping_var <- full_formula[length(full_formula)]

  # Then add the grouping var to list of all variables, and calculate formula length (for later use, to iterate)
  all_vars[length(all_vars) + 1] <- grouping_var
  formula_length <- length(all_vars) # this returns the number of elements in the all_vars list

  # ii) Create vectors to fill
  l1_vars_A <- c()
  l2_vars_A <- c()

  # iii) Sort variables into l1_vars_A and l2_vars_A

  # (a) Pull your temp dataset to work with

  temp_data <- data %>%
    dplyr::select(tidyselect::all_of(all_vars)) # this pulls a df of all variables and values used in the modelA

  # (b) Isolate the largest group from the dataframe, which you'll use to test variances of variables to sort into l1 and l2 lists
  number <- temp_data %>%
    dplyr::group_by_at(formula_length) %>% #group by the clustering variable, which is the last variable in the df (by virtue of how formula(modelA) works, where it pulls out the formula, and the last variable is on the other side of the |, i.e., the clustering variable)
    dplyr::count() %>%
    dplyr::ungroup() %>% #have to ungroup because otherwise top_n will return n rows from each group, rather than n groups
    dplyr::top_n(1) %>% # returns the ID and N of the largest group
    dplyr::pull(1) # returns the number of the group that you'll use for your variance check

  # (c) Filter temp_data by the number you extracted in 3b
  temp_data_number <- temp_data %>%
    dplyr::ungroup() %>%
    dplyr::filter(temp_data[formula_length] == number) #temp_data[formula_length] is the column that holds the clustering variable

  # (d) Iterate through temp_data_number, calculating the variance for each variable in all_vars, and then sorting by whether variance is 0 (l2) or non-zero (l1)
  x <- 2 # setting a counter overall, starting at 2 to skip the outcome variable (which is otherwise var1 in l1_vars_A)
  l1_counter <- 1 # setting a counter for adding to l1_vars_A list
  l2_counter <- 1 # setting a counter for adding to l2_vars_A list
  while (x < formula_length) {
    if (var(temp_data_number[x]) == 0) {
      l2_vars_A[l2_counter] <- names(temp_data_number[x])
      l2_counter <- l2_counter + 1
    } else {
      l1_vars_A[l1_counter] <- names(temp_data_number[x])
      l1_counter <- l1_counter + 1
    }
    x <- x + 1 # iterate the counter
  }

  # Step 4: pull variable names for L1 predictors with random slopes into a variable called random_slope_vars_A

  temp_cov_list <- nlme::ranef(modelA)

  # determine where to start pulling from the temp_cov_list, depending on whether you need to bypass the `(Intercept)` column
  if (has_intercept == 1) {
    running_count <- 2
  } else {
    running_count <- 1
  }

  random_slope_vars_A <- c()
  x <- 1 # counter for indexing in random_slope_vars_A

  # running count less than or equal to list length, so it traverses the entire list (doesn't leave last element off)
  while (running_count <= length(temp_cov_list)) {
    random_slope_vars_A[x] <- names(temp_cov_list[running_count])
    x <- x + 1
    running_count <- running_count + 1
  }

  # Step 5: determine value of centeredwithincluster

  for(var in l1_vars_A) {

    # Sum the l1 column at hand (var in l1_vars_A)
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
      break # so break out of the for loop because if at least one L1 var is not CWC, then the variables will need to be centered by the r2mlm_ function
    }
  }

  # Step 5: pull column numbers for _covs variables
  # 5a) within_covs (l1 variables)
  # for (each value in l1_vars_A list) {match(value, names(data))}
  within_A <- c()
  i = 0
  for (var in l1_vars_A) {
    i = i + 1
    tmp <- match(var, names(data))
    within_A[[i]] <- tmp
  }

  # 5b) pull column numbers for between_covs (l2 variables)
  between_A <- c()
  i = 1
  for (var in l2_vars_A) {
    tmp <- match(var, names(data))
    between_A[[i]] <- tmp
    i = i + 1
  }

  # 5c) pull column numbers for random_covs (l1 variables with random slopes)
  random_A <- c()
  i = 1
  for (var in random_slope_vars_A) {
    tmp <- match(var, names(data))
    random_A[[i]] <- tmp
    i = i + 1
  }

  # Step 6: pull gamma values (fixed slopes)
  # 6a) gamma_w, fixed slopes for L1 variables (from l1_vars_A list)
  gammaw_A <- c()
  i = 1
  for (var in l1_vars_A) {
    gammaw_A[[i]] <- nlme::fixef(modelA)[var]
    i = i + 1
  }

  # 6b) gamma_b, intercept value if hasintercept = TRUE, and fixed slopes for L2 variables (from between list)
  gammab_A <- c()
  if (has_intercept == TRUE) {
    gammab_A[[1]] <- nlme::fixef(modelA)[1]
    i = 2
  } else {
    i = 1
  }
  for (var in l2_vars_A) {
    gammab_A[[i]] <- nlme::fixef(modelA)[var]
    i = i + 1
  }

  # Step 7: Tau matrix

  tau_A <- nlme::getVarCov(modelA)

  # Step 8: sigma^2 value, Rij

  sigma2_A <- modelA$sigma


  # EXTRACT FOR MODEL B

  if ((terms(modelB) %>% attr("intercept")) == 1) {
    has_intercept = TRUE
  } else {
    has_intercept = FALSE
  }

  # Step 3: Pull l1 var names into a variable called l1_vars_B, and vice versa for l2 var names

  # i) Pull all variable names from the formula
  all_vars <- all.vars(formula(modelB))

  # in nlme, formula(modelB) doesn't return grouping var, but we'll need that later on, so let's grab it here
  full_formula <- all.vars(asOneFormula(modelB))
  grouping_var <- full_formula[length(full_formula)]

  # Then add the grouping var to list of all variables, and calculate formula length (for later use, to iterate)
  all_vars[length(all_vars) + 1] <- grouping_var
  formula_length <- length(all_vars) # this returns the number of elements in the all_vars list

  # ii) Create vectors to fill
  l1_vars_B <- c()
  l2_vars_B <- c()

  # iii) Sort variables into l1_vars_B and l2_vars_B

  # (a) Pull your temp dataset to work with

  temp_data <- data %>%
    dplyr::select(tidyselect::all_of(all_vars)) # this pulls a df of all variables and values used in the modelB

  # (b) Isolate the largest group from the dataframe, which you'll use to test variances of variables to sort into l1 and l2 lists
  number <- temp_data %>%
    dplyr::group_by_at(formula_length) %>% #group by the clustering variable, which is the last variable in the df (by virtue of how formula(modelB) works, where it pulls out the formula, and the last variable is on the other side of the |, i.e., the clustering variable)
    dplyr::count() %>%
    dplyr::ungroup() %>% #have to ungroup because otherwise top_n will return n rows from each group, rather than n groups
    dplyr::top_n(1) %>% # returns the ID and N of the largest group
    dplyr::pull(1) # returns the number of the group that you'll use for your variance check

  # (c) Filter temp_data by the number you extracted in 3b
  temp_data_number <- temp_data %>%
    dplyr::ungroup() %>%
    dplyr::filter(temp_data[formula_length] == number) #temp_data[formula_length] is the column that holds the clustering variable

  # (d) Iterate through temp_data_number, calculating the variance for each variable in all_vars, and then sorting by whether variance is 0 (l2) or non-zero (l1)
  x <- 2 # setting a counter overall, starting at 2 to skip the outcome variable (which is otherwise var1 in l1_vars_B)
  l1_counter <- 1 # setting a counter for adding to l1_vars_B list
  l2_counter <- 1 # setting a counter for adding to l2_vars_B list
  while (x < formula_length) {
    if (var(temp_data_number[x]) == 0) {
      l2_vars_B[l2_counter] <- names(temp_data_number[x])
      l2_counter <- l2_counter + 1
    } else {
      l1_vars_B[l1_counter] <- names(temp_data_number[x])
      l1_counter <- l1_counter + 1
    }
    x <- x + 1 # iterate the counter
  }

  # Step 4: pull variable names for L1 predictors with random slopes into a variable called random_slope_vars_B

  temp_cov_list <- nlme::ranef(modelB)

  # determine where to start pulling from the temp_cov_list, depending on whether you need to bypass the `(Intercept)` column
  if (has_intercept == 1) {
    running_count <- 2
  } else {
    running_count <- 1
  }

  random_slope_vars_B <- c()
  x <- 1 # counter for indexing in random_slope_vars_B

  # running count less than or equal to list length, so it traverses the entire list (doesn't leave last element off)
  while (running_count <= length(temp_cov_list)) {
    random_slope_vars_B[x] <- names(temp_cov_list[running_count])
    x <- x + 1
    running_count <- running_count + 1
  }

  # Step 5: determine value of centeredwithincluster

  for(var in l1_vars_B) {

    # Sum the l1 column at hand (var in l1_vars_B)
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
      break # so break out of the for loop because if at least one L1 var is not CWC, then the variables will need to be centered by the r2mlm_ function
    }
  }

  # Step 5: pull column numbers for _covs variables
  # 5a) within_covs (l1 variables)
  # for (each value in l1_vars_B list) {match(value, names(data))}
  within_B <- c()
  i = 0
  for (var in l1_vars_B) {
    i = i + 1
    tmp <- match(var, names(data))
    within_B[[i]] <- tmp
  }

  # 5b) pull column numbers for between_covs (l2 variables)
  between_B <- c()
  i = 1
  for (var in l2_vars_B) {
    tmp <- match(var, names(data))
    between_B[[i]] <- tmp
    i = i + 1
  }

  # 5c) pull column numbers for random_covs (l1 variables with random slopes)
  random_B <- c()
  i = 1
  for (var in random_slope_vars_B) {
    tmp <- match(var, names(data))
    random_B[[i]] <- tmp
    i = i + 1
  }

  # Step 6: pull gamma values (fixed slopes)
  # 6a) gamma_w, fixed slopes for L1 variables (from l1_vars_B list)
  gammaw_B <- c()
  i = 1
  for (var in l1_vars_B) {
    gammaw_B[[i]] <- fixef(modelB)[var]
    i = i + 1
  }

  # 6b) gamma_b, intercept value if hasintercept = TRUE, and fixed slopes for L2 variables (from between list)
  gammab_B <- c()
  if (has_intercept == TRUE) {
    gammab_B[[1]] <- nlme::fixef(modelB)[1]
    i = 2
  } else {
    i = 1
  }
  for (var in l2_vars_B) {
    gammab_B[[i]] <- nlme::fixef(modelB)[var]
    i = i + 1
  }

  # Step 7: Tau matrix

  tau_B <- nlme::getVarCov(modelB)

  # Step 8: sigma^2 value, Rij

  sigma2_B <- modelB$sigma

  # Step 9: input everything into r2mlm_

  r2mlm_comp_manual(as.data.frame(data), within_covs_modA = within_A, between_covs_modA = between_A, random_covs_modA = random_A, gamma_w_modA = gammaw_A, gamma_b_modA = gammab_A, Tau_modA = tau_A, sigma2_modA = sigma2_A, within_covs_modB = within_B, between_covs_modB = between_B, random_covs_modB = random_B, gamma_w_modB = gammaw_B, gamma_b_modB = gammab_B, Tau_modB = tau_B, sigma2_modB = sigma2_B)

}
