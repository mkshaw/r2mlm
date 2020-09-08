#' Compute R-squared differences between two multilevel models, automatically
#' inputting parameter estimates.
#'
#' \code{r2mlm_comp} reads in two multilevel models (MLMs) (generated using
#' \code{\link[lme4]{lmer}} or \code{\link[nlme]{nlme}}) under comparison
#' (designated Model A and Model B), and outputs all R-squared measures in the
#' Rights and Sterba (2019) framework for both models, as well as R-squared
#' differences between the two models. Definitions of these R-squared difference
#' measures are provided in Rights & Sterba (2020) Table 1; importantly, to
#' detect the impact of a specific kind of term (e.g., the kind of term added to
#' Model A to form Model B), a particular target single-source R-squared
#' difference measure from this framework is used. For instructions on how to
#' identify which target single-source R-squared difference measure to interpret
#' to detect the impact of which kind of term that distinguishes Model A from B,
#' see Rights and Sterba (2020) Table 2. Additionally, this function produces
#' side-by-side graphical comparisons of the R-squared measures for Model A vs.
#' Model B that can be used to visualize changes in each measure across models.
#' This function assumes all level-1 predictors are cluster-mean-centered, for
#' reasons described in Rights & Sterba (2020). Any number of level-1 and/or
#' level-2 predictors is supported and any of the level-1 predictors can have
#' random slopes. This function can be used with either the hierarchical or the
#' simultaneous model-building approach described in Rights and Sterba (2020).
#' This function can be used with either nested or non-nested model comparisons
#' (in which R-squared estimates for Model A are subtracted from those for Model
#' B).
#'
#' Assumes that both models are fit with lmer or both models are fit with nlme.
#'
#' @param modelA,modelB Models generated using \code{\link[lme4]{lmer}} or
#'   \code{\link[nlme]{nlme}}.
#'
#' @return If the inputs are valid models, then the output will be a list and
#'   associated graphical representation of R-squared decompositions. If the
#'   models are not valid, the function will return an error prompting the user
#'   to input valid models.
#'
#' @examples
#' # Using lme4 for your model
#' # The "bobyqa" optimizer is required for these particular models to converge
#'
#' # Model A, no "salary" components included
#'
#' modelA_lme4 <- lmer(satisfaction ~ 1 + control_c + control_m + s_t_ratio + (1
#' + control_c | schoolID), data = teachsat, REML = TRUE, control =
#' lmerControl(optimizer = "bobyqa"))
#'
#' # Model B, full model with "salary" components included
#'
#' modelB_lme4 <- lmer(satisfaction ~ 1 + salary_c + control_c + salary_m +
#' control_m + s_t_ratio + (1 + salary_c + control_c | schoolID), data =
#' teachsat, REML = TRUE, control = lmerControl(optimizer = "bobyqa"))
#'
#' r2mlm_comp(modelA_lme4, modelB_lme4)
#'
#' # Using nlme for your model
#'
#' # Model A, no "salary" components included
#'
#' modelA_nlme <- lme(satisfaction ~ 1 + control_c + control_m + s_t_ratio,
#'                   random = ~ 1 + control_c | schoolID,
#'                   data = teachsat,
#'                   method = "REML",
#'                   control = lmeControl(optimizer = "bobyqa"))
#'
#' # Model B, full model with "salary" components included
#'
#' modelB_nlme <- lme(satisfaction ~ 1 + salary_c + control_c + salary_m +
#'                   control_m + s_t_ratio,
#'                   random = ~ 1 + salary_c + control_c | schoolID,
#'                   data = teachsat,
#'                   method = "REML",
#'                   control = lmeControl(optimizer = "bobyqa"))
#'
#' r2mlm_comp(modelA_nlme, modelB_nlme)
#'
#' @seealso \href{https://doi.org/10.1037/met0000184}{Rights, J. D., & Sterba,
#'   S. K. (2019). Quantifying explained variance in multilevel models: An
#'   integrative framework for defining R-squared measures. Psychological
#'   Methods, 24(3), 309–338.}
#' @seealso \href{https://doi.org/10.1080/00273171.2019.1660605}{Rights, J. D., &
#'   Sterba, S. K. (2020). New recommendations on the use of R-squared
#'   differences in multilevel model comparisons. Multivariate Behavioral
#'   Research.}
#'
#' @family r2mlm model comparison functions
#'
#' @importFrom lme4 fortify.merMod ranef fixef VarCorr getME
#' @importFrom nlme asOneFormula
#' @importFrom magrittr %>%
#' @importFrom stats terms formula model.frame
#' @importFrom stringr str_split_fixed
#' @importFrom rlang := .data
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

  # Step 0: throw error if model contains higher-order terms
  temp_formula <- formula(modelA)
  grepl_array <- grepl("I(", temp_formula, fixed = TRUE)

  for (bool in grepl_array) {
    if (bool == TRUE) {
      stop("Higher-order terms created with I() syntax are not currently accepted. To include a higher-order term such as x^2 or x^3, you must manually include them as separate columns in your dataset.")
    }
  }

  # r2mlm_ wrapper sub-functions, but for modelA

  # Step 1: pull data

  data <- fortify.merMod(modelA)

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
  l1_vars_A <- c()
  l2_vars_A <- c()

  # iv) Sort variables into l1_vars_A and l2_vars_A

  # (a) Pull your temp dataset to work with

  temp_data <- model.frame(modelA) # this pulls a df of all variables and values used in the modelA

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
      dplyr::select(cluster_variable, variable) %>%
      dplyr::group_map(~ var(.))

    # variable to track variance
    variance_tracker <- 0

    # add up the variance from each cluster
    for (i in t) {
      variance_tracker <- variance_tracker + i
    }

    # if the sum of variance is 0, then each cluster has 0 variance, so it's an L2 variable
    if (variance_tracker == 0) {
      l2_vars_A[l2_counter] <- variable
      l2_counter <- l2_counter + 1
    } else {
      l1_vars_A[l1_counter] <- variable
      l1_counter <- l1_counter + 1
    }

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

  interaction_vars_A <- c()
  x <- 1
  for (term in attr(terms(modelA), "term.labels")) {

    if (grepl(":", term) == TRUE) {
      interaction_vars_A[x] <- term
      x <- x + 1
    }

  }

  data <- data %>% dplyr::ungroup() # need to ungroup in order to create new columns

  for (whole in interaction_vars_A) {

    half1 <- str_split_fixed(whole, ":", 2)[1]
    half2 <- str_split_fixed(whole, ":", 2)[2]

    if (!is.na(match(half1, l2_vars_A)) && !is.na(match(half2, l2_vars_A))) {
      l2_vars_A[l2_counter] <- whole
      l2_counter <- l2_counter + 1
    } else {
      l1_vars_A[l1_counter] <- whole
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

  if (is.null(l1_vars_A)) {
    centeredwithincluster <- TRUE
  } else {
    for (variable in l1_vars_A) {

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

  # Step 5: pull column numbers for _covs variables
  # 5a) within_covs (l1 variables)
  # for (each value in l1_vars_A list) {match(value, names(data))}
  within_A <- c()
  i = 0
  for (var in l1_vars_A) {
    i = i + 1
    tmp <- match(var, names(data))
    within_A[i] <- tmp
  }

  # 5b) pull column numbers for between_covs (l2 variables)
  between_A <- c()
  i = 1
  for (var in l2_vars_A) {
    tmp <- match(var, names(data))
    between_A[i] <- tmp
    i = i + 1
  }

  # 5c) pull column numbers for random_covs (l1 variables with random slopes)
  random_A <- c()
  i = 1
  for (var in random_slope_vars_A) {
    tmp <- match(var, names(data))
    random_A[i] <- tmp
    i = i + 1
  }

  # Step 6: pull gamma values (fixed slopes)
  # 6a) gamma_w, fixed slopes for L1 variables (from l1_vars_A list)
  gammaw_A <- c()
  i = 1
  for (var in l1_vars_A) {
    gammaw_A[i] <- lme4::fixef(modelA)[var]
    i = i + 1
  }

  # 6b) gamma_b, intercept value if hasintercept = TRUE, and fixed slopes for L2 variables (from between list)
  gammab_A <- c()
  if (has_intercept == TRUE) {
    gammab_A[1] <- lme4::fixef(modelA)[1]
    i = 2
  } else {
    i = 1
  }
  for (var in l2_vars_A) {
    gammab_A[i] <- lme4::fixef(modelA)[var]
    i = i + 1
  }

  # Step 7: Tau matrix, results from VarCorr(modelA)

  vcov <- lme4::VarCorr(modelA)
  tau_A <- as.matrix(Matrix::bdiag(vcov))

  # Step 8: sigma^2 value, Rij

  sigma2_A <- lme4::getME(modelA, "sigma")^2





  # r2mlm_ sub-functions, for model B

  # Step 0: throw error if model contains higher-order terms
  temp_formula <- formula(modelB)
  grepl_array <- grepl("I(", temp_formula, fixed = TRUE)

  for (bool in grepl_array) {
    if (bool == TRUE) {
      stop("Higher-order terms created with I() syntax are not currently accepted. To include a higher-order term such as x^2 or x^3, you must manually include them as separate columns in your dataset.")
    }
  }


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
  l1_vars_B <- c()
  l2_vars_B <- c()

  # iv) Sort variables into l1_vars_B and l2_vars_B

  # (a) Pull your temp dataset to work with

  temp_data <- model.frame(modelB) # this pulls a df of all variables and values used in the modelB

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
      dplyr::select(cluster_variable, variable) %>%
      dplyr::group_map(~ var(.))

    # variable to track variance
    variance_tracker <- 0

    # add up the variance from each cluster
    for (i in t) {
      variance_tracker <- variance_tracker + i
    }

    # if the sum of variance is 0, then each cluster has 0 variance, so it's an L2 variable
    if (variance_tracker == 0) {
      l2_vars_B[l2_counter] <- variable
      l2_counter <- l2_counter + 1
    } else {
      l1_vars_B[l1_counter] <- variable
      l1_counter <- l1_counter + 1
    }

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

  interaction_vars_B <- c()
  x <- 1
  for (term in attr(terms(modelB), "term.labels")) {

    if (grepl(":", term) == TRUE) {
      interaction_vars_B[x] <- term
      x <- x + 1
    }

  }

  data <- data %>% dplyr::ungroup() # need to ungroup in order to create new columns

  for (whole in interaction_vars_B) {

    half1 <- str_split_fixed(whole, ":", 2)[1]
    half2 <- str_split_fixed(whole, ":", 2)[2]

    if (!is.na(match(half1, l2_vars_B)) && !is.na(match(half2, l2_vars_B))) {
      l2_vars_B[l2_counter] <- whole
      l2_counter <- l2_counter + 1
    } else {
      l1_vars_B[l1_counter] <- whole
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

  if (is.null(l1_vars_B)) {
    centeredwithincluster <- TRUE
  } else {
    for (variable in l1_vars_B) {

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

  # Step 5: pull column numbers for _covs variables
  # 5a) within_covs (l1 variables)
  # for (each value in l1_vars_B list) {match(value, names(data))}
  within_B <- c()
  i = 0
  for (var in l1_vars_B) {
    i = i + 1
    tmp <- match(var, names(data))
    within_B[i] <- tmp
  }

  # 5b) pull column numbers for between_covs (l2 variables)
  between_B <- c()
  i = 1
  for (var in l2_vars_B) {
    tmp <- match(var, names(data))
    between_B[i] <- tmp
    i = i + 1
  }

  # 5c) pull column numbers for random_covs (l1 variables with random slopes)
  random_B <- c()
  i = 1
  for (var in random_slope_vars_B) {
    tmp <- match(var, names(data))
    random_B[i] <- tmp
    i = i + 1
  }

  # Step 6: pull gamma values (fixed slopes)
  # 6a) gamma_w, fixed slopes for L1 variables (from l1_vars_B list)
  gammaw_B <- c()
  i = 1
  for (var in l1_vars_B) {
    gammaw_B[i] <- lme4::fixef(modelB)[var]
    i = i + 1
  }

  # 6b) gamma_b, intercept value if hasintercept = TRUE, and fixed slopes for L2 variables (from between list)
  gammab_B <- c()
  if (has_intercept == TRUE) {
    gammab_B[1] <- lme4::fixef(modelB)[1]
    i = 2
  } else {
    i = 1
  }
  for (var in l2_vars_B) {
    gammab_B[i] <- lme4::fixef(modelB)[var]
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

  # Step 0: throw error if model contains higher-order terms
  temp_formula <- formula(modelA)
  grepl_array <- grepl("I(", temp_formula, fixed = TRUE)

  for (bool in grepl_array) {
    if (bool == TRUE) {
      stop("Higher-order terms created with I() syntax are not currently accepted. To include a higher-order term such as x^2 or x^3, you must manually include them as separate columns in your dataset.")
    }
  }


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

  # in nlme, formula(model) doesn't return grouping var, but we'll need that later on, so let's grab it here
  full_formula <- all.vars(asOneFormula(modelA))
  cluster_variable <- full_formula[length(full_formula)]

  # Then add the grouping var to list of all variables, and calculate formula length (for later use, to iterate)
  all_vars[length(all_vars) + 1] <- cluster_variable
  formula_length <- length(all_vars) # this returns the number of elements in the all_vars list

  # ii) determine whether data is appropriate format. Only the cluster variable can be a factor, for now

  # b) If any of those variables is non-numeric, then throw an error

  all_vars_except_cluster <- all.vars(formula(modelA))

  for (var in all_vars_except_cluster) {

    if (!(class(data[[var]]) == "integer") && !(class(data[[var]]) == "numeric")) {
      stop("Your data must be numeric. Only the cluster variable can be a factor.")
    }

  }

  # iii) Create vectors to fill
  l1_vars_A <- c()
  l2_vars_A <- c()

  # iii) Sort variables into l1_vars_A and l2_vars_A

  # (a) Pull your temp dataset to work with

  temp_data <- data %>%
    dplyr::select(tidyselect::all_of(all_vars)) # this pulls a df of all variables and values used in the modelA

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
      dplyr::select(cluster_variable, variable) %>%
      dplyr::group_map(~ var(.))

    # variable to track variance
    variance_tracker <- 0

    # add up the variance from each cluster
    for (i in t) {
      variance_tracker <- variance_tracker + i
    }

    # if the sum of variance is 0, then each cluster has 0 variance, so it's an L2 variable
    if (variance_tracker == 0) {
      l2_vars_A[l2_counter] <- variable
      l2_counter <- l2_counter + 1
    } else {
      l1_vars_A[l1_counter] <- variable
      l1_counter <- l1_counter + 1
    }

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

  interaction_vars_A <- c()
  x <- 1
  for (term in attr(terms(modelA), "term.labels")) {

    if (grepl(":", term) == TRUE) {
      interaction_vars_A[x] <- term
      x <- x + 1
    }

  }

  data <- data %>% dplyr::ungroup() # need to ungroup in order to create new columns

  for (whole in interaction_vars_A) {

    half1 <- str_split_fixed(whole, ":", 2)[1]
    half2 <- str_split_fixed(whole, ":", 2)[2]

    if (!is.na(match(half1, l2_vars_A)) && !is.na(match(half2, l2_vars_A))) {
      l2_vars_A[l2_counter] <- whole
      l2_counter <- l2_counter + 1
    } else {
      l1_vars_A[l1_counter] <- whole
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

  if (is.null(l1_vars_A)) {
    centeredwithincluster <- TRUE
  } else {
    for (variable in l1_vars_A) {

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
  # for (each value in l1_vars_A list) {match(value, names(data))}
  within_A <- c()
  i = 0
  for (var in l1_vars_A) {
    i = i + 1
    tmp <- match(var, names(data))
    within_A[i] <- tmp
  }

  # 6b) pull column numbers for between_covs (l2 variables)
  between_A <- c()
  i = 1
  for (var in l2_vars_A) {
    tmp <- match(var, names(data))
    between_A[i] <- tmp
    i = i + 1
  }

  # 6c) pull column numbers for random_covs (l1 variables with random slopes)
  random_A <- c()
  i = 1
  for (var in random_slope_vars_A) {
    tmp <- match(var, names(data))
    random_A[i] <- tmp
    i = i + 1
  }

  # Step 7: pull gamma values (fixed slopes)
  # 7a) gamma_w, fixed slopes for L1 variables (from l1_vars_A list)
  gammaw_A <- c()
  i = 1
  for (var in l1_vars_A) {
    gammaw_A[i] <- nlme::fixef(modelA)[var]
    i = i + 1
  }

  # 7b) gamma_b, intercept value if hasintercept = TRUE, and fixed slopes for L2 variables (from between list)
  gammab_A <- c()
  if (has_intercept == TRUE) {
    gammab_A[1] <- nlme::fixef(modelA)[1]
    i = 2
  } else {
    i = 1
  }
  for (var in l2_vars_A) {
    gammab_A[i] <- nlme::fixef(modelA)[var]
    i = i + 1
  }

  # Step 8: Tau matrix

  tau_A <- nlme::getVarCov(modelA)

  # Step 9: sigma^2 value, Rij

  sigma2_A <- modelA$sigma^2


  # EXTRACT FOR MODEL B

  # Step 0: throw error if model contains higher-order terms
  temp_formula <- formula(modelB)
  grepl_array <- grepl("I(", temp_formula, fixed = TRUE)

  for (bool in grepl_array) {
    if (bool == TRUE) {
      stop("Higher-order terms created with I() syntax are not currently accepted. To include a higher-order term such as x^2 or x^3, you must manually include them as separate columns in your dataset.")
    }
  }

  # Step 1: determine intercept

  if ((terms(modelB) %>% attr("intercept")) == 1) {
    has_intercept = TRUE
  } else {
    has_intercept = FALSE
  }

  # Step 3: Pull l1 var names into a variable called l1_vars_B, and vice versa for l2 var names

  # i) Pull all variable names from the formula
  all_vars <- all.vars(formula(modelB))

  # in nlme, formula(model) doesn't return grouping var, but we'll need that later on, so let's grab it here
  full_formula <- all.vars(asOneFormula(modelB))
  cluster_variable <- full_formula[length(full_formula)]

  # Then add the grouping var to list of all variables, and calculate formula length (for later use, to iterate)
  all_vars[length(all_vars) + 1] <- cluster_variable
  formula_length <- length(all_vars) # this returns the number of elements in the all_vars list

  # ii) determine whether data is appropriate format. Only the cluster variable can be a factor, for now

  # b) If any of those variables is non-numeric, then throw an error

  all_vars_except_cluster <- all.vars(formula(modelB))

  for (var in all_vars_except_cluster) {

    if (!(class(data[[var]]) == "integer") && !(class(data[[var]]) == "numeric")) {
      stop("Your data must be numeric. Only the cluster variable can be a factor.")
    }

  }

  # iii) Create vectors to fill
  l1_vars_B <- c()
  l2_vars_B <- c()

  # iii) Sort variables into l1_vars_B and l2_vars_B

  # (a) Pull your temp dataset to work with

  temp_data <- data %>%
    dplyr::select(tidyselect::all_of(all_vars)) # this pulls a df of all variables and values used in the modelB

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
      dplyr::select(cluster_variable, variable) %>%
      dplyr::group_map(~ var(.))

    # variable to track variance
    variance_tracker <- 0

    # add up the variance from each cluster
    for (i in t) {
      variance_tracker <- variance_tracker + i
    }

    # if the sum of variance is 0, then each cluster has 0 variance, so it's an L2 variable
    if (variance_tracker == 0) {
      l2_vars_B[l2_counter] <- variable
      l2_counter <- l2_counter + 1
    } else {
      l1_vars_B[l1_counter] <- variable
      l1_counter <- l1_counter + 1
    }

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

  interaction_vars_B <- c()
  x <- 1
  for (term in attr(terms(modelB), "term.labels")) {

    if (grepl(":", term) == TRUE) {
      interaction_vars_B[x] <- term
      x <- x + 1
    }

  }

  data <- data %>% dplyr::ungroup() # need to ungroup in order to create new columns

  for (whole in interaction_vars_B) {

    half1 <- str_split_fixed(whole, ":", 2)[1]
    half2 <- str_split_fixed(whole, ":", 2)[2]

    if (!is.na(match(half1, l2_vars_B)) && !is.na(match(half2, l2_vars_B))) {
      l2_vars_B[l2_counter] <- whole
      l2_counter <- l2_counter + 1
    } else {
      l1_vars_B[l1_counter] <- whole
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

  if (is.null(l1_vars_B)) {
    centeredwithincluster <- TRUE
  } else {
    for (variable in l1_vars_B) {

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
  # for (each value in l1_vars_B list) {match(value, names(data))}
  within_B <- c()
  i = 0
  for (var in l1_vars_B) {
    i = i + 1
    tmp <- match(var, names(data))
    within_B[i] <- tmp
  }

  # 6b) pull column numbers for between_covs (l2 variables)
  between_B <- c()
  i = 1
  for (var in l2_vars_B) {
    tmp <- match(var, names(data))
    between_B[i] <- tmp
    i = i + 1
  }

  # 6c) pull column numbers for random_covs (l1 variables with random slopes)
  random_B <- c()
  i = 1
  for (var in random_slope_vars_B) {
    tmp <- match(var, names(data))
    random_B[i] <- tmp
    i = i + 1
  }

  # Step 7: pull gamma values (fixed slopes)
  # 7a) gamma_w, fixed slopes for L1 variables (from l1_vars_B list)
  gammaw_B <- c()
  i = 1
  for (var in l1_vars_B) {
    gammaw_B[i] <- fixef(modelB)[var]
    i = i + 1
  }

  # 7b) gamma_b, intercept value if hasintercept = TRUE, and fixed slopes for L2 variables (from between list)
  gammab_B <- c()
  if (has_intercept == TRUE) {
    gammab_B[1] <- nlme::fixef(modelB)[1]
    i = 2
  } else {
    i = 1
  }
  for (var in l2_vars_B) {
    gammab_B[i] <- nlme::fixef(modelB)[var]
    i = i + 1
  }

  # Step 8: Tau matrix

  tau_B <- nlme::getVarCov(modelB)

  # Step 9: sigma^2 value, Rij

  sigma2_B <- modelB$sigma^2

  # Step 10: input everything into r2mlm_

  r2mlm_comp_manual(as.data.frame(data), within_covs_modA = within_A, between_covs_modA = between_A, random_covs_modA = random_A, gamma_w_modA = gammaw_A, gamma_b_modA = gammab_A, Tau_modA = tau_A, sigma2_modA = sigma2_A, within_covs_modB = within_B, between_covs_modB = between_B, random_covs_modB = random_B, gamma_w_modB = gammaw_B, gamma_b_modB = gammab_B, Tau_modB = tau_B, sigma2_modB = sigma2_B)

}
