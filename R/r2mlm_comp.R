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
#' @param data Dataset with rows denoting observations and columns denoting
#'   variables.
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
#'                   control = lmeControl(opt = "optim"))
#'
#' # Model B, full model with "salary" components included
#'
#' modelB_nlme <- lme(satisfaction ~ 1 + salary_c + control_c + salary_m +
#'                   control_m + s_t_ratio,
#'                   random = ~ 1 + salary_c + control_c | schoolID,
#'                   data = teachsat,
#'                   method = "REML",
#'                   control = lmeControl(opt = "optim"))
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
#' @importFrom stats terms formula model.frame na.omit
#' @importFrom stringr str_split_fixed
#' @importFrom rlang := .data
#'
#' @export

# 1 r2mlm_comp_wrapper -----------------------------------------------------

r2mlm_comp <- function(modelA, modelB, data = NULL) {

  # throw error if model contains higher-order terms
  temp_formula <- formula(modelA)
  grepl_array <- grepl("I(", temp_formula, fixed = TRUE)

  for (bool in grepl_array) {
    if (bool == TRUE) {
      stop("Error: r2mlm does not allow for models fit using the I() function; user must thus manually include any desired transformed predictor variables such as x^2 or x^3 as separate columns in dataset.")
    }
  }

  temp_formula <- formula(modelB)
  grepl_array <- grepl("I(", temp_formula, fixed = TRUE)

  for (bool in grepl_array) {
    if (bool == TRUE) {
      stop("Error: r2mlm does not allow for models fit using the I() function; user must thus manually include any desired transformed predictor variables such as x^2 or x^3 as separate columns in dataset.")
    }
  }

  # make sure models are same type, call helper functions
  if (typeof(modelA) != typeof(modelB)) {
    stop("Your models must be of the same type, either both lme4 or both nlme models.")
  }

  if (typeof(modelA) == "list") {
    r2mlm_comp_nlme(modelA, modelB, data)
  } else if (typeof(modelA) == "S4") {
    r2mlm_comp_lmer(modelA, modelB, data)
  } else {
    stop("Error. You must input models generated using either the lme4 or nlme package.")
  }

}

# 2 r2mlm_comp_lmer --------------------------------------------------------

r2mlm_comp_lmer <- function(modelA, modelB, data) {

  # r2mlm_ wrapper sub-functions, but for modelA

  # Step 1) check if model has_intercept

  if ((terms(modelA) %>% attr("intercept")) == 1) {
    has_intercept = TRUE
  } else {
    has_intercept = FALSE
  }

  # Step 2) Pull all variable names from the formula
  all_variables <- all.vars(formula(modelA))
  cluster_variable <- all_variables[length(all_variables)] # pull cluster, we'll need it later

  # Step 3a) pull and prepare data

  if(is.null(data)) {
    data <- check_hierarchical(modelA, modelB, "lme4", cluster_variable)
  } else {
    # get interaction variables from both models, create list to add to dataframe
    interaction_vars_A <- get_interaction_vars(modelA)
    interaction_vars_B <- get_interaction_vars(modelB)

    interaction_vars <- unique(append(interaction_vars_A, interaction_vars_B)) # merge lists and use unique() to avoid duplicates

    # get data (when is.null(data), these steps are happening in check_hierarchical -> prepare_data)
    data <- na.omit(data) %>%
      add_interaction_vars_to_data(interaction_vars, .) %>%
      group_data(cluster_variable, .)

  }

  # Step 3b) determine whether data is appropriate format. Only the cluster variable can be a factor, for now
  # a) Pull all variables except for cluster

  outcome_and_predictors <- all_variables[1:(length(all_variables) - 1)]

  # b) If any of those variables is non-numeric, then throw an error

  for (variable in outcome_and_predictors) {

    if (!(class(data[[variable]]) == "integer") && !(class(data[[variable]]) == "numeric")) {
      stop("Your data must be numeric. Only the cluster variable can be a factor.")
    }

  }

  # Step 4) Fill l1 and l2 vectors

  # * Step 4a) Define variables you'll be sorting
  # take the outcome out of the predictors with [2:length(outcome_and_predictors)], then add interaction terms

  # sometimes the outcome_and_predictors is only an outcome variable (for the null model). If that's the case, then
  #     predictors is null, just call get_interaction_vars just in case

  if (length(outcome_and_predictors) == 1) {
    predictors <- get_interaction_vars(modelA)
  } else {
    predictors <- append(outcome_and_predictors[2:length(outcome_and_predictors)], get_interaction_vars(modelA))
  }

  # * Step 4b) Create and fill vectors
  l1_vars_A <- sort_variables(data, predictors, cluster_variable)$l1_vars
  l2_vars_A <- sort_variables(data, predictors, cluster_variable)$l2_vars

  # Step 5) pull variable names for L1 predictors with random slopes into a variable called random_slope_vars

  random_slope_vars_A <- get_random_slope_vars(modelA, has_intercept, "lme4")

  # Step 6) determine value of centeredwithincluster

  if (is.null(l1_vars_A)) {
    centeredwithincluster <- TRUE
  } else {
    centeredwithincluster <- get_cwc(l1_vars_A, cluster_variable, data)
  }

  # Step 7) pull column numbers for _covs variables
  # 7a) within_covs (l1 variables)
  # for (each value in l1_vars list) {match(value, names(data))}

  within_A <- get_covs(l1_vars_A, data)

  # 7b) pull column numbers for between_covs (l2 variables)

  between_A <- get_covs(l2_vars_A, data)

  # 7c) pull column numbers for random_covs (l1 variables with random slopes)

  random_A <- get_covs(random_slope_vars_A, data)

  # Step 8) pull gamma values (fixed slopes)
  # 8a) gamma_w, fixed slopes for L1 variables (from l1_vars list)
  gammaw_A <- c()
  i = 1
  for (var in l1_vars_A) {
    gammaw_A[i] <- lme4::fixef(modelA)[var]
    i = i + 1
  }

  # 8b) gamma_b, intercept value if hasintercept = TRUE, and fixed slopes for L2 variables (from between list)
  gammab_A <- c()
  if (has_intercept == TRUE) {
    gammab_A[1] <- lme4::fixef(modelA)[1]
    i = 2
  } else {
    i = 1
  }
  for (variable in l2_vars_A) {
    gammab_A[i] <- lme4::fixef(modelA)[variable]
    i = i + 1
  }

  # Step 9) Tau matrix, results from VarCorr(modelA)

  vcov <- lme4::VarCorr(modelA)
  tau_A <- as.matrix(Matrix::bdiag(vcov))

  # Step 10) sigma^2 value, Rij

  sigma2_A <- lme4::getME(modelA, "sigma")^2





  # r2mlm_ sub-functions, for model B

  # Step 1) check if model has_intercept

  if ((terms(modelB) %>% attr("intercept")) == 1) {
    has_intercept = TRUE
  } else {
    has_intercept = FALSE
  }

  # Step 2) Pull all variable names from the formula
  all_variables <- all.vars(formula(modelB))
  cluster_variable <- all_variables[length(all_variables)] # pull cluster, we'll need it later

  # Step 3) determine whether data is appropriate format. Only the cluster variable can be a factor, for now
  # a) Pull all variables except for cluster

  outcome_and_predictors <- all_variables[1:(length(all_variables) - 1)]

  # b) If any of those variables is non-numeric, then throw an error

  for (variable in outcome_and_predictors) {

    if (!(class(data[[variable]]) == "integer") && !(class(data[[variable]]) == "numeric")) {
      stop("Your data must be numeric. Only the cluster variable can be a factor.")
    }

  }

  # Step 4) Fill l1 and l2 vectors

  # * Step 4a) Define variables you'll be sorting
  # take the outcome out of the predictors with [2:length(outcome_and_predictors)], then add interaction terms

  # sometimes the outcome_and_predictors is only an outcome variable (for the null model). If that's the case, then
  #     predictors is null, just call get_interaction_vars just in case

  if (length(outcome_and_predictors) == 1) {
    predictors <- get_interaction_vars(modelB)
  } else {
    predictors <- append(outcome_and_predictors[2:length(outcome_and_predictors)], get_interaction_vars(modelB))
  }

  # * Step 4b) Create and fill vectors
  l1_vars_B <- sort_variables(data, predictors, cluster_variable)$l1_vars
  l2_vars_B <- sort_variables(data, predictors, cluster_variable)$l2_vars

  # Step 5) pull variable names for L1 predictors with random slopes into a variable called random_slope_vars

  random_slope_vars_B <- get_random_slope_vars(modelB, has_intercept, "lme4")

  # Step 6) determine value of centeredwithincluster

  if (is.null(l1_vars_B)) {
    centeredwithincluster <- TRUE
  } else {
    centeredwithincluster <- get_cwc(l1_vars_B, cluster_variable, data)
  }

  # Step 7) pull column numbers for _covs variables
  # 7a) within_covs (l1 variables)
  # for (each value in l1_vars list) {match(value, names(data))}

  within_B <- get_covs(l1_vars_B, data)

  # 7b) pull column numbers for between_covs (l2 variables)

  between_B <- get_covs(l2_vars_B, data)

  # 7c) pull column numbers for random_covs (l1 variables with random slopes)

  random_B <- get_covs(random_slope_vars_B, data)

  # Step 8) pull gamma values (fixed slopes)
  # 8a) gamma_w, fixed slopes for L1 variables (from l1_vars list)
  gammaw_B <- c()
  i = 1
  for (var in l1_vars_B) {
    gammaw_B[i] <- lme4::fixef(modelB)[var]
    i = i + 1
  }

  # 8b) gamma_b, intercept value if hasintercept = TRUE, and fixed slopes for L2 variables (from between list)
  gammab_B <- c()
  if (has_intercept == TRUE) {
    gammab_B[1] <- lme4::fixef(modelB)[1]
    i = 2
  } else {
    i = 1
  }
  for (variable in l2_vars_B) {
    gammab_B[i] <- lme4::fixef(modelB)[variable]
    i = i + 1
  }

  # Step 9) Tau matrix, results from VarCorr(modelB)

  vcov <- lme4::VarCorr(modelB)
  tau_B <- as.matrix(Matrix::bdiag(vcov))

  # Step 10) sigma^2 value, Rij

  sigma2_B <- lme4::getME(modelB, "sigma")^2


  r2mlm_comp_manual(as.data.frame(data), within_covs_modA = within_A, between_covs_modA = between_A, random_covs_modA = random_A, gamma_w_modA = gammaw_A, gamma_b_modA = gammab_A, Tau_modA = tau_A, sigma2_modA = sigma2_A, within_covs_modB = within_B, between_covs_modB = between_B, random_covs_modB = random_B, gamma_w_modB = gammaw_B, gamma_b_modB = gammab_B, Tau_modB = tau_B, sigma2_modB = sigma2_B)

}

# 3 r2mlm_comp_nlme --------------------------------------------------------

r2mlm_comp_nlme <- function(modelA, modelB, data) {

  # EXTRACT FOR MODEL A

  # Step 1) check if modelA has_intercept

  if ((terms(modelA) %>% attr("intercept")) == 1) {
    has_intercept = TRUE
  } else {
    has_intercept = FALSE
  }

  # Step 2) Pull all variable names from the formula
  all_variables <- all.vars(formula(modelA))
  cluster_variable <- nlme::getGroups(modelA) %>% attr("label") # in nlme, formula(model) doesn't return grouping var, but we'll need that later on, so let's grab it here

  # Then add the grouping var to list of all variables, and calculate formula length (for later use, to iterate)
  all_variables[length(all_variables) + 1] <- cluster_variable
  formula_length <- length(all_variables) # this returns the number of elements in the all_variables list

  # Step 3a) pull and prepare data

  if(is.null(data)) {
    data <- check_hierarchical(modelA, modelB, "nlme", cluster_variable)
  } else {
    # get interaction variables from both models, create list to add to dataframe
    interaction_vars_A <- get_interaction_vars(modelA)
    interaction_vars_B <- get_interaction_vars(modelB)

    interaction_vars <- unique(append(interaction_vars_A, interaction_vars_B)) # merge lists and use unique() to avoid duplicates

    # get data (when is.null(data), these steps are happening in check_hierarchical -> prepare_data)
    data <- na.omit(data) %>%
      add_interaction_vars_to_data(interaction_vars, .) %>%
      group_data(cluster_variable, .)

  }

  # * Step 3b) determine whether data is appropriate format. Only the cluster variable can be a factor, for now

  outcome_and_predictors <- all.vars(formula(modelA))

  for (variable in outcome_and_predictors) {

    if (!(class(data[[variable]]) == "integer") && !(class(data[[variable]]) == "numeric")) {
      stop("Your data must be numeric. Only the cluster variable can be a factor.")
    }

  }

  # Step 4) Fill l1 and l2 vectors

  # * Step 4a) Define variables you'll be sorting
  # sometimes the outcome_and_predictors is only an outcome variable (for the null model). If that's the case, then
  #     predictors is null, just call get_interaction_vars just in case

  if (length(outcome_and_predictors) == 1) {
    predictors <- get_interaction_vars(modelA)
  } else {
    predictors <- append(outcome_and_predictors[2:length(outcome_and_predictors)], get_interaction_vars(modelA))
  }

  # * Step 4b) Create and fill vectors
  l1_vars_A <- sort_variables(data, predictors, cluster_variable)$l1_vars
  l2_vars_A <- sort_variables(data, predictors, cluster_variable)$l2_vars

  # Step 5) pull variable names for L1 predictors with random slopes into a variable called random_slope_vars

  random_slope_vars <- get_random_slope_vars(modelA, has_intercept, "nlme")

  # Step 6) determine value of centeredwithincluster

  if (is.null(l1_vars_A)) {
    centeredwithincluster <- TRUE
  } else {
    centeredwithincluster <- get_cwc(l1_vars_A, cluster_variable, data)
  }

  # Step 7) pull column numbers for _covs variables
  # 7a) within_covs (l1 variables)
  # for (each value in l1_vars list) {match(value, names(data))}
  within_A <- get_covs(l1_vars_A, data)

  # 7b) pull column numbers for between_covs (l2 variables)
  between_A <- get_covs(l2_vars_A, data)

  # 7c) pull column numbers for random_covs (l1 variables with random slopes)
  random_A <- get_covs(random_slope_vars, data)

  # Step 8) pull gamma values (fixed slopes)
  # 8a) gamma_w, fixed slopes for L1 variables (from l1_vars list)
  gammaw_A <- c()
  i = 1
  for (variable in l1_vars_A) {
    gammaw_A[i] <- nlme::fixef(modelA)[variable]
    i = i + 1
  }

  # 8b) gamma_b, intercept value if hasintercept = TRUE, and fixed slopes for L2 variables (from between list)
  gammab_A <- c()
  if (has_intercept == TRUE) {
    gammab_A[1] <- nlme::fixef(modelA)[1]
    i = 2
  } else {
    i = 1
  }
  for (variable in l2_vars_A) {
    gammab_A[i] <- nlme::fixef(modelA)[variable]
    i = i + 1
  }

  # Step 9) Tau matrix

  tau_A <- nlme::getVarCov(modelA)

  # Step 10) sigma^2 value, Rij

  sigma2_A <- modelA$sigma^2









  # EXTRACT FOR MODEL B

  # Step 1) check if modelB has_intercept

  if ((terms(modelB) %>% attr("intercept")) == 1) {
    has_intercept = TRUE
  } else {
    has_intercept = FALSE
  }

  # Step 2) Pull all variable names from the formula
  all_variables <- all.vars(formula(modelB))
  cluster_variable <- nlme::getGroups(modelB) %>% attr("label") # in nlme, formula(model) doesn't return grouping var, but we'll need that later on, so let's grab it here

  # Then add the grouping var to list of all variables, and calculate formula length (for later use, to iterate)
  all_variables[length(all_variables) + 1] <- cluster_variable
  formula_length <- length(all_variables) # this returns the number of elements in the all_variables list

  # * Step 3b) determine whether data is appropriate format. Only the cluster variable can be a factor, for now

  outcome_and_predictors <- all.vars(formula(modelB))

  for (variable in outcome_and_predictors) {

    if (!(class(data[[variable]]) == "integer") && !(class(data[[variable]]) == "numeric")) {
      stop("Your data must be numeric. Only the cluster variable can be a factor.")
    }

  }

  # Step 4) Fill l1 and l2 vectors

  # * Step 4a) Define variables you'll be sorting
  # sometimes the outcome_and_predictors is only an outcome variable (for the null model). If that's the case, then
  #     predictors is null, just call get_interaction_vars just in case

  if (length(outcome_and_predictors) == 1) {
    predictors <- get_interaction_vars(modelB)
  } else {
    predictors <- append(outcome_and_predictors[2:length(outcome_and_predictors)], get_interaction_vars(modelB))
  }

  # * Step 4b) Create and fill vectors
  l1_vars_B <- sort_variables(data, predictors, cluster_variable)$l1_vars
  l2_vars_B <- sort_variables(data, predictors, cluster_variable)$l2_vars

  # Step 5) pull variable names for L1 predictors with random slopes into a variable called random_slope_vars

  random_slope_vars <- get_random_slope_vars(modelB, has_intercept, "nlme")

  # Step 6) determine value of centeredwithincluster

  if (is.null(l1_vars_B)) {
    centeredwithincluster <- TRUE
  } else {
    centeredwithincluster <- get_cwc(l1_vars_B, cluster_variable, data)
  }

  # Step 7) pull column numbers for _covs variables
  # 7a) within_covs (l1 variables)
  # for (each value in l1_vars list) {match(value, names(data))}
  within_B <- get_covs(l1_vars_B, data)

  # 7b) pull column numbers for between_covs (l2 variables)
  between_B <- get_covs(l2_vars_B, data)

  # 7c) pull column numbers for random_covs (l1 variables with random slopes)
  random_B <- get_covs(random_slope_vars, data)

  # Step 8) pull gamma values (fixed slopes)
  # 8a) gamma_w, fixed slopes for L1 variables (from l1_vars list)
  gammaw_B <- c()
  i = 1
  for (variable in l1_vars_B) {
    gammaw_B[i] <- nlme::fixef(modelB)[variable]
    i = i + 1
  }

  # 8b) gamma_b, intercept value if hasintercept = TRUE, and fixed slopes for L2 variables (from between list)
  gammab_B <- c()
  if (has_intercept == TRUE) {
    gammab_B[1] <- nlme::fixef(modelB)[1]
    i = 2
  } else {
    i = 1
  }
  for (variable in l2_vars_B) {
    gammab_B[i] <- nlme::fixef(modelB)[variable]
    i = i + 1
  }

  # Step 9) Tau matrix

  tau_B <- nlme::getVarCov(modelB)

  # Step 10) sigma^2 value, Rij

  sigma2_B <- modelB$sigma^2






  # Step 11) input everything into r2mlm_comp_manual

  r2mlm_comp_manual(as.data.frame(data), within_covs_modA = within_A, between_covs_modA = between_A, random_covs_modA = random_A, gamma_w_modA = gammaw_A, gamma_b_modA = gammab_A, Tau_modA = tau_A, sigma2_modA = sigma2_A, within_covs_modB = within_B, between_covs_modB = between_B, random_covs_modB = random_B, gamma_w_modB = gammaw_B, gamma_b_modB = gammab_B, Tau_modB = tau_B, sigma2_modB = sigma2_B)

}
