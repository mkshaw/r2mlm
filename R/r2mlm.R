#' Compute R-squared values for multilevel models, automatically inputting
#' parameter estimates.
#'
#' \code{r2mlm} reads in a multilevel model (MLM) object generated using
#' \code{\link[lme4]{lmer}} or \code{\link[nlme]{nlme}}, and outputs all
#' relevant R-squared measures from the Rights and Sterba (2019) framework of
#' multilevel model R-squared measures, which can be visualized together as a
#' set using the outputted bar chart decompositions of outcome variance. That is,
#' when predictors are cluster-mean-centered, all R-squared measures from Rights
#' & Sterba (2019) Table 1 and decompositions from Rights & Sterba (2019) Figure
#' 1 are outputted. When predictors are not cluster-mean-centered, the total
#' R-squared measures from Rights & Sterba (2019) Table 5, as well as bar chart
#' decompositions are outputted. Any number of level-1 and/or level-2 predictors
#' is supported. Any of the level-1 predictors can have random slopes.
#'
#' \code{r2mlm} first determines whether a given model was generated using
#' \code{\link[lme4]{lmer}} or \code{\link[nlme]{nlme}}, then passes the model
#' to helper functions that pull the raw data and parameter estimates from the
#' model, and pass that information to \code{\link{r2mlm_manual}}.
#'
#' Previous MLM literature has offered two perspectives on how to treat variance
#' attributable to random intercepts and slopes, called the “marginal” and
#' “conditional” approaches (e.g., Edwards et al., 2008; Orelien & Edwards,
#' 2008; Vonesh & Chinchilli, 1997; Wang & Schaalje, 2009; Xu, 2003). In the
#' marginal approach, all variance attributable to predictor and cluster mean
#' variation is treated as unexplained. In the conditional approach, all
#' variance attributable to predictors and cluster mean variation is treated as
#' explained. This package offers researchers access to both the marginal and
#' conditional approaches. There are 5 marginal measures: f1_total, f2_total,
#' f_total, f1_within, and f2_between. The other 7 measures are conditional:
#' v_total, m_total, fv_total, fvm_total, v_within, f1v_within, and m_between.
#'
#' @param model A model generated using \code{\link[lme4]{lmer}} or
#'   \code{\link[nlme]{nlme}}. Note that models using \code{lmer} must specify
#'   random effects at the end of the model, like so: \code{outcome ~ 1 +
#'   fixed_effects + (random_effects | cluster_variable)}. Anything else (e.g.,
#'   \code{outcome ~ 1 + (random_effects | cluster_variable) + fixed_effects})
#'   will not work.
#' @param bargraph Optional bar graph output, default is TRUE.
#'
#' @return If the input is a valid model, then the output will be a list and
#'   associated graphical representation of R-squared decompositions. If the
#'   model is not valid, it will return an error prompting the user to input a
#'   valid model.
#'
#' @examples
#' # Using lme4 for your model
#'
#' # The "bobyqa" optimizer is required for this particular model to converge
#'
#' model_lme4 <- lmer(satisfaction ~ 1 + salary_c + control_c + salary_m + control_m +
#' s_t_ratio + (1 + salary_c + control_c| schoolID), data = teachsat, REML =
#' TRUE, control = lmerControl(optimizer = "bobyqa"))
#'
#' r2mlm(model_lme4)
#'
#' # Using nlme for your model
#'
#' model_nlme <- lme(satisfaction ~ 1 + salary_c + control_c + salary_m +
#'                   control_m + s_t_ratio,
#'                   random = ~ 1 + salary_c + control_c | schoolID,
#'                   data = teachsat,
#'                   method = "REML",
#'                   control = lmeControl(opt = "optim"))
#'
#' r2mlm(model_nlme)
#'
#' @seealso Rights, J. D., & Sterba, S. K. (2019). Quantifying explained
#'   variance in multilevel models: An integrative framework for defining
#'   R-squared measures. Psychological Methods, 24(3), 309–338.
#'   <doi:10.1037/met0000184>
#'
#' @family r2mlm single model functions
#'
#' @importFrom lme4 ranef fixef VarCorr getME
#' @importFrom nlme asOneFormula
#' @importFrom magrittr %>%
#' @importFrom stats terms formula model.frame
#' @importFrom stringr str_split_fixed
#' @importFrom rlang := .data
#'
#'
#' @export

# 1 r2mlm Master Function -------------------------------------------------

r2mlm <- function(model, bargraph = TRUE) {

  # throw error if model contains higher-order terms
  temp_formula <- formula(model)
  grepl_array <- grepl("I(", temp_formula, fixed = TRUE)

  for (bool in grepl_array) {
    if (bool == TRUE) {
      stop("Error: r2mlm does not allow for models fit using the I() function; user must thus manually include any desired transformed predictor variables such as x^2 or x^3 as separate columns in dataset.")
    }
  }

  # call appropriate r2mlm helper function
  if (typeof(model) == "list") {
    r2mlm_nlme(model, bargraph)
  } else if (typeof(model) == "S4") {
    r2mlm_lmer(model, bargraph)
  } else {
    stop("You must input a model generated using either the lme4 or nlme package.")
  }

}

# 2 r2mlm_lmer helper function --------------------------------------------

r2mlm_lmer <- function(model, bargraph) {

  # Step 1) check if model has_intercept

  if ((terms(model) %>% attr("intercept")) == 1) {
    has_intercept = TRUE
  } else {
    has_intercept = FALSE
  }

  # Step 2) Pull all variable names from the formula
  all_variables <- all.vars(formula(model))
  cluster_variable <- all_variables[length(all_variables)] # pull cluster, we'll need it later

  # Step 3a) pull and prepare data

  data <- prepare_data(model, "lme4", cluster_variable)

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
    predictors <- get_interaction_vars(model)
  } else {
    predictors <- append(outcome_and_predictors[2:length(outcome_and_predictors)], get_interaction_vars(model))
  }

  # * Step 4b) Create and fill vectors
  l1_vars <- sort_variables(data, predictors, cluster_variable)$l1_vars
  l2_vars <- sort_variables(data, predictors, cluster_variable)$l2_vars

  # Step 5) pull variable names for L1 predictors with random slopes into a variable called random_slope_vars

  random_slope_vars <- get_random_slope_vars(model, has_intercept, "lme4")

  # Step 6) determine value of centeredwithincluster

  if (is.null(l1_vars)) {
    centeredwithincluster <- TRUE
  } else {
    centeredwithincluster <- get_cwc(l1_vars, cluster_variable, data)
  }

  # Step 7) pull column numbers for _covs variables
  # 7a) within_covs (l1 variables)
  # for (each value in l1_vars list) {match(value, names(data))}

  within <- get_covs(l1_vars, data)

  # 7b) pull column numbers for between_covs (l2 variables)

  between <- get_covs(l2_vars, data)

  # 7c) pull column numbers for random_covs (l1 variables with random slopes)

  random <- get_covs(random_slope_vars, data)

  # Step 8) pull gamma values (fixed slopes)
  # 8a) gamma_w, fixed slopes for L1 variables (from l1_vars list)
  gammaw <- c()
  i = 1
  for (variable in l1_vars) {
    gammaw[i] <- fixef(model)[variable]
    i = i + 1
  }

  # 8b) gamma_b, intercept value if hasintercept = TRUE, and fixed slopes for L2 variables (from between list)
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

  # Step 9) Tau matrix, results from VarCorr(model)

  vcov <- VarCorr(model)
  tau <- as.matrix(Matrix::bdiag(vcov))

  # Step 10) sigma^2 value, Rij

  sigma2 <- getME(model, "sigma")^2

  # Step 11) input everything into r2MLM

  r2mlm_manual(as.data.frame(data), within_covs = within, between_covs = between, random_covs = random, gamma_w = gammaw, gamma_b = gammab, Tau = tau, sigma2 = sigma2, has_intercept = has_intercept, clustermeancentered = centeredwithincluster, bargraph = bargraph)

}

# 3 r2mlm_nlme helper function --------------------------------------------

r2mlm_nlme <- function(model, bargraph) {

  # Step 1) check if model has_intercept

  if ((terms(model) %>% attr("intercept")) == 1) {
    has_intercept = TRUE
  } else {
    has_intercept = FALSE
  }

  # Step 2) Pull all variable names from the formula
  all_variables <- all.vars(formula(model))
  cluster_variable <- nlme::getGroups(model) %>% attr("label") # in nlme, formula(model) doesn't return grouping var, but we'll need that later on, so let's grab it here

  # Then add the grouping var to list of all variables, and calculate formula length (for later use, to iterate)
  all_variables[length(all_variables) + 1] <- cluster_variable
  formula_length <- length(all_variables) # this returns the number of elements in the all_vars list TODO remove this

  # Step 3a) pull and prepare data
  data <- prepare_data(model, "nlme", cluster_variable)

  # * Step 3b) determine whether data is appropriate format. Only the cluster variable can be a factor, for now

  outcome_and_predictors <-  all_variables[1:length(all_variables) - 1] # Pull all variables except for cluster

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
    predictors <- get_interaction_vars(model)
  } else {
    predictors <- append(outcome_and_predictors[2:length(outcome_and_predictors)], get_interaction_vars(model))
  }

  # * Step 4b) Create and fill vectors
  l1_vars <- sort_variables(data, predictors, cluster_variable)$l1_vars
  l2_vars <- sort_variables(data, predictors, cluster_variable)$l2_vars

  # Step 5) pull variable names for L1 predictors with random slopes into a variable called random_slope_vars

  random_slope_vars <- get_random_slope_vars(model, has_intercept, "nlme")

  # Step 6) determine value of centeredwithincluster

  if (is.null(l1_vars)) {
    centeredwithincluster <- TRUE
  } else {
    centeredwithincluster <- get_cwc(l1_vars, cluster_variable, data)
  }

  # Step 7) pull column numbers for _covs variables
  # 7a) within_covs (l1 variables)
  # for (each value in l1_vars list) {match(value, names(data))}
  within <- get_covs(l1_vars, data)

  # 7b) pull column numbers for between_covs (l2 variables)
  between <- get_covs(l2_vars, data)

  # 7c) pull column numbers for random_covs (l1 variables with random slopes)
  random <- get_covs(random_slope_vars, data)

  # Step 8) pull gamma values (fixed slopes)
  # 8a) gamma_w, fixed slopes for L1 variables (from l1_vars list)
  gammaw <- c()
  i = 1
  for (variable in l1_vars) {
    gammaw[i] <- nlme::fixef(model)[variable]
    i = i + 1
  }

  # 8b) gamma_b, intercept value if hasintercept = TRUE, and fixed slopes for L2 variables (from between list)
  gammab <- c()
  if (has_intercept == TRUE) {
    gammab[1] <- nlme::fixef(model)[1]
    i = 2
  } else {
    i = 1
  }
  for (variable in l2_vars) {
    gammab[i] <- nlme::fixef(model)[variable]
    i = i + 1
  }

  # Step 9) Tau matrix

  tau <- nlme::getVarCov(model)

  # Step 10) sigma^2 value, Rij

  sigma2 <- model$sigma^2

  # Step 11) input everything into r2mlm

  r2mlm_manual(as.data.frame(data), within_covs = within, between_covs = between, random_covs = random, gamma_w = gammaw, gamma_b = gammab, Tau = tau, sigma2 = sigma2, has_intercept = has_intercept, clustermeancentered = centeredwithincluster, bargraph = bargraph)

}


