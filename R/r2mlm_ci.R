#' Compute confidence intervals for R-squared in multilevel models.
#'
#' \code{r2mlm_ci} reads in a multilevel model (MLM) object generated using
#' \code{\link[lme4]{lmer}} or \code{\link[nlme]{nlme}} along with bootstrap
#' specifications and returns the upper and lower bounds of confidence intervals
#' for all R-squared measures available in the \code{r2mlm} framework.
#'
#' Note that bootstrapping confidence intervals for many R-squared values at
#' once is computationally intensive, and as a result runs somewhat slowly. For
#' this reason, the progress bar displays by default.
#'
#' @param model A model generated using \code{\link[lme4]{lmer}} or
#'   \code{\link[nlme]{nlme}}. Note that models using \code{lmer} must specify
#'   random effects at the end of the model, like so: \code{outcome ~ 1 +
#'   fixed_effects + (random_effects | cluster_variable)}. Anything else (e.g.,
#'   \code{outcome ~ 1 + (random_effects | cluster_variable) + fixed_effects})
#'   will not work.
#' @param nsim The number of bootstrapping iterations to use for the
#'   bootstrapped sampling distribution. Common values are 500 and 1000.
#' @param boottype A character vector for the type of bootstrapping to perform.
#'   Options are parametric and residual. Parametric bootstrapping assumes
#'   normally distributed residuals, whereas residual does not.
#' @param confinttype A character vector for the type of confidence interval to
#'   calculate. Options are norm (for normal), basic, and perc (for percentile).
#' @param level The desired confidence level, defaults to 0.95, yielding a 95%
#'   confidence interval.
#' @param progress TRUE/FALSE for printing progress bar or not, defaults to
#'   TRUE.
#'
#' @return If the input is a valid model, then the output will be a list of
#'   R-squared confidence intervals for all 12 measures estimated by the r2mlm
#'   function.
#'
#' @examples
#'
#' \dontrun{
#' # The "bobyqa" optimizer is required for this particular model to converge
#'
#' model_lme4 <- lmer(satisfaction ~ 1 + salary_c + control_c + salary_m + control_m +
#' s_t_ratio + (1 + salary_c + control_c| schoolID), data = teachsat, REML =
#' TRUE, control = lmerControl(optimizer = "bobyqa"))
#'
#' r2mlm_ci(model = model_lme4,
#'          nsim = 100,
#'          boottype = c("residual"),
#'          confinttype = c("perc"),
#'          level = 0.95,
#'          progress = TRUE)
#' }
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
#' @importFrom stats terms formula model.frame confint
#' @importFrom stringr str_split_fixed
#' @importFrom rlang := .data
#' @importFrom methods is
#'
#' @export

r2mlm_ci <- function(model, nsim, boottype, confinttype, level = 0.95, progress = T) {

  # function to calculate r2s
  r2s <- function(.) {
    r2mlm(., bargraph = F)$R2s[c(1:8, 10, 13, 16, 18)]
  }

  boo <- bootmlm::bootstrap_mer(x = model,
                                FUN = r2s,
                                nsim = nsim,
                                type = boottype,
                                .progress = progress)

  tryCatch(
    expr = {
      confints <- confint(boo,
                          type = confinttype,
                          level = level)
      rownames(confints) <- c("f1_total", "f2_total", "v_total", "m_total", "f_total", "fv_total", "fvm_total", "f1_within", "v_within", "fv_within", "f2_between", "m_between")
      return(confints)
    },
    error = function(e) {
      # proper dimensions for t if not centered
      boo$t0 <- boo$t0[1:5] # only 5 R2s
      boo$t <- boo$t[1:nsim, 1:5] # matrix of 5 R2 values for number of simulations
      confints <- confint(boo,
                          type = confinttype,
                          level = level)
      rownames(confints) <- c("f_total", "v_total", "m_total", "fv_total", "fvm_total")
      return(confints)
    }
  )


}

