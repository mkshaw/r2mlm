#' Compute R-squared values for longitudinal multilevel models, manually
#' inputting parameter estimates.
#'
#' \code{r2mlm_long_manual} takes as input raw data and  multilevel model (MLM)
#' parameter estimates and outputs all relevant R-squared measures as well as an
#' accompanying bar chart. This function extends the \code{r2mlm_manual}
#' function by allowing researchers to input heteroscedastic variance estimates,
#' and by providing level-specific measures for non-cluster-mean-centered
#' models.
#'
#' This function reads in raw data as well as parameter estimates from the
#' researcher’s previously fit longitudinal growth model (hence, any software
#' program can have been used to fit the researcher’s longitudinal growth model
#' prior to the use of this R function, so long as parameter estimates from the
#' fitted model are recorded; note that this function accommodates
#' non-longitudnal models as well). This function then outputs R-squared
#' measures as well as variance decompositions and associated bar charts
#' outlined in Rights & Sterba (2021). This function allows researchers to input
#' heteroscedastic residual variance by including multiple estimates, for
#' example, corresponding to individual timepoints. Users need not specify if
#' predictors are person-mean-centered or not—the function will automatically
#' output total, within-person, and between-person variance attributable to each
#' potential source of explained variance (f1, f2, v1, v2, and m). Note,
#' however, that the interpretations of these sources differ for
#' person-mean-centered versus non-person-mean-centered models and that variance
#' attributable to v2 will necessarily be 0 for person-mean-centered models.
#'
#' @param data Dataset with rows denoting observations and columns denoting
#'   variables
#' @param covs list of predictors in the dataset that have fixed components of
#'   slopes included in the model (if none, set to NULL)
#' @param random_covs list of predictors in the dataset that have random
#'   components of slopes included in the model (if none, set to NULL)
#' @param clusterID variable name in dataset corresponding to cluster (e.g.,
#'   person) identification
#' @param gammas vector containing estimated fixed components of all slopes,
#'   listed in the order specified in covs (if none, set to NULL)
#' @param Tau random effect covariance matrix; the first row and the first
#'   column denote the intercept variance and covariances and each subsequent
#'   row/column denotes a given random slope’s variance and covariances (to be
#'   entered in the order listed by random_covs)
#' @param sigma2 level-1 residual variance; can be entered as a single number,
#'   or as a set of numbers, for example corresponding to different residual
#'   variances at individual timepoints; if entered as a set of numbers,
#'   function will assume equal weights and take the raw average of these to
#'   estimate the expectation of the error variance
#' @param bargraph Optional bar graph output, default is TRUE.
#'
#' @return If the input is valid, then the output will be a list and associated
#'   graphical representation of R-squared decompositions. If the input is not
#'   valid, it will return an error.
#'
#' @seealso Rights, J. D., & Sterba, S. K. (2021). Effect size measures for
#'   longitudinal growth analyses: Extending a framework of multilevel model
#'   R-squareds to accommodate heteroscedasticity, autocorrelation,
#'   nonlinearity, and alternative centering strategies. New Directions for
#'   Child and Adolescent Development, 2021, 65– 110. <doi:10.1002/cad.20387>
#'
#' @family r2mlm single model functions
#'
#' @importFrom rockchalk gmc

r2mlm_long_manual <- function(data, covs, random_covs, clusterID,
                       gammas, Tau, sigma2, bargraph = TRUE) {

    if (is.null(covs) == FALSE) {
      centered_data <- gmc(data, covs, clusterID)
      phi_w <- var(centered_data[, c(paste0(covs, "_dev"))])
      phi_b <- var(centered_data[, c(paste0(covs, "_mn"))])
      gammas <- matrix(c(gammas), ncol = 1)
      f1 <- t(gammas) %*% phi_w %*% gammas
      f2 <- t(gammas) %*% phi_b %*% gammas
    }
    else{
      f1 <- 0
      f2 <- 0
    }

    if (is.null(random_covs) == FALSE) {
      centered_data_rand <- gmc(data, random_covs, clusterID)
      Sig_w <- var(centered_data_rand[, c(paste0(random_covs, "_dev"))])
      Sig_b <- var(centered_data_rand[, c(paste0(random_covs, "_mn"))])
      m_mat <-
        matrix(c(colMeans(cbind(1, data[, c(random_covs)]))), ncol = 1)
      v1 <- sum(diag(Tau[2:nrow(Tau), 2:nrow(Tau)] %*% Sig_w))
      v2 <- sum(diag(Tau[2:nrow(Tau), 2:nrow(Tau)] %*% Sig_b))
    }
    else{
      v1 <- 0
      v2 <- 0
      m_mat <- 1
    }

    m <- t(m_mat) %*% Tau %*% m_mat

    sigma <- mean(sigma2)

    #decompositions

    decomp_fixed_within <- f1 / sum(f1, f2, v1, v2, m, sigma)
    decomp_fixed_between <- f2 / sum(f1, f2, v1, v2, m, sigma)
    decomp_varslopes_within <- v1 / sum(f1, f2, v1, v2, m, sigma)
    decomp_varslopes_between <- v2 / sum(f1, f2, v1, v2, m, sigma)
    decomp_varmeans <- m / sum(f1, f2, v1, v2, m, sigma)
    decomp_sigma <- sigma / sum(f1, f2, v1, v2, m, sigma)

    decomp_fixed_within_w <- f1 / sum(f1, v1, sigma)
    decomp_fixed_between_b <- f2 / sum(f2, v2, m)
    decomp_varslopes_within_w <- v1 / sum(f1, v1, sigma)
    decomp_varslopes_between_b <- v2 / sum(f2, v2, m)
    decomp_varmeans_b <- m / sum(f2, v2, m)
    decomp_sigma_w <- sigma / sum(f1, v1, sigma)

    #barchart
    if (bargraph == TRUE) {

      contributions_stacked <-
        matrix(
          c(
            decomp_fixed_within, decomp_fixed_between, decomp_varslopes_within,
            decomp_varslopes_between, decomp_varmeans, decomp_sigma,
            decomp_fixed_within_w, 0, decomp_varslopes_within_w,
            0, 0, decomp_sigma_w,
            0, decomp_fixed_between_b, 0,
            decomp_varslopes_between_b, decomp_varmeans_b, 0
          ),
          6, 3
        )
      colnames(contributions_stacked) <- c("total", "within", "between")
      rownames(contributions_stacked) <- c(
        "fixed slopes (within)",
        "fixed slopes (between)",
        "slope variation (within)",
        "slope variation (between)",
        "intercept variation (between)",
        "residual (within)"
      )


      barplot(
        contributions_stacked,
        main = "Decomposition",
        horiz = FALSE,
        ylim = c(0, 1),
        col = c(
          "darkred",
          "steelblue",
          "darkred",
          "steelblue",
          "midnightblue",
          "white"
        ),
        ylab = "proportion of variance",
        density = c(NA, NA, 30, 40, 40, NA),
        angle = c(0, 45, 0, 90, 135, 0),
        xlim = c(0, 1.5),
        width = c(.3, .3)
      )
      legend(
        1.1,
        .65,
        legend = rownames(contributions_stacked),
        fill = c(
          "darkred",
          "steelblue",
          "darkred",
          "steelblue",
          "midnightblue",
          "white"
        ),
        cex = .7,
        pt.cex = 1,
        xpd = T,
        density = c(NA, NA, 30, 40, 40, NA),
        angle = c(0, 45, 0, 90, 135, 0)
      )
    }

    #create tables for output

    decomp_table <-
      matrix(
        c(
          decomp_fixed_within, decomp_fixed_between, decomp_varslopes_within,
          decomp_varslopes_between, decomp_varmeans, decomp_sigma,
          decomp_fixed_within_w, "NA", decomp_varslopes_within_w,
          "NA", "NA", decomp_sigma_w,
          "NA", decomp_fixed_between_b, "NA",
          decomp_varslopes_between_b, decomp_varmeans_b, "NA"
        ),
        6,
        3
      )
    colnames(decomp_table) <- c("total", "within", "between")
    rownames(decomp_table) <- c(
      "fixed slopes (within)",
      "fixed slopes (between)",
      "slope variation (within)",
      "slope variation (between)",
      "intercept variation (between)",
      "residual (within)"
    )

    R2_table <-
      matrix(
        c(
          decomp_fixed_within, decomp_fixed_between, decomp_varslopes_within,
          decomp_varslopes_between, decomp_varmeans, decomp_fixed_within + decomp_fixed_between,
          decomp_fixed_within + decomp_fixed_between + decomp_varslopes_within + decomp_varslopes_between, decomp_fixed_within + decomp_fixed_between + decomp_varslopes_within + decomp_varslopes_between + decomp_varmeans, decomp_fixed_within_w,
          "NA", decomp_varslopes_within_w, "NA",
          "NA", "NA", decomp_fixed_within_w + decomp_varslopes_within_w,
          "NA", "NA", decomp_fixed_between_b,
          "NA", decomp_varslopes_between_b, decomp_varmeans_b,
          "NA", decomp_fixed_between_b + decomp_varslopes_between_b, "NA"
        ),
        8,
        3
      )

    colnames(R2_table) <- c("total", "within", "between")
    rownames(R2_table) <- c("f1", "f2", "v1", "v2", "m", "f", "fv", "fvm")

    Output <- list(noquote(decomp_table), noquote(R2_table))
    names(Output) <- c("Decompositions", "R2s")

    return(Output)
  }
