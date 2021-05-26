#' Compute R-squared values for three-level multilevel models, manually
#' inputting parameter estimates.
#'
#' \code{r2mlm3_manual} takes as input raw data and three-level multilevel model
#' (MLM) parameter estimates and outputs all relevant R-squared measures as well
#' as an accompanying bar chart.
#'
#' This function can also accommodate two-level models. To input results
#' for two-level models, set the following arguments equal to NULL: l3_covs,
#' random_covs13, random_covs23, gamma_3, Tau13, Tau23.
#'
#' @param data Dataset with rows denoting observations and columns denoting
#'   variables
#' @param l1_covs Vector of numbers (or variable names) corresponding to the
#'   columns in the dataset of the level-1 predictors used in the MLM (if none
#'   used, set to NULL)
#' @param l2_covs Vector of numbers (or variable names) corresponding to the
#'   columns in the dataset of the level-2 predictors used in the MLM (if none
#'   used, set to NULL)
#' @param l3_covs Vector of numbers (or variable names) corresponding to the
#'   columns in the dataset of the level-3 predictors used in the MLM (if none
#'   used, set to NULL)
#' @param random_covs12 Vector of numbers (or variable names) corresponding to
#'   the columns in the dataset of the level-1 predictors that have random
#'   slopes across level-2 units in the MLM (if no such random slopes, set to
#'   NULL)
#' @param random_covs13 Vector of numbers (or variable names) corresponding to
#'   the columns in the dataset of the level-1 predictors that have random
#'   slopes across level-3 units in the MLM (if no such random slopes, set to
#'   NULL)
#' @param random_covs23 Vector of numbers (or variable names) corresponding to
#'   the columns in the dataset of the level-2 predictors that have random
#'   slopes across level-3 units in the MLM (if no such random slopes, set to
#'   NULL)
#' @param gamma_1 Vector of fixed slope estimates for all level-1 predictors, to
#'   be entered in the order of the predictors listed by l1_covs (if none, set
#'   to NULL)
#' @param gamma_2 Vector of fixed slope estimates for all level-2 predictors, to
#'   be entered in the order of the predictors listed by l2_covs (if none, set
#'   to NULL)
#' @param gamma_3 Vector of fixed slope estimates for all level-3 predictors, to
#'   be entered in the order of the predictors listed by l3_covs (if none, set
#'   to NULL)
#' @param Tau12 For cluster-mean-centered model results (set to NULL if entering
#'   non-cluster-mean-centered model results), this is the random effect
#'   covariance matrix with the first row/column denoting the intercept variance
#'   and covariances across level-2 units and each subsequent row/column denotes
#'   a given level-1 predictor’s random slope variance and covariances across
#'   level-2 units (to be entered in the order listed by random_covs12; if none,
#'   set to NULL)
#' @param Tau13 For cluster-mean-centered model results (set to NULL if entering
#'   non-cluster-mean-centered model results), this is the random effect
#'   covariance matrix with the first row/column denoting the intercept variance
#'   and covariances across level-3 units and each subsequent row/column denotes
#'   a given level-1 predictor’s random slope variance and covariances across
#'   level-3 units (to be entered in the order listed by random_covs13; if none,
#'   set to NULL)
#' @param Tau23 For cluster-mean-centered model results (set to NULL if entering
#'   non-cluster-mean-centered model results), this is the random effect
#'   covariance matrix with each row/column denoting a given level-2 predictor’s
#'   random slope variance and covariances across level-3 units (to be entered
#'   in the order listed by random_covs23; if none, set to NULL)
#' @param sigma2 Level-1 residual variance
#' @param clustermeancentered By default, this argument is set to TRUE,
#'   indicating that cluster-mean-centered model results are being inputted.
#'   When instead entering non-cluster-mean-centered model results, set this
#'   argument to FALSE. Additionally, for non-cluster-mean-centered model
#'   results, random effect variances/covariances are to be entered in arguments
#'   Tau2_noncmc and Tau3_noncmc (defined below), rather than in the Tau12,
#'   Tau13, and Tau23 arguments used for cluster-mean-centered model results.
#'   Additionally, when entering non-cluster-mean-centered model results, user
#'   must specify l2clusterID_noncmc and l3clusterID_noncmc (neither of which
#'   are necessary for cluster-mean-centered model results). Function input is
#'   otherwise the same for cluster-mean-centered and non-cluster-mean-centered
#'   model results.
#' @param Tau2_noncmc For non-cluster-mean-centered model results, this is the
#'   level-2 random effect covariance matrix; the first row/column denotes the
#'   intercept variance and covariances across level-2 units and each subsequent
#'   row/column denotes a given predictor’s random slope variance and
#'   covariances across level-2 units (to be entered in the order listed by
#'   randomcovsl2_noncmc; by default, this argument is set to NULL)
#' @param Tau3_noncmc For non-cluster-mean-centered model results, this is the
#'   level-3 random effect covariance matrix; the first row/column denotes the
#'   intercept variance and covariances across level-3 units and each subsequent
#'   row/column denotes a given predictor’s random slope variance and
#'   covariances across level-3 units (to be entered in the order listed by
#'   randomcovsl2_noncmc; by default, this argument is set to NULL)
#' @param l2clusterID_noncmc For non-cluster-mean-centered model results, this
#'   is the number (or variable name) corresponding to the column in the dataset
#'   containing the level-2 cluster identification (function assumes that each
#'   level-2 cluster ID is unique; by default, this argument is set to NULL)
#' @param l3clusterID_noncmc For non-cluster-mean-centered model results, this
#'   is the number (or variable name) corresponding to the column in the dataset
#'   containing the level-3 cluster identification (function assumes that each
#'   level-3 cluster ID is unique; by default, this argument is set to NULL)
#' @param bargraph Optional bar graph output, default is TRUE.
#'
#' @return If the input is valid, then the output will be a list and associated
#'   graphical representation of R-squared decompositions. If the input is not
#'   valid, it will return an error.
#'
#' @family r2mlm single model functions
#'
#' @importFrom stats var
#' @importFrom graphics barplot legend
#'
#' @export

r2mlm3_manual <-
  function(data, l1_covs, l2_covs, l3_covs, random_covs12, random_covs13,
           random_covs23, gamma_1, gamma_2, gamma_3, Tau12, Tau13, Tau23,
           sigma2, clustermeancentered = TRUE, Tau2_noncmc = NULL,
           Tau3_noncmc = NULL, l2clusterID_noncmc = NULL, l3clusterID_noncmc = NULL,
           bargraph = TRUE) {

    if (clustermeancentered == TRUE) {
      ##compute phis
      if (is.null(l1_covs) == T) {
        phi_1 <- 0
        gamma_1 <- 0
      } else{
        phi_1 <- var(data[, l1_covs], na.rm = T)
      }
      if (is.null(l2_covs) == T) {
        phi_2 <- 0
        gamma_2 <- 0
      } else{
        phi_2 <- var(data[, l2_covs], na.rm = T)
      }
      if (is.null(l3_covs) == T) {
        phi_3 <- 0
        gamma_3 <- 0
      } else{
        phi_3 <- var(data[, l3_covs], na.rm = T)
      }

      ##compute variances/covariances of random components

      if (is.null(Tau13) == TRUE)
        Tau13 <- matrix(c(0), 1, 1)

      if (is.null(random_covs12) == F) {
        var_randomcovs12 <- var(cbind(1, data[, c(random_covs12)]), na.rm = T)
        psi12 <- matrix(c(diag(Tau12)), ncol = 1)
        kappa12 <- matrix(c(Tau12[lower.tri(Tau12) == TRUE]), ncol = 1)
        v12 <- matrix(c(diag(var_randomcovs12)), ncol = 1)
        r12 <-
          matrix(c(var_randomcovs12[lower.tri(var_randomcovs12) == TRUE]), ncol =
                   1)
      } else{
        var_randomcovs12 <- 0
        psi12 <- 0
        kappa12 <- 0
        v12 <- 0
        r12 <- 0
      }

      if (is.null(random_covs13) == F) {
        var_randomcovs13 <- var(cbind(1, data[, c(random_covs13)]), na.rm = T)
        psi13 <- matrix(c(diag(Tau13)), ncol = 1)
        kappa13 <- matrix(c(Tau13[lower.tri(Tau13) == TRUE]), ncol = 1)
        v13 <- matrix(c(diag(var_randomcovs13)), ncol = 1)
        r13 <-
          matrix(c(var_randomcovs13[lower.tri(var_randomcovs13) == TRUE]), ncol =
                   1)
      } else{
        var_randomcovs13 <- 0
        psi13 <- 0
        kappa13 <- 0
        v13 <- 0
        r13 <- 0
      }

      if (is.null(random_covs23) == F) {
        var_randomcovs23 <- var(cbind(data[, c(random_covs23)]), na.rm = T)
        psi23 <- matrix(c(diag(Tau23)), ncol = 1)
        kappa23 <- matrix(c(Tau23[lower.tri(Tau23) == TRUE]), ncol = 1)
        v23 <- matrix(c(diag(var_randomcovs23)), ncol = 1)
        r23 <-
          matrix(c(var_randomcovs23[lower.tri(var_randomcovs23) == TRUE]), ncol =
                   1)
      } else{
        var_randomcovs23 <- 0
        psi23 <- 0
        kappa23 <- 0
        v23 <- 0
        r23 <- 0
      }

      ##compute level-specific and total variance

      l1var <-
        (t(gamma_1) %*% phi_1 %*% gamma_1) + (t(v12) %*% psi12 + 2 * (t(r12) %*%
                                                                        kappa12)) + (t(v13) %*% psi13 + 2 * (t(r13) %*% kappa13)) + sigma2

      l2var <-
        (t(gamma_2) %*% phi_2 %*% gamma_2) + (t(v23) %*% psi23 + 2 * (t(r23) %*% kappa23)) + Tau12[1, 1]

      l3var <- (t(gamma_3) %*% phi_3 %*% gamma_3) + Tau13[1, 1]

      totalvar <- l1var + l2var + l3var

      #total R2 measures

      R2_f1_t <- (t(gamma_1) %*% phi_1 %*% gamma_1) / totalvar
      R2_f2_t <- (t(gamma_2) %*% phi_2 %*% gamma_2) / totalvar
      R2_f3_t <- (t(gamma_3) %*% phi_3 %*% gamma_3) / totalvar
      R2_v12_t <- (t(v12) %*% psi12 + 2 * (t(r12) %*% kappa12)) / totalvar
      R2_v13_t <- (t(v13) %*% psi13 + 2 * (t(r13) %*% kappa13)) / totalvar
      R2_v23_t <- (t(v23) %*% psi23 + 2 * (t(r23) %*% kappa23)) / totalvar
      R2_m2_t <- Tau12[1, 1] / totalvar
      R2_m3_t <- Tau13[1, 1] / totalvar

      #l1 R2 measures
      R2_f1_1 <- (t(gamma_1) %*% phi_1 %*% gamma_1) / l1var
      R2_v12_1 <- (t(v12) %*% psi12 + 2 * (t(r12) %*% kappa12)) / l1var
      R2_v13_1 <- (t(v13) %*% psi13 + 2 * (t(r13) %*% kappa13)) / l1var

      #l2 R2 measures
      R2_f2_2 <- (t(gamma_2) %*% phi_2 %*% gamma_2) / l2var
      R2_v23_2 <- (t(v23) %*% psi23 + 2 * (t(r23) %*% kappa23)) / l2var
      R2_m2_2 <- Tau12[1, 1] / l2var

      #l3 R2 measures
      R2_f3_3 <-  (t(gamma_3) %*% phi_3 %*% gamma_3) / l3var
      R2_m3_3 <- Tau13[1, 1] / l3var

      R2_table <-
        matrix(
          c(
            R2_f1_t, R2_f2_t, R2_f3_t, R2_v12_t,
            R2_v13_t, R2_v23_t, R2_m2_t, R2_m3_t,
            R2_f1_1, "NA", "NA", R2_v12_1,
            R2_v13_1, "NA", "NA", "NA",
            "NA", R2_f2_2, "NA", "NA",
            "NA", R2_v23_2, R2_m2_2, "NA",
            "NA", "NA", R2_f3_3, "NA",
            "NA", "NA", "NA", R2_m3_3
          )
          ,
          ncol = 4
        )
      rownames(R2_table) <-
        c("f1", "f2", "f3", "v12", "v13", "v23", "m2", "m3")
      colnames(R2_table) <- c("total", "l1", "l2", "l3")

      ##barchart
      if (bargraph == TRUE) {
        contributions_stacked <-
          matrix(
            c(R2_f1_t, R2_f2_t, R2_f3_t, R2_v12_t,
              R2_v13_t, R2_v23_t, R2_m2_t, R2_m3_t,
              sigma2 / totalvar, R2_f1_1, 0, 0,
              R2_v12_1, R2_v13_1, 0, 0,
              0, sigma2 / l1var, 0, R2_f2_2,
              0, 0, 0, R2_v23_2,
              R2_m2_2, 0, 0, 0,
              0, R2_f3_3, 0, 0,
              0, 0, R2_m3_3, 0
            ),
            nrow = 9,
            ncol = 4
          )
        colnames(contributions_stacked) <-
          c("total", "level-1", "level-2", "level-3")
        rownames(contributions_stacked) <- c("f1", "f2", "f3",
                                             "v12", "v13", "v23",
                                             "m2", "m3", "resid")

        barplot(
          contributions_stacked,
          main = "Decomposition",
          horiz = FALSE,
          ylim = c(0, 1),
          col = c(
            "darkred",
            "steelblue",
            "darkgoldenrod1",
            "darkred",
            "darkred",
            "midnightblue",
            "midnightblue",
            "darkgoldenrod1",
            "white"
          ),
          ylab = "proportion of variance",
          density = c(NA, NA, NA, 20, 40, 20, 40, 20, NA),
          angle = c(0, 0, 0, 45, 45, 135, 135, 0, 0)
        )
        legend(
          4.89,
          .7,
          title = "source",
          legend = rownames(contributions_stacked),
          fill = c(
            "darkred",
            "steelblue",
            "darkgoldenrod1",
            "darkred",
            "darkred",
            "midnightblue",
            "midnightblue",
            "darkgoldenrod1",
            "white"
          ),
          cex = .7,
          pt.cex = 1,
          xpd = TRUE,
          density = c(NA, NA, NA, 30, 50, 30, 50, 30, NA),
          angle = c(0, 0, 0, 45, 45, 135, 135, 0, 0)
        )
      }

      Output <- list(noquote(R2_table))
      names(Output) <- c("R2s")
    }

    if (clustermeancentered == FALSE) {
      allfixedcovs_noncmc <- c(l1_covs, l2_covs, l3_covs)
      allgamma_noncmc <- c(gamma_1, gamma_2, gamma_3)
      randomcovsl2_noncmc <- random_covs12
      randomcovsl3_noncmc <- c(random_covs13, random_covs23)

      ##create version of dataset that decomposes predictors into level-specific portions

      data_temp <- data

      if (is.null(allfixedcovs_noncmc) == F) {
        data_fixed_centered <-
          matrix(NA, nrow(data), length(allfixedcovs_noncmc) * 3)
        colnames(data_fixed_centered) <-
          c(paste0("p", seq(length(
            allfixedcovs_noncmc
          )), "_l1"),
          paste0("p", seq(length(
            allfixedcovs_noncmc
          )), "_l2"),
          paste0("p", seq(length(
            allfixedcovs_noncmc
          )), "_l3"))


        if (length(l2_covs) > 0) {
          data_fixed_centered[, (length(l1_covs) + 1):(length(allfixedcovs_noncmc))] <-
            0
        }

        if (length(l3_covs) > 0) {
          data_fixed_centered[, (length(allfixedcovs_noncmc) + length(l1_covs) + length(l2_covs) +
                                   1):(2 * length(allfixedcovs_noncmc))] <- 0
          data_fixed_centered[, (length(l1_covs) + 1):(length(allfixedcovs_noncmc))] <-
            0
        }
      }

      if (is.null(randomcovsl2_noncmc) == F) {
        data_randoml2_centered <-
          matrix(NA, nrow(data), length(randomcovsl2_noncmc) * 3)
        colnames(data_randoml2_centered) <-
          c(paste0("p", seq(length(
            randomcovsl2_noncmc
          )), "_l1"),
          paste0("p", seq(length(
            randomcovsl2_noncmc
          )), "_l2"),
          paste0("p", seq(length(
            randomcovsl2_noncmc
          )), "_l3"))
      }

      if (is.null(randomcovsl3_noncmc) == F) {
        data_randoml3_centered <-
          matrix(NA, nrow(data), length(randomcovsl3_noncmc) * 3)
        colnames(data_randoml3_centered) <-
          c(paste0("p", seq(length(
            randomcovsl3_noncmc
          )), "_l1"),
          paste0("p", seq(length(
            randomcovsl3_noncmc
          )), "_l2"),
          paste0("p", seq(length(
            randomcovsl3_noncmc
          )), "_l3"))
      }


      ##loop through level-2 clusters
      for (clus2 in seq(max(as.numeric(data[, l2clusterID_noncmc])))) {
        for (i in seq(nrow(data))) {
          if (data[i, l2clusterID_noncmc] == clus2) {
            if (is.null(allfixedcovs_noncmc) == F) {
              for (ncov in seq(length(allfixedcovs_noncmc))) {
                ##compute level-1 portion of predictor
                data_fixed_centered[i, paste0("p", ncov, "_l1")] <-
                  data[i, allfixedcovs_noncmc[ncov]] - mean(data_temp[which(data_temp[, l2clusterID_noncmc] ==
                                                                              clus2), allfixedcovs_noncmc[ncov]], na.rm = T)
                ##compute level-2 portion of predictor
                data_fixed_centered[i, paste0("p", ncov, "_l2")] <-
                  mean(data_temp[which(data_temp[, l2clusterID_noncmc] == clus2), allfixedcovs_noncmc[ncov]], na.rm =
                         T)
              }
            }

            if (is.null(randomcovsl2_noncmc) == F) {
              for (ncov in seq(length(randomcovsl2_noncmc))) {
                ##compute level-1 portion of predictor
                data_randoml2_centered[i, paste0("p", ncov, "_l1")] <-
                  data[i, randomcovsl2_noncmc[ncov]] - mean(data_temp[which(data_temp[, l2clusterID_noncmc] ==
                                                                              clus2), randomcovsl2_noncmc[ncov]], na.rm = T)
                ##compute level-2 portion of predictor
                data_randoml2_centered[i, paste0("p", ncov, "_l2")] <-
                  mean(data_temp[which(data_temp[, l2clusterID_noncmc] == clus2), randomcovsl2_noncmc[ncov]], na.rm =
                         T)
              }
            }

            if (is.null(randomcovsl3_noncmc) == F) {
              for (ncov in seq(length(randomcovsl3_noncmc))) {
                ##compute level-1 portion of predictor
                data_randoml3_centered[i, paste0("p", ncov, "_l1")] <-
                  data[i, randomcovsl3_noncmc[ncov]] - mean(data_temp[which(data_temp[, l2clusterID_noncmc] ==
                                                                              clus2), randomcovsl3_noncmc[ncov]], na.rm = T)
                ##compute level-2 portion of predictor
                data_randoml3_centered[i, paste0("p", ncov, "_l2")] <-
                  mean(data_temp[which(data_temp[, l2clusterID_noncmc] == clus2), randomcovsl3_noncmc[ncov]], na.rm =
                         T)
              }
            }
          }
        }
      }

      ##loop through level-3 clusters
      for (clus3 in seq(max(as.numeric(data[, l3clusterID_noncmc])))) {
        for (i in seq(nrow(data))) {
          if (data[i, "schoolid"] == clus3) {
            if (is.null(allfixedcovs_noncmc) == F) {
              for (ncov in seq(length(allfixedcovs_noncmc))) {
                ##compute level-3 portion of predictor
                data_fixed_centered[i, paste0("p", ncov, "_l3")] <-
                  mean(data_temp[which(data_temp[, l3clusterID_noncmc] == clus3), allfixedcovs_noncmc[ncov]], na.rm =
                         T)
              }
            }

            if (is.null(randomcovsl2_noncmc) == F) {
              for (ncov in seq(length(randomcovsl2_noncmc))) {
                ##compute level-3 portion of predictor
                data_randoml2_centered[i, paste0("p", ncov, "_l3")] <-
                  mean(data_temp[which(data_temp[, l3clusterID_noncmc] == clus3), randomcovsl2_noncmc[ncov]], na.rm =
                         T)
              }
            }

            if (is.null(randomcovsl3_noncmc) == F) {
              for (ncov in seq(length(randomcovsl3_noncmc))) {
                ##compute level-3 portion of predictor
                data_randoml3_centered[i, paste0("p", ncov, "_l3")] <-
                  mean(data_temp[which(data_temp[, l3clusterID_noncmc] == clus3), randomcovsl3_noncmc[ncov]], na.rm =
                         T)
              }
            }

          }
        }
      }

      if (is.null(allfixedcovs_noncmc) == F) {
        for (ncov in seq(length(allfixedcovs_noncmc))) {
          ##cluster-mean-center level-2 predictor
          data_fixed_centered[, paste0("p", ncov, "_l2")] <-
            data_fixed_centered[, paste0("p", ncov, "_l2")] - data_fixed_centered[, paste0("p", ncov, "_l3")]
        }
      }

      if (is.null(randomcovsl2_noncmc) == F) {
        for (ncov in seq(length(randomcovsl2_noncmc))) {
          ##cluster-mean-center level-2 predictor
          data_randoml2_centered[, paste0("p", ncov, "_l2")] <-
            data_randoml2_centered[, paste0("p", ncov, "_l2")] - data_randoml2_centered[, paste0("p", ncov, "_l3")]
        }
      }

      if (is.null(randomcovsl3_noncmc) == F) {
        for (ncov in seq(length(randomcovsl3_noncmc))) {
          ##cluster-mean-center level-2 predictor
          data_randoml3_centered[, paste0("p", ncov, "_l2")] <-
            data_randoml3_centered[, paste0("p", ncov, "_l2")] - data_randoml3_centered[, paste0("p", ncov, "_l3")]
        }
      }

      ##make lower-level portions of cluster-level predictors equal to 0
      if (is.null(allfixedcovs_noncmc) == F) {
        if (length(l2_covs) > 0) {
          data_fixed_centered[, (length(l1_covs) + 1):(length(allfixedcovs_noncmc))] <-
            0
        }

        if (length(l3_covs) > 0) {
          data_fixed_centered[, (length(allfixedcovs_noncmc) + length(l1_covs) + length(l2_covs) +
                                   1):(2 * length(allfixedcovs_noncmc))] <- 0
          data_fixed_centered[, (length(l1_covs) + 1):(length(allfixedcovs_noncmc))] <-
            0
        }
      }


      if (is.null(allfixedcovs_noncmc) == F) {
        ##compute phis
        phi1 <-
          var(data_fixed_centered[, seq(length(allfixedcovs_noncmc))], na.rm = T)
        phi2 <-
          var(data_fixed_centered[, c(length(allfixedcovs_noncmc) + seq(length(allfixedcovs_noncmc)))], na.rm =
                T)
        phi3 <-
          var(data_fixed_centered[, c(2 * length(allfixedcovs_noncmc) + seq(length(allfixedcovs_noncmc)))], na.rm =
                T)

        ##compute variance explained by fs
        f1 <- t(allgamma_noncmc) %*% phi1 %*% allgamma_noncmc
        f2 <- t(allgamma_noncmc) %*% phi2 %*% allgamma_noncmc
        f3 <- t(allgamma_noncmc) %*% phi3 %*% allgamma_noncmc
      } else{
        f1 <- f2 <- f3 <- 0
      }

      ##make lower-level portions of cluster-level predictors equal to 0
      if (is.null(randomcovsl3_noncmc) == F) {
        if (length(random_covs23) > 0) {
          data_randoml3_centered[, (length(random_covs13) + 1):(length(randomcovsl3_noncmc))] <-
            0
        }
      }

      if (is.null(randomcovsl2_noncmc) == F) {
        ##compute Sigmas
        Sigma12 <-
          var(cbind(0, data_randoml2_centered[, seq(length(randomcovsl2_noncmc))]), na.rm =
                T)
        Sigma22 <-
          var(cbind(0, data_randoml2_centered[, c(length(randomcovsl2_noncmc) + seq(length(randomcovsl2_noncmc)))]), na.rm =
                T)
        Sigma32 <-
          var(cbind(0, data_randoml2_centered[, c(2 * length(randomcovsl2_noncmc) +
                                                    seq(length(randomcovsl2_noncmc)))]), na.rm = T)

        ##compute variance explained by vs
        v12 <- as.numeric(sum(diag(Sigma12 %*% Tau2_noncmc)))
        v22 <- as.numeric(sum(diag(Sigma22 %*% Tau2_noncmc)))
        v32 <- as.numeric(sum(diag(Sigma32 %*% Tau2_noncmc)))

        ##compute m vectors

        m_mat12 <-
          matrix(c(as.numeric(colMeans(
            cbind(0, data_randoml2_centered[, seq(length(randomcovsl2_noncmc))]), na.rm =
              T
          ))), ncol = 1)
        m_mat22 <-
          matrix(c(as.numeric(colMeans(
            cbind(0, data_randoml2_centered[, c(length(randomcovsl2_noncmc) + seq(length(randomcovsl2_noncmc)))]), na.rm =
              T
          ))), ncol = 1)
        m_mat32 <-
          matrix(c(as.numeric(colMeans(
            cbind(1, data_randoml2_centered[, c(2 * length(randomcovsl2_noncmc) + seq(length(randomcovsl2_noncmc)))]), na.rm =
              T
          ))), ncol = 1)

        ##compute variance explained by m

        m2 <-
          (
            t(m_mat12) %*% Tau2_noncmc %*% m_mat12 + t(m_mat22) %*% Tau2_noncmc %*% m_mat22 +
              t(m_mat32) %*% Tau2_noncmc %*% m_mat32
          )

      } else{
        v12 <- v22 <- v32 <- 0
        m2 <- Tau2_noncmc[1, 1]
      }

      if (is.null(randomcovsl3_noncmc) == F) {
        ##compute Sigmas
        Sigma13 <-
          var(cbind(0, data_randoml3_centered[, seq(length(randomcovsl3_noncmc))]), na.rm =
                T)
        Sigma23 <-
          var(cbind(0, data_randoml3_centered[, c(length(randomcovsl3_noncmc) + seq(length(randomcovsl3_noncmc)))]), na.rm =
                T)
        Sigma33 <-
          var(cbind(0, data_randoml3_centered[, c(2 * length(randomcovsl3_noncmc) +
                                                    seq(length(randomcovsl3_noncmc)))]), na.rm = T)

        ##compute variance explained by vs
        v13 <- as.numeric(sum(diag(Sigma13 %*% Tau3_noncmc)))
        v23 <- as.numeric(sum(diag(Sigma23 %*% Tau3_noncmc)))
        v33 <- as.numeric(sum(diag(Sigma33 %*% Tau3_noncmc)))

        ##compute m vectors
        m_mat13 <-
          matrix(c(as.numeric(colMeans(
            cbind(0, data_randoml3_centered[, seq(length(randomcovsl3_noncmc))]), na.rm =
              T
          ))), ncol = 1)
        m_mat23 <-
          matrix(c(as.numeric(colMeans(
            cbind(0, data_randoml3_centered[, c(length(randomcovsl3_noncmc) + seq(length(randomcovsl3_noncmc)))]), na.rm =
              T
          ))), ncol = 1)
        m_mat33 <-
          matrix(c(as.numeric(colMeans(
            cbind(1, data_randoml3_centered[, c(2 * length(randomcovsl3_noncmc) + seq(length(randomcovsl3_noncmc)))]), na.rm =
              T
          ))), ncol = 1)

        ##compute variance explained by m
        m3 <-
          (
            t(m_mat13) %*% Tau3_noncmc %*% m_mat13 + t(m_mat23) %*% Tau3_noncmc %*% m_mat23 +
              t(m_mat33) %*% Tau3_noncmc %*% m_mat33
          )

      } else{
        v13 <- v23 <- v33 <- 0
        m3 <- Tau3_noncmc[1, 1]
      }

      ##compute R2 measures

      totvar <- sum(f1, f2, f3, v12, v13, v22, v23, v32, v33, m2, m3, sigma2)
      l1var <- sum(f1, v12, v13, sigma2)
      l2var <- sum(f2, v22, v23, m2)
      l3var <- sum(f3, v32, v33, m3)

      R2_f1_t <- f1 / totvar
      R2_f2_t <- f2 / totvar
      R2_f3_t <- f3 / totvar
      R2_v12_t <- v12 / totvar
      R2_v13_t <- v13 / totvar
      R2_v22_t <- v22 / totvar
      R2_v23_t <- v23 / totvar
      R2_v32_t <- v32 / totvar
      R2_v33_t <- v33 / totvar
      R2_m2_t <- m2 / totvar
      R2_m3_t <- m3 / totvar

      R2_f1_1 <- f1 / l1var
      R2_f2_2 <- f2 / l2var
      R2_f3_3 <- f3 / l3var
      R2_v12_1 <- v12 / l1var
      R2_v13_1 <- v13 / l1var
      R2_v22_2 <- v22 / l2var
      R2_v23_2 <- v23 / l2var
      R2_v32_3 <- v32 / l3var
      R2_v33_3 <- v33 / l3var
      R2_m2_2 <- m2 / l2var
      R2_m3_3 <- m3 / l3var

      ##create table of measures

      R2_table <-
        matrix(
          c(
            R2_f1_t, R2_f2_t, R2_f3_t, R2_v12_t,
            R2_v13_t, R2_v22_t, R2_v23_t, R2_v32_t,
            R2_v33_t, R2_m2_t, R2_m3_t, R2_f1_1,
            "NA", "NA", R2_v12_1, R2_v13_1,
            "NA", "NA", "NA", "NA",
            "NA", "NA", "NA", R2_f2_2,
            "NA", "NA", "NA", R2_v22_2,
            R2_v23_2, "NA", "NA", R2_m2_2,
            "NA", "NA", "NA", R2_f3_3,
            "NA", "NA", "NA", "NA",
            R2_v32_3, R2_v33_3, "NA", R2_m3_3
          ),
          ncol = 4
        )
      rownames(R2_table) <-
        c("f1",
          "f2",
          "f3",
          "v12",
          "v13",
          "v22",
          "v23",
          "v32",
          "v33",
          "m2",
          "m3")
      colnames(R2_table) <- c("total", "l1", "l2", "l3")

      ##barchart

      contributions_stacked <-
        matrix(
          c(
            R2_f1_t, R2_f2_t, R2_f3_t, R2_v12_t,
            R2_v13_t, R2_v22_t, R2_v23_t, R2_v32_t,
            R2_v33_t, R2_m2_t, R2_m3_t, sigma2 / totvar,
            R2_f1_1, 0, 0, R2_v12_1,
            R2_v13_1, 0, 0, 0,
            0, 0, 0, sigma2 / l1var,
            0, R2_f2_2, 0, 0,
            0, R2_v22_2, R2_v23_2, 0,
            0, R2_m2_2, 0, 0,
            0, 0, R2_f3_3, 0,
            0, 0, 0, R2_v32_3,
            R2_v33_3, 0, R2_m3_3, 0
          ),
          nrow = 12,
          ncol = 4
        )
      colnames(contributions_stacked) <-
        c("total", "level-1", "level-2", "level-3")
      rownames(contributions_stacked) <- c("f1",
                                           "f2",
                                           "f3",
                                           "v12",
                                           "v13",
                                           "v22",
                                           "v23",
                                           "v32",
                                           "v33",
                                           "m2",
                                           "m3",
                                           "resid")

      barplot(
        contributions_stacked,
        main = "Decomposition",
        horiz = FALSE,
        ylim = c(0, 1),
        col = c(
          "darkred",
          "steelblue",
          "darkgoldenrod1",
          "darkred",
          "darkred",
          "midnightblue",
          "midnightblue",
          "darkgoldenrod1",
          "darkgoldenrod1",
          "midnightblue",
          "darkgoldenrod1",
          "white"
        ),
        ylab = "proportion of variance",
        density = c(NA, NA, NA, 20, 40, 40, 20, 40, 20, 40, 20, NA),
        angle = c(0, 0, 0, 45, 45, 90, 135, 90, 90, 135, 0, 0)
      )
      legend(
        4.89,
        .7,
        title = "source",
        legend = rownames(contributions_stacked),
        fill = c(
          "darkred",
          "steelblue",
          "darkgoldenrod1",
          "darkred",
          "darkred",
          "midnightblue",
          "midnightblue",
          "darkgoldenrod1",
          "darkgoldenrod1",
          "midnightblue",
          "darkgoldenrod1",
          "white"
        ),
        cex = .7,
        pt.cex = 1,
        xpd = TRUE,
        density = c(NA, NA, NA, 30, 50, 50, 30, 50, 30, 50, 30, NA),
        angle = c(0, 0, 0, 45, 45, 90, 135, 90, 90, 135, 0, 0)
      )


      Output <- list(noquote(R2_table))
      names(Output) <- c("R2s")
    }

    return(Output)
  }
