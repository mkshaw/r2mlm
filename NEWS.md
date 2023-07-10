# r2mlm 0.3.4

## Minor Edits

* Replaces rockchalk dependency with misty v. 0.4.12.

# r2mlm 0.3.3

## Minor Edits

* Fixes typo in r2mlm3_manual (#59).
* Fixes typo in r2mlm3_manual (#62).

# r2mlm 0.3.2

## Major Changes

* Output now returns as numeric rather than characters. (#55)

## Minor Edits

* Removes broomExtra dependency. (#52, #57)
* Changes how variable types are checked from if() to is().

# r2mlm 0.3.1

## Major Changes
* Exported `r2mlm_long_manual` to be user-facing.

## Minor Edits
* Updated `r2mlm_manual` and `r2mlm_comp_manual` documentation to reflect changes to `teachsat` dataset implemented in version 0.3.0. (#53)

# r2mlm 0.3.0

## Major Changes
* Adds two manual functions: one for 3-level models (r2mlm3_manual) and one for models with heteroscedasticity, autocorrelation, nonlinearity, and non-centered-within-cluster models (r2mlm_long_manual)
* Bar graph output is now optional. The default behaviour is to output bar graphs, but if you don't want graphical output, the argument is `bargraph = FALSE`. For example, `r2mlm(model, bargraph = FALSE)`. (Issue #46)

## Bug Fixes
* To test whether clusters are mean-centered, the code computes cluster means for all level-1 variables, sees if the means are roughly zero (< .0000001), and if yes then it assigns `clustermeancentered = TRUE`. This update changes the code to test whether the *absolute value* of the means are roughly zero, to address the case in which a cluster has a negative non-zero mean (that would otherwise mistakenly be assigned to `clustermeancentered = TRUE` because the negative number is less than 0.0000001). (Issue #41)
* Fixes an issue where models with non-cwc interaction terms were returning results as though they were centered-within-cluster. r2mlm returns non-cwc results, r2mlm_comp breaks. (Issue #42)
* Fixed an error thrown if certain groups only have one unit: "Error in if (variance_tracker == 0) { : missing value where TRUE/FALSE needed." Fixed this (#44).

## Minor Edits
* Changed simulated data

# r2mlm 0.2.0

## Major Changes
* Can now accept data with missing points, handles it with listwise deletion via broomExtra::augment(model). (#23, #29)
* Related to accepting missing data (#23, #29), this update changes `r2mlm_comp()` to accept optional data argument. You can now call `r2mlm_comp(modelA, modelB)` or `r2mlm_comp(modelA, modelB, data)`. If data is provided, the function will use that data. If data is not provided and models are hierarchically nested, the function will extract data automatically. If data is not provided and models are not hierarchically nested, the function will throw an error asking users to input data.

## Bug Fixes
* Bug fix: when groups of 1 exist, variance was returning as NA, generating "Error in if (variance_tracker == 0) { : missing value where TRUE/FALSE needed" (#26)

## Minor Edits
* Fixed typo in `r2mlm_manual()` documentation (#33)
* Updates documentation of `r2mlm()` and `r2mlm_comp()` to note that models run in `lme4` must be formatted with random effects at the end of the formula. (#30)
* Refactored to increase modularity, adding files: utils.R, prepare_data.R
