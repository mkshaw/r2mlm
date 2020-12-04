# r2mlm v 0.1.1

* Fixed typo in `r2mlm_manual()` documentation (#33)
* Refactored to increase modularity, adding files: utils.R, prepare_data.R
* Can now accept data with missing points, handles it with listwise deletion via broomExtra::augment(model). (#23, #29)
* Related to accepting missing data (#23, #29), this update changes `r2mlm_comp()` to accept optional data argument. You can now call `r2mlm_comp(modelA, modelB)` or `r2mlm_comp(modelA, modelB, data)`. If data is provided, the function will use that data. If data is not provided and models are hierarchically nested, the function will extract data automatically. If data is not provided and models are not hierarchically nested, the function will throw an error asking users to input data.
* Bug fix: when groups of 1 exist, variance was returning as NA, generating "Error in if (variance_tracker == 0) { : missing value where TRUE/FALSE needed" (#26)
* Updates documentation of `r2mlm()` and `r2mlm_comp()` to note that lme4 must be formatted with random effects at the end of the formula. (#30)
