# r2mlm v 0.1.1

* Fixed typo in `r2mlm_manual()` documentation (#33)
* Refactored to increase modularity, adding files: utils.R, prepare_data.R
* Can now accept data with missing points, handles it with listwise deletion via broom::augment(model). (#23, #29)
* Bug fix: when groups of 1 exist, variance was returning as NA, generating "Error in if (variance_tracker == 0) { : missing value where TRUE/FALSE needed" (#26)
