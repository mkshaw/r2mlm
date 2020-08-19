
<!-- README.md is generated from README.Rmd. Please edit that file -->

# r2mlm

<!-- badges: start -->

[![CRAN\_Status\_Badge](https://img.shields.io/cran/v/r2mlm)](https://cran.r-project.org/package=r2mlm)
<!-- badges: end -->

r2mlm generates R-squared measures for multilevel models, based on the
framework described in [Rights & Sterba
(2019)](https://psycnet.apa.org/record/2018-32178-001). You provide a
multilevel model generated using
[lme4](https://cran.r-project.org/web/packages/lme4/index.html) or
[nlme](https://cran.r-project.org/web/packages/nlme/index.html), and
r2mlm outputs both total- and level-specific variance explained.

## Installation

You can install the released version of r2mlm from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("r2mlm")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mkshaw/r2mlm")
```

## Example

As an example, say you have a dataset of teachers clustered within
schools. Your model has teacher job satisfaction as the outcome,
predicted by teacher salary (within-school, level 1 variable) and
student-teacher ratio (between-school, level 2 variable).

``` r
library(r2mlm)
#> Loading required package: lme4
#> Loading required package: Matrix
#> Loading required package: nlme
#> 
#> Attaching package: 'nlme'
#> The following object is masked from 'package:lme4':
#> 
#>     lmList

# Generate the model, in this case with lme4:
model <- lmer(satisfaction ~ 1 + salary_c + s_t_ratio + (1 | schoolID), data = teachsat, REML = TRUE)

# Generate R-squared measures for that model:
r2mlm(model)
```

<img src="man/figures/README-example-1.png" width="100%" />

    #> $Decompositions
    #>                 total              within            between          
    #> fixed, within   0.0812033364578788 0.141733976636916 NA               
    #> fixed, between  0.0644784051104297 NA                0.150977764921643
    #> slope variation 0                  0                 NA               
    #> mean variation  0.362593787565715  NA                0.849022235078357
    #> sigma2          0.491724470865976  0.858266023363084 NA               
    #> 
    #> $R2s
    #>     total              within            between          
    #> f1  0.0812033364578788 0.141733976636916 NA               
    #> f2  0.0644784051104297 NA                0.150977764921643
    #> v   0                  0                 NA               
    #> m   0.362593787565715  NA                0.849022235078357
    #> f   0.145681741568309  NA                NA               
    #> fv  0.145681741568309  0.141733976636916 NA               
    #> fvm 0.508275529134024  NA                NA

## Package Functions

There are two main functions currently available in r2mlm:

1.  `r2mlm()`, for computing variance explained for a single model.
2.  `r2mlm_comp()`, for comparing variance explained between two
    different models.

In some cases, you might run a model that will not converge in `nlme` or
`lme4`, but will converge in another software (e.g., HLM, Mplus). If you
run your model and have the associated output, you can manually generate
R-squared measures using `r2mlm_manual()` or `r2mlm_comp_manual()`. For
manual entry details, see the help pages:

``` r
?r2mm_manual()
?r2mlm_comp_manual()
```

## Framework Assumptions

This framework of variance explained assumes the following:

1.  Two-level multilevel linear models (3 or more levels not supported
    at this time);
2.  Normal outcome variable;
3.  Homoscedastic residual variances at both level 1 and level 2;
4.  Level-1 predictors are centered within cluster.

## Getting Help
