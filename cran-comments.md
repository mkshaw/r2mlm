## Test environments
* local macOS install, R 3.6.3
* ubuntu 16.04 (on travis-ci), R 4.0.2
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs.

There were 3 NOTEs:

* running examples for arch 'i386' ... [30s] NOTE
Examples with CPU (user + system) or elapsed time > 10s
            user system elapsed
r2mlm_comp 12.89   0.08   13.14

* running examples for arch 'x64' ... [30s] NOTE
Examples with CPU (user + system) or elapsed time > 10s
            user system elapsed
r2mlm_comp 13.42      0   13.58

The examples compare results for two models, which is somewhat computationally intensive.

* Found the following URLs which should use \doi (with the DOI name only):
  File 'r2mlm.Rd':
    https://doi.org/10.1037/met0000184
  File 'r2mlm_comp.Rd':
    https://doi.org/10.1037/met0000184
    https://doi.org/10.1080/00273171.2019.1660605
  File 'r2mlm_comp_manual.Rd':
    https://doi.org/10.1037/met0000184
    https://doi.org/10.1080/00273171.2019.1660605
  File 'r2mlm_manual.Rd':
    https://doi.org/10.1037/met0000184
    
Those URLs are hyperlinks created using href, which requires links as input.

## Downstream dependencies
There are currently no downstream dependencies for this package.
