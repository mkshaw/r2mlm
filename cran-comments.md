## Test environments
* local macOS install, R 3.6.3
* ubuntu 16.04 (on travis-ci), R 4.0.2
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs.

There were 2 NOTEs:

* running examples for arch 'i386' ... [30s] NOTE
Examples with CPU (user + system) or elapsed time > 10s
            user system elapsed
r2mlm_comp 12.89   0.08   13.14

* running examples for arch 'x64' ... [30s] NOTE
Examples with CPU (user + system) or elapsed time > 10s
            user system elapsed
r2mlm_comp 13.42      0   13.58

The examples compare results for two models, which is somewhat computationally intensive.

## Downstream dependencies
There are currently no downstream dependencies for this package.
