## Resubmission

## Test environments
* local macOS install, R 3.6.3
* ubuntu 16.04 (on travis-ci), R 4.0.2
* win-builder (devel and release)

## R CMD check results

There was 1 NOTE:
Found the following (possibly) invalid URLs:
  URL: https://doi.org/10.1037/met0000184
    From: README.md
    Status: 400
    Message: Bad Request

Found the following (possibly) invalid DOIs:
  DOI: 10.1037/met0000184
    From: DESCRIPTION
    Status: Bad Request
    Message: 400

The DOI works for me when I paste it into a browser, and is properly enclosed in angle brackets in the DESCRIPTION file.

There were no ERRORs, WARNINGs.

## Downstream dependencies
There are currently no downstream dependencies for this package.
