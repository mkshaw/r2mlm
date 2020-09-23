## Resubmission

This is a resubmission. In this version I have:

* Fixed README URLs to be in canonical form
* Added DOIs to references in the Description field of the DESCRIPTION file

## Test environments
* local macOS install, R 3.6.3
* ubuntu 16.04 (on travis-ci), R 4.0.0
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs.

There were 3 NOTEs:

* checking for future file timestamps (710ms)
   unable to verify current time
   
Based on my Google searching, I believe this is an issue with worldclockapi.com.

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Mairead Shaw <mairead.shaw@mail.mcgill.ca>'

New submission

This is my first CRAN submission.

* Possibly mis-spelled words in DESCRIPTION:
  
    Sterba (12:16, 21:46)
    
Sterba is the last name of one of the authors of this package and a paper referenced in the documentation.

## Downstream dependencies
There are currently no downstream dependencies for this package.
