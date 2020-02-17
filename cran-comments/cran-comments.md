## Resubmission

This is a resubmission. In this version I have changed the following:

* Added a reference web site to DESCRIPTION file
* Replaced T by TRUE and F by FALSE
* Removed getwd() and added user input for specifying directories
* Added \value to all .Rd files
* Information messages written to the console are now generated using message()
* data folder was added, containing data to enable examples to run
* on.exit() is used to ensure par() parameters are reset on exit
* \donttest{} was removed, functioning examples were created for selected functions

Additional changes can be found in the NEWS.md file.

## Test environments
* local Windows 10 install, R 3.5.1
* Windows, using devtools::check_win_devel()
* Linux, using devtools::check_rhub()

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

On NOTE was produced:
 checking CRAN incoming feasibility ... NOTE
  Maintainer: ‘Matthias Klein <klein.663@osu.edu>’

It seems that this is not a note I could possibly fix, but rather an internal note.

## Downstream dependencies
There are currently no downstream dependencies for this package.


## Previous cran-comments

## Test environments
* local Windows 10 install, R 3.5.1
* Windows, using devtools::check_win_devel(), devtools::check_win_release(), and devtools::check_win_oldrelease()
* Linux, using devtools::check_rhub()

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

On NOTE was produced:
 checking CRAN incoming feasibility ... NOTE
  Maintainer: ‘Matthias Klein <klein.663@osu.edu>’

It seems that this is not a note I could possibly fix, but rather an internal note.

## Downstream dependencies
There are currently no downstream dependencies for this package.

