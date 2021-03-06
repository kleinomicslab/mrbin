% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrbin.R
\name{mrbinrun}
\alias{mrbinrun}
\title{A function performing all data read and processing steps.}
\usage{
mrbinrun()
}
\value{
{None}
}
\description{
This function reads parameters from the global variable mrbin.env$mrbinparam and
performs the following operations:
Reading NMR files, creating bins, removing solvent area, removing additional
user-defined areas, summing up bins that contain unstable peaks such as
citric acid, removes noise bins, crops HSQC spectra to the diagonal area,
performs PQN scaling, replaces negative values, log transforms and displays a
PCA plot. Parameters are then saved in a text file. These can be recreated
using recreatemrbin().
}
\examples{
setParam(parameters=list(dimension="2D",binwidth2D=0.05,binheight=3,PQNScaling="No",
         fixNegatives="No",logTrafo="No",signal_to_noise2D=20,
         NMRfolders=c(system.file("extdata/1/12/pdata/10",package="mrbin"),
                      system.file("extdata/2/12/pdata/10",package="mrbin"))))
mrbinrun()
}
