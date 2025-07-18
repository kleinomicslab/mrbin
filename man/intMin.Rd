% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrbin.R
\name{intMin}
\alias{intMin}
\title{A function for changing plotNMR plots.}
\usage{
intMin(dimension = "1D", refreshPlot = TRUE, value = NULL)
}
\arguments{
\item{dimension}{Dimension to use. Defaults to "1D"}

\item{refreshPlot}{Refresh plot automatically. Defaults to TRUE}

\item{value}{Set exact value. Defaults to NULL}
}
\value{
{None}
}
\description{
This function decreases the intensity of the current NMR spectrum plot.
}
\examples{
resetEnv()#clean up previous data from the package environment 
mrbin(silent=TRUE,setDefault=TRUE,parameters=list(dimension="1D",binwidth1D=.1,
         PQNScaling="No",noiseRemoval="No",trimZeros="No",
         fixNegatives="No",logTrafo="No",PCA="No",verbose=TRUE,
         NMRvendor="mrbin",
         example=TRUE,#only used for the package examples
         NMRfolders=system.file("extdata/1.mr1",package="mrbin")))
plotNMR()
intMin()
}
