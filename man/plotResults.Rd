% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrbin.R
\name{plotResults}
\alias{plotResults}
\title{A function for plotting quality indicators, including PCA plots.}
\usage{
plotResults(mrbinResults, defineGroups = TRUE, process = TRUE, silent = FALSE)
}
\arguments{
\item{mrbinResults}{An mrbin object}

\item{defineGroups}{Should group membership be highlighted in PCA?}

\item{process}{If set to FALSE, the file name will be extended by "Raw" to indicate that data has not been processed yet}

\item{silent}{If set to TRUE, plots will be saved but not shown for the binning step for speed purposes}
}
\value{
{None}
}
\description{
This function plots boxplots (bin-wise and sample-wise) as visual quality indicators.
It also performs PCA, then plots PC1 and PC2 and loading plots.
}
\examples{
resetEnv()#clean up previous data from the package environment 
results<-mrbin(silent=TRUE,setDefault=FALSE,parameters=list(dimension="2D",
    binRegion=c(8,1,15,140),binwidth2D=0.2,binheight=4,solventRegion=c(5.5,4.2),
    PQNScaling="No",noiseRemoval="Yes",trimZeros="Yes",
    fixNegatives="No",logTrafo="No",PCA="No",signal_to_noise2D=10,
    example=TRUE,#only used for the package examples
    NMRfolders=c(system.file("extdata/1.mr2",package="mrbin"),
                 system.file("extdata/2.mr2",package="mrbin")),
    NMRvendor="mrbin"))
plotResults(results)
}
