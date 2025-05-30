% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrbin.R
\name{addToPlot}
\alias{addToPlot}
\title{A function for adding NMR spectra to the plot list.}
\usage{
addToPlot(
  folder = NULL,
  dimension = "1D",
  NMRvendor = "Bruker",
  useAsNames = "Folder names",
  add = TRUE,
  omitCurrent = FALSE
)
}
\arguments{
\item{folder}{Defines the exact NMR data folder. If NULL, mrbin parameter set is used}

\item{dimension}{Defines the data dimension, "1D" or "2D". Only used if not NULL}

\item{NMRvendor}{Defines the NMR manufacturer, default is "Bruker"}

\item{useAsNames}{How should sample names be generated}

\item{add}{Add spectra to existing list, or replace existing spectra. Default is TRUE}

\item{omitCurrent}{Omit the "current spectrum" spot and start filling the additional lists immediately. Default is FALSE}
}
\value{
{none}
}
\description{
This function adds a spectrum to the plot list.
}
\examples{
\donttest{ addToPlot() }
}
