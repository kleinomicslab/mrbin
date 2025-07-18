% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrbin.R
\name{mrplot}
\alias{mrplot}
\title{A function for plotting NMR spectra.}
\usage{
mrplot(
  hideMenu = FALSE,
  folders = NULL,
  NMRvendor = "Bruker",
  dimensions = NULL,
  intensity1D = NULL,
  zoom = NULL,
  color = NULL,
  background = NULL,
  lwd = 1,
  plotTitle = "",
  showNames = "Spectrum titles",
  graphics = TRUE,
  highlight = NULL,
  binlist = NULL,
  annotate = NULL,
  metaboliteIdentities = NULL,
  annotateColors = c("black", "red", "orange3", "yellow4", "green3", "blue", "purple",
    "violet", "brown4", "chartreuse4", "blue4", "deeppink", "orangered", "olivedrab",
    "cadetblue", "tomato3"),
  annotateAngles = c(35, -35, 20, -20, 45, -45, 60, -60, 75, -75, 15, -15),
  hideExcludedAnnotations = FALSE,
  ...
)
}
\arguments{
\item{hideMenu}{Do not show the menu. Defaults to FALSE}

\item{folders}{Optional vector of folder names of spectra to load. Defaults to NULL}

\item{NMRvendor}{Optional value naming the format of spectra to load. Defaults to Bruker}

\item{dimensions}{Optional vector dimensions of spectra to load. Defaults to NULL}

\item{intensity1D}{Optional value of initial 1D intensity. Defaults to NULL}

\item{zoom}{Optional vector of initial zoom area. Defaults to NULL}

\item{color}{Defines the color of the spectrum plot. If NULL, a rainbow theme is used for 2D NMR}

\item{background}{Background color, defaults to NULL (no background fill, usually results in a white background)}

\item{lwd}{Line width, defaults to 1}

\item{plotTitle}{Plot title, defaults to "" (empty)}

\item{showNames}{Display list of spectrum titles in plot, defaults to "Spectrum titles". Other options are "" and "Folder names"}

\item{graphics}{Controls whether pop-up windows are shown for selections. Defaults to TRUE.}

\item{highlight}{A vector of up to 2 frequencies that will be highlighted in the plot. If 2 values are provided the distance in Hz is shown as well. Defaults to NULL.}

\item{binlist}{Optional: A vector containing bin names as they are generated by mrbin. These bins will be marked by rectangles in the plot. This could be useful for metabolite identification when having a list of significantly changing signals. Default is NULL.}

\item{annotate}{Should peak annotation regions be shown?}

\item{metaboliteIdentities}{Optional: A file path or 4-column matrix where each row belongs to one unique metabolite signal (left, right, top, bottom borders). Row names are metabolite names. For a file, this needs to be the file path for a .csv file containing such a matrix, the first columns containing metabolite names and the first row being a header. Each row belongs to one unique metabolite signal (left, right, top, bottom borders). Row names are metabolite names.}

\item{annotateColors}{Colors for annotation boxes}

\item{annotateAngles}{Angles for shading of annotation boxes}

\item{hideExcludedAnnotations}{Should excluded peak annotation regions be hidden?}

\item{...}{Additional graphical parameters that will be passed to the functions plot, lines, and/or contour}
}
\value{
{None}
}
\description{
This function plots NMR spectra. A menu of commands is displayed to edit the
plot view and add spectra. Multiple spectra will be overlaid, and if both
1D and 2D spectra are selected, they are shown in two plots with matched ranges.
}
\examples{
resetEnv()#clean up previous data from the package environment 
metaboliteIdentities=matrix(c(1.346,1.324,21,23,1,1,
                              4.12,4.1,70.8578,71.653,0,1,
                              3.052,3.043,30.5,33.5,1,1,
                              4.066,4.059,57,59.5,1,0,
                              2.582,2.479,46,49,1,1,
                              2.737,2.634,46,49,1,1),
                   ncol=6,byrow=TRUE)
rownames(metaboliteIdentities)=c("Lactate","Lactate","Creatinine","Creatinine","Citrate","Citrate")
colnames(metaboliteIdentities)=c("left","right","top","bottom","usePeak1D","usePeak2D")
mrplot(folders=c(system.file("extdata/1.mr2",package="mrbin"),
                 system.file("extdata/1.mr1",package="mrbin"),
                 system.file("extdata/2.mr1",package="mrbin"),
                 system.file("extdata/3.mr1",package="mrbin")),
       NMRvendor="mrbin",
       example=TRUE,#only used for the package examples
       dimensions=c("2D","1D","1D","1D"),zoom=c(2.8,2.4,20,60),
       highlight=c(2.564,2.537),
       binlist=c("2.725,2.675","2.575,2.525"),
       annotate=TRUE,metaboliteIdentities=metaboliteIdentities,
       plotTitle="Significant Bins",intensity1D=24,hideMenu=TRUE)
}
