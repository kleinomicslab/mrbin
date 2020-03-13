#mrbin - Collection of R functions for analyzing NMR metabolomics data.
#Written by Matthias Klein, The Ohio State University
#
#Package: mrbin
#Title: Magnetic Resonance Binning, Integration and Normalization
#Version: 1.3.0.9000
#Authors@R:
#    person(given = "Matthias",
#           family = "Klein",
#           role = c("aut", "cre"),
#           email = "klein.663@osu.edu",
#           comment = c(ORCID = "0000-0001-7455-5381"))
#Description: Nuclear Magnetic Resonance is widely used in the Life Sciences.
#    This package converts 1D or 2D data into a matrix of values
#    suitable for further data analysis and performs basic processing steps in a
#    reproducible way. Negative values, a common issue in these data, are replaced
#    by positive values. All used parameters are stored in a readable text file
#    and can be restored from that file to enable exact reproduction of the data
#    at a later time.
#License: GPL-3

#' @importFrom grDevices colorRamp heat.colors rainbow rgb devAskNewPage dev.copy dev.off pdf
#' @importFrom graphics axis contour hist legend lines par plot text boxplot points rect polygon box
#' @importFrom stats heatmap median prcomp quantile sd
#' @importFrom utils flush.console select.list write.csv
NULL

#' A function replacing negative values.
#'
#' This function replaces (column-wise) negative values by a small positive
#' number. The number is calculated as an affine transformation to the range of
#' the lowest positive number to 0,01*the lowest positive number (of this
#' column). Ranks stay unchanged. Positive numbers are not altered.
#' If sample-wise noise levels are available, the median noise level of samples
#' with negative values is calculated and replaces the lowest positive number in
#' case it is smaller. If no noise data is available, the 1% percentile of all
#' positive values in the data set is used as an estimate.
#' It is recommended to us this function AFTER noise removal and other data
#' clean-up methods, as it may alter (reduce) the noise level.
#' If no NMR data and noise levels are provided as arguments, the function will
#' use NMR data and noise levels from the global variables mrbin.env$bins and
#' mrbin.env$mrbinTMP.
#' @param NMRdata A matrix containing NMR data. Columns=frequencies,rows=samples
#' @param noiseLevels A vector
#' @return NMRdata An invisible matrix containing NMR data without negative values.
#' @export
#' @examples
#'  Example<-mrbin(silent=TRUE,
#'                    parameters=list(verbose=FALSE,dimension="1D",PQNScaling="No",
#'                    binwidth1D=0.005,signal_to_noise1D=1,PCA="No",binRegion=c(9.5,7.5,10,156),
#'                    saveFiles="No",referenceScaling="No",noiseRemoval="No",
#'                    fixNegatives="No",logTrafo="No",noiseThreshold=.05,
#'                    NMRfolders=c(system.file("extdata/2/10/pdata/10",package="mrbin"),
#'                               system.file("extdata/3/10/pdata/10",package="mrbin"))
#'                    ))
#'  sum(Example$bins<=0)
#'  exampleNMRpositive<-atnv(NMRdata=Example$bins, noiseLevels=Example$parameters$noise_level)
#'  sum(exampleNMRpositive<=0)

atnv<-function(NMRdata=NULL,noiseLevels=NULL){
 if(!is.null(mrbin.env$bins)|!is.null(NMRdata)){
     if(is.null(NMRdata)){
       if("bins"%in%ls(envir=mrbin.env)&"mrbinparam"%in%ls(envir=mrbin.env)){
         NMRdata <- mrbin.env$bins
         noiseLevels <- mrbin.env$mrbinparam$noise_level
       }
     }
     if(is.null(noiseLevels)){
             noiseTMP<-sort(NMRdata[NMRdata>0])[ceiling(.01*length(NMRdata[NMRdata>0]))]
     }
     for(i in 1:ncol(NMRdata)){
        negatives<-NMRdata[,i]<=0

        if(!is.null(noiseLevels)){#If noise levels are available, restrict range to below noise
           if(length(negatives)>0){
             noiseTMP<-stats::median(mrbin.env$mrbinparam$noise_level[negatives])
           } else {
             noiseTMP<-stats::median(mrbin.env$mrbinparam$noise_level)
           }
        }
        if(sum(negatives)>0){
            minTMP<-min(NMRdata[negatives,i])#select lowest bin
            maxTMP<-min(noiseTMP,min(NMRdata[!negatives,i]))#select lowest bin above 0
            NMRdata[negatives,i]<-(NMRdata[negatives,i]+(maxTMP-minTMP))/
                                          (maxTMP-minTMP)*(maxTMP*.99)+maxTMP*.01
        }
     }
     if("bins"%in%ls(envir=mrbin.env)&"mrbinparam"%in%ls(envir=mrbin.env)){
          if(nrow(mrbin.env$bins)==1){
            mrbin.env$bins<-matrix(NMRdata,nrow=1)
            rownames(mrbin.env$bins)<-rownames(NMRdata)
            colnames(mrbin.env$bins)<-colnames(NMRdata)
          } else {
            mrbin.env$bins<-NMRdata
          }
     }
     invisible(NMRdata)
 }
}

#' A function executed when loading this package
#'
#' This function resets the parameter variables.
#' @param  libname Library name
#' @param  pkgname Package name
#' @return {None}
#' @export
#' @examples
#' \donttest{ .onLoad() }

.onLoad <- function(libname, pkgname){
    assign("mrbin.env",new.env(emptyenv()),parent.env(environment()))
    resetEnv()
}

#' A parameter resetting function
#'
#' This function resets the parameter variables.
#' @return {None}
#' @export
#' @examples
#' resetEnv()

resetEnv<-function(){
    if(!exists("mrbin.env", mode="environment")) .onLoad()
    assign("bins",NULL,mrbin.env)
    assign("paramChangeFlag",FALSE,mrbin.env)
    assign("mrbinTMP",list(
               mrbinversion="1.3.0.9000",
               binsRaw=NULL,
               medians=NULL,
               noise_level_Raw=NULL,
               noise_level_Raw_TMP=NULL,
               meanNumberOfPointsPerBin=NULL,
               meanNumberOfPointsPerBin_TMP=NULL,
               PCA=NULL,
               binTMP=NULL,
               binNames=NULL,
               nbins1=NULL,
               nbins2=NULL,
               nbins=NULL,
               binRegions=matrix(ncol=4,nrow=0,dimnames=list(NULL,c("left","right","top","bottom"))),
               currentFolder=NULL,
               currentSpectrum=NULL,
               currentSpectrumName=NULL,
               currentSpectrumFolderName=NULL,
               currentSpectrumEXPNO=NULL,
               currentSpectrumFolderName_EXPNO=NULL,
               currentSpectrumTitle=NULL
    ),mrbin.env)
    assign("requiredParam1D",c(
               "binwidth1D","reference1D","signal_to_noise1D","noiseRange1d",
               "dimension","binMethod","binRegion","specialBinList","referenceScaling",
               "removeSolvent","removeAreas","sumBins","noiseRemoval","PQNScaling",
               "fixNegatives","logTrafo","defineGroups","PCA","solventRegion",
               "removeAreaList","sumBinList","noiseThreshold","PQNshowHist",
               "PQNminimumFeatures","PCAtitlelength","createBins","useAsNames","saveFiles",
               "verbose","Factors","NMRfolders"
               ),mrbin.env)
    assign("requiredParam2D",c(
               "binwidth2D","binheight","cropHSQC","reference2D","signal_to_noise2D",
               "noiseRange2d","croptopRight","croptopLeft","cropbottomRight","cropbottomLeft",
               "PQNsugarArea",
               "dimension","binMethod","binRegion","specialBinList","referenceScaling",
               "removeSolvent","removeAreas","sumBins","noiseRemoval","PQNScaling",
               "fixNegatives","logTrafo","defineGroups","PCA","solventRegion",
               "removeAreaList","sumBinList","noiseThreshold","PQNshowHist",
               "PQNminimumFeatures","PCAtitlelength","createBins","useAsNames","saveFiles",
               "verbose","Factors","NMRfolders"
               ),mrbin.env)
    assign("mrbinparam", list(
               dimension="1D",
               binMethod="Rectangular bins",#"Custom bin list"
               binwidth1D=.01,
               binwidth2D=.02,
               binheight=1,
               binRegion=c(9.5,0.5,10,156),
               specialBinList=matrix(ncol=4,nrow=0,dimnames=list(NULL,c("left","right","top","bottom"))),
               referenceScaling="Yes",
               removeSolvent="Yes",
               removeAreas="No",
               sumBins="No",
               noiseRemoval="Yes",
               cropHSQC="Yes",
               PQNScaling="Yes",
               fixNegatives="Yes",
               logTrafo="Yes",
               defineGroups="Yes",
               PCA="Yes",
               reference1D=c(.03,-0.03),
               reference2D=c(.04,-0.04,-2,2),
               solventRegion=c(4.95,4.65),
               removeAreaList=matrix(ncol=4,nrow=0,dimnames=list(NULL,c("left","right","top","bottom"))),
               sumBinList=matrix(ncol=4,nrow=0,dimnames=list(NULL,c("left","right","top","bottom"))),
               noiseThreshold=0.75,
               signal_to_noise1D=25,
               signal_to_noise2D=6,
               noiseRange2d=c(3.3,2.3,90,110),
               noiseRange1d=c(10,9.5),
               croptopRight=c(0,-1.50),#only 2D, defines edge points of the cropped area
               croptopLeft=c(0,3.5),
               cropbottomRight=c(160,6),
               cropbottomLeft=c(160,10),
               PQNsugarArea=c(5.4,3.35,72,100),#exclude most glucose to reduce impact
               PQNshowHist=FALSE,#show histograms of quotients
               PQNminimumFeatures=40,#Top number of features to proceed
               PCAtitlelength=8,
               createBins="Yes",
               useAsNames="Folder names",#"Folder names and EXPNO", "Spectrum titles"
               NMRvendor="Bruker",#NMR vendor. Currently, only Bruker data is supported.
               mrbinversionTracker=mrbin.env$mrbinTMP$mrbinversion,
               saveFiles="No",
               verbose=TRUE,
               outputFileName=NULL,
               NMRfolders=NULL,
               Factors=NULL,
               medians=NULL,
               noise_level=NULL,
               numberOfFeaturesRaw=NULL,
               numberOfFeaturesAfterRemovingSolvent=NULL,
               numberOfFeaturesAfterRemovingAreas=NULL,
               numberOfFeaturesAfterSummingBins=NULL,
               numberOfFeaturesAfterNoiseRemoval=NULL,
               numberOfFeaturesAfterCropping=NULL
               ),mrbin.env)
    assign("mrbinparam_copy",mrbin.env$mrbinparam,mrbin.env)
    assign("mrbinplot",list(
               lowestContour=.01,
               plotRegion=NULL,
               intensityScale=1,
               nContours=60,
               heatmap=FALSE),mrbin.env)
}

#' A function setting the parameters and performing binning and data processing
#'
#' This function guides the user through the set-up of parameters, starts binning
#' and performs the chosen data processing steps
#' If a list of parameters is provided and silent is set to TRUE, no user input
#' is requested and binning and data processing are performed silently.
#' @param parameters Optional: A list of parameters to be used. If omitted, the user will be asked through a series of question to set the parameters.
#' @param silent If TRUE, the user will be asked no questions and binning and data analysis will run according to the current parameters. Defaults to FALSE.
#' @param setDefault If TRUE, all current parameters will be replaced by the default parameters (before loading any provided parameters sets). Defaults to FALSE.
#' @return An invisible list containing bins (data after processing), parameters, and factors
#' @export
#' @examples
#' # Let the user set parameters interactively
#' \donttest{ results <- mrbin() }
#' # Set parameters in command line.
#' mrbinExample<-mrbin(silent=TRUE,setDefault=TRUE,parameters=list(dimension="1D",
#'                 binwidth1D=0.05,signal_to_noise1D=50,
#'                 NMRfolders=c(system.file("extdata/1/10/pdata/10",package="mrbin"),
#'                             system.file("extdata/2/10/pdata/10",package="mrbin"),
#'                             system.file("extdata/3/10/pdata/10",package="mrbin"),
#'                             system.file("extdata/4/10/pdata/10",package="mrbin")),
#'                 Factors=factor(c("Group A","Group A","Group A","Group B","Group B"))))

mrbin<-function(silent=FALSE,setDefault=FALSE,parameters=NULL){
  if(!exists("mrbin.env", mode="environment")) .onLoad()
  if(setDefault) resetEnv()
  if(!is.null(parameters)){
      setParam(parameters)
  }
  if(mrbin.env$mrbinparam$verbose){
       message(paste("\nmrbin version ",mrbin.env$mrbinTMP$mrbinversion,"\n",sep=""))
  }
  stopTMP<-FALSE
  selectionRepeat<-""
  if(silent) startmrbin<-"Start binning now"
  #Create bin list?
  if(!silent){
   selection<-utils::select.list(c("Yes","No"),preselect="Yes",
              title="mrbin: Create bins?",graphics=TRUE)
     if(!selection=="Yes")         stopTMP<-TRUE
     if(selection=="Yes"&!stopTMP){
       selectStep<-1
       lastStepDone<-FALSE
       while(!lastStepDone&!stopTMP){
         if(selectStep==1){
           #Select new parameters?
           selectionNewTMP<-c("Edit parameters","Reload from file")
           if(!is.null(mrbin.env$mrbinparam$NMRfolders)) selectionNewTMP<-c(selectionNewTMP,"Use current parameters")
           selectionRepeat<-utils::select.list(selectionNewTMP,preselect="Edit parameters",
                                              title="Edit parameters or use existing?",graphics=TRUE)
           if(length(selectionRepeat)==0|selectionRepeat=="") stopTMP<-TRUE
           if(selectionRepeat=="Reload from file"&!stopTMP){
             recreatemrbin()
             selectionRepeat2<-utils::select.list(c("Edit parameters","Use parameters from file without changes",
                                              "Go back"),
                                              preselect="Edit parameters",
                                              title="Edit parameters or use as is?",graphics=TRUE)
             if(length(selectionRepeat2)==0|selectionRepeat=="") stopTMP<-TRUE
             if(selectionRepeat2=="Go back"&!stopTMP) selectStep<-selectStep-2
             if(selectionRepeat2=="Edit parameters"&!stopTMP) selectionRepeat<-"Edit parameters"
             if(selectionRepeat2=="Use parameters from file without changes"&!stopTMP) selectionRepeat<-"Use current parameters"
           }
           if(!stopTMP) selectStep<-selectStep+1
         }
         if(selectionRepeat=="Use current parameters"&!stopTMP){
           selectStep<-18
         }
         if(selectionRepeat=="Edit parameters"&!stopTMP){
           if(selectStep==2){
             #1D or 2D data?
             dimension<-utils::select.list(c("1D","2D","Go back"),
                                     preselect=mrbin.env$mrbinparam$dimension,
                                     title="1D or 2D spectra?",graphics=TRUE)
             if(length(dimension)==0|dimension=="") stopTMP<-TRUE
             if(!stopTMP&!dimension=="Go back"){
               if(!identical(mrbin.env$mrbinparam$dimension,dimension)) mrbin.env$paramChangeFlag<-TRUE
               if(dimension%in%c("1D","2D")) mrbin.env$mrbinparam$dimension<-dimension
             }
           if(!stopTMP&dimension=="Go back") selectStep<-selectStep-2
           if(!stopTMP) selectStep<-selectStep+1
           }
           if(selectStep==3){
             if(!stopTMP){
               #Select folders
               if(!is.null(mrbin.env$mrbinparam$NMRfolders)){
                    selectionFolders<-utils::select.list(c("Yes","No",
                                      "Remove spectra from previous list","Go back"),
                                      preselect="Yes",
                                      title="Use previous spectra list?",graphics=TRUE)
                    if(length(selectionFolders)==0|selectionFolders=="") stopTMP<-TRUE
                    if(!stopTMP){
                      if(selectionFolders=="No"){
                        selectionFolders<-selectFolders()
                        if(selectionFolders=="stop")  stopTMP<-TRUE
                      }
                    }
                    if(!stopTMP){
                      if(selectionFolders=="Remove spectra from previous list")  removeSpectrum()
                    }
               } else {
                    selectionFolders<-selectFolders()
                    if(selectionFolders=="stop")  stopTMP<-TRUE
               }
             }
             if(!stopTMP&selectionFolders=="Go back") selectStep<-selectStep-2
             if(!stopTMP) selectStep<-selectStep+1
           }
           if(selectStep==4){
              if(!stopTMP){
                  mrbin.env$mrbinTMP$currentFolder<-mrbin.env$mrbinparam$NMRfolders[1]
                  readNMR()
              }
              if(!stopTMP){
                if(dimension=="1D") dimlength<-2
                if(dimension=="2D") dimlength<-4
                #Use rectangluar bins or use special bin list, e.g. for lipids
                binMethodpreSelect<-mrbin.env$mrbinparam$binMethod
                if(binMethodpreSelect=="Custom bin list") binMethodpreSelect<-"User defined bin list"
                binMethod<-utils::select.list(c("Rectangular bins","User defined bin list","Go back"),
                           preselect=binMethodpreSelect,
                           ,title ="Binning method: ",graphics=TRUE)
                if(length(binMethod)==0|binMethod=="") stopTMP<-TRUE
                if(!stopTMP){
                  if(!binMethod=="Go back"){
                    if(binMethod=="User defined bin list") binMethod<-"Custom bin list"
                    if(!identical(mrbin.env$mrbinparam$binMethod,binMethod)) mrbin.env$paramChangeFlag<-TRUE
                    mrbin.env$mrbinparam$binMethod<-binMethod
                    #Bin region
                    adjRegion<-""
                    if(mrbin.env$mrbinparam$binMethod=="Rectangular bins"){
                      accept<-FALSE
                      while(!accept&!stopTMP){
                        binRegionText<-paste(paste(c("left=","ppm, right=","ppm, top=","ppm, bottom=")[1:dimlength],
                                          mrbin.env$mrbinparam$binRegion[1:dimlength],collapse="",sep=""),"ppm",sep="")
                        plotNMR(region="all",rectangleRegions=matrix(mrbin.env$mrbinparam$binRegion,ncol=4),color="black",
                                manualScale=FALSE,
                                plotTitle=paste("Bin region\n",binRegionText,
                                sep=""))
                        adjRegion<-utils::select.list(c(paste("Keep: ",binRegionText,collapse=""),
                                   "Change..."),preselect=paste("Use ",binRegionText,collapse=""),
                                   title ="Bin region [ppm]: ",graphics=TRUE)
                        if(length(adjRegion)==0|adjRegion=="") stopTMP<-TRUE
                        if(!stopTMP){
                          if(adjRegion=="Change..."&!stopTMP){
                            regionTMP<-readline(prompt=paste("New left border, press enter to keep ",
                                      mrbin.env$mrbinparam$binRegion[1],": ",sep=""))
                            if(!regionTMP==""){
                                mrbin.env$paramChangeFlag<-TRUE
                                mrbin.env$mrbinparam$binRegion[1]<-as.numeric(regionTMP)
                            }
                            regionTMP<-readline(prompt=paste("New right border, press enter to keep ",
                                      mrbin.env$mrbinparam$binRegion[2],": ",sep=""))
                            if(!regionTMP==""){
                                mrbin.env$paramChangeFlag<-TRUE
                                mrbin.env$mrbinparam$binRegion[2]<-as.numeric(regionTMP)
                            }
                            if(mrbin.env$mrbinparam$dimension=="2D"){
                                regionTMP<-readline(prompt=paste("New top border, press enter to keep ",
                                          mrbin.env$mrbinparam$binRegion[3],": ",sep=""))
                                if(!regionTMP=="") {
                                    mrbin.env$paramChangeFlag<-TRUE
                                    mrbin.env$mrbinparam$binRegion[3]<-as.numeric(regionTMP)
                                }
                                regionTMP<-readline(prompt=paste("New bottom border, press enter to keep ",
                                          mrbin.env$mrbinparam$binRegion[4],": ",sep=""))
                                if(!regionTMP=="") {
                                    mrbin.env$paramChangeFlag<-TRUE
                                    mrbin.env$mrbinparam$binRegion[4]<-as.numeric(regionTMP)
                                }
                              }
                          } else {
                             accept<-TRUE
                          }
                        }
                      }
                    }
                  } else {
                     selectStep<-selectStep-2
                  }
                if(!stopTMP) selectStep<-selectStep+1
                }
              }
            }
            if(selectStep==5){ #Define bins
              #Bin width and height
              adjbinRegion<-""
              addbinRegion<-""
              if(mrbin.env$mrbinparam$dimension=="1D"&!stopTMP&mrbin.env$mrbinparam$binMethod=="Rectangular bins"){
                  accept<-FALSE
                  widthAdjust<-""
                  while(!accept&!stopTMP&!widthAdjust=="Go back"){
                    plotNMR(region=c(1.5,1.15,16.5,27),
                            rectangleRegions=matrix(c(1.35+as.numeric(mrbin.env$mrbinparam$binwidth1D),
                                                    1.35,21,21+1),ncol=4),
                            color="black",
                            manualScale=FALSE,
                            plotTitle=paste("Bin size\nwidth=",mrbin.env$mrbinparam$binwidth1D,"ppm",
                            sep=""))
                    binWidthTitle<-paste("Keep: width=",mrbin.env$mrbinparam$binwidth1D,"ppm",sep="")
                    widthAdjust<-utils::select.list(c(binWidthTitle,"Change...","Go back"),
                                 preselect=binWidthTitle,
                                 title ="Bin width [ppm]: ",graphics=TRUE)
                    if(length(widthAdjust)==0|widthAdjust=="") stopTMP<-TRUE
                    if(widthAdjust=="Change..."){
                         widthTMP<-readline(prompt=paste("New 1D bin width, press enter to keep ",
                                   mrbin.env$mrbinparam$binwidth1D,": ",sep=""))
                         if(!widthTMP=="") {
                             mrbin.env$paramChangeFlag<-TRUE
                             mrbin.env$mrbinparam$binwidth1D<-as.numeric(widthTMP)
                         }
                    }
                    if(widthAdjust==binWidthTitle) accept<-TRUE
                  }
                  if(widthAdjust=="Go back"){
                     selectStep<-selectStep-2
                  }
                  if(!stopTMP) selectStep<-selectStep+1
              }
              if(mrbin.env$mrbinparam$dimension=="2D"&!stopTMP&mrbin.env$mrbinparam$binMethod=="Rectangular bins"){
                  accept<-FALSE
                  widthAdjust<-""
                  #heightAdjust<-""
                  while(!accept&!stopTMP&!widthAdjust=="Go back"){
                    plotNMR(region=c(2,1,16.5,27),
                            rectangleRegions=matrix(c(1.5+as.numeric(mrbin.env$mrbinparam$binwidth2D),
                                                    1.5,21,21+mrbin.env$mrbinparam$binheight),ncol=4),
                            color="black",manualScale=FALSE,
                            plotTitle=paste("Bin size\nwidth=",mrbin.env$mrbinparam$binwidth2D,
                                      "ppm, height=",mrbin.env$mrbinparam$binheight,"ppm",sep=""))
                    currentBinSize<-paste("Keep: width=",mrbin.env$mrbinparam$binwidth2D,
                                      "ppm, height=",mrbin.env$mrbinparam$binheight,"ppm",sep="")
                    widthAdjust<-utils::select.list(c(currentBinSize,"Change...","Go back"),
                                 preselect=currentBinSize,
                                 title ="Bin size [ppm]: ",graphics=TRUE)
                    if(length(widthAdjust)==0|widthAdjust=="") stopTMP<-TRUE
                    if(widthAdjust==currentBinSize) accept<-TRUE
                    if(widthAdjust=="Change..."){
                         widthTMP<-readline(prompt=paste("New 2D bin width, press enter to keep ",
                                   mrbin.env$mrbinparam$binwidth2D,": ",sep=""))
                         if(!widthTMP=="") {
                             mrbin.env$paramChangeFlag<-TRUE
                             mrbin.env$mrbinparam$binwidth2D<-as.numeric(widthTMP)
                         }
                         heightTMP<-readline(prompt=paste("New 2D bin height, press enter to keep ",
                                    mrbin.env$mrbinparam$binheight,": ",sep=""))
                         if(!heightTMP=="") {
                           mrbin.env$paramChangeFlag<-TRUE
                           mrbin.env$mrbinparam$binheight<-as.numeric(heightTMP)
                         }
                    }
                  }
                  if(widthAdjust=="Go back"){
                     selectStep<-selectStep-2
                  }
                  if(!stopTMP) selectStep<-selectStep+1
              }
              #Set custom bin list
              if(!stopTMP&mrbin.env$mrbinparam$binMethod=="Custom bin list"){
                adjbinRegion<-""
                if(nrow(mrbin.env$mrbinparam$specialBinList)>1){
                  for(ibinRegions in 1:nrow(mrbin.env$mrbinparam$specialBinList)){
                     if(!stopTMP&!adjbinRegion=="Go back"){
                      mean1<-mean(mrbin.env$mrbinparam$specialBinList[ibinRegions,1:2])
                      range1<-mrbin.env$mrbinparam$specialBinList[ibinRegions,1]-mrbin.env$mrbinparam$specialBinList[ibinRegions,2]
                      mean2<-mean(mrbin.env$mrbinparam$specialBinList[ibinRegions,3:4])
                      range2<-mrbin.env$mrbinparam$specialBinList[ibinRegions,4]-mrbin.env$mrbinparam$specialBinList[ibinRegions,3]
                      if(mrbin.env$mrbinparam$dimension=="1D"){
                        plotNMR(region=c(mean1+4*range1,mean1-4*range1,mean2-4*range2,mean2+4*range2),
                                rectangleRegions=matrix(c(mrbin.env$mrbinparam$specialBinList[ibinRegions,1],
                                                        mrbin.env$mrbinparam$specialBinList[ibinRegions,2],0,2),ncol=4),
                                color="black",
                                manualScale=FALSE,
                                plotTitle=paste("Custom bins\nleft=",mrbin.env$mrbinparam$specialBinList[ibinRegions,1],
                                          "ppm, right=",mrbin.env$mrbinparam$specialBinList[ibinRegions,2],"ppm",sep=""))
                      }
                      if(mrbin.env$mrbinparam$dimension=="2D"){
                        plotNMR(region=c(mean1+4*range1,mean1-4*range1,mean2-4*range2,mean2+4*range2),
                                rectangleRegions=matrix(mrbin.env$mrbinparam$specialBinList[ibinRegions,1:4],ncol=4),
                                color="black",
                                manualScale=FALSE,
                                plotTitle=paste("Custom bins\nleft=",mrbin.env$mrbinparam$specialBinList[ibinRegions,1],
                                          "ppm, right=",mrbin.env$mrbinparam$specialBinList[ibinRegions,2],"ppm",sep=""))
                      }

                      adjbinRegion<-utils::select.list(c(paste(mrbin.env$mrbinparam$specialBinList[ibinRegions,1:dimlength],
                                 collapse=","),"Change...","Go back"),
                                 preselect=paste(mrbin.env$mrbinparam$specialBinList[ibinRegions,1:dimlength],collapse=","),
                                 title =paste("Change region?",rownames(mrbin.env$mrbinparam$specialBinList)[ibinRegions]),graphics=TRUE)
                      if(length(adjbinRegion)==0|adjbinRegion=="") stopTMP<-TRUE
                     }
                    if(adjbinRegion=="Change..."&!stopTMP){
                      #nameTMP<-readline(prompt=paste("New name, press enter to keep ",
                      #          rownames(mrbin.env$mrbinparam$specialBinList)[ibinRegions],": ",sep=""))
                      #if(!nameTMP=="") {
                      #       mrbin.env$paramChangeFlag<-TRUE
                      #       rownames(mrbin.env$mrbinparam$specialBinList)[ibinRegions]<-nameTMP
                      #}
                      regionTMP<-readline(prompt=paste("New left border, press enter to keep ",
                                mrbin.env$mrbinparam$specialBinList[ibinRegions,1],": ",sep=""))
                      if(!regionTMP=="") {
                             mrbin.env$paramChangeFlag<-TRUE
                             mrbin.env$mrbinparam$specialBinList[ibinRegions,1]<-as.numeric(regionTMP)
                      }
                      regionTMP<-readline(prompt=paste("New right border, press enter to keep ",
                                mrbin.env$mrbinparam$specialBinList[ibinRegions,2],": ",sep=""))
                      if(!regionTMP=="") {
                             mrbin.env$paramChangeFlag<-TRUE
                             mrbin.env$mrbinparam$specialBinList[ibinRegions,2]<-as.numeric(regionTMP)
                      }
                      if(mrbin.env$mrbinparam$dimension=="2D"&!stopTMP){
                        regionTMP<-readline(prompt=paste("New top border, press enter to keep ",
                                  mrbin.env$mrbinparam$specialBinList[ibinRegions,3],": ",sep=""))
                        if(!regionTMP=="") {
                               mrbin.env$paramChangeFlag<-TRUE
                               mrbin.env$mrbinparam$specialBinList[ibinRegions,3]<-as.numeric(regionTMP)
                        }
                        regionTMP<-readline(prompt=paste("New bottom border, press enter to keep ",
                                  mrbin.env$mrbinparam$specialBinList[ibinRegions,4],": ",sep=""))
                        if(!regionTMP=="") {
                               mrbin.env$paramChangeFlag<-TRUE
                               mrbin.env$mrbinparam$specialBinList[ibinRegions,4]<-as.numeric(regionTMP)
                        }

                      }
                    }
                  }
                }
                addBinFlag<-TRUE
                addbinRegion<-""
                if(!stopTMP&!adjbinRegion=="Go back"){
                   while(addBinFlag&!addbinRegion=="Go back"){
                      addbinRegion<-utils::select.list(c("Yes","No","Go back"),preselect="No",
                                   title ="Add additional bin?",graphics=TRUE)
                      if(length(addbinRegion)==0|adjbinRegion=="") stopTMP<-TRUE
                      if(!stopTMP){
                        if(addbinRegion=="No") addBinFlag<-FALSE
                        if(addbinRegion=="Yes"){
                          mrbin.env$paramChangeFlag<-TRUE
                          mrbin.env$mrbinparam$specialBinList<-rbind(mrbin.env$mrbinparam$specialBinList,c(0,0,0,160))
                          ibinRegions<-nrow(mrbin.env$mrbinparam$specialBinList)
                          #nameTMP<-readline(prompt="New bin name: ")
                          #if(!nameTMP=="") {
                          #rownames(mrbin.env$mrbinparam$specialBinList)[ibinRegions]<-nameTMP
                          #}
                          regionTMP<-readline(prompt=paste("Left border, press enter to keep ",
                                    mrbin.env$mrbinparam$specialBinList[ibinRegions,1],": ",sep=""))
                          if(!regionTMP=="") {
                                 mrbin.env$paramChangeFlag<-TRUE
                                 mrbin.env$mrbinparam$specialBinList[ibinRegions,1]<-as.numeric(regionTMP)
                          }
                          regionTMP<-readline(prompt=paste("Right border, press enter to keep ",
                                    mrbin.env$mrbinparam$specialBinList[ibinRegions,2],": ",sep=""))
                          if(!regionTMP=="") {
                                 mrbin.env$paramChangeFlag<-TRUE
                                 mrbin.env$mrbinparam$specialBinList[ibinRegions,2]<-as.numeric(regionTMP)
                          }
                          if(mrbin.env$mrbinparam$dimension=="2D"&!stopTMP){
                            regionTMP<-readline(prompt=paste("Top border, press enter to keep ",
                                      mrbin.env$mrbinparam$specialBinList[ibinRegions,3],": ",sep=""))
                            if(!regionTMP=="") {
                                   mrbin.env$paramChangeFlag<-TRUE
                                   mrbin.env$mrbinparam$specialBinList[ibinRegions,3]<-as.numeric(regionTMP)
                            }
                            regionTMP<-readline(prompt=paste("Bottom border, press enter to keep ",
                                      mrbin.env$mrbinparam$specialBinList[ibinRegions,4],": ",sep=""))
                            if(!regionTMP=="") {
                                   mrbin.env$paramChangeFlag<-TRUE
                                   mrbin.env$mrbinparam$specialBinList[ibinRegions,4]<-as.numeric(regionTMP)
                            }

                          }
                       }
                     }
                    }
                  }
                  if(!stopTMP&!adjbinRegion=="Go back"&!addbinRegion=="Go back"){
                    mrbin.env$mrbinTMP$binRegions<-mrbin.env$mrbinparam$specialBinList
                  }
                  if(adjbinRegion=="Go back"|addbinRegion=="Go back"){
                     selectStep<-selectStep-2
                  }
                  if(!stopTMP) selectStep<-selectStep+1
              }
            }
            if(selectStep==6){
              #Scale to reference
              if(!stopTMP){
                adjRegion<-""
                referenceScaling<-utils::select.list(c("Yes","No","Go back"),
                                         preselect=mrbin.env$mrbinparam$referenceScaling,
                                              title = "Scale to reference signal?",graphics=TRUE)
                if(length(referenceScaling)==0|referenceScaling=="") stopTMP<-TRUE
                if(!stopTMP&!referenceScaling=="Go back"){
                  if(!identical(mrbin.env$mrbinparam$referenceScaling,referenceScaling)) mrbin.env$paramChangeFlag<-TRUE
                  mrbin.env$mrbinparam$referenceScaling<-referenceScaling
                  if(mrbin.env$mrbinparam$referenceScaling=="Yes"){
                    if(mrbin.env$mrbinparam$dimension=="1D"){
                    accept<-FALSE
                    adjRegion<-""
                    while(!accept&!stopTMP&!adjRegion=="Go back"){
                      mean1<-mean(mrbin.env$mrbinparam$reference1D)
                      range1<-max(mrbin.env$mrbinparam$reference1D)-min(mrbin.env$mrbinparam$reference1D)
                      plotNMR(region=c(mean1+4*range1,mean1-4*range1,-10,10),
                              rectangleRegions=matrix(c(mrbin.env$mrbinparam$reference1D[1],
                                                      mrbin.env$mrbinparam$reference1D[2],0,2),ncol=4),
                              color="black",
                              manualScale=FALSE,
                              plotTitle=paste("Reference region\nleft=",mrbin.env$mrbinparam$reference1D[1],
                                        "ppm, right=",mrbin.env$mrbinparam$reference1D[2],"ppm",sep=""))
                      RefRegionTitle<-paste("Keep: left=",mrbin.env$mrbinparam$reference1D[1],
                                            "ppm, right=",mrbin.env$mrbinparam$reference1D[2],"ppm",sep="")
                      adjRegion<-utils::select.list(c(RefRegionTitle,
                                 "Change...","Go back"),
                                 preselect=RefRegionTitle,title ="Reference region [ppm]: ",graphics=TRUE)
                      if(length(adjRegion)==0|adjRegion=="") stopTMP<-TRUE
                      if(adjRegion==RefRegionTitle) accept<-TRUE
                      if(adjRegion=="Change..."){
                        regionTMP<-readline(prompt=paste("New left border, press enter to keep ",
                                  mrbin.env$mrbinparam$reference1D[1],": ",sep=""))
                        if(!regionTMP=="") {
                               mrbin.env$paramChangeFlag<-TRUE
                               mrbin.env$mrbinparam$reference1D[1]<-as.numeric(regionTMP)
                        }
                        regionTMP<-readline(prompt=paste("New right border, press enter to keep ",
                                  mrbin.env$mrbinparam$reference1D[2],": ",sep=""))
                        if(!regionTMP=="") {
                               mrbin.env$paramChangeFlag<-TRUE
                               mrbin.env$mrbinparam$reference1D[2]<-as.numeric(regionTMP)
                        }
                      }
                      }
                    }
                    if(mrbin.env$mrbinparam$dimension=="2D"){
                      accept<-FALSE
                      adjRegion<-""
                      while(!accept&!stopTMP&!adjRegion=="Go back"){
                        mean1<-mean(mrbin.env$mrbinparam$reference2D[1:2])
                        range1<-max(mrbin.env$mrbinparam$reference2D[1:2])-min(mrbin.env$mrbinparam$reference2D[1:2])
                        mean2<-mean(mrbin.env$mrbinparam$reference2D[3:4])
                        range2<-max(mrbin.env$mrbinparam$reference2D[3:4])-min(mrbin.env$mrbinparam$reference2D[3:4])
                        plotNMR(region=c(mean1+4*range1,mean1-4*range1,
                                         mean2-4*range2,mean2+4*range2),
                                rectangleRegions=matrix(mrbin.env$mrbinparam$reference2D,ncol=4),
                                color="black",
                                manualScale=FALSE,
                                plotTitle=paste("Reference region\nleft=",mrbin.env$mrbinparam$reference2D[1],
                                          "ppm, right=",mrbin.env$mrbinparam$reference2D[2],
                                          "ppm, top=",mrbin.env$mrbinparam$reference2D[3],
                                          "ppm, bottom=",mrbin.env$mrbinparam$reference2D[4],"ppm",sep=""))
                        RefRegionTitle<-paste("Keep: left=",mrbin.env$mrbinparam$reference2D[1],
                                              "ppm, right=",mrbin.env$mrbinparam$reference2D[2],
                                              "ppm, top=",mrbin.env$mrbinparam$reference2D[3],
                                              "ppm, bottom=",mrbin.env$mrbinparam$reference2D[4],"ppm",sep="")
                        adjRegion<-utils::select.list(c(RefRegionTitle,
                                   "Change...","Go back"),
                                   preselect=RefRegionTitle,title ="Reference region [ppm]: ",graphics=TRUE)
                        if(length(adjRegion)==0|adjRegion=="") stopTMP<-TRUE
                        if(adjRegion==RefRegionTitle) accept<-TRUE
                        if(adjRegion=="Change..."){
                          regionTMP<-readline(prompt=paste("New left border, press enter to keep ",
                                    mrbin.env$mrbinparam$reference2D[1],": ",sep=""))
                          if(!regionTMP=="") {
                                 mrbin.env$paramChangeFlag<-TRUE
                                 mrbin.env$mrbinparam$reference2D[1]<-as.numeric(regionTMP)
                          }
                          regionTMP<-readline(prompt=paste("New right border, press enter to keep ",
                                    mrbin.env$mrbinparam$reference2D[2],": ",sep=""))
                          if(!regionTMP=="") {
                                 mrbin.env$paramChangeFlag<-TRUE
                                 mrbin.env$mrbinparam$reference2D[2]<-as.numeric(regionTMP)
                          }
                          regionTMP<-readline(prompt=paste("New top border, press enter to keep ",
                                    mrbin.env$mrbinparam$reference2D[3],": ",sep=""))
                          if(!regionTMP=="") {
                                 mrbin.env$paramChangeFlag<-TRUE
                                 mrbin.env$mrbinparam$reference2D[3]<-as.numeric(regionTMP)
                          }
                          regionTMP<-readline(prompt=paste("New bottom border, press enter to keep ",
                                    mrbin.env$mrbinparam$reference2D[4],": ",sep=""))
                          if(!regionTMP=="") {
                                 mrbin.env$paramChangeFlag<-TRUE
                                 mrbin.env$mrbinparam$reference2D[4]<-as.numeric(regionTMP)
                          }
                        }
                      }
                    }
                  }
                }
                if(referenceScaling=="Go back"|adjRegion=="Go back"){
                   selectStep<-selectStep-2
                }
                if(!stopTMP) selectStep<-selectStep+1
              }
            }
            if(selectStep==7){
              #Remove solvent
              if(!stopTMP){
                adjRegion<-""
                removeSolvent<-utils::select.list(c("Yes","No","Go back"),
                                         preselect=mrbin.env$mrbinparam$removeSolvent,
                                         title = "Remove solvent area?",graphics=TRUE)
                if(length(removeSolvent)==0|removeSolvent=="") stopTMP<-TRUE
                if(!stopTMP&!removeSolvent=="Go back"){
                  mrbin.env$mrbinparam$removeSolvent<-removeSolvent
                  if(mrbin.env$mrbinparam$removeSolvent=="Yes"){
                    accept<-FALSE
                    adjRegion<-""
                    while(!accept&!stopTMP&!adjRegion=="Go back"){
                      mean1<-mean(mrbin.env$mrbinparam$solventRegion[1:2])
                      range1<-max(mrbin.env$mrbinparam$solventRegion[1:2])-min(mrbin.env$mrbinparam$solventRegion[1:2])
                      plotNMR(region=c(mean1+8*range1,mean1-8*range1,-10,160),
                              rectangleRegions=matrix(c(mrbin.env$mrbinparam$solventRegion[1],
                                                      mrbin.env$mrbinparam$solventRegion[2],-1000,1000),ncol=4),
                              color="black",
                              manualScale=FALSE,
                              plotTitle=paste("Solvent region\nleft=",mrbin.env$mrbinparam$solventRegion[1],
                                        "ppm, right=",mrbin.env$mrbinparam$solventRegion[2],"ppm",sep=""))
                      SolventRegionTitle<-paste("Keep: left=",mrbin.env$mrbinparam$solventRegion[1],
                                            "ppm, right=",mrbin.env$mrbinparam$solventRegion[2],"ppm",sep="")
                      adjRegion<-utils::select.list(c(SolventRegionTitle,
                                 "Change...","Go back"),
                                 preselect=SolventRegionTitle,title ="Solvent region to be removed: ",graphics=TRUE)
                      if(length(adjRegion)==0|adjRegion=="") stopTMP<-TRUE
                      if(adjRegion==SolventRegionTitle) accept=TRUE
                      if(adjRegion=="Change..."&!stopTMP){
                        regionTMP<-readline(prompt=paste("New left border, press enter to keep ",
                                  mrbin.env$mrbinparam$solventRegion[1],": ",sep=""))
                        if(!regionTMP=="") {
                               mrbin.env$mrbinparam$solventRegion[1]<-as.numeric(regionTMP)
                        }
                        regionTMP<-readline(prompt=paste("New right border, press enter to keep ",
                                  mrbin.env$mrbinparam$solventRegion[2],": ",sep=""))
                        if(!regionTMP=="") {
                               mrbin.env$mrbinparam$solventRegion[2]<-as.numeric(regionTMP)
                        }
                      }
                    }
                  }
                }
                if(removeSolvent=="Go back"|adjRegion=="Go back"){
                   selectStep<-selectStep-2
                }
                if(!stopTMP) selectStep<-selectStep+1
              }
            }
            if(selectStep==8){
            #Remove additional areas
              if(!stopTMP){
                removeAreaListTMP<-""
                removeAreas<-utils::select.list(c("Yes","No","Go back"),preselect=mrbin.env$mrbinparam$removeAreas,
                                         title = "Remove additional areas?",graphics=TRUE)
                if(length(removeAreas)==0|removeAreas=="") stopTMP<-TRUE
                if(!stopTMP&!removeAreas=="Go back"){
                  mrbin.env$mrbinparam$removeAreas<-removeAreas
                  if(mrbin.env$mrbinparam$removeAreas=="Yes"){
                    addAreasFlag<-TRUE
                    if(nrow(mrbin.env$mrbinparam$removeAreaList)>0){
                        addAreasFlag<-FALSE
                        removeAreaListTMP<-utils::select.list(c("Keep","New","Add to existing list","Go back"),
                                           preselect="Keep",
                                           title = "Use previous area list or define new?",graphics=TRUE)
                        if(length(removeAreaListTMP)==0|removeAreaListTMP=="") stopTMP<-TRUE
                        if(!removeAreaListTMP=="Keep"&!stopTMP&!removeAreaListTMP=="Go back"){
                          addAreasFlag<-TRUE
                          if(removeAreaListTMP=="New"&!stopTMP){
                            mrbin.env$mrbinparam$removeAreaList<-matrix(ncol=4,nrow=0,dimnames=list(NULL,c("left","right","top","bottom")))
                          }
                        }
                    }
                    if(!stopTMP){
                      iaddAreas<-nrow(mrbin.env$mrbinparam$removeAreaList)+1
                    }
                    while(addAreasFlag&!stopTMP){
                      mrbin.env$mrbinparam$removeAreaList<-rbind(mrbin.env$mrbinparam$removeAreaList,c(0,0,0,0))
                      mrbin.env$mrbinparam$removeAreaList[iaddAreas,1]<-as.numeric(readline(prompt="Left border: "))
                      mrbin.env$mrbinparam$removeAreaList[iaddAreas,2]<-as.numeric(readline(prompt="Right border: "))
                      if(mrbin.env$mrbinparam$dimension=="2D"&!stopTMP){
                        mrbin.env$mrbinparam$removeAreaList[iaddAreas,3]<-as.numeric(readline(prompt="Top border: "))
                        mrbin.env$mrbinparam$removeAreaList[iaddAreas,4]<-as.numeric(readline(prompt="Bottom border: "))
                      }
                      iaddAreas<-iaddAreas+1
                      keepAdding<-utils::select.list(c("No","Yes"),preselect="No",multiple=FALSE,
                                         title = "Add additional areas?",graphics=TRUE)
                      if(length(keepAdding)==0|keepAdding=="") stopTMP<-TRUE
                      if(keepAdding=="No"&!stopTMP)  addAreasFlag<-FALSE
                    }
                  }
                }
                if(removeAreas=="Go back"|removeAreaListTMP=="Go back"){
                   selectStep<-selectStep-2
                }
                if(!stopTMP) selectStep<-selectStep+1
              }
            }
            if(selectStep==9){
              #Sum up bins of unstable peaks
              if(!stopTMP){
                sumBinListTMP<-""
                sumBins<-utils::select.list(c("Sum bins of unstable peaks (e.g. citrate)","No","Go back"),
                                     preselect=mrbin.env$mrbinparam$sumBins,
                                     title = "Sum bins of unstable peaks?",graphics=TRUE)
                if(length(sumBins)==0|sumBins=="") stopTMP<-TRUE
                if(!stopTMP&!sumBins=="Go back"){
                  if(sumBins=="Sum bins of unstable peaks (e.g. citrate)"){
                    mrbin.env$mrbinparam$sumBins<-"Yes"
                  } else {
                    mrbin.env$mrbinparam$sumBins<-sumBins
                  }
                  if(mrbin.env$mrbinparam$sumBins=="Yes"&!stopTMP){
                      addAreasFlag<-TRUE
                      if(nrow(mrbin.env$mrbinparam$sumBinList)>0&!stopTMP){
                          addAreasFlag<-FALSE
                          sumBinListTMP<-utils::select.list(c("Keep","New","Add to existing list","Go back"),
                                         preselect="Keep",
                                         title = "Use previous area list or define new?",graphics=TRUE)
                          if(length(sumBinListTMP)==0|sumBinListTMP=="") stopTMP<-TRUE
                          if(!sumBinListTMP=="Keep"&!stopTMP&!sumBinListTMP=="Go back"){
                            addAreasFlag<-TRUE
                            if(!sumBinListTMP=="New"&!stopTMP){
                                mrbin.env$mrbinparam$sumBinList<-matrix(ncol=4,nrow=0,dimnames=list(NULL,c("left","right","top","bottom")))
                            }
                          }
                      }
                      if(!stopTMP&!sumBinListTMP=="Go back"){
                        iaddAreas<-nrow(mrbin.env$mrbinparam$sumBinList)+1
                      }
                      while(addAreasFlag&!stopTMP&!sumBinListTMP=="Go back"){
                        mrbin.env$mrbinparam$sumBinList<-rbind(mrbin.env$mrbinparam$sumBinList,c(0,0,0,0))
                        mrbin.env$mrbinparam$sumBinList[iaddAreas,1]<-as.numeric(readline(prompt="Left border: "))
                        mrbin.env$mrbinparam$sumBinList[iaddAreas,2]<-as.numeric(readline(prompt="Right border: "))
                        if(mrbin.env$mrbinparam$dimension=="2D"&!stopTMP){
                            mrbin.env$mrbinparam$sumBinList[iaddAreas,3]<-as.numeric(readline(prompt="Top border: "))
                            mrbin.env$mrbinparam$sumBinList[iaddAreas,4]<-as.numeric(readline(prompt="Bottom border: "))
                        }
                        iaddAreas<-iaddAreas+1
                        keepAdding<-utils::select.list(c("No","Yes"),preselect="No",title = "Add additional areas?",graphics=TRUE)
                        if(length(keepAdding)==0|keepAdding=="") stopTMP<-TRUE
                        if(keepAdding=="No"&!stopTMP)  addAreasFlag<-FALSE
                      }
                  }
                }
                if(sumBins=="Go back"|sumBinListTMP=="Go back"){
                   selectStep<-selectStep-2
                }
                if(!stopTMP) selectStep<-selectStep+1
              }
            }
            if(selectStep==10){
              #Remove noise
              if(!stopTMP){
                SNRTMP<-""
                noiseTMP<-""
                noiseRemoval<-utils::select.list(c("Yes","No","Go back"),
                                         preselect=mrbin.env$mrbinparam$noiseRemoval,
                                         title="Remove noise?",graphics=TRUE)
                if(length(noiseRemoval)==0|noiseRemoval=="") stopTMP<-TRUE
                if(!stopTMP&!noiseRemoval=="Go back"){
                  mrbin.env$mrbinparam$noiseRemoval<-noiseRemoval
                  if(mrbin.env$mrbinparam$noiseRemoval=="Yes"&!stopTMP){
                    if(mrbin.env$mrbinparam$dimension=="1D"&!stopTMP){
                      SNRTMP<-utils::select.list(c(mrbin.env$mrbinparam$signal_to_noise1D,"Change...","Go back"),
                              preselect=
                              as.character(mrbin.env$mrbinparam$signal_to_noise1D),title="Signal-to-noise ratio (SNR):",graphics=TRUE)
                      if(length(SNRTMP)==0|SNRTMP=="") stopTMP<-TRUE
                      if(SNRTMP=="Change...") SNRTMP<-readline(prompt=paste(
                                           "New 1D signal to noise ratio, press enter to keep ",
                                           mrbin.env$mrbinparam$signal_to_noise1D,": ",sep=""))
                      if(!SNRTMP==""&!stopTMP&!SNRTMP=="Go back") {
                          mrbin.env$mrbinparam$signal_to_noise1D<-as.numeric(SNRTMP)
                      }
                    }
                    if(mrbin.env$mrbinparam$dimension=="2D"&!stopTMP){
                      SNRTMP<-utils::select.list(c(mrbin.env$mrbinparam$signal_to_noise2D,"Change...","Go back"),
                              preselect=
                              as.character(mrbin.env$mrbinparam$signal_to_noise2D),title="Signal-to-noise ratio (SNR):",graphics=TRUE)
                      if(length(SNRTMP)==0|SNRTMP=="") stopTMP<-TRUE
                      if(SNRTMP=="Change...") SNRTMP<-readline(prompt=paste(
                                           "New 2D signal to noise ratio, press enter to keep ",
                                           mrbin.env$mrbinparam$signal_to_noise2D,": ",sep=""))
                      if(!SNRTMP==""&!stopTMP&!SNRTMP=="Go back") {
                           mrbin.env$mrbinparam$signal_to_noise2D<-as.numeric(SNRTMP)
                      }
                    }
                    if(!stopTMP&!SNRTMP=="Go back"){
                      noiseTMP<-utils::select.list(unique(c(as.character(mrbin.env$mrbinparam$noiseThreshold),
                                             "0.2","0.75","0.05","Change...","Go back")),
                                             preselect=as.character(mrbin.env$mrbinparam$noiseThreshold),
                                             title="Minimum percentage > noise",graphics=TRUE)
                      if(length(noiseTMP)==0|noiseTMP=="") stopTMP<-TRUE
                      if(!stopTMP&!noiseTMP=="Go back"){
                        if(noiseTMP=="Change..."&!stopTMP){
                              noiseTMP<-readline(prompt=paste("New noise threshold, press enter to keep ",
                                        mrbin.env$mrbinparam$noiseThreshold,": ",sep=""))
                        }
                        if(!noiseTMP==""&!stopTMP&!noiseTMP=="Go back")  {
                          mrbin.env$mrbinparam$noiseThreshold<-as.numeric(noiseTMP)
                        }
                      }
                    }
                  }
                }
                if(noiseRemoval=="Go back"|SNRTMP=="Go back"|noiseTMP=="Go back"){
                   selectStep<-selectStep-2
                }
                if(!stopTMP) selectStep<-selectStep+1
              }
            }
            if(selectStep==11){
              #Crop HSQCs
              if(mrbin.env$mrbinparam$dimension=="2D"&!stopTMP){
                 plotNMR(region="all",
                        polygonRegion=matrix(c(mrbin.env$mrbinparam$croptopRight,
                                      mrbin.env$mrbinparam$croptopLeft,
                                      mrbin.env$mrbinparam$cropbottomLeft,
                                      mrbin.env$mrbinparam$cropbottomRight),
                                      ncol=2,byrow=TRUE),
                        color="black",
                        manualScale=FALSE,
                        plotTitle=paste("Crop spectrum to diagonal",sep=""))
                 cropHSQC<-utils::select.list(c("Yes","No","Go back"),
                                       preselect=mrbin.env$mrbinparam$cropHSQC,
                                       title="Crop HSQCs?",graphics=TRUE)
                 if(length(cropHSQC)==0|cropHSQC=="") stopTMP<-TRUE
                 if(!stopTMP&!cropHSQC=="Go back"){
                   mrbin.env$mrbinparam$cropHSQC<-cropHSQC
                 }
                if(cropHSQC=="Go back"){
                   selectStep<-selectStep-2
                }
              }
              if(!stopTMP) selectStep<-selectStep+1
            }
            if(selectStep==12){
              #PQN scaling?
              if(!stopTMP){
                PQNScaling<-utils::select.list(c("Yes","No","Go back"),
                                        preselect=mrbin.env$mrbinparam$PQNScaling,
                                        title = "PQN normalization?",graphics=TRUE)
                if(length(PQNScaling)==0|PQNScaling=="") stopTMP<-TRUE
                if(!stopTMP&!PQNScaling=="Go back"){
                  mrbin.env$mrbinparam$PQNScaling<-PQNScaling
                }
                if(!stopTMP&PQNScaling=="Go back"){
                   selectStep<-selectStep-2
                }
                if(!stopTMP) selectStep<-selectStep+1
              }
            }
            if(selectStep==13){
              #Replace negative values
              if(!stopTMP){
                fixNegatives<-utils::select.list(c("Yes","No","Go back"),
                                          preselect=mrbin.env$mrbinparam$fixNegatives,
                                          title="Replace negative values?",graphics=TRUE)
                if(length(fixNegatives)==0|fixNegatives=="") stopTMP<-TRUE
                if(!stopTMP&!fixNegatives=="Go back"){
                  mrbin.env$mrbinparam$fixNegatives<-fixNegatives
                }
                if(!stopTMP&fixNegatives=="Go back"){
                   selectStep<-selectStep-2
                }
                if(!stopTMP) selectStep<-selectStep+1
              }
            }
            if(selectStep==14){
              #Log scaling?
              if(!stopTMP){
                logTrafo<-utils::select.list(c("Yes","No","Go back"),
                                      preselect=mrbin.env$mrbinparam$logTrafo,
                                      title="Log transformation?",graphics=TRUE)
                if(length(logTrafo)==0|logTrafo=="") stopTMP<-TRUE
                if(!stopTMP&!logTrafo=="Go back"){
                  mrbin.env$mrbinparam$logTrafo<- logTrafo
                }
                if(!stopTMP&logTrafo=="Go back"){
                   selectStep<-selectStep-2
                }
                if(!stopTMP) selectStep<-selectStep+1
              }
            }
            if(selectStep==15){
              #Define sample names
              if(!stopTMP){
                NamesDictTMP<-c("Folder names","Spectrum titles","Folder names and EXPNO")
                names(NamesDictTMP)<-paste(NamesDictTMP," (\"",c(mrbin.env$mrbinTMP$currentSpectrumFolderName,
                                     mrbin.env$mrbinTMP$currentSpectrumTitle,
                                     mrbin.env$mrbinTMP$currentSpectrumFolderName_EXPNO),
                                     "\", ...)",sep="")
                NamesDictTMP2<-names(NamesDictTMP)#paste(NamesDictTMP," (",c(mrbin.env$mrbinTMP$currentSpectrumFolderName,
                               #      mrbin.env$mrbinTMP$currentSpectrumTitle,
                               #      mrbin.env$mrbinTMP$currentSpectrumFolderName_EXPNO),")",sep="")
                names(NamesDictTMP2)<-NamesDictTMP
                useAsNames<-utils::select.list(c(names(NamesDictTMP),"Go back"),
                                          preselect=NamesDictTMP2[mrbin.env$mrbinparam$useAsNames],
                                          title = "Create sample names from",graphics=TRUE)
                if(length(useAsNames)==0|useAsNames=="") stopTMP<-TRUE
                if(!stopTMP&!useAsNames=="Go back"){
                  mrbin.env$mrbinparam$useAsNames<-NamesDictTMP[useAsNames]
                }
                if(!stopTMP&useAsNames=="Go back"){
                   selectStep<-selectStep-2
                }
                if(!stopTMP) selectStep<-selectStep+1
              }
            }
            if(selectStep==16){
              #Plot results
              if(!stopTMP){
                PCAtitlelength<-""
                PCA<-utils::select.list(c("Yes","No","Go back"),
                                  preselect = mrbin.env$mrbinparam$PCA,
                                  title = "Create result plot?",graphics=TRUE)
                if(length(PCA)==0|PCA=="") stopTMP<-TRUE
                if(!stopTMP&!PCA=="Go back"){
                  mrbin.env$mrbinparam$PCA<-PCA
                  if(!stopTMP&mrbin.env$mrbinparam$PCA=="Yes"){
                    currentPCAtitlelength<-as.character(mrbin.env$mrbinparam$PCAtitlelength)
                    if(mrbin.env$mrbinparam$useAsNames=="Spectrum titles") Title<-mrbin.env$mrbinTMP$currentSpectrumTitle
                    if(mrbin.env$mrbinparam$useAsNames=="Folder names") Title<-mrbin.env$mrbinTMP$currentSpectrumFolderName
                    if(mrbin.env$mrbinparam$useAsNames=="Folder names and EXPNO") Title<-mrbin.env$mrbinTMP$currentSpectrumFolderName_EXPNO
                    TitleListTMP<-unique(c(4,6,8,500,currentPCAtitlelength))
                    names(TitleListTMP)<-unique(c(4,6,8,500,currentPCAtitlelength))
                    TitleListTMPDict<-as.character(TitleListTMP)
                    names(TitleListTMPDict)<-TitleListTMP
                    TitleListTMPDict[TitleListTMPDict=="500"] <-"All"
                    TitleListTMP2<-c(TitleListTMPDict,"Custom...","Go back")
                    names_TitleListTMP2<-NULL
                    for(i_TitleListTMP2 in 1:length(TitleListTMPDict)){
                      names_TitleListTMP2<-c(names_TitleListTMP2,paste(TitleListTMPDict[i_TitleListTMP2]," letters (\"",
                                           substr(Title,1,as.numeric(TitleListTMP[i_TitleListTMP2])),"\", ...)",sep=""))
                    }
                    names_TitleListTMP2<-c(names_TitleListTMP2,"Custom...","Go back")
                    names(TitleListTMP2)<-names_TitleListTMP2
                    TitleListTMP4<-TitleListTMP2
                    TitleListTMP4[TitleListTMP4=="All"]<-"500"
                    TitleListTMP3<-c(paste(TitleListTMPDict," (",substr(Title,1,TitleListTMP),")",sep=""),"Custom...","Go back")
                    names(TitleListTMP3)<-TitleListTMP2
                    PCAtitlelength<-utils::select.list(names(TitleListTMP2),
                                      preselect = names(TitleListTMP2)[TitleListTMP4==as.character(mrbin.env$mrbinparam$PCAtitlelength)],
                                      title = "Crop titles for plot?",graphics=TRUE)
                    if(length(PCAtitlelength)==0|PCAtitlelength=="") stopTMP<-TRUE
                    if(!stopTMP&!PCAtitlelength=="Go back"){
                      if(TitleListTMP2[PCAtitlelength]=="All"){
                        mrbin.env$mrbinparam$PCAtitlelength<-500
                      } else {
                        if(PCAtitlelength=="Custom..."){
                            PCAtitlelengthTMP<-readline(prompt=paste("New title length, press enter to keep ",
                                mrbin.env$mrbinparam$PCAtitlelength,": ",sep=""))
                            if(!PCAtitlelengthTMP=="") {
                                mrbin.env$mrbinparam$PCAtitlelength<-as.numeric(PCAtitlelengthTMP)
                            } else {
                                #PCAtitlelength<-as.character(mrbin.env$mrbinparam$PCAtitlelength)
                            }
                        } else {
                           mrbin.env$mrbinparam$PCAtitlelength<-as.numeric(TitleListTMP2[PCAtitlelength])
                        }
                      }
                    }
                  }
                }
                if(!stopTMP&(PCA=="Go back"|PCAtitlelength=="Go back")){
                   selectStep<-selectStep-2
                }
                if(!stopTMP) selectStep<-selectStep+1
              }
            }
            if(selectStep==17){
              #Define groups
              if(!stopTMP){
                selectionFactors<-""
                defineGroupsTMP<-mrbin.env$mrbinparam$defineGroups
                if(defineGroupsTMP=="Create groups for PCA plot") defineGroupsTMP<-"Yes"
                defineGroups<-utils::select.list(c("Create groups for PCA plot","No","Go back"),
                                          preselect=defineGroupsTMP,
                                          title = "Define group members for PCA?",graphics=TRUE)
                if(length(defineGroups)==0|defineGroups=="") stopTMP<-TRUE
                if(!stopTMP&!defineGroups=="Go back"){
                  if(defineGroups=="Create groups for PCA plot"){
                    mrbin.env$mrbinparam$defineGroups<-"Yes"
                  } else {
                    mrbin.env$mrbinparam$defineGroups<-defineGroups
                  }
                  if(mrbin.env$mrbinparam$defineGroups=="Yes"){
                    if(!is.null(mrbin.env$mrbinparam$Factors)){
                       if(length(mrbin.env$mrbinparam$Factors)==length(mrbin.env$mrbinparam$NMRfolders)){
                         yesOptionTMP<-paste("Use previous factor list (",
                                             paste(mrbin.env$mrbinparam$Factors[1:max(3,length(mrbin.env$mrbinparam$Factors))],
                                                   sep=", ",collapse=", "),
                                             ", ...)",sep="")
                         selectionFactors<-utils::select.list(c(yesOptionTMP,"No","Go back"),
                                           preselect=yesOptionTMP,
                                           title = "Use previous factor list?",graphics=TRUE)
                         if(length(selectionFactors)==0|selectionFactors=="") stopTMP<-TRUE
                         if(!stopTMP&selectionFactors=="No")  setFactors()
                       } else {
                         setFactors()
                       }
                    } else {
                      setFactors()
                    }
                  }
                }
                if(!stopTMP&(defineGroups=="Go back"|selectionFactors=="Go back")){
                   selectStep<-selectStep-2
                }
                if(!stopTMP) selectStep<-selectStep+1
              }
            }
     }
     if(selectStep==18){
       #Save output files to hard drive?
       if(!stopTMP){
         saveFilesTMP<-utils::select.list(c("Yes","No","Go back"),
                     preselect=mrbin.env$mrbinparam$saveFiles,
                     title ="Save output to disk?",graphics=TRUE)
         if(length(saveFilesTMP)==0|saveFilesTMP=="") stopTMP<-TRUE
         if(!stopTMP&!saveFilesTMP=="Go back"){
          mrbin.env$mrbinparam$saveFiles<-saveFilesTMP
          if(mrbin.env$mrbinparam$saveFiles=="Yes"&!stopTMP){
            enterFoldersTMP<-readline(prompt="Enter starting folder path. (Examples: Windows: \"C:\\\", Apple: \"/\") : ")
            if(enterFoldersTMP=="") stopTMP<-TRUE
            if(!stopTMP){
             parentFolder<-gsub('\\\\',"/",enterFoldersTMP)
             #Browse
             selectFlag<-0
             while(selectFlag<1){
               folderListFull<-list.dirs(path=parentFolder,recursive = FALSE,full.names=TRUE)
               folderList<-list.dirs(path=parentFolder,recursive = FALSE,full.names=FALSE)
               if(length(strsplit(parentFolder,split="/")[[1]])>1){
                  folderListTMP<-c("..                                                                            ",
                                  folderList)
               } else {
                  folderListTMP<-folderList
               }
               selectFolders<-utils::select.list(folderListTMP,preselect=NULL,multiple=TRUE,
                                title = "Go to folder, then click OK",graphics=TRUE)
               if(!length(selectFolders)==1) {
                  selectFolders<-parentFolder
                  selectList<-parentFolder
                  selectFlag<-1
               }
               if(length(selectFolders)==1&selectFlag<1){
                    if(selectFolders=="..                                                                            "){
                         if(length(strsplit(parentFolder,split="/")[[1]])>1){
                            selectList<-paste(rev(rev(strsplit(parentFolder,split="/")[[1]])[-1]),sep="",collapse="/")
                         } else {
                            selectList<-parentFolder
                         }
                    } else {
                         selectList<-folderListFull[which(folderList%in%selectFolders)]
                    }
               }
               parentFolder<-selectList
             }
             filenameTMP<-utils::select.list(c(paste("mrbin_",gsub(":","-",gsub(" ","_",Sys.time())),
                         sep=""),"Change..."),
                         title ="Output file name: ",graphics=TRUE)
             if(length(filenameTMP)==0|filenameTMP=="") stopTMP<-TRUE
             if(!stopTMP){
              if(filenameTMP=="Change..."&!stopTMP){
                filenameTMP<-readline(prompt=paste("New file name, press enter to use ",
                          paste("mrbin_",gsub(":","-",gsub(" ","_",Sys.time())),sep=""),": \n",sep=""))
                #if(!filenameTMP=="") filenameTMP
                if(filenameTMP=="") filenameTMP<-paste("mrbin_",gsub(":","-",gsub(" ","_",Sys.time())),
                           sep="")
               }
               mrbin.env$mrbinparam$outputFileName<-gsub("//","/",paste(
                                parentFolder,"/",filenameTMP,sep=""))
             }
            }
           }
        }
        if(!stopTMP&saveFilesTMP=="Go back"){
           selectStep<-selectStep-2
        }
        if(!stopTMP) selectStep<-selectStep+1
      }
     }
     if(selectStep==19){
       if(!stopTMP){
         keepDataTMP<-""
         startmrbin<-utils::select.list(c("Start binning now","I'll do it later","Go back"),
                    preselect="Start binning now",
                    title = "You're all set. Start binning?",graphics=TRUE)
         if(length(startmrbin)==0|startmrbin=="") stopTMP<-TRUE
         if(!stopTMP&!startmrbin=="Go back"){
           if(startmrbin=="Start binning now"&!stopTMP){
             if(!is.null(mrbin.env$mrbinTMP$binsRaw)){
               if(!mrbin.env$paramChangeFlag){
                 if(nrow(mrbin.env$mrbinTMP$binsRaw)==length(mrbin.env$mrbinparam$NMRfolders)) {
                   createBinNumbers()
                   if(mrbin.env$mrbinTMP$nbins==ncol(mrbin.env$mrbinTMP$binsRaw)){
                     keepDataTMP<-utils::select.list(c("Calculate new bin data (recommended)",
                                  "Use existing bin data (may cause inconsistent data)",
                                  "Go back"),
                                  preselect="Create new bin data (recommended)",title =
                                  "Previous data available",graphics=TRUE)
                     if(length(keepDataTMP)==0|keepDataTMP=="") stopTMP<-TRUE
                     if(keepDataTMP=="Use existing bin data (may cause inconsistent data)"&!stopTMP){
                          mrbin.env$bins<-mrbin.env$mrbinTMP$binsRaw
                          mrbin.env$mrbinparam$createBins<-"No"
                          createBinRegions()
                     }
                   }
                 }
               }
             }
           }
         }
         if(!stopTMP&(startmrbin=="Go back"|keepDataTMP=="Go back")){
           selectStep<-selectStep-2
         }
         if(!stopTMP) selectStep<-selectStep+1
       }
     }
     if(selectStep==20){
       if(!stopTMP)  lastStepDone<-TRUE
     }
    }
   }
 }
 if(!stopTMP){
   if(startmrbin=="Start binning now"){
     mrbinrun()
     invisible(list(bins=mrbin.env$bins,factors=mrbin.env$mrbinparam$Factors,parameters=mrbin.env$mrbinparam))
   }
 }
}

#' A function performing all data read and processing steps.
#'
#' This function reads parameters from the global variable mrbin.env$mrbinparam and
#' performs the following operations:
#' Reading NMR files, creating bins, removing solvent area, removing additional
#' user-defined areas, summing up bins that contain unstable peaks such as
#' citric acid, removes noise bins, crops HSQC spectra to the diagonal area,
#' performs PQN scaling, replaces negative values, log transforms and displays a
#' PCA plot. Parameters are then saved in a text file. These can be recreated
#' using recreatemrbin().
#' @return {None}
#' @export
#' @examples
#' setParam(parameters=list(dimension="2D",binwidth2D=0.05,binheight=3,PQNScaling="No",
#'          fixNegatives="No",logTrafo="No",
#'          NMRfolders=c(system.file("extdata/1/12/pdata/10",package="mrbin"),
#'                       system.file("extdata/2/12/pdata/10",package="mrbin"))))
#' mrbinrun()

mrbinrun<-function(){
  if(!exists("mrbin.env", mode="environment")) .onLoad()
  if(!is.null(mrbin.env$mrbinparam$NMRfolders)){
    #if(!"mrbinparam"%in%ls(envir = mrbin.env))    resetEnv()
    if(mrbin.env$mrbinparam$createBins=="Yes") createBinRegions()
    if(mrbin.env$mrbinparam$removeSolvent=="Yes") removeSolvent()
    if(mrbin.env$mrbinparam$removeAreas=="Yes") removeAreas()
    if(mrbin.env$mrbinparam$sumBins=="Yes") sumBins()
    if(mrbin.env$mrbinparam$cropHSQC=="Yes"&mrbin.env$mrbinparam$dimension=="2D") cropNMR()
    if(mrbin.env$mrbinparam$createBins=="Yes") binMultiNMR()
    if(mrbin.env$mrbinparam$noiseRemoval=="Yes") removeNoise()
    if(mrbin.env$mrbinparam$createBins=="Yes") createBinNames()
    if(mrbin.env$mrbinparam$PQNScaling=="Yes") PQNScaling()
    if(mrbin.env$mrbinparam$fixNegatives=="Yes") atnv()
    if(mrbin.env$mrbinparam$logTrafo=="Yes") logTrafo()
    if(mrbin.env$mrbinparam$saveFiles=="Yes"){
      dput(mrbin.env$mrbinparam, file = paste(mrbin.env$mrbinparam$outputFileName,".txt",sep=""))
      utils::write.csv(mrbin.env$bins, file = paste(mrbin.env$mrbinparam$outputFileName,"bins.csv",sep=""))
    }
    if(mrbin.env$mrbinparam$PCA=="Yes"){
      plotResults()
      if(mrbin.env$mrbinparam$saveFiles=="Yes"){
        grDevices::dev.copy(pdf,paste(mrbin.env$mrbinparam$outputFileName,"plot.pdf",sep=""))
        grDevices::dev.off()
      }
    }
     mrbin.env$paramChangeFlag<-FALSE
     mrbin.env$mrbinparam$createBins<-"Yes"
     resultOutputTMP<-c("\nNumber of spectra: ",nrow(mrbin.env$bins),"\n",
         "Number of bins at start: ",mrbin.env$mrbinparam$numberOfFeaturesRaw,"\n")
     if(!is.null(mrbin.env$mrbinparam$numberOfFeaturesAfterRemovingSolvent)&mrbin.env$mrbinparam$removeSolvent=="Yes"){
          resultOutputTMP<-c(resultOutputTMP,
             "Number of bins after removing solvent: ",mrbin.env$mrbinparam$numberOfFeaturesAfterRemovingSolvent,"\n")
     }
     if(!is.null(mrbin.env$mrbinparam$numberOfFeaturesAfterRemovingAreas)&mrbin.env$mrbinparam$removeAreas=="Yes"){
          resultOutputTMP<-c(resultOutputTMP,
             "Number of bins after removing areas: ",mrbin.env$mrbinparam$numberOfFeaturesAfterRemovingAreas,"\n")
     }
     if(!is.null(mrbin.env$mrbinparam$numberOfFeaturesAfterSummingBins)&mrbin.env$mrbinparam$sumBins=="Yes"){
          resultOutputTMP<-c(resultOutputTMP,
             "Number of bins after summing bins: ",mrbin.env$mrbinparam$numberOfFeaturesAfterSummingBins,"\n")
     }
     if(!is.null(mrbin.env$mrbinparam$numberOfFeaturesAfterCropping)&mrbin.env$mrbinparam$cropHSQC=="Yes"&mrbin.env$mrbinparam$dimension=="2D"){
          resultOutputTMP<-c(resultOutputTMP,
             "Number of bins after cropping: ",mrbin.env$mrbinparam$numberOfFeaturesAfterCropping,"\n")
     }
     if(!is.null(mrbin.env$mrbinparam$numberOfFeaturesAfterNoiseRemoval)&mrbin.env$mrbinparam$noiseRemoval=="Yes"){
          resultOutputTMP<-c(resultOutputTMP,
             "Number of bins after noise removal: ",mrbin.env$mrbinparam$numberOfFeaturesAfterNoiseRemoval,"\n")
     }
     resultOutputTMP<-paste(resultOutputTMP,sep="")
     if(mrbin.env$mrbinparam$verbose){
       printParameters()
       message(resultOutputTMP)
     }
  }
}

#' A function for printing parameters to the screen.
#'
#' This function reads parameters from the global variable mrbin.env$mrbinparam and
#' prints the required R code for creating a data set to the screen.
#' @return {None}
#' @export
#' @examples
#' printParameters()

printParameters<-function(){
  if(!exists("mrbin.env", mode="environment")) .onLoad()
  if(mrbin.env$mrbinparam$verbose){
    if(mrbin.env$mrbinparam$dimension=="1D") paramNamesTMP<-mrbin.env$requiredParam1D
    if(mrbin.env$mrbinparam$dimension=="2D") paramNamesTMP<-mrbin.env$requiredParam2D
      printTMPline<-paste("results<-mrbin(silent=FALSE,parameters=list(")
      printTMP<-NULL
      counter<-0
      for(i in paramNamesTMP){
        counter<-counter+1
        if(is.character(mrbin.env$mrbinparam[[i]])){
          sepSymbol<-"\""
        } else {
          sepSymbol<-""
        }
        returnSymbol<-"\n  "
        #if(length(mrbin.env$mrbinparam[[i]])<17&!i=="NMRfolders"){
        if(!i=="NMRfolders"){
          returnSymbol<-""
        } #else {
        #}
        vectorSymbol1<-""
        vectorSymbol2<-""
        if(is.matrix(mrbin.env$mrbinparam[[i]])){
          if(nrow(mrbin.env$mrbinparam[[i]])==0){
            vectorSymbol1<-"matrix(nrow=0"
            vectorSymbol2<-paste(",ncol=",ncol(mrbin.env$mrbinparam[[i]]),")",sep="")
          } else {
            vectorSymbol1<-"matrix(c("
            vectorSymbol2<-paste("),ncol=",ncol(mrbin.env$mrbinparam[[i]]),")",sep="")
          }
        }
        if(length(mrbin.env$mrbinparam[[i]])>1){
          if(is.vector(mrbin.env$mrbinparam[[i]])){
            vectorSymbol1<-"c("
            vectorSymbol2<-")"
          }
          if(is.factor(mrbin.env$mrbinparam[[i]])){
            vectorSymbol1<-"factor(c("
            vectorSymbol2<-"))"
            sepSymbol<-"\""
          }
        }
        sepTMP<-paste(sepSymbol,",",returnSymbol,sepSymbol,sep="")
        if(is.null(mrbin.env$mrbinparam[[i]])){
            valueTMP<-"NULL"
        } else {
          if(is.factor(mrbin.env$mrbinparam[[i]])){
            valueTMP<-paste(as.character(mrbin.env$mrbinparam[[i]]),sep=sepTMP,collapse=sepTMP)
          } else {
            valueTMP<-paste(mrbin.env$mrbinparam[[i]],sep=sepTMP,collapse=sepTMP)
          }
        }
        if(counter<length(paramNamesTMP)){
          TMPline<-paste(i,"=",vectorSymbol1,sepSymbol,valueTMP,
                      sepSymbol,vectorSymbol2,",",sep="")
          if((nchar(printTMPline)+nchar(TMPline))<=79){
            printTMPline<-paste(printTMPline,TMPline,sep="")
          } else {
             printTMP<-paste(printTMP,printTMPline,"\n ",sep="")
              printTMPline<-TMPline
          }
        } else {
          TMPline<-paste(i,"=",vectorSymbol1,sepSymbol,valueTMP,
                      sepSymbol,vectorSymbol2,"",sep="")
          if((nchar(printTMPline)+nchar(TMPline))<=77){
            printTMP<-paste(printTMP,printTMPline,TMPline,"))\n",sep="")
          } else {
             printTMP<-paste(printTMP,printTMPline,"\n ",sep="")
             printTMP<-paste(printTMP,TMPline,"))\n",sep="")
          }

        }
      }
      #printTMP<-paste(printTMP,"))\n",sep="")
      cat("\nTo recreate this data set, use the following code:\n\n##################################################\n\n")
      cat(printTMP)
      cat("\n##################################################\n\nTo recreate this data set, use the code above this line.\n")
  }
}


#' A function recreating parameters from previous runs.
#'
#' This function reads parameters from a text file that was created during a
#' previous run or mrbin(). After reading, the data can be recreated using
#' mrbin(). File names in mrbin$param might need to be updated.
#' using recreatemrbin().
#' @param filename File path/name of the mrbin parameter file to be loaded
#' @return {None}
#' @export
#' @examples
#' ## Lets the user browse for the parameter file
#' \donttest{ recreatemrbin() }
#' # Insert full folder path and file name
#' recreatemrbin(system.file("extdata/mrbin.txt",package="mrbin"))

recreatemrbin<-function(filename=NULL){
  if(!exists("mrbin.env", mode="environment")) .onLoad()
  if(is.null(filename)){
   selectFlag<-0
   enterFolders<-utils::select.list(c("Browse...","Enter full file path manually"),
                 preselect="Browse...",title="Select parameter file:",graphics=TRUE)
   enterFoldersTMP<-readline(prompt="Enter file path: ")
   if(!enterFoldersTMP==""){
    parentFolder<-gsub('\\\\',"/",enterFoldersTMP)
    if(enterFolders=="Enter file path manually"){
        selectFolders<-parentFolder
        selectList<-parentFolder
        selectFlag<-1
    }
    while(selectFlag<1){
     folderListFull1<-list.dirs(path=parentFolder,recursive = FALSE,full.names=TRUE)
     folderList1<-list.dirs(path=parentFolder,recursive = FALSE,full.names=FALSE)
     folderListFull<-c(folderListFull1,list.files(path = parentFolder, pattern = ".txt",recursive = FALSE,full.names=TRUE))
     folderList<-c(folderList1,list.files(path = parentFolder, pattern = ".txt",recursive = FALSE,full.names=FALSE))
     if(length(strsplit(parentFolder,split="/")[[1]])>1){
        folderListTMP<-c("..                                                                            ",folderList)
     } else {
        folderListTMP<-folderList
     }
     selectFolders<-utils::select.list(folderListTMP,preselect=NULL,multiple=FALSE,
                      title = paste("Browse to parameter file:",parentFolder),graphics=TRUE)
     if(selectFolders%in%c("..                                                                            ",folderList1)){
          if(selectFolders=="..                                                                            "){
               if(length(strsplit(parentFolder,split="/")[[1]])>1){
                  parentFolder<-paste(rev(rev(strsplit(parentFolder,split="/")[[1]])[-1]),sep="",collapse="/")
               }
          } else {
               parentFolder<-folderListFull[which(folderList%in%selectFolders)]
          }
      } else {
           selectFlag<-1
           selectList<-folderListFull[which(folderList%in%selectFolders)]
      }
     }
     filename<-selectList
  }
 }
 if(!is.null(filename)){
    setParam(dget(filename))
 } else {
    warning("No parameter file name specified. No file loaded.")
 }
}

#' A function setting parameters and checking for consistency.
#'
#' This function set parameters and checks parameters for consistency.
#' @param parameters List of parameters to be set
#' @return {None}
#' @export
#' @examples
#' setParam(parameters=list(dimension="1D"))

setParam<-function(parameters=NULL){
  if(!is.null(parameters)){
    mrbin.env$mrbinparam_copy<-mrbin.env$mrbinparam
    if("dimension"%in%names(parameters)){
      dimTMP<-parameters$dimension
    } else {
      dimTMP<-mrbin.env$mrbinparam$dimension
    }
    if(dimTMP=="1D") diffSet3<-setdiff(mrbin.env$requiredParam1D,names(parameters))
    if(dimTMP=="2D") diffSet3<-setdiff(mrbin.env$requiredParam2D,names(parameters))
    diffSet2<-setdiff(names(parameters),names(mrbin.env$mrbinparam_copy))
    intersectSet<-intersect(names(parameters),names(mrbin.env$mrbinparam_copy))
    if(length(diffSet2)>0){
       warning(paste("Unexpected parameters: ",
           paste(diffSet2,sep=", ", collapse=", "),"\n",
           "These parameters are not used. Potentially they were created in a different mrbin version.",
           sep=""))
    }
    if(length(diffSet3)>0){
       if(mrbin.env$mrbinparam$verbose) message(paste("Current values are used for missing parameters: ",
           paste(diffSet3,sep=", ", collapse=", "),sep=""))
    }
    if(length(intersectSet)>0){
       for(iintersectSet in intersectSet){
            mrbin.env$mrbinparam[[iintersectSet]]<-parameters[[iintersectSet]]
       }
    }
    if(!mrbin.env$mrbinTMP$mrbinversion==mrbin.env$mrbinparam$mrbinversionTracker){
       warning(paste("Imported file was created using another mrbin version: ",mrbin.env$mrbinparam$mrbinversionTracker,
           ".\n For exact reproduction of results, please get the old version at: kleinomicslab.com",
           sep=""))
    }
  } else {
      stop("setParam: No parameters were provided.")
  }
}

#' A function for log transforming data.
#'
#' This function simply log transforms. Will not work with negative data.
#' @return {None}
#' @export
#' @examples
#' mrbinExample<-mrbin(silent=TRUE,setDefault=TRUE,parameters=list(dimension="1D", logTrafo="No",
#'                     binwidth1D=0.05,signal_to_noise1D=50, verbose=FALSE, PCA="No",
#'                     NMRfolders=c(system.file("extdata/1/10/pdata/10",package="mrbin"),
#'                                 system.file("extdata/2/10/pdata/10",package="mrbin"))))
#' logTrafo()

logTrafo<-function(){
 if(!is.null(mrbin.env$bins)){
     if(sum(mrbin.env$bins<=0)>0){
         stop("Log transform does not work with negative values.")
     } else {
        if(nrow(mrbin.env$bins)==1){
          rownamesTMP<-rownames(mrbin.env$bins)
          colnamesTMP<-colnames(mrbin.env$bins)
          mrbin.env$bins<-matrix(log(mrbin.env$bins),nrow=1)
          rownames(mrbin.env$bins)<-rownamesTMP
          colnames(mrbin.env$bins)<-colnamesTMP
        } else {
         mrbin.env$bins<-log(mrbin.env$bins)
        }
     }
 }
}

#' A function for interactively setting the current spectrum.
#'
#' This function lets the user pick a spectrum from the list of spectra
#' analysis. This function is meant only for use within the mrbin function.
#' @return {None}
#' @export
#' @examples
#' \donttest{ setCurrentSpectrum() }

setCurrentSpectrum<-function(){
       newCurrent<-utils::select.list(mrbin.env$mrbinparam$NMRfolders,preselect=mrbin.env$mrbinTMP$currentFolder,
                                          title="Select new current spectrum.",graphics=TRUE)
       if(!newCurrent==""){
            mrbin.env$mrbinTMP$currentFolder<-newCurrent
            readNMR()
       }
}

#' A function for removing a spectrum.
#'
#' This function lets the user pick spectra from a list for removal from data
#' analysis. This function is meant only for use within the mrbin function.
#' @return {None}
#' @export
#' @examples
#' \donttest{ removeSpectrum() }

removeSpectrum<-function(){
 if(!is.null(mrbin.env$mrbinparam$NMRfolders)){
    listTMP<-utils::select.list(mrbin.env$mrbinparam$NMRfolders,preselect = NULL, multiple = TRUE,title ="Select spectra to be removed",graphics=TRUE)
    if(length(listTMP)>0){
      if(!listTMP==""){
         if(length(mrbin.env$mrbinparam$Factors)==length(mrbin.env$mrbinparam$NMRfolders)){
           mrbin.env$mrbinparam$Factors<-mrbin.env$mrbinparam$Factors[-which(mrbin.env$mrbinparam$NMRfolders%in%listTMP)]
         }
         if(!is.null(mrbin.env$bins)){
           if(nrow(mrbin.env$bins)==length(mrbin.env$mrbinparam$NMRfolders)){
            if((nrow(mrbin.env$bins)-length(listTMP))==1){
              rownamesTMP<-rownames(mrbin.env$bins)[-which(rownames(mrbin.env$bins)%in%listTMP)]
              colnamesTMP<-colnames(mrbin.env$bins)
              mrbin.env$bins<-matrix(mrbin.env$bins[-which(rownames(mrbin.env$bins)%in%listTMP),],nrow=1)
              rownames(mrbin.env$bins)<-rownamesTMP
              colnames(mrbin.env$bins)<-colnamesTMP
            } else {
               mrbin.env$bins<-mrbin.env$bins[-which(rownames(mrbin.env$bins)%in%listTMP),]
            }
           }
         }
         mrbin.env$mrbinparam$NMRfolders<-mrbin.env$mrbinparam$NMRfolders[-which(mrbin.env$mrbinparam$NMRfolders%in%listTMP)]
      }
    }
 }
}

#' A function for setting group members.
#'
#' This function lets the user pick samples from a list to assign them to
#' groups. This function is meant only for use within the mrbin function.
#' @return {None}
#' @export
#' @examples
#' \donttest{ setFactors() }

setFactors<-function(){
   if(length(mrbin.env$mrbinparam$NMRfolders)>0){
      Factors<-rep("Group 0",length(mrbin.env$mrbinparam$NMRfolders))
      names(Factors)<-mrbin.env$mrbinparam$NMRfolders
      flag<-TRUE
      i<-0
      while(flag){
        i<-i+1
        listTMP<-utils::select.list(names(Factors),preselect = NULL, multiple = TRUE,title ="Please select group members",graphics=TRUE)
        groupNameTMP<-utils::select.list(c(paste("Group",i),"Enter new name"),preselect=paste("Group",i),
                  multiple=FALSE,title ="Group name?",graphics=TRUE)
        if(groupNameTMP=="Enter new name"){
          groupNameTMP<-readline(prompt=paste("New group name, press enter to use \"Group ",i,"\": ",sep=""))
          if(groupNameTMP=="") groupNameTMP<-paste("Group ",i,sep="")
        }
        Factors[listTMP]<-groupNameTMP
        select<-utils::select.list(c("Yes","No"),preselect = "No",multiple = FALSE,
            title = paste("Define additional groups?",sep=""),graphics=TRUE)
        if(select=="No") flag<-FALSE
      }
      Factors<-as.factor(Factors)
      mrbin.env$mrbinparam$Factors <- Factors
   }
}

#' A function for selecting NMR data folders.
#'
#' This function selects the correct folder selection function for the vendor
#' (currently only Bruker). This function is meant only for use within the
#' mrbin function.
#' @return An invisible list of folder names, or "Go back" or "stop"
#' @export
#' @examples
#' \donttest{ selectFolders() }

selectFolders<-function(){#Select NMR spectral folders
      selectionFolders<-""
      if(mrbin.env$mrbinparam$NMRvendor=="Bruker"){
          selectionFolders<-selectBrukerFolders()
      }  else {
          stop(paste("No folder selection function defined for vendor ",mrbin.env$mrbinparam$NMRvendor,".\n",sep=""))
      }
      invisible(selectionFolders)
}

#' A function for selecting Bruker NMR data folders.
#'
#' This function lets the user set NMR data folders interactively (for Bruker data). This function
#' is meant only for use within the mrbin function.
#' @return An invisible list of folder names, or "Go back" or "stop"
#' @export
#' @examples
#' \donttest{ selectBrukerFolders() }

selectBrukerFolders<-function(){#Select Bruker NMR spectral folders
  selectionFolders<-""
  NMRfoldersTMP<-NULL
  datanameDict<-c("1r","2rr")
  names(datanameDict)<-c("1D","2D")
  datanameTmp<-datanameDict[mrbin.env$mrbinparam$dimension]
  singleFolderFlag<-FALSE
  enterFolders<-utils::select.list(c("Browse...",#"Enter parent folder path manually",
                 "Go back"),
                 preselect="Browse...",title="Set NMR parent folder:",graphics=TRUE)
  if(enterFolders==""|length(enterFolders)==0){
       selectionFolders<-"stop"
  } else {
    if(enterFolders=="Go back"){
       selectionFolders<-"Go back"
    } else {
      folderPrompt<-"Enter folder path: "
      if(enterFolders=="Browse..."){
        folderPrompt<-"Enter starting folder path. (Examples: Windows: \"C:\\\", Apple: \"/\") : "
      }
      enterFoldersTMP<-readline(prompt=folderPrompt)
      if(!enterFoldersTMP==""){
       parentFolder<-gsub('\\\\',"/",enterFoldersTMP)
       #if(enterFolders=="Enter parent folder path manually"){
       #     selectFolders<-parentFolder
       #     selectList<-parentFolder
       #     singleFolderFlag<-TRUE
       #}
       if(enterFolders=="Browse..."){
       if(!singleFolderFlag) message("Hint: When reaching the NMR parent folder, click OK WITHOUT selecting\nany folder.")
       selectFlag<-0
       while(selectFlag<1){
         folderListFull<-list.dirs(path=parentFolder,recursive = FALSE,full.names=TRUE)
         folderList<-list.dirs(path=parentFolder,recursive = FALSE,full.names=FALSE)
         if(length(strsplit(parentFolder,split="/")[[1]])>1){
            folderListTMP<-c("..                                                                            ",folderList)
         } else {
            folderListTMP<-folderList
         }
         if(!singleFolderFlag){
             selectFolders<-utils::select.list(folderListTMP,preselect=NULL,multiple=TRUE,
                              title = "Go to parent folder, then click OK",graphics=TRUE)
         }
         if(!length(selectFolders)==1) {
            selectFolders<-parentFolder
            selectList<-parentFolder
            singleFolderFlag<-TRUE
         }
         if(length(selectFolders)==1&!singleFolderFlag){
              if(selectFolders=="..                                                                            "){
                   if(length(strsplit(parentFolder,split="/")[[1]])>1){
                      selectList<-paste(rev(rev(strsplit(parentFolder,split="/")[[1]])[-1]),sep="",collapse="/")
                   } else {
                      selectList<-parentFolder
                   }
              } else {
                   selectList<-folderListFull[which(folderList%in%selectFolders)]
              }
          }
          #}
          if(singleFolderFlag){
             if(mrbin.env$mrbinparam$verbose) message("Looking for NMR data in folder...")
             utils::flush.console()
             singleFolderFlag<-FALSE
             spectrum_path_list<-NULL
             subdirsTmp0<-list.dirs(path = selectList,recursive = FALSE,full.names=FALSE)
             spectrum_proc_path<-NULL
             if(length(subdirsTmp0)>0){
               #Look for Bruker NMR data in folder.
               #if(sum(sapply(suppressWarnings(as.numeric(subdirsTmp0)),is.numeric))>0){ #Look for folders "1", "2" etc (EXPNO)
               #    spectrum_path2<-paste(selectList,"/",subdirsTmp0[which(sapply(suppressWarnings(as.numeric(subdirsTmp0)),is.numeric))],sep="")
               if(sum(!is.na(suppressWarnings(as.numeric(subdirsTmp0))))>0){ #Look for folders "1", "2" etc (EXPNO)
                   spectrum_path2<-paste(selectList,"/",subdirsTmp0[!is.na(suppressWarnings(as.numeric(subdirsTmp0)))],sep="")
                   for(i in spectrum_path2){
                      if(dir.exists(paste(i,"/pdata",sep=""))){
                          subdirsTmp2<-list.dirs(path = paste(i,"/pdata",sep=""),recursive = FALSE,full.names=FALSE)
                          if(sum(!is.na(suppressWarnings(as.numeric(subdirsTmp2))))>0){ #Look for folders "1", "2" etc (PROCNO)
                                spectrum_path3<-paste(i,"/pdata/",subdirsTmp2[!is.na(suppressWarnings(as.numeric(subdirsTmp2)))],sep="")
                                spectrum_path_list<-c(spectrum_path_list,spectrum_path3)
                          }
                      }
                   }
                 if(length(spectrum_path_list)>0){
                   for(i in 1:length(spectrum_path_list)){
                       if(datanameTmp%in%list.files(spectrum_path_list[i])){
                           spectrum_proc_path<-c(spectrum_proc_path,spectrum_path_list[i])
                       }
                    }
                  spectrum_path_list<-NULL
                 }
               }
               for(isubDirs in 1:length(subdirsTmp0)){#Look for Bruker NMR data in subfolders.
                 subdirsTmp<-list.dirs(path = paste(selectList,"/",subdirsTmp0[isubDirs],sep=""),recursive = FALSE,full.names=FALSE)
                 if(length(subdirsTmp)>0){
                   if(sum(!is.na(suppressWarnings(as.numeric(subdirsTmp))))>0){ #Look for folders "1", "2" etc (EXPNO)
                       spectrum_path2<-paste(selectList,"/",subdirsTmp0[isubDirs],"/",subdirsTmp[!is.na(suppressWarnings(as.numeric(subdirsTmp)))],sep="")
                       for(i in spectrum_path2){
                          if(dir.exists(paste(i,"/pdata",sep=""))){
                              subdirsTmp2<-list.dirs(path = paste(i,"/pdata",sep=""),recursive = FALSE,full.names=FALSE)
                              if(sum(!is.na(suppressWarnings(as.numeric(subdirsTmp2))))>0){ #Look for folders "1", "2" etc (PROCNO)
                                    spectrum_path3<-paste(i,"/pdata/",subdirsTmp2[!is.na(suppressWarnings(as.numeric(subdirsTmp2)))],sep="")
                                    spectrum_path_list<-c(spectrum_path_list,spectrum_path3)
                              }
                          }
                       }
                   }
                 }
               }
               if(length(spectrum_path_list)>0){
                 for(i in 1:length(spectrum_path_list)){
                     if(datanameTmp%in%list.files(spectrum_path_list[i])){
                         spectrum_proc_path<-c(spectrum_proc_path,spectrum_path_list[i])
                     }
                  }
               }
             }
             if(is.null(spectrum_proc_path)){
                message(paste("No Bruker NMR data found in folder."))
                noDataFoundExit<-utils::select.list(c("Continue browsing","Exit"),
                                preselect="Continue browsing",
                                title = "No NMR data found in folder.",graphics=TRUE)
                if(noDataFoundExit=="Continue browsing"){
                      selectFlag<-0
                } else {
                      selectFlag<-1
                }
             } else {
               spectrum_proc_pathTMPName<-NULL
               spectrum_proc_pathTMPEXPNO<-NULL
               for(i_findDuplicates in 1:length(spectrum_proc_path)){
                 spectrum_proc_pathTMPName<-c(spectrum_proc_pathTMPName,rev(strsplit(spectrum_proc_path,split="/")[[i_findDuplicates]])[4])
                 spectrum_proc_pathTMPEXPNO<-c(spectrum_proc_pathTMPEXPNO,rev(strsplit(spectrum_proc_path,split="/")[[i_findDuplicates]])[3])
               }
               if(sum(duplicated(spectrum_proc_pathTMPName))>0){
                 for(i_findDuplicates in unique(spectrum_proc_pathTMPName[which(duplicated(spectrum_proc_pathTMPName))])){
                   spectrum_proc_path[which(spectrum_proc_pathTMPName == i_findDuplicates)]<-
                       spectrum_proc_path[which(spectrum_proc_pathTMPName == i_findDuplicates)][
                       order(as.numeric(spectrum_proc_pathTMPEXPNO[which(spectrum_proc_pathTMPName == i_findDuplicates)]))]
                 }
               }
               spectrum_proc_pathTMP<-c("Select all",spectrum_proc_path,"Select all")
               NMRfoldersTMP<-utils::select.list(spectrum_proc_pathTMP,preselect=NULL,multiple=TRUE,
                         title = "Select data sets",graphics=TRUE)
               if("Select all"%in%NMRfoldersTMP) NMRfoldersTMP<-spectrum_proc_path
               #NMRfolders<-c(NMRfolders,NMRfoldersTMP)
               if(!is.null(NMRfoldersTMP)) mrbin.env$mrbinparam$NMRfolders<-unique(c(mrbin.env$mrbinparam$NMRfolders,NMRfoldersTMP))
               if(mrbin.env$mrbinparam$verbose) cat(paste("Adding to list:\n",paste(NMRfoldersTMP,"\n",
                                                    sep="",collapse="")))
               addSpectrumTMP<-TRUE
               while(addSpectrumTMP){
                 yesorno<-utils::select.list(c("Keep current spectra list","Add additional spectra","Remove spectra from list"),
                           preselect="Keep current spectra list",multiple=FALSE,title="Add additional spectra?",graphics=TRUE)
                 if(yesorno==""){
                     addSpectrumTMP<-FALSE
                 }
                 if(yesorno=="Add additional spectra"){
                     selectFlag<-0
                     addSpectrumTMP<-FALSE
                 }
                 if(yesorno=="Remove spectra from list") {
                     removeSpectrum()
                     addSpectrumTMP<-TRUE
                 }
                 if(yesorno=="Keep current spectra list") {
                     selectFlag<-1
                     addSpectrumTMP<-FALSE
                 }
               }
             }
         } else {
            parentFolder<-selectList
         }
       }
       }
      }
      #if(!is.null(NMRfolders)) mrbin.env$mrbinparam$NMRfolders<-unique(c(mrbin.env$mrbinparam$NMRfolders,NMRfolders))
    }
  }
  invisible(selectionFolders)
}

#' A function for binning multiple NMR spectra.
#'
#' This function creates bins for each spectrum in mrbin.env$mrbinparam$NMRfolders and
#' saves the bins to mrbin.env$bins. Should only be run within the mrbin function.
#' @return {None}
#' @export
#' @examples
#' \donttest{ binMultiNMR() }

binMultiNMR<-function(){
 if(!is.null(mrbin.env$mrbinparam$NMRfolders)){
    mrbin.env$bins<-NULL
    mrbin.env$mrbinTMP$binNames<-NULL
    mrbin.env$mrbinTMP$binTMP<-NULL
    #Open and bin all spectra
    mrbin.env$mrbinTMP$binsRaw<-NULL
    mrbin.env$mrbinTMP$noise_level_Raw<-rep(NA,length(mrbin.env$mrbinparam$NMRfolders))
    mrbin.env$mrbinTMP$meanNumberOfPointsPerBin<-rep(NA,length(mrbin.env$mrbinparam$NMRfolders))
    if(mrbin.env$mrbinparam$verbose){
      cat("Binning spectrum: ")
      utils::flush.console()
    }
    for(i in 1:length(mrbin.env$mrbinparam$NMRfolders)){
        if(mrbin.env$mrbinparam$verbose) cat(paste(i," ",sep=""))
        utils::flush.console()
        mrbin.env$mrbinTMP$currentFolder<-mrbin.env$mrbinparam$NMRfolders[i]
        readNMR()
        calculateNoise()
        mrbin.env$mrbinTMP$noise_level_Raw[i]<-mrbin.env$mrbinTMP$noise_level_Raw_TMP
        mrbin.env$mrbinTMP$binTMP<-NULL
        binSingleNMR()
        mrbin.env$mrbinTMP$meanNumberOfPointsPerBin[i]<-mrbin.env$mrbinTMP$meanNumberOfPointsPerBin_TMP
        if(is.null(mrbin.env$mrbinTMP$binsRaw)){
            mrbin.env$mrbinTMP$binsRaw<-matrix(rep(0,length(mrbin.env$mrbinTMP$binTMP)*
                                          length(mrbin.env$mrbinparam$NMRfolders)),
                                          nrow=length(mrbin.env$mrbinparam$NMRfolders))
            colnames(mrbin.env$mrbinTMP$binsRaw)<-names(mrbin.env$mrbinTMP$binTMP)
            rownames(mrbin.env$mrbinTMP$binsRaw)<-paste("TemporaryRowName_",1:length(mrbin.env$mrbinparam$NMRfolders),sep="")
        }
        if(is.null(mrbin.env$mrbinTMP$binTMP)){
            stop("Error while binning spectrum.")
        } else {
            mrbin.env$mrbinTMP$binsRaw[i,]<-mrbin.env$mrbinTMP$binTMP
            currentSpectrumNameTMP<-mrbin.env$mrbinTMP$currentSpectrumName
            i_currentSpectrumNameTMP<-2
            while(currentSpectrumNameTMP%in%rownames(mrbin.env$mrbinTMP$binsRaw)){
                              currentSpectrumNameTMP<-paste(mrbin.env$mrbinTMP$currentSpectrumName,".",i_currentSpectrumNameTMP,sep="")
               i_currentSpectrumNameTMP<-i_currentSpectrumNameTMP+1
            }
            if(i_currentSpectrumNameTMP>2) warning(paste("Renamed duplicate spectrum title (",
                                           currentSpectrumNameTMP,"). Consider using an alternative naming method.\n",
                                           sep=""))
            rownames(mrbin.env$mrbinTMP$binsRaw)[i]<-currentSpectrumNameTMP
            mrbin.env$mrbinTMP$currentSpectrumName<-currentSpectrumNameTMP
        }
    }
    cat("Done.\n")
    utils::flush.console()
    mrbin.env$bins<-mrbin.env$mrbinTMP$binsRaw
    #mrbin.env$mrbinparam$numberOfFeaturesRaw<-ncol(mrbin.env$bins)
    mrbin.env$mrbinparam$noise_level<-mrbin.env$mrbinTMP$noise_level_Raw*(mrbin.env$mrbinTMP$meanNumberOfPointsPerBin^(-.5))
 }
}

#' A function for creating bin regions.
#'
#' This function creates regions for the bins. This function is
#' meant only for use within the mrbin function.
#' @return {None}
#' @export
#' @examples
#' \donttest{ createBinRegions() }

createBinRegions<-function(){
   createBinNumbers()
   if(mrbin.env$mrbinparam$binMethod=="Rectangular bins"){

       mrbin.env$mrbinTMP$binRegions<-matrix(ncol=4,
                                        nrow=mrbin.env$mrbinTMP$nbins,
                                        dimnames=list(NULL,c("left","right","top","bottom")))
       if(mrbin.env$mrbinparam$dimension=="1D"){
              mrbin.env$mrbinTMP$binRegions[,1]<-mrbin.env$mrbinparam$binRegion[1]-((1:mrbin.env$mrbinTMP$nbins)-1)*mrbin.env$mrbinparam$binwidth1D
              mrbin.env$mrbinTMP$binRegions[,2]<-mrbin.env$mrbinparam$binRegion[1]-((1:mrbin.env$mrbinTMP$nbins))*mrbin.env$mrbinparam$binwidth1D
              #rownames(mrbin.env$mrbinTMP$binRegions)<-apply(mrbin.env$mrbinTMP$binRegions[,1:2],1,mean)
              #rownames(mrbin.env$mrbinTMP$binRegions)<-apply(mrbin.env$mrbinTMP$binRegions[,1:2],1,paste,collapse=",")
       }
       if(mrbin.env$mrbinparam$dimension=="2D"){
              mrbin.env$mrbinTMP$binRegions[,1]<-rep(mrbin.env$mrbinparam$binRegion[1]-((1:mrbin.env$mrbinTMP$nbins2)-1)*mrbin.env$mrbinparam$binwidth2D,
                                                 mrbin.env$mrbinTMP$nbins1)
              mrbin.env$mrbinTMP$binRegions[,2]<-rep(mrbin.env$mrbinparam$binRegion[1]-((1:mrbin.env$mrbinTMP$nbins2))*mrbin.env$mrbinparam$binwidth2D,
                                                 mrbin.env$mrbinTMP$nbins1)
              mrbin.env$mrbinTMP$binRegions[,3]<-sort(rep(mrbin.env$mrbinparam$binRegion[4]-((1:mrbin.env$mrbinTMP$nbins1))*mrbin.env$mrbinparam$binheight,
                                                 mrbin.env$mrbinTMP$nbins2),decreasing=TRUE)
              mrbin.env$mrbinTMP$binRegions[,4]<-sort(rep(mrbin.env$mrbinparam$binRegion[4]-((1:mrbin.env$mrbinTMP$nbins1)-1)*mrbin.env$mrbinparam$binheight,
                                                 mrbin.env$mrbinTMP$nbins2),decreasing=TRUE)
              #rownames(mrbin.env$mrbinTMP$binRegions)<-paste(apply(mrbin.env$mrbinTMP$binRegions[,3:4],1,mean),
              #                                          apply(mrbin.env$mrbinTMP$binRegions[,1:2],1,mean),sep=",")
              #rownames(mrbin.env$mrbinTMP$binRegions)<-apply(mrbin.env$mrbinTMP$binRegions,1,paste,collapse=",")

       }
       #mrbin.env$mrbinTMP$binNames<-rownames(mrbin.env$mrbinTMP$binRegions)
  }
  mrbin.env$mrbinparam$numberOfFeaturesRaw<-nrow(mrbin.env$mrbinTMP$binRegions)
}

#' A function for creating bin titles.
#'
#' This function creates titles for the bins to represent their ppm range. This function is
#' meant only for use within the mrbin function.
#' @return {None}
#' @export
#' @examples
#' \donttest{ createBinNames() }

createBinNames<-function(){
   if(mrbin.env$mrbinparam$dimension=="1D"){
          rownames(mrbin.env$mrbinTMP$binRegions)<-apply(mrbin.env$mrbinTMP$binRegions[,1:2],1,paste,collapse=",")
   }
   if(mrbin.env$mrbinparam$dimension=="2D"){
          rownames(mrbin.env$mrbinTMP$binRegions)<-apply(mrbin.env$mrbinTMP$binRegions,1,paste,collapse=",")
   }
   mrbin.env$mrbinTMP$binNames<-rownames(mrbin.env$mrbinTMP$binRegions)
   colnames(mrbin.env$bins)<-rownames(mrbin.env$mrbinTMP$binRegions)
}

#' A function for creating bin numbers.
#'
#' This function calculates numbers of bins from the chosen parameters. This function is
#' meant only for use within the mrbin function.
#' @return {None}
#' @export
#' @examples
#' \donttest{ createBinNumbers() }

createBinNumbers<-function(){
   if(mrbin.env$mrbinparam$binMethod=="Custom bin list"){
          mrbin.env$mrbinTMP$nbins<-nrow(mrbin.env$mrbinTMP$binRegions)
   }
   if(mrbin.env$mrbinparam$binMethod=="Rectangular bins"){
       if(mrbin.env$mrbinparam$dimension=="2D"){
          mrbin.env$mrbinTMP$nbins2<-ceiling((mrbin.env$mrbinparam$binRegion[1]-mrbin.env$mrbinparam$binRegion[2])/mrbin.env$mrbinparam$binwidth2D)
          mrbin.env$mrbinTMP$nbins1<-ceiling((mrbin.env$mrbinparam$binRegion[4]-mrbin.env$mrbinparam$binRegion[3])/mrbin.env$mrbinparam$binheight)
          mrbin.env$mrbinTMP$nbins<-mrbin.env$mrbinTMP$nbins2*mrbin.env$mrbinTMP$nbins1
      }
      if(mrbin.env$mrbinparam$dimension=="1D"){
          mrbin.env$mrbinTMP$nbins<-ceiling((mrbin.env$mrbinparam$binRegion[1]-mrbin.env$mrbinparam$binRegion[2])/mrbin.env$mrbinparam$binwidth1D)
      }
   }
}

#' A function for binning a single NMR spectrum.
#'
#' This function creates bins for the current spectrum. This function is
#' meant only for use within the mrbin function.
#' @return {None}
#' @export
#' @examples
#' \donttest{ binSingleNMR() }

binSingleNMR<-function(){
 if(!is.null(mrbin.env$mrbinTMP$currentFolder)){
   warningFlag<-FALSE
    numberOfPointsPerBin<-NULL
   if(mrbin.env$mrbinparam$dimension=="2D"){#2d spectra
      #if(is.null(mrbin.env$mrbinTMP$binNames))          createBinRegions()
      #Create index of signals in each bin
      NMRspectrumRownames<-as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrum))
      NMRspectrumColnames<-as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrum))
      #counter<-1
      mrbin.env$mrbinTMP$binTMP<-rep(0,nrow(mrbin.env$mrbinTMP$binRegions))
      names(mrbin.env$mrbinTMP$binTMP)<-rownames(mrbin.env$mrbinTMP$binRegions)
      for(ibinTMP in 1:nrow(mrbin.env$mrbinTMP$binRegions)){
        numberOfPointsPerBinTMP<-sum(NMRspectrumRownames<=mrbin.env$mrbinTMP$binRegions[ibinTMP,4]&
                                         NMRspectrumRownames>mrbin.env$mrbinTMP$binRegions[ibinTMP,3])*
                                         sum(NMRspectrumColnames<=mrbin.env$mrbinTMP$binRegions[ibinTMP,1]&
                                         NMRspectrumColnames>mrbin.env$mrbinTMP$binRegions[ibinTMP,2])
        if(numberOfPointsPerBinTMP>0){
          numberOfPointsPerBin<-c(numberOfPointsPerBin,numberOfPointsPerBinTMP)
          mrbin.env$mrbinTMP$binTMP[ibinTMP]<-sum(mrbin.env$mrbinTMP$currentSpectrum[
                                         NMRspectrumRownames<=mrbin.env$mrbinTMP$binRegions[ibinTMP,4]&
                                         NMRspectrumRownames>mrbin.env$mrbinTMP$binRegions[ibinTMP,3],
                                         NMRspectrumColnames<=mrbin.env$mrbinTMP$binRegions[ibinTMP,1]&
                                         NMRspectrumColnames>mrbin.env$mrbinTMP$binRegions[ibinTMP,2]])/
                                        (sum(NMRspectrumRownames<=mrbin.env$mrbinTMP$binRegions[ibinTMP,4]&
                                         NMRspectrumRownames>mrbin.env$mrbinTMP$binRegions[ibinTMP,3])*
                                         sum(NMRspectrumColnames<=mrbin.env$mrbinTMP$binRegions[ibinTMP,1]&
                                         NMRspectrumColnames>mrbin.env$mrbinTMP$binRegions[ibinTMP,2]))
        } else {
          warningFlag<-TRUE
        }
     }
   } else {#1D spectra
       #if(is.null(mrbin.env$mrbinTMP$binNames)) createBinRegions()
       mrbin.env$mrbinTMP$binTMP<-rep(0,nrow(mrbin.env$mrbinTMP$binRegions))
       names(mrbin.env$mrbinTMP$binTMP)<-rownames(mrbin.env$mrbinTMP$binRegions)
       #counter<-1
       NMRspectrumNames<-as.numeric(names(mrbin.env$mrbinTMP$currentSpectrum))
     for(ibinTMP in 1:nrow(mrbin.env$mrbinTMP$binRegions)){
        numberOfPointsPerBinTMP<-sum(NMRspectrumNames<=mrbin.env$mrbinTMP$binRegions[ibinTMP,1]&
          NMRspectrumNames>mrbin.env$mrbinTMP$binRegions[ibinTMP,2])
        if(numberOfPointsPerBinTMP>0){
           numberOfPointsPerBin<-c(numberOfPointsPerBin,numberOfPointsPerBinTMP)
           mrbin.env$mrbinTMP$binTMP[ibinTMP]<-sum(mrbin.env$mrbinTMP$currentSpectrum[
                                   NMRspectrumNames<=mrbin.env$mrbinTMP$binRegions[ibinTMP,1]&
                                   NMRspectrumNames>mrbin.env$mrbinTMP$binRegions[ibinTMP,2]])/
                                  (sum(NMRspectrumNames<=mrbin.env$mrbinTMP$binRegions[ibinTMP,1]&
                                   NMRspectrumNames>mrbin.env$mrbinTMP$binRegions[ibinTMP,2]))
        } else {
          warningFlag<-TRUE
        }
     }
    }
    mrbin.env$mrbinTMP$meanNumberOfPointsPerBin_TMP<-stats::median(numberOfPointsPerBin)
    if(warningFlag){
      warning(paste("Binning region may be larger than total spectrum size for spectrum ",mrbin.env$mrbinTMP$currentFolder,sep=""))
    }
 }
}

#' A function for reading NMR spectra.
#'
#' This function picks the correct NMR reading function, based on vendor. This function is
#' meant only for use within the mrbin function.
#' @return {None}
#' @export
#' @examples
#' \donttest{ readNMR() }

readNMR<-function(){#Read NMR spectral data
 if(!is.null(mrbin.env$mrbinTMP$currentFolder)){
  if(mrbin.env$mrbinparam$NMRvendor=="Bruker"){
      readBruker()
  }  else {
      stop(paste("No data import function defined for vendor ",mrbin.env$mrbinparam$NMRvendor,".\n",sep=""))
  }
  if(mrbin.env$mrbinparam$referenceScaling=="Yes") referenceScaling()
 }
}

#' A function for reading Bruker NMR spectra.
#'
#' This function reads Bruker NMR data.
#' @param folder Defines the exact NMR data folder. If NULL, mrbin parameter set is used
#' @param dimension Defines the data dimension, "1D" or "2D". Only used if "folder" is not NULL
#' @return An (invisible) object containing spectral data
#' @export
#' @examples
#' exampleData<-readBruker(folder=system.file("extdata/1/10/pdata/10",package="mrbin"),
#'                         dimension="1D")

readBruker<-function(folder=NULL,dimension=NULL){#Read Bruker NMR spectral data
 datanameDict<-c("1r","2rr")
 names(datanameDict)<-c("1D","2D")
 if(is.null(folder)){
   spectrum_proc_path<-gsub('\\\\',"/",mrbin.env$mrbinTMP$currentFolder)
   datanameTmp<-datanameDict[mrbin.env$mrbinparam$dimension]
 } else {
   spectrum_proc_path<-folder
   datanameTmp<-datanameDict[dimension]
 }
 if(!is.null(spectrum_proc_path)){
     BYTORDP_Dict<-c("little","big")
     names(BYTORDP_Dict)<-c(0,1)
     TITLE<-scan(file=paste(spectrum_proc_path,"/title",sep=""),what="character",sep="\n",quiet=TRUE)[1]
     proc<-scan(file=paste(spectrum_proc_path,"/procs",sep=""),what="character",sep="\n",quiet=TRUE)
     proc<-gsub("#","",proc)
     proc<-gsub("\\$"," ",proc)
     SI2<-as.numeric(strsplit(proc[grep(" SI=",proc)],split="= ")[[1]][2])
     BYTORDP<-BYTORDP_Dict[strsplit(proc[grep(" BYTORDP=",proc)],split="= ")[[1]][2]]
     NC_proc<-as.numeric(strsplit(proc[grep(" NC_proc=",proc)],split="= ")[[1]][2])
     XDIM2<-as.numeric(strsplit(proc[grep(" XDIM=",proc)],split="= ")[[1]][2])
     OFFSET2<-as.numeric(strsplit(proc[grep(" OFFSET=",proc)],split="= ")[[1]][2])#ppm
     SF2<-as.numeric(strsplit(proc[grep(" SF=",proc)],split="= ")[[1]][2])#MHz
     SW_p2<-as.numeric(strsplit(proc[grep(" SW_p=",proc)],split="= ")[[1]][2])#Hz
     rightlimit2<-OFFSET2-SW_p2/SF2
     leftlimit2<-OFFSET2
     frequencynames2<-OFFSET2-(0:(SI2-1))*SW_p2/SF2/SI2
     n<-SI2
     if(datanameTmp=="2rr"){#Only 2D
         proc2<-scan(file=paste(spectrum_proc_path,"/proc2s",sep=""),what="character",sep="\n",quiet=TRUE)
         proc2<-gsub("#","",proc2)
         proc2<-gsub("\\$"," ",proc2)
         SI1<-as.numeric(strsplit(proc2[grep(" SI=",proc2)],split="= ")[[1]][2])
         XDIM1<-as.numeric(strsplit(proc2[grep(" XDIM=",proc2)],split="= ")[[1]][2])
         OFFSET1<-as.numeric(strsplit(proc2[grep(" OFFSET=",proc2)],split="= ")[[1]][2])#ppm
         SF1<-as.numeric(strsplit(proc2[grep(" SF=",proc2)],split="= ")[[1]][2])#MHz
         SW_p1<-as.numeric(strsplit(proc2[grep(" SW_p=",proc2)],split="= ")[[1]][2])#Hz
         n<-SI2*SI1
         rightlimit1<-OFFSET1-SW_p1/SF1
         leftlimit1<-OFFSET1
         frequencynames1<-OFFSET1-(0:(SI1-1))*SW_p1/SF1/SI1
     }
     currentSpectrum<-readBin(paste(spectrum_proc_path,"/",datanameTmp,sep=""),
                          what='integer',size=4,n=n,endian=BYTORDP)
     currentSpectrum<-currentSpectrum*2^NC_proc
     if(datanameTmp=="2rr"){
        currentSpectrumTMP<-currentSpectrum
        currentSpectrum<-matrix(currentSpectrum,ncol=SI2,byrow =TRUE)#ncol=SI2)
        #Needs to be sorted,XDIM2: number of columns/block, XDIM1: rows/block
        counter<-0
        for(j in 1:(nrow(currentSpectrum)/XDIM1)){
            for(i in 1:(ncol(currentSpectrum)/XDIM2)){
                currentSpectrum[(j-1)*XDIM1+(1:XDIM1),(i-1)*XDIM2+(1:XDIM2)]<-
                   matrix(currentSpectrumTMP[counter*XDIM2*XDIM1 +(1:(XDIM2*XDIM1))],ncol=XDIM2,byrow=TRUE)
                counter<-counter+1
            }
        }
        rownames(currentSpectrum)<-frequencynames1
        colnames(currentSpectrum)<-frequencynames2
     } else {
       names(currentSpectrum)<-frequencynames2
     }
  if(is.null(folder)){
    mrbin.env$mrbinTMP$currentSpectrum<-currentSpectrum
    mrbin.env$mrbinTMP$currentSpectrumTitle<-TITLE
    mrbin.env$mrbinTMP$currentSpectrumFolderName<-rev(strsplit(mrbin.env$mrbinTMP$currentFolder,"/")[[1]])[4]
    mrbin.env$mrbinTMP$currentSpectrumEXPNO<-rev(strsplit(mrbin.env$mrbinTMP$currentFolder,"/")[[1]])[3]
    mrbin.env$mrbinTMP$currentSpectrumFolderName_EXPNO<-paste(
                 mrbin.env$mrbinTMP$currentSpectrumFolderName,
                 paste(c("_","0","0","0","0")[1:max(1,5-nchar(mrbin.env$mrbinTMP$currentSpectrumEXPNO))],sep="",collapse=""),
                 mrbin.env$mrbinTMP$currentSpectrumEXPNO,sep="")
    if(mrbin.env$mrbinparam$useAsNames=="Spectrum titles")    mrbin.env$mrbinTMP$currentSpectrumName<-TITLE
    if(mrbin.env$mrbinparam$useAsNames=="Folder names")    mrbin.env$mrbinTMP$currentSpectrumName<-mrbin.env$mrbinTMP$currentSpectrumFolderName
    if(mrbin.env$mrbinparam$useAsNames=="Folder names and EXPNO")    mrbin.env$mrbinTMP$currentSpectrumName<-mrbin.env$mrbinTMP$currentSpectrumFolderName_EXPNO
  }
  invisible(currentSpectrum)
 }
}

#' A function for scaling to the reference area.
#'
#' This function scales NMR data to the reference area.
#' @return {None}
#' @export
#' @examples
#' resetEnv()
#' mrbin(silent=TRUE,setDefault=TRUE,parameters=list(dimension="1D",
#'          referenceScaling="No",binwidth1D=0.05,PQNScaling="No",PCA="No",
#'          NMRfolders=system.file("extdata/1/10/pdata/10",package="mrbin")))
#' referenceScaling()

referenceScaling<-function(){
 if(!is.null(mrbin.env$bins)){
    if(mrbin.env$mrbinparam$dimension=="1D"){
         mrbin.env$mrbinTMP$currentSpectrum<-mrbin.env$mrbinTMP$currentSpectrum/mean(mrbin.env$mrbinTMP$currentSpectrum[
                                         which(as.numeric(names(mrbin.env$mrbinTMP$currentSpectrum))<=
                                         mrbin.env$mrbinparam$reference1D[1]&
                                         as.numeric(names(mrbin.env$mrbinTMP$currentSpectrum))>
                                         mrbin.env$mrbinparam$reference1D[2])],na.rm=TRUE)
    }
    if(mrbin.env$mrbinparam$dimension=="2D"){
         selectedTMP1<-which(as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrum))<=
                         mrbin.env$mrbinparam$reference2D[4]&as.numeric(rownames(
                         mrbin.env$mrbinTMP$currentSpectrum))>mrbin.env$mrbinparam$reference2D[3])
         selectedTMP2<-which(as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrum))<=
                         mrbin.env$mrbinparam$reference2D[1]&
                         as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrum))>
                         mrbin.env$mrbinparam$reference2D[2])

         mrbin.env$mrbinTMP$currentSpectrum<-mrbin.env$mrbinTMP$currentSpectrum/mean(mrbin.env$mrbinTMP$currentSpectrum[
                                     selectedTMP1,selectedTMP2],na.rm=TRUE)
    }
 }
}

#' A function for summing bins.
#'
#' This function sums up bins. The sums are saved to the middle (median) bin of
#' the original area. All other bins of the area are removed then. This is handy
#' for signals that are know to vary between spectra due to pH or salt content,
#' such as citric acid. Intended for use within the mrbin function.
#' @return {None}
#' @export
#' @examples
#' \donttest{ sumBins() }

sumBins<-function(){#sum up regions with shifting peaks and remove remaining bins
 #if(!is.null(mrbin.env$bins)){
  if(nrow(mrbin.env$mrbinparam$sumBinList)>0){
    for(i in 1:nrow(mrbin.env$mrbinparam$sumBinList)){
       limits<-mrbin.env$mrbinparam$sumBinList[i,]
       if(mrbin.env$mrbinparam$dimension=="1D"){
          TMP<-apply(mrbin.env$mrbinTMP$binRegions[,1:2],1,mean)>limits[2]& apply(mrbin.env$mrbinTMP$binRegions[,1:2],1,mean)<limits[1]
          if(sum(TMP)>0){
              i_TMP<-quantile(x=1:sum(TMP), probs = .5,type=3)#define "middle" bin. This one will be kept
              #mrbin.env$bins[,which(TMP)[i_TMP]]<-apply(mrbin.env$bins[,which(TMP)],1,sum)
              #if(nrow(mrbin.env$bins)==1){
              #  rownamesTMP<-rownames(mrbin.env$bins)
              #  colnamesTMP<-colnames(mrbin.env$bins)[-which(TMP)[-i_TMP]]
              #  mrbin.env$bins<-matrix(mrbin.env$bins[,-which(TMP)[-i_TMP]],nrow=1)
              #  rownames(mrbin.env$bins)<-rownamesTMP
              #  colnames(mrbin.env$bins)<-colnamesTMP
              #} else {
              #  mrbin.env$bins<-mrbin.env$bins[,-which(TMP)[-i_TMP]]
              #}
              mrbin.env$mrbinTMP$binRegions[which(TMP)[i_TMP],1]<-max(mrbin.env$mrbinTMP$binRegions[which(TMP),1])
              mrbin.env$mrbinTMP$binRegions[which(TMP)[i_TMP],2]<-min(mrbin.env$mrbinTMP$binRegions[which(TMP),2])
              mrbin.env$mrbinTMP$binRegions<-mrbin.env$mrbinTMP$binRegions[-which(TMP)[-i_TMP],]
          }
       } else {#2D limits=c(4.04,4.08,58,60)
           NMRdataNames<-cbind(apply(mrbin.env$mrbinTMP$binRegions[,3:4],1,mean),apply(mrbin.env$mrbinTMP$binRegions[,1:2],1,mean))
           TMP<-NMRdataNames[,2]>limits[2]&NMRdataNames[,2]<limits[1]&
                  NMRdataNames[,1]>limits[3]& NMRdataNames[,1]<limits[4]
           if(sum(TMP)>0){
              i_TMP<-quantile(x=1:sum(TMP), probs = .5,type=3)#define "middle" bin. This one will be kept
              #mrbin.env$bins[,which(TMP)[i_TMP]]<-apply(mrbin.env$bins[,which(TMP)],1,sum)
              #if(nrow(mrbin.env$bins)==1){
              #  rownamesTMP<-rownames(mrbin.env$bins)
              #  colnamesTMP<-colnames(mrbin.env$bins)[-which(TMP)[-i_TMP]]
              #  mrbin.env$bins<-matrix(mrbin.env$bins[,-which(TMP)[-i_TMP]],nrow=1)
              #  rownames(mrbin.env$bins)<-rownamesTMP
              #  colnames(mrbin.env$bins)<-colnamesTMP
              #} else {
              #  mrbin.env$bins<-mrbin.env$bins[,-which(TMP)[-i_TMP]]
              #}
              mrbin.env$mrbinTMP$binRegions[which(TMP)[i_TMP],1]<-max(mrbin.env$mrbinTMP$binRegions[which(TMP),1])
              mrbin.env$mrbinTMP$binRegions[which(TMP)[i_TMP],2]<-min(mrbin.env$mrbinTMP$binRegions[which(TMP),2])
              mrbin.env$mrbinTMP$binRegions[which(TMP)[i_TMP],3]<-min(mrbin.env$mrbinTMP$binRegions[which(TMP),3])
              mrbin.env$mrbinTMP$binRegions[which(TMP)[i_TMP],4]<-max(mrbin.env$mrbinTMP$binRegions[which(TMP),4])
              mrbin.env$mrbinTMP$binRegions<-mrbin.env$mrbinTMP$binRegions[-which(TMP)[-i_TMP],]
           }
       }
    }
  }
  mrbin.env$mrbinparam$numberOfFeaturesAfterSummingBins<-nrow(mrbin.env$mrbinTMP$binRegions)
 #}
}

#' A function for removing the solvent region.
#'
#' This function removes the solvent region. Should only be run within the mrbin function.
#' @return {None}
#' @export
#' @examples
#' \donttest{ removeSolvent() }

removeSolvent<-function(){
 #if(!is.null(mrbin.env$bins)){
   solventTMP<-(mrbin.env$mrbinparam$solventRegion[2]>=mrbin.env$mrbinTMP$binRegions[,2]&
               mrbin.env$mrbinparam$solventRegion[2]<=mrbin.env$mrbinTMP$binRegions[,1])|
               (mrbin.env$mrbinparam$solventRegion[1]>=mrbin.env$mrbinTMP$binRegions[,2]&
               mrbin.env$mrbinparam$solventRegion[1]<=mrbin.env$mrbinTMP$binRegions[,1])|
               (mrbin.env$mrbinTMP$binRegions[,2]>=mrbin.env$mrbinparam$solventRegion[2]&
               mrbin.env$mrbinTMP$binRegions[,2]<=mrbin.env$mrbinparam$solventRegion[1])|
               (mrbin.env$mrbinTMP$binRegions[,1]>=mrbin.env$mrbinparam$solventRegion[2]&
               mrbin.env$mrbinTMP$binRegions[,1]<=mrbin.env$mrbinparam$solventRegion[1])
   #if(mrbin.env$mrbinparam$dimension=="1D"){
       #NMRdataNames<-apply(mrbin.env$mrbinTMP$binRegions[,1:2],1,mean)
   #}
   #if(mrbin.env$mrbinparam$dimension=="2D"){
       #NMRdataNames<-cbind(apply(mrbin.env$mrbinTMP$binRegions[,3:4],1,mean),apply(mrbin.env$mrbinTMP$binRegions[,1:2],1,mean))
       #solventTMP<-NMRdataNames[,2]>mrbin.env$mrbinparam$solventRegion[2]&NMRdataNames[,2]<mrbin.env$mrbinparam$solventRegion[1]
   #}
   if(sum(solventTMP)>0){
      #if(nrow(mrbin.env$bins)==1){
      #  rownamesTMP<-rownames(mrbin.env$bins)
      #  colnamesTMP<-colnames(mrbin.env$bins)[-which(solventTMP)]
      #  mrbin.env$bins<-matrix(mrbin.env$bins[,-which(solventTMP)],nrow=1)
      #  rownames(mrbin.env$bins)<-rownamesTMP
      #  colnames(mrbin.env$bins)<-colnamesTMP
      #} else {
      #  mrbin.env$bins<-mrbin.env$bins[,-which(solventTMP)]
      #}
      mrbin.env$mrbinTMP$binRegions<-mrbin.env$mrbinTMP$binRegions[-which(solventTMP),]
   }
   mrbin.env$mrbinparam$numberOfFeaturesAfterRemovingSolvent<-nrow(mrbin.env$mrbinTMP$binRegions)
 #}
}

#' A function for removing additional regions.
#'
#' This function removes additional regions. This can be useful when some areas
#' are visibly affected by spectral artifacts. Should only be run from within the mrbin function.
#' @return {None}
#' @export
#' @examples
#' \donttest{ removeAreas() }

removeAreas<-function(){#limits=c(4.75,4.95,-10,160)
 #if(!is.null(mrbin.env$bins)){
  if(nrow(mrbin.env$mrbinparam$removeAreaList)>0){
     if(mrbin.env$mrbinparam$dimension=="1D")    NMRdataNames<-apply(mrbin.env$mrbinTMP$binRegions[,1:2],1,mean)
     if(mrbin.env$mrbinparam$dimension=="2D"){
         NMRdataNames<-cbind(apply(mrbin.env$mrbinTMP$binRegions[,3:4],1,mean),apply(mrbin.env$mrbinTMP$binRegions[,1:2],1,mean))
     }
     removeTMP<-NULL
     for(i in 1:nrow(mrbin.env$mrbinparam$removeAreaList)){
         limits<-mrbin.env$mrbinparam$removeAreaList[i,]
         if(mrbin.env$mrbinparam$dimension=="1D"){
             removeTMP2<-NMRdataNames>limits[2]&NMRdataNames<limits[1]
         }
         if(mrbin.env$mrbinparam$dimension=="2D"){
             removeTMP2<-NMRdataNames[,2]>limits[2]&NMRdataNames[,2]<limits[1]&
                NMRdataNames[,1]>limits[3]& NMRdataNames[,1]<limits[4]
         }
         if(sum(removeTMP2)>0)        removeTMP<-c(removeTMP,which(removeTMP2))
     }
     if(!is.null(removeTMP)){
         removeTMP<-unique(removeTMP)
         #if(nrow(mrbin.env$bins)==1){
         #   rownamesTMP<-rownames(mrbin.env$bins)
         #   colnamesTMP<-colnames(mrbin.env$bins)[-removeTMP]
         #   mrbin.env$bins<-matrix(mrbin.env$bins[,-removeTMP],nrow=1)
         #   rownames(mrbin.env$bins)<-rownamesTMP
         #   colnames(mrbin.env$bins)<-colnamesTMP
         #} else {
         #   mrbin.env$bins<-mrbin.env$bins[,-removeTMP]
         #}
         mrbin.env$mrbinTMP$binRegions<-mrbin.env$mrbinTMP$binRegions[-removeTMP,]
     }
  }
  mrbin.env$mrbinparam$numberOfFeaturesAfterRemovingAreas<-nrow(mrbin.env$mrbinTMP$binRegions)
 #}
}

#' A function for calculating noise levels.
#'
#' This function calculates noise levels for the current spectrum. Only for use
#' within the mrbin function.
#' @return {None}
#' @export
#' @examples
#' \donttest{ calculateNoise() }

calculateNoise<-function(){
 if(!is.null(mrbin.env$mrbinTMP$currentSpectrum)){
  #if(mrbin.env$mrbinparam$dimension=="2D"){
  #     NMRdataNames<-cbind(apply(mrbin.env$mrbinTMP$binRegions[,3:4],1,mean),apply(mrbin.env$mrbinTMP$binRegions[,1:2],1,mean))
  #}
  #noise_level<-rep(0,length(mrbin.env$mrbinparam$NMRfolders))
  #names(noise_level)<-rownames(mrbin.env$bins)
  #for(j in 1:nrow(mrbin.env$bins)){
    if(mrbin.env$mrbinparam$dimension=="1D"){
         noise_level<-stats::sd(mrbin.env$mrbinTMP$currentSpectrum[
                             which(as.numeric(names(mrbin.env$mrbinTMP$currentSpectrum))<=max(mrbin.env$mrbinparam$noiseRange1d[1:2])&
                                   as.numeric(names(mrbin.env$mrbinTMP$currentSpectrum))>=min(mrbin.env$mrbinparam$noiseRange1d[1:2]))])
    }
    if(mrbin.env$mrbinparam$dimension=="2D"){
         noise_level<-stats::sd(mrbin.env$mrbinTMP$currentSpectrum[
                              which(as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrum))>=min(mrbin.env$mrbinparam$noiseRange2d[3:4])&
                                as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrum))<=max(mrbin.env$mrbinparam$noiseRange2d[3:4])),
                              which(as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrum))<=max(mrbin.env$mrbinparam$noiseRange2d[1:2])&
                                as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrum))>=min(mrbin.env$mrbinparam$noiseRange2d[1:2]))])
    }
  #}
  mrbin.env$mrbinTMP$noise_level_Raw_TMP<-noise_level
 }
}

#' A function for removing bins below noise level.
#'
#' This function checks for each bin (column) whether its level is below the
#' individual noise level times the signal-to-noise ratio. If less than the
#' defined threshold level are above noise*SNR, the whole bin is removed.
#' @return {None}
#' @export
#' @examples
#' mrbin(silent=TRUE,setDefault=TRUE,parameters=list(dimension="1D",
#'                     binwidth1D=0.05,noiseRemoval="No",PQNScaling="No",
#'                     fixNegatives="No",logTrafo="No",PCA="No",verbose=FALSE,
#'                     NMRfolders=c(system.file("extdata/1/10/pdata/10",package="mrbin"),
#'                                 system.file("extdata/2/10/pdata/10",package="mrbin"),
#'                                 system.file("extdata/3/10/pdata/10",package="mrbin"))))
#' removeNoise()

removeNoise<-function(){#remove noise peaks
 if(!is.null(mrbin.env$bins)){
    minimumNumber<-max(1,floor(mrbin.env$mrbinparam$noiseThreshold*nrow(mrbin.env$bins)))
    colnames_NMRdata_no_noise<-NULL
    if(mrbin.env$mrbinparam$dimension=="2D"){
         NMRdataNames<-cbind(apply(mrbin.env$mrbinTMP$binRegions[,3:4],1,mean),apply(mrbin.env$mrbinTMP$binRegions[,1:2],1,mean))
         SNR<-mrbin.env$mrbinparam$signal_to_noise2D
    } else {
         SNR<-mrbin.env$mrbinparam$signal_to_noise1D
    }
    #calculateNoise()
    for(i in 1:ncol(mrbin.env$bins)){#Keep only bins where at least X spectra are > SNR
          if(sum(mrbin.env$bins[,i]>mrbin.env$mrbinparam$noise_level*SNR)>=minimumNumber){
              colnames_NMRdata_no_noise<-c(colnames_NMRdata_no_noise,i)
          }
    }
    if(!is.null(colnames_NMRdata_no_noise)){
        if(nrow(mrbin.env$bins)==1){
            rownamesTMP<-rownames(mrbin.env$bins)
            colnamesTMP<-colnames(mrbin.env$bins)[colnames_NMRdata_no_noise]
            mrbin.env$bins<-matrix(mrbin.env$bins[,colnames_NMRdata_no_noise],nrow=1)
            rownames(mrbin.env$bins)<-rownamesTMP
            colnames(mrbin.env$bins)<-colnamesTMP
        } else {
          mrbin.env$bins<-mrbin.env$bins[,colnames_NMRdata_no_noise]
        }

        mrbin.env$mrbinTMP$binRegions<-mrbin.env$mrbinTMP$binRegions[colnames_NMRdata_no_noise,]
    } else {
        stop("No bins above noise level. Noise removal stopped.\n")
    }
  mrbin.env$mrbinparam$numberOfFeaturesAfterNoiseRemoval<-ncol(mrbin.env$bins)
 }
}

#' A function for cropping HSQC spectra.
#'
#' This function crops HSQC spectra to the region along the diagonal to remove
#' uninformative signals. Will work only for 1H-13C HSQC spectra.
#' @param plot Should a plot of the bins before and after cropping be shown? Defaults to FALSE.
#' @return {None}
#' @export
#' @examples
#' setParam(parameters=list(dimension="2D",binwidth2D=1,binheight=4,cropHSQC="No",PCA="No",
#'          PQNScaling="No",noiseRemoval="No",removeSolvent="No",
#'          NMRfolders=c(system.file("extdata/1/12/pdata/10",package="mrbin"))))
#' mrbinrun()
#' cropNMR()

cropNMR<-function(plot=FALSE){
 #if(!is.null(mrbin.env$bins)){
 if(mrbin.env$mrbinparam$dimension=="2D"){
   selectedCols<-NULL
   coordTmp2<-cbind(apply(mrbin.env$mrbinTMP$binRegions[,3:4],1,mean),apply(mrbin.env$mrbinTMP$binRegions[,1:2],1,mean))
   for(j in 1:nrow(mrbin.env$mrbinTMP$binRegions)){
      coordTmp<-coordTmp2[j,]#1=C,2=H
      if(((coordTmp[2]-mrbin.env$mrbinparam$croptopLeft[2])*(mrbin.env$mrbinparam$cropbottomLeft[1]-mrbin.env$mrbinparam$croptopLeft[1])-
          (coordTmp[1]-mrbin.env$mrbinparam$croptopLeft[1])*(mrbin.env$mrbinparam$cropbottomLeft[2]-mrbin.env$mrbinparam$croptopLeft[2]))<0&
         ((coordTmp[2]-mrbin.env$mrbinparam$croptopRight[2])*(mrbin.env$mrbinparam$cropbottomRight[1]-mrbin.env$mrbinparam$croptopRight[1])-
          (coordTmp[1]-mrbin.env$mrbinparam$croptopRight[1])*(mrbin.env$mrbinparam$cropbottomRight[2]-mrbin.env$mrbinparam$croptopRight[2]))>0
        ){#along diagonal? outer product
          selectedCols<-c(selectedCols,j)
      }
  }
  if(plot){
      graphics::plot(t(matrix(as.numeric(unlist(strsplit(colnames(mrbin.env$bins),","))),
           nrow=2))[,c(2,1)],
           xlim=c(10,0),ylim=c(160,0),xlab="1H",ylab="13C",pch=20,col="black")
      graphics::points(t(matrix(as.numeric(unlist(strsplit(colnames(mrbin.env$bins[,selectedCols]),","))),
           nrow=2))[,c(2,1)],col="red",pch=20)
  }
  #if(nrow(mrbin.env$bins)==1){
  #        rownamesTMP<-rownames(mrbin.env$bins)
  #        colnamesTMP<-colnames(mrbin.env$bins)[selectedCols]
  #        mrbin.env$bins<-matrix(mrbin.env$bins[,selectedCols],nrow=1)
  #        rownames(mrbin.env$bins)<-rownamesTMP
  #        colnames(mrbin.env$bins)<-colnamesTMP
  #} else {
  #  mrbin.env$bins<-mrbin.env$bins[,selectedCols]
  #}
  mrbin.env$mrbinTMP$binRegions<-mrbin.env$mrbinTMP$binRegions[selectedCols,]
  mrbin.env$mrbinparam$numberOfFeaturesAfterCropping<-nrow(mrbin.env$mrbinTMP$binRegions)
  }
 #}
}

#' A function for PQN scaling.
#'
#' This function performs PQN scaling. To further exclude unreliable noise, only
#' the most intense signals are used. For 1H-13C HSQC spectra, most of the sugar
#' regions are excluded to avoid a dominating effect of the multiple sugar
#' signals.
#' @return {None}
#' @export
#' @examples
#' mrbinExample<-mrbin(silent=TRUE,setDefault=TRUE,parameters=list(dimension="1D",
#'                     binwidth1D=0.05,PQNScaling="No",PCA="No",
#'                     NMRfolders=c(system.file("extdata/1/10/pdata/10",package="mrbin"),
#'                                 system.file("extdata/2/10/pdata/10",package="mrbin"),
#'                                 system.file("extdata/3/10/pdata/10",package="mrbin"))))
#' PQNScaling()

PQNScaling<-function(){#Scale to PQN
 if(!is.null(mrbin.env$bins)){
  if(nrow(mrbin.env$bins)>1){
    #Create synthetic median spectrum by averaging all spectra
    NMRdataTmp<-rbind(mrbin.env$bins,apply(mrbin.env$bins,2,mean))
    rownames(NMRdataTmp)[nrow(NMRdataTmp)]<-"Median"
    medianSample<-nrow(NMRdataTmp)
    #For HSQC spetra: Remove most sugar signals to get a better fold change estimate
    if(mrbin.env$mrbinparam$dimension == "2D" & mrbin.env$mrbinparam$cropHSQC=="Yes") {
        selectedCols<-NULL
        for(j in 1:ncol(NMRdataTmp)){
            coordTmp<-cbind(apply(mrbin.env$mrbinTMP$binRegions[,3:4],1,mean),apply(mrbin.env$mrbinTMP$binRegions[,1:2],1,mean))[j,]#1=C,2=H
                   if(!(coordTmp[2]>mrbin.env$mrbinparam$PQNsugarArea[2]&coordTmp[2]<mrbin.env$mrbinparam$PQNsugarArea[1]&
                        coordTmp[1]>mrbin.env$mrbinparam$PQNsugarArea[3]&coordTmp[1]<mrbin.env$mrbinparam$PQNsugarArea[4])){
                            selectedCols<-c(selectedCols,j)
                   }
        }
        NMRdataTmp2<-NMRdataTmp[,selectedCols]
    } else { #1D spectra
      NMRdataTmp2<-NMRdataTmp
    }
    #Calculate fold changes versus reference sample
    NMRdataTmp_scaledMedian<-NMRdataTmp
    PQNminimumFeatures<-min(max(mrbin.env$mrbinparam$PQNminimumFeatures,floor(.25*ncol(NMRdataTmp2))),ncol(NMRdataTmp2))
    #if(mrbin.env$mrbinparam$PQNminimumFeatures>ncol(NMRdataTmp2)) mrbin.env$mrbinparam$PQNminimumFeatures<-ncol(NMRdataTmp2)
    medianFoldChanges<-rep(0,nrow(NMRdataTmp_scaledMedian))
    names(medianFoldChanges)<-rownames(NMRdataTmp_scaledMedian)
    if(mrbin.env$mrbinparam$PQNshowHist){#Plot distribution of fold changes per sample
      oldpar<-graphics::par("mar","mfrow")
      on.exit(graphics::par(oldpar))
      graphics::par(mfrow=c(ceiling(sqrt(nrow(NMRdataTmp2))),ceiling(sqrt(nrow(NMRdataTmp2)))))
      graphics::par(mar=c(.1,.1,3,.1))
    }
    for(i in 1:nrow(NMRdataTmp2)){#scale to spectral area, use only points above X% quantile for better reliability
        overlapAboveNoise<- sort(NMRdataTmp2[i,],index.return=TRUE,decreasing=TRUE)$ix[1:PQNminimumFeatures]#$ix returns the index for matrix data
        medianFoldChanges[i]<-stats::median(NMRdataTmp2[i,overlapAboveNoise]/
                                     NMRdataTmp2[medianSample,overlapAboveNoise])
        NMRdataTmp_scaledMedian[i,]<-NMRdataTmp[i,]/medianFoldChanges[i]
      if(mrbin.env$mrbinparam$PQNshowHist){#Plot distribution of fold changes per sample
            graphics::hist(NMRdataTmp2[i,overlapAboveNoise]/NMRdataTmp2[medianSample,overlapAboveNoise],breaks=60,
            main=rownames(NMRdataTmp2)[i],xlab="",ylab="")
            graphics::lines(rep(medianFoldChanges[i],2),c(0,20000),col='red')
      }
    }
    #Remove reference sample from list
    NMRdataTmp_scaledMedian<-NMRdataTmp_scaledMedian[-nrow(NMRdataTmp_scaledMedian),]
    if(nrow(mrbin.env$bins)==1){
            rownamesTMP<-rownames(mrbin.env$bins)
            colnamesTMP<-colnames(mrbin.env$bins)
            mrbin.env$bins<-matrix(NMRdataTmp_scaledMedian,nrow=1)
            rownames(mrbin.env$bins)<-rownamesTMP
            colnames(mrbin.env$bins)<-colnamesTMP
    } else {
      mrbin.env$bins<-NMRdataTmp_scaledMedian
    }
    mrbin.env$mrbinparam$medians<-medianFoldChanges
   } else {
      warning("Too few samples to perform PQN normalization.")
   }
 }
}

#' A function for plotting the current bin positions
#'
#' This function plots the current binned spectrum.
#' @param showtitle If TRUE, the spectrum title is shown. Defaults to TRUE.
#' @return {None}
#' @export
#' @examples
#' mrbinExample<-mrbin(silent=TRUE,setDefault=TRUE,parameters=list(dimension="1D",
#'                     binwidth1D=0.05,PQNScaling="No",PCA="No",
#'                     NMRfolders=c(system.file("extdata/1/10/pdata/10",package="mrbin"))))
#' plotBins()

plotBins<-function(showtitle=TRUE){#Plot bin positions
 if(!is.null(mrbin.env$bins)){
   devAskNewPage(ask = FALSE)
   if(mrbin.env$mrbinparam$dimension=="2D"){
      NMRdata_select<-mrbin.env$bins[mrbin.env$mrbinTMP$currentSpectrumName,]
      NMRdata_select<-NMRdata_select[order(abs(NMRdata_select))]
      if(showtitle) titleTMP<-mrbin.env$mrbinTMP$currentSpectrumName
      if(!showtitle) titleTMP<-""
      allBins2<-cbind(apply(mrbin.env$mrbinTMP$binRegions[,3:4],1,mean),apply(mrbin.env$mrbinTMP$binRegions[,1:2],1,mean))
      if(showtitle){
        colorRampHSQC<-grDevices::colorRamp(c("blue","green","orange","red","red"))
        colorTMP<-grDevices::rgb(colorRampHSQC(NMRdata_select/max(abs(NMRdata_select))/2+.5),maxColorValue=255)
      } else {
        colorTMP<-"white"
      }
      graphics::plot(allBins2[,2],allBins2[,1],ylim=mrbin.env$mrbinparam$binRegion[4:3],xlim=mrbin.env$mrbinparam$binRegion[1:2],
           cex=.5,pch=15,main=titleTMP,col=colorTMP, ask=FALSE)
      graphics::rect(xleft=mrbin.env$mrbinTMP$binRegions[,1], ybottom=mrbin.env$mrbinTMP$binRegions[,4],
                     xright=mrbin.env$mrbinTMP$binRegions[,2], ytop=mrbin.env$mrbinTMP$binRegions[,3], col = "green", border = NA)
      graphics::box()
      utils::flush.console()
   }
   if(mrbin.env$mrbinparam$dimension=="1D") {
      colorRampHSQC<-grDevices::colorRamp(c("blue","green","orange","red","red"))
      NMRdata_select<-mrbin.env$bins[mrbin.env$mrbinTMP$currentSpectrumName,]
      NMRdata_select<-NMRdata_select[order(abs(NMRdata_select))]
      if(showtitle) titleTMP<-mrbin.env$mrbinTMP$currentSpectrumName
      if(!showtitle) titleTMP<-""
      graphics::plot(NULL,NULL,ylim=c(0,1),xlim=mrbin.env$mrbinparam$binRegion[1:2],yaxt="n",
           main=titleTMP, ask=FALSE)
      graphics::rect(xleft=mrbin.env$mrbinTMP$binRegions[,1], ybottom=0, xright=mrbin.env$mrbinTMP$binRegions[,2], ytop=2.0, col = "green", border = NA)
      graphics::lines(as.numeric(names(mrbin.env$mrbinTMP$currentSpectrum))
           #apply(mrbin.env$mrbinTMP$binRegions[,1:2],1,mean)
           ,mrbin.env$mrbinTMP$currentSpectrum/(sort(mrbin.env$mrbinTMP$currentSpectrum)[ceiling(.9999*length(mrbin.env$mrbinTMP$currentSpectrum))]),
           col="black")
      graphics::box()
      utils::flush.console()
   }
 }
}

#' A function for plotting quality indicators, including PCA plots.
#'
#' This function plots boxplots (bin-wise and sample-wise) as visual quality indicators. It also performs PCA, then plots PC1 and PC2 and loading plots.
#' @return {None}
#' @export
#' @examples
#' mrbinExample<-mrbin(silent=TRUE,setDefault=TRUE,parameters=list(dimension="1D",
#'                     binwidth1D=0.05,PCA="No",
#'                     NMRfolders=c(system.file("extdata/1/10/pdata/10",package="mrbin"),
#'                                 system.file("extdata/2/10/pdata/10",package="mrbin"),
#'                                 system.file("extdata/3/10/pdata/10",package="mrbin"))))
#' plotResults()

plotResults<-function(){
 if(!is.null(mrbin.env$bins)){
    if(is.null(mrbin.env$mrbinparam$Factors)){
      FactorsTMP<-factor(rep("Group 0",nrow(mrbin.env$bins)))
    } else {
      FactorsTMP<-mrbin.env$mrbinparam$Factors
    }
    colorPalette<-grDevices::rainbow(length(levels(FactorsTMP)))
    mrbin.env$mrbinTMP$PCA<-stats::prcomp(mrbin.env$bins)
    oldpar<-graphics::par("mar","mfrow","mgp")
    on.exit(graphics::par(oldpar))
    devAskNewPage(ask = FALSE)
    graphics::par(mfrow=c(2,3),mar=c(3.1,2,2,0.5))
    plotBins(showtitle=FALSE)
    axisCex1<-.1
    if(ncol(mrbin.env$bins)<100) axisCex1<-.25
    if(ncol(mrbin.env$bins)<25) axisCex1<-.5
    if(ncol(mrbin.env$bins)<15) axisCex1<-1
    axisCex2<-.1
    if(nrow(mrbin.env$bins)<100) axisCex2<-.25
    if(nrow(mrbin.env$bins)<25) axisCex2<-.5
    if(nrow(mrbin.env$bins)<15) axisCex2<-1
    graphics::boxplot(mrbin.env$bins,main="",xlab="Bins",ylab="",boxwex=1,ask=FALSE,xaxt="n")
    graphics::axis(1,las=2,cex.axis=axisCex1,at=1:ncol(mrbin.env$bins),labels=colnames(mrbin.env$bins))
    graphics::boxplot(t(mrbin.env$bins),main="",xlab="Samples",ylab="",boxwex=1,ask=FALSE,xaxt="n")
    graphics::axis(1,las=2,cex.axis=axisCex2,at=1:nrow(mrbin.env$bins),labels=substr(rownames(mrbin.env$bins),1,mrbin.env$mrbinparam$PCAtitlelength))
    utils::flush.console()
    if(nrow(mrbin.env$bins)>1){
      graphics::par(mgp=c(1,1,0))
      graphics::plot(mrbin.env$mrbinTMP$PCA$rotation,pch=16,cex=.75,main="PCA Loadings Plot", ask=FALSE, xaxt='n', yaxt='n')
      graphics::text(mrbin.env$mrbinTMP$PCA$rotation,labels=substr(rownames(mrbin.env$mrbinTMP$PCA$rotation),1,12),pos=4,cex=1)
      numlevels<-NULL
      if(mrbin.env$mrbinparam$defineGroups=="Yes"){
        for(i in 1:nlevels(FactorsTMP)) numlevels<-c(numlevels,as.numeric(
                         FactorsTMP[which(FactorsTMP==levels(FactorsTMP)[i])][1]))
      }
      addTMP1<-.1*(max(mrbin.env$mrbinTMP$PCA$x[,1])-min(mrbin.env$mrbinTMP$PCA$x[,1]))
      addTMP2<-.1*(max(mrbin.env$mrbinTMP$PCA$x[,2])-min(mrbin.env$mrbinTMP$PCA$x[,2]))
      if(mrbin.env$mrbinparam$defineGroups=="Yes"){
        PCAFactors<-as.numeric(FactorsTMP)
      } else {
        PCAFactors<-rep(1,nrow(mrbin.env$bins))
      }
      graphics::plot(mrbin.env$mrbinTMP$PCA$x[,1],mrbin.env$mrbinTMP$PCA$x[,2],
           xlim=c(min(mrbin.env$mrbinTMP$PCA$x[,1])-addTMP1,max(mrbin.env$mrbinTMP$PCA$x[,1])+addTMP1),
           ylim=c(min(mrbin.env$mrbinTMP$PCA$x[,2])-addTMP2,max(mrbin.env$mrbinTMP$PCA$x[,2])+addTMP2),
           pch=PCAFactors+14, xaxt='n', yaxt='n',
           col=colorPalette[PCAFactors],
           main="PCA",
           xlab=paste("PC1 (",round(100*mrbin.env$mrbinTMP$PCA$sdev[1]/sum(mrbin.env$mrbinTMP$PCA$sdev),1),"%)",sep=""),
           ylab=paste("PC2 (",round(100*mrbin.env$mrbinTMP$PCA$sdev[2]/sum(mrbin.env$mrbinTMP$PCA$sdev),1),"%)",sep="")
           ,cex=.75, ask=FALSE)
      graphics::text(mrbin.env$mrbinTMP$PCA$x,labels=paste(substr(rownames(mrbin.env$mrbinTMP$PCA$x),1,mrbin.env$mrbinparam$PCAtitlelength)),pos=3,cex=1,
             col=colorPalette[PCAFactors])
      graphics::plot(NULL,NULL,xaxt="n",yaxt="n",ask=FALSE,axes=FALSE,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",main="")
      graphics::par(mar=c(4.2,0,2.8,0.5))
      if(mrbin.env$mrbinparam$defineGroups=="Yes") graphics::legend("left", legend=levels(FactorsTMP),
              col=colorPalette[numlevels],pch=numlevels+14,cex=1,bg="white")
      utils::flush.console()
    } else {
       warning("Too few samples to perform PCA.")
    }
 }
}

#' A function for plotting NMR spectra.
#'
#' This function plots the current NMR spectrum. If no parameters are provided, parameters
#' are read from the mrbin.env environment variables, set by mrbin.
#' To change the plot, use zoom(),
#' zoomIn(), zoomOut(), intPlus(), intMin(), left(), right().
#' For 2D data use additionally: contMin(), contPlus(), up(), down()
#' @param region A vector defining the plot region (left, right, top, bottom)
#' @param rectangleRegions A 4-column matrix defining areas where to plot rectangles
#' @param rectangleColors Define colors for the rectangles
#' @param polygonRegion Defines 4 corners of a polygon to be plotted
#' @param color Defines the color of the spectrum plot. If NULL, a rainbow theme is used for 2D NMR
#' @param add If TRUE, additional spectrum plots are overlaid with the current plot
#' @param manualScale If TRUE, scaling factor is taken from environment variables.
#' @param plotTitle Defines the main title of the plot
#' @return {None}
#' @export
#' @examples
#' mrbin(silent=TRUE,setDefault=TRUE,parameters=list(dimension="1D",binwidth1D=.1,
#'          PQNScaling="No",
#'          fixNegatives="No",logTrafo="No",PCA="No",verbose=FALSE,
#'          NMRfolders=system.file("extdata/1/10/pdata/10",package="mrbin")))
#' plotNMR()

plotNMR<-function(region=NULL,rectangleRegions=NULL,
                   rectangleColors=c("green","orange","blue","red","yellow","gray","purple"),
                   polygonRegion=NULL,
                   color=NULL,add=FALSE,
                   manualScale=TRUE,plotTitle=""){
 if(!is.null(mrbin.env$mrbinTMP$currentSpectrum)){
   devAskNewPage(ask = FALSE)
   if(is.null(region)){
       if(is.null(mrbin.env$mrbinplot$plotRegion)){
            if(is.matrix(mrbin.env$mrbinTMP$currentSpectrum)){
                mrbin.env$mrbinplot$plotRegion<-
                          c(max(as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrum))),
                            min(as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrum))),
                            min(as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrum))),
                            max(as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrum))))
            } else {
                mrbin.env$mrbinplot$plotRegion<-
                         c(max(as.numeric(names(mrbin.env$mrbinTMP$currentSpectrum))),
                         min(as.numeric(names(mrbin.env$mrbinTMP$currentSpectrum))),
                         -10,160)
            }
       }
       region<-mrbin.env$mrbinplot$plotRegion
   }
   if(length(region)==1){
     if(region=="all"){
          if(is.matrix(mrbin.env$mrbinTMP$currentSpectrum)){
                region<-c(max(as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrum))),
                            min(as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrum))),
                            min(as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrum))),
                            max(as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrum))))
          } else {
                region<-c(max(as.numeric(names(mrbin.env$mrbinTMP$currentSpectrum))),
                         min(as.numeric(names(mrbin.env$mrbinTMP$currentSpectrum))),
                         -10,160)
          }
     }
   }
   if(is.matrix(mrbin.env$mrbinTMP$currentSpectrum)){#2D spectra
      if(is.null(color)) color<-rev(grDevices::rainbow(mrbin.env$mrbinplot$nContours))
      spectrumTMP<-mrbin.env$mrbinTMP$currentSpectrum[which(as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrum))<
                           region[4]&
                          as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrum))>=region[3]),
                         which(as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrum))>=region[2]&
                         as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrum))<region[1])
                         ]
      if(manualScale){
        if(sum(spectrumTMP<(mrbin.env$mrbinplot$lowestContour*max(mrbin.env$mrbinTMP$currentSpectrum)))>0){
             spectrumTMP[spectrumTMP<(mrbin.env$mrbinplot$lowestContour*max(mrbin.env$mrbinTMP$currentSpectrum))]<-0
        }
      } else {
        spectrumTMP[spectrumTMP<=0]<-1e-8
        spectrumTMP<-log(100000*spectrumTMP/max(spectrumTMP))
      }
      displaySize<-512#Reduce resolution for faster plotting of high-res spectra
      if(nrow(spectrumTMP)>(2*displaySize)){
         sizeRegion1<- max(2,floor(nrow(spectrumTMP)/displaySize))
         nRegion1<-max(2,floor(nrow(spectrumTMP)/sizeRegion1))
         spectrumTMP2<-matrix(rep(0,nRegion1*ncol(spectrumTMP)),nrow=nRegion1)
         bottomTMP<-mean(as.numeric(rownames(spectrumTMP)[1:sizeRegion1]))
         topTMP<-mean(as.numeric(rownames(spectrumTMP)[((nRegion1-1)*sizeRegion1)+1:sizeRegion1]))
         for(i in 1:nRegion1) {
            spectrumTMP2[i,]<-apply(spectrumTMP[((i-1)*sizeRegion1)+1:sizeRegion1,],2,max)
         }
         rownames(spectrumTMP2)<-bottomTMP+1:nRegion1*(topTMP-bottomTMP)/nRegion1
         colnames(spectrumTMP2)<-colnames(spectrumTMP)
         spectrumTMP<-spectrumTMP2
      }
      if(ncol(spectrumTMP)>(2*displaySize)){
         sizeRegion2<- max(2,floor(ncol(spectrumTMP)/displaySize))
         nRegion2<-max(2,floor(ncol(spectrumTMP)/sizeRegion2))
         spectrumTMP2<-matrix(rep(0,nRegion2*nrow(spectrumTMP)),ncol=nRegion2)
         leftTMP<-mean(as.numeric(colnames(spectrumTMP)[1:sizeRegion2]))
         rightTMP<-mean(as.numeric(colnames(spectrumTMP)[((nRegion2-1)*sizeRegion2)+1:sizeRegion2]))
         for(i in 1:nRegion2) {
            spectrumTMP2[,i]<-apply(spectrumTMP[,((i-1)*sizeRegion2)+1:sizeRegion2],1,max)
         }
         colnames(spectrumTMP2)<-leftTMP+1:nRegion2*(rightTMP-leftTMP)/nRegion2
         rownames(spectrumTMP2)<-rownames(spectrumTMP)
         spectrumTMP<-spectrumTMP2
      }
      if(!manualScale){
        while(sum(spectrumTMP>0)>.0025*length(spectrumTMP)){
            spectrumTMP[spectrumTMP>0][spectrumTMP[spectrumTMP>0]<sort(
                       spectrumTMP[spectrumTMP>0])[floor(.5*length(spectrumTMP[spectrumTMP>0]))]]<-0
        }
      }
      options(max.contour.segments=1000)
      if(!add){
          graphics::plot(NULL,NULL,
            type="l",xlim=c(-region[1],-region[2]),
            ylim=c(-region[4],-region[3]),main=plotTitle,xaxt="n",yaxt="n",
            xlab="Chemical shift [ppm]",ylab="Chemical shift [ppm]")
          magnitude2<-10^round(log(max(as.numeric(colnames(spectrumTMP)))-
                      min(as.numeric(colnames(spectrumTMP))),base=10))/10
          magnitude1<-10^round(log(max(as.numeric(rownames(spectrumTMP)))-
                       min(as.numeric(rownames(spectrumTMP))),base=10))/10
          graphics::axis(2,
               at=-(0:100*magnitude1+floor(min(as.numeric(rownames(spectrumTMP)))/magnitude1)*magnitude1)
               ,labels=(0:100*magnitude1+floor(min(as.numeric(rownames(spectrumTMP)))/magnitude1)*magnitude1)
               )
          graphics::axis(1,
               at=(0:100*magnitude2+floor(min(-as.numeric(colnames(spectrumTMP)))/magnitude2)*magnitude2)
               ,labels=-(0:100*magnitude2+floor(min(-as.numeric(colnames(spectrumTMP)))/magnitude2)*magnitude2)
               )
      }
      if(!is.null(rectangleRegions)){
          graphics::rect(xleft=-rectangleRegions[,1], ybottom=-rectangleRegions[,4],
                         xright=-rectangleRegions[,2], ytop=-rectangleRegions[,3],
                         col = rectangleColors, border = NA)
          graphics::box()
      }
      if(!is.null(polygonRegion)){
          graphics::polygon(-1*polygonRegion[,2],-1*polygonRegion[,1],
                    #c(0,10,11,1),c(6e+5,0,1000,6e+5)
                    col=rectangleColors, border = NA)
          graphics::box()
      }
      suppressWarnings(graphics::contour(x = -as.numeric(colnames(spectrumTMP)),
        y = -as.numeric(rownames(spectrumTMP)),
        z = t(spectrumTMP)*mrbin.env$mrbinplot$intensityScale,
        levels =  (mrbin.env$mrbinplot$lowestContour+
                      .8*(1:mrbin.env$mrbinplot$nContours-1)/
                      (mrbin.env$mrbinplot$nContours)*(1-
                      mrbin.env$mrbinplot$lowestContour))* max(spectrumTMP)
        ,drawlabels=FALSE,col=color,lwd=1,add=TRUE))
      utils::flush.console()
   } else {  #1D
      if(is.null(color)) color<-"black"
      if(manualScale){
        ymin<-min(mrbin.env$mrbinTMP$currentSpectrum)
        ymax<-max(mrbin.env$mrbinTMP$currentSpectrum)/mrbin.env$mrbinplot$intensityScale
      } else {
        ymin<-0
        spectrumTMP<-mrbin.env$mrbinTMP$currentSpectrum[as.numeric(names(mrbin.env$mrbinTMP$currentSpectrum<region[1]))&
                   as.numeric(names(mrbin.env$mrbinTMP$currentSpectrum))>region[2]]
        ymax<-sort(spectrumTMP,decreasing=TRUE)[ceiling(.01*length(spectrumTMP))]
      }
      if(!add){
          graphics::plot(NULL,NULL,
            type="l",xlim=c(region[1],region[2]),
            ylim=c(ymin,ymax),main=plotTitle,
            xlab="Chemical shift [ppm]",ylab="Intensity")
      }
      if(!is.null(rectangleRegions)){
          graphics::rect(xleft=rectangleRegions[,1], ybottom=0, xright=rectangleRegions[,2], ytop=ymax*2,
                    col = rectangleColors, border = NA)
          graphics::box()
      }
      graphics::lines(as.numeric(names(mrbin.env$mrbinTMP$currentSpectrum)),mrbin.env$mrbinTMP$currentSpectrum,
                col=color)
      utils::flush.console()
    }
 }
}

#' A function for changing plotNMR plots.
#'
#' This function increases the intensity of the current NMR spectrum plot.
#' @return {None}
#' @export
#' @examples
#' mrbin(silent=TRUE,setDefault=TRUE,parameters=list(dimension="1D",binwidth1D=.1,
#'          PQNScaling="No",
#'          fixNegatives="No",logTrafo="No",PCA="No",verbose=FALSE,
#'          NMRfolders=system.file("extdata/1/10/pdata/10",package="mrbin")))
#' plotNMR()
#' intPlus()

intPlus<-function(){#increase plot intensity
 if(!is.null(mrbin.env$mrbinTMP$currentSpectrum)){
   mrbin.env$mrbinplot$intensityScale<-mrbin.env$mrbinplot$intensityScale*2
   plotNMR()
 }
}

#' A function for changing plotNMR plots.
#'
#' This function decreases the intensity of the current NMR spectrum plot.
#' @return {None}
#' @export
#' @examples
#' mrbin(silent=TRUE,setDefault=TRUE,parameters=list(dimension="1D",binwidth1D=.1,
#'          PQNScaling="No",
#'          fixNegatives="No",logTrafo="No",PCA="No",verbose=FALSE,
#'          NMRfolders=system.file("extdata/1/10/pdata/10",package="mrbin")))
#' plotNMR()
#' intMin()

intMin<-function(){#decrease plot intensity
 if(!is.null(mrbin.env$mrbinTMP$currentSpectrum)){
   mrbin.env$mrbinplot$intensityScale<-mrbin.env$mrbinplot$intensityScale*0.5
   plotNMR()
 }
}

#' A function for changing plotNMR plots.
#'
#' This function increases the minimum contour level of the current 2D NMR
#' spectrum plot.
#' @return {None}
#' @export
#' @examples
#' mrbin(silent=TRUE,parameters=list(dimension="2D",binwidth2D=0.5,
#'          binheight=20,PQNScaling="No",referenceScaling="No",binRegion=c(8,2,20,140),
#'          fixNegatives="No",logTrafo="No",PCA="No",verbose=FALSE,saveFiles="No",
#'          NMRfolders=c(system.file("extdata/1/12/pdata/10",package="mrbin"))))
#' plotNMR()
#' contPlus()

contPlus<-function(){#decrease plot intensity
 if(!is.null(mrbin.env$mrbinTMP$currentSpectrum)){
   mrbin.env$mrbinplot$lowestContour<-mrbin.env$mrbinplot$lowestContour*1.5
   plotNMR()
 }
}

#' A function for changing plotNMR plots.
#'
#' This function decreases the minimum contour level of the current 2D NMR
#' spectrum plot.
#' @return {None}
#' @export
#' @examples
#' mrbin(silent=TRUE,parameters=list(dimension="2D",binwidth2D=0.5,
#'          binheight=20,PQNScaling="No",referenceScaling="No",binRegion=c(8,2,20,140),
#'          fixNegatives="No",logTrafo="No",PCA="No",verbose=FALSE,saveFiles="No",
#'          NMRfolders=c(system.file("extdata/1/12/pdata/10",package="mrbin"))))
#' plotNMR()
#' contMin()

contMin<-function(){#decrease plot intensity
 if(!is.null(mrbin.env$mrbinTMP$currentSpectrum)){
   mrbin.env$mrbinplot$lowestContour<-mrbin.env$mrbinplot$lowestContour*0.75
   plotNMR()
 }
}

#' A function for changing plotNMR plots.
#'
#' This function changes the plot region of the current NMR plot.
#' @param left New left boundary
#' @param right New right boundary
#' @param top New top boundary
#' @param bottom New bottom boundary
#' @return {None}
#' @export
#' @examples
#' # Ask for user input
#' \donttest{ zoom() }
#' mrbin(silent=TRUE,setDefault=TRUE,parameters=list(dimension="1D",binwidth1D=.1,
#'          PQNScaling="No",
#'          fixNegatives="No",logTrafo="No",PCA="No",verbose=FALSE,
#'          NMRfolders=system.file("extdata/1/10/pdata/10",package="mrbin")))
#' plotNMR()
#' zoom(left=4.6,right=2,top=10,bottom=150)

zoom<-function(left=NULL,right=NULL,top=NULL,bottom=NULL){
 if(!is.null(mrbin.env$mrbinTMP$currentSpectrum)){
   if(is.null(left)) stop("Please set left limit\n")
   if(is.null(right)) stop("Please set right limit\n")
   if(is.matrix(mrbin.env$mrbinTMP$currentSpectrum)){
      if(is.null(top)) stop("Please set top limit\n")
      if(is.null(bottom)) stop("Please set bottom limit\n")
   }
   mrbin.env$mrbinplot$plotRegion[1]<-left
   mrbin.env$mrbinplot$plotRegion[2]<-right
   mrbin.env$mrbinplot$plotRegion[3]<-top
   mrbin.env$mrbinplot$plotRegion[4]<-bottom
   plotNMR()
 }
}

#' A function for changing plotNMR plots.
#'
#' This function zooms into the plot region of the current NMR plot.
#' @return {None}
#' @export
#' @examples
#' mrbin(silent=TRUE,setDefault=TRUE,parameters=list(dimension="1D",binwidth1D=.1,
#'          PQNScaling="No",
#'          fixNegatives="No",logTrafo="No",PCA="No",verbose=FALSE,
#'          NMRfolders=system.file("extdata/1/10/pdata/10",package="mrbin")))
#' plotNMR()
#' zoomIn()

zoomIn<-function(){#Zoom into NMR spectrum plot
 if(!is.null(mrbin.env$mrbinTMP$currentSpectrum)){
   if(mrbin.env$mrbinparam$dimension=="1D"){
       leftMax<-max(as.numeric(names(mrbin.env$mrbinTMP$currentSpectrum)))
       rightMax<-min(as.numeric(names(mrbin.env$mrbinTMP$currentSpectrum)))
       topMax<--10
       bottomMax<-160
   }
   if(mrbin.env$mrbinparam$dimension=="2D"){
       leftMax<-max(as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrum)))
       rightMax<-min(as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrum)))
       topMax<-min(as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrum)))
       bottomMax<-max(as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrum)))
   }
   mrbin.env$mrbinplot$plotRegion[1]<-min(leftMax,
                             mrbin.env$mrbinplot$plotRegion[1]-(mrbin.env$mrbinplot$plotRegion[1]-mrbin.env$mrbinplot$plotRegion[2])/4)
   mrbin.env$mrbinplot$plotRegion[2]<-max(rightMax,
                             mrbin.env$mrbinplot$plotRegion[1]-(mrbin.env$mrbinplot$plotRegion[1]-mrbin.env$mrbinplot$plotRegion[2])*3/4)
   mrbin.env$mrbinplot$plotRegion[3]<-max(topMax,
                             mrbin.env$mrbinplot$plotRegion[4]-(mrbin.env$mrbinplot$plotRegion[4]-mrbin.env$mrbinplot$plotRegion[3])*3/4)
   mrbin.env$mrbinplot$plotRegion[4]<-min(bottomMax,
                             mrbin.env$mrbinplot$plotRegion[4]-(mrbin.env$mrbinplot$plotRegion[4]-mrbin.env$mrbinplot$plotRegion[3])/4)
   plotNMR()
 }
}

#' A function for changing plotNMR plots.
#'
#' This function zooms out from the plot region of the current NMR plot.
#' @return {None}
#' @export
#' @examples
#' mrbin(silent=TRUE,setDefault=TRUE,parameters=list(dimension="1D",binwidth1D=.1,
#'          PQNScaling="No",
#'          fixNegatives="No",logTrafo="No",PCA="No",verbose=FALSE,
#'          NMRfolders=system.file("extdata/1/10/pdata/10",package="mrbin")))
#' plotNMR()
#' zoomIn()
#' zoomOut()

zoomOut<-function(){#Zoom out from NMR spectrum plot
 if(!is.null(mrbin.env$mrbinTMP$currentSpectrum)){
   if(mrbin.env$mrbinparam$dimension=="1D"){
       leftMax<-max(as.numeric(names(mrbin.env$mrbinTMP$currentSpectrum)))
       rightMax<-min(as.numeric(names(mrbin.env$mrbinTMP$currentSpectrum)))
       topMax<--10
       bottomMax<-160
   }
   if(mrbin.env$mrbinparam$dimension=="2D"){
       leftMax<-max(as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrum)))
       rightMax<-min(as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrum)))
       topMax<-min(as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrum)))
       bottomMax<-max(as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrum)))
   }
   mrbin.env$mrbinplot$plotRegion[1]<-min(leftMax,
                             mrbin.env$mrbinplot$plotRegion[1]+(mrbin.env$mrbinplot$plotRegion[1]-mrbin.env$mrbinplot$plotRegion[2])/2)
   mrbin.env$mrbinplot$plotRegion[2]<-max(rightMax,
                             mrbin.env$mrbinplot$plotRegion[1]-(mrbin.env$mrbinplot$plotRegion[1]-mrbin.env$mrbinplot$plotRegion[2])*3/2)
   mrbin.env$mrbinplot$plotRegion[3]<-max(topMax,
                             mrbin.env$mrbinplot$plotRegion[4]-(mrbin.env$mrbinplot$plotRegion[4]-mrbin.env$mrbinplot$plotRegion[3])*3/2)
   mrbin.env$mrbinplot$plotRegion[4]<-min(bottomMax,
                             mrbin.env$mrbinplot$plotRegion[4]+(mrbin.env$mrbinplot$plotRegion[4]-mrbin.env$mrbinplot$plotRegion[3])/2)
   plotNMR()
 }
}

#' A function for changing plotNMR plots.
#'
#' This function moves left the plot region of the current NMR plot.
#' @return {None}
#' @export
#' @examples
#' mrbin(silent=TRUE,parameters=list(dimension="1D",binwidth1D=.5,
#'          PQNScaling="No",saveFiles="No",referenceScaling="No",
#'          fixNegatives="No",logTrafo="No",PCA="No",verbose=FALSE,
#'          NMRfolders=system.file("extdata/1/10/pdata/10",package="mrbin")))
#' plotNMR()
#' zoomIn()
#' left()

left<-function(){
 if(!is.null(mrbin.env$mrbinTMP$currentSpectrum)){
   if(mrbin.env$mrbinparam$dimension=="1D"){
       leftMax<-max(as.numeric(names(mrbin.env$mrbinTMP$currentSpectrum)))
       rightMax<-min(as.numeric(names(mrbin.env$mrbinTMP$currentSpectrum)))
   }
   if(mrbin.env$mrbinparam$dimension=="2D"){
       leftMax<-max(as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrum)))
       rightMax<-min(as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrum)))
   }
   mrbin.env$mrbinplot$plotRegion[1]<-min(leftMax,
                             mrbin.env$mrbinplot$plotRegion[1]+(mrbin.env$mrbinplot$plotRegion[1]-mrbin.env$mrbinplot$plotRegion[2])*.1)
   mrbin.env$mrbinplot$plotRegion[2]<-max(rightMax,
                             mrbin.env$mrbinplot$plotRegion[2]+(mrbin.env$mrbinplot$plotRegion[1]-mrbin.env$mrbinplot$plotRegion[2])*.1)
   plotNMR()
 }
}

#' A function for changing plotNMR plots.
#'
#' This function moves right the plot region of the current NMR plot.
#' @return {None}
#' @export
#' @examples
#' mrbin(silent=TRUE,parameters=list(dimension="1D",binwidth1D=.5,
#'          PQNScaling="No",saveFiles="No",referenceScaling="No",
#'          fixNegatives="No",logTrafo="No",PCA="No",verbose=FALSE,
#'          NMRfolders=system.file("extdata/1/10/pdata/10",package="mrbin")))
#' plotNMR()
#' zoomIn()
#' right()

right<-function(){
 if(!is.null(mrbin.env$mrbinTMP$currentSpectrum)){
   if(mrbin.env$mrbinparam$dimension=="1D"){
       leftMax<-max(as.numeric(names(mrbin.env$mrbinTMP$currentSpectrum)))
       rightMax<-min(as.numeric(names(mrbin.env$mrbinTMP$currentSpectrum)))
   }
   if(mrbin.env$mrbinparam$dimension=="2D"){
       leftMax<-max(as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrum)))
       rightMax<-min(as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrum)))
   }
   mrbin.env$mrbinplot$plotRegion[1]<-min(leftMax,
                             mrbin.env$mrbinplot$plotRegion[1]-(mrbin.env$mrbinplot$plotRegion[1]-mrbin.env$mrbinplot$plotRegion[2])*.1)
   mrbin.env$mrbinplot$plotRegion[2]<-max(rightMax,
                             mrbin.env$mrbinplot$plotRegion[2]-(mrbin.env$mrbinplot$plotRegion[1]-mrbin.env$mrbinplot$plotRegion[2])*.1)
   plotNMR()
 }
}

#' A function for changing plotNMR plots.
#'
#' This function moves down the plot region of the current NMR plot (only 2D).
#' @return {None}
#' @export
#' @examples
#' mrbin(silent=TRUE,parameters=list(dimension="2D",binwidth2D=0.5,
#'          binheight=20,PQNScaling="No",referenceScaling="No",binRegion=c(8,2,20,140),
#'          fixNegatives="No",logTrafo="No",PCA="No",verbose=FALSE,saveFiles="No",
#'          NMRfolders=c(system.file("extdata/1/12/pdata/10",package="mrbin"))))
#' plotNMR()
#' zoomIn()
#' down()

down<-function(){
 if(!is.null(mrbin.env$mrbinTMP$currentSpectrum)){
   if(mrbin.env$mrbinparam$dimension=="2D"){
       topMax<-min(as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrum)))
       bottomMax<-max(as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrum)))
       mrbin.env$mrbinplot$plotRegion[3]<-max(topMax,
                                 mrbin.env$mrbinplot$plotRegion[3]-(mrbin.env$mrbinplot$plotRegion[3]-mrbin.env$mrbinplot$plotRegion[4])*.1)
       mrbin.env$mrbinplot$plotRegion[4]<-min(bottomMax,
                                 mrbin.env$mrbinplot$plotRegion[4]-(mrbin.env$mrbinplot$plotRegion[3]-mrbin.env$mrbinplot$plotRegion[4])*.1)
       plotNMR()
   }
 }
}

#' A function for changing plotNMR plots.
#'
#' This function moves up the plot region of the current NMR plot (only 2D).
#' @return {None}
#' @export
#' @examples
#' mrbin(silent=TRUE,parameters=list(dimension="2D",binwidth2D=0.5,
#'          binheight=20,PQNScaling="No",referenceScaling="No",binRegion=c(8,2,20,140),
#'          fixNegatives="No",logTrafo="No",PCA="No",verbose=FALSE,saveFiles="No",
#'          NMRfolders=c(system.file("extdata/1/12/pdata/10",package="mrbin"))))
#' plotNMR()
#' zoomIn()
#' up()

up<-function(){
 if(!is.null(mrbin.env$mrbinTMP$currentSpectrum)){
   if(mrbin.env$mrbinparam$dimension=="2D"){
       topMax<-min(as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrum)))
       bottomMax<-max(as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrum)))
       mrbin.env$mrbinplot$plotRegion[3]<-max(topMax,
                                 mrbin.env$mrbinplot$plotRegion[3]+(mrbin.env$mrbinplot$plotRegion[3]-mrbin.env$mrbinplot$plotRegion[4])*.1)
       mrbin.env$mrbinplot$plotRegion[4]<-min(bottomMax,
                                 mrbin.env$mrbinplot$plotRegion[4]+(mrbin.env$mrbinplot$plotRegion[3]-mrbin.env$mrbinplot$plotRegion[4])*.1)
       plotNMR()
   }
 }
}

#' A function for saving the package environment.
#'
#' This function returns a list of all objects of the current package environment.
#' This may be helpful for debugging or for accessing NMR spectral data and the raw bin data.
#' @return A list containing all objects from the local package environment.
#' @export
#' @examples
#' tempList<-getEnv()

getEnv<-function(){
  mrbin.envCopy <- mget(ls(envir = mrbin.env), mrbin.env)
  return(mrbin.envCopy)
}

#' A function for changing and adding variables in the package environment.
#'
#' This function can change variables in the current package environment.
#' This may be helpful for debugging or for some plotting functions.
#' @param variableList A list containing all objects to be saved in the local package environment.
#' @return {None}
#' @export
#' @examples
#' \donttest{ putToEnv(list(bins=NULL)) }

putToEnv<-function(variableList){
  if(!is.list(variableList)){
    warning("Parameter is not a list. No action performed.")
  } else {
   for(i in 1:length(variableList)){
      assign(names(variableList)[i],variableList[[i]],envir =mrbin.env)
   }
  }
}
