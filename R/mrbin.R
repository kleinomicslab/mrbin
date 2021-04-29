#mrbin - Collection of R functions for analyzing NMR metabolomics data.
#Written by Matthias Klein, The Ohio State University
#
#Package: mrbin
#Title: Magnetic Resonance Binning, Integration and Normalization
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
#'                    parameters=list(verbose=TRUE,dimension="1D",PQNScaling="No",
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
         noiseLevels <- apply(mrbin.env$mrbinTMP$noise_level,1,median)
       }
     }
     if(is.null(noiseLevels)){
         noiseTMP<-sort(NMRdata[NMRdata>0])[ceiling(.01*length(NMRdata[
                        NMRdata>0]))]
     }
     for(i in 1:ncol(NMRdata)){
        negatives<-NMRdata[,i]<=0
        if(sum(negatives)>0){
          if(!is.null(noiseLevels)){#If noise levels are available, restrict range to below noise
             #if(sum(negatives)>0){
               noiseTMP<-stats::median(noiseLevels[negatives])
             #} else {
             #  noiseTMP<-stats::median(noiseLevels)
             #}
          }
            minTMP<-min(NMRdata[negatives,i])#select lowest bin
            if(sum(!negatives)>0){
              maxTMP<-min(noiseTMP,min(NMRdata[!negatives,i]))#select lowest bin above 0
            } else {
              maxTMP<-noiseTMP
            }
            NMRdata[negatives,i]<-(NMRdata[negatives,i]+(maxTMP-minTMP))/
                                          (maxTMP-minTMP)*(maxTMP*.99)+
                                           maxTMP*.01
        }
     }
     if("bins"%in%ls(envir=mrbin.env)&"mrbinparam"%in%ls(envir=mrbin.env)){
       if(!is.null(mrbin.env$bins)){
          if(nrow(mrbin.env$bins)==1){
            mrbin.env$bins<-matrix(NMRdata,nrow=1)
            rownames(mrbin.env$bins)<-rownames(NMRdata)
            colnames(mrbin.env$bins)<-colnames(NMRdata)
          } else {
            mrbin.env$bins<-NMRdata
          }
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
#' @keywords internal
#' @noRd
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
               mrbinversion=as.character(utils::packageVersion("mrbin")),
               binsRaw=NULL,
               medians=NULL,
               noise_level_TMP=NULL,
               noise_level=NULL,
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
               currentSpectrumOriginal=NULL,
               currentSpectrumName=NULL,
               currentSpectrumFolderName=NULL,
               currentSpectrumEXPNO=NULL,
               currentSpectrumFolderName_EXPNO=NULL,
               currentSpectrumTitle=NULL,
               i=1
    ),mrbin.env)
    assign("requiredParam",c(
               "dimension","binMethod","binRegion","specialBinList","referenceScaling",
               "removeSolvent","removeAreas","sumBins","noiseRemoval","PQNScaling",
               "PQNIgnoreSugarArea","PQNsugarArea","fixNegatives","logTrafo","tryParallel",
               "defineGroups","PCA","solventRegion","noiseThreshold","trimZeros",
               "PQNminimumFeatures","PCAtitlelength","createBins","useAsNames","saveFiles",
               "outputFileName","verbose","removeAreaList","sumBinList","Factors","NMRfolders"
               ),mrbin.env)
    assign("requiredParam1D",c(
               "binwidth1D","reference1D","signal_to_noise1D","noiseRange1d",
               mrbin.env$requiredParam
               ),mrbin.env)
    assign("requiredParam2D",c(
               "binwidth2D","binheight","cropHSQC","reference2D","signal_to_noise2D",
               "noiseRange2d","croptopRight","croptopLeft","cropbottomRight","cropbottomLeft",
               mrbin.env$requiredParam
               ),mrbin.env)
    assign("mrbinparam", list(
               dimension="1D",
               binMethod="Rectangular bins",#"Custom bin list"
               binwidth1D=.01,
               binwidth2D=.02,
               binheight=1,
               binRegion=c(9.5,0.5,10,156),
               specialBinList=NULL,
               referenceScaling="Yes",
               removeSolvent="Yes",
               removeAreas="No",
               sumBins="No",
               trimZeros="Yes",
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
               removeAreaList=NULL,
               sumBinList=NULL,
               showSpectrumPreview="Yes",
               noiseThreshold=0.75,
               signal_to_noise1D=25,
               signal_to_noise2D=6,
               noiseRange2d=c(3.3,2.3,90,110),
               noiseRange1d=c(10,9.5),
               croptopRight=c(0,-1.50),#only 2D, defines edge points of the cropped area
               croptopLeft=c(0,3.5),
               cropbottomRight=c(160,6),
               cropbottomLeft=c(160,10),
               PQNIgnoreSugarArea="Yes",#exclude most glucose to reduce impact
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
               noise_level_Raw=NULL,
               noise_level_adjusted=NULL,
               numberOfFeaturesRaw=NULL,
               numberOfFeaturesAfterRemovingSolvent=NULL,
               numberOfFeaturesAfterRemovingAreas=NULL,
               numberOfFeaturesAfterSummingBins=NULL,
               numberOfFeaturesAfterTrimmingZeros=NULL,
               numberOfFeaturesAfterNoiseRemoval=NULL,
               numberOfFeaturesAfterCropping=NULL,
               tryParallel=TRUE,
               warningMessages=NULL
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
#' and performs the chosen data processing steps.
#' If a list of parameters is provided and silent is set to TRUE, no user input
#' is requested and binning and data processing are performed silently.
#' @param parameters Optional: A list of parameters, see examples for details. If omitted, the user will be asked through a series of question to set the parameters.
#' @param silent If TRUE, the user will be asked no questions and binning and data analysis will run according to the current parameters. Defaults to FALSE.
#' @param setDefault If TRUE, all current parameters will be replaced by the default parameters (before loading any provided parameters sets). Defaults to FALSE.
#' @return An invisible list containing bins (data after processing), parameters, and factors
#' @export
#' @examples
#' # Set parameters in command line.
#' mrbinExample<-mrbin(silent=TRUE,setDefault=TRUE,parameters=list(dimension="1D",
#'                 binwidth1D=0.05,signal_to_noise1D=100,
#'                 NMRfolders=c(system.file("extdata/1/10/pdata/10",package="mrbin"),
#'                             system.file("extdata/2/10/pdata/10",package="mrbin"),
#'                             system.file("extdata/3/10/pdata/10",package="mrbin")),
#'                 Factors=factor(c("Group A","Group A","Group B"))))

mrbin<-function(silent=FALSE,setDefault=FALSE,parameters=NULL){
  if(!exists("mrbin.env", mode="environment")) .onLoad()
  if(setDefault) resetEnv()
  if(!is.null(parameters)){
      setParam(parameters)
  }
  #if(mrbin.env$mrbinparam$verbose){
  #     message(paste("\nmrbin version ",mrbin.env$mrbinTMP$mrbinversion,"\n",sep=""), appendLF = FALSE)
  #}
  stopTMP<-FALSE
  selectionRepeat<-""
  if(silent) startmrbin<-"Start binning now"
  #Create bin list?
  if(!silent){
   #On Apple or Mac computer, display a hint for installing Quartz
   if(Sys.info()['sysname']=='Darwin'){
     if(mrbin.env$mrbinparam$verbose){
       message("Hint: If you see text lists instead of dialog boxes, please install xquartz \nfrom https://www.xquartz.org")
       utils::flush.console()
     }
   }
   selectStep<--2
   lastStepDone<-FALSE
   while(!lastStepDone&!stopTMP){
     if(selectStep==-2){
       selectTMPNo<-"Do not show verbose hints and results"
       selectTMPYes<-"Show verbose hints and results (recommended)"
       preselectTMP<-selectTMPNo
       if(mrbin.env$mrbinparam$verbose) preselectTMP<-selectTMPYes
       selection<-utils::select.list(c(selectTMPYes,selectTMPNo),preselect=preselectTMP,
                  title="Show verbose hints?",graphics=TRUE)
       if(length(selection)==0|selection=="") stopTMP<-TRUE
       if(!stopTMP){
        if(selection==selectTMPYes) mrbin.env$mrbinparam$verbose<-TRUE
        if(selection==selectTMPNo) mrbin.env$mrbinparam$verbose<-FALSE
       }
       if(!stopTMP) selectStep<-selectStep+1
     }
     if(selectStep==-1){
       selectTMPNo<-"Do not use parallel computing"
       selectTMPYes<-"Use parallel package for speed"
       preselectTMP<-selectTMPNo
       if(mrbin.env$mrbinparam$tryParallel) preselectTMP<-selectTMPYes
       selection<-utils::select.list(c(selectTMPYes,selectTMPNo),preselect=preselectTMP,
                  title="Try parallel computing for speed?",graphics=TRUE)
       if(length(selection)==0|selection=="") stopTMP<-TRUE
       if(!stopTMP){
        if(selection==selectTMPYes) mrbin.env$mrbinparam$tryParallel<-TRUE
        if(selection==selectTMPNo) mrbin.env$mrbinparam$tryParallel<-FALSE
       }
       if(!stopTMP) selectStep<-selectStep+1
     }
     if(selectStep==0){
       selectTMPNo<-"Do not show previews (e.g. for slow hardware)"
       selectTMPYes<-"Show spectrum previews (recommended)"
       preselectTMP<-selectTMPNo
       if(mrbin.env$mrbinparam$showSpectrumPreview=="Yes") preselectTMP<-selectTMPYes
       selection<-utils::select.list(c(selectTMPYes,selectTMPNo,"Go back"),preselect=preselectTMP,
                  title="Show spectrum previews?",graphics=TRUE)
       if(length(selection)==0|selection=="") stopTMP<-TRUE
       if(!stopTMP){
          if(selection==selectTMPYes) mrbin.env$mrbinparam$showSpectrumPreview<-"Yes"
          if(selection==selectTMPNo) mrbin.env$mrbinparam$showSpectrumPreview<-"No"
       }
       if(!stopTMP&selection=="Go back") selectStep<-selectStep-2
       if(!stopTMP) selectStep<-selectStep+1
     }
       #selectStep<-1
       #lastStepDone<-FALSE
       #while(!lastStepDone&!stopTMP){
         if(selectStep==1){
           #Select new parameters?
           selectionNewTMP<-c("Edit parameters","Reload from file")
           if(!is.null(mrbin.env$mrbinparam$NMRfolders)) selectionNewTMP<-c(selectionNewTMP,"Use current parameters")
           selectionRepeat<-utils::select.list(c(selectionNewTMP,"Go back"),preselect="Edit parameters",
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
           if(!stopTMP&selectionRepeat=="Go back") selectStep<-selectStep-2
           if(!stopTMP) selectStep<-selectStep+1
         }
         if(selectionRepeat=="Use current parameters"&!stopTMP){
           selectStep<-19
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
               if(dimension%in%c("1D","2D")){
                 mrbin.env$mrbinparam$dimension<-dimension
                 if(dimension=="1D"){
                   dimlength<-2
                 }
                 if(dimension=="2D"){
                   dimlength<-4
                 }

               }
             }
           if(!stopTMP&dimension=="Go back") selectStep<-selectStep-2
           if(!stopTMP) selectStep<-selectStep+1
           }
           if(selectStep==3){
             if(!stopTMP){
               #Select folders
               addFoldersTMP<-""
               if(!is.null(mrbin.env$mrbinparam$NMRfolders)){
                    selectionFoldersYes<-paste("Keep previous list (",length(mrbin.env$mrbinparam$NMRfolders),
                                         " spectra)",sep="")
                    selectionFolders<-utils::select.list(c("Create new spectra list",selectionFoldersYes,
                                      "Add or remove spectra from previous list","Go back"),
                                      preselect=selectionFoldersYes,
                                      title="Use previous spectra list?",graphics=TRUE)
                    if(length(selectionFolders)==0|selectionFolders=="") stopTMP<-TRUE
                    if(!stopTMP){
                      if(selectionFolders=="Create new spectra list"){
                        selectionFolders<-selectFolders()
                        if(selectionFolders=="stop")  stopTMP<-TRUE
                      }
                    }
                    if(!stopTMP){
                      if(selectionFolders=="Add or remove spectra from previous list"){
                        removeSpectrum()
                        addFoldersTMP<-utils::select.list(c("Add spectra to list","Keep list",
                                      "Go back"),
                                      preselect="Keep list",
                                      title="Add spectra to list?",graphics=TRUE)
                        if(length(addFoldersTMP)==0|addFoldersTMP=="") stopTMP<-TRUE
                        if(!stopTMP){
                          if(addFoldersTMP=="Add spectra to list"){
                            selectionFolders<-selectFolders(keep=TRUE)
                            if(selectionFolders=="stop")  stopTMP<-TRUE
                          }
                        }
                      }
                    }
               } else {
                    selectionFolders<-selectFolders()
                    if(selectionFolders=="stop")  stopTMP<-TRUE
               }
             }
             if(!stopTMP&(selectionFolders=="Go back"|addFoldersTMP=="Go back")) selectStep<-selectStep-2
             if(!stopTMP) selectStep<-selectStep+1
           }
           if(selectStep==4){
              if(!stopTMP){
                  mrbin.env$mrbinTMP$currentFolder<-mrbin.env$mrbinparam$NMRfolders[1]
                  if(mrbin.env$mrbinparam$showSpectrumPreview=="Yes"|mrbin.env$mrbinparam$PCA=="Yes"){
                    readNMR2()
                  } else {
                    readNMR2(onlyTitles=TRUE)
                  }
              }
              if(!stopTMP){
                #Use rectangluar bins or use special bin list, e.g. for lipids
                binMethodpreSelect<-mrbin.env$mrbinparam$binMethod
                userDefBinTMP<-"User defined bin list, e.g. for lipid analysis"
                if(binMethodpreSelect=="Custom bin list") binMethodpreSelect<-userDefBinTMP
                binMethod<-utils::select.list(c("Rectangular bins",userDefBinTMP,"Go back"),
                           preselect=binMethodpreSelect,
                           ,title ="Binning method: ",graphics=TRUE)
                if(length(binMethod)==0|binMethod=="") stopTMP<-TRUE
                if(!stopTMP){
                  if(!binMethod=="Go back"){
                    if(binMethod==userDefBinTMP) binMethod<-"Custom bin list"
                    if(!identical(mrbin.env$mrbinparam$binMethod,binMethod)) mrbin.env$paramChangeFlag<-TRUE
                    mrbin.env$mrbinparam$binMethod<-binMethod
                    #Bin region
                    adjRegion<-""
                    if(mrbin.env$mrbinparam$binMethod=="Rectangular bins"){
                      if(mrbin.env$mrbinparam$verbose){
                        message("Hint: Include all visible peaks, excluding reference")
                        utils::flush.console()
                      }
                      accept<-FALSE
                      while(!accept&!stopTMP){
                        binRegionText<-paste(paste(c("left=","ppm, right=","ppm, top=","ppm, bottom=")[1:dimlength],
                                          mrbin.env$mrbinparam$binRegion[1:dimlength],collapse="",sep=""),"ppm",sep="")
                        if(mrbin.env$mrbinparam$showSpectrumPreview=="Yes") plotNMR(region="all",rectangleRegions=matrix(mrbin.env$mrbinparam$binRegion,ncol=4),color="black",
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
                  if(mrbin.env$mrbinparam$verbose){
                    message("Hint: Should exceed peak size, and include at least 3 data points")
                    utils::flush.console()
                  }
                  accept<-FALSE
                  widthAdjust<-""
                  while(!accept&!stopTMP&!widthAdjust=="Go back"){
                    if(mrbin.env$mrbinparam$showSpectrumPreview=="Yes") plotNMR(region=c(1.5,1.15,16.5,27),
                            rectangleRegions=matrix(c(1.35+as.numeric(mrbin.env$mrbinparam$binwidth1D),
                                                    1.35,21,21+1),ncol=4),
                            color="black", showGrid=TRUE,
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
                  if(mrbin.env$mrbinparam$verbose){
                    message("Hint: Should exceed peak size, and include at least 3x3 data points")
                    utils::flush.console()
                  }
                  accept<-FALSE
                  widthAdjust<-""
                  #heightAdjust<-""
                  while(!accept&!stopTMP&!widthAdjust=="Go back"){
                    if(mrbin.env$mrbinparam$showSpectrumPreview=="Yes") plotNMR(region=c(2,1,16.5,27),
                            rectangleRegions=matrix(c(1.5+as.numeric(mrbin.env$mrbinparam$binwidth2D),
                                                    1.5,21,21+mrbin.env$mrbinparam$binheight),ncol=4),
                            color="black",manualScale=FALSE, showGrid=TRUE,
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
                if(!is.null(mrbin.env$mrbinparam$specialBinList)){
                  if(nrow(mrbin.env$mrbinparam$specialBinList)==0) mrbin.env$mrbinparam$specialBinList<-NULL
                }
                adjbinRegionSelect<-""
                adjbinRegionAccept<-""
                if(!is.null(mrbin.env$mrbinparam$specialBinList)){
                  if(nrow(mrbin.env$mrbinparam$specialBinList)==1){
                    specialBinList_s<-""
                    specialBinList_dots<-""
                  } else {
                    specialBinList_s<-"s"
                    specialBinList_dots<-", ..."
                    }
                  keepbinRegionText<-paste(paste(c("left=","ppm, right=","ppm, top=","ppm, bottom=")[1:dimlength],
                              mrbin.env$mrbinparam$specialBinList[1:dimlength],collapse="",sep=""),"ppm",sep="")
                  keepbinRegionYes<-paste("Keep previous bin list (",nrow(mrbin.env$mrbinparam$specialBinList),
                                               " bin",specialBinList_s,", ",keepbinRegionText,specialBinList_dots,")",sep="")
                  editbinRegionYes<-"Edit previous bin list"
                  preselectbinRegion<-keepbinRegionYes
                  keepbinRegionIndex<-c(1,2,3,4)
                } else {
                  preselectbinRegion<-"Create new bin list"
                  keepbinRegionYes<-"Keep previous bin list"
                  editbinRegionYes<-"Edit previous bin list"
                  keepbinRegionIndex<-c(1,4)
                }
                if(mrbin.env$mrbinparam$showSpectrumPreview=="Yes") plotNMR(region="all",rectangleRegions=mrbin.env$mrbinparam$specialBinList,color="black",
                        manualScale=FALSE,rectangleColors="green",
                        plotTitle=paste("Bin regions\n",
                        sep=""))
                adjbinRegionSelect<-utils::select.list(c("Create new bin list",keepbinRegionYes,
                                  editbinRegionYes,"Go back")[keepbinRegionIndex],
                                  preselect=preselectbinRegion,
                                  title="Create new bin list?",graphics=TRUE)
                if(length(adjbinRegionSelect)==0|adjbinRegionSelect=="") stopTMP<-TRUE
                if(!stopTMP){
                  if(adjbinRegionSelect=="Create new bin list"){
                    mrbin.env$mrbinparam$specialBinList<-NULL
                  }
                }
                #if(!stopTMP&!is.null(mrbin.env$mrbinparam$specialBinList)){
                if(!stopTMP){
                  if(adjbinRegionSelect=="Create new bin list"|adjbinRegionSelect==editbinRegionYes){
                    if(mrbin.env$mrbinparam$verbose){
                      message("Hint: Should exceed peak size, and include at least 3 data points")
                      utils::flush.console()
                    }
                    #for(ibinRegions in 1:nrow(mrbin.env$mrbinparam$specialBinList)){
                    if(is.null(mrbin.env$mrbinparam$specialBinList)){
                            mrbin.env$mrbinparam$specialBinList<-matrix(ncol=4,nrow=0,dimnames=list(NULL,c("left","right","top","bottom")))
                    }
                    ibinRegions <- 1
                    adjbinRegion<-""
                    addbinRegion<-""
                    while(ibinRegions <= (nrow(mrbin.env$mrbinparam$specialBinList)+1)&!stopTMP&
                          !adjbinRegion=="Go back"&!addbinRegion=="Go back"&!addbinRegion=="No"){
                      if(!stopTMP&!adjbinRegion=="Go back"){
                        if(ibinRegions>nrow(mrbin.env$mrbinparam$specialBinList)){
                          addbinRegion<-utils::select.list(c("Yes","No","Go back"),preselect="Yes",
                                     title ="Add a new bin?",graphics=TRUE)
                          if(length(addbinRegion)==0|addbinRegion==""){
                            stopTMP<-TRUE
                            addbinRegion<-""
                          }
                          if(!stopTMP){
                            if(addbinRegion=="Yes"){
                              mrbin.env$paramChangeFlag<-TRUE
                              mrbin.env$mrbinparam$specialBinList<-rbind(mrbin.env$mrbinparam$specialBinList,c(0,0,0,0))
                              if(nrow(mrbin.env$mrbinparam$specialBinList)==1){
                                rownames(mrbin.env$mrbinparam$specialBinList)<-""
                              } else {
                                rownames(mrbin.env$mrbinparam$specialBinList)[ibinRegions]<-""
                              }

                            }
                          }
                        }
                        if(!stopTMP&!adjbinRegion=="Go back"&!addbinRegion=="Go back"&!addbinRegion=="No"){
                          mean1<-mean(mrbin.env$mrbinparam$specialBinList[ibinRegions,1:2])
                          range1<-mrbin.env$mrbinparam$specialBinList[ibinRegions,1]-mrbin.env$mrbinparam$specialBinList[ibinRegions,2]
                          mean2<-mean(mrbin.env$mrbinparam$specialBinList[ibinRegions,3:4])
                          range2<-mrbin.env$mrbinparam$specialBinList[ibinRegions,4]-mrbin.env$mrbinparam$specialBinList[ibinRegions,3]
                          regionTMP<-c(mean1+3.5*range1,mean1-3.5*range1,mean2-3.5*range2,mean2+3.5*range2)
                          showGridTMP<-TRUE
                          if(sum(mrbin.env$mrbinparam$specialBinList[ibinRegions,]==0)==4){
                            regionTMP<-"all"
                            showGridTMP<-FALSE
                          }
                          if(mrbin.env$mrbinparam$dimension=="1D"){
                            if(!range1==0){
                              if(mrbin.env$mrbinparam$showSpectrumPreview=="Yes") plotNMR(region=regionTMP,
                                      rectangleRegions=matrix(c(mrbin.env$mrbinparam$specialBinList[ibinRegions,1],
                                                              mrbin.env$mrbinparam$specialBinList[ibinRegions,2],0,2),ncol=4),
                                      color="black", showGrid=showGridTMP,manualScale=FALSE,
                                      plotTitle=paste("Custom bins\nleft=",mrbin.env$mrbinparam$specialBinList[ibinRegions,1],
                                                "ppm, right=",mrbin.env$mrbinparam$specialBinList[ibinRegions,2],"ppm",sep=""))
                            }
                          }
                          if(mrbin.env$mrbinparam$dimension=="2D"){
                            if(!range1==0&!range2==0){
                              if(mrbin.env$mrbinparam$showSpectrumPreview=="Yes") plotNMR(region=regionTMP,
                                      rectangleRegions=matrix(mrbin.env$mrbinparam$specialBinList[ibinRegions,1:4],ncol=4),
                                      color="black",manualScale=FALSE, showGrid=showGridTMP,
                                      plotTitle=paste("Custom bins\nleft=",mrbin.env$mrbinparam$specialBinList[ibinRegions,1],
                                                "ppm, right=",mrbin.env$mrbinparam$specialBinList[ibinRegions,2],"ppm",sep=""))
                            }
                          }
                          adjbinRegionAccept<-paste(paste(c("Keep left=","ppm, right=","ppm, top=","ppm, bottom=")[1:dimlength],
                                            mrbin.env$mrbinparam$specialBinList[ibinRegions,1:dimlength],collapse="",sep=""),"ppm",sep="")
                          if(sum(mrbin.env$mrbinparam$specialBinList[ibinRegions,]==0)==4){
                            adjbinRegion<-"Change..."
                          } else {
                            if(rownames(mrbin.env$mrbinparam$specialBinList)[ibinRegions]==""){
                              binTitleTMP<-""
                            } else {
                              binTitleTMP<-paste(" (\"",rownames(mrbin.env$mrbinparam$specialBinList)[ibinRegions],"\")",sep="")
                            }
                            adjbinRegion<-utils::select.list(c(adjbinRegionAccept,"Change...","Remove bin","Go back"),
                                     preselect=adjbinRegionAccept,
                                     title =paste("Edit bin ",ibinRegions,binTitleTMP,"?",sep=""),graphics=TRUE)
                          }
                          if(length(adjbinRegion)==0|adjbinRegion=="") stopTMP<-TRUE
                        }
                      }
                      if(adjbinRegion=="Change..."&!stopTMP&!adjbinRegion=="Go back"){

                        if(rownames(mrbin.env$mrbinparam$specialBinList)[ibinRegions]==""){
                          promptTMP<-paste("New bin name, press enter for no name: ",sep="")
                        } else {
                          promptTMP<-paste("New bin name, press enter to keep ",
                                    rownames(mrbin.env$mrbinparam$specialBinList)[ibinRegions],": ",sep="")
                        }
                        nameTMP<-readline(prompt=promptTMP)
                        if(!nameTMP=="") {
                               mrbin.env$paramChangeFlag<-TRUE
                               rownames(mrbin.env$mrbinparam$specialBinList)[ibinRegions]<-nameTMP
                        }
                        regionTMP<-readline(prompt=paste("Bin ",ibinRegions,": left border, press enter to keep ",
                                  mrbin.env$mrbinparam$specialBinList[ibinRegions,1],": ",sep=""))
                        if(!regionTMP=="") {
                               mrbin.env$paramChangeFlag<-TRUE
                               mrbin.env$mrbinparam$specialBinList[ibinRegions,1]<-as.numeric(regionTMP)
                        }
                        regionTMP<-readline(prompt=paste("Bin ",ibinRegions,": right border, press enter to keep ",
                                  mrbin.env$mrbinparam$specialBinList[ibinRegions,2],": ",sep=""))
                        if(!regionTMP=="") {
                               mrbin.env$paramChangeFlag<-TRUE
                               mrbin.env$mrbinparam$specialBinList[ibinRegions,2]<-as.numeric(regionTMP)
                        }
                        if(mrbin.env$mrbinparam$specialBinList[ibinRegions,1]<mrbin.env$mrbinparam$specialBinList[ibinRegions,2]){
                          TMP<-mrbin.env$mrbinparam$specialBinList[ibinRegions,1]
                          mrbin.env$mrbinparam$specialBinList[ibinRegions,1]<-mrbin.env$mrbinparam$specialBinList[ibinRegions,2]
                          mrbin.env$mrbinparam$specialBinList[ibinRegions,2]<-TMP
                        }
                        if(mrbin.env$mrbinparam$dimension=="2D"&!stopTMP){
                          regionTMP<-readline(prompt=paste("Bin ",ibinRegions,": top border, press enter to keep ",
                                    mrbin.env$mrbinparam$specialBinList[ibinRegions,3],": ",sep=""))
                          if(!regionTMP=="") {
                                 mrbin.env$paramChangeFlag<-TRUE
                                 mrbin.env$mrbinparam$specialBinList[ibinRegions,3]<-as.numeric(regionTMP)
                          }
                          regionTMP<-readline(prompt=paste("Bin ",ibinRegions,": bottom border, press enter to keep ",
                                    mrbin.env$mrbinparam$specialBinList[ibinRegions,4],": ",sep=""))
                          if(!regionTMP=="") {
                                 mrbin.env$paramChangeFlag<-TRUE
                                 mrbin.env$mrbinparam$specialBinList[ibinRegions,4]<-as.numeric(regionTMP)
                          }
                          if(mrbin.env$mrbinparam$specialBinList[ibinRegions,4]<mrbin.env$mrbinparam$specialBinList[ibinRegions,3]){
                            TMP<-mrbin.env$mrbinparam$specialBinList[ibinRegions,3]
                            mrbin.env$mrbinparam$specialBinList[ibinRegions,3]<-mrbin.env$mrbinparam$specialBinList[ibinRegions,4]
                            mrbin.env$mrbinparam$specialBinList[ibinRegions,4]<-TMP
                          }
                        }
                      }
                      if(adjbinRegion=="Remove bin"&!stopTMP&!adjbinRegion=="Go back"){
                        mrbin.env$mrbinparam$specialBinList<-mrbin.env$mrbinparam$specialBinList[-ibinRegions,,drop=FALSE]
                      }
                      if(adjbinRegion==adjbinRegionAccept&!stopTMP&!adjbinRegion=="Go back"){
                        ibinRegions <- ibinRegions+1
                      }
                    }
                    if(nrow(mrbin.env$mrbinparam$specialBinList)==0){
                      mrbin.env$mrbinparam$specialBinList<-NULL
                      adjbinRegion<-"Go back"
                    }
                  }
                }
                if(!stopTMP&!addbinRegion=="Go back"&!adjbinRegion=="Go back"){
                  if(!is.null(mrbin.env$mrbinparam$specialBinList)){
                    #if(nrow(mrbin.env$mrbinparam$specialBinList)>0){
                    #  mrbin.env$mrbinTMP$binRegions<-mrbin.env$mrbinparam$specialBinList
                    #  if(!is.matrix(mrbin.env$mrbinTMP$binRegions)) mrbin.env$mrbinTMP$binRegions<-matrix(mrbin.env$mrbinTMP$binRegions,ncol=4)
                    #}
                  } else {
                    addbinRegion<-"Go back"
                  }
                }
                if(adjbinRegion=="Go back"|addbinRegion=="Go back"|adjbinRegionSelect=="Go back"){
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
                    if(mrbin.env$mrbinparam$verbose){
                      message("Hint: Include reference signal with some margin")
                      utils::flush.console()
                    }
                    if(mrbin.env$mrbinparam$dimension=="1D"){
                    accept<-FALSE
                    adjRegion<-""
                    while(!accept&!stopTMP&!adjRegion=="Go back"){
                      mean1<-mean(mrbin.env$mrbinparam$reference1D)
                      range1<-max(mrbin.env$mrbinparam$reference1D)-min(mrbin.env$mrbinparam$reference1D)
                      if(mrbin.env$mrbinparam$showSpectrumPreview=="Yes") plotNMR(region=c(mean1+4*range1,mean1-4*range1,-10,10),
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
                        if(mrbin.env$mrbinparam$showSpectrumPreview=="Yes") plotNMR(region=c(mean1+4*range1,mean1-4*range1,
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
                    if(mrbin.env$mrbinparam$verbose){
                      message("Hint: Include solvent signal and some margin")
                      utils::flush.console()
                    }
                    accept<-FALSE
                    adjRegion<-""
                    while(!accept&!stopTMP&!adjRegion=="Go back"){
                      mean1<-mean(mrbin.env$mrbinparam$solventRegion[1:2])
                      range1<-max(mrbin.env$mrbinparam$solventRegion[1:2])-min(mrbin.env$mrbinparam$solventRegion[1:2])
                      if(mrbin.env$mrbinparam$showSpectrumPreview=="Yes") plotNMR(region=c(mean1+6*range1,mean1-6*range1,-10,160),
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
                adjbinRegion<-""
                addbinRegion<-""
                if(mrbin.env$mrbinparam$verbose){
                  message("Hint: Remove solvents, signals from blank samples, spectral artifacts")
                  utils::flush.console()
                }
                removeAreas<-utils::select.list(c("Yes","No","Go back"),preselect=mrbin.env$mrbinparam$removeAreas,
                                         title = "Remove additional areas?",graphics=TRUE)
                if(length(removeAreas)==0|removeAreas=="") stopTMP<-TRUE
                if(!stopTMP&!removeAreas=="Go back"){
                  mrbin.env$mrbinparam$removeAreas<-removeAreas
                  if(mrbin.env$mrbinparam$removeAreas=="Yes"){
                    addAreasFlag<-TRUE
                    if(!is.null( mrbin.env$mrbinparam$removeAreaList)){
                      if(nrow(mrbin.env$mrbinparam$removeAreaList)==0){
                        mrbin.env$mrbinparam$removeAreaList<-NULL
                      }
                    }
                    if(!is.null(mrbin.env$mrbinparam$removeAreaList)){
                      if(nrow(mrbin.env$mrbinparam$removeAreaList)>0){
                          addAreasFlag<-FALSE
                          if(nrow(mrbin.env$mrbinparam$removeAreaList)==1){
                            regions_s<-""
                            regions_dots<-""
                          }
                          if(nrow(mrbin.env$mrbinparam$removeAreaList)>1){
                            regions_s<-"s"
                            regions_dots<-", ..."
                          }
                          preselectKeepTMP<-paste("Keep current list (",nrow(mrbin.env$mrbinparam$removeAreaList)," region",regions_s,", ",
                                            paste(c("left=","ppm, right=","ppm, top=","ppm ,bottom=")[1:dimlength],
                                              mrbin.env$mrbinparam$removeAreaList[1,1:dimlength],
                                              sep="",collapse=""),
                                            "ppm",regions_dots,")",sep="")
                          preselectKeepTMPYes<-preselectKeepTMP
                          keepbinRegionIndex<-c(1,2,3,4)
                      }
                    } else {
                       preselectKeepTMP<-"Keep current list"
                       preselectKeepTMPYes<-"Create new list"
                       keepbinRegionIndex<-c(1,4)
                    }
                    if(mrbin.env$mrbinparam$showSpectrumPreview=="Yes") plotNMR(region="all",rectangleRegions=mrbin.env$mrbinparam$removeAreaList,color="black",
                            manualScale=FALSE,rectangleColors="green",
                            plotTitle=paste("Removed areas\n",
                            sep=""))
                    removeAreaListTMP<-utils::select.list(c("Create new list",preselectKeepTMP,"Edit current list","Go back")[keepbinRegionIndex],
                                       preselect=preselectKeepTMPYes,
                                       title = "Use previous area list or define new?",graphics=TRUE)
                    if(length(removeAreaListTMP)==0|removeAreaListTMP=="") stopTMP<-TRUE
                    if(!removeAreaListTMP==preselectKeepTMP&!stopTMP&!removeAreaListTMP=="Go back"){
                      addAreasFlag<-TRUE
                      if(removeAreaListTMP=="Create new list"&!stopTMP){
                        mrbin.env$mrbinparam$removeAreaList<-matrix(ncol=4,nrow=0,dimnames=list(NULL,c("left","right","top","bottom")))
                      }
                    }
                    #if(!stopTMP){
                    #  iaddAreas<-nrow(mrbin.env$mrbinparam$removeAreaList)+1
                    #}
                    if(!stopTMP){
                      if(removeAreaListTMP=="Create new list"|removeAreaListTMP=="Edit current list"){
                        #for(ibinRegions in 1:nrow(mrbin.env$mrbinparam$specialBinList)){
                        if(is.null(mrbin.env$mrbinparam$removeAreaList)){
                                mrbin.env$mrbinparam$removeAreaList<-matrix(ncol=4,nrow=0,dimnames=list(NULL,c("left","right","top","bottom")))
                        }
                        ibinRegions <- 1
                        adjbinRegion<-""
                        addbinRegion<-""
                        adjbinRegionAccept<-""
                        #while(addAreasFlag&!stopTMP){
                        while(ibinRegions <= (nrow(mrbin.env$mrbinparam$removeAreaList)+1)&!stopTMP&
                              !adjbinRegion=="Go back"&!addbinRegion=="Go back"&!addbinRegion=="No"){
                          if(!stopTMP&!adjbinRegion=="Go back"){
                            if(ibinRegions>nrow(mrbin.env$mrbinparam$removeAreaList)){
                              addbinRegion<-utils::select.list(c("Yes","No","Go back"),preselect="No",
                                         title ="Add a new region?",graphics=TRUE)
                              if(length(addbinRegion)==0|addbinRegion==""){
                                stopTMP<-TRUE
                                addbinRegion<-""
                              }
                              if(!stopTMP){
                                if(addbinRegion=="Yes"){
                                  mrbin.env$paramChangeFlag<-TRUE
                                  mrbin.env$mrbinparam$removeAreaList<-rbind(mrbin.env$mrbinparam$removeAreaList,c(0,0,0,0))
                                }
                              }
                            }
                            if(!stopTMP&!adjbinRegion=="Go back"&!addbinRegion=="Go back"&!addbinRegion=="No"){
                              mean1<-mean(mrbin.env$mrbinparam$removeAreaList[ibinRegions,1:2])
                              range1<-mrbin.env$mrbinparam$removeAreaList[ibinRegions,1]-mrbin.env$mrbinparam$removeAreaList[ibinRegions,2]
                              mean2<-mean(mrbin.env$mrbinparam$removeAreaList[ibinRegions,3:4])
                              range2<-mrbin.env$mrbinparam$removeAreaList[ibinRegions,4]-mrbin.env$mrbinparam$removeAreaList[ibinRegions,3]
                              regionTMP<-c(mean1+4*range1,mean1-4*range1,mean2-4*range2,mean2+4*range2)
                              if(sum(mrbin.env$mrbinparam$removeAreaList[ibinRegions,]==0)==4) regionTMP<-"all"
                              if(mrbin.env$mrbinparam$dimension=="1D"){
                                if(mrbin.env$mrbinparam$showSpectrumPreview=="Yes") plotNMR(region=regionTMP,
                                        rectangleRegions=matrix(c(mrbin.env$mrbinparam$removeAreaList[ibinRegions,1],
                                                                mrbin.env$mrbinparam$removeAreaList[ibinRegions,2],0,2),ncol=4),
                                        color="black",
                                        manualScale=FALSE,
                                        plotTitle=paste("Remove area\nleft=",mrbin.env$mrbinparam$removeAreaList[ibinRegions,1],
                                                  "ppm, right=",mrbin.env$mrbinparam$removeAreaList[ibinRegions,2],"ppm",sep=""))
                              }
                              if(mrbin.env$mrbinparam$dimension=="2D"){
                                if(mrbin.env$mrbinparam$showSpectrumPreview=="Yes") plotNMR(region=regionTMP,
                                        rectangleRegions=matrix(mrbin.env$mrbinparam$removeAreaList[ibinRegions,1:4],ncol=4),
                                        color="black",
                                        manualScale=FALSE,
                                        plotTitle=paste("Remove area\nleft=",mrbin.env$mrbinparam$removeAreaList[ibinRegions,1],
                                                  "ppm, right=",mrbin.env$mrbinparam$removeAreaList[ibinRegions,2],"ppm",sep=""))
                              }
                              adjbinRegionAccept<-paste(paste(c("Keep left=","ppm, right=","ppm, top=","ppm, bottom=")[1:dimlength],
                                                mrbin.env$mrbinparam$removeAreaList[ibinRegions,1:dimlength],collapse="",sep=""),"ppm",sep="")
                              if(sum(mrbin.env$mrbinparam$removeAreaList[ibinRegions,]==0)==4){
                                adjbinRegion<-"Change..."
                              } else {
                                adjbinRegion<-utils::select.list(c(adjbinRegionAccept,"Change...","Remove entry","Go back"),
                                         preselect=adjbinRegionAccept,
                                         title =paste("Keep region ",ibinRegions,"?",sep=""),graphics=TRUE)
                              }
                              if(length(adjbinRegion)==0|adjbinRegion=="") stopTMP<-TRUE
                            }
                          }
                          if(!stopTMP){
                          if(adjbinRegion=="Change..."&!adjbinRegion=="Go back"){
                            regionTMP<-readline(prompt=paste("Region ",ibinRegions,": left border, press enter to keep ",
                                      mrbin.env$mrbinparam$removeAreaList[ibinRegions,1],": ",sep=""))
                            if(!regionTMP=="") {
                                   mrbin.env$paramChangeFlag<-TRUE
                                   mrbin.env$mrbinparam$removeAreaList[ibinRegions,1]<-as.numeric(regionTMP)
                            }
                            regionTMP<-readline(prompt=paste("Region ",ibinRegions,": right border, press enter to keep ",
                                      mrbin.env$mrbinparam$removeAreaList[ibinRegions,2],": ",sep=""))
                            if(!regionTMP=="") {
                                   mrbin.env$paramChangeFlag<-TRUE
                                   mrbin.env$mrbinparam$removeAreaList[ibinRegions,2]<-as.numeric(regionTMP)
                            }
                            if(mrbin.env$mrbinparam$removeAreaList[ibinRegions,1]<mrbin.env$mrbinparam$removeAreaList[ibinRegions,2]){
                              TMP<-mrbin.env$mrbinparam$removeAreaList[ibinRegions,1]
                              mrbin.env$mrbinparam$removeAreaList[ibinRegions,1]<-mrbin.env$mrbinparam$removeAreaList[ibinRegions,2]
                              mrbin.env$mrbinparam$removeAreaList[ibinRegions,2]<-TMP
                            }
                          if(mrbin.env$mrbinparam$dimension=="2D"&!stopTMP){
                            regionTMP<-readline(prompt=paste("Region ",ibinRegions,": top border, press enter to keep ",
                                      mrbin.env$mrbinparam$removeAreaList[ibinRegions,3],": ",sep=""))
                            if(!regionTMP=="") {
                                   mrbin.env$paramChangeFlag<-TRUE
                                   mrbin.env$mrbinparam$removeAreaList[ibinRegions,3]<-as.numeric(regionTMP)
                            }
                            regionTMP<-readline(prompt=paste("Region ",ibinRegions,": bottom border, press enter to keep ",
                                      mrbin.env$mrbinparam$removeAreaList[ibinRegions,4],": ",sep=""))
                            if(!regionTMP=="") {
                                   mrbin.env$paramChangeFlag<-TRUE
                                   mrbin.env$mrbinparam$removeAreaList[ibinRegions,4]<-as.numeric(regionTMP)
                            }
                          if(mrbin.env$mrbinparam$removeAreaList[ibinRegions,4]<mrbin.env$mrbinparam$removeAreaList[ibinRegions,3]){
                            TMP<-mrbin.env$mrbinparam$removeAreaList[ibinRegions,3]
                            mrbin.env$mrbinparam$removeAreaList[ibinRegions,3]<-mrbin.env$mrbinparam$removeAreaList[ibinRegions,4]
                            mrbin.env$mrbinparam$removeAreaList[ibinRegions,4]<-TMP
                          }
                         }
                        }
                        if(adjbinRegion=="Remove entry"&!stopTMP&!adjbinRegion=="Go back"){
                          mrbin.env$mrbinparam$removeAreaList<-mrbin.env$mrbinparam$removeAreaList[-ibinRegions,,drop=FALSE]
                        }
                        }
                        if(!stopTMP){
                          if(adjbinRegion==adjbinRegionAccept){
                            ibinRegions <- ibinRegions+1
                          }
                        }
                      }
                    }
                  }
                }
               }
              }
              if(!is.null(mrbin.env$mrbinparam$removeAreaList)){
                if(nrow(mrbin.env$mrbinparam$removeAreaList)==0){
                  mrbin.env$mrbinparam$removeAreaList<-NULL
                  adjbinRegion<-"Go back"
                }
              }
              if(adjbinRegion=="Go back"|addbinRegion=="Go back"|removeAreaListTMP=="Go back"|removeAreas=="Go back"){
                 selectStep<-selectStep-2
              }
              if(!stopTMP) selectStep<-selectStep+1
            }
            if(selectStep==9){
              #Sum up bins of unstable peaks
              if(!stopTMP){
                if(mrbin.env$mrbinparam$verbose){
                  message("Hint: Signals that differ in chemical shift from sample to sample")
                  utils::flush.console()
                }
                sumBinListTMP<-""
                removeAreaListTMP<-""
                adjbinRegion<-""
                addbinRegion<-""
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
                      if(!is.null( mrbin.env$mrbinparam$sumBinList)){
                        if(nrow(mrbin.env$mrbinparam$sumBinList)==0){
                          mrbin.env$mrbinparam$sumBinList<-NULL
                        }
                      }
                      if(!is.null( mrbin.env$mrbinparam$sumBinList)){
                        if(nrow(mrbin.env$mrbinparam$sumBinList)>0&!stopTMP){
                            addAreasFlag<-FALSE
                            if(nrow(mrbin.env$mrbinparam$sumBinList)==1){
                              regions_s<-""
                              regions_dots<-""
                            }
                            if(nrow(mrbin.env$mrbinparam$sumBinList)>1){
                              regions_s<-"s"
                              regions_dots<-", ..."
                            }
                            preselectKeepTMP<-paste("Keep current list (",nrow(mrbin.env$mrbinparam$sumBinList)," region",regions_s,", ",
                                              paste(c("left=","ppm, right=","ppm, top=","ppm ,bottom=")[1:dimlength],
                                                mrbin.env$mrbinparam$sumBinList[1,1:dimlength],
                                                sep="",collapse=""),
                                              "ppm",regions_dots,")",sep="")
                            preselectKeepTMPYes<-preselectKeepTMP
                            keepbinRegionIndex<-c(1,2,3,4)
                        }
                      } else {
                        preselectKeepTMP<-"Keep current list"
                        preselectKeepTMPYes<-"Create new list"
                        keepbinRegionIndex<-c(1,4)
                      }
                      if(mrbin.env$mrbinparam$showSpectrumPreview=="Yes") plotNMR(region="all",rectangleRegions=mrbin.env$mrbinparam$sumBinList,color="black",
                              manualScale=FALSE,rectangleColors="green",
                              plotTitle=paste("Summed areas\n",
                              sep=""))
                      sumBinListTMP<-utils::select.list(c("Create new list",preselectKeepTMP,"Edit current list","Go back")[keepbinRegionIndex],
                                     preselect=preselectKeepTMPYes,
                                     title = "Use previous area list or define new?",graphics=TRUE)
                      if(length(sumBinListTMP)==0|sumBinListTMP=="") stopTMP<-TRUE
                      if(!sumBinListTMP==preselectKeepTMP&!stopTMP&!sumBinListTMP=="Go back"){
                        addAreasFlag<-TRUE
                        if(sumBinListTMP=="Create new list"&!stopTMP){
                            mrbin.env$mrbinparam$sumBinList<-matrix(ncol=4,nrow=0,dimnames=list(NULL,c("left","right","top","bottom")))
                        }
                      }
                    #}
                    #if(!stopTMP&!sumBinListTMP=="Go back"){
                    #  iaddAreas<-nrow(mrbin.env$mrbinparam$sumBinList)+1
                    #}
                    if(!stopTMP){
                      if(sumBinListTMP=="Create new list"|sumBinListTMP=="Edit current list"){
                        #for(ibinRegions in 1:nrow(mrbin.env$mrbinparam$specialBinList)){
                        if(is.null(mrbin.env$mrbinparam$sumBinList)){
                                mrbin.env$mrbinparam$sumBinList<-matrix(ncol=4,nrow=0,dimnames=list(NULL,c("left","right","top","bottom")))
                        }
                        ibinRegions <- 1
                        adjbinRegion<-""
                        addbinRegion<-""
                        adjbinRegionAccept<-""
                        #while(addAreasFlag&!stopTMP){
                        while(ibinRegions <= (nrow(mrbin.env$mrbinparam$sumBinList)+1)&!stopTMP&
                              !adjbinRegion=="Go back"&!addbinRegion=="Go back"&!addbinRegion=="No"){
                          if(!stopTMP&!adjbinRegion=="Go back"){
                            if(ibinRegions>nrow(mrbin.env$mrbinparam$sumBinList)){
                              addbinRegion<-utils::select.list(c("Yes","No","Go back"),preselect="No",
                                         title ="Add a new region?",graphics=TRUE)
                              if(length(addbinRegion)==0|addbinRegion==""){
                                stopTMP<-TRUE
                                addbinRegion<-""
                              }
                              if(!stopTMP){
                                if(addbinRegion=="Yes"){
                                  mrbin.env$paramChangeFlag<-TRUE
                                  mrbin.env$mrbinparam$sumBinList<-rbind(mrbin.env$mrbinparam$sumBinList,c(0,0,0,0))
                                }
                              }
                            }
                            if(!stopTMP&!adjbinRegion=="Go back"&!addbinRegion=="Go back"&!addbinRegion=="No"){
                              mean1<-mean(mrbin.env$mrbinparam$sumBinList[ibinRegions,1:2])
                              range1<-mrbin.env$mrbinparam$sumBinList[ibinRegions,1]-mrbin.env$mrbinparam$sumBinList[ibinRegions,2]
                              mean2<-mean(mrbin.env$mrbinparam$sumBinList[ibinRegions,3:4])
                              range2<-mrbin.env$mrbinparam$sumBinList[ibinRegions,4]-mrbin.env$mrbinparam$sumBinList[ibinRegions,3]
                              regionTMP<-c(mean1+4*range1,mean1-4*range1,mean2-4*range2,mean2+4*range2)
                              if(sum(mrbin.env$mrbinparam$sumBinList[ibinRegions,]==0)==4) regionTMP<-"all"
                              if(mrbin.env$mrbinparam$dimension=="1D"){
                                if(mrbin.env$mrbinparam$showSpectrumPreview=="Yes") plotNMR(region=regionTMP,
                                        rectangleRegions=matrix(c(mrbin.env$mrbinparam$sumBinList[ibinRegions,1],
                                                                mrbin.env$mrbinparam$sumBinList[ibinRegions,2],0,2),ncol=4),
                                        color="black",
                                        manualScale=FALSE,
                                        plotTitle=paste("Sum area\nleft=",mrbin.env$mrbinparam$sumBinList[ibinRegions,1],
                                                  "ppm, right=",mrbin.env$mrbinparam$sumBinList[ibinRegions,2],"ppm",sep=""))
                              }
                              if(mrbin.env$mrbinparam$dimension=="2D"){
                                if(mrbin.env$mrbinparam$showSpectrumPreview=="Yes") plotNMR(region=regionTMP,
                                        rectangleRegions=matrix(mrbin.env$mrbinparam$sumBinList[ibinRegions,1:4],ncol=4),
                                        color="black",
                                        manualScale=FALSE,
                                        plotTitle=paste("Sum area\nleft=",mrbin.env$mrbinparam$sumBinList[ibinRegions,1],
                                                  "ppm, right=",mrbin.env$mrbinparam$sumBinList[ibinRegions,2],"ppm",sep=""))
                              }
                              adjbinRegionAccept<-paste(paste(c("Keep left=","ppm, right=","ppm, top=","ppm, bottom=")[1:dimlength],
                                                mrbin.env$mrbinparam$sumBinList[ibinRegions,1:dimlength],collapse="",sep=""),"ppm",sep="")
                              if(sum(mrbin.env$mrbinparam$sumBinList[ibinRegions,]==0)==4){
                                adjbinRegion<-"Change..."
                              } else {
                                adjbinRegion<-utils::select.list(c(adjbinRegionAccept,"Change...","Remove","Go back"),
                                         preselect=adjbinRegionAccept,
                                         title =paste("Keep region ",ibinRegions,"?",sep=""),graphics=TRUE)
                              }
                              if(length(adjbinRegion)==0|adjbinRegion=="") stopTMP<-TRUE
                            }
                          }
                          if(adjbinRegion=="Change..."&!stopTMP&!adjbinRegion=="Go back"){
                            regionTMP<-readline(prompt=paste("Region ",ibinRegions,": left border, press enter to keep ",
                                      mrbin.env$mrbinparam$sumBinList[ibinRegions,1],": ",sep=""))
                            if(!regionTMP=="") {
                                   mrbin.env$paramChangeFlag<-TRUE
                                   mrbin.env$mrbinparam$sumBinList[ibinRegions,1]<-as.numeric(regionTMP)
                            }
                            regionTMP<-readline(prompt=paste("Region ",ibinRegions,": right border, press enter to keep ",
                                      mrbin.env$mrbinparam$sumBinList[ibinRegions,2],": ",sep=""))
                            if(!regionTMP=="") {
                                   mrbin.env$paramChangeFlag<-TRUE
                                   mrbin.env$mrbinparam$sumBinList[ibinRegions,2]<-as.numeric(regionTMP)
                            }
                            if(mrbin.env$mrbinparam$sumBinList[ibinRegions,1]<mrbin.env$mrbinparam$sumBinList[ibinRegions,2]){
                              TMP<-mrbin.env$mrbinparam$sumBinList[ibinRegions,1]
                              mrbin.env$mrbinparam$sumBinList[ibinRegions,1]<-mrbin.env$mrbinparam$sumBinList[ibinRegions,2]
                              mrbin.env$mrbinparam$sumBinList[ibinRegions,2]<-TMP
                            }
                          if(mrbin.env$mrbinparam$dimension=="2D"&!stopTMP){
                            regionTMP<-readline(prompt=paste("Region ",ibinRegions,": top border, press enter to keep ",
                                      mrbin.env$mrbinparam$sumBinList[ibinRegions,3],": ",sep=""))
                            if(!regionTMP=="") {
                                   mrbin.env$paramChangeFlag<-TRUE
                                   mrbin.env$mrbinparam$sumBinList[ibinRegions,3]<-as.numeric(regionTMP)
                            }
                            regionTMP<-readline(prompt=paste("Region ",ibinRegions,": bottom border, press enter to keep ",
                                      mrbin.env$mrbinparam$sumBinList[ibinRegions,4],": ",sep=""))
                            if(!regionTMP=="") {
                                   mrbin.env$paramChangeFlag<-TRUE
                                   mrbin.env$mrbinparam$sumBinList[ibinRegions,4]<-as.numeric(regionTMP)
                            }
                          if(mrbin.env$mrbinparam$sumBinList[ibinRegions,4]<mrbin.env$mrbinparam$sumBinList[ibinRegions,3]){
                            TMP<-mrbin.env$mrbinparam$sumBinList[ibinRegions,3]
                            mrbin.env$mrbinparam$sumBinList[ibinRegions,3]<-mrbin.env$mrbinparam$sumBinList[ibinRegions,4]
                            mrbin.env$mrbinparam$sumBinList[ibinRegions,4]<-TMP
                          }
                         }
                        }
                        if(adjbinRegion=="Remove"&!stopTMP&!adjbinRegion=="Go back"){
                          mrbin.env$mrbinparam$sumBinList<-mrbin.env$mrbinparam$sumBinList[-ibinRegions,,drop=FALSE]
                        }
                        if(adjbinRegion==adjbinRegionAccept&!stopTMP){
                          ibinRegions <- ibinRegions+1
                        }
                      }
                    }
                  }
                if(is.null(mrbin.env$mrbinparam$sumBinList)){
                  adjbinRegion<-"Go back"
                } else {
                  if(nrow(mrbin.env$mrbinparam$sumBinList)==0){
                    mrbin.env$mrbinparam$sumBinList<-NULL
                    adjbinRegion<-"Go back"
                  }
                }
                }
                }
                if(adjbinRegion=="Go back"|addbinRegion=="Go back"|sumBins=="Go back"|sumBinListTMP=="Go back"){
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
                adjNoiseRegion<-""
                if(mrbin.env$mrbinparam$verbose){
                  message("Hint: Recommended to increase statistical power")
                  utils::flush.console()
                }
                noiseRemoval<-utils::select.list(c("Yes","No","Go back"),
                                         preselect=mrbin.env$mrbinparam$noiseRemoval,
                                         title="Remove noise?",graphics=TRUE)
                if(length(noiseRemoval)==0|noiseRemoval=="") stopTMP<-TRUE
                if(!stopTMP&!noiseRemoval=="Go back"){
                  mrbin.env$mrbinparam$noiseRemoval<-noiseRemoval
                  if(mrbin.env$mrbinparam$noiseRemoval=="Yes"&!stopTMP){
                   adjNoiseRegionAcceptFlag<-FALSE
                   adjNoiseRegion<-""
                   while(!adjNoiseRegionAcceptFlag&!stopTMP&!adjNoiseRegion=="Go back"){
                    if(mrbin.env$mrbinparam$verbose){
                      message("Hint: Region should have no real signals")
                      utils::flush.console()
                    }
                    if(mrbin.env$mrbinparam$dimension=="1D"){
                        if(mrbin.env$mrbinparam$showSpectrumPreview=="Yes") plotNMR(#region=regionTMP,
                                rectangleRegions=matrix(c(mrbin.env$mrbinparam$noiseRange1d,-1000,1000),ncol=4),
                                color="black",
                                manualScale=FALSE,
                                plotTitle=paste("Noise region\nleft=",mrbin.env$mrbinparam$noiseRange1d[1],
                                          "ppm, right=",mrbin.env$mrbinparam$noiseRange1d[2],"ppm",sep=""))
                        adjNoiseRegionAccept<-paste(paste(c("Keep left=","ppm, right=","ppm, top=","ppm, bottom=")[1:dimlength],
                                          mrbin.env$mrbinparam$noiseRange1d,collapse="",sep=""),"ppm",sep="")
                      } else {
                        if(mrbin.env$mrbinparam$showSpectrumPreview=="Yes") plotNMR(#region=regionTMP,
                                rectangleRegions=matrix(mrbin.env$mrbinparam$noiseRange2d,ncol=4),
                                color="black",
                                manualScale=FALSE,
                                plotTitle=paste("Noise region\nleft=",mrbin.env$mrbinparam$noiseRange2d[1],
                                          "ppm, right=",mrbin.env$mrbinparam$noiseRange2d[2],"ppm, top=",mrbin.env$mrbinparam$noiseRange2d[3],
                                          "ppm, bottom=",mrbin.env$mrbinparam$noiseRange2d[4],"ppm",sep=""))

                        adjNoiseRegionAccept<-paste(paste(c("Keep left=","ppm, right=","ppm, top=","ppm, bottom=")[1:dimlength],
                                          mrbin.env$mrbinparam$noiseRange2d,collapse="",sep=""),"ppm",sep="")
                      }
                    adjNoiseRegion<-utils::select.list(c(adjNoiseRegionAccept,"Change...","Go back"),
                             preselect=adjNoiseRegionAccept,
                             title ="Keep noise region?",graphics=TRUE)
                    if(length(adjNoiseRegion)==0|adjNoiseRegion=="") stopTMP<-TRUE
                    if(!stopTMP){
                          if(adjNoiseRegion==adjNoiseRegionAccept) adjNoiseRegionAcceptFlag<-TRUE
                          if(adjNoiseRegion=="Change..."&!adjNoiseRegion=="Go back"){
                           if(mrbin.env$mrbinparam$dimension=="1D"){
                            regionTMP<-readline(prompt=paste("Noise: left border, press enter to keep ",
                                      mrbin.env$mrbinparam$noiseRange1d[1],": ",sep=""))
                            if(!regionTMP=="") {
                                   mrbin.env$paramChangeFlag<-TRUE
                                   mrbin.env$mrbinparam$noiseRange1d[1]<-as.numeric(regionTMP)
                            }
                            regionTMP<-readline(prompt=paste("Noise: right border, press enter to keep ",
                                      mrbin.env$mrbinparam$noiseRange1d[2],": ",sep=""))
                            if(!regionTMP=="") {
                                   mrbin.env$paramChangeFlag<-TRUE
                                   mrbin.env$mrbinparam$noiseRange1d[2]<-as.numeric(regionTMP)
                            }
                            if(mrbin.env$mrbinparam$noiseRange1d[1]<mrbin.env$mrbinparam$noiseRange1d[2]){
                              TMP<-mrbin.env$mrbinparam$noiseRange1d[1]
                              mrbin.env$mrbinparam$noiseRange1d[1]<-mrbin.env$mrbinparam$noiseRange1d[2]
                              mrbin.env$mrbinparam$noiseRange1d[2]<-TMP
                            }
                           }
                          if(mrbin.env$mrbinparam$dimension=="2D"&!stopTMP){
                            regionTMP<-readline(prompt=paste("Noise: left border, press enter to keep ",
                                      mrbin.env$mrbinparam$noiseRange2d[1],": ",sep=""))
                            if(!regionTMP=="") {
                                   mrbin.env$paramChangeFlag<-TRUE
                                   mrbin.env$mrbinparam$noiseRange2d[1]<-as.numeric(regionTMP)
                            }
                            regionTMP<-readline(prompt=paste("Noise: right border, press enter to keep ",
                                      mrbin.env$mrbinparam$noiseRange2d[2],": ",sep=""))
                            if(!regionTMP=="") {
                                   mrbin.env$paramChangeFlag<-TRUE
                                   mrbin.env$mrbinparam$noiseRange2d[2]<-as.numeric(regionTMP)
                            }
                            if(mrbin.env$mrbinparam$noiseRange2d[1]<mrbin.env$mrbinparam$noiseRange2d[2]){
                              TMP<-mrbin.env$mrbinparam$noiseRange2d[1]
                              mrbin.env$mrbinparam$noiseRange2d[1]<-mrbin.env$mrbinparam$noiseRange2d[2]
                              mrbin.env$mrbinparam$noiseRange2d[2]<-TMP
                            }
                            regionTMP<-readline(prompt=paste("Noise: top border, press enter to keep ",
                                      mrbin.env$mrbinparam$noiseRange2d[3],": ",sep=""))
                            if(!regionTMP=="") {
                                   mrbin.env$paramChangeFlag<-TRUE
                                   mrbin.env$mrbinparam$noiseRange2d[3]<-as.numeric(regionTMP)
                            }
                            regionTMP<-readline(prompt=paste("Noise: bottom border, press enter to keep ",
                                      mrbin.env$mrbinparam$noiseRange2d[4],": ",sep=""))
                            if(!regionTMP=="") {
                                   mrbin.env$paramChangeFlag<-TRUE
                                   mrbin.env$mrbinparam$noiseRange2d[4]<-as.numeric(regionTMP)
                            }
                          if(mrbin.env$mrbinparam$noiseRange2d[4]<mrbin.env$mrbinparam$noiseRange2d[3]){
                            TMP<-mrbin.env$mrbinparam$noiseRange2d[3]
                            mrbin.env$mrbinparam$noiseRange2d[3]<-mrbin.env$mrbinparam$noiseRange2d[4]
                            mrbin.env$mrbinparam$noiseRange2d[4]<-TMP
                          }
                         }
                        }
                    }
                    }
                    if(!stopTMP&mrbin.env$mrbinparam$verbose){
                      message("Hint: Increase to remove more noise, decrease to allow low intensity signals")
                      utils::flush.console()
                    }
                    if(mrbin.env$mrbinparam$dimension=="1D"&!stopTMP&!adjNoiseRegion=="Go back"){
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
                    if(mrbin.env$mrbinparam$dimension=="2D"&!stopTMP&!adjNoiseRegion=="Go back"){
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
                    if(!stopTMP&mrbin.env$mrbinparam$verbose){
                      message("Hint: 0.75 if metabolites are always present (serum), 0.2 if \n metabolites are absent in some samples (urine, cell culture)")
                      utils::flush.console()
                    }
                    if(!stopTMP&!SNRTMP=="Go back"&!adjNoiseRegion=="Go back"){
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
                if(noiseRemoval=="Go back"|SNRTMP=="Go back"|noiseTMP=="Go back"|adjNoiseRegion=="Go back"){
                   selectStep<-selectStep-2
                }
                if(!stopTMP) selectStep<-selectStep+1
              }
            }
            if(selectStep==11){
              #Crop HSQCs
              if(mrbin.env$mrbinparam$dimension=="2D"&!stopTMP){
                 if(!stopTMP&mrbin.env$mrbinparam$verbose){
                   message("Hint: This removes additional noise")
                   utils::flush.console()
                 }
                 if(mrbin.env$mrbinparam$showSpectrumPreview=="Yes") plotNMR(region="all",
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
              #Trim zeros
              if(!stopTMP){
                if(mrbin.env$mrbinparam$verbose){
                  message("Hint: Recommended after removing solvents or noise")
                  utils::flush.console()
                }
                trimZeros<-utils::select.list(c("Yes","No","Go back"),
                                         preselect=mrbin.env$mrbinparam$trimZeros,
                                         title="Trim zero-value bins?",graphics=TRUE)
                if(length(trimZeros)==0|trimZeros=="") stopTMP<-TRUE
                if(!stopTMP&!trimZeros=="Go back"){
                  mrbin.env$mrbinparam$trimZeros<-trimZeros
                }
                if(trimZeros=="Go back"&!stopTMP){
                   if(mrbin.env$mrbinparam$dimension=="2D") selectStep<-selectStep-2
                   if(mrbin.env$mrbinparam$dimension=="1D") selectStep<-selectStep-3
                }
                if(!stopTMP) selectStep<-selectStep+1
              }
            }
            if(selectStep==13){
              #PQN scaling?
              if(!stopTMP){
                PQNScalingIgnoreSugar<-""
                if(!stopTMP&mrbin.env$mrbinparam$verbose){
                  message("Hint: Recommended for urine and tissue extracts")
                  utils::flush.console()
                }
                PQNScaling<-utils::select.list(c("Yes","No","Go back"),
                                        preselect=mrbin.env$mrbinparam$PQNScaling,
                                        title = "PQN normalization?",graphics=TRUE)
                if(length(PQNScaling)==0|PQNScaling=="") stopTMP<-TRUE
                if(!stopTMP&!PQNScaling=="Go back"){
                  mrbin.env$mrbinparam$PQNScaling<-PQNScaling
                  if(mrbin.env$mrbinparam$PQNScaling=="Yes"){
                    if(mrbin.env$mrbinparam$verbose){
                      message("Hint: Improves PQN but works only for 1H and 1H-13C spectra")
                      utils::flush.console()
                    }
                    PQNScalingIgnoreSugar<-utils::select.list(c("Yes","No","Go back"),
                                            preselect=mrbin.env$mrbinparam$PQNIgnoreSugarArea,
                                            title = "Ignore glucose for PQN?",graphics=TRUE)
                    if(length(PQNScalingIgnoreSugar)==0|PQNScalingIgnoreSugar=="") stopTMP<-TRUE
                    if(!stopTMP&!PQNScalingIgnoreSugar=="Go back"){
                      mrbin.env$mrbinparam$PQNScaling<-PQNScaling
                    }
                  }
                }
                if(!stopTMP&(PQNScaling=="Go back"|PQNScalingIgnoreSugar=="Go back")){
                   selectStep<-selectStep-2
                }
                if(!stopTMP) selectStep<-selectStep+1
              }
            }
            if(selectStep==14){
              if(!stopTMP&mrbin.env$mrbinparam$verbose){
                message("Hint: Recommended in most cases")
                utils::flush.console()
              }
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
            if(selectStep==15){
              #Log scaling?
              if(!stopTMP){
                if(!stopTMP&mrbin.env$mrbinparam$verbose){
                  message("Hint: Recommended if data needs to be normal or fold changes are of interest")
                  utils::flush.console()
                }
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
            if(selectStep==16){
              #Define sample names
              if(!stopTMP){
                if(!stopTMP&mrbin.env$mrbinparam$verbose){
                  message("Hint: If only EXPNO differs choose Folder names and EXPNO")
                  utils::flush.console()
                }
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
            if(selectStep==17){
              #Plot results
              if(!stopTMP){
                if(!stopTMP&mrbin.env$mrbinparam$verbose){
                  message("Hint: Recommended for quality control")
                  utils::flush.console()
                }
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
                    if(!stopTMP&mrbin.env$mrbinparam$verbose){
                      message("Hint: Recommended for a nicer plot, make sure names are unique")
                      utils::flush.console()
                    }
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
            if(selectStep==18){
              #Define groups
              if(!stopTMP){
                if(!stopTMP&mrbin.env$mrbinparam$verbose){
                  message("Hint: Define treatment groups for color-coded PCA")
                  utils::flush.console()
                }
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
                                             paste(mrbin.env$mrbinparam$Factors[1:min(3,length(mrbin.env$mrbinparam$Factors))],
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
     if(selectStep==19){
       #Save output files to hard drive?
       if(!stopTMP){
         saveFilesTMP2<-"Select new folder and file name"
         #if(!stopTMP&mrbin.env$mrbinparam$verbose) message("Hint: Recommended if you want to reuse the data later")
         saveFilesTMP<-utils::select.list(c("Yes","No","Go back"),
                     preselect=mrbin.env$mrbinparam$saveFiles,
                     title ="Save output to disk?",graphics=TRUE)
         if(length(saveFilesTMP)==0|saveFilesTMP=="") stopTMP<-TRUE
         if(!stopTMP&!saveFilesTMP=="Go back"){
          mrbin.env$mrbinparam$saveFiles<-saveFilesTMP
          if(mrbin.env$mrbinparam$saveFiles=="Yes"&!stopTMP){
            if(!is.null(mrbin.env$mrbinparam$outputFileName)){
             keepFileTMP<-paste("Keep ",mrbin.env$mrbinparam$outputFileName,sep="")
             saveFilesTMP2<-utils::select.list(c(keepFileTMP,"Select new folder and file name","Go back"),
                         preselect=keepFileTMP,
                         title ="Keep file name and folder?",graphics=TRUE)
             if(length(saveFilesTMP2)==0|saveFilesTMP=="") stopTMP<-TRUE
            }
            if(!stopTMP&saveFilesTMP2=="Select new folder and file name"){
               enterFoldersTMP<-readline(prompt="Enter starting folder path. (Examples: Windows: \"C:\\\", Apple: \"/\") : ")
               if(enterFoldersTMP=="") saveFilesTMP2<-"Go back"
               if(!saveFilesTMP2=="Go back"){
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
        }
        if(!stopTMP&(saveFilesTMP=="Go back"|saveFilesTMP2=="Go back")){
           selectStep<-selectStep-2
        }
        if(!stopTMP) selectStep<-selectStep+1
      }
     }
     if(selectStep==20){
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
                   #createBinNumbers()
                   #createBinRegions()
                   #if(mrbin.env$mrbinTMP$nbins==ncol(mrbin.env$mrbinTMP$binsRaw)){
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
                   #}
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
     if(selectStep==21){
       if(!stopTMP)  lastStepDone<-TRUE
     }
    }
   }
 #}
 if(!stopTMP){
   if(mrbin.env$mrbinparam$verbose){
       printParameters()
   }
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
#'          fixNegatives="No",logTrafo="No",signal_to_noise2D=20,
#'          NMRfolders=c(system.file("extdata/1/12/pdata/10",package="mrbin"),
#'                       system.file("extdata/2/12/pdata/10",package="mrbin"))))
#' mrbinrun()

mrbinrun<-function(){
  if(!exists("mrbin.env", mode="environment")) .onLoad()
  if(!is.null(mrbin.env$mrbinparam$NMRfolders)){
    #if(mrbin.env$mrbinparam$PCA=="Yes"){
    mrbin.env$mrbinTMP$currentFolder<-mrbin.env$mrbinparam$NMRfolders[1]
    readNMR2()
    #}    #if(!"mrbinparam"%in%ls(envir = mrbin.env))    resetEnv()
    #if(mrbin.env$mrbinparam$createBins=="Yes")
    createBinNumbers()
    #if(mrbin.env$mrbinparam$binMethod=="Rectangular bins"&mrbin.env$mrbinparam$createBins=="Yes") createBinRegions()
    #if(mrbin.env$mrbinparam$createBins=="Yes")
    createBinRegions()
    mrbin.env$mrbinparam$numberOfFeaturesRaw<-nrow(mrbin.env$mrbinTMP$binRegions)
    if(mrbin.env$mrbinparam$removeSolvent=="Yes") removeSolvent()
    if(mrbin.env$mrbinparam$removeAreas=="Yes") removeAreas()
    if(mrbin.env$mrbinparam$sumBins=="Yes") sumBins()
    if(mrbin.env$mrbinparam$cropHSQC=="Yes"&mrbin.env$mrbinparam$dimension=="2D") cropNMR()
    if(mrbin.env$mrbinparam$createBins=="Yes") binMultiNMR()
    if(mrbin.env$mrbinparam$trimZeros=="Yes") trimZeros()
    if(mrbin.env$mrbinparam$noiseRemoval=="Yes") removeNoise()
    #if(mrbin.env$mrbinparam$createBins=="Yes")
    createBinNames()
    if(mrbin.env$mrbinparam$fixNegatives=="Yes") atnv()
    if(mrbin.env$mrbinparam$PQNScaling=="Yes") PQNScaling()
    if(mrbin.env$mrbinparam$logTrafo=="Yes") logTrafo()
    #Sort bins
    if(nrow(mrbin.env$bins)>1){
      binRegionsTMP<-mrbin.env$mrbinTMP$binRegions
      mrbin.env$bins<-mrbin.env$bins[,rev(order(binRegionsTMP[,3])),drop=FALSE]
      binRegionsTMP<-binRegionsTMP[rev(order(binRegionsTMP[,3])),,drop=FALSE]
      mrbin.env$bins<-mrbin.env$bins[,rev(order(binRegionsTMP[,1])),drop=FALSE]
      #binRegionsTMP<-binRegionsTMP[,order(binRegionsTMP[,1])]
    }
    if(mrbin.env$mrbinparam$saveFiles=="Yes"){
      dput(mrbin.env$mrbinparam, file = paste(mrbin.env$mrbinparam$outputFileName,".txt",sep=""))
      utils::write.csv(mrbin.env$bins, file = paste(mrbin.env$mrbinparam$outputFileName,"bins.csv",sep=""))
    }
    if(mrbin.env$mrbinparam$PCA=="Yes"){
      readNMR2()
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
     if(!is.null(mrbin.env$mrbinparam$numberOfFeaturesAfterTrimmingZeros)&mrbin.env$mrbinparam$trimZeros=="Yes"){
          resultOutputTMP<-c(resultOutputTMP,
             "Number of bins after trimming zero-value bins: ",mrbin.env$mrbinparam$numberOfFeaturesAfterTrimmingZeros,"\n")
     }
     if(!is.null(mrbin.env$mrbinparam$numberOfFeaturesAfterNoiseRemoval)&mrbin.env$mrbinparam$noiseRemoval=="Yes"){
          resultOutputTMP<-c(resultOutputTMP,
             "Number of bins after noise removal: ",mrbin.env$mrbinparam$numberOfFeaturesAfterNoiseRemoval,"\n")
     }
     resultOutputTMP<-paste(resultOutputTMP,sep="")
     if(mrbin.env$mrbinparam$verbose){
       #printParameters()
       message(resultOutputTMP, appendLF = FALSE)
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
            vectorSymbol1<-"matrix(c(\n  "
            if(is.null(rownames(mrbin.env$mrbinparam[[i]]))){
              vectorSymbol2<-paste("\n  ),ncol=",ncol(mrbin.env$mrbinparam[[i]]),",byrow=TRUE)",sep="")
            } else {
              vectorSymbol2<-paste("\n  ),ncol=",ncol(mrbin.env$mrbinparam[[i]]),
                             ",dimnames=list(c(\n  \"",
                             paste(rownames(mrbin.env$mrbinparam[[i]]),sep="\",\"",collapse="\",\""),
                             "\"\n  ),NULL),byrow=TRUE",
                             ")",sep="")
            }
          }
        }
        if(length(mrbin.env$mrbinparam[[i]])>1){
          if(is.vector(mrbin.env$mrbinparam[[i]])){
            vectorSymbol1<-"c("
            vectorSymbol2<-")"
            if(i=="NMRfolders"){
              vectorSymbol1<-"c(\n  "
              vectorSymbol2<-"\n )"
            }
          }
        }
        if(is.factor(mrbin.env$mrbinparam[[i]])){
          vectorSymbol1<-"factor(c("
          vectorSymbol2<-"))"
          if(nchar(paste(mrbin.env$mrbinparam[[i]],sep=",",collapse=","))>60){
            vectorSymbol1<-"factor(c(\n  "
            vectorSymbol2<-"\n  ))"
          }
          sepSymbol<-"\""
        }
        sepTMP<-paste(sepSymbol,",",returnSymbol,sepSymbol,sep="")
        sepTMPReturn<-paste(sepSymbol,",","\n  ",sepSymbol,sep="")
        if(is.null(mrbin.env$mrbinparam[[i]])){
            valueTMP<-"NULL"
        } else {
          if(is.factor(mrbin.env$mrbinparam[[i]])){
            lengthCutOff2<-8
            if(length(mrbin.env$mrbinparam[[i]])<=lengthCutOff2){
              valueTMP<-paste(as.character(mrbin.env$mrbinparam[[i]]),sep=sepTMP,collapse=sepTMP)
            } else {
              valueTMP<-NULL
              for(iFactorLengthTMP in 1:ceiling(length(mrbin.env$mrbinparam[[i]])/lengthCutOff2)){
                sepFactorTMPTMP<-sepTMPReturn
                if(iFactorLengthTMP==ceiling(length(mrbin.env$mrbinparam[[i]])/lengthCutOff2)){
                  sepFactorTMPTMP<-""
                }
                valueTMP<-paste(valueTMP,
                          paste(as.character(mrbin.env$mrbinparam[[i]])[((iFactorLengthTMP-1)*lengthCutOff2+1):
                                 min(length(mrbin.env$mrbinparam[[i]]),(iFactorLengthTMP-1)*lengthCutOff2+lengthCutOff2)],
                            sep=sepTMP,collapse=sepTMP),sepFactorTMPTMP,
                          sep="",collapse="")
              }
            }
          } else {
            lengthCutOff<-12
            if(length(mrbin.env$mrbinparam[[i]])<=lengthCutOff){
              valueTMP<-paste(t(mrbin.env$mrbinparam[[i]]),sep=sepTMP,collapse=sepTMP)
            } else {
              valueTMP<-NULL
              for(iFactorLengthTMP in 1:ceiling(length(mrbin.env$mrbinparam[[i]])/lengthCutOff)){
                sepFactorTMPTMP<-sepTMPReturn
                if(iFactorLengthTMP==ceiling(length(mrbin.env$mrbinparam[[i]])/lengthCutOff)){
                  sepFactorTMPTMP<-""
                }
                valueTMP<-paste(valueTMP,
                          paste(t(mrbin.env$mrbinparam[[i]])[((iFactorLengthTMP-1)*lengthCutOff+1):
                                 min(length(mrbin.env$mrbinparam[[i]]),(iFactorLengthTMP-1)*lengthCutOff+lengthCutOff)],
                            sep=sepTMP,collapse=sepTMP),sepFactorTMPTMP,
                          sep="",collapse="")
              }
            }
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
      message("\nTo recreate this data set, use the following code:\n\n##################################################\n\n", appendLF = FALSE)
      #Necessary as message() cuts of strings after around 8187 characters (Windows 10)
      messageMaxLength<-2000
      for(iLengthMessage in 1:ceiling(nchar(printTMP)/messageMaxLength)){
        message(substr(printTMP,
                      (iLengthMessage-1)*messageMaxLength+1,
                      (iLengthMessage-1)*messageMaxLength+messageMaxLength),
                appendLF = FALSE)
      }
      message("\n##################################################\n\nTo recreate this data set, use the code above this line.\n", appendLF = FALSE)
      if(mrbin.env$mrbinparam$saveFiles=="Yes"){
        message("To load the saved data set, use the following code:\n\n", appendLF = FALSE)
        message(" data<-read.csv(\n  \"",mrbin.env$mrbinparam$outputFileName,
                "bins.csv\",\n", appendLF = FALSE)
        message("  check.names = FALSE,row.names=1)\n\n", appendLF = FALSE)
      }
      utils::flush.console()
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
    diffSet2<-setdiff(names(parameters)[!sapply(parameters, is.null)],names(mrbin.env$mrbinparam_copy))
    intersectSet<-intersect(names(parameters),names(mrbin.env$mrbinparam_copy))
    if(length(diffSet2)>0){
       warning(paste("Unexpected parameters: ",
           paste(diffSet2,sep=", ", collapse=", "),"\n",
           "These parameters are not used. Potentially they were created in a different mrbin version.",
           sep=""))
    }
    if(length(diffSet3)>0){
       if(mrbin.env$mrbinparam$verbose) message(paste("Current values are used for missing parameters: ",
           paste(diffSet3,sep=", ", collapse=", "),"\n",sep=""), appendLF = FALSE)
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
#'                     binwidth1D=0.05,signal_to_noise1D=50, verbose=TRUE, PCA="No",
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
#' @param  spectrumNumber If provided, this number will be used; defaults to NULL
#' @return {None}
#' @export
#' @examples
#' \donttest{ setCurrentSpectrum(spectrumNumber=1) }

setCurrentSpectrum<-function(spectrumNumber=NULL){
     if(is.null(spectrumNumber)){
       newCurrent<-utils::select.list(mrbin.env$mrbinparam$NMRfolders,preselect=mrbin.env$mrbinTMP$currentFolder,
                                          title="Select new current spectrum.",graphics=TRUE)
     } else {
       if(length(mrbin.env$mrbinparam$NMRfolders)>=spectrumNumber){
         newCurrent<- mrbin.env$mrbinparam$NMRfolders[spectrumNumber]
       } else {
         newCurrent<-""
       }
     }
     if(!newCurrent==""){
          mrbin.env$mrbinTMP$currentFolder<-newCurrent
          readNMR2()
     }
}

#' A function for removing a spectrum.
#'
#' This function lets the user pick spectra from a list for removal from data
#' analysis. This function is meant only for use within the mrbin function.
#' @return {None}
#' @keywords internal
#' @noRd
#' @examples
#' \donttest{ removeSpectrum() }

removeSpectrum<-function(){
 if(!is.null(mrbin.env$mrbinparam$NMRfolders)){
    listTMP<-utils::select.list(mrbin.env$mrbinparam$NMRfolders,preselect = NULL, multiple = TRUE,title ="Remove spectra, cancel to keep all",graphics=TRUE)
    continueTMP<-FALSE
    if(length(listTMP)>0){
      if(length(listTMP)>1){
        continueTMP<-TRUE
      } else {
        if(!listTMP==""){
          continueTMP<-TRUE
        }
      }
      if(continueTMP){
         if(length(mrbin.env$mrbinparam$Factors)==length(mrbin.env$mrbinparam$NMRfolders)){
           mrbin.env$mrbinparam$Factors<-mrbin.env$mrbinparam$Factors[-which(mrbin.env$mrbinparam$NMRfolders%in%listTMP),drop=FALSE]
         }
         if(!is.null(mrbin.env$bins)){
           if(nrow(mrbin.env$bins)==length(mrbin.env$mrbinparam$NMRfolders)){
            #if((nrow(mrbin.env$bins)-length(listTMP))==1){
               mrbin.env$bins<-mrbin.env$bins[-which(rownames(mrbin.env$bins)%in%listTMP),,drop=FALSE]
              #rownamesTMP<-rownames(mrbin.env$bins)[-which(rownames(mrbin.env$bins)%in%listTMP)]
              #colnamesTMP<-colnames(mrbin.env$bins)
              #mrbin.env$bins<-matrix(mrbin.env$bins[-which(rownames(mrbin.env$bins)%in%listTMP),],nrow=1)
              #rownames(mrbin.env$bins)<-rownamesTMP
              #colnames(mrbin.env$bins)<-colnamesTMP
            #} else {
            #   mrbin.env$bins<-mrbin.env$bins[-which(rownames(mrbin.env$bins)%in%listTMP),,drop=FALSE]
            #}
           }
         }
         mrbin.env$mrbinparam$NMRfolders<-mrbin.env$mrbinparam$NMRfolders[-which(mrbin.env$mrbinparam$NMRfolders%in%listTMP),drop=FALSE]
      }
    }
 }
}

#' A function for setting group members.
#'
#' This function lets the user pick samples from a list to assign them to
#' groups. This function is meant only for use within the mrbin function.
#' @return {None}
#' @keywords internal
#' @noRd
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
#' @param keep Keep and add to current list of spectra, or create an all new list
#' @return An invisible list of folder names, or "Go back" or "stop"
#' @keywords internal
#' @noRd
#' @examples
#' \donttest{ selectFolders() }

selectFolders<-function(keep=FALSE){#Select NMR spectral folders
      selectionFolders<-""
      if(mrbin.env$mrbinparam$NMRvendor=="Bruker"){
          selectionFolders<-selectBrukerFolders(keep)
      }  else {
          stop(paste("No folder selection function defined for vendor ",mrbin.env$mrbinparam$NMRvendor,".\n",sep=""))
      }
      invisible(selectionFolders)
}

#' A function for selecting Bruker NMR data folders.
#'
#' This function lets the user set NMR data folders interactively (for Bruker data). This function
#' is meant only for use within the mrbin function.
#' @param keep Keep and add to current list of spectra, or create an all new list
#' @return An invisible list of folder names, or "Go back" or "stop"
#' @keywords internal
#' @noRd
#' @examples
#' \donttest{ selectBrukerFolders() }

selectBrukerFolders<-function(keep=FALSE){#Select Bruker NMR spectral folders
  selectionFolders<-""
  if(!keep){
    mrbin.env$mrbinparam$NMRfolders<-NULL
  }
  NMRfoldersTMP<-NULL
  datanameDict<-c("1r","2rr")
  names(datanameDict)<-c("1D","2D")
  datanameTmp<-datanameDict[mrbin.env$mrbinparam$dimension]
  singleFolderFlag<-FALSE
  enterFolders<-utils::select.list(c("Browse...",#"Enter parent folder path manually",
                 "Go back"),
                 preselect="Browse...",title="Set NMR parent folder:",graphics=TRUE)
  if(enterFolders==""|length(enterFolders)==0){
       if(is.null(mrbin.env$mrbinparam$NMRfolders)) selectionFolders<-"stop"
  } else {
    if(enterFolders=="Go back"){
       selectionFolders<-"Go back"
    } else {
      folderPrompt<-"Enter folder path: "
      if(enterFolders=="Browse..."){
        folderPrompt<-"Enter starting folder path (example: Windows: \"C:\\\", Apple: \"/\"): "
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
       if(!singleFolderFlag) message("After reaching the NMR parent folder, click OK WITHOUT selecting a folder.\n", appendLF = FALSE)
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
                       list.filesTMP<-list.files(spectrum_path_list[i])
                       if(datanameTmp%in%list.filesTMP&"title"%in%list.filesTMP){
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
                     list.filesTMP<-list.files(spectrum_path_list[i])
                     if(datanameTmp%in%list.filesTMP&"title"%in%list.filesTMP){
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
               #if(mrbin.env$mrbinparam$verbose) message(paste("Adding to list:\n",paste(NMRfoldersTMP,"\n",
               #                                     sep="",collapse="")), appendLF = FALSE)
               addSpectrumTMP<-TRUE
               while(addSpectrumTMP){
                 yesornoPreSelect<-paste("Keep current spectra list (",length(mrbin.env$mrbinparam$NMRfolders)," spectra)",sep="")
                 yesorno<-utils::select.list(c(yesornoPreSelect,"Add additional spectra","Remove spectra from list"),
                           preselect=yesornoPreSelect,multiple=FALSE,title="Add additional spectra?",graphics=TRUE)
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
                 if(yesorno==yesornoPreSelect) {
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
      } else {
        selectionFolders<-"Go back"
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
#' @keywords internal
#' @noRd
#' @examples
#' \donttest{ binMultiNMR() }

binMultiNMR<-function(){
 if(!is.null(mrbin.env$mrbinparam$NMRfolders)){
    mrbin.env$bins<-NULL
    mrbin.env$mrbinTMP$binNames<-NULL
    #mrbin.env$mrbinTMP$binTMP<-NULL
    #Open and bin all spectra
    #mrbin.env$mrbinTMP$binsRaw<-NULL
    #mrbin.env$mrbinparam$noise_level_Raw<-NULL
    #mrbin.env$mrbinTMP$noise_level<-NULL
    #mrbin.env$mrbinTMP$meanNumberOfPointsPerBin<-NULL
    if(mrbin.env$mrbinparam$verbose){
      message("Binning spectra... ", appendLF = FALSE)
      utils::flush.console()
    }
    #listProgressTMP2<-(1:10)/10
    #Before binning first spectrum
    mrbin.env$mrbinparam$noise_level_Raw<-rep(NA,length(mrbin.env$mrbinparam$NMRfolders))
    mrbin.env$mrbinTMP$noise_level<-matrix(rep(NA,length(mrbin.env$mrbinparam$NMRfolders)*nrow(mrbin.env$mrbinTMP$binRegions)),ncol=nrow(mrbin.env$mrbinTMP$binRegions))
    mrbin.env$mrbinTMP$meanNumberOfPointsPerBin<-matrix(rep(NA,length(mrbin.env$mrbinparam$NMRfolders)*nrow(mrbin.env$mrbinTMP$binRegions)),ncol=nrow(mrbin.env$mrbinTMP$binRegions))
    mrbin.env$mrbinTMP$binTMP<-NULL
    mrbin.env$mrbinTMP$currentFolder<-mrbin.env$mrbinparam$NMRfolders[1]
    useParallel<-FALSE
    if(mrbin.env$mrbinparam$tryParallel){
      if(requireNamespace("parallel",quietly=TRUE)){
        try(cluster<-parallel::makeCluster(parallel::detectCores()),silent=TRUE)
        #Test cluster formation
        test<-try(parallel::clusterEvalQ(cluster,1+1),silent=TRUE)
        if(is.list(test)){
         useParallel<-TRUE
        } else {
          try(parallel::stopCluster(cluster),silent=TRUE)
        }
      } else {
         warning("Package parallel not found, using regular (slower) mode.")
      }
    }
    if(useParallel){
      #try(parallel::clusterExport(cluster, "mrbin.env"),silent=TRUE)
      try(
        parallel::clusterExport(cluster, c(
          "readNMR","readBruker","referenceScaling","removeSolvent2",
          "removeAreas2","binSingleNMR","calculateNoise",
          "checkBaseline"))
      ,silent=TRUE)
      try(
        binData<-parallel::parLapply(cluster,
          mrbin.env$mrbinparam$NMRfolders,binMultiNMR2,
          dimension=mrbin.env$mrbinparam$dimension,
          binRegions=mrbin.env$mrbinTMP$binRegions,
          referenceScaling=mrbin.env$mrbinparam$referenceScaling,
          removeSolvent=mrbin.env$mrbinparam$removeSolvent,
          removeAreas=mrbin.env$mrbinparam$removeAreas,
          reference1D=mrbin.env$mrbinparam$reference1D,
          reference2D=mrbin.env$mrbinparam$reference2D,
          solventRegion=mrbin.env$mrbinparam$solventRegion,
          removeAreaList=mrbin.env$mrbinparam$removeAreaList,
          NMRvendor=mrbin.env$mrbinparam$NMRvendor,
          noiseRange1d=mrbin.env$mrbinparam$noiseRange1d,
          noiseRange2d=mrbin.env$mrbinparam$noiseRange2d,
          binMethod=mrbin.env$mrbinparam$binMethod,
          useAsNames=mrbin.env$mrbinparam$useAsNames
        )
      ,silent=TRUE)
      try(parallel::stopCluster(cluster),silent=TRUE)
      if(exists("binData")){
        if(!is.list(binData)){
           useParallel<-FALSE
           warning("Parallel computing not successful, using regular mode.")
          }
      } else {
         useParallel<-FALSE
         warning("Parallel computing did not succeed, using regular mode.")
      }
    }
    if(!useParallel){
      binData<-lapply(
        mrbin.env$mrbinparam$NMRfolders,binMultiNMR2,
        dimension=mrbin.env$mrbinparam$dimension,
        binRegions=mrbin.env$mrbinTMP$binRegions,
        referenceScaling=mrbin.env$mrbinparam$referenceScaling,
        removeSolvent=mrbin.env$mrbinparam$removeSolvent,
        removeAreas=mrbin.env$mrbinparam$removeAreas,
        reference1D=mrbin.env$mrbinparam$reference1D,
        reference2D=mrbin.env$mrbinparam$reference2D,
        solventRegion=mrbin.env$mrbinparam$solventRegion,
        removeAreaList=mrbin.env$mrbinparam$removeAreaList,
        NMRvendor=mrbin.env$mrbinparam$NMRvendor,
        noiseRange1d=mrbin.env$mrbinparam$noiseRange1d,
        noiseRange2d=mrbin.env$mrbinparam$noiseRange2d,
        binMethod=mrbin.env$mrbinparam$binMethod,
        useAsNames=mrbin.env$mrbinparam$useAsNames
        )
    }
    if(mrbin.env$mrbinparam$verbose) message("done.\n", appendLF = FALSE)
    utils::flush.console()
    mrbin.env$mrbinTMP$binsRaw<-matrix(rep(0,nrow(mrbin.env$mrbinTMP$binRegions)*#length(mrbin.env$mrbinTMP$binTMP)*
                                      length(mrbin.env$mrbinparam$NMRfolders)),
                                      nrow=length(mrbin.env$mrbinparam$NMRfolders))
    currentSpectrumNameTMP<-paste("TemporaryRowName_",1:length(mrbin.env$mrbinparam$NMRfolders),sep="")
    for(ibinData in 1:length(binData)){
      #for(iWarning in 1:length(binData[[ibinData]]$warningMessage)){
        if(!is.null(binData[[ibinData]]$warningMessage)){
          warning(binData[[ibinData]]$warningMessage)
          #save warning messages to parameters
           mrbin.env$mrbinparam$warningMessages<-
              c(mrbin.env$mrbinparam$warningMessages,
                binData[[ibinData]]$warningMessage)
        }
      #}
      mrbin.env$mrbinTMP$binsRaw[ibinData,]<-binData[[ibinData]]$binTMP
      mrbin.env$mrbinTMP$meanNumberOfPointsPerBin[ibinData,]<-binData[[ibinData]]$meanNumberOfPointsPerBin_TMP#mrbin.env$mrbinTMP$meanNumberOfPointsPerBin_TMP
      mrbin.env$mrbinparam$noise_level_Raw[ibinData]<-binData[[ibinData]]$noise_level_Raw_TMP
      mrbin.env$mrbinTMP$noise_level[ibinData,]<-binData[[ibinData]]$noise_level_TMP
      mrbin.env$mrbinparam$noise_level_adjusted[ibinData]<-median(binData[[ibinData]]$noise_level_TMP)
      currentSpectrumNameTMP[ibinData]<-binData[[ibinData]]$currentSpectrumName
    }
    i_currentSpectrumNameTMP<-2
    while(sum(duplicated(currentSpectrumNameTMP))>0){
       currentSpectrumNameTMP[duplicated(currentSpectrumNameTMP)]<-paste(
            currentSpectrumNameTMP[duplicated(currentSpectrumNameTMP)],".",
            i_currentSpectrumNameTMP,sep="")
       i_currentSpectrumNameTMP<-i_currentSpectrumNameTMP+1
    }
    if(i_currentSpectrumNameTMP>2){
      warning(paste(
        "Renamed duplicate spectrum titles. Please use a different naming method.",
        sep=""))
    }
    rownames(mrbin.env$mrbinTMP$binsRaw)<-currentSpectrumNameTMP
    colnames(mrbin.env$mrbinTMP$binsRaw)<-names(binData[[1]]$binTMP)#names(mrbin.env$mrbinTMP$binTMP)
    #rownames(mrbin.env$mrbinTMP$binsRaw)<-paste("TemporaryRowName_",1:length(mrbin.env$mrbinparam$NMRfolders),sep="")

    mrbin.env$bins<-mrbin.env$mrbinTMP$binsRaw
    #mrbin.env$mrbinparam$numberOfFeaturesRaw<-ncol(mrbin.env$bins)
 }
}

#' A function for creating bin regions.
#'
#' This function creates regions for the bins. This function is
#' meant only for use within the mrbin function.
#' @return {None}
#' @keywords internal
#' @noRd
#' @examples
#' \donttest{ createBinRegions() }

createBinRegions<-function(){
   #createBinNumbers()
   if(mrbin.env$mrbinparam$binMethod=="Rectangular bins"){
       mrbin.env$mrbinTMP$binRegions<-matrix(ncol=4,
                                        nrow=mrbin.env$mrbinTMP$nbins,
                                        dimnames=list(NULL,c("left","right","top","bottom")))
       if(mrbin.env$mrbinparam$dimension=="1D"){
          decimalDigits<-max(nchar(strsplit(as.character(mrbin.env$mrbinparam$binRegion[1]),"[.]")[[1]][2]),
                             nchar(strsplit(as.character(mrbin.env$mrbinparam$binwidth1D),"[.]")[[1]][2]),
                             0,na.rm=TRUE)+4
          mrbin.env$mrbinTMP$binRegions[,1]<-round(mrbin.env$mrbinparam$binRegion[1]-((1:mrbin.env$mrbinTMP$nbins)-1)*mrbin.env$mrbinparam$binwidth1D,decimalDigits)
          mrbin.env$mrbinTMP$binRegions[,2]<-round(mrbin.env$mrbinparam$binRegion[1]-((1:mrbin.env$mrbinTMP$nbins))*mrbin.env$mrbinparam$binwidth1D,decimalDigits)
       }
       if(mrbin.env$mrbinparam$dimension=="2D"){
          decimalDigits1<-max(nchar(strsplit(as.character(mrbin.env$mrbinparam$binRegion[1]),"[.]")[[1]][2]),
                              nchar(strsplit(as.character(mrbin.env$mrbinparam$binwidth2D),"[.]")[[1]][2]),
                              0,na.rm=TRUE)+2
          decimalDigits2<-max(nchar(strsplit(as.character(mrbin.env$mrbinparam$binRegion[4]),"[.]")[[1]][2]),
                              nchar(strsplit(as.character(mrbin.env$mrbinparam$binheight),"[.]")[[1]][2]),
                              0,na.rm=TRUE)+2
          mrbin.env$mrbinTMP$binRegions[,1]<-round(rep(mrbin.env$mrbinparam$binRegion[1]-((1:mrbin.env$mrbinTMP$nbins2)-1)*mrbin.env$mrbinparam$binwidth2D,
                                             mrbin.env$mrbinTMP$nbins1),decimalDigits1)
          mrbin.env$mrbinTMP$binRegions[,2]<-round(rep(mrbin.env$mrbinparam$binRegion[1]-((1:mrbin.env$mrbinTMP$nbins2))*mrbin.env$mrbinparam$binwidth2D,
                                             mrbin.env$mrbinTMP$nbins1),decimalDigits1)
          mrbin.env$mrbinTMP$binRegions[,3]<-round(sort(rep(mrbin.env$mrbinparam$binRegion[4]-((1:mrbin.env$mrbinTMP$nbins1))*mrbin.env$mrbinparam$binheight,
                                             mrbin.env$mrbinTMP$nbins2),decreasing=TRUE),decimalDigits2)
          mrbin.env$mrbinTMP$binRegions[,4]<-round(sort(rep(mrbin.env$mrbinparam$binRegion[4]-((1:mrbin.env$mrbinTMP$nbins1)-1)*mrbin.env$mrbinparam$binheight,
                                             mrbin.env$mrbinTMP$nbins2),decreasing=TRUE),decimalDigits2)
       }
  }
  if(mrbin.env$mrbinparam$binMethod=="Custom bin list"){
    mrbin.env$mrbinTMP$binRegions<-mrbin.env$mrbinparam$specialBinList
    if(!is.matrix(mrbin.env$mrbinTMP$binRegions)) mrbin.env$mrbinTMP$binRegions<-matrix(mrbin.env$mrbinTMP$binRegions,ncol=4)
  }
}

#' A function for creating bin titles.
#'
#' This function creates titles for the bins to represent their ppm range. This function is
#' meant only for use within the mrbin function.
#' @return {None}
#' @keywords internal
#' @noRd
#' @examples
#' \donttest{ createBinNames() }

createBinNames<-function(){
   rownamesSepSignTMP<-""
   if(!is.null(rownames(mrbin.env$mrbinTMP$binRegions))){
     if(sum(rownames(mrbin.env$mrbinTMP$binRegions)=="")<nrow(mrbin.env$mrbinTMP$binRegions)){#For custom bin lists
        rownamesSepSignTMP<-";"
     }
   }
   if(mrbin.env$mrbinparam$dimension=="1D"){
     if(nrow(mrbin.env$mrbinTMP$binRegions)>1){
          namesTMP<-apply(mrbin.env$mrbinTMP$binRegions[,1:2],1,paste,collapse=",")
          rownames(mrbin.env$mrbinTMP$binRegions)<-paste(rownames(mrbin.env$mrbinTMP$binRegions),namesTMP,
                                                         sep=rownamesSepSignTMP)
     } else {
          namesTMP<-paste(mrbin.env$mrbinTMP$binRegions[,1:2],collapse=",")
          rownames(mrbin.env$mrbinTMP$binRegions)<-paste(rownames(mrbin.env$mrbinTMP$binRegions),
                                                         namesTMP,
                                                         sep=rownamesSepSignTMP)
     }
   }
   if(mrbin.env$mrbinparam$dimension=="2D"){
     if(nrow(mrbin.env$mrbinTMP$binRegions)>1){
          namesTMP<-apply(mrbin.env$mrbinTMP$binRegions,1,paste,collapse=",")
          rownames(mrbin.env$mrbinTMP$binRegions)<-paste(rownames(mrbin.env$mrbinTMP$binRegions),
                                                         namesTMP,
                                                         sep=rownamesSepSignTMP)
     } else {
          namesTMP<-paste(mrbin.env$mrbinTMP$binRegions,collapse=",")
          rownames(mrbin.env$mrbinTMP$binRegions)<-paste(rownames(mrbin.env$mrbinTMP$binRegions),
                                                         namesTMP,
                                                         sep=rownamesSepSignTMP)
     }
   }
   mrbin.env$mrbinTMP$binNames<-rownames(mrbin.env$mrbinTMP$binRegions)
   colnames(mrbin.env$bins)<-rownames(mrbin.env$mrbinTMP$binRegions)
}

#' A function for creating bin numbers.
#'
#' This function calculates numbers of bins from the chosen parameters. This function is
#' meant only for use within the mrbin function.
#' @return {None}
#' @keywords internal
#' @noRd
#' @examples
#' \donttest{ createBinNumbers() }

createBinNumbers<-function(){
   if(mrbin.env$mrbinparam$binMethod=="Custom bin list"){
          mrbin.env$mrbinTMP$nbins<-nrow(mrbin.env$mrbinparam$specialBinList)
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



#' A function for reading NMR spectra.
#'
#' This function picks the correct NMR reading function, based on vendor. This function is
#' meant only for use within the mrbin function.
#' @param onlyTitles Read only spectrum titles, but no data. Defaults to FALSE
#' @return {none}
#' @keywords internal
#' @noRd
#' @examples
#' \donttest{ readNMR2() }

readNMR2<-function(onlyTitles=FALSE){#Read NMR spectral data
 #if(!is.null(mrbin.env$mrbinTMP$currentFolder)){
   NMRdataList<-readNMR(onlyTitles=onlyTitles,
           folder=mrbin.env$mrbinTMP$currentFolder,
           dimension=mrbin.env$mrbinparam$dimension,
           NMRvendor=mrbin.env$mrbinparam$NMRvendor,
           useAsNames=mrbin.env$mrbinparam$useAsNames)
   if(!onlyTitles){
     mrbin.env$mrbinTMP$currentSpectrum<-NMRdataList$currentSpectrum
     mrbin.env$mrbinTMP$currentSpectrumOriginal<-NMRdataList$currentSpectrum
   }
   mrbin.env$mrbinTMP$currentSpectrumTitle<-NMRdataList$currentSpectrumTitle
   mrbin.env$mrbinTMP$currentSpectrumFolderName<-NMRdataList$currentSpectrumFolderName
   mrbin.env$mrbinTMP$currentSpectrumEXPNO<-NMRdataList$currentSpectrumEXPNO
   mrbin.env$mrbinTMP$currentSpectrumFolderName_EXPNO<-NMRdataList$currentSpectrumFolderName_EXPNO
   mrbin.env$mrbinTMP$currentSpectrumName<-NMRdataList$currentSpectrumName
}


#' A function for summing bins.
#'
#' This function sums up bins. The sums are saved to the middle (median) bin of
#' the original area. All other bins of the area are removed then. This is handy
#' for signals that are know to vary between spectra due to pH or salt content,
#' such as citric acid. Intended for use within the mrbin function.
#' @return {None}
#' @keywords internal
#' @noRd
#' @examples
#' \donttest{ sumBins() }

sumBins<-function(){#sum up regions with shifting peaks and remove remaining bins
 #if(!is.null(mrbin.env$bins)){
  if(nrow(mrbin.env$mrbinparam$sumBinList)>0&nrow(mrbin.env$mrbinTMP$binRegions)>1){
    for(i in 1:nrow(mrbin.env$mrbinparam$sumBinList)){
       limits<-mrbin.env$mrbinparam$sumBinList[i,]
       if(mrbin.env$mrbinparam$dimension=="1D"){
          #Find partially overlapping bins
          TMP_left<-mrbin.env$mrbinTMP$binRegions[,2]<limits[1]&mrbin.env$mrbinTMP$binRegions[,1]>=limits[1]
          mrbin.env$mrbinTMP$binRegions[TMP_left,2]<-limits[1]
          TMP_right<-mrbin.env$mrbinTMP$binRegions[,1]>limits[2]&mrbin.env$mrbinTMP$binRegions[,2]<=limits[2]
          mrbin.env$mrbinTMP$binRegions[TMP_right,1]<-limits[2]
          #Find completely overlapping bins
          TMP<-mrbin.env$mrbinTMP$binRegions[,2]>limits[2]&mrbin.env$mrbinTMP$binRegions[,1]<limits[1]
          if(sum(TMP)>0){
              #i_TMP<-quantile(x=1:sum(TMP), probs = .5,type=3)#define "middle" bin. This one will be kept
              #i_TMP2<-which(TMP)[i_TMP]
              #mrbin.env$mrbinTMP$binRegions[i_TMP2,1]<-limits[1]
              #mrbin.env$mrbinTMP$binRegions[i_TMP2,2]<-limits[2]
              #if(sum(TMP)>1)
              mrbin.env$mrbinTMP$binRegions<-mrbin.env$mrbinTMP$binRegions[!TMP,,drop=FALSE]
          } #else {
              #mrbin.env$mrbinTMP$binRegions<-rbind(c(0,0,0,0),mrbin.env$mrbinTMP$binRegions)
              #i_TMP2<-1#nrow(mrbin.env$mrbinTMP$binRegions)
              #mrbin.env$mrbinTMP$binRegions[i_TMP2,1]<-limits[1]
              #mrbin.env$mrbinTMP$binRegions[i_TMP2,2]<-limits[2]
          #}
          #if(!is.matrix(mrbin.env$mrbinTMP$binRegions)) mrbin.env$mrbinTMP$binRegions<-matrix(mrbin.env$mrbinTMP$binRegions,ncol=4)

       } else {#2D limits=c(4.04,4.08,58,60)
           #NMRdataNames<-cbind(apply(mrbin.env$mrbinTMP$binRegions[,3:4],1,mean),apply(mrbin.env$mrbinTMP$binRegions[,1:2],1,mean))
           #Find partially overlapping bins ("corners" are an issue here)
           TMP_left<-mrbin.env$mrbinTMP$binRegions[,2]<limits[1]&mrbin.env$mrbinTMP$binRegions[,1]>limits[1]&
                     mrbin.env$mrbinTMP$binRegions[,3]>=limits[3]&mrbin.env$mrbinTMP$binRegions[,4]<=limits[4]
           mrbin.env$mrbinTMP$binRegions[TMP_left,2]<-limits[1]
           TMP_right<-mrbin.env$mrbinTMP$binRegions[,1]>limits[2]&mrbin.env$mrbinTMP$binRegions[,2]<limits[2]&
                     mrbin.env$mrbinTMP$binRegions[,3]>=limits[3]&mrbin.env$mrbinTMP$binRegions[,4]<=limits[4]
           mrbin.env$mrbinTMP$binRegions[TMP_right,1]<-limits[2]
           TMP_top<-mrbin.env$mrbinTMP$binRegions[,1]<=limits[1]&mrbin.env$mrbinTMP$binRegions[,2]>=limits[2]&
                     mrbin.env$mrbinTMP$binRegions[,3]<limits[3]&mrbin.env$mrbinTMP$binRegions[,4]>limits[3]
           mrbin.env$mrbinTMP$binRegions[TMP_top,3]<-limits[3]
           TMP_bottom<-mrbin.env$mrbinTMP$binRegions[,1]<=limits[1]&mrbin.env$mrbinTMP$binRegions[,2]>=limits[2]&
                     mrbin.env$mrbinTMP$binRegions[,4]>limits[4]&mrbin.env$mrbinTMP$binRegions[,3]<limits[4]
           mrbin.env$mrbinTMP$binRegions[TMP_bottom,3]<-limits[4]
           #Find completely overlapping bins
           TMP<-mrbin.env$mrbinTMP$binRegions[,2]>limits[2]&mrbin.env$mrbinTMP$binRegions[,1]<limits[1]&
             mrbin.env$mrbinTMP$binRegions[,3]>limits[3]& mrbin.env$mrbinTMP$binRegions[,4]<limits[4]
           if(sum(TMP)>0){
              #i_TMP<-quantile(x=1:sum(TMP), probs = .5,type=3)#define "middle" bin. This one will be kept
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
              #mrbin.env$mrbinTMP$binRegions[which(TMP)[i_TMP],]<-limits
              #if(sum(TMP)>1)
              mrbin.env$mrbinTMP$binRegions<-mrbin.env$mrbinTMP$binRegions[!TMP,,drop=FALSE]
           } #else {
           #   mrbin.env$mrbinTMP$binRegions<-rbind(mrbin.env$mrbinTMP$binRegions,c(0,0,0,0))
           #   i_TMP2<-nrow(mrbin.env$mrbinTMP$binRegions)
           #   mrbin.env$mrbinTMP$binRegions[i_TMP2,]<-limits
           #}
          # if(!is.matrix(mrbin.env$mrbinTMP$binRegions)) mrbin.env$mrbinTMP$binRegions<-matrix(mrbin.env$mrbinTMP$binRegions,ncol=4)
       }
       #New bins are added on top of list to ensure they see all data points before they are set to NA
       mrbin.env$mrbinTMP$binRegions<-rbind(c(0,0,0,0),mrbin.env$mrbinTMP$binRegions)
       mrbin.env$mrbinTMP$binRegions[1,]<-limits
    }
  }
  mrbin.env$mrbinparam$numberOfFeaturesAfterSummingBins<-nrow(mrbin.env$mrbinTMP$binRegions)
 #}
}

#' A function for removing the solvent region from binned data.
#'
#' This function removes the solvent region. Should only be run within the mrbin function.
#' @return {None}
#' @keywords internal
#' @noRd
#' @examples
#' \donttest{ removeSolvent() }

removeSolvent<-function(){
   solventTMP<-mrbin.env$mrbinTMP$binRegions[,2]>mrbin.env$mrbinparam$solventRegion[2]&
               mrbin.env$mrbinTMP$binRegions[,1]<mrbin.env$mrbinparam$solventRegion[1]
   if(sum(solventTMP)>0){
      mrbin.env$mrbinTMP$binRegions<-mrbin.env$mrbinTMP$binRegions[-which(solventTMP),,drop=FALSE]
      if(!is.matrix(mrbin.env$mrbinTMP$binRegions)) mrbin.env$mrbinTMP$binRegions<-matrix(mrbin.env$mrbinTMP$binRegions,ncol=4)
   }
   mrbin.env$mrbinparam$numberOfFeaturesAfterRemovingSolvent<-nrow(mrbin.env$mrbinTMP$binRegions)
 #}
}


#' A function for removing additional regions from binned data.
#'
#' This function removes additional regions. This can be useful when some areas
#' are visibly affected by spectral artifacts. Should only be run from within the mrbin function.
#' @return {None}
#' @keywords internal
#' @noRd
#' @examples
#' \donttest{ removeAreas() }

removeAreas<-function(){#limits=c(4.75,4.95,-10,160)
  if(nrow(mrbin.env$mrbinparam$removeAreaList)>0){
     removeTMP<-NULL
     for(i in 1:nrow(mrbin.env$mrbinparam$removeAreaList)){
         limits<-mrbin.env$mrbinparam$removeAreaList[i,]
         if(mrbin.env$mrbinparam$dimension=="1D"){
             removeTMP2<-mrbin.env$mrbinTMP$binRegions[,2]>limits[2]&mrbin.env$mrbinTMP$binRegions[,1]<limits[1]
         }
         if(mrbin.env$mrbinparam$dimension=="2D"){
             removeTMP2<-mrbin.env$mrbinTMP$binRegions[,2]>limits[2]&mrbin.env$mrbinTMP$binRegions[,1]<limits[1]&
                mrbin.env$mrbinTMP$binRegions[,3]>limits[3]& mrbin.env$mrbinTMP$binRegions[,4]<limits[4]
         }
         if(sum(removeTMP2)>0)        removeTMP<-c(removeTMP,which(removeTMP2))
     }
     if(!is.null(removeTMP)){
         removeTMP<-unique(removeTMP)
         mrbin.env$mrbinTMP$binRegions<-mrbin.env$mrbinTMP$binRegions[-removeTMP,,drop=FALSE]
         if(!is.matrix(mrbin.env$mrbinTMP$binRegions)) mrbin.env$mrbinTMP$binRegions<-matrix(mrbin.env$mrbinTMP$binRegions,ncol=4)
     }
  }
  mrbin.env$mrbinparam$numberOfFeaturesAfterRemovingAreas<-nrow(mrbin.env$mrbinTMP$binRegions)
 #}
}


#' A function for trimming zero-values bins.
#'
#' This function removes zero-values bins. These might be created during removal of
#' solvent and additional areas, or at the edges of the spectrum. Only for internal use.
#' @return {None}
#' @keywords internal
#' @noRd
#' @examples
#' \donttest{ trimZeros() }

trimZeros<-function(){
   mrbin.env$binsTEST<-mrbin.env$bins
  if(nrow(mrbin.env$bins)>1){
    TMP<-apply(mrbin.env$bins==0,2,sum)/nrow(mrbin.env$bins)<.75
    if(sum(TMP)>0){
      mrbin.env$bins<-mrbin.env$bins[,which(TMP),drop=FALSE]
      mrbin.env$mrbinTMP$binRegions<-mrbin.env$mrbinTMP$binRegions[which(TMP),,drop=FALSE]
    }
  } else {
    #if(sum(sum(mrbin.env$bins==0)/nrow(mrbin.env$bins)<.75)>0){
    #  mrbin.env$bins<-mrbin.env$bins[,which(sum(mrbin.env$bins==0)/nrow(mrbin.env$bins)<.75),drop=FALSE]
    #}
    warning("Too few samples for zero trimming.")
  }
   mrbin.env$mrbinparam$numberOfFeaturesAfterTrimmingZeros<-ncol(mrbin.env$bins)
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
#'                     fixNegatives="No",logTrafo="No",PCA="No",verbose=TRUE,
#'                     NMRfolders=c(system.file("extdata/1/10/pdata/10",package="mrbin"),
#'                                 system.file("extdata/2/10/pdata/10",package="mrbin"),
#'                                 system.file("extdata/3/10/pdata/10",package="mrbin"))))
#' removeNoise()

removeNoise<-function(){#remove noise peaks
 if(nrow(mrbin.env$mrbinTMP$binRegions)>1){
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
          if(sum(mrbin.env$bins[,i]>(mrbin.env$mrbinTMP$noise_level[,i]*SNR))>=minimumNumber){
              colnames_NMRdata_no_noise<-c(colnames_NMRdata_no_noise,i)
          }
    }
    if(!is.null(colnames_NMRdata_no_noise)){
        #if(nrow(mrbin.env$bins)==1){
        #    rownamesTMP<-rownames(mrbin.env$bins)
        #    colnamesTMP<-colnames(mrbin.env$bins)[colnames_NMRdata_no_noise]
        #    mrbin.env$bins<-matrix(mrbin.env$bins[,colnames_NMRdata_no_noise,drop=FALSE],nrow=1)
        #    rownames(mrbin.env$bins)<-rownamesTMP
        #    colnames(mrbin.env$bins)<-colnamesTMP
        #} else {
          mrbin.env$bins<-mrbin.env$bins[,colnames_NMRdata_no_noise,drop=FALSE]
        #}
        mrbin.env$mrbinTMP$binRegions<-mrbin.env$mrbinTMP$binRegions[colnames_NMRdata_no_noise,,drop=FALSE]
        if(!is.matrix(mrbin.env$mrbinTMP$binRegions)) mrbin.env$mrbinTMP$binRegions<-matrix(mrbin.env$mrbinTMP$binRegions,ncol=4)
    } else {
        warning("No bins above noise level. Noise removal stopped.")
    }
  mrbin.env$mrbinparam$numberOfFeaturesAfterNoiseRemoval<-ncol(mrbin.env$bins)
 } else {
   warning("Too few bins for noise removal. Noise removal stopped.")
 }
}


#' A function for checking for baseline distortions.
#'
#' This function checks for each spectrum whether the median intensity in the
#' noise region is further than 10 standard deviation from zero. In this
#' case, a warning is displayed.
#' @param NMRdata Spectral data
#' @param dimension Dimension
#' @param currentSpectrumName Spectrum name
#' @param noiseRange1d Noise range
#' @param noiseRange2d Noise range
#' @return An (invisible) object containing the warning message
#' @keywords internal
#' @noRd
#' @examples
#' \donttest{ checkBaseline()  }

checkBaseline<-function(NMRdata=NULL,dimension="1D",currentSpectrumName=NULL,
   noiseRange1d=NULL,noiseRange2d=NULL){#check noise region for baseline distortions
  if(!is.null(NMRdata)){
    if(dimension=="1D"){
         baseline_level<-stats::median(NMRdata[
                 which(as.numeric(names(NMRdata))<=max(noiseRange1d[1:2])&
                      as.numeric(names(NMRdata))>=min(noiseRange1d[1:2]))])
         sd_level<-stats::sd(NMRdata[
                 which(as.numeric(names(NMRdata))<=max(noiseRange1d[1:2])&
                      as.numeric(names(NMRdata))>=min(noiseRange1d[1:2]))])
    }
    if(dimension=="2D"){
         baseline_level<-stats::median(NMRdata[
               which(as.numeric(rownames(NMRdata))>=min(noiseRange2d[3:4])&
                 as.numeric(rownames(NMRdata))<=max(noiseRange2d[3:4])),
               which(as.numeric(colnames(NMRdata))<=max(noiseRange2d[1:2])&
                 as.numeric(colnames(NMRdata))>=min(noiseRange2d[1:2]))])
         sd_level<-stats::sd(NMRdata[
               which(as.numeric(rownames(NMRdata))>=min(noiseRange2d[3:4])&
                 as.numeric(rownames(NMRdata))<=max(noiseRange2d[3:4])),
               which(as.numeric(colnames(NMRdata))<=max(noiseRange2d[1:2])&
                 as.numeric(colnames(NMRdata))>=min(noiseRange2d[1:2]))])
    }
    warningMessage<-NULL
    if(abs(baseline_level)/sd_level>10){
      warningMessage<-paste("Baseline is not flat in noise region for sample:\n",
                  currentSpectrumName,
                  "\nPlease check phase and baseline. Results may be corrupted.",
                  sep="")
    }
    invisible(warningMessage)
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
#' Example<-mrbin(silent=TRUE,
#'          parameters=list(dimension="2D",binwidth2D=1,binheight=4,cropHSQC="No",PCA="No",
#'          PQNScaling="No",noiseRemoval="No",removeSolvent="No",verbose=TRUE,
#'          NMRfolders=c(system.file("extdata/1/12/pdata/10",package="mrbin"))))
#' cropNMR()

cropNMR<-function(plot=FALSE){
if(nrow(mrbin.env$mrbinTMP$binRegions)>1){
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
  mrbin.env$mrbinTMP$binRegions<-mrbin.env$mrbinTMP$binRegions[selectedCols,,drop=FALSE]
  if(!is.matrix(mrbin.env$mrbinTMP$binRegions)) mrbin.env$mrbinTMP$binRegions<-matrix(mrbin.env$mrbinTMP$binRegions,ncol=4)
  mrbin.env$mrbinparam$numberOfFeaturesAfterCropping<-nrow(mrbin.env$mrbinTMP$binRegions)
  }
 } else {
   warning("Too few bins for cropping. Cropping stopped.")
 }

}

#' A function for PQN scaling.
#'
#' This function performs PQN scaling. To further exclude unreliable noise, only
#' the most intense signals are used. For 1H and 1H-13C HSQC spectra, most of
#' the sugar regions can be excluded to avoid a dominating effect of the
#' multiple glucose signals.
#' @param NMRdata A matrix containing NMR data. Columns=frequencies,rows=samples
#' @param ignoreGlucose A character value ("Yes" or "No")
#' @param dimension A character value ("1D" or "2D")
#' @param ppmNames A character value ("borders" or "mean")
#' @param sugarArea A numeric vector defining the the borders of glucose area
#' @param minimumFeatures A numeric value defining minimum feature number used
#' @param showHist A logical value, default is FALSE
#' @return NMRdata An invisible matrix containing scaled NMR data.
#' @export
#' @examples
#' mrbinExample<-mrbin(silent=TRUE,setDefault=TRUE,parameters=list(dimension="1D",
#'                     binwidth1D=0.05,PQNScaling="No",PCA="No",
#'                     NMRfolders=c(system.file("extdata/1/10/pdata/10",package="mrbin"),
#'                                 system.file("extdata/2/10/pdata/10",package="mrbin"),
#'                                 system.file("extdata/3/10/pdata/10",package="mrbin"))))
#' PQNScaling()

PQNScaling<-function(NMRdata=NULL,ignoreGlucose="Yes",dimension="1D",
                     ppmNames="borders",sugarArea=c(5.4,3.35,72,100),
                     minimumFeatures=40,showHist=FALSE){
  dataProvidedtoFunction<-TRUE
  if(is.null(NMRdata)){
    dataProvidedtoFunction<-FALSE
    NMRdata<-mrbin.env$bins
    ignoreGlucose<-mrbin.env$mrbinparam$PQNIgnoreSugarArea
    dimension<-mrbin.env$mrbinparam$dimension
    ppmNames<-"borders"
    sugarArea<-mrbin.env$mrbinparam$PQNsugarArea
    minimumFeatures<-mrbin.env$mrbinparam$PQNminimumFeatures
    showHist<-mrbin.env$mrbinparam$PQNshowHist
  } else {
    NMRdata<-as.matrix(NMRdata)
  }
  if(ncol(NMRdata)>1){
    if(nrow(NMRdata)>1){
      #Create synthetic median spectrum by averaging all spectra
      NMRdataTmp<-rbind(NMRdata,apply(NMRdata,2,mean))
      rownames(NMRdataTmp)[nrow(NMRdataTmp)]<-"Median"
      medianSample<-nrow(NMRdataTmp)
      #Remove most sugar signals to get a better fold change estimate
      if(ignoreGlucose=="Yes") {
        if(dimension == "2D" ) {
          selectedCols<-NULL
          if(ppmNames=="borders"){
            if(dataProvidedtoFunction){
              colnamesTMP<-matrix(as.numeric(unlist(strsplit(colnames(
                         NMRdataTmp),","))),ncol=4,byrow=TRUE)
              coordTmpAll<-cbind(apply(colnamesTMP[,3:4],1,mean),
                  apply(colnamesTMP[,1:2],1,mean))

            } else {
              coordTmpAll<-cbind(apply(mrbin.env$mrbinTMP$binRegions[,3:4],1,mean),
                  apply(mrbin.env$mrbinTMP$binRegions[,1:2],1,mean))
            }
          }
          if(ppmNames=="mean"){
              colnamesTMP<-matrix(as.numeric(unlist(strsplit(colnames(
                            NMRdataTmp),","))),ncol=2,byrow=TRUE)
              coordTmpAll<-colnamesTMP
          }
          for(j in 1:ncol(NMRdataTmp)){
            coordTmp<-coordTmpAll[j,]#1=C,2=H
            if(!(coordTmp[2]>sugarArea[2]&coordTmp[2]<sugarArea[1]&
              coordTmp[1]>sugarArea[3]&coordTmp[1]<sugarArea[4])){
                selectedCols<-c(selectedCols,j)
            }
          }
          NMRdataTmp2<-NMRdataTmp[,selectedCols,drop=FALSE]
      } else { #1D spectra
          selectedCols<-NULL
          if(ppmNames=="borders"){
            if(dataProvidedtoFunction){
              colnamesTMP<-matrix(as.numeric(unlist(strsplit(colnames(
                            NMRdataTmp),","))),ncol=2,byrow=TRUE)
              coordTmpAll<-apply(colnamesTMP[,1:2],1,mean)
            } else {
              coordTmpAll<-apply(mrbin.env$mrbinTMP$binRegions[,1:2],1,mean)
            }
          }
          if(ppmNames=="mean"){
              coordTmpAll<-as.numeric(colnames(NMRdataTmp))
          }
          for(j in 1:ncol(NMRdataTmp)){
            coordTmp<-coordTmpAll[j]#
            if(!(coordTmp>sugarArea[2]&
                 coordTmp<sugarArea[1])){
                   selectedCols<-c(selectedCols,j)
            }
          }
          NMRdataTmp2<-NMRdataTmp[,selectedCols,drop=FALSE]
      }
    }
    #Calculate fold changes versus reference sample
    NMRdataTmp_scaledMedian<-NMRdataTmp
    PQNminimumFeatures<-min(max(minimumFeatures,floor(.25*ncol(NMRdataTmp2))),ncol(NMRdataTmp2))
    medianFoldChanges<-rep(0,nrow(NMRdataTmp_scaledMedian))
    names(medianFoldChanges)<-rownames(NMRdataTmp_scaledMedian)
    if(showHist){#Plot distribution of fold changes per sample
      oldpar<-graphics::par("mar","mfrow")
      on.exit(graphics::par(oldpar))
      graphics::par(mfrow=c(ceiling(sqrt(nrow(NMRdataTmp2))),ceiling(sqrt(nrow(NMRdataTmp2)))))
      graphics::par(mar=c(.1,.1,3,.1))
    }
    for(i in 1:nrow(NMRdataTmp2)){#scale to spectral area, use only points above X% quantile for better reliability
        overlapAboveNoise<- sort(NMRdataTmp2[i,],index.return=TRUE,decreasing=TRUE)$ix[1:PQNminimumFeatures]#$ix returns the index for matrix data
        #cat(paste(PQNminimumFeatures,i,nrow(NMRdataTmp2),ncol(NMRdataTmp2),"\n"))
        #cat(is.list(NMRdataTmp2),is.list(NMRdataTmp),is.list(NMRdata),"\n")
        #cat(NMRdataTmp2[i,1:5],"\n")
        #cat(sort(NMRdataTmp2[i,],index.return=TRUE,decreasing=TRUE)$ix[1:5],"\n")
        #cat(overlapAboveNoise,"\n")
        medianFoldChanges[i]<-stats::median(NMRdataTmp2[i,overlapAboveNoise]/
                                     NMRdataTmp2[medianSample,overlapAboveNoise])
        NMRdataTmp_scaledMedian[i,]<-NMRdataTmp[i,]/medianFoldChanges[i]
      if(showHist){#Plot distribution of fold changes per sample
            graphics::hist(NMRdataTmp2[i,overlapAboveNoise]/NMRdataTmp2[medianSample,overlapAboveNoise],breaks=60,
            main=rownames(NMRdataTmp2)[i],xlab="",ylab="")
            graphics::lines(rep(medianFoldChanges[i],2),c(0,20000),col='red')
      }
    }
    #Remove reference sample from list
    NMRdataTmp_scaledMedian<-NMRdataTmp_scaledMedian[-nrow(NMRdataTmp_scaledMedian),,drop=FALSE]
    if(!dataProvidedtoFunction){
      mrbin.env$bins<-NMRdataTmp_scaledMedian
      mrbin.env$mrbinparam$medians<-medianFoldChanges
    }
   } else {
      warning("Too few samples to perform PQN normalization.")
   }
 } else {
   warning("Too few bins to perform PQN scaling.")
 }
 invisible(NMRdataTmp_scaledMedian)
}


#' A function for plotting quality indicators, including PCA plots.
#'
#' This function plots boxplots (bin-wise and sample-wise) as visual quality indicators. It also performs PCA, then plots PC1 and PC2 and loading plots.
#' @return {None}
#' @export
#' @examples
#' mrbinExample<-mrbin(silent=TRUE,setDefault=FALSE,parameters=list(dimension="2D",
#'                     binwidth2D=0.05,binheight=3,
#'                     PQNScaling="No",noiseRemoval="Yes",trimZeros="Yes",
#'                     fixNegatives="No",logTrafo="No",PCA="No",
#'                     NMRfolders=c(system.file("extdata/1/12/pdata/10",package="mrbin"),
#'                                 system.file("extdata/2/12/pdata/10",package="mrbin"),
#'                                 system.file("extdata/3/12/pdata/10",package="mrbin"))))
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
    if(mrbin.env$mrbinparam$dimension=="1D"){
      plotNMR(rectangleRegions=mrbin.env$mrbinTMP$binRegions,
          color="black",rectangleColors="green", rectangleFront=FALSE,
          manualScale=FALSE,
          plotTitle="")
    }
    if(mrbin.env$mrbinparam$dimension=="2D"){
      plotNMR(rectangleRegions=mrbin.env$mrbinTMP$binRegions,
          color="darkgray",rectangleColors="blue", rectangleFront=TRUE,
          manualScale=FALSE,
          plotTitle="")
    }
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
      numlevels<-NULL
      if(mrbin.env$mrbinparam$defineGroups=="Yes"){
        for(i in 1:nlevels(FactorsTMP)) numlevels<-c(numlevels,as.numeric(
                         FactorsTMP[which(FactorsTMP==levels(FactorsTMP)[i])][1]))
      }
      addTMP1<-.15*(max(mrbin.env$mrbinTMP$PCA$x[,1])-min(mrbin.env$mrbinTMP$PCA$x[,1]))
      xlimTMP<-c(min(mrbin.env$mrbinTMP$PCA$x[,1])-addTMP1,max(mrbin.env$mrbinTMP$PCA$x[,1])+addTMP1)
      addTMP1rot<-.25*(max(mrbin.env$mrbinTMP$PCA$rotation[,1])-min(mrbin.env$mrbinTMP$PCA$rotation[,1]))
      xlimTMProt<-c(min(mrbin.env$mrbinTMP$PCA$rotation[,1])-addTMP1rot,max(mrbin.env$mrbinTMP$PCA$rotation[,1])+addTMP1rot)
      xlabTMP<-paste("PC1 (",round(100*(mrbin.env$mrbinTMP$PCA$sdev[1]^2)/sum(mrbin.env$mrbinTMP$PCA$sdev^2),1),"%)",sep="")
      if(ncol(mrbin.env$mrbinTMP$PCA$x)>1){
        addTMP2<-.15*(max(mrbin.env$mrbinTMP$PCA$x[,2])-min(mrbin.env$mrbinTMP$PCA$x[,2]))
        ylimTMP<-c(min(mrbin.env$mrbinTMP$PCA$x[,2])-addTMP2,max(mrbin.env$mrbinTMP$PCA$x[,2])+addTMP2)
        ylabTMP<-paste("PC2 (",round(100*(mrbin.env$mrbinTMP$PCA$sdev[2]^2)/sum(mrbin.env$mrbinTMP$PCA$sdev^2),1),"%)",sep="")
        addTMP2rot<-.25*(max(mrbin.env$mrbinTMP$PCA$rotation[,2])-min(mrbin.env$mrbinTMP$PCA$rotation[,2]))
        ylimTMProt<-c(min(mrbin.env$mrbinTMP$PCA$rotation[,2])-addTMP2rot,max(mrbin.env$mrbinTMP$PCA$rotation[,2])+addTMP2rot)
      } else {
        ylimTMP<-NULL
        ylabTMP<-""
        ylimTMProt<-NULL
      }
      if(mrbin.env$mrbinparam$defineGroups=="Yes"){
        PCAFactors<-as.numeric(FactorsTMP)
      } else {
        PCAFactors<-rep(1,nrow(mrbin.env$bins))
      }
      graphics::plot(mrbin.env$mrbinTMP$PCA$rotation,xlim=xlimTMProt,ylim=ylimTMProt,xlab="PC1",ylab="PC2",
                     pch=16,cex=.75,main="PCA Loadings Plot", ask=FALSE, xaxt='n', yaxt='n')
      graphics::text(mrbin.env$mrbinTMP$PCA$rotation,labels=#substr(
                     rownames(mrbin.env$mrbinTMP$PCA$rotation)#,1,12)
                     ,pos=4,cex=1)
      graphics::plot(mrbin.env$mrbinTMP$PCA$x#[,1],mrbin.env$mrbinTMP$PCA$x[,2],
           ,xlim=xlimTMP,
           ylim=ylimTMP,
           pch=PCAFactors+14, xaxt='n', yaxt='n',
           col=colorPalette[PCAFactors],
           main="PCA",
           xlab=xlabTMP,
           ylab=ylabTMP,
           cex=.75, ask=FALSE)
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
#' @param rectangleFront Plot rectangles in front of spectrum rather than in background (only 2D)
#' @param polygonRegion Defines 4 corners of a polygon to be plotted
#' @param color Defines the color of the spectrum plot. If NULL, a rainbow theme is used for 2D NMR
#' @param add If TRUE, additional spectrum plots are overlaid with the current plot
#' @param manualScale If TRUE, scaling factor is taken from environment variables
#' @param plotTitle Defines the main title of the plot
#' @param showGrid Shows a grid of data points. Defaults to FALSE
#' @return {None}
#' @export
#' @examples
#' mrbin(silent=TRUE,setDefault=TRUE,parameters=list(dimension="1D",binwidth1D=.1,
#'          PQNScaling="No",noiseRemoval="No",trimZeros="No",
#'          fixNegatives="No",logTrafo="No",PCA="No",verbose=TRUE,
#'          NMRfolders=system.file("extdata/1/10/pdata/10",package="mrbin")))
#' plotNMR()

plotNMR<-function(region=NULL,rectangleRegions=NULL,
                   rectangleColors=c("green","orange","blue","red","yellow","gray","purple"),
                   rectangleFront=FALSE,
                   polygonRegion=NULL,
                   color=NULL,add=FALSE,showGrid=FALSE,
                   manualScale=TRUE,plotTitle=""){
 if(!is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
   devAskNewPage(ask = FALSE)
   if(is.null(region)){
       if(is.null(mrbin.env$mrbinplot$plotRegion)){
            if(is.matrix(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
                mrbin.env$mrbinplot$plotRegion<-
                          c(max(as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrumOriginal))),
                            min(as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrumOriginal))),
                            min(as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrumOriginal))),
                            max(as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrumOriginal))))
            } else {
                mrbin.env$mrbinplot$plotRegion<-
                         c(max(as.numeric(names(mrbin.env$mrbinTMP$currentSpectrumOriginal))),
                         min(as.numeric(names(mrbin.env$mrbinTMP$currentSpectrumOriginal))),
                         -10,160)
            }
       }
       region<-mrbin.env$mrbinplot$plotRegion
   }
   if(length(region)==1){
     if(region=="all"){
          if(is.matrix(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
                region<-c(max(as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrumOriginal))),
                            min(as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrumOriginal))),
                            min(as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrumOriginal))),
                            max(as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrumOriginal))))
          } else {
                region<-c(max(as.numeric(names(mrbin.env$mrbinTMP$currentSpectrumOriginal))),
                         min(as.numeric(names(mrbin.env$mrbinTMP$currentSpectrumOriginal))),
                         -10,160)
          }
     }
   }
   if(is.matrix(mrbin.env$mrbinTMP$currentSpectrumOriginal)){#2D spectra
      if(is.null(color)) color<-rev(grDevices::rainbow(mrbin.env$mrbinplot$nContours))
      spectrumTMP<-mrbin.env$mrbinTMP$currentSpectrumOriginal[which(as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrumOriginal))<
                           region[4]&
                          as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrumOriginal))>=region[3]),
                         which(as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrumOriginal))>=region[2]&
                         as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrumOriginal))<region[1])
                         ]
      if(manualScale){
        if(sum(spectrumTMP<(mrbin.env$mrbinplot$lowestContour*max(mrbin.env$mrbinTMP$currentSpectrumOriginal)))>0){
             spectrumTMP[spectrumTMP<(mrbin.env$mrbinplot$lowestContour*max(mrbin.env$mrbinTMP$currentSpectrumOriginal))]<-0
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
      if(!is.null(rectangleRegions)&!rectangleFront){
          graphics::rect(xleft=-rectangleRegions[,1], ybottom=-rectangleRegions[,4],
                         xright=-rectangleRegions[,2], ytop=-rectangleRegions[,3],
                         col = rectangleColors, border = NA)
          graphics::box()
      }
      if(showGrid){
          GridxCoordTMP<-as.numeric(colnames(spectrumTMP))
          GridxCoordTMP<-GridxCoordTMP[GridxCoordTMP<region[1]&GridxCoordTMP>region[2]]
          GridxCoordTMP<--GridxCoordTMP
          GridyCoordTMP<-as.numeric(rownames(spectrumTMP))
          GridyCoordTMP<-GridyCoordTMP[GridyCoordTMP<region[4]&GridyCoordTMP>region[3]]
          GridyCoordTMP<--GridyCoordTMP
          GridxCoord<-rep(GridxCoordTMP,length(GridyCoordTMP))
          GridyCoord<-sort(rep(GridyCoordTMP,length(GridxCoordTMP)))
          graphics::points(GridxCoord,GridyCoord,
                         col = "darkgray",pch="+",cex=.6)
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
      if(!is.null(rectangleRegions)&rectangleFront){
          #graphics::rect(xleft=-rectangleRegions[,1], ybottom=-rectangleRegions[,4],
          #               xright=-rectangleRegions[,2], ytop=-rectangleRegions[,3],
          #               col = rectangleColors, border = NA)
          graphics::rect(xleft=-rectangleRegions[,1], ybottom=-rectangleRegions[,4],
                         xright=-rectangleRegions[,2], ytop=-rectangleRegions[,3],
                         col = rectangleColors, border = rectangleColors)
          graphics::box()
      }
      utils::flush.console()
   } else {  #1D
      if(is.null(color)) color<-"black"
      if(manualScale){
        ymin<-min(mrbin.env$mrbinTMP$currentSpectrumOriginal)
        ymax<-max(mrbin.env$mrbinTMP$currentSpectrumOriginal)/mrbin.env$mrbinplot$intensityScale
      } else {
        ymin<-0
        spectrumTMP<-mrbin.env$mrbinTMP$currentSpectrumOriginal[as.numeric(names(mrbin.env$mrbinTMP$currentSpectrumOriginal<region[1]))&
                   as.numeric(names(mrbin.env$mrbinTMP$currentSpectrumOriginal))>region[2]]
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
      graphics::lines(as.numeric(names(mrbin.env$mrbinTMP$currentSpectrumOriginal)),mrbin.env$mrbinTMP$currentSpectrumOriginal,
                col=color)
      if(showGrid){
          GridxCoord<-as.numeric(names(mrbin.env$mrbinTMP$currentSpectrumOriginal))
          Gridx<-mrbin.env$mrbinTMP$currentSpectrumOriginal[GridxCoord<region[1]&GridxCoord>region[2]]
          GridxCoord<-GridxCoord[GridxCoord<region[1]&GridxCoord>region[2]]
          graphics::points(GridxCoord,Gridx,
                         col = "darkgray",pch="+",cex=.85)
      }
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
#'          PQNScaling="No",noiseRemoval="No",trimZeros="No",
#'          fixNegatives="No",logTrafo="No",PCA="No",verbose=TRUE,
#'          NMRfolders=system.file("extdata/1/10/pdata/10",package="mrbin")))
#' plotNMR()
#' intPlus()

intPlus<-function(){#increase plot intensity
 if(!is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
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
#'          PQNScaling="No",noiseRemoval="No",trimZeros="No",
#'          fixNegatives="No",logTrafo="No",PCA="No",verbose=TRUE,
#'          NMRfolders=system.file("extdata/1/10/pdata/10",package="mrbin")))
#' plotNMR()
#' intMin()

intMin<-function(){#decrease plot intensity
 if(!is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
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
#' readBruker(folder=system.file("extdata/1/12/pdata/10",package="mrbin"),dimension="2D")
#' plotNMR()
#' contPlus()

contPlus<-function(){#decrease plot intensity
 if(!is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
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
#'          binheight=3,PQNScaling="No",referenceScaling="No",binRegion=c(4,3,60,65),
#'          noiseRemoval="No",trimZeros="No",cropHSQC="No",
#'          fixNegatives="No",logTrafo="No",PCA="No",verbose=TRUE,saveFiles="No",
#'          NMRfolders=c(system.file("extdata/1/12/pdata/10",package="mrbin"))))
#' plotNMR()
#' contMin()

contMin<-function(){#decrease plot intensity
 if(!is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
   mrbin.env$mrbinplot$lowestContour<-mrbin.env$mrbinplot$lowestContour*0.75
   plotNMR()
 }
}

#' A function for changing plotNMR plots.
#'
#' This function changes the plot region of the current NMR plot. Can be called with
#' no arguments: zoom(). In this case the user will be asked for manual input.
#' @param left New left boundary
#' @param right New right boundary
#' @param top New top boundary
#' @param bottom New bottom boundary
#' @return {None}
#' @export
#' @examples
#' mrbin(silent=TRUE,setDefault=TRUE,parameters=list(dimension="1D",binwidth1D=.1,
#'          PQNScaling="No",noiseRemoval="No",trimZeros="No",
#'          fixNegatives="No",logTrafo="No",PCA="No",verbose=TRUE,
#'          NMRfolders=system.file("extdata/1/10/pdata/10",package="mrbin")))
#' plotNMR()
#' zoom(left=4.6,right=2,top=10,bottom=150)

zoom<-function(left=NULL,right=NULL,top=NULL,bottom=NULL){
 if(!is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
   if(is.null(left)) stop("Please set left limit\n")
   if(is.null(right)) stop("Please set right limit\n")
   if(is.matrix(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
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
#'          PQNScaling="No",noiseRemoval="No",trimZeros="No",
#'          fixNegatives="No",logTrafo="No",PCA="No",verbose=TRUE,
#'          NMRfolders=system.file("extdata/1/10/pdata/10",package="mrbin")))
#' plotNMR()
#' zoomIn()

zoomIn<-function(){#Zoom into NMR spectrum plot
 if(!is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
   if(mrbin.env$mrbinparam$dimension=="1D"){
       leftMax<-max(as.numeric(names(mrbin.env$mrbinTMP$currentSpectrumOriginal)))
       rightMax<-min(as.numeric(names(mrbin.env$mrbinTMP$currentSpectrumOriginal)))
       topMax<--10
       bottomMax<-160
   }
   if(mrbin.env$mrbinparam$dimension=="2D"){
       leftMax<-max(as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrumOriginal)))
       rightMax<-min(as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrumOriginal)))
       topMax<-min(as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrumOriginal)))
       bottomMax<-max(as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrumOriginal)))
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
#'          PQNScaling="No",noiseRemoval="No",trimZeros="No",
#'          fixNegatives="No",logTrafo="No",PCA="No",verbose=TRUE,
#'          NMRfolders=system.file("extdata/1/10/pdata/10",package="mrbin")))
#' plotNMR()
#' zoomIn()
#' zoomOut()

zoomOut<-function(){#Zoom out from NMR spectrum plot
 if(!is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
   if(mrbin.env$mrbinparam$dimension=="1D"){
       leftMax<-max(as.numeric(names(mrbin.env$mrbinTMP$currentSpectrumOriginal)))
       rightMax<-min(as.numeric(names(mrbin.env$mrbinTMP$currentSpectrumOriginal)))
       topMax<--10
       bottomMax<-160
   }
   if(mrbin.env$mrbinparam$dimension=="2D"){
       leftMax<-max(as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrumOriginal)))
       rightMax<-min(as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrumOriginal)))
       topMax<-min(as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrumOriginal)))
       bottomMax<-max(as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrumOriginal)))
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
#'          noiseRemoval="No",trimZeros="No",
#'          PQNScaling="No",saveFiles="No",referenceScaling="No",
#'          fixNegatives="No",logTrafo="No",PCA="No",verbose=TRUE,
#'          NMRfolders=system.file("extdata/1/10/pdata/10",package="mrbin")))
#' plotNMR()
#' zoomIn()
#' left()

left<-function(){
 if(!is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
   if(mrbin.env$mrbinparam$dimension=="1D"){
       leftMax<-max(as.numeric(names(mrbin.env$mrbinTMP$currentSpectrumOriginal)))
       rightMax<-min(as.numeric(names(mrbin.env$mrbinTMP$currentSpectrumOriginal)))
   }
   if(mrbin.env$mrbinparam$dimension=="2D"){
       leftMax<-max(as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrumOriginal)))
       rightMax<-min(as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrumOriginal)))
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
#'          noiseRemoval="No",trimZeros="No",
#'          PQNScaling="No",saveFiles="No",referenceScaling="No",
#'          fixNegatives="No",logTrafo="No",PCA="No",verbose=TRUE,
#'          NMRfolders=system.file("extdata/1/10/pdata/10",package="mrbin")))
#' plotNMR()
#' zoomIn()
#' right()

right<-function(){
 if(!is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
   if(mrbin.env$mrbinparam$dimension=="1D"){
       leftMax<-max(as.numeric(names(mrbin.env$mrbinTMP$currentSpectrumOriginal)))
       rightMax<-min(as.numeric(names(mrbin.env$mrbinTMP$currentSpectrumOriginal)))
   }
   if(mrbin.env$mrbinparam$dimension=="2D"){
       leftMax<-max(as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrumOriginal)))
       rightMax<-min(as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrumOriginal)))
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
#'          binheight=3,PQNScaling="No",referenceScaling="No",binRegion=c(4,3,60,65),
#'          noiseRemoval="No",trimZeros="No",cropHSQC="No",
#'          fixNegatives="No",logTrafo="No",PCA="No",verbose=TRUE,saveFiles="No",
#'          NMRfolders=c(system.file("extdata/1/12/pdata/10",package="mrbin"))))
#' plotNMR()
#' zoomIn()
#' down()

down<-function(){
 if(!is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
   if(mrbin.env$mrbinparam$dimension=="2D"){
       topMax<-min(as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrumOriginal)))
       bottomMax<-max(as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrumOriginal)))
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
#'          binheight=3,PQNScaling="No",referenceScaling="No",binRegion=c(4,3,60,65),
#'          noiseRemoval="No",trimZeros="No",cropHSQC="No",
#'          fixNegatives="No",logTrafo="No",PCA="No",verbose=TRUE,saveFiles="No",
#'          NMRfolders=c(system.file("extdata/1/12/pdata/10",package="mrbin"))))
#' plotNMR()
#' zoomIn()
#' up()

up<-function(){
 if(!is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
   if(mrbin.env$mrbinparam$dimension=="2D"){
       topMax<-min(as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrumOriginal)))
       bottomMax<-max(as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrumOriginal)))
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

#' A function for binning NMR spectra.
#'
#' This function performs binning of a selected spectrum. This function is
#' meant only for use within the mrbin function.
#' @param folder Defines the exact NMR data folder. If NULL, mrbin parameter set is used
#' @param dimension Dimension
#' @param binRegions Bin regions
#' @param referenceScaling "Yes" or "No"
#' @param removeSolvent "Yes" or "No"
#' @param removeAreas "Yes" or "No"
#' @param reference1D Reference region
#' @param reference2D Reference region
#' @param solventRegion Solvent Region
#' @param removeAreaList Regions to be removed
#' @param NMRvendor Defines the NMR manufacturer, default is "Bruker"
#' @param noiseRange1d Noise range
#' @param noiseRange2d Noise range
#' @param binMethod Binning method
#' @param useAsNames How should sample names be generated
#' @return An (invisible) list containing binned data and related data
#' @keywords internal
#' @noRd
#' @examples
#' \donttest{ binMultiNMR2() }

binMultiNMR2<-function(folder=NULL,dimension="1D",
   binRegions=NULL,referenceScaling="No",removeSolvent="No",
   removeAreas="No",reference1D=NULL,reference2D=NULL,solventRegion=NULL,
   removeAreaList=NULL,NMRvendor="Bruker",noiseRange1d=NULL,
   noiseRange2d=NULL,binMethod="Rectangular bins",
   useAsNames="Spectrum titles"
   ){#Bin NMR spectral data
  if(!is.null(folder)){
    warningMessage<-NULL
    NMRdataList<-readNMR(folder=folder,dimension=dimension,
                  NMRvendor=NMRvendor,useAsNames=useAsNames)
    NMRdata<-NMRdataList$currentSpectrum
    NMRdataOriginal<-NMRdata
    scalingFactor<-1
    if(referenceScaling=="Yes"){
      referenceScalingList<-referenceScaling(NMRdata=NMRdata,
          reference1D=reference1D,reference2D=reference2D,dimension=dimension)
      NMRdata<-referenceScalingList$scaledSpectrum
      scalingFactor<-referenceScalingList$scalingFactor
      if(scalingFactor<0){
        scalingFactor<-abs(scalingFactor)
        warningMessage<-#c(warningMessage,
          paste(
          "Reference signal is negative for sample:\n",
          NMRdataList$currentSpectrumName,
          "\nPlease check phase and baseline. Results may be corrupted.",
          sep="")#)
        #warning(warningMessage)
      }
    }
    if(removeSolvent=="Yes") NMRdata<-removeSolvent2(NMRdata=NMRdata,
               dimension=dimension,solventRegion=solventRegion)
    if(removeAreas=="Yes") NMRdata<-removeAreas2(NMRdata=NMRdata,
               dimension=dimension,removeAreaList=removeAreaList)
    #mrbin.env$mrbinTMP$binTMP<-NULL
    binData<-binSingleNMR(currentSpectrum=NMRdata,dimension=dimension,
              binRegions=binRegions,binMethod=binMethod)
    noiseData<-calculateNoise(NMRdata=NMRdataOriginal,
               pointsPerBin=binData$pointsPerBin,dimension=dimension,
               noiseRange1d=noiseRange1d,noiseRange2d=noiseRange2d
               )#list
    if(referenceScaling=="Yes"){
      if(scalingFactor<(3*noiseData$noise_level)){
        warningMessage<-c(warningMessage,paste(
          "Reference signal is very low for sample:\n",
          NMRdataList$currentSpectrumName,
          "\nPlease check if reference peak is at 0ppm. Results may be corrupted.",
          sep=""))
      }
    }
    warningMessageTMP<-
       checkBaseline(NMRdata=NMRdataOriginal,dimension=dimension,
               currentSpectrumName=NMRdataList$currentSpectrumName,
               noiseRange1d=noiseRange1d,noiseRange2d=noiseRange2d)
    if(!is.null(warningMessageTMP)){
      warningMessage<-paste(warningMessage,warningMessageTMP,sep="\n")
    }
    #noiseDataScaled<-calculateNoise(NMRdata=NMRdata,
    #           pointsPerBin=binData$pointsPerBin,dimension=dimension,
    #           noiseRange1d=noiseRange1d,noiseRange2d=noiseRange2d
    #           )#list
    invisible(list(binTMP=binData$binTMP,
       meanNumberOfPointsPerBin_TMP=binData$pointsPerBin,#mrbin.env$mrbinTMP$meanNumberOfPointsPerBin_TMP
       noise_level_Raw_TMP=noiseData$noise_level,
       noise_level_TMP=noiseData$noise_level_TMP/scalingFactor,#noiseDataScaled$noise_level_TMP,
       currentSpectrumName=NMRdataList$currentSpectrumName,
       warningMessage=warningMessage
       ))
  }
}


#' A function for reading NMR spectra.
#'
#' This function picks the correct NMR reading function, based on vendor. This function is
#' meant only for use within the mrbin function.
#' @param folder Defines the exact NMR data folder. If NULL, mrbin parameter set is used
#' @param dimension Defines the data dimension, "1D" or "2D". Only used if not NULL
#' @param onlyTitles Read only spectrum titles, but no data. Defaults to FALSE
#' @param NMRvendor Defines the NMR manufacturer, default is "Bruker"
#' @param useAsNames How should sample names be generated
#' @return An (invisible) list containing spectral data and the spectrum name
#' @keywords internal
#' @noRd
#' @examples
#' \donttest{ readNMR() }

readNMR<-function(folder=NULL,dimension=NULL,onlyTitles=FALSE,
          NMRvendor="Bruker",useAsNames="Spectrum titles"){#Read NMR spectral data
 #if(!is.null(mrbin.env$mrbinTMP$currentFolder)){
  if(NMRvendor=="Bruker"){
      currentSpectrum<-readBruker(folder=folder,dimension=dimension,
                      onlyTitles=onlyTitles,useAsNames=useAsNames)
  }  else {
      stop(paste("No data import function defined for vendor ",NMRvendor,".\n",sep=""))
  }
  invisible(currentSpectrum)
 #}
}

#' A function for reading Bruker NMR spectra.
#'
#' This function reads Bruker NMR data. 1D and 2D data are supported.
#' @param folder Defines the exact NMR data folder. If NULL, mrbin parameter set is used
#' @param dimension Defines the data dimension, "1D" or "2D". Only used if not NULL
#' @param onlyTitles Read only spectrum titles, but no data. Defaults to FALSE
#' @param useAsNames How should sample names be generated
#' @return An (invisible) list containing spectral data and the spectrum name
#' @export
#' @examples
#' exampleData<-readBruker(folder=system.file("extdata/1/10/pdata/10",package="mrbin"),
#'                         dimension="1D")

readBruker<-function(folder=NULL,dimension=NULL,onlyTitles=FALSE,
  useAsNames="Spectrum titles"){#Read Bruker NMR spectral data
 datanameDict<-c("1r","2rr")
 names(datanameDict)<-c("1D","2D")
 #if(is.null(folder)){
 #  spectrum_proc_path<-gsub('\\\\',"/",mrbin.env$mrbinTMP$currentFolder)
 #} else {
   spectrum_proc_path<-folder
 #}
 #if(is.null(dimension)){
 #  datanameTmp<-datanameDict[mrbin.env$mrbinparam$dimension]
 #} else {
   datanameTmp<-datanameDict[dimension]
 #}
 if(!is.null(spectrum_proc_path)){
   BYTORDP_Dict<-c("little","big")
   names(BYTORDP_Dict)<-c(0,1)
   TITLE<-scan(file=paste(spectrum_proc_path,"/title",sep=""),what="character",sep="\n",quiet=TRUE)[1]
   if(!onlyTitles){
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
   } else {
     currentSpectrum <- NULL
   }
   currentSpectrumTitle<-TITLE
   currentSpectrumFolderName<-rev(strsplit(spectrum_proc_path,"/")[[1]])[4]
   currentSpectrumEXPNO<-rev(strsplit(spectrum_proc_path,"/")[[1]])[3]
   currentSpectrumFolderName_EXPNO<-paste(
               currentSpectrumFolderName,
               paste(c("_","0","0","0","0")[1:max(1,5-nchar(currentSpectrumEXPNO))],sep="",collapse=""),
               currentSpectrumEXPNO,sep="")
   if(useAsNames=="Spectrum titles")    titleFinal<-currentSpectrumTitle
   if(useAsNames=="Folder names")    titleFinal<-currentSpectrumFolderName
   if(useAsNames=="Folder names and EXPNO")    titleFinal<-currentSpectrumFolderName_EXPNO
   #if(is.null(folder)){
   #  if(!onlyTitles){
   #    mrbin.env$mrbinTMP$currentSpectrum<-currentSpectrum
   #    mrbin.env$mrbinTMP$currentSpectrumOriginal<-currentSpectrum
   #  }
   #  mrbin.env$mrbinTMP$currentSpectrumTitle<-currentSpectrumTitle
   #  mrbin.env$mrbinTMP$currentSpectrumFolderName<-currentSpectrumFolderName
   #  mrbin.env$mrbinTMP$currentSpectrumEXPNO<-currentSpectrumEXPNO
   #  mrbin.env$mrbinTMP$currentSpectrumFolderName_EXPNO<-currentSpectrumFolderName_EXPNO
   #  mrbin.env$mrbinTMP$currentSpectrumName<-titleFinal
   #}
   invisible(list(currentSpectrum=currentSpectrum,currentSpectrumName=titleFinal,
           currentSpectrumTitle=currentSpectrumTitle,
           currentSpectrumFolderName=currentSpectrumFolderName,
           currentSpectrumEXPNO=currentSpectrumEXPNO,
           currentSpectrumFolderName_EXPNO=currentSpectrumFolderName_EXPNO))
 }
}

#' A function for scaling to the reference area.
#'
#' This function scales NMR data to the reference area.
#' @param NMRdata Spectral data
#' @param reference1D Reference region
#' @param reference2D Reference region
#' @param dimension Dimension
#' @return An (invisible) list of spectral data and the scaling factor
#' @keywords internal
#' @noRd
#' @examples
#' \donttest{
#' resetEnv()
#' mrbin(silent=TRUE,setDefault=TRUE,parameters=list(dimension="1D",
#'          referenceScaling="No",binwidth1D=0.05,PQNScaling="No",PCA="No",
#'          NMRfolders=system.file("extdata/1/10/pdata/10",package="mrbin")))
#' referenceScaling()
#' }

referenceScaling<-function(NMRdata=NULL,reference1D=NULL,reference2D=NULL,dimension="1D"){
  dataProvided<-TRUE
  #if(is.null(NMRdata)){
  #  NMRdata<-mrbin.env$mrbinTMP$currentSpectrum
  #  dataProvided<-FALSE
  #}
  if(dimension=="1D"){
       scalingFactor<-mean(NMRdata[which(as.numeric(names(NMRdata))<=
                           reference1D[1]&as.numeric(names(NMRdata))>
                                       reference1D[2])],na.rm=TRUE)
       scaledSpectrum<-NMRdata/scalingFactor
  }
  if(dimension=="2D"){
       selectedTMP1<-which(as.numeric(rownames(NMRdata))<=
                       reference2D[4]&as.numeric(rownames(
                       NMRdata))>reference2D[3])
       selectedTMP2<-which(as.numeric(colnames(NMRdata))<=
                       reference2D[1]&
                       as.numeric(colnames(NMRdata))>
                       reference2D[2])
       scalingFactor<-mean(NMRdata[selectedTMP1,selectedTMP2],na.rm=TRUE)
       scaledSpectrum<-NMRdata/scalingFactor
  }
  #if(!dataProvided){
  #  mrbin.env$mrbinTMP$currentSpectrum<-scaledSpectrum
  #}
  invisible(list(scaledSpectrum=scaledSpectrum,scalingFactor=scalingFactor))
}

#' A function for removing the solvent region from raw data.
#'
#' This function removes the solvent region. Should only be run within the mrbin function.
#' @param NMRdata Spectral data
#' @param dimension Dimension
#' @param solventRegion Solvent Region
#' @return An (invisible) object containing spectral data
#' @keywords internal
#' @noRd
#' @examples
#' \donttest{ removeSolvent2() }

removeSolvent2<-function(NMRdata=NULL,dimension="1D",solventRegion=NULL){
  #dataProvided<-TRUE
  #if(is.null(NMRdata)){
  #  NMRdata<-mrbin.env$mrbinTMP$currentSpectrum
  #  dataProvided<-FALSE
  #}
   if(dimension=="1D"){
       newNMRdata<-NMRdata[which(
                      as.numeric(names(NMRdata))>solventRegion[1]|
                      as.numeric(names(NMRdata))<solventRegion[2]
                      )]
   }
   if(dimension=="2D"){
       newNMRdata<-NMRdata[,which(
                      as.numeric(colnames(NMRdata))>solventRegion[1]|
                      as.numeric(colnames(NMRdata))<solventRegion[2]
                      )]
   }
  #if(!dataProvided){
  #  mrbin.env$mrbinTMP$currentSpectrum<-newNMRdata
  #}
  invisible(newNMRdata)
}

#' A function for removing additional regions from raw data.
#'
#' This function removes additional regions. This can be useful when some areas
#' are visibly affected by spectral artifacts. Should only be run from within the mrbin function.
#' @param NMRdata Spectral data
#' @param dimension Dimension
#' @param removeAreaList Regions to be removed
#' @return An (invisible) object containing spectral data
#' @keywords internal
#' @noRd
#' @examples
#' \donttest{ removeAreas2() }

removeAreas2<-function(NMRdata=NULL,dimension="1D",removeAreaList=NULL){
  if(nrow(removeAreaList)>0){
    #dataProvided<-TRUE
    #if(is.null(NMRdata)){
    #  NMRdata<-mrbin.env$mrbinTMP$currentSpectrum
    #  dataProvided<-FALSE
    #}
     for(i in 1:nrow(removeAreaList)){
       if(dimension=="1D"){
           newNMRdata<-NMRdata[which(
                          as.numeric(names(NMRdata))>removeAreaList[i,1]|
                          as.numeric(names(NMRdata))<removeAreaList[i,2]
                          )]
       }
       if(dimension=="2D"){
           newNMRdata<-NMRdata
           newNMRdata[which(
                          as.numeric(rownames(newNMRdata))<=removeAreaList[i,4]&
                          as.numeric(rownames(newNMRdata))>=removeAreaList[i,3]
                          ),which(
                          as.numeric(colnames(newNMRdata))<=removeAreaList[i,1]&
                          as.numeric(colnames(newNMRdata))>=removeAreaList[i,2]
                          )]<-NA
           #Does not work - Cannot remove a few data points from the middle of a matrix
           #mrbin.env$mrbinTMP$currentSpectrum<-mrbin.env$mrbinTMP$currentSpectrum[which(
           #               as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrum))>mrbin.env$mrbinparam$removeAreaList[i,4]|
           #               as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrum))<mrbin.env$mrbinparam$removeAreaList[i,3]
           #               ),which(
           #               as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrum))>mrbin.env$mrbinparam$removeAreaList[i,1]|
           #               as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrum))<mrbin.env$mrbinparam$removeAreaList[i,2]
           #               )]
       }
     }
    #if(!dataProvided){
    #  mrbin.env$mrbinTMP$currentSpectrum<-newNMRdata
    #}
    invisible(newNMRdata)
  }
}

#' A function for binning a single NMR spectrum.
#'
#' This function creates bins for the current spectrum. This function is
#' meant only for use within the mrbin function.
#' @param NMRdata Spectral data
#' @param dimension Dimension
#' @param binRegions Bin regions
#' @param binMethod Binning method
#' @return An (invisible) list containing binned data
#' @keywords internal
#' @noRd
#' @examples
#' \donttest{ binSingleNMR() }

binSingleNMR<-function(currentSpectrum=NULL,dimension="1D",
  binRegions=NULL,binMethod="Rectangular bins"){
  #dataProvided<-TRUE
  #if(is.null(currentSpectrum)){
  #  currentSpectrum<-mrbin.env$mrbinTMP$currentSpectrum
  #  dataProvided<-FALSE
  #}
  warningFlag<-FALSE
  numberOfPointsPerBin<-NULL
  binTMP<-rep(0,nrow(binRegions))
  names(binTMP)<-rownames(binRegions)
  if(dimension=="2D"){#2d spectra
    #Create index of signals in each bin
    NMRspectrumRownames<-as.numeric(rownames(currentSpectrum))
    NMRspectrumColnames<-as.numeric(colnames(currentSpectrum))
    #counter<-1
    for(ibinTMP in 1:nrow(binRegions)){
      rowsTMP<-NMRspectrumRownames<=binRegions[ibinTMP,4]&
                            NMRspectrumRownames>binRegions[ibinTMP,3]
      colsTMP<-NMRspectrumColnames<=binRegions[ibinTMP,1]&
                    NMRspectrumColnames>binRegions[ibinTMP,2]
      numberOfPointsPerBinTMP<-(sum(rowsTMP)*sum(colsTMP))-sum(is.na(
        currentSpectrum[rowsTMP,colsTMP]))
      numberOfPointsPerBin<-c(numberOfPointsPerBin,numberOfPointsPerBinTMP)
      if(numberOfPointsPerBinTMP>0){
        binTMP[ibinTMP]<-sum(currentSpectrum[rowsTMP,colsTMP],na.rm=TRUE)/
                         numberOfPointsPerBinTMP
        #Set to NA: Make sure each point is counted only once for rectangular bins. For custom bins lists, double counting may be on purpose
        if(binMethod=="Rectangular bins") currentSpectrum[rowsTMP,colsTMP]<-NA
      }
   }
 } else {#1D spectra
   NMRspectrumNames<-as.numeric(names(currentSpectrum))
   if(TRUE){
     for(ibinTMP in 1:nrow(binRegions)){
        indexTMP<-NMRspectrumNames<=binRegions[ibinTMP,1]&
                                      NMRspectrumNames>binRegions[ibinTMP,2]
        numberOfPointsPerBinTMP<-sum(indexTMP)
        numberOfPointsPerBin<-c(numberOfPointsPerBin,numberOfPointsPerBinTMP)
        if(numberOfPointsPerBinTMP>0){
           binTMP[ibinTMP]<-sum(currentSpectrum[
                            indexTMP])/numberOfPointsPerBinTMP
           #Set to NA: Make sure each point is counted only once for rectangular bins. For custom bins lists, double counting may be on purpose
           #if(binMethod=="Rectangular bins") currentSpectrum[NMRspectrumNames<=binRegions[ibinTMP,1]&
           #                        NMRspectrumNames>binRegions[ibinTMP,2]]<-NA
        }
     }
   }
   if(FALSE){
     #avoid time-consuming loop - but this takes far too much memory in
     #parallel mode. Without parallel, it might save only little time
     indexMatrix<-matrix(rep(NMRspectrumNames,nrow(binRegions)),
                         nrow=nrow(binRegions),byrow=TRUE)
     valueMatrix<-matrix(rep(currentSpectrum,nrow(binRegions)),
                         nrow=nrow(binRegions),byrow=TRUE)
     indexMatrix2<-indexMatrix<=binRegions[,1]&
                   indexMatrix>binRegions[,2]
     numberOfPointsPerBin<-apply(indexMatrix2,1,sum)
     binTMP<-apply(valueMatrix*indexMatrix2,1,sum)/numberOfPointsPerBin
     names(binTMP)<-rownames(binRegions)
     binTMP[is.na(binTMP)]<-0
   }
  }
  if(is.null(numberOfPointsPerBin)){
    warning("No bin contains any data for a spectrum.")
  } #else {
    #if(!dataProvided){
    #  mrbin.env$mrbinTMP$meanNumberOfPointsPerBin_TMP<-numberOfPointsPerBin
    #  mrbin.env$mrbinTMP$binTMP<-binTMP
    #}
  #}
  invisible(list(binTMP=binTMP,pointsPerBin=numberOfPointsPerBin))
}

#' A function for calculating noise levels.
#'
#' This function calculates noise levels for the current spectrum. Only for use
#' within the mrbin function.
#' @param NMRdata Spectral data
#' @param pointsPerBin Mean number of data points per bin
#' @param dimension Dimension
#' @param noiseRange1d Noise range
#' @param noiseRange2d Noise range
#' @return An (invisible) object containing noise level
#' @keywords internal
#' @noRd
#' @examples
#' \donttest{ calculateNoise() }

calculateNoise<-function(NMRdata=NULL,pointsPerBin=NULL,dimension="1D",
   noiseRange1d=NULL,noiseRange2d=NULL){
  #dataProvided<-TRUE
  #if(is.null(NMRdata)){
  #  NMRdata<-mrbin.env$mrbinTMP$currentSpectrumOriginal
  #  dataProvided<-FALSE
  #}
  #if(is.null(pointsPerBin)){
  #  pointsPerBin<-mrbin.env$mrbinTMP$meanNumberOfPointsPerBin_TMP
  #}
  if(!is.null(NMRdata)){
    if(dimension=="1D"){
         noise_level<-stats::sd(NMRdata[
                 which(as.numeric(names(NMRdata))<=max(noiseRange1d[1:2])&
                      as.numeric(names(NMRdata))>=min(noiseRange1d[1:2]))])
    }
    if(dimension=="2D"){
         noise_level<-stats::sd(NMRdata[
               which(as.numeric(rownames(NMRdata))>=min(noiseRange2d[3:4])&
                 as.numeric(rownames(NMRdata))<=max(noiseRange2d[3:4])),
               which(as.numeric(colnames(NMRdata))<=max(noiseRange2d[1:2])&
                 as.numeric(colnames(NMRdata))>=min(noiseRange2d[1:2]))])
    }
    noise_level_TMP<-noise_level*(pointsPerBin^(-.5))
    #if(!dataProvided){
    #  mrbin.env$mrbinTMP$noise_level_Raw_TMP<-noise_level
    #  mrbin.env$mrbinTMP$noise_level_TMP<-noise_level_TMP
    #}
    invisible(list(noise_level=noise_level,noise_level_TMP=noise_level_TMP))
 }
}
