#mrbin - Collection of R functions for analyzing NMR metabolomics data.
#Written by Matthias Klein, The Ohio State University
#
#Package: mrbin
#Title: Binning, Scaling and Normalization of NMR Data
#Version: 1.1.0
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

#' @importFrom grDevices colorRamp heat.colors rainbow rgb devAskNewPage
#' @importFrom graphics axis contour hist legend lines par plot text boxplot points rect
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
#' @return NMRdata A matrix containing NMR data without negative values.
#' @export
#' @examples
#' ## Use own data
#' # atnv(NMRdataMatrix,noiseLevelVector)
#' ## Use current mrbin data. Requires data loaded using mrbin()
#' # atnv()

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
     return(NMRdata)
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
#' ## Don't run
#' # .onLoad()

.onLoad <- function(libname, pkgname){
    assign("mrbin.env",new.env(emptyenv()),parent.env(environment()))
    resetEnv()
}

#' A parameter resetting function
#'
#' This function resets the parameter variables.
#' @param
#' @return {None}
#' @export
#' @examples
#' ## Don't run
#' # resetEnv()

resetEnv<-function(){
    if(!exists("mrbin.env", mode="environment")) .onLoad()
    assign("bins",NULL,mrbin.env)
    assign("paramChangeFlag",FALSE,mrbin.env)
    assign("mrbinTMP",list(
               mrbinversion="1.1.0",
               binsRaw=NULL,
               medians=NULL,
               PCA=NULL,
               binTMP=NULL,
               binNames=NULL,
               nbins1=NULL,
               nbins2=NULL,
               nbins=NULL,
               binRegions=matrix(ncol=4,nrow=0,dimnames=list(NULL,c("left","right","top","bottom"))),
               currentFolder=NULL,
               currentSpectrum=NULL,
               currentSpectrumName=NULL
    ),mrbin.env)
    assign("mrbinparam", list(
               dimension="2D",
               binMethod="Rectangular bins",#"Special bin list"
               binwidth1D=.003,
               binwidth2D=.02,
               binheight=1,
               binRegion=c(9.5,0.5,10,156),
               specialBinList=matrix(c(
                               5.45,5.2,0,160,
                               2.9,2.74,0,160,
                               2.44,2.2,0,160,
                               2.44,2.3326,0,160,
                               2.3326,2.2,0,160,
                               2.14,1.93,0,160,
                               1.41,1.2,0,160,
                               0.94,0.8,0,160,
                               4.325,4.26,0,160
                               ),ncol=4,byrow=TRUE,dimnames=list(c(
                               "-CH=CH- Methene",
                               "=CH-CH2-CH= Diallylic",
                               "COO-CH2-CH2- Methylene_to_carboxyl",
                               "COO-CH2-CH2- Methylene_to_carboxyl_Free Fatty Acids",
                               "COO-CH2-CH2- Methylene_to_carboxyl_Triacylglycerides, Phospholipids",
                               "-CH2-CH=CH- Allylic",
                               "-CH2- Methylene",
                               "-CH3 Methyl",
                               "Glycerol"
                               ),
                                c("left","right","top","bottom"))),
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
               signal_to_noise1D=10,
               signal_to_noise2D=3,
               noiseRange2d=c(2.3,3.3,90,110),
               noiseRange1d=c(9.4,10),
               croptopRight=c(0,-1.50),#only 2D, defines edge points of the cropped area
               croptopLeft=c(0,3.5),
               cropbottomRight=c(160,6),
               cropbottomLeft=c(160,10),
               PQNsugarArea=c(5.4,3.35,72,100),#exclude most glucose to reduce impact
               PQNshowHist=FALSE,#show histograms of quotients
               PQNminimumFeatures=40,#Top number of features to proceed
               PCAtitlelength=8,
               createBins="Yes",
               readNMR="Yes",
               useAsNames="Folder names",#Use spectrum titles as sample names or folder names will be used
               NMRvendor="Bruker",#NMR vendor. Currently, only Bruker data is supported.
               mrbinversionTracker=mrbin.env$mrbinTMP$mrbinversion,
               saveFiles="No",
               readFromFile=TRUE,#FALSE only for testing example
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
               plotRegion=c(9.5,.5,-1,156),
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
#' ## Let the user set parameters interactively
#' # results <- mrbin()
#' ## Set all parameters in the console and run without user intervention
#' # results <- mrbin(silent=TRUE,
#' #              setDefault=FALSE,
#' #              parameters=list(
#' #              dimension="1D",
#' #              binMethod="Rectangular bins",
#' #              binwidth1D=.01,
#' #              referenceScaling="Yes",
#' #              removeSolvent="Yes",
#' #              removeAreas="No",
#' #              sumBins="No",
#' #              noiseRemoval="Yes",
#' #              PQNScaling="Yes",
#' #              fixNegatives="Yes",
#' #              logTrafo="Yes",
#' #              saveFiles="Yes",
#' #              outputFileName="mrbin_test_results",
#' #              NMRfolders=c("C:/NMR/Sample01/1/pdata/1",
#' #                           "C:/NMR/Sample19/1/pdata/1",
#' #                           "C:/NMR/Sample61/1/pdata/1")
#' #      ))
#' ## Test example, uses variables instead of files for speed
#' results <- mrbin(silent=TRUE,
#'                  parameters=list(dimension="1D",
#'                         binwidth1D=0.05,
#'                         signal_to_noise1D=35,
#'                         PCAtitlelength=20,
#'                         readFromFile = FALSE,#Only for testing, should be TRUE in all other cases
#'                         NMRfolders = c("spectrum1",
#'                                        "spectrum2",
#'                                        "spectrum3")#Usually this should be folder names
#'                         )
#'         )


mrbin<-function(silent=FALSE,setDefault=FALSE,parameters=NULL){
  if(!exists("mrbin.env", mode="environment")) .onLoad()
  message(paste("\nmrbin version ",mrbin.env$mrbinTMP$mrbinversion,"\n",
       "Binning, Scaling and Normalization of NMR Data\n",sep=""))
  if(setDefault) resetEnv()
  if(!is.null(parameters)){
      setParam(parameters)
  }
  stopTMP<-FALSE
  selectionRepeat<-""
  if(silent) startmrbin<-"Start binning now"
  #Create bin list?
  if(!silent){
   selection<-utils::select.list(c("Yes","No"),preselect="Yes",
              title="mrbin: Create bin list?",graphics=TRUE)
     if(!selection=="Yes")         stopTMP<-TRUE
     if(selection=="Yes"&!stopTMP){
       #Select new parameters?
       selectionRepeat<-utils::select.list(c("Set new","Use current","Reload from file"),preselect="Set new",
                                          title="Set parameters or use existing parameters?",graphics=TRUE)
       if(length(selectionRepeat)==0|selectionRepeat=="") stopTMP<-TRUE
       if(selectionRepeat=="Reload from file"&!stopTMP){
         recreatemrbin()
         selectionRepeat2<-utils::select.list(c("Adjust","Use as is"),preselect="Adjust",
                                          title="Adjust parameters or use as is?",graphics=TRUE)
         if(length(selectionRepeat2)==0|selectionRepeat=="") stopTMP<-TRUE
         if(selectionRepeat2=="Adjust"&!stopTMP) selectionRepeat<-"Set new"
       }
       if(selectionRepeat=="Set new"&!stopTMP){
         #1D or 2D data?
         dimension<-utils::select.list(c("1D","2D"),preselect=mrbin.env$mrbinparam$dimension,
                                 title="1D or 2D spectra?",graphics=TRUE)
         if(length(dimension)==0|dimension=="") stopTMP<-TRUE
          if(!stopTMP){
            if(!identical(mrbin.env$mrbinparam$dimension,dimension)) mrbin.env$paramChangeFlag<-TRUE
            mrbin.env$mrbinparam$dimension<-dimension
            if(dimension=="1D") dimlength<-2
            if(dimension=="2D") dimlength<-4
            #Use rectangluar bins or use special bin list, e.g. for lipids
            binMethod<-utils::select.list(c("Rectangular bins","Special bin list"),preselect=mrbin.env$mrbinparam$binMethod,
                       ,title ="Binning method: ",graphics=TRUE)
            if(length(binMethod)==0|binMethod=="") stopTMP<-TRUE
            if(!stopTMP){
              if(!identical(mrbin.env$mrbinparam$binMethod,binMethod)) mrbin.env$paramChangeFlag<-TRUE
              mrbin.env$mrbinparam$binMethod<-binMethod
              #Bin region
              adjRegion<-""
              if(mrbin.env$mrbinparam$binMethod=="Rectangular bins"){
                adjRegion<-utils::select.list(c(paste(mrbin.env$mrbinparam$binRegion[1:dimlength],
                           collapse=","),"Change..."),preselect=paste(mrbin.env$mrbinparam$binRegion[1:dimlength],
                           collapse=","),title ="Bin region [ppm]: ",graphics=TRUE)
                if(length(adjRegion)==0|adjRegion=="") stopTMP<-TRUE
              }
            }
          }
          if(!stopTMP&mrbin.env$mrbinparam$binMethod=="Rectangular bins"){
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
          }
          }
          #Bin width and height
          if(mrbin.env$mrbinparam$dimension=="1D"&!stopTMP&mrbin.env$mrbinparam$binMethod=="Rectangular bins"){
              widthAdjust<-utils::select.list(c(mrbin.env$mrbinparam$binwidth1D,"Change..."),preselect=
                           as.character(mrbin.env$mrbinparam$binwidth1D),title ="Bin width [ppm]: ",graphics=TRUE)
              if(length(widthAdjust)==0|widthAdjust=="") stopTMP<-TRUE
              if(widthAdjust=="Change..."){
                   widthTMP<-readline(prompt=paste("New 1D bin width, press enter to keep ",
                             mrbin.env$mrbinparam$binwidth1D,": ",sep=""))
                   if(!widthTMP=="") {
                       mrbin.env$paramChangeFlag<-TRUE
                       mrbin.env$mrbinparam$binwidth1D<-as.numeric(widthTMP)
                   }
              }
          }
          if(mrbin.env$mrbinparam$dimension=="2D"&!stopTMP&mrbin.env$mrbinparam$binMethod=="Rectangular bins"){
              widthAdjust<-utils::select.list(c(mrbin.env$mrbinparam$binwidth2D,"Change..."),preselect=
                           as.character(mrbin.env$mrbinparam$binwidth2D),title ="Bin width [ppm]: ",graphics=TRUE)
              if(length(widthAdjust)==0|widthAdjust=="") stopTMP<-TRUE
              if(widthAdjust=="Change..."){
                   widthTMP<-readline(prompt=paste("New 2D bin width, press enter to keep ",
                             mrbin.env$mrbinparam$binwidth2D,": ",sep=""))
                   if(!widthTMP=="") {
                       mrbin.env$paramChangeFlag<-TRUE
                       mrbin.env$mrbinparam$binwidth2D<-as.numeric(widthTMP)
                   }
              }
              if(!stopTMP){
                heightAdjust<-utils::select.list(c(mrbin.env$mrbinparam$binheight,"Change..."),preselect=
                              as.character(mrbin.env$mrbinparam$binheight),title ="Bin height [ppm]: ",graphics=TRUE)
                if(length(heightAdjust)==0|heightAdjust=="") stopTMP<-TRUE
                if(heightAdjust=="Change..."){
                     heightTMP<-readline(prompt=paste("New 2D bin height, press enter to keep ",
                                mrbin.env$mrbinparam$binheight,": ",sep=""))
                     if(!heightTMP=="") {
                       mrbin.env$paramChangeFlag<-TRUE
                       mrbin.env$mrbinparam$binheight<-as.numeric(heightTMP)
                     }
                }
              }
          }
          #Set special bin list
          if(!stopTMP&mrbin.env$mrbinparam$binMethod=="Special bin list"){
              for(ibinRegions in 1:nrow(mrbin.env$mrbinparam$specialBinList)){
               if(!stopTMP){
                adjbinRegion<-utils::select.list(c(paste(mrbin.env$mrbinparam$specialBinList[ibinRegions,1:dimlength],collapse=","),"Change..."),
                           preselect=paste(mrbin.env$mrbinparam$specialBinList[ibinRegions,1:dimlength],collapse=","),
                           title =paste("Change region?",rownames(mrbin.env$mrbinparam$specialBinList)[ibinRegions]),graphics=TRUE)
                if(length(adjbinRegion)==0|adjbinRegion=="") stopTMP<-TRUE
               }
              if(adjbinRegion=="Change..."&!stopTMP){
                nameTMP<-readline(prompt=paste("New name, press enter to keep ",
                          rownames(mrbin.env$mrbinparam$specialBinList)[ibinRegions],": ",sep=""))
                if(!nameTMP=="") {
                       mrbin.env$paramChangeFlag<-TRUE
                       rownames(mrbin.env$mrbinparam$specialBinList)[ibinRegions]<-nameTMP
                }
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
              addBinFlag<-TRUE
              if(!stopTMP){
               while(addBinFlag){
                  addbinRegion<-utils::select.list(c("Yes","No"),preselect="No",
                               title ="Add additional bin?",graphics=TRUE)
                  if(length(addbinRegion)==0|adjbinRegion=="") stopTMP<-TRUE
                  if(addbinRegion=="No"&!stopTMP) addBinFlag<-FALSE
                  if(addbinRegion=="Yes"&!stopTMP){
                    mrbin.env$paramChangeFlag<-TRUE
                    mrbin.env$mrbinparam$specialBinList<-rbind(mrbin.env$mrbinparam$specialBinList,c(0,0,0,160))
                    ibinRegions<-nrow(mrbin.env$mrbinparam$specialBinList)
                    nameTMP<-readline(prompt="New bin name: ")
                    #if(!nameTMP=="") {
                    rownames(mrbin.env$mrbinparam$specialBinList)[ibinRegions]<-nameTMP
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
              if(!stopTMP){
              mrbin.env$mrbinTMP$binRegions<-mrbin.env$mrbinparam$specialBinList
              }
          }
          #Scale to reference
          if(!stopTMP){
            referenceScaling<-utils::select.list(c("Yes","No"),preselect=mrbin.env$mrbinparam$referenceScaling,
                                          title = "Scale to reference signal?",graphics=TRUE)
            if(length(referenceScaling)==0|referenceScaling=="") stopTMP<-TRUE
          }
          if(!stopTMP){
            if(!identical(mrbin.env$mrbinparam$referenceScaling,referenceScaling)) mrbin.env$paramChangeFlag<-TRUE
            mrbin.env$mrbinparam$referenceScaling<-referenceScaling
          }
          if(mrbin.env$mrbinparam$referenceScaling=="Yes"&!stopTMP){
            if(mrbin.env$mrbinparam$dimension=="1D"&!stopTMP){
              adjRegion<-utils::select.list(c(paste(mrbin.env$mrbinparam$reference1D,collapse=","),"Change..."),
                         preselect=paste(mrbin.env$mrbinparam$reference1D,collapse=","),title ="Reference region [ppm]: ",graphics=TRUE)
              if(length(adjRegion)==0|adjRegion=="") stopTMP<-TRUE
              if(adjRegion=="Change..."&!stopTMP){
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
            if(mrbin.env$mrbinparam$dimension=="2D"&!stopTMP){
              adjRegion<-utils::select.list(c(paste(mrbin.env$mrbinparam$reference2D,collapse=","),"Change..."),
                         preselect=paste(mrbin.env$mrbinparam$reference2D,collapse=","),title ="Reference region [ppm]: ",graphics=TRUE)
              if(length(adjRegion)==0|adjRegion=="") stopTMP<-TRUE
              if(adjRegion=="Change..."&!stopTMP){
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
          #Remove solvent
          if(!stopTMP){
            removeSolvent<-utils::select.list(c("Yes","No"),preselect=mrbin.env$mrbinparam$removeSolvent,
                                     title = "Remove solvent area?",graphics=TRUE)
            if(length(removeSolvent)==0|removeSolvent=="") stopTMP<-TRUE
          }
          if(!stopTMP){
            mrbin.env$mrbinparam$removeSolvent<-removeSolvent
          }
          if(mrbin.env$mrbinparam$removeSolvent=="Yes"&!stopTMP){
              adjRegion<-utils::select.list(c(paste(mrbin.env$mrbinparam$solventRegion,collapse=","),"Change..."),
                         preselect=paste(mrbin.env$mrbinparam$solventRegion,collapse=","),title ="Solvent region [ppm]: ",graphics=TRUE)
              if(length(adjRegion)==0|adjRegion=="") stopTMP<-TRUE
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
          #Remove additional areas
          if(!stopTMP){
            removeAreas<-utils::select.list(c("Yes","No"),preselect=mrbin.env$mrbinparam$removeAreas,
                                     title = "Remove additional areas?",graphics=TRUE)
            if(length(removeAreas)==0|removeAreas=="") stopTMP<-TRUE
          }
          if(!stopTMP){
            mrbin.env$mrbinparam$removeAreas<-removeAreas
          }
          if(mrbin.env$mrbinparam$removeAreas=="Yes"&!stopTMP){
            addAreasFlag<-TRUE
            if(nrow(mrbin.env$mrbinparam$removeAreaList)>0&!stopTMP){
                addAreasFlag<-FALSE
                removeAreaListTMP<-utils::select.list(c("Keep","New","Add to existing list"),preselect="Keep",
                                   title = "Use previous area list or define new?",graphics=TRUE)
                if(length(removeAreaListTMP)==0|removeAreaListTMP=="") stopTMP<-TRUE
                if(!removeAreaListTMP=="Keep"&!stopTMP){
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
              keepAdding<-utils::select.list(c("No","Yes"),preselect="No",multiple=FALSE,title = "Add more areas?",graphics=TRUE)
              if(length(keepAdding)==0|keepAdding=="") stopTMP<-TRUE
              if(keepAdding=="No"&!stopTMP)  addAreasFlag<-FALSE
            }
          }
          #Sum up bins of unstable peaks
          if(!stopTMP){
            sumBins<-utils::select.list(c("Yes","No"),preselect=mrbin.env$mrbinparam$sumBins,
                                 title = "Sum bins of unstable peaks?",graphics=TRUE)
            if(length(sumBins)==0|sumBins=="") stopTMP<-TRUE
          }
          if(!stopTMP){
            mrbin.env$mrbinparam$sumBins<-sumBins
          }
          if(mrbin.env$mrbinparam$sumBins=="Yes"&!stopTMP){
              addAreasFlag<-TRUE
              if(nrow(mrbin.env$mrbinparam$sumBinList)>0&!stopTMP){
                  addAreasFlag<-FALSE
                  sumBinListTMP<-utils::select.list(c("Keep","New","Add to existing list"),preselect="Keep",
                                 title = "Use previous area list or define new?",graphics=TRUE)
                  if(length(sumBinListTMP)==0|sumBinListTMP=="") stopTMP<-TRUE
                  if(!sumBinListTMP=="Keep"&!stopTMP){
                    addAreasFlag<-TRUE
                    if(!sumBinListTMP=="New"&!stopTMP){
                        mrbin.env$mrbinparam$sumBinList<-matrix(ncol=4,nrow=0,dimnames=list(NULL,c("left","right","top","bottom")))
                    }
                  }
              }
              if(!stopTMP){
                iaddAreas<-nrow(mrbin.env$mrbinparam$sumBinList)+1
              }
              while(addAreasFlag&!stopTMP){
                mrbin.env$mrbinparam$sumBinList<-rbind(mrbin.env$mrbinparam$sumBinList,c(0,0,0,0))
                mrbin.env$mrbinparam$sumBinList[iaddAreas,1]<-as.numeric(readline(prompt="Left border: "))
                mrbin.env$mrbinparam$sumBinList[iaddAreas,2]<-as.numeric(readline(prompt="Right border: "))
                if(mrbin.env$mrbinparam$dimension=="2D"&!stopTMP){
                    mrbin.env$mrbinparam$sumBinList[iaddAreas,3]<-as.numeric(readline(prompt="Top border: "))
                    mrbin.env$mrbinparam$sumBinList[iaddAreas,4]<-as.numeric(readline(prompt="Bottom border: "))
                }
                iaddAreas<-iaddAreas+1
                keepAdding<-utils::select.list(c("No","Yes"),preselect="No",title = "Add more areas?",graphics=TRUE)
                if(length(keepAdding)==0|keepAdding=="") stopTMP<-TRUE
                if(keepAdding=="No"&!stopTMP)  addAreasFlag<-FALSE
              }
          }
          #Remove noise
          if(!stopTMP){
            noiseRemoval<-utils::select.list(c("Yes","No"),
                                     preselect=mrbin.env$mrbinparam$noiseRemoval,title="Remove noise?",graphics=TRUE)
            if(length(noiseRemoval)==0|noiseRemoval=="") stopTMP<-TRUE
          }
          if(!stopTMP){
            mrbin.env$mrbinparam$noiseRemoval<-noiseRemoval
          }
          if(mrbin.env$mrbinparam$noiseRemoval=="Yes"&!stopTMP){
            if(mrbin.env$mrbinparam$dimension=="1D"&!stopTMP){
              SNRTMP<-utils::select.list(c(mrbin.env$mrbinparam$signal_to_noise1D,"Change..."),preselect=
                      as.character(mrbin.env$mrbinparam$signal_to_noise1D),title="Signal-to-noise ratio (SNR):",graphics=TRUE)
              if(length(SNRTMP)==0|SNRTMP=="") stopTMP<-TRUE
              if(SNRTMP=="Change..."&!stopTMP) SNRTMP<-readline(prompt=paste(
                                   "New 1D signal to noise ratio, press enter to keep ",
                                   mrbin.env$mrbinparam$signal_to_noise1D,": ",sep=""))
              if(!SNRTMP==""&!stopTMP) {
                  mrbin.env$mrbinparam$signal_to_noise1D<-as.numeric(SNRTMP)
              }
            }
            if(mrbin.env$mrbinparam$dimension=="2D"&!stopTMP){
              SNRTMP<-utils::select.list(c(mrbin.env$mrbinparam$signal_to_noise2D,"Change..."),preselect=
                      as.character(mrbin.env$mrbinparam$signal_to_noise2D),title="Signal-to-noise ratio (SNR):",graphics=TRUE)
              if(length(SNRTMP)==0|SNRTMP=="") stopTMP<-TRUE
              if(SNRTMP=="Change..."&!stopTMP) SNRTMP<-readline(prompt=paste(
                                   "New 2D signal to noise ratio, press enter to keep ",
                                   mrbin.env$mrbinparam$signal_to_noise2D,": ",sep=""))
              if(!SNRTMP==""&!stopTMP) {
                   mrbin.env$mrbinparam$signal_to_noise2D<-as.numeric(SNRTMP)
              }
            }
            if(!stopTMP){
              noiseTMP<-utils::select.list(unique(c(as.character(mrbin.env$mrbinparam$noiseThreshold),"0.2","0.75","0.05","Change...")),
                                     preselect=as.character(mrbin.env$mrbinparam$noiseThreshold),title="Minimum ratio > SNR",graphics=TRUE)
              if(length(noiseTMP)==0|noiseTMP=="") stopTMP<-TRUE
            }
            if(!stopTMP){
              if(noiseTMP=="Change..."&!stopTMP){
                    noiseTMP<-readline(prompt=paste("New noise threshold, press enter to keep ",
                              mrbin.env$mrbinparam$noiseThreshold,": ",sep=""))
                    if(!noiseTMP=="")  {
                       mrbin.env$mrbinparam$noiseThreshold<-as.numeric(noiseTMP)
                    }
              } else {
                mrbin.env$mrbinparam$noiseThreshold<-as.numeric(noiseTMP)
              }
            }
          }
          #Crop HSQCs
          if(mrbin.env$mrbinparam$dimension=="2D"&!stopTMP){
               cropHSQC<-utils::select.list(c("Yes","No"),preselect=mrbin.env$mrbinparam$cropHSQC,
                                     title="Crop H-C HSQCs?",graphics=TRUE)
               if(length(cropHSQC)==0|cropHSQC=="") stopTMP<-TRUE
               if(!stopTMP){
                 mrbin.env$mrbinparam$cropHSQC<-cropHSQC
               }
          }
          #PQN scaling?
          if(!stopTMP){
            PQNScaling<-utils::select.list(c("Yes","No"),preselect=mrbin.env$mrbinparam$PQNScaling,
                                    title = "PQN normalization?",graphics=TRUE)
            if(length(PQNScaling)==0|PQNScaling=="") stopTMP<-TRUE
          }
          if(!stopTMP){
            mrbin.env$mrbinparam$PQNScaling<-PQNScaling
          }
          #Replace negative values
          if(!stopTMP){
            fixNegatives<-utils::select.list(c("Yes","No"),preselect=mrbin.env$mrbinparam$fixNegatives,
                                      title="Replace negative values?",graphics=TRUE)
            if(length(fixNegatives)==0|fixNegatives=="") stopTMP<-TRUE
          }
          if(!stopTMP){
            mrbin.env$mrbinparam$fixNegatives<-fixNegatives
          }
          #Log scaling?
          if(!stopTMP){
            logTrafo<-utils::select.list(c("Yes","No"),preselect=mrbin.env$mrbinparam$logTrafo,
                                  title="Log transformation?",graphics=TRUE)
            if(length(logTrafo)==0|logTrafo=="") stopTMP<-TRUE
          }
          if(!stopTMP){
            mrbin.env$mrbinparam$logTrafo<- logTrafo
          }
          #PCA
          if(!stopTMP){
            PCA<-utils::select.list(c("Yes","No"),preselect = mrbin.env$mrbinparam$PCA,
                              title = "Create result plot?",graphics=TRUE)
            if(length(PCA)==0|PCA=="") stopTMP<-TRUE
          }
          if(!stopTMP){
            mrbin.env$mrbinparam$PCA<-PCA
          }
          if(!stopTMP&mrbin.env$mrbinparam$PCA=="Yes"){
            PCAtitlelength<-utils::select.list(c(sort(unique(c(2,4,6,8,
                              mrbin.env$mrbinparam$PCAtitlelength))),
                              "Use whole title","Custom..."),
                              preselect = as.character(mrbin.env$mrbinparam$PCAtitlelength),
                              title = "Crop titles for plot?",graphics=TRUE)
            if(length(PCAtitlelength)==0|PCAtitlelength=="") stopTMP<-TRUE
          }
          if(!stopTMP&mrbin.env$mrbinparam$PCA=="Yes"){
                 if(PCAtitlelength=="Use whole title") PCAtitlelength<-"500"
          }
          if(!stopTMP&mrbin.env$mrbinparam$PCA=="Yes"){
                  if(PCAtitlelength=="Custom..."){
                      PCAtitlelengthTMP<-readline(prompt=paste("New title length, press enter to keep ",
                          mrbin.env$mrbinparam$PCAtitlelength,": ",sep=""))
                      if(!PCAtitlelengthTMP=="") {
                          PCAtitlelength<-PCAtitlelengthTMP
                      } else {
                          PCAtitlelength<-as.character(mrbin.env$mrbinparam$PCAtitlelength)
                      }
                  }
          }
          if(!stopTMP&mrbin.env$mrbinparam$PCA=="Yes"){
                 mrbin.env$mrbinparam$PCAtitlelength<-as.numeric(PCAtitlelength)
          }
          #Select folders
          if(!stopTMP){
            if(!is.null(mrbin.env$mrbinparam$NMRfolders)){
                selectionFolders<-utils::select.list(c("Yes","No","Remove spectra from previous list"),preselect="Yes",
                                  title="Use previous spectra list?",graphics=TRUE)
                if(length(selectionFolders)==0|selectionFolders=="") stopTMP<-TRUE
                if(!stopTMP){
                  if(selectionFolders=="No")  selectFolders()
                }
                if(!stopTMP){
                  if(selectionFolders=="Remove spectra from previous list")  removeSpectrum()
                }
            } else {
                selectFolders()
            }
          }
          #Define groups
          if(!stopTMP){
            defineGroups<-utils::select.list(c("Yes","No"),preselect=mrbin.env$mrbinparam$defineGroups,
                                      title = "Define group members?",graphics=TRUE)
            if(length(defineGroups)==0|defineGroups=="") stopTMP<-TRUE
          }
          if(!stopTMP){
            mrbin.env$mrbinparam$defineGroups<-defineGroups
          }
          if(!stopTMP){
            if(mrbin.env$mrbinparam$defineGroups=="Yes"&!stopTMP){
                if(!is.null(mrbin.env$mrbinparam$Factors)&!stopTMP){
                   if(length(mrbin.env$mrbinparam$Factors)==length(mrbin.env$mrbinparam$NMRfolders)&!stopTMP){
                     selectionFactors<-utils::select.list(c("Yes","No"),preselect="Yes",
                                       title = "Use previous factor list?",graphics=TRUE)
                     if(length(selectionFactors)==0|selectionFactors=="") stopTMP<-TRUE
                     if(!stopTMP){
                       if(selectionFactors=="No")  setFactors()
                     }
                   }
                } else {
                     setFactors()
                }
            }
          }
          #Define sample names
          if(!stopTMP){
            useAsNames<-utils::select.list(c("Folder names","Spectrum titles"),
                                      preselect=mrbin.env$mrbinparam$useAsNames,
                                      title = "Create sample names from",graphics=TRUE)
            if(length(useAsNames)==0|useAsNames=="") stopTMP<-TRUE
          }
          if(!stopTMP){
            mrbin.env$mrbinparam$useAsNames<-useAsNames
          }
   }
   #Save output files to hard drive?
   if(!stopTMP){
     saveFilesTMP<-utils::select.list(c("Yes","No"),preselect="Yes",
                 title ="Save output files to hard drive?",graphics=TRUE)
     if(length(saveFilesTMP)==0|saveFilesTMP=="") stopTMP<-TRUE
   }
   if(!stopTMP){
    mrbin.env$mrbinparam$saveFiles<-saveFilesTMP
    if(mrbin.env$mrbinparam$saveFiles=="Yes"&!stopTMP){
      enterFoldersTMP<-readline(prompt="Enter output folder path: ")
      if(enterFoldersTMP=="") stopTMP<-TRUE
      if(!stopTMP){
       parentFolder<-gsub('\\\\',"/",enterFoldersTMP)
       filenameTMP<-utils::select.list(c(paste("mrbin_",gsub(":","-",gsub(" ","_",Sys.time())),
                   sep=""),"Change..."),
                   title ="Output file name: ",graphics=TRUE)
       if(length(filenameTMP)==0|filenameTMP=="") stopTMP<-TRUE
       if(!stopTMP){
        if(filenameTMP=="Change..."&!stopTMP){
          filenameTMP<-readline(prompt=paste("New file name, press enter to use ",
                    paste("mrbin_",gsub(":","-",gsub(" ","_",Sys.time())),sep=""),": \n",sep=""))
          if(!filenameTMP=="") filenameTMP
          if(filenameTMP=="") filenameTMP<-paste("mrbin_",gsub(":","-",gsub(" ","_",Sys.time())),
                     sep="")
         }
         mrbin.env$mrbinparam$outputFileName<-gsub("//","/",paste(
                          parentFolder,"/",filenameTMP,sep=""))
       }
      }
     }
   }
   if(!stopTMP){
     startmrbin<-utils::select.list(c("Start binning now","I'll do it later"),
                preselect="Start binning now",title = "You're all set. Start binning?",graphics=TRUE)
     if(length(startmrbin)==0|startmrbin=="") stopTMP<-TRUE
   }
   if(!stopTMP){
     if(startmrbin=="Start binning now"&!stopTMP){
         if(!is.null(mrbin.env$mrbinTMP$binsRaw)){
             if(!mrbin.env$paramChangeFlag){
                 if(nrow(mrbin.env$mrbinTMP$binsRaw)==length(mrbin.env$mrbinparam$NMRfolders)) {
                      createBinNumbers()
                      if(mrbin.env$mrbinTMP$nbins==ncol(mrbin.env$mrbinTMP$binsRaw)){
                           keepDataTMP<-utils::select.list(c("Create new bin data (recommended)",
                                        "Use existing bin data (may cause inconsistent data)"),
                                        preselect="Create new bin data (recommended)",title =
                                        "Raw bin data exists in memory.",graphics=TRUE)
                           if(length(keepDataTMP)==0|keepDataTMP=="") stopTMP<-TRUE
                           if(keepDataTMP=="Use existing bin data (may cause inconsistent data)"&!stopTMP){
                                mrbin.env$bins<-mrbin.env$mrbinTMP$binsRaw
                                mrbin.env$mrbinparam$createBins<-"No"
                           }

                      }
                 }
             }
         }
     }
   }
  }
 }
 if(!stopTMP){
   if(startmrbin=="Start binning now"){
     mrbinrun()
     mrbin.env$paramChangeFlag<-FALSE
     mrbin.env$mrbinparam$createBins<-"Yes"
     resultOutputTMP<-c("\nNumber of spectra: ",nrow(mrbin.env$bins),"\n",
         "Number of bins at start: ",mrbin.env$mrbinparam$numberOfFeaturesRaw,"\n")
     if(!is.null(mrbin.env$mrbinparam$numberOfFeaturesAfterRemovingSolvent)){
          resultOutputTMP<-c(resultOutputTMP,
             "Number of bins after removing solvent: ",mrbin.env$mrbinparam$numberOfFeaturesAfterRemovingSolvent,"\n")
     }
     if(!is.null(mrbin.env$mrbinparam$numberOfFeaturesAfterRemovingAreas)){
          resultOutputTMP<-c(resultOutputTMP,
             "Number of bins after removing areas: ",mrbin.env$mrbinparam$numberOfFeaturesAfterRemovingAreas,"\n")
     }
     if(!is.null(mrbin.env$mrbinparam$numberOfFeaturesAfterSummingBins)){
          resultOutputTMP<-c(resultOutputTMP,
             "Number of bins after summing bins: ",mrbin.env$mrbinparam$numberOfFeaturesAfterSummingBins,"\n")
     }
     if(!is.null(mrbin.env$mrbinparam$numberOfFeaturesAfterNoiseRemoval)){
          resultOutputTMP<-c(resultOutputTMP,
             "Number of bins after noise removal: ",mrbin.env$mrbinparam$numberOfFeaturesAfterNoiseRemoval,"\n")
     }
     if(!is.null(mrbin.env$mrbinparam$numberOfFeaturesAfterCropping)){
          resultOutputTMP<-c(resultOutputTMP,
             "Number of bins after cropping: ",mrbin.env$mrbinparam$numberOfFeaturesAfterCropping,"\n")
     }
     resultOutputTMP<-paste( resultOutputTMP,sep="")
     message(resultOutputTMP)
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
#' @param
#' @return {None}
#' @export
#' @examples
#' ## Only works after setting all parameters
#' # mrbinrun()

mrbinrun<-function(){
  if(!exists("mrbin.env", mode="environment")) .onLoad()
  if(!is.null(mrbin.env$mrbinparam$NMRfolders)){
    #if(!"mrbinparam"%in%ls(envir = mrbin.env))    resetEnv()
    if(mrbin.env$mrbinparam$createBins=="Yes") binMultiNMR()
    if(mrbin.env$mrbinparam$removeSolvent=="Yes") removeSolvent()
    if(mrbin.env$mrbinparam$removeAreas=="Yes") removeAreas()
    if(mrbin.env$mrbinparam$sumBins=="Yes") sumBins()
    if(mrbin.env$mrbinparam$noiseRemoval=="Yes") removeNoise()
    if(mrbin.env$mrbinparam$cropHSQC=="Yes"&mrbin.env$mrbinparam$dimension=="2D") cropNMR()
    if(mrbin.env$mrbinparam$PQNScaling=="Yes") PQNScaling()
    if(mrbin.env$mrbinparam$fixNegatives=="Yes") atnv()
    if(mrbin.env$mrbinparam$logTrafo=="Yes") logTrafo()
    if(mrbin.env$mrbinparam$saveFiles=="Yes"){
      dput(mrbin.env$mrbinparam, file = paste(mrbin.env$mrbinparam$outputFileName,".txt",sep=""))
      utils::write.csv(mrbin.env$bins, file = paste(mrbin.env$mrbinparam$outputFileName,"bins.csv",sep=""))
    }
    if(mrbin.env$mrbinparam$PCA=="Yes") plotResults()
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
#' # recreatemrbin()
#' ## Insert full folder path and file name
#' # recreatemrbin("mrbin_01-03-20_12-00-10.txt")


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
     #if(!singleFolderFlag){
         selectFolders<-utils::select.list(folderListTMP,preselect=NULL,multiple=FALSE,
                          title = paste("Browse to parameter file:",parentFolder),graphics=TRUE)

     #}
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
    message("No filename was specified. No file loaded.")
 }
}

#' A function setting parameters and checking for consistency.
#'
#' This function set parameters and checks parameters for consistency.
#' @param parameters List of parameters to be set
#' @return {None}
#' @export
#' @examples
#' ## Provide the parameters you like to set
#' # setParam(list(dimension="1D"))

setParam<-function(parameters=NULL){
  if(!is.null(parameters)){
    mrbin.env$mrbinparam_copy<-mrbin.env$mrbinparam
    mrbin.env$mrbinparam<-parameters
    diffSet<-setdiff(names(mrbin.env$mrbinparam_copy),names(mrbin.env$mrbinparam))
    diffSet3<-setdiff(diffSet,c("mrbinversionTracker", "medians", "noise_level",
               "numberOfFeaturesRaw", "numberOfFeaturesAfterRemovingSolvent",
               "numberOfFeaturesAfterRemovingAreas","numberOfFeaturesAfterSummingBins",
               "numberOfFeaturesAfterNoiseRemoval", "numberOfFeaturesAfterCropping"))
    diffSet2<-setdiff(names(mrbin.env$mrbinparam),names(mrbin.env$mrbinparam_copy))
    if(length(diffSet2)>0){
       message(paste("Unexpected parameters: ",
           paste(diffSet2,sep=", ", collapse=", "),"\n",
           "These parameters will not be used. Potentially they were created in a different mrbin version.\n",
           "Check for new versions at: kleinomicslab.com", sep=""))
    }
    if(length(diffSet3)>0){
       message(paste("Previous or default values are used for missing parameters: ",
           paste(diffSet,sep=", ", collapse=", "),sep=""))
    }
    if(length(diffSet)>0){
       for(idiffSet in diffSet){
            mrbin.env$mrbinparam[[idiffSet]]<-mrbin.env$mrbinparam_copy[[idiffSet]]
       }
    }
    if(!mrbin.env$mrbinTMP$mrbinversion==mrbin.env$mrbinparam$mrbinversionTracker){
       message(paste("Imported file was created using the older mrbin version ",mrbin.env$mrbinparam$mrbinversionTracker,
           ".\n For exact reproduction of results, please get the old version at: kleinomicslab.com",
           sep=""))
    } else {
        if("NMRfolders"%in%names(parameters)){
          message(paste("mrbin.env$mrbinparam$NMRfolders may need adjustment if folders were renamed.",
           sep=""))
        }
    }
  } else {
      message("No parameters were provided.")
  }
}


#' A function for log transforming data.
#'
#' This function simply log transforms. Will not work with negative data.
#' @param
#' @return {None}
#' @export
#' @examples
#' ## Will work only after creating bins
#' # logTrafo()


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

#' A function for setting the current spectrum.
#'
#' This function lets the user pick a spectrum from the list of spectra
#' analysis.
#' @param
#' @return {None}
#' @export
#' @examples
#' ## Will work only after setting a list of spectra
#' # setCurrentSpectrum()

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
#' analysis.
#' @param
#' @return {None}
#' @export
#' @examples
#' ## Will work only after creating a spectrum list
#' # removeSpectrum()

removeSpectrum<-function(){
 if(!is.null(mrbin.env$mrbinparam$NMRfolders)){
    listTMP<-utils::select.list(mrbin.env$mrbinparam$NMRfolders,preselect = NULL, multiple = TRUE,title ="Select spectra to be removed",graphics=TRUE)
    if(length(listTMP)>0){
      if(!listTMP==""){
         if(length(mrbin.env$mrbinparam$Factors)==length(mrbin.env$mrbinparam$NMRfolders)){
           mrbin.env$mrbinparam$Factors<-mrbin.env$mrbinparam$Factors[-which(mrbin.env$mrbinparam$NMRfolders%in%listTMP)]
         }
         if(nrow(mrbin.env$bins)==length(mrbin.env$mrbinparam$NMRfolders)){
          if((nrow(mrbin.env$bins)-listTMP)==1){
            rownamesTMP<-rownames(mrbin.env$bins)[-which(rownames(mrbin.env$bins)%in%listTMP)]
            colnamesTMP<-colnames(mrbin.env$bins)
            mrbin.env$bins<-matrix(mrbin.env$bins[-which(rownames(mrbin.env$bins)%in%listTMP),],nrow=1)
            rownames(mrbin.env$bins)<-rownamesTMP
            colnames(mrbin.env$bins)<-colnamesTMP
          } else {
             mrbin.env$bins<-mrbin.env$bins[-which(rownames(mrbin.env$bins)%in%listTMP),]
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
#' groups.
#' @param
#' @return {None}
#' @export
#' @examples
#' ## Will work only after setting a sectrum list
#' # setFactors()

setFactors<-function(){
   if(length(mrbin.env$mrbinparam$NMRfolders)>0){
      Factors<-rep("Group 0",length(mrbin.env$mrbinparam$NMRfolders))
      names(Factors)<-mrbin.env$mrbinparam$NMRfolders
      flag<-TRUE
      i<-0
      while(flag){
        i<-i+1
        listTMP<-utils::select.list(names(Factors),preselect = NULL, multiple = TRUE,title ="Please select group members")
        groupNameTMP<-utils::select.list(c(paste("Group",i),"Enter new name"),preselect=paste("Group",i),
                  multiple=FALSE,title ="Group name?",graphics=TRUE)
        if(groupNameTMP=="Enter new name"){
          groupNameTMP<-readline(prompt=paste("New group name (enter to use \"Group ",i,"\"): ",sep=""))
          if(groupNameTMP=="") groupNameTMP<-paste("Group ",i,sep="")
        }
        Factors[listTMP]<-groupNameTMP
        select<-utils::select.list(c("Yes","No"),preselect = "No",multiple = FALSE,
            title = paste("Define more groups?",sep=""),graphics=TRUE)
        if(select=="No") flag<-FALSE
      }
      Factors<-as.factor(Factors)
      mrbin.env$mrbinparam$Factors <- Factors
   }
}

#' A function for selecting NMR data folders.
#'
#' This function selects the correct folder selection function for the vendor
#' (currently only Bruker).
#' @param
#' @return {None}
#' @export
#' @examples
#' ## Requires user input
#' # selectFolders()

selectFolders<-function(){#Select NMR spectral folders
      if(mrbin.env$mrbinparam$NMRvendor=="Bruker"){
          selectBrukerFolders()
      }  else {
          stop(paste("No folder selection function defined for vendor ",mrbin.env$mrbinparam$NMRvendor,".\n",sep=""))
      }
}

#' A function for selecting Bruker NMR data folders.
#'
#' This function lets the user set NMR data folders (for Bruker data).
#' @param
#' @return {None}
#' @export
#' @examples
#' ## Requires user input
#' # selectBrukerFolders()

selectBrukerFolders<-function(){#Select Bruker NMR spectral folders
  NMRfolders<-NULL
  datanameDict<-c("1r","2rr")
  names(datanameDict)<-c("1D","2D")
  datanameTmp<-datanameDict[mrbin.env$mrbinparam$dimension]
  singleFolderFlag<-FALSE
  enterFolders<-utils::select.list(c("Browse...","Enter parent folder path manually"),
                 preselect="Browse...",title="Set NMR parent folder:",graphics=TRUE)
  enterFoldersTMP<-readline(prompt="Enter folder path: ")
  if(!enterFoldersTMP==""){
   parentFolder<-gsub('\\\\',"/",enterFoldersTMP)
   if(enterFolders=="Enter parent folder path manually"){
        selectFolders<-parentFolder
        selectList<-parentFolder
        singleFolderFlag<-TRUE
   }
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
                          title = paste("Browse to NMR parent folder:",parentFolder),graphics=TRUE)
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
      if(singleFolderFlag){
         message("Looking for NMR data in folder...")
         utils::flush.console()
         singleFolderFlag<-FALSE
         spectrum_path_list<-NULL
         subdirsTmp0<-list.dirs(path = selectList,recursive = FALSE,full.names=FALSE)
         if(length(subdirsTmp0)>0){
           for(isubDirs in 1:length(subdirsTmp0)){#Look for Bruker NMR data in folder.
             subdirsTmp<-list.dirs(path = paste(selectList,"/",subdirsTmp0[isubDirs],sep=""),recursive = FALSE,full.names=FALSE)
             if(length(subdirsTmp)>0){
               if(sum(sapply(suppressWarnings(as.numeric(subdirsTmp)),is.numeric))>0){ #Look for folders "1", "2" etc (EXPNO)
                   spectrum_path2<-paste(selectList,"/",subdirsTmp0[isubDirs],"/",subdirsTmp[which(sapply(suppressWarnings(as.numeric(subdirsTmp)),is.numeric))],sep="")
                   for(i in spectrum_path2){
                      if(dir.exists(paste(i,"/pdata",sep=""))){
                          subdirsTmp2<-list.dirs(path = paste(i,"/pdata",sep=""),recursive = FALSE,full.names=FALSE)
                          if(sum(sapply(suppressWarnings(as.numeric(subdirsTmp2)),is.numeric))>0){ #Look for folders "1", "2" etc (PROCNO)
                                spectrum_path3<-paste(i,"/pdata/",subdirsTmp2[which(sapply(suppressWarnings(as.numeric(subdirsTmp2)),is.numeric))],sep="")
                                spectrum_path_list<-c(spectrum_path_list,spectrum_path3)
                          }
                      }
                   }
               }
             }
             spectrum_proc_path<-NULL
             if(length(spectrum_path_list)>0){
               for(i in 1:length(spectrum_path_list)){
                   if(datanameTmp%in%list.files(spectrum_path_list[i])){
                       spectrum_proc_path<-c(spectrum_proc_path,spectrum_path_list[i])
                   }
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
           NMRfolders<-c(NMRfolders,
                     utils::select.list(spectrum_proc_path,preselect=NULL,multiple=TRUE,
                     title = "Select data sets",graphics=TRUE))
           yesorno<-utils::select.list(c("No","Yes"),preselect="No",multiple=FALSE,title="Add more spectra?",graphics=TRUE)
           if(yesorno=="Yes"){
               selectFlag<-0
           } else {
               selectFlag<-1
           }
         }
     } else {
        parentFolder<-selectList
     }
   }
   }
  }
  if(!is.null(NMRfolders)) mrbin.env$mrbinparam$NMRfolders<-NMRfolders
}

#' A function for binning multiple NMR spectra.
#'
#' This function creates bins for each spectrum in mrbin.env$mrbinparam$NMRfolders and
#' saves the bins to mrbin.env$bins.
#' @param
#' @return {None}
#' @export
#' @examples
#' ## Requires spectrum list
#' # binMultiNMR()

binMultiNMR<-function(){
 if(!is.null(mrbin.env$mrbinparam$NMRfolders)){
    mrbin.env$bins<-NULL
    mrbin.env$mrbinTMP$binNames<-NULL
    mrbin.env$mrbinTMP$binTMP<-NULL
    #Open and bin all spectra
    if(mrbin.env$mrbinparam$readNMR=="Yes"){
      mrbin.env$mrbinTMP$binsRaw<-NULL
      for(i in 1:length(mrbin.env$mrbinparam$NMRfolders)){
          message(paste("Binning spectrum ",i,sep=""))
          utils::flush.console()
          mrbin.env$mrbinTMP$currentFolder<-mrbin.env$mrbinparam$NMRfolders[i]
          readNMR()
          mrbin.env$mrbinTMP$binTMP<-NULL
          binSingleNMR()
          if(is.null(mrbin.env$mrbinTMP$binsRaw)){
              mrbin.env$mrbinTMP$binsRaw<-matrix(rep(0,length(mrbin.env$mrbinTMP$binTMP)*
                                            length(mrbin.env$mrbinparam$NMRfolders)),
                                            nrow=length(mrbin.env$mrbinparam$NMRfolders))
              colnames(mrbin.env$mrbinTMP$binsRaw)<-names(mrbin.env$mrbinTMP$binTMP)
              rownames(mrbin.env$mrbinTMP$binsRaw)<-1:length(mrbin.env$mrbinparam$NMRfolders)
          }
          if(is.null(mrbin.env$mrbinTMP$binTMP)){
              stop("Error while binning spectrum.")
          } else {
              mrbin.env$mrbinTMP$binsRaw[i,]<-mrbin.env$mrbinTMP$binTMP
              rownames(mrbin.env$mrbinTMP$binsRaw)[i]<-mrbin.env$mrbinTMP$currentSpectrumName#rev(strsplit(mrbin.env$mrbinTMP$currentFolder, "/")[[1]])[1]
          }
      }
    }
    mrbin.env$bins<-mrbin.env$mrbinTMP$binsRaw
    mrbin.env$mrbinparam$numberOfFeaturesRaw<-ncol(mrbin.env$bins)
 }
}

#' A function for creating bin titles.
#'
#' This function creates titles for the bins to represent their average ppm
#' value.
#' @param
#' @return {None}
#' @export
#' @examples
#' ## Requires data
#' # createBinNames()

createBinNames<-function(){
   createBinNumbers()

   if(mrbin.env$mrbinparam$binMethod=="Special bin list"){

   }
   if(mrbin.env$mrbinparam$binMethod=="Rectangular bins"){

       mrbin.env$mrbinTMP$binRegions<-matrix(ncol=4,
                                        nrow=mrbin.env$mrbinTMP$nbins,
                                        dimnames=list(NULL,c("left","right","top","bottom")))
       if(mrbin.env$mrbinparam$dimension=="1D"){
              mrbin.env$mrbinTMP$binRegions[,1]<-mrbin.env$mrbinparam$binRegion[1]-((1:mrbin.env$mrbinTMP$nbins)-1)*mrbin.env$mrbinparam$binwidth1D
              mrbin.env$mrbinTMP$binRegions[,2]<-mrbin.env$mrbinparam$binRegion[1]-((1:mrbin.env$mrbinTMP$nbins))*mrbin.env$mrbinparam$binwidth1D
              rownames(mrbin.env$mrbinTMP$binRegions)<-apply(mrbin.env$mrbinTMP$binRegions[,1:2],1,mean)
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
              rownames(mrbin.env$mrbinTMP$binRegions)<-paste(apply(mrbin.env$mrbinTMP$binRegions[,3:4],1,mean),
                                                        apply(mrbin.env$mrbinTMP$binRegions[,1:2],1,mean),sep=",")

       }
       mrbin.env$mrbinTMP$binNames<-rownames(mrbin.env$mrbinTMP$binRegions)
  }
}
#' A function for creating bin numbers.
#'
#' This function calculates numbers of bins from the chosen parameters
#' value.
#' @param
#' @return {None}
#' @export
#' @examples
#' ##Requires data
#' #createBinNumbers()

createBinNumbers<-function(){
   if(mrbin.env$mrbinparam$binMethod=="Special bin list"){
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
#' This function creates bins for the current spectrum.
#' @param
#' @return {None}
#' @export
#' @examples
#' ## Requires data
#' # binSingleNMR()

binSingleNMR<-function(){
 if(!is.null(mrbin.env$mrbinTMP$currentFolder)){
   if(mrbin.env$mrbinparam$dimension=="2D"){#2d spectra
      if(is.null(mrbin.env$mrbinTMP$binNames))          createBinNames()
      #Create index of signals in each bin
      NMRspectrumRownames<-as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrum))
      NMRspectrumColnames<-as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrum))
      #counter<-1
      mrbin.env$mrbinTMP$binTMP<-rep(0,nrow(mrbin.env$mrbinTMP$binRegions))
      names(mrbin.env$mrbinTMP$binTMP)<-rownames(mrbin.env$mrbinTMP$binRegions)
     for(ibinTMP in 1:nrow(mrbin.env$mrbinTMP$binRegions)){
          mrbin.env$mrbinTMP$binTMP[ibinTMP]<-sum(mrbin.env$mrbinTMP$currentSpectrum[
                                         NMRspectrumRownames<=mrbin.env$mrbinTMP$binRegions[ibinTMP,4]&
                                         NMRspectrumRownames>mrbin.env$mrbinTMP$binRegions[ibinTMP,3],
                                         NMRspectrumColnames<=mrbin.env$mrbinTMP$binRegions[ibinTMP,1]&
                                         NMRspectrumColnames>mrbin.env$mrbinTMP$binRegions[ibinTMP,2]])/
                                        (sum(NMRspectrumRownames<=mrbin.env$mrbinTMP$binRegions[ibinTMP,4]&
                                         NMRspectrumRownames>mrbin.env$mrbinTMP$binRegions[ibinTMP,3])*
                                         sum(NMRspectrumColnames<=mrbin.env$mrbinTMP$binRegions[ibinTMP,1]&
                                         NMRspectrumColnames>mrbin.env$mrbinTMP$binRegions[ibinTMP,2]))
     }
   } else {#1D spectra
       if(is.null(mrbin.env$mrbinTMP$binNames)) createBinNames()
       mrbin.env$mrbinTMP$binTMP<-rep(0,nrow(mrbin.env$mrbinTMP$binRegions))
       names(mrbin.env$mrbinTMP$binTMP)<-rownames(mrbin.env$mrbinTMP$binRegions)
       #counter<-1
       NMRspectrumNames<-as.numeric(names(mrbin.env$mrbinTMP$currentSpectrum))
     for(ibinTMP in 1:nrow(mrbin.env$mrbinTMP$binRegions)){
          mrbin.env$mrbinTMP$binTMP[ibinTMP]<-sum(mrbin.env$mrbinTMP$currentSpectrum[
                                         NMRspectrumNames<=mrbin.env$mrbinTMP$binRegions[ibinTMP,1]&
                                         NMRspectrumNames>mrbin.env$mrbinTMP$binRegions[ibinTMP,2]])/
                                        (sum(NMRspectrumNames<=mrbin.env$mrbinTMP$binRegions[ibinTMP,1]&
                                         NMRspectrumNames>mrbin.env$mrbinTMP$binRegions[ibinTMP,2]))
     }
    }
    if(sum(mrbin.env$mrbinTMP$binTMP==0)>0){
        message("Warning: Binning region may be larger than total spectrum size.")
    }
 }
}

#' A function for reading NMR spectra.
#'
#' This function picks the correct NMR reading function, based on vendor.
#' @param
#' @return {None}
#' @export
#' @examples
#' ## Requires data
#' # readNMR()

readNMR<-function(){#Read NMR spectral data
 if(!is.null(mrbin.env$mrbinTMP$currentFolder)){
   if(mrbin.env$mrbinparam$readFromFile){
      if(mrbin.env$mrbinparam$NMRvendor=="Bruker"){
          readBruker()
      }  else {
          stop(paste("No data import function defined for vendor ",mrbin.env$mrbinparam$NMRvendor,".\n",sep=""))
      }
   } else { #From variable, so far only used in example code
         mrbin.env$mrbinTMP$currentSpectrum<-eval(parse(text = mrbin.env$mrbinTMP$currentFolder))#Turn character into variable name
         mrbin.env$mrbinTMP$currentSpectrumName<-mrbin.env$mrbinTMP$currentFolder
   }
  if(mrbin.env$mrbinparam$referenceScaling=="Yes") referenceScaling()
 }
}

#' A function for reading Bruker NMR spectra.
#'
#' This function reads Bruker NMR data.
#' @param
#' @return {None}
#' @export
#' @examples
#' ## Requires folder names
#' # readBruker()

readBruker<-function(){#Read Bruker NMR spectral data
 if(!is.null(mrbin.env$mrbinTMP$currentFolder)){
     datanameDict<-c("1r","2rr")
     names(datanameDict)<-c("1D","2D")
     datanameTmp<-datanameDict[mrbin.env$mrbinparam$dimension]
     spectrum_proc_path<-mrbin.env$mrbinTMP$currentFolder
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
  mrbin.env$mrbinTMP$currentSpectrum<-currentSpectrum
  if(mrbin.env$mrbinparam$useAsNames=="Spectrum titles")    mrbin.env$mrbinTMP$currentSpectrumName<-TITLE
  if(mrbin.env$mrbinparam$useAsNames=="Folder names")    mrbin.env$mrbinTMP$currentSpectrumName<-rev(strsplit(mrbin.env$mrbinTMP$currentFolder,"/")[[1]])[4]
 }
}

#' A function for scaling to the reference area.
#'
#' This function scales NMR data to the reference area.
#' @param
#' @return {None}
#' @export
#' @examples
#' ## Requires data
#' # referenceScaling()

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

#' A function for summing up bins.
#'
#' This function sums up bins. The sums are saved to the middle (median) bin of
#' the original area. All other bins of the area are removed then. This is handy
#' for signals that are know to vary between spectra due to pH or salt content,
#' such as citric acid.
#' @param
#' @return {None}
#' @export
#' @examples
#' ## Requires data
#' # sumBins()

sumBins<-function(){#sum up regions with shifting peaks and remove remaining bins
 if(!is.null(mrbin.env$bins)){
  if(nrow(mrbin.env$mrbinparam$sumBinList)>0){
    for(i in 1:nrow(mrbin.env$mrbinparam$sumBinList)){
       limits<-mrbin.env$mrbinparam$sumBinList[i,]
       if(mrbin.env$mrbinparam$dimension=="1D"){
          TMP<-apply(mrbin.env$mrbinTMP$binRegions[,1:2],1,mean)>limits[2]& apply(mrbin.env$mrbinTMP$binRegions[,1:2],1,mean)<limits[1]
          if(sum(TMP)>0){
              i_TMP<-quantile(x=1:sum(TMP), probs = .5,type=3)#define "middle" bin. This one will be kept
              mrbin.env$bins[,which(TMP)[i_TMP]]<-apply(mrbin.env$bins[,which(TMP)],1,sum)
              if(nrow(mrbin.env$bins)==1){
                rownamesTMP<-rownames(mrbin.env$bins)
                colnamesTMP<-colnames(mrbin.env$bins)[-which(TMP)[-i_TMP]]
                mrbin.env$bins<-matrix(mrbin.env$bins[,-which(TMP)[-i_TMP]],nrow=1)
                rownames(mrbin.env$bins)<-rownamesTMP
                colnames(mrbin.env$bins)<-colnamesTMP
              } else {
                mrbin.env$bins<-mrbin.env$bins[,-which(TMP)[-i_TMP]]
              }
              mrbin.env$mrbinTMP$binRegions<-mrbin.env$mrbinTMP$binRegions[-which(TMP)[-i_TMP],]
          }
       } else {#2D limits=c(4.04,4.08,58,60)
           NMRdataNames<-cbind(apply(mrbin.env$mrbinTMP$binRegions[,3:4],1,mean),apply(mrbin.env$mrbinTMP$binRegions[,1:2],1,mean)) #t(matrix(as.numeric(unlist(strsplit(colnames(mrbin.env$bins),","))),nrow=2))
           TMP<-NMRdataNames[,2]>limits[2]&NMRdataNames[,2]<limits[1]&
                  NMRdataNames[,1]>limits[3]& NMRdataNames[,1]<limits[4]
           if(sum(TMP)>0){
              i_TMP<-quantile(x=1:sum(TMP), probs = .5,type=3)#define "middle" bin. This one will be kept
              mrbin.env$bins[,which(TMP)[i_TMP]]<-apply(mrbin.env$bins[,which(TMP)],1,sum)
              if(nrow(mrbin.env$bins)==1){
                rownamesTMP<-rownames(mrbin.env$bins)
                colnamesTMP<-colnames(mrbin.env$bins)[-which(TMP)[-i_TMP]]
                mrbin.env$bins<-matrix(mrbin.env$bins[,-which(TMP)[-i_TMP]],nrow=1)
                rownames(mrbin.env$bins)<-rownamesTMP
                colnames(mrbin.env$bins)<-colnamesTMP
              } else {
                mrbin.env$bins<-mrbin.env$bins[,-which(TMP)[-i_TMP]]
              }
              mrbin.env$mrbinTMP$binRegions<-mrbin.env$mrbinTMP$binRegions[-which(TMP)[-i_TMP],]
           }
       }
    }
  }
  mrbin.env$mrbinparam$numberOfFeaturesAfterSummingBins<-ncol(mrbin.env$bins)
 }
}

#' A function for removing the solvent region.
#'
#' This function removes the solvent region.
#' @param
#' @return {None}
#' @export
#' @examples
#' ## Requires data
#' # removeSolvent()

removeSolvent<-function(){
 if(!is.null(mrbin.env$bins)){
   if(mrbin.env$mrbinparam$dimension=="1D"){
       NMRdataNames<-apply(mrbin.env$mrbinTMP$binRegions[,1:2],1,mean)
       solventTMP<-NMRdataNames>mrbin.env$mrbinparam$solventRegion[2]&NMRdataNames<mrbin.env$mrbinparam$solventRegion[1]
   }
   if(mrbin.env$mrbinparam$dimension=="2D"){
       NMRdataNames<-cbind(apply(mrbin.env$mrbinTMP$binRegions[,3:4],1,mean),apply(mrbin.env$mrbinTMP$binRegions[,1:2],1,mean))
       solventTMP<-NMRdataNames[,2]>mrbin.env$mrbinparam$solventRegion[2]&NMRdataNames[,2]<mrbin.env$mrbinparam$solventRegion[1]
   }
   if(sum(solventTMP)>0){
      if(nrow(mrbin.env$bins)==1){
        rownamesTMP<-rownames(mrbin.env$bins)
        colnamesTMP<-colnames(mrbin.env$bins)[-which(solventTMP)]
        mrbin.env$bins<-matrix(mrbin.env$bins[,-which(solventTMP)],nrow=1)
        rownames(mrbin.env$bins)<-rownamesTMP
        colnames(mrbin.env$bins)<-colnamesTMP
      } else {
        mrbin.env$bins<-mrbin.env$bins[,-which(solventTMP)]
      }
      mrbin.env$mrbinTMP$binRegions<-mrbin.env$mrbinTMP$binRegions[-which(solventTMP),]
   }
   mrbin.env$mrbinparam$numberOfFeaturesAfterRemovingSolvent<-ncol(mrbin.env$bins)
 }
}

#' A function for removing additional regions.
#'
#' This function removes additional regions. This can be useful when some areas
#' are visibly affected by spectral artifacts.
#' @param
#' @return {None}
#' @export
#' @examples
#' ## Requires data
#' # removeAreas()

removeAreas<-function(){#limits=c(4.75,4.95,-10,160)
 if(!is.null(mrbin.env$bins)){
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
         if(nrow(mrbin.env$bins)==1){
            rownamesTMP<-rownames(mrbin.env$bins)
            colnamesTMP<-colnames(mrbin.env$bins)[-removeTMP]
            mrbin.env$bins<-matrix(mrbin.env$bins[,-removeTMP],nrow=1)
            rownames(mrbin.env$bins)<-rownamesTMP
            colnames(mrbin.env$bins)<-colnamesTMP
         } else {
            mrbin.env$bins<-mrbin.env$bins[,-removeTMP]
         }
         mrbin.env$mrbinTMP$binRegions<-mrbin.env$mrbinTMP$binRegions[-removeTMP,]
     }
  }
  mrbin.env$mrbinparam$numberOfFeaturesAfterRemovingAreas<-ncol(mrbin.env$bins)
 }
}

#' A function for calculating noise levels.
#'
#' This function calculates noise levels.
#' @param
#' @return {None}
#' @export
#' @examples
#' ## Requires data
#' # calculateNoise()

calculateNoise<-function(){
 if(!is.null(mrbin.env$bins)){
  if(mrbin.env$mrbinparam$dimension=="2D"){
       NMRdataNames<-cbind(apply(mrbin.env$mrbinTMP$binRegions[,3:4],1,mean),apply(mrbin.env$mrbinTMP$binRegions[,1:2],1,mean))
  }
  noise_level<-rep(0,nrow(mrbin.env$bins))
  names(noise_level)<-rownames(mrbin.env$bins)
  for(j in 1:nrow(mrbin.env$bins)){
    if(mrbin.env$mrbinparam$dimension=="1D"){
         noise_level[j]<-stats::sd(mrbin.env$bins[j,which(apply(mrbin.env$mrbinTMP$binRegions[,1:2],1,mean)>mrbin.env$mrbinparam$noiseRange1d[1]&
                                        apply(mrbin.env$mrbinTMP$binRegions[,1:2],1,mean)<mrbin.env$mrbinparam$noiseRange1d[2])])
    }
    if(mrbin.env$mrbinparam$dimension=="2D"){
         noise_level[j]<-stats::sd(mrbin.env$bins[j,which(NMRdataNames[,2]>mrbin.env$mrbinparam$noiseRange2d[1]&
                                NMRdataNames[,2]<mrbin.env$mrbinparam$noiseRange2d[2]&
                                NMRdataNames[,1]>mrbin.env$mrbinparam$noiseRange2d[3]&
                                NMRdataNames[,1]<mrbin.env$mrbinparam$noiseRange2d[4])])
    }
  }
  mrbin.env$mrbinparam$noise_level<-noise_level
 }
}

#' A function for removing bins below noise level.
#'
#' This function checks for each bin (column) whether its level is below the
#' individual noise level times the signal-to-noise ratio. If less than the
#' defined threshold level are above noise*SNR, the whole bin is removed.
#' @param
#' @return {None}
#' @export
#' @examples
#' ## Requires data
#' # removeNoise()

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
    calculateNoise()
    for(i in 1:ncol(mrbin.env$bins)){#Keep only bins where at least X spectra are > SNR
          if(sum(mrbin.env$bins[,i]>mrbin.env$mrbinparam$noise_level*SNR)>=minimumNumber){
              colnames_NMRdata_no_noise<-c(colnames_NMRdata_no_noise,i)#colnames(mrbin.env$bins)[i])
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
#' ## Requires data
#' # cropNMR()

cropNMR<-function(plot=FALSE){
 if(!is.null(mrbin.env$bins)){
   if(mrbin.env$mrbinparam$dimension=="2D"){
     selectedCols<-NULL
     for(j in 1:ncol(mrbin.env$bins)){
        coordTmp<-cbind(apply(mrbin.env$mrbinTMP$binRegions[,3:4],1,mean),apply(mrbin.env$mrbinTMP$binRegions[,1:2],1,mean))[j,]#1=C,2=H
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
    if(nrow(mrbin.env$bins)==1){
            rownamesTMP<-rownames(mrbin.env$bins)
            colnamesTMP<-colnames(mrbin.env$bins)[selectedCols]
            mrbin.env$bins<-matrix(mrbin.env$bins[,selectedCols],nrow=1)
            rownames(mrbin.env$bins)<-rownamesTMP
            colnames(mrbin.env$bins)<-colnamesTMP
    } else {
      mrbin.env$bins<-mrbin.env$bins[,selectedCols]
    }
  }
  mrbin.env$mrbinparam$numberOfFeaturesAfterCropping<-ncol(mrbin.env$bins)
 }
}

#' A function for PQN scaling.
#'
#' This function performs PQN scaling. To further exclude unreliable noise, only
#' the most intense signals are used. For 1H-13C HSQC spectra, most of the sugar
#' regions are excluded to avoid a dominating effect of the multiple sugar
#' signals.
#' @param
#' @return {None}
#' @export
#' @examples
#' ## Requires data
#' # PQNScaling()

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
      message("Too few samples to perform PQN normalization.")
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
#' ## Requires data
#' # plotBins()

plotBins<-function(showtitle=TRUE){#Plot bin positions
 if(!is.null(mrbin.env$bins)){
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
        colorTMP<-"green"
      }
      graphics::plot(allBins2[,2],allBins2[,1],ylim=mrbin.env$mrbinparam$binRegion[4:3],xlim=mrbin.env$mrbinparam$binRegion[1:2],
           cex=.5,pch=15,main=titleTMP,col=colorTMP, ask=FALSE)
   }
   if(mrbin.env$mrbinparam$dimension=="1D") {
      colorRampHSQC<-grDevices::colorRamp(c("blue","green","orange","red","red"))
      NMRdata_select<-mrbin.env$bins[mrbin.env$mrbinTMP$currentSpectrumName,]
      NMRdata_select<-NMRdata_select[order(abs(NMRdata_select))]
      if(showtitle) titleTMP<-mrbin.env$mrbinTMP$currentSpectrumName
      if(!showtitle) titleTMP<-""
      graphics::plot(NULL,NULL,ylim=c(0,1),xlim=mrbin.env$mrbinparam$binRegion[1:2],yaxt="n",
           main=titleTMP, ask=FALSE)
      graphics::rect(xleft=mrbin.env$mrbinTMP$binRegions[,1], ybottom=0, xright=mrbin.env$mrbinTMP$binRegions[,2], ytop=1.0, col = "green", border = NA)
      graphics::lines(as.numeric(names(mrbin.env$mrbinTMP$currentSpectrum)),
           mrbin.env$mrbinTMP$currentSpectrum/(sort(mrbin.env$mrbinTMP$currentSpectrum)[ceiling(.9999*length(mrbin.env$mrbinTMP$currentSpectrum))]),
           col="black")

      #graphics::text(5,120,mrbin.env$mrbinTMP$currentSpectrumName,col='gray23'       )
   }
 }
}

#' A function for plotting quality indicators, including PCA plots.
#'
#' This function plots boxplots (bin-wise and sample-wise) as visual quality indicators. It also performs PCA, then plots PC1 and PC2 and loading plots.
#' @param
#' @return {None}
#' @export
#' @examples
#' ## Requires data
#' # plotResults()

plotResults<-function(){
 if(!is.null(mrbin.env$bins)){
    if(is.null(mrbin.env$mrbinparam$Factors)) mrbin.env$mrbinparam$Factors<-factor(rep("Group 0",nrow(mrbin.env$bins)))
    colorPalette<-grDevices::rainbow(length(levels(mrbin.env$mrbinparam$Factors)))
    mrbin.env$mrbinTMP$PCA<-stats::prcomp(mrbin.env$bins)
    oldpar<-graphics::par("mar","mfrow","xpd")
    on.exit(graphics::par(oldpar))
    devAskNewPage(ask = FALSE)
    graphics::par(mfrow=c(2,3),mar=c(3.1,2,1.0,0.5))
    plotBins(showtitle=FALSE)
    graphics::boxplot(mrbin.env$bins,main="",xlab="Bins",ylab="",boxwex=1,cex.lab=.25, ask=FALSE)
    graphics::boxplot(t(mrbin.env$bins),main="",xlab="Samples",ylab="",boxwex=1,cex.lab=.25, ask=FALSE)
    if(nrow(mrbin.env$bins)>1){
      graphics::par(xpd=TRUE,mar=c(4.2,4.1,2.8,0.5))
      graphics::plot(mrbin.env$mrbinTMP$PCA$x[,1],mrbin.env$mrbinTMP$PCA$x[,2],
           pch=as.numeric(mrbin.env$mrbinparam$Factors)+14,
           col=colorPalette[as.numeric(mrbin.env$mrbinparam$Factors)],
           main="PCA",
           xlab=paste("PC1 (",round(100*mrbin.env$mrbinTMP$PCA$sdev[1]/sum(mrbin.env$mrbinTMP$PCA$sdev),1),"%)",sep=""),
           ylab=paste("PC2 (",round(100*mrbin.env$mrbinTMP$PCA$sdev[2]/sum(mrbin.env$mrbinTMP$PCA$sdev),1),"%)",sep="")
           ,cex=.75, ask=FALSE)
      numlevels<-NULL
      for(i in 1:nlevels(mrbin.env$mrbinparam$Factors)) numlevels<-c(numlevels,as.numeric(
                         mrbin.env$mrbinparam$Factors[which(mrbin.env$mrbinparam$Factors==levels(mrbin.env$mrbinparam$Factors)[i])][1]))
      graphics::legend("bottomleft", inset=c(-.075,-.1),#c(-0.1*(max(mrbin.env$mrbinTMP$PCA$x[,1])-min(mrbin.env$mrbinTMP$PCA$x[,1])),-0.9*(max(mrbin.env$mrbinTMP$PCA$x[,2])-min(mrbin.env$mrbinTMP$PCA$x[,2]))),
              legend=levels(mrbin.env$mrbinparam$Factors),
              col=colorPalette[numlevels],
              pch=numlevels+14,
              cex=1,bg="white")
      graphics::text(mrbin.env$mrbinTMP$PCA$x,labels=paste(substr(rownames(mrbin.env$mrbinTMP$PCA$x),1,mrbin.env$mrbinparam$PCAtitlelength)),pos=3,cex=1,
             col=colorPalette[as.numeric(mrbin.env$mrbinparam$Factors)])
      graphics::par(xpd=FALSE)
      graphics::plot(mrbin.env$mrbinTMP$PCA$rotation,pch=16,cex=.75,main="Loadings Plot", ask=FALSE)
      graphics::text(mrbin.env$mrbinTMP$PCA$rotation,labels=substr(rownames(mrbin.env$mrbinTMP$PCA$rotation),1,12),pos=4,cex=1)
      graphics::plot(NULL,NULL,xaxt="n",yaxt="n",ask=FALSE,axes=FALSE,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",main="")
    } else {
       message("Too few samples to perform PCA.")
    }
 }
}

#' A function for plotting NMR spectra.
#'
#' This function plots the current NMR spectrum. To change the plot, use zoom(),
#' zoomIn(), zoomOut(), intPlus(), intMin(), left(), right().
#' For 2D data use additionally: contMin(), contPlus(), up(), down()
#' @param
#' @return {None}
#' @export
#' @examples
#' ## Requires data
#' # plotNMR()

plotNMR<-function(){
 if(!is.null(mrbin.env$mrbinTMP$currentSpectrum)){
   if(is.matrix(mrbin.env$mrbinTMP$currentSpectrum)){#2D spectra
      spectrumTMP<-mrbin.env$mrbinTMP$currentSpectrum[which(as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrum))<
                           mrbin.env$mrbinplot$plotRegion[4]&
                          as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrum))>=mrbin.env$mrbinplot$plotRegion[3]),
                         which(as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrum))>=mrbin.env$mrbinplot$plotRegion[2]&
                         as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrum))<mrbin.env$mrbinplot$plotRegion[1])
                         ]
      if(sum(spectrumTMP<(mrbin.env$mrbinplot$lowestContour*max(mrbin.env$mrbinTMP$currentSpectrum)))>0){
           spectrumTMP[spectrumTMP<(mrbin.env$mrbinplot$lowestContour*max(mrbin.env$mrbinTMP$currentSpectrum))]<-0
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
      if(mrbin.env$mrbinplot$heatmap){
          stats::heatmap(spectrumTMP, Rowv = NA, Colv = NA,col =rev(grDevices::heat.colors(128)),
               scale="none",
               labRow="", labCol="")
      } else {
        options(max.contour.segments=1000)
        graphics::contour(x = -as.numeric(colnames(spectrumTMP)),
          y = -as.numeric(rownames(spectrumTMP)),
          z = t(spectrumTMP)*mrbin.env$mrbinplot$intensityScale,
          levels =  (mrbin.env$mrbinplot$lowestContour+
                        .8*(1:mrbin.env$mrbinplot$nContours-1)/
                        (mrbin.env$mrbinplot$nContours)*(1-
                        mrbin.env$mrbinplot$lowestContour))* max(spectrumTMP)
          ,drawlabels =FALSE
          ,col = rev(grDevices::rainbow(mrbin.env$mrbinplot$nContours))
          ,xaxt="none",yaxt="none"
          ,lwd=1,
          main=mrbin.env$mrbinTMP$currentSpectrumName
          )
        magnitude2<-10^round(log(max(as.numeric(colnames(spectrumTMP)))-
                    min(as.numeric(colnames(spectrumTMP))),base=10))/10
        graphics::axis(1,
             at=(0:100*magnitude2+floor(min(-as.numeric(colnames(spectrumTMP)))/magnitude2)*magnitude2)
             ,labels=-(0:100*magnitude2+floor(min(-as.numeric(colnames(spectrumTMP)))/magnitude2)*magnitude2)
             )
        magnitude1<-10^round(log(max(as.numeric(rownames(spectrumTMP)))-
                     min(as.numeric(rownames(spectrumTMP))),base=10))/10
        graphics::axis(2,
             at=-(0:100*magnitude1+floor(min(as.numeric(rownames(spectrumTMP)))/magnitude1)*magnitude1)
             ,labels=(0:100*magnitude1+floor(min(as.numeric(rownames(spectrumTMP)))/magnitude1)*magnitude1)
             )
      }
   } else {  #1D
      graphics::plot(as.numeric(names(mrbin.env$mrbinTMP$currentSpectrum)),mrbin.env$mrbinTMP$currentSpectrum,
           type="l",xlim=c(mrbin.env$mrbinplot$plotRegion[1],mrbin.env$mrbinplot$plotRegion[2]),
           ylim=c(min(as.numeric(names(mrbin.env$mrbinTMP$currentSpectrum))),
                   max(as.numeric(names(mrbin.env$mrbinTMP$currentSpectrum))))/mrbin.env$mrbinplot$intensityScale,
           xlab="Chemical shift [ppm]",ylab="Intensity")
    }
 }
}

#' A function for changing plotNMR plots.
#'
#' This function increases the intensity of the current NMR spectrum plot.
#' @param
#' @return {None}
#' @export
#' @examples
#' ## Requires data
#' # intPlus()

intPlus<-function(){#increase plot intensity
 if(!is.null(mrbin.env$mrbinTMP$currentSpectrum)){
   mrbin.env$mrbinplot$intensityScale<-mrbin.env$mrbinplot$intensityScale*2
   plotNMR()
 }
}

#' A function for changing plotNMR plots.
#'
#' This function decreases the intensity of the current NMR spectrum plot.
#' @param
#' @return {None}
#' @export
#' @examples
#' ## Requires data
#' # intMin()

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
#' @param
#' @return {None}
#' @export
#' @examples
#' ## Requires data
#' # contPlus()

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
#' @param
#' @return {None}
#' @export
#' @examples
#' ## Requires data
#' # contMin()

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
#' ## Requires data
#' # zoom(left=4.6,right=2,top=20,bottom=100)
#' ## Ask for user input
#' # zoom()

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
#' @param
#' @return {None}
#' @export
#' @examples
#' ## Requires data
#' # zoomIn()

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
#' @param
#' @return {None}
#' @export
#' @examples
#' ## Requires data
#' # zoomOut()

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
#' @param
#' @return {None}
#' @export
#' @examples
#' ## Requires data
#' # left()

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
#' @param
#' @return {None}
#' @export
#' @examples
#' ## Requires data
#' # right()

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
#' @param
#' @return {None}
#' @export
#' @examples
#' ## Requires data
#' # down()

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
#' @param
#' @return {None}
#' @export
#' @examples
#' ## Requires data
#' # up()

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
#' @param
#' @return A list containing all objects from the local package environment.
#' @export
#' @examples
#' ## Requires data
#' # tempList<-getEnv()

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
#' ## Requires data
#' # putToEnv(list(mrbinparam=tempList$mrbinparam))

putToEnv<-function(variableList){
  if(!is.list(variableList)){
    message("Parameter is not a list. No action performed.")
  } else {
   for(i in 1:length(variableList)){
      assign(names(variableList)[i],variableList[[i]],envir =mrbin.env)
   }
  }
}
