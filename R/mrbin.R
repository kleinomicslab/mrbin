#mrbin - Collection of R funtions for analyzing NMR spectral data.
#Written by Matthias Klein, The Ohio State University
#Jan 30 2020
#Package: mrbin
#Title: Binning, Scaling and Normalization of NMR Data
#Version: 0.94
#Authors@R:
#    person(given = "Matthias",
#           family = "Klein",
#           role = c("aut", "cre"),
#           email = "klein.663@osu.edu",
#           comment = c(ORCID = "0000-0001-7455-5381"))
#Description: This package can create bin lists from Bruker
#    NMR spectra and perform basic processing steps.
#    The function negTrafo() replaces negative values by
#    small values in the noise range, preserving the order
#    of values. mrbinGUI() lets the user set parameters,
#    mrbin() performs the follwoing operations: Reading NMR
#    files, creating bins, removing solvent area, removing
#    additional user-defined areas, summing up bins that
#    contain instable peaks such as citric acid, removes
#    noise bins, crops HSQC spectra to the diagonal area,
#    performs PQN scaling, replaces negative values, log
#    transforms and displays a PCA plot. Parameters are
#    saved in a file and can be restored using
#    recreatemrbin().
#License: GPL-3


#' A parameter resetting function
#'
#' This function resets the parameter variables.
#' @param
#' @keywords
#' @export
#' @examples
#' resetEnv()

resetEnv<-function(){

    mrbinbins<<-NULL

    mrbinTMP <<- list(
               mrbinversion="0.94",
               binsRaw=NULL,
               medians=NULL,
               PCA=NULL,
               binTMP=NULL,
               binNames=NULL,
               nbins1=NULL,
               nbins2=NULL,
               nbins=NULL,
               currentFolder=NULL,
               currentSpectrum=NULL,
               currentSpectrumName=NULL
    )
    mrbinparam <<- list(
               dimension="2D",
               binwidth1D=.003,
               binwidth2D=.02,
               binheight=1,
               binRegion=c(9.5,0.5,10,156),
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
               PQNshowHist=F,#show histograms of quotients
               PQNminimumFeatures=40,#Top number of features to proceed
               PCAtitlelength=6,
               createBins="Yes",
               readNMR="Yes",
               useAsNames="Folder names",#Use spectrum titles as sample names or folder names will be used
               NMRvendor="Bruker",#NMR vendor. Currently, only Bruker data is supported.
               mrbinversionTracker=mrbinTMP$mrbinversion,
               outputFileName="mrbinparam.txt",
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
               )
    mrbinplot <<- list(
               lowestContour=.01,
               plotRegion=c(9.5,.5,-1,156),
               intensityScale=1,
               nContours=60,
               heatmap=F
    )
}

#' A function setting the parameters
#'
#' This function guides the user through the set-up of parameters and starts mrbin().
#' @param
#' @keywords
#' @export
#' @examples
#' mrbinGUI()

mrbinGUI<-function(){
   if(!"mrbinparam"%in%ls(envir = .GlobalEnv))    resetEnv()
   selection<-select.list(c("Yes","No"),preselect="Yes",
              title=paste("mrbin ",mrbinTMP$mrbinversion,". Create bin list?",sep=""))
     if(length(selection)==0|selection=="") invisible()
     if(selection=="Yes"){
       selectionRepeat<-select.list(c("Set new","Use current"),preselect="Set new",
                        title="Set parameters or use existing parameters?")
       if(length(selectionRepeat)==0|selectionRepeat=="") stop("User canceled function.")
       if(selectionRepeat=="Set new"){
          dimension<-select.list(c("1D","2D"),preselect=mrbinparam$dimension,
                                 title="1D or 2D spectra?")
          if(length(dimension)==0|dimension=="") stop("User canceled function.")
          mrbinparam$dimension<<-dimension
          if(dimension=="1D") dimlength<-2
          if(dimension=="2D") dimlength<-4
          adjRegion<-select.list(c(paste(mrbinparam$binRegion[1:dimlength],collapse=","),"Change..."),
                     preselect=paste(mrbinparam$binRegion,collapse=" "),title ="Bin region [ppm]: ")
          if(length(adjRegion)==0|adjRegion=="") stop("User canceled function.")
          if(adjRegion=="Change..."){
            regionTMP<-readline(prompt=paste("New left border, press enter to keep ",
                      mrbinparam$binRegion[1],": ",sep=""))
            if(!regionTMP=="") mrbinparam$binRegion[1]<<-as.numeric(regionTMP)
            regionTMP<-readline(prompt=paste("New right border, press enter to keep ",
                      mrbinparam$binRegion[2],": ",sep=""))
            if(!regionTMP=="") mrbinparam$binRegion[2]<<-as.numeric(regionTMP)
              if(mrbinparam$dimension=="2D"){
                regionTMP<-readline(prompt=paste("New top border, press enter to keep ",
                          mrbinparam$binRegion[3],": ",sep=""))
                if(!regionTMP=="") mrbinparam$binRegion[3]<<-as.numeric(regionTMP)
                regionTMP<-readline(prompt=paste("New bottom border, press enter to keep ",
                          mrbinparam$binRegion[4],": ",sep=""))
                if(!regionTMP=="") mrbinparam$binRegion[4]<<-as.numeric(regionTMP)
              }
          }
          if(mrbinparam$dimension=="1D"){
              widthAdjust<-select.list(c(mrbinparam$binwidth1D,"Change..."),preselect=
                           as.character(mrbinparam$binwidth1D),title ="Bin width [ppm]: ")
              if(length(widthAdjust)==0|widthAdjust=="") stop("User canceled function.")
              if(widthAdjust=="Change..."){
                   widthTMP<-readline(prompt=paste("New 1D bin width, press enter to keep ",
                             mrbinparam$binwidth1D,": ",sep=""))
                   if(!widthTMP=="") mrbinparam$binwidth1D<<-as.numeric(widthTMP)
              }
          }
          if(mrbinparam$dimension=="2D"){
              widthAdjust<-select.list(c(mrbinparam$binwidth2D,"Change..."),preselect=
                           as.character(mrbinparam$binwidth2D),title ="Bin width [ppm]: ")
              if(length(widthAdjust)==0|widthAdjust=="") stop("User canceled function.")
              if(widthAdjust=="Change..."){
                   widthTMP<-readline(prompt=paste("New 2D bin width, press enter to keep ",
                             mrbinparam$binwidth2D,": ",sep=""))
                   if(!widthTMP=="") mrbinparam$binwidth2D<<-as.numeric(widthTMP)
              }
              heightAdjust<-select.list(c(mrbinparam$binheight,"Change..."),preselect=
                            as.character(mrbinparam$binheight),title ="Bin height [ppm]: ")
              if(length(heightAdjust)==0|heightAdjust=="") stop("User canceled function.")
              if(heightAdjust=="Change..."){
                   heightTMP<-readline(prompt=paste("New 2D bin height, press enter to keep ",
                              mrbinparam$binheight,": ",sep=""))
                   if(!heightTMP=="") mrbinparam$binheight<<-as.numeric(heightTMP)
              }
          }
          referenceScaling<-select.list(c("Yes","No"),preselect=mrbinparam$referenceScaling,
                                        title = "Scale to reference signal?")
          if(length(referenceScaling)==0|referenceScaling=="") stop("User canceled function.")
          mrbinparam$referenceScaling<<-referenceScaling
          if(mrbinparam$referenceScaling=="Yes"){
            if(mrbinparam$dimension=="1D"){
              adjRegion<-select.list(c(paste(mrbinparam$reference1D,collapse=" "),"Change..."),
                         preselect=paste(mrbinparam$reference1D,collapse=" "),title ="Reference region [ppm]: ")
              if(length(adjRegion)==0|adjRegion=="") stop("User canceled function.")
              if(adjRegion=="Change..."){
                regionTMP<-readline(prompt=paste("New left border, press enter to keep ",
                          mrbinparam$reference1D[1],": ",sep=""))
                if(!regionTMP=="") mrbinparam$reference1D[1]<<-as.numeric(regionTMP)
                regionTMP<-readline(prompt=paste("New right border, press enter to keep ",
                          mrbinparam$reference1D[2],": ",sep=""))
                if(!regionTMP=="") mrbinparam$reference1D[2]<<-as.numeric(regionTMP)
              }
            }
            if(mrbinparam$dimension=="2D"){
              adjRegion<-select.list(c(paste(mrbinparam$reference2D,collapse=" "),"Change..."),
                         preselect=paste(mrbinparam$reference2D,collapse=" "),title ="Reference region [ppm]: ")
              if(length(adjRegion)==0|adjRegion=="") stop("User canceled function.")
              if(adjRegion=="Change..."){
              regionTMP<-readline(prompt=paste("New left border, press enter to keep ",
                        mrbinparam$reference2D[1],": ",sep=""))
              if(!regionTMP=="") mrbinparam$reference2D[1]<<-as.numeric(regionTMP)
              regionTMP<-readline(prompt=paste("New right border, press enter to keep ",
                        mrbinparam$reference2D[2],": ",sep=""))
              if(!regionTMP=="") mrbinparam$reference2D[2]<<-as.numeric(regionTMP)
              regionTMP<-readline(prompt=paste("New top border, press enter to keep ",
                        mrbinparam$reference2D[3],": ",sep=""))
              if(!regionTMP=="") mrbinparam$reference2D[3]<<-as.numeric(regionTMP)
              regionTMP<-readline(prompt=paste("New bottom border, press enter to keep ",
                        mrbinparam$reference2D[4],": ",sep=""))
              if(!regionTMP=="") mrbinparam$reference2D[4]<<-as.numeric(regionTMP)
              }
            }
          }
          removeSolvent<-select.list(c("Yes","No"),preselect=mrbinparam$removeSolvent,
                                   title = "Remove solvent area?")
          if(length(removeSolvent)==0|removeSolvent=="") stop("User canceled function.")
          mrbinparam$removeSolvent<<-removeSolvent
          if(mrbinparam$removeSolvent=="Yes"){
              adjRegion<-select.list(c(paste(mrbinparam$solventRegion,collapse=" "),"Change..."),
                         preselect=paste(mrbinparam$solventRegion,collapse=" "),title ="Solvent region [ppm]: ")
              if(length(adjRegion)==0|adjRegion=="") stop("User canceled function.")
              if(adjRegion=="Change..."){
                regionTMP<-readline(prompt=paste("New left border, press enter to keep ",
                          mrbinparam$solventRegion[1],": ",sep=""))
                if(!regionTMP=="") mrbinparam$solventRegion[1]<<-as.numeric(regionTMP)
                regionTMP<-readline(prompt=paste("New right border, press enter to keep ",
                          mrbinparam$solventRegion[2],": ",sep=""))
                if(!regionTMP=="") mrbinparam$solventRegion[2]<<-as.numeric(regionTMP)
              }
          }
          removeAreas<-select.list(c("Yes","No"),preselect=mrbinparam$removeAreas,
                                   title = "Remove additional areas?")
          if(length(removeAreas)==0|removeAreas=="") stop("User canceled function.")
          mrbinparam$removeAreas<<-removeAreas
          if(mrbinparam$removeAreas=="Yes"){
            addAreasFlag<-T
            if(nrow(mrbinparam$removeAreaList)>0){
                addAreasFlag<-F
                removeAreaListTMP<-select.list(c("Keep","New","Add to existing list"),preselect="Keep",
                                   title = "Use previous area list or define new?")
                if(length(removeAreaListTMP)==0|removeAreaListTMP=="") stop("User canceled function.")
                if(!removeAreaListTMP=="Keep"){
                  addAreasFlag<-T
                  if(removeAreaListTMP=="New"){
                    mrbinparam$removeAreaList<<-matrix(ncol=4,nrow=0,dimnames=list(NULL,c("left","right","top","bottom")))
                  }
                }
            }
            iaddAreas<-nrow(mrbinparam$removeAreaList)+1
            while(addAreasFlag){
              mrbinparam$removeAreaList<<-rbind(mrbinparam$removeAreaList,c(0,0,0,0))
              mrbinparam$removeAreaList[iaddAreas,1]<<-as.numeric(readline(prompt="Left border: "))
              mrbinparam$removeAreaList[iaddAreas,2]<<-as.numeric(readline(prompt="Right border: "))
              if(mrbinparam$dimension=="2D"){
                mrbinparam$removeAreaList[iaddAreas,3]<<-as.numeric(readline(prompt="Top border: "))
                mrbinparam$removeAreaList[iaddAreas,4]<<-as.numeric(readline(prompt="Bottom border: "))
              }
              iaddAreas<-iaddAreas+1
              keepAdding<-select.list(c("No","Yes"),preselect="No",multiple=F,title = "Add more areas?")
              if(length(keepAdding)==0|keepAdding=="") stop("User canceled function.")
              if(keepAdding=="No")  addAreasFlag<-F
            }
          }
          sumBins<-select.list(c("Yes","No"),preselect=mrbinparam$sumBins,
                               title = "Sum bins of unstable peaks?")
          if(length(sumBins)==0|sumBins=="") stop("User canceled function.")
          mrbinparam$sumBins<<-sumBins
          if(mrbinparam$sumBins=="Yes"){
              addAreasFlag<-T
              if(nrow(mrbinparam$sumBinList)>0){
                  addAreasFlag<-F
                  sumBinListTMP<-select.list(c("Keep","New","Add to existing list"),preselect="Keep",
                                 title = "Use previous area list or define new?")
                  if(length(sumBinListTMP)==0|sumBinListTMP=="") stop("User canceled function.")
                  if(!sumBinListTMP=="Keep"){
                    addAreasFlag<-T
                    if(!sumBinListTMP=="New"){
                        mrbinparam$sumBinList<<-matrix(ncol=4,nrow=0,dimnames=list(NULL,c("left","right","top","bottom")))
                    }
                  }
              }
              iaddAreas<-nrow(mrbinparam$sumBinList)+1
              while(addAreasFlag){
                mrbinparam$sumBinList<<-rbind(mrbinparam$sumBinList,c(0,0,0,0))
                mrbinparam$sumBinList[iaddAreas,1]<<-as.numeric(readline(prompt="Left border: "))
                mrbinparam$sumBinList[iaddAreas,2]<<-as.numeric(readline(prompt="Right border: "))
                if(mrbinparam$dimension=="2D"){
                    mrbinparam$sumBinList[iaddAreas,3]<<-as.numeric(readline(prompt="Top border: "))
                    mrbinparam$sumBinList[iaddAreas,4]<<-as.numeric(readline(prompt="Bottom border: "))
                }
                iaddAreas<-iaddAreas+1
                keepAdding<-select.list(c("No","Yes"),preselect="No",title = "Add more areas?")
                if(length(keepAdding)==0|keepAdding=="") stop("User canceled function.")
                if(keepAdding=="No")  addAreasFlag<-F
              }
          }
          noiseRemoval<-select.list(c("Yes","No"),
                                   preselect=mrbinparam$noiseRemoval,title="Remove noise?")
          if(length(noiseRemoval)==0|noiseRemoval==""){
           cat("test")
           stop("User canceled function.")
          }
          mrbinparam$noiseRemoval<<-noiseRemoval
          if(mrbinparam$noiseRemoval=="Yes"){
            if(mrbinparam$dimension=="1D"){
              SNRTMP<-select.list(c(mrbinparam$signal_to_noise1D,"Change..."),preselect=
                      as.character(mrbinparam$signal_to_noise1D),title="Signal-to-noise ratio (SNR):")
              if(length(SNRTMP)==0|SNRTMP=="") stop("User canceled function.")
              if(SNRTMP=="Change...") SNRTMP<-readline(prompt=paste(
                                   "New 1D signal to noise ratio, press enter to keep ",
                                   mrbinparam$signal_to_noise1D,": ",sep=""))
              if(!SNRTMP=="") mrbinparam$signal_to_noise1D<<-as.numeric(SNRTMP)
            }
            if(mrbinparam$dimension=="2D"){
              SNRTMP<-select.list(c(mrbinparam$signal_to_noise2D,"Change..."),preselect=
                      as.character(mrbinparam$signal_to_noise2D),title="Signal-to-noise ratio (SNR):")
              if(length(SNRTMP)==0|SNRTMP=="") stop("User canceled function.")
              if(SNRTMP=="Change...") SNRTMP<-readline(prompt=paste(
                                   "New 2D signal to noise ratio, press enter to keep ",
                                   mrbinparam$signal_to_noise2D,": ",sep=""))
              if(!SNRTMP=="") mrbinparam$signal_to_noise2D<<-as.numeric(SNRTMP)
            }
             noiseTMP<-select.list(unique(c(as.character(mrbinparam$noiseThreshold),"0.2","0.75","0.05","Change...")),
                                   preselect=as.character(mrbinparam$noiseThreshold),title="Minimum ratio > SNR")
             if(length(noiseTMP)==0|noiseTMP=="") stop("User canceled function.")
             if(noiseTMP=="Change..."){
                  noiseTMP<-readline(prompt=paste("New noise threshold, press enter to keep ",
                            mrbinparam$noiseThreshold,": ",sep=""))
                  if(!noiseTMP=="") mrbinparam$noiseThreshold<<-as.numeric(noiseTMP)
             } else {
              mrbinparam$noiseThreshold<<-as.numeric(noiseTMP)
             }
          }
          if(mrbinparam$dimension=="2D"){
               cropHSQC<-select.list(c("Yes","No"),
                                                               preselect=mrbinparam$cropHSQC,
                                                               title="Crop H-C HSQCs?")
               if(length(cropHSQC)==0|cropHSQC=="") stop("User canceled function.")
               mrbinparam$cropHSQC<<-cropHSQC
          }
          fixNegatives<-select.list(c("Yes","No"),preselect=mrbinparam$fixNegatives,
                                    title="Replace negative values?")
          if(length(fixNegatives)==0|fixNegatives=="") stop("User canceled function.")
          mrbinparam$fixNegatives<<-fixNegatives
          #Log scaling?
          logTrafo<-select.list(c("Yes","No"),preselect=mrbinparam$logTrafo,
                                title="Log transformation?")
          if(length(logTrafo)==0|logTrafo=="") stop("User canceled function.")
          mrbinparam$logTrafo<<- logTrafo
          #PQN scaling?
          PQNScaling<-select.list(c("Yes","No"),preselect=mrbinparam$PQNScaling,
                                  title = "PQN normalization?")
          if(length(PQNScaling)==0|PQNScaling=="") stop("User canceled function.")
          mrbinparam$PQNScaling<<-PQNScaling
          #PCA
          mrbinparam$PCA<<-select.list(c("Yes","No"),preselect = mrbinparam$PCA,
                            title = "Create PCA plot?")
          #Select folders
          if(!is.null(mrbinparam$NMRfolders)){
              selectionFolders<-select.list(c("Yes","No"),preselect="Yes",
                                title="Use previous spectra list?")
              if(length(selectionFolders)==0|selectionFolders=="") stop("User canceled function.")
              if(selectionFolders=="No")  selectFolders()
          } else {
              selectFolders()
          }
          #Define groups
          defineGroups<-select.list(c("Yes","No"),preselect=mrbinparam$defineGroups,
                                    title = "Define group members?")
          if(length(defineGroups)==0|defineGroups=="") stop("User canceled function.")
          mrbinparam$defineGroups<<-defineGroups
          if(mrbinparam$defineGroups=="Yes"){
              if(!is.null(mrbinparam$Factors)){
                 if(length(mrbinparam$Factors)==length(mrbinparam$NMRfolders)){
                   selectionFactors<-select.list(c("Yes","No"),preselect="Yes",
                                     title = "Use previous factor list?")
                   if(length(selectionFactors)==0|selectionFactors=="") stop("User canceled function.")
                   if(selectionFactors=="No")  setFactors()
                 }
              } else {
                   setFactors()
              }
          }
          #Define sample names
          useAsNames<-select.list(c("Folder names","Spectrum titles"),
                                    preselect=mrbinparam$useAsNames,
                                    title = "Create sample names from")
          if(length(useAsNames)==0|useAsNames=="") stop("User canceled function.")
          mrbinparam$useAsNames<<-useAsNames
          #Define file name
          filenameTMP<-select.list(c(paste("mrbinparam_",gsub(":","-",gsub(" ","_",Sys.time())),".txt",
                       sep=""),"Change..."),#preselect=paste(mrbinparam$solventRegion,collapse=" "),
                       title ="Output file name: ")
          if(length(filenameTMP)==0|filenameTMP=="") stop("User canceled function.")
          if(filenameTMP=="Change..."){
            filenameTMP<-readline(prompt=paste("New file name, press enter to use ",
                      paste("mrbinparam_",gsub(":","-",gsub(" ","_",Sys.time())),".txt",sep=""),": \n",sep=""))
            if(!filenameTMP=="") mrbinparam$outputFileName<<-filenameTMP
            if(filenameTMP=="") mrbinparam$outputFileName<<-paste("mrbinparam_",gsub(":","-",gsub(" ","_",Sys.time())),".txt",
                       sep="")


          } else {
              mrbinparam$outputFileName<<-filenameTMP
          }
   }
   startmrbin<-select.list(c("Start binning now","I'll do it later, using mrbin()"),
              preselect="Start binning now",title = "You're all set. Start binning?")
   if(length(startmrbin)==0|startmrbin=="") stop("User canceled function.")
   if(startmrbin=="Start binning now")  mrbin()
  }
}

#' A function performing all data read and processing steps.
#'
#' This function reads parameters from the global variable mrbinparam and
#' performs the follwoing operations:
#' Reading NMR files, creating bins, removing solvent area, removing additional
#' user-defined areas, summing up bins that contain instable peaks such as
#' citric acid, removes noise bins, crops HSQC spectra to the diagonal area,
#' performs PQN scaling, replaces negative values, log transforms and displays a
#' PCA plot. Parameters are then saved in a text file. These can be recreated
#' using recreatemrbin().
#' @param
#' @keywords
#' @export
#' @examples
#' mrbin()

mrbin<-function(){
    if(!"mrbinparam"%in%ls(envir = .GlobalEnv))    resetEnv()
    if(mrbinparam$createBins=="Yes") binMultiNMR()
    if(mrbinparam$removeSolvent=="Yes") removeSolvent()
    if(mrbinparam$removeAreas=="Yes") removeAreas()
    if(mrbinparam$sumBins=="Yes") sumBins()
    if(mrbinparam$noiseRemoval=="Yes") removeNoise()
    if(mrbinparam$cropHSQC=="Yes"&mrbinparam$dimension=="2D") cropNMR()
    if(mrbinparam$PQNScaling=="Yes") PQNScaling()
    if(mrbinparam$fixNegatives=="Yes") negTrafo()
    if(mrbinparam$logTrafo=="Yes") logTrafo()
    if(mrbinparam$PCA=="Yes") PCA()
    dput(mrbinparam, file = mrbinparam$outputFileName)
    write.csv(mrbinbins, file = "mrbinbins.csv")
}

#' A function recreating parameters from previous runs.
#'
#' This function reads parameters from a text file taht was created during a
#' previous run or mrbin(). After reading, the data can be recreated using
#' mrbin(). File names in mrbin$param might need to be updated.
#' using recreatemrbin().
#' @param file File path/name of the mrbin parameter file to be loaded
#' @keywords
#' @export
#' @examples
#' recreatemrbin()

recreatemrbin<-function(file){
    if(!"mrbinparam"%in%ls(envir = .GlobalEnv))    resetEnv()
    mrbinparam<-dget(file)
    if(!mrbinTMP$mrbinversion==mrbinparam$mrbinversion){
       cat(paste("Imported file was created using the older mrbin version ",mrbinparam$mrbinversion,
           ". For exact reproduction of results, please get the old version at: kleinomicslab.com\n",
           sep=""))
    } else {
        cat(paste("Parameter file was succesfully loaded. Most likely, folder names will need to be updated.\n",
           "Please update the object mrbinparam$NMRfolders if necessary. To recreate the data, enter mrbin().\n",
           sep=""))
    }
}

#' A function replacing negative values.
#'
#' This function replaces (column-wise) negative values by a small positive
#' number. The number is calculated as a linear transformation to the range of
#' the lowest positive number to 0,01*the lowest positive number (of this
#' column). Ranks stay unchanged. Positive numbers are not altered.
#' If sample-wise noise levels are avaliable, the median nosie level of samples
#' with negative values is calculated and replaces the lowest positive number in
#' case it is smaller. If no noise data is available, the 1% percentile of all
#' positive values in the data set is used as an estimate.
#' It is recommended to us this function AFTER noise removal and other data
#' clean-up methods, as it may alter (reduce) the noise level.
#' If no NMR data and noise levels are provided as arguments, the funtion will
#' use NMR data and noise levels from the global variables mrbinbins and
#' mrbinTMP.
#' @param NMRdata A matrix containing NMR data. Columns=frequencies,rows=samples
#' @param noiseLevels A vector
#' @return NMRdata A matrix containing NMR data without negative values.
#' @keywords
#' @export
#' @examples negTrafo()
#' @examples negTrafo(NMRdataMatrix,noiseLevelVector)
#' negTrafo()

negTrafo<-function(NMRdata=NULL,noiseLevels=NULL){
     if(is.null(NMRdata)){
       if("mrbinbins"%in%ls(envir=.GlobalEnv)&"mrbinparam"%in%ls(envir=.GlobalEnv)){
         NMRdata <- mrbinbins
         noiseLevels <- mrbinparam$noise_level
       }
     }
     if(is.null(noiseLevels)){
             noiseTMP<-quantile(x=NMRdata[NMRdata>0], probs = .01,type=3)#min(NMRdata[NMRdata>0])
     }
     for(i in 1:ncol(NMRdata)){
        negatives<-NMRdata[,i]<=0
        if(!is.null(noiseLevels)){#If noise levels are available, restrict range to below noise
             noiseTMP<-median(mrbinparam$noise_level[negatives])
        }
        if(sum(negatives)>0){
            minTMP<-min(NMRdata[negatives,i])#select lowest bin
            maxTMP<-min(noiseTMP,min(NMRdata[!negatives,i]))#select lowest bin above 0

            #mrbinbins[negatives,i]<<-(NMRdata[negatives,i]+(maxTMP-minTMP))/
            #                              (maxTMP-minTMP)*(maxTMP*.99)+maxTMP*.01
            NMRdata[negatives,i]<-(NMRdata[negatives,i]+(maxTMP-minTMP))/
                                          (maxTMP-minTMP)*(maxTMP*.99)+maxTMP*.01
        }
     }
     if("mrbinbins"%in%ls(envir=.GlobalEnv)&"mrbinparam"%in%ls(envir=.GlobalEnv)){
          mrbinbins<<-NMRdata
     }
     return(NMRdata)
}

#' A function for log transforming data.
#'
#' This function simply log transforms. Will not work with negative data.
#' @param
#' @keywords
#' @export
#' @examples
#' logTrafo()

logTrafo<-function(){
     if(sum(mrbinbins<=0)>0){
         stop("Log transform does not work with negative values.")
     } else {
         mrbinbins<<-log(mrbinbins)
     }
}

#' A function for removing a spectrum.
#'
#' This function lets the user pick spectra from a list for removal from data
#' analysis.
#' @param
#' @keywords
#' @export
#' @examples
#' removeSpectrum()

removeSpectrum<-function(){
    listTMP<-select.list(rownames(mrbinbins),preselect = NULL, multiple = T,title ="Select spectra to be removed")
    if(length(listTMP)>0){
         mrbinbins<<-mrbinbins[-which(rownames(mrbinbins)%in%listTMP)]
    }
}

#' A function for setting group members.
#'
#' This function lets the user pick samples from a list to assign them to
#' groups.
#' @param
#' @keywords
#' @export
#' @examples
#' setFactors()

setFactors<-function(){
   if(length(mrbinparam$NMRfolders)>0){
      Factors<-rep("Group 0",length(mrbinparam$NMRfolders))
      names(Factors)<-mrbinparam$NMRfolders
      flag<-T
      i<-0
      while(flag){
        i<-i+1
        listTMP<-select.list(names(Factors),preselect = NULL, multiple = T,title ="Please select group members")
        groupNameTMP<-select.list(c(paste("Group",i),"Enter new name"),preselect=paste("Group",i),
                  multiple=F,title ="Group name?")
        if(groupNameTMP=="Enter new name"){
          groupNameTMP<-readline(prompt=paste("New group name (enter to use \"Group ",i,"\"): ",sep=""))
          if(groupNameTMP=="") groupNameTMP<-paste("Group ",i,sep="")
        }
        Factors[listTMP]<-groupNameTMP
        select<-select.list(c("Yes","No"),preselect = "No",multiple = F,
            title = paste("Define more groups?",sep=""))
        if(select=="No") flag<-F
      }
      Factors<-as.factor(Factors)
      mrbinparam$Factors <<- Factors
   }
}

#' A function for selecting NMR data folders.
#'
#' This function selects the correct folder selection function for the vendor
#' (currently only Bruker).
#' @param
#' @keywords
#' @export
#' @examples
#' selectFolders()

selectFolders<-function(){#Select NMR spectral folders
      if(mrbinparam$NMRvendor=="Bruker"){
          selectBrukerFolders()
      }  else {
          stop(paste("No folder selection function defined for vendor ",mrbinparam$NMRvendor,".\n",sep=""))
      }
}

#' A function for selecting Bruker NMR data folders.
#'
#' This function lets the user set NMR data folders (for Bruker data).
#' @param
#' @keywords
#' @export
#' @examples
#' selectBrukerFolders()

selectBrukerFolders<-function(){#Select Bruker NMR spectral folders
   NMRfolders<-NULL
   parentFolder<-getwd()
   datanameDict<-c("1r","2rr")
   names(datanameDict)<-c("1D","2D")
   datanameTmp<-datanameDict[mrbinparam$dimension]
   singleFolderFlag<-F
   enterFolders<-select.list(c("Browse...","Enter parent folder path manually"),
                 preselect="Browse...",title="Set NMR parent folder:")
   if(enterFolders=="Enter parent folder path manually"){
      enterFoldersTMP<-readline(prompt="Enter parent folder path, press enter to browse instead: ")
      if(!enterFoldersTMP==""){
        parentFolder<-gsub('\\\\',"/",enterFoldersTMP)
        selectFolders<-parentFolder
        selectList<-parentFolder
        singleFolderFlag<-T
      }
   }
   selectFlag<-0
   while(selectFlag<1){
     folderListFull<-list.dirs(path=parentFolder,recursive = F,full.names=T)
     folderList<-list.dirs(path=parentFolder,recursive = F,full.names=F)
     if(length(strsplit(parentFolder,split="/")[[1]])>1){
        folderListTMP<-c("..                                                                            ",folderList)
     } else {
        folderListTMP<-folderList
     }
     if(!singleFolderFlag){
         selectFolders<-select.list(folderListTMP,preselect=NULL,multiple=T,
                          title = paste("Browse to NMR parent folder:",parentFolder))
     }
     if(!length(selectFolders)==1) {
        selectFolders<-parentFolder
        selectList<-parentFolder
        singleFolderFlag<-T
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
      #if(length(selectFolders)>1){
      #    cat("\n Please select only single folders in this step to browse to the NMR parent folder.\n")
      #}
      if(singleFolderFlag){
         cat("Looking for NMR data in folder...\n")
         singleFolderFlag<-F
         spectrum_path_list<-NULL
         subdirsTmp0<-list.dirs(path = selectList,recursive = F,full.names=F)
         if(length(subdirsTmp0)>0){
           for(isubDirs in 1:length(subdirsTmp0)){#Look for Bruker NMR data in folder.
             subdirsTmp<-list.dirs(path = paste(selectList,"/",subdirsTmp0[isubDirs],sep=""),recursive = F,full.names=F)
             if(length(subdirsTmp)>0){
               if(sum(sapply(as.numeric(subdirsTmp),is.numeric))>0){ #Look for folders "1", "2" etc (EXPNO)
                   spectrum_path2<-paste(selectList,"/",subdirsTmp0[isubDirs],"/",subdirsTmp[which(sapply(as.numeric(subdirsTmp),is.numeric))],sep="")
                   for(i in spectrum_path2){
                      if(dir.exists(paste(i,"/pdata",sep=""))){
                          subdirsTmp2<-list.dirs(path = paste(i,"/pdata",sep=""),recursive = F,full.names=F)
                          if(sum(sapply(as.numeric(subdirsTmp2),is.numeric))>0){ #Look for folders "1", "2" etc (PROCNO)
                                spectrum_path3<-paste(i,"/pdata/",subdirsTmp2[which(sapply(as.numeric(subdirsTmp2),is.numeric))],sep="")
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
            cat(paste("No Bruker NMR data found in folder.\n"))
            noDataFoundExit<-select.list(c("Continue brwosing","Exit"),
                            preselect="Continue brwosing",
                            title = "No NMR data found in folder.")
            if(noDataFoundExit=="Continue browsing"){
                  selectFlag<-0
            } else {
                  selectFlag<-1
            }
         } else {
           NMRfolders<-c(NMRfolders,
                     select.list(spectrum_proc_path,preselect=NULL,multiple=T,
                     title = "Select data sets"))
           yesorno<-select.list(c("No","Yes"),preselect="No",multiple=F,title="Add more spectra?")
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
   mrbinparam$NMRfolders<<-NMRfolders
}

#' A function for binning multiple NMR spectra.
#'
#' This function creates bins for each spectrum in mrbinparam$NMRfolders and
#' saves the bins to mrbinbins.
#' @param
#' @keywords
#' @export
#' @examples
#' binMultiNMR()

binMultiNMR<-function(){
    mrbinbins<<-NULL
    mrbinTMP$binNames<<-NULL
    mrbinTMP$binTMP<<-NULL
    #Open and bin all spectra
    if(mrbinparam$readNMR=="Yes"){
      mrbinbinsRaw<<-NULL
      for(i in 1:length(mrbinparam$NMRfolders)){
          cat("Binning spectrum",i,"\n")
          flush.console()
          mrbinTMP$currentFolder<<-mrbinparam$NMRfolders[i]
          readNMR()
          mrbinTMP$binTMP<<-NULL
          binSingleNMR()
          if(is.null(mrbinbinsRaw)){
              mrbinbinsRaw<<-matrix(rep(0,length(mrbinTMP$binTMP)*length(mrbinparam$NMRfolders)),nrow=length(mrbinparam$NMRfolders))
              colnames(mrbinbinsRaw)<<-names(mrbinTMP$binTMP)
              rownames(mrbinbinsRaw)<<-1:length(mrbinparam$NMRfolders)
          }
          if(is.null(mrbinTMP$binTMP)){
              stop("Error while binning spectrum.")
          } else {
              mrbinbinsRaw[i,]<<-mrbinTMP$binTMP
              rownames(mrbinbinsRaw)[i]<<-mrbinTMP$currentSpectrumName#rev(strsplit(mrbinTMP$currentFolder, "/")[[1]])[1]
          }
      }
    }
    mrbinbins<<-mrbinbinsRaw
    mrbinparam$numberOfFeaturesRaw<<-ncol(mrbinbins)
}

#' A function for creating bin titles.
#'
#' This function creates titles for the bins to represent their average ppm
#' value.
#' @param
#' @keywords
#' @export
#' @examples
#' createBinNames()

createBinNames<-function(){
   if(mrbinparam$dimension=="2D"){
      mrbinTMP$nbins2<<-ceiling((mrbinparam$binRegion[1]-mrbinparam$binRegion[2])/mrbinparam$binwidth2D)
      mrbinTMP$nbins1<<-ceiling((mrbinparam$binRegion[4]-mrbinparam$binRegion[3])/mrbinparam$binheight)
      mrbinTMP$nbins<<-mrbinTMP$nbins2*mrbinTMP$nbins1
      binNames2<-mrbinparam$binRegion[1]+mrbinparam$binwidth2D/2-(1:mrbinTMP$nbins2)*mrbinparam$binwidth2D
      binNames1<-mrbinparam$binRegion[4]+mrbinparam$binheight/2-(1:mrbinTMP$nbins1)*mrbinparam$binheight
      mrbinTMP$binNames<<-paste(sort(rep(binNames1,mrbinTMP$nbins2),decreasing=T),
                         rep(binNames2,mrbinTMP$nbins1),sep=",")
  }
  if(mrbinparam$dimension=="1D"){
      mrbinTMP$nbins<<-ceiling((mrbinparam$binRegion[1]-mrbinparam$binRegion[2])/mrbinparam$binwidth1D)
      mrbinTMP$binNames<<-mrbinparam$binRegion[1]+mrbinparam$binwidth1D/2-(1:mrbinTMP$nbins)*mrbinparam$binwidth1D
  }
}

#' A function for binning a sinlge NMR spectrum.
#'
#' This function creates bins for the current spectrum.
#' @param
#' @keywords
#' @export
#' @examples
#' binSingleNMR()

binSingleNMR<-function(){
   if(mrbinparam$dimension=="2D"){#2d spectra
      if(is.null(mrbinTMP$binNames))          createBinNames()
      #Create index of signals in each bin
      NMRspectrumRownames<-as.numeric(rownames(mrbinTMP$currentSpectrum))
      NMRspectrumColnames<-as.numeric(colnames(mrbinTMP$currentSpectrum))
      counter<-1
      mrbinTMP$binTMP<<-rep(0,length(mrbinTMP$binNames))
      names(mrbinTMP$binTMP)<<-mrbinTMP$binNames
      which_j<-matrix(rep(F,nrow(mrbinTMP$currentSpectrum)*mrbinTMP$nbins1),ncol=mrbinTMP$nbins1)
      which_i<-matrix(rep(F,ncol(mrbinTMP$currentSpectrum)*mrbinTMP$nbins2),nrow=mrbinTMP$nbins2)
      for(j in 1:mrbinTMP$nbins1){
           which_j[,j]<-(NMRspectrumRownames>(mrbinparam$binRegion[4]-j*mrbinparam$binheight)&
                          NMRspectrumRownames<=(mrbinparam$binRegion[4]-(j-1)*mrbinparam$binheight))
      }
      for(i in 1:mrbinTMP$nbins2){
           which_i[i,]<-(NMRspectrumColnames>(mrbinparam$binRegion[1]-i*mrbinparam$binwidth2D)&
                          NMRspectrumColnames<=(mrbinparam$binRegion[1]-(i-1)*mrbinparam$binwidth2D))
      }
      for(j in 1:mrbinTMP$nbins1){
          for(i in 1:mrbinTMP$nbins2){
             mrbinTMP$binTMP[counter]<<-sum(mrbinTMP$currentSpectrum[which_j[,j],which_i[i,]])/
                                        (sum(which_j[,j])*sum(which_i[i,]))
             counter<-counter+1
         }
      }
   } else {#1D spectra
       if(is.null(mrbinTMP$binNames)) createBinNames()
       mrbinTMP$binTMP<<-rep(0,length(mrbinTMP$binNames))
       names(mrbinTMP$binTMP)<<-mrbinTMP$binNames
       counter<-1
       NMRspectrumNames<-as.numeric(names(mrbinTMP$currentSpectrum))
       for(i in 1:mrbinTMP$nbins){
               indexTMP<-which(NMRspectrumNames>(mrbinparam$binRegion[1]-i*mrbinparam$binwidth1D)&
                            NMRspectrumNames<=(mrbinparam$binRegion[1]-(i-1)*mrbinparam$binwidth1D))
               mrbinTMP$binTMP[counter]<<-sum(mrbinTMP$currentSpectrum[indexTMP])/length(indexTMP)
               counter<-counter+1
       }
    }
    if(sum(mrbinTMP$binTMP==0)>0){
        cat("Warning: Binning region may be larger than total spectrum size.\n")
    }
}

#' A function for reading NMR spectra.
#'
#' This function picks the correct NMR reading function, based on vendor.
#' @param
#' @keywords
#' @export
#' @examples
#' readNMR()

readNMR<-function(){#Read NMR spectral data
      if(mrbinparam$NMRvendor=="Bruker"){
          readBruker()
      }  else {
          stop(paste("No data import function defined for vendor ",mrbinparam$NMRvendor,".\n",sep=""))
      }
}

#' A function for reading Bruker NMR spectra.
#'
#' This function reads Bruker NMR data.
#' @param
#' @keywords
#' @export
#' @examples
#' readBruker()

readBruker<-function(){#Read Bruker NMR spectral data
     datanameDict<-c("1r","2rr")
     names(datanameDict)<-c("1D","2D")
     datanameTmp<-datanameDict[mrbinparam$dimension]
     spectrum_proc_path<-mrbinTMP$currentFolder
     BYTORDP_Dict<-c("little","big")
     names(BYTORDP_Dict)<-c(0,1)
     TITLE<-scan(file=paste(spectrum_proc_path,"/title",sep=""),what="character",sep="\n",quiet=T)[1]
     proc<-scan(file=paste(spectrum_proc_path,"/procs",sep=""),what="character",sep="\n",quiet=T)
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
         proc2<-scan(file=paste(spectrum_proc_path,"/proc2s",sep=""),what="character",sep="\n",quiet=T)
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
        currentSpectrum<-matrix(currentSpectrum,ncol=SI2,byrow =T)#ncol=SI2)
        #Needs to be sorted,XDIM2: number of columns/block, XDIM1: rows/block
        counter<-0
        for(j in 1:(nrow(currentSpectrum)/XDIM1)){
            for(i in 1:(ncol(currentSpectrum)/XDIM2)){
                currentSpectrum[(j-1)*XDIM1+(1:XDIM1),(i-1)*XDIM2+(1:XDIM2)]<-
                   matrix(currentSpectrumTMP[counter*XDIM2*XDIM1 +(1:(XDIM2*XDIM1))],ncol=XDIM2,byrow=T)
                counter<-counter+1
            }
        }
        rownames(currentSpectrum)<-frequencynames1
        colnames(currentSpectrum)<-frequencynames2
     } else {
       names(currentSpectrum)<-frequencynames2
       #if(XDIM2>1){cat("Spectrum might be distorted (XDIM>1)\n")}
     }
  mrbinTMP$currentSpectrum<<-currentSpectrum
  if(mrbinparam$useAsNames=="Spectrum titles")    mrbinTMP$currentSpectrumName<<-TITLE
  if(mrbinparam$useAsNames=="Folder names")    mrbinTMP$currentSpectrumName<<-rev(strsplit(mrbinTMP$currentFolder,"/")[[1]])[4]

  if(mrbinparam$referenceScaling=="Yes") referenceScaling()
}

#' A function for scaling to the reference area.
#'
#' This function scales NMR data to the reference area.
#' @param
#' @keywords
#' @export
#' @examples
#' referenceScaling()

referenceScaling<-function(){
    if(mrbinparam$dimension=="1D"){
         mrbinTMP$currentSpectrum<<-mrbinTMP$currentSpectrum/mean(mrbinTMP$currentSpectrum[
                                         which(as.numeric(names(mrbinTMP$currentSpectrum))<=
                                         mrbinparam$reference1D[1]&
                                         as.numeric(names(mrbinTMP$currentSpectrum))>
                                         mrbinparam$reference1D[2])],na.rm=T)
    }
    if(mrbinparam$dimension=="2D"){
         selectedTMP1<-which(as.numeric(rownames(mrbinTMP$currentSpectrum))<=
                         mrbinparam$reference2D[4]&as.numeric(rownames(
                         mrbinTMP$currentSpectrum))>mrbinparam$reference2D[3])
         selectedTMP2<-which(as.numeric(colnames(mrbinTMP$currentSpectrum))<=
                         mrbinparam$reference2D[1]&
                         as.numeric(colnames(mrbinTMP$currentSpectrum))>
                         mrbinparam$reference2D[2])

         mrbinTMP$currentSpectrum<<-mrbinTMP$currentSpectrum/mean(mrbinTMP$currentSpectrum[
                                     selectedTMP1,selectedTMP2],na.rm=T)
    }
}

#' A function for summing up bins.
#'
#' This function sums up bins. The sums are saved to the middle (median) bin of
#' the original area. All other bins of the area are removed then. Thsi is handy
#' for signals that are know to vary between spectra due to pH or salt content,
#' such as citric acid.
#' @param
#' @keywords
#' @export
#' @examples
#' sumBins()

sumBins<-function(){#sum up regions with shifting peaks and remove remaining bins
  if(nrow(mrbinparam$sumBinList)>0){
    for(i in 1:nrow(mrbinparam$sumBinList)){
       limits<-mrbinparam$sumBinList[i,]
       if(mrbinparam$dimension=="1D"){
          TMP<-as.numeric(colnames(mrbinbins))>limits[2]& as.numeric(colnames(mrbinbins))<limits[1]
          if(sum(TMP)>0){
              i_TMP<-quantile(x=1:sum(TMP), probs = .5,type=3)#define "middle" bin. This one will be kept
              mrbinbins[,which(TMP)[i_TMP]]<<-apply(mrbinbins[,which(TMP)],1,sum)
              mrbinbins<<-mrbinbins[,-which(TMP)[-i_TMP]]
          }
       } else {#2D limits=c(4.04,4.08,58,60)
           NMRdataNames<-t(matrix(as.numeric(unlist(strsplit(colnames(mrbinbins),","))),nrow=2))
           TMP<-NMRdataNames[,2]>limits[2]&NMRdataNames[,2]<limits[1]&
                  NMRdataNames[,1]>limits[3]& NMRdataNames[,1]<limits[4]
           if(sum(TMP)>0){
              i_TMP<-quantile(x=1:sum(TMP), probs = .5,type=3)#define "middle" bin. This one will be kept
              mrbinbins[,which(TMP)[i_TMP]]<<-apply(mrbinbins[,which(TMP)],1,sum)
              mrbinbins<<-mrbinbins[,-which(TMP)[-i_TMP]]
           }
       }
    }
  }
  mrbinparam$numberOfFeaturesAfterSummingBins<<-ncol(mrbinbins)
}

#' A function for removing the solvent region.
#'
#' This function removes the solvent region.
#' @param
#' @keywords
#' @export
#' @examples
#' removeSolvent()

removeSolvent<-function(){
   if(mrbinparam$dimension=="1D"){
       NMRdataNames<-as.numeric(colnames(mrbinbins))
       solventTMP<-NMRdataNames>mrbinparam$solventRegion[2]&NMRdataNames<mrbinparam$solventRegion[1]
   }
   if(mrbinparam$dimension=="2D"){
       NMRdataNames<-t(matrix(as.numeric(unlist(strsplit(colnames(mrbinbins),","))),nrow=2))
       solventTMP<-NMRdataNames[,2]>mrbinparam$solventRegion[2]&NMRdataNames[,2]<mrbinparam$solventRegion[1]
   }
   if(sum(solventTMP)>0)     mrbinbins<<-mrbinbins[,-which(solventTMP)]
   mrbinparam$numberOfFeaturesAfterRemovingSolvent<<-ncol(mrbinbins)
}

#' A function for removing additional regions.
#'
#' This function removes additional regions. This can be useful when some areas
#' are visibly affected by spectral artifacts.
#' @param
#' @keywords
#' @export
#' @examples
#' removeAreas()

removeAreas<-function(){#limits=c(4.75,4.95,-10,160)
  if(nrow(mrbinparam$removeAreaList)>0){
     if(mrbinparam$dimension=="1D")    NMRdataNames<-as.numeric(colnames(mrbinbins))
     if(mrbinparam$dimension=="2D"){
         NMRdataNames<-t(matrix(as.numeric(unlist(strsplit(colnames(mrbinbins),","))),nrow=2))
     }
     removeTMP<-NULL
     for(i in 1:nrow(mrbinparam$removeAreaList)){
         limits<-mrbinparam$removeAreaList[i,]
         if(mrbinparam$dimension=="1D"){
             removeTMP2<-NMRdataNames>limits[2]&NMRdataNames<limits[1]
         }
         if(mrbinparam$dimension=="2D"){
             removeTMP2<-NMRdataNames[,2]>limits[2]&NMRdataNames[,2]<limits[1]&
                NMRdataNames[,1]>limits[3]& NMRdataNames[,1]<limits[4]
         }
         if(sum(removeTMP2)>0)        removeTMP<-c(removeTMP,which(removeTMP2))
     }
     if(!is.null(removeTMP)){
         removeTMP<-unique(removeTMP)
         mrbinbins<<-mrbinbins[,-removeTMP]
     }
  }
  mrbinparam$numberOfFeaturesAfterRemovingAreas<<-ncol(mrbinbins)
}

#' A function for calculating noise levels.
#'
#' This function calculates noise levels.
#' @param
#' @keywords
#' @export
#' @examples
#' calculateNoise()

calculateNoise<-function(){
  if(mrbinparam$dimension=="2D"){
       NMRdataNames<-t(matrix(as.numeric(unlist(strsplit(colnames(mrbinbins),","))),nrow=2))
  }
  noise_level<-rep(0,nrow(mrbinbins))
  names(noise_level)<-rownames(mrbinbins)
  for(j in 1:nrow(mrbinbins)){
    if(mrbinparam$dimension=="1D"){
         noise_level[j]<-sd(mrbinbins[j,which(as.numeric(colnames(mrbinbins))>mrbinparam$noiseRange1d[1]&
                                        as.numeric(colnames(mrbinbins))<mrbinparam$noiseRange1d[2])])
    }
    if(mrbinparam$dimension=="2D"){
         noise_level[j]<-sd(mrbinbins[j,which(NMRdataNames[,2]>mrbinparam$noiseRange2d[1]&
                                NMRdataNames[,2]<mrbinparam$noiseRange2d[2]&
                                NMRdataNames[,1]>mrbinparam$noiseRange2d[3]&
                                NMRdataNames[,1]<mrbinparam$noiseRange2d[4])])
    }
  }
  mrbinparam$noise_level<<-noise_level
}

#' A function for removing bins below noise level.
#'
#' This function checks for each bin (column) whether its level is below the
#' individual noise level times the signal-to-noise ratio. If less than the
#' defined threshold level are above noise*SNR, the whole bin is removed.
#' @param
#' @keywords
#' @export
#' @examples
#' removeNoise()

removeNoise<-function(){#remove noise peaks
    minimumNumber<-max(1,ceiling(mrbinparam$noiseThreshold*nrow(mrbinbins)))
    colnames_NMRdata_no_noise<-NULL
    if(mrbinparam$dimension=="2D"){
         NMRdataNames<-t(matrix(as.numeric(unlist(strsplit(colnames(mrbinbins),","))),nrow=2))
         SNR<-mrbinparam$signal_to_noise2D
    } else {
         SNR<-mrbinparam$signal_to_noise1D
    }
    calculateNoise()
    for(i in 1:ncol(mrbinbins)){#Keep only bins where at least X spectra are > SNR
          if(sum(mrbinbins[,i]>mrbinparam$noise_level*SNR)>=minimumNumber){
              colnames_NMRdata_no_noise<-c(colnames_NMRdata_no_noise,i)#colnames(mrbinbins)[i])
          }
    }
    if(!is.null(colnames_NMRdata_no_noise)){
        mrbinbins<<-mrbinbins[,colnames_NMRdata_no_noise]
    } else {
        stop("No bins above noise level. Noise removal stopped.\n")
    }
  mrbinparam$numberOfFeaturesAfterNoiseRemoval<<-ncol(mrbinbins)
}

#' A function for cropping HSQC spectra.
#'
#' This function crops HSQC spectra to the region along the diagonal to remove
#' uninformative signals. Will work only for 1H-13C HSQC spectra.
#' @param
#' @keywords
#' @export
#' @examples
#' cropNMR()

cropNMR<-function(){
   if(mrbinparam$dimension=="2D"){
     selectedCols<-NULL
     for(j in 1:ncol(mrbinbins)){
        coordTmp<-as.numeric(strsplit(colnames(mrbinbins)[j],",")[[1]])#1=C,2=H
        if(((coordTmp[2]-mrbinparam$croptopLeft[2])*(mrbinparam$cropbottomLeft[1]-mrbinparam$croptopLeft[1])-
            (coordTmp[1]-mrbinparam$croptopLeft[1])*(mrbinparam$cropbottomLeft[2]-mrbinparam$croptopLeft[2]))<0&
           ((coordTmp[2]-mrbinparam$croptopRight[2])*(mrbinparam$cropbottomRight[1]-mrbinparam$croptopRight[1])-
            (coordTmp[1]-mrbinparam$croptopRight[1])*(mrbinparam$cropbottomRight[2]-mrbinparam$croptopRight[2]))>0
          ){#along diagonal? outer product
            selectedCols<-c(selectedCols,j)
        }
    }
    mrbinbins<<-mrbinbins[,selectedCols]
    #plot(t(matrix(as.numeric(unlist(strsplit(colnames(NMRdataTmp),","))),
    #     nrow=2))[,c(2,1)],
    #     xlim=c(10,0),ylim=c(160,0),xlab="1H",ylab="13C",pch=20,col="black")
    #points(t(matrix(as.numeric(unlist(strsplit(colnames(NMRdataTmp2),","))),
    #     nrow=2))[,c(2,1)],col="red",pch=20)
  }
  mrbinparam$numberOfFeaturesAfterCropping<<-ncol(mrbinbins)
}

#' A function for PQN scaling.
#'
#' This function performs PQN scaling. To further exlude unreliable noise, only
#' the most intense signals are used. For 1H-13C HSQC spectra, most of the sugar
#' regions are excluded to avoid a dominating ffect of the multiple sugar
#' signals.
#' @param
#' @keywords
#' @export
#' @examples
#' PQNScaling()

PQNScaling<-function(){#Scale to PQN
 if(ncol(mrbinbins)>1){
  #Create synthetic median spectrum by averaging all spectra
  NMRdataTmp<-rbind(mrbinbins,apply(mrbinbins,2,mean))
  rownames(NMRdataTmp)[nrow(NMRdataTmp)]<-"Median"
  medianSample<-nrow(NMRdataTmp)
  #For HSQC spetra: Remove most sugar signals to get a better fold change estimate
  if(mrbinparam$dimension == "2D" & mrbinparam$cropHSQC=="Yes") {
      selectedCols<-NULL
      for(j in 1:ncol(NMRdataTmp)){
          coordTmp<-as.numeric(strsplit(colnames(NMRdataTmp)[j],",")[[1]])#1=C,2=H
                 if(!(coordTmp[2]>mrbinparam$PQNsugarArea[2]&coordTmp[2]<mrbinparam$PQNsugarArea[1]&
                      coordTmp[1]>mrbinparam$PQNsugarArea[3]&coordTmp[1]<mrbinparam$PQNsugarArea[4])){
                          selectedCols<-c(selectedCols,j)
                 }
      }
      NMRdataTmp2<-NMRdataTmp[,selectedCols]
  } else { #1D spectra
    NMRdataTmp2<-NMRdataTmp
  }
  #Calculate fold changes versus reference sample
  NMRdataTmp_scaledMedian<-NMRdataTmp
  if(mrbinparam$PQNminimumFeatures>ncol(NMRdataTmp2)) mrbinparam$PQNminimumFeatures<-ncol(NMRdataTmp2)
  medianFoldChanges<-rep(0,nrow(NMRdataTmp_scaledMedian))
  names(medianFoldChanges)<-rownames(NMRdataTmp_scaledMedian)
  if(mrbinparam$PQNshowHist){#Plot distribution of fold changes per sample
    par(mfrow=c(ceiling(sqrt(nrow(NMRdataTmp2))),ceiling(sqrt(nrow(NMRdataTmp2)))))
    par(mar=c(.1,.1,3,.1))
  }
  for(i in 1:nrow(NMRdataTmp2)){#scale to spectral area, use only points above X% quantile for better reliability
      overlapAboveNoise<- sort(NMRdataTmp2[i,],index.return=T,decreasing=T)$ix[1:mrbinparam$PQNminimumFeatures]
      medianFoldChanges[i]<-median(NMRdataTmp2[i,overlapAboveNoise]/
                                   NMRdataTmp2[medianSample,overlapAboveNoise])
      NMRdataTmp_scaledMedian[i,]<-NMRdataTmp[i,]/medianFoldChanges[i]
    if(mrbinparam$PQNshowHist){#Plot distribution of fold changes per sample
          hist(NMRdataTmp2[i,overlapAboveNoise]/NMRdataTmp2[medianSample,overlapAboveNoise],breaks=60,
          main=rownames(NMRdataTmp2)[i],xlab="",ylab="")
          lines(rep(medianFoldChanges[i],2),c(0,20000),col='red')
    }
  }
  #Remove reference sample from list
  NMRdataTmp_scaledMedian<-NMRdataTmp_scaledMedian[-nrow(NMRdataTmp_scaledMedian),]
  mrbinbins<<-NMRdataTmp_scaledMedian
  mrbinparam$medians<<-medianFoldChanges
 } else {
    stop("Too few samples to perform PQN normalization.\n")
 }
}

#' A function for plotting the current 2D bin data.
#'
#' This function plots the current binned spectrum. Works only for 2D data.
#' @param
#' @keywords
#' @export
#' @examples
#' plotBins()

plotBins<-function(){#Plot 2D NMR bin data
    colorRampHSQC<-colorRamp(c("blue","green","orange","red","red"),space="rgb")
    NMRdata_select<-mrbinbins[mrbinTMP$currentSpectrumName,]
    NMRdata_select<-NMRdata_select[order(abs(NMRdata_select))]
    allBins2<-t(matrix(as.numeric(unlist(strsplit(names(NMRdata_select),","))),nrow=2))
    #par(bg = 'black', fg='gray23')
    plot(allBins2[,2],allBins2[,1],ylim=c(156,-5),xlim=c(11,-1),
         cex=.5,pch=15
         ,main=mrbinTMP$currentSpectrumName
         ,col=rgb(colorRampHSQC(NMRdata_select/max(abs(NMRdata_select))/2+.5),maxColorValue=255))
    #text(5,120,mrbinTMP$currentSpectrumName,col='gray23'       )
}

#' A function for performing PCA.
#'
#' This function performs PCA, then plots PC1 and PC2 and loading plots.
#' @param
#' @keywords
#' @export
#' @examples
#' PCA()

PCA<-function(){
 if(ncol(mrbinbins)>1){
      if(is.null(mrbinparam$Factors)) mrbinparam$Factors<<-factor(rep("Group 0",nrow(mrbinbins)))
      colorPalette<-rainbow(length(levels(mrbinparam$Factors)))
      mrbinTMP$PCA<<-prcomp(mrbinbins)
      par(mfrow=c(1,2),xpd=T,mar=c(8.1,4.1,4.1,2.1))
      plot(mrbinTMP$PCA$x[,1],mrbinTMP$PCA$x[,2],
           pch=as.numeric(mrbinparam$Factors)+14,
           col=colorPalette[as.numeric(mrbinparam$Factors)],
           main="PCA",
           xlab=paste("PC1 (",round(100*mrbinTMP$PCA$sdev[1]/sum(mrbinTMP$PCA$sdev),1),"%)",sep=""),
           ylab=paste("PC2 (",round(100*mrbinTMP$PCA$sdev[2]/sum(mrbinTMP$PCA$sdev),1),"%)",sep="")
           ,cex=.75
           )
      numlevels<-NULL
      for(i in 1:nlevels(mrbinparam$Factors)) numlevels<-c(numlevels,as.numeric(
                         mrbinparam$Factors[which(mrbinparam$Factors==levels(mrbinparam$Factors)[i])][1]))
      legend("bottomleft", inset=c(-.1,-.25),#c(-0.1*(max(mrbinTMP$PCA$x[,1])-min(mrbinTMP$PCA$x[,1])),-0.9*(max(mrbinTMP$PCA$x[,2])-min(mrbinTMP$PCA$x[,2]))),
              legend=levels(mrbinparam$Factors),
              col=colorPalette[numlevels],
              pch=numlevels+14,
              cex=.75)
      text(mrbinTMP$PCA$x,labels=paste(substr(rownames(mrbinTMP$PCA$x),1,mrbinparam$PCAtitlelength)),pos=3,cex=.5,
             col=colorPalette[as.numeric(mrbinparam$Factors)])
      par(xpd=F)
      plot(mrbinTMP$PCA$rotation,pch=16,cex=.5,main="Loadings Plot")
      text(mrbinTMP$PCA$rotation,labels=substr(rownames(mrbinTMP$PCA$rotation),1,12),pos=4,cex=.5)
 } else {
    stop("Too few samples to perform PQN normalization.\n")
 }
}

#' A function for plotting NMR spectra.
#'
#' This function plots the current NMR spectrum. To change the plot, use zoom(),
#' zoomIn(), zoomOut(), intPlus(), intMin(), left(), right().
#' For 2D data use additionally: contMin(), contPlus(), up(), down()
#' @param
#' @keywords
#' @export
#' @examples
#' plotNMR()

plotNMR<-function(){
   if(is.matrix(mrbinTMP$currentSpectrum)){#2D spectra
      spectrumTMP<-mrbinTMP$currentSpectrum[which(as.numeric(rownames(mrbinTMP$currentSpectrum))<
                           mrbinplot$plotRegion[4]&
                          as.numeric(rownames(mrbinTMP$currentSpectrum))>=mrbinplot$plotRegion[3]),
                         which(as.numeric(colnames(mrbinTMP$currentSpectrum))>=mrbinplot$plotRegion[2]&
                         as.numeric(colnames(mrbinTMP$currentSpectrum))<mrbinplot$plotRegion[1])
                         ]
      if(sum(spectrumTMP<(mrbinplot$lowestContour*max(mrbinTMP$currentSpectrum)))>0){
           spectrumTMP[spectrumTMP<(mrbinplot$lowestContour*max(mrbinTMP$currentSpectrum))]<-0
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
      if(mrbinplot$heatmap){
          heatmap(spectrumTMP, Rowv = NA, Colv = NA,col =rev(heat.colors(128)),
               scale="none",
               labRow="", labCol="")
      } else {
        options(max.contour.segments=1000)
        contour(x = -as.numeric(colnames(spectrumTMP)),
          y = -as.numeric(rownames(spectrumTMP)),
          z = t(spectrumTMP)*mrbinplot$intensityScale,
          levels =  (mrbinplot$lowestContour+
                        .8*(1:mrbinplot$nContours-1)/
                        (mrbinplot$nContours)*(1-
                        mrbinplot$lowestContour))* max(spectrumTMP)
          ,drawlabels =F
          ,col = rev(rainbow(mrbinplot$nContours))
          ,xaxt="none",yaxt="none"
          ,lwd=1,
          main=mrbinTMP$currentSpectrumName
          )
        magnitude2<-10^round(log(max(as.numeric(colnames(spectrumTMP)))-
                    min(as.numeric(colnames(spectrumTMP))),base=10))/10
        axis(1,
             at=(0:100*magnitude2+floor(min(-as.numeric(colnames(spectrumTMP)))/magnitude2)*magnitude2)
             ,labels=-(0:100*magnitude2+floor(min(-as.numeric(colnames(spectrumTMP)))/magnitude2)*magnitude2)
             )
        magnitude1<-10^round(log(max(as.numeric(rownames(spectrumTMP)))-
                     min(as.numeric(rownames(spectrumTMP))),base=10))/10
        axis(2,
             at=-(0:100*magnitude1+floor(min(as.numeric(rownames(spectrumTMP)))/magnitude1)*magnitude1)
             ,labels=(0:100*magnitude1+floor(min(as.numeric(rownames(spectrumTMP)))/magnitude1)*magnitude1)
             )
      }
   } else {  #1D
      plot(as.numeric(names(mrbinTMP$currentSpectrum)),mrbinTMP$currentSpectrum,
           type="l",xlim=c(mrbinplot$plotRegion[1],mrbinplot$plotRegion[2]),
           ylim=c(min(as.numeric(names(mrbinTMP$currentSpectrum))),
                   max(as.numeric(names(mrbinTMP$currentSpectrum))))/mrbinplot$intensityScale,
           xlab="Chemical shift [ppm]",ylab="Intensity")
    }
}

#' A function for changing plotNMR plots.
#'
#' This function increases the intensity of the current NMR spectrum plot.
#' @param
#' @keywords
#' @export
#' @examples
#' intPlus()

intPlus<-function(){#increase plot intensity
   mrbinplot$intensityScale<<-mrbinplot$intensityScale*2
   plotNMR()
}

#' A function for changing plotNMR plots.
#'
#' This function decreases the intensity of the current NMR spectrum plot.
#' @param
#' @keywords
#' @export
#' @examples
#' intMin()

intMin<-function(){#decrease plot intensity
   mrbinplot$intensityScale<<-mrbinplot$intensityScale*0.5
   plotNMR()
}

#' A function for changing plotNMR plots.
#'
#' This function increases the minimum contour level of the current 2D NMR
#' spectrum plot.
#' @param
#' @keywords
#' @export
#' @examples
#' contPlus()

contPlus<-function(){#decrease plot intensity
   mrbinplot$lowestContour<<-mrbinplot$lowestContour*1.5
   plotNMR()
}

#' A function for changing plotNMR plots.
#'
#' This function decreases the minimum contour level of the current 2D NMR
#' spectrum plot.
#' @param
#' @keywords
#' @export
#' @examples
#' contMin()

contMin<-function(){#decrease plot intensity
   mrbinplot$lowestContour<<-mrbinplot$lowestContour*0.75
   plotNMR()
}

#' A function for changing plotNMR plots.
#'
#' This function changes the plot region of the current NMR plot.
#' @param left New left boundary
#' @param right New right boundary
#' @param top New top boundary
#' @param bottom New bottom boundary
#' @keywords
#' @export
#' @examples
#' zoom()

zoom<-function(left=NULL,right=NULL,top=NULL,bottom=NULL){
   if(is.null(left)) stop("Please set left limit\n")
   if(is.null(right)) stop("Please set right limit\n")
   if(is.matrix(mrbinTMP$currentSpectrum)){
      if(is.null(top)) stop("Please set top limit\n")
      if(is.null(bottom)) stop("Please set bottom limit\n")
   }
   mrbinplot$plotRegion[1]<<-left
   mrbinplot$plotRegion[2]<<-right
   mrbinplot$plotRegion[3]<<-top
   mrbinplot$plotRegion[4]<<-bottom
   plotNMR()
}

#' A function for changing plotNMR plots.
#'
#' This function zooms into the plot region of the current NMR plot.
#' @param
#' @keywords
#' @export
#' @examples
#' zoomIn()

zoomIn<-function(){#Zoom into NMR spectrum plot
   if(mrbinparam$dimension=="1D"){
       leftMax<-max(as.numeric(names(mrbinTMP$currentSpectrum)))
       rightMax<-min(as.numeric(names(mrbinTMP$currentSpectrum)))
       topMax<--10
       bottomMax<-160
   }
   if(mrbinparam$dimension=="2D"){
       leftMax<-max(as.numeric(colnames(mrbinTMP$currentSpectrum)))
       rightMax<-min(as.numeric(colnames(mrbinTMP$currentSpectrum)))
       topMax<-min(as.numeric(rownames(mrbinTMP$currentSpectrum)))
       bottomMax<-max(as.numeric(rownames(mrbinTMP$currentSpectrum)))
   }
   mrbinplot$plotRegion[1]<<-min(leftMax,
                             mrbinplot$plotRegion[1]-(mrbinplot$plotRegion[1]-mrbinplot$plotRegion[2])/4)
   mrbinplot$plotRegion[2]<<-max(rightMax,
                             mrbinplot$plotRegion[1]-(mrbinplot$plotRegion[1]-mrbinplot$plotRegion[2])*3/4)
   mrbinplot$plotRegion[3]<<-max(topMax,
                             mrbinplot$plotRegion[4]-(mrbinplot$plotRegion[4]-mrbinplot$plotRegion[3])*3/4)
   mrbinplot$plotRegion[4]<<-min(bottomMax,
                             mrbinplot$plotRegion[4]-(mrbinplot$plotRegion[4]-mrbinplot$plotRegion[3])/4)
   plotNMR()
}

#' A function for changing plotNMR plots.
#'
#' This function zooms out from the plot region of the current NMR plot.
#' @param
#' @keywords
#' @export
#' @examples
#' zoomOut()

zoomOut<-function(){#Zoom out from NMR spectrum plot
   if(mrbinparam$dimension=="1D"){
       leftMax<-max(as.numeric(names(mrbinTMP$currentSpectrum)))
       rightMax<-min(as.numeric(names(mrbinTMP$currentSpectrum)))
       topMax<--10
       bottomMax<-160
   }
   if(mrbinparam$dimension=="2D"){
       leftMax<-max(as.numeric(colnames(mrbinTMP$currentSpectrum)))
       rightMax<-min(as.numeric(colnames(mrbinTMP$currentSpectrum)))
       topMax<-min(as.numeric(rownames(mrbinTMP$currentSpectrum)))
       bottomMax<-max(as.numeric(rownames(mrbinTMP$currentSpectrum)))
   }
   mrbinplot$plotRegion[1]<<-min(leftMax,
                             mrbinplot$plotRegion[1]+(mrbinplot$plotRegion[1]-mrbinplot$plotRegion[2])/2)
   mrbinplot$plotRegion[2]<<-max(rightMax,
                             mrbinplot$plotRegion[1]-(mrbinplot$plotRegion[1]-mrbinplot$plotRegion[2])*3/2)
   mrbinplot$plotRegion[3]<<-max(topMax,
                             mrbinplot$plotRegion[4]-(mrbinplot$plotRegion[4]-mrbinplot$plotRegion[3])*3/2)
   mrbinplot$plotRegion[4]<<-min(bottomMax,
                             mrbinplot$plotRegion[4]+(mrbinplot$plotRegion[4]-mrbinplot$plotRegion[3])/2)
   plotNMR()
}

#' A function for changing plotNMR plots.
#'
#' This function moves left the plot region of the current NMR plot.
#' @param
#' @keywords
#' @export
#' @examples
#' left()

left<-function(){
   if(mrbinparam$dimension=="1D"){
       leftMax<-max(as.numeric(names(mrbinTMP$currentSpectrum)))
       rightMax<-min(as.numeric(names(mrbinTMP$currentSpectrum)))
   }
   if(mrbinparam$dimension=="2D"){
       leftMax<-max(as.numeric(colnames(mrbinTMP$currentSpectrum)))
       rightMax<-min(as.numeric(colnames(mrbinTMP$currentSpectrum)))
   }
   mrbinplot$plotRegion[1]<<-min(leftMax,
                             mrbinplot$plotRegion[1]+(mrbinplot$plotRegion[1]-mrbinplot$plotRegion[2])*.1)
   mrbinplot$plotRegion[2]<<-max(rightMax,
                             mrbinplot$plotRegion[2]+(mrbinplot$plotRegion[1]-mrbinplot$plotRegion[2])*.1)
   plotNMR()
}

#' A function for changing plotNMR plots.
#'
#' This function moves right the plot region of the current NMR plot.
#' @param
#' @keywords
#' @export
#' @examples
#' right()

right<-function(){
   if(mrbinparam$dimension=="1D"){
       leftMax<-max(as.numeric(names(mrbinTMP$currentSpectrum)))
       rightMax<-min(as.numeric(names(mrbinTMP$currentSpectrum)))
   }
   if(mrbinparam$dimension=="2D"){
       leftMax<-max(as.numeric(colnames(mrbinTMP$currentSpectrum)))
       rightMax<-min(as.numeric(colnames(mrbinTMP$currentSpectrum)))
   }
   mrbinplot$plotRegion[1]<<-min(leftMax,
                             mrbinplot$plotRegion[1]-(mrbinplot$plotRegion[1]-mrbinplot$plotRegion[2])*.1)
   mrbinplot$plotRegion[2]<<-max(rightMax,
                             mrbinplot$plotRegion[2]-(mrbinplot$plotRegion[1]-mrbinplot$plotRegion[2])*.1)
   plotNMR()
}

#' A function for changing plotNMR plots.
#'
#' This function moves down the plot region of the current NMR plot (only 2D).
#' @param
#' @keywords
#' @export
#' @examples
#' down()

down<-function(){
   if(mrbinparam$dimension=="2D"){
       topMax<-min(as.numeric(rownames(mrbinTMP$currentSpectrum)))
       bottomMax<-max(as.numeric(rownames(mrbinTMP$currentSpectrum)))
       mrbinplot$plotRegion[3]<<-max(topMax,
                                 mrbinplot$plotRegion[3]-(mrbinplot$plotRegion[3]-mrbinplot$plotRegion[4])*.1)
       mrbinplot$plotRegion[4]<<-min(bottomMax,
                                 mrbinplot$plotRegion[4]-(mrbinplot$plotRegion[3]-mrbinplot$plotRegion[4])*.1)
       plotNMR()
   }
}

#' A function for changing plotNMR plots.
#'
#' This function moves up the plot region of the current NMR plot (only 2D).
#' @param
#' @keywords
#' @export
#' @examples
#' up()

up<-function(){
   if(mrbinparam$dimension=="2D"){
       topMax<-min(as.numeric(rownames(mrbinTMP$currentSpectrum)))
       bottomMax<-max(as.numeric(rownames(mrbinTMP$currentSpectrum)))
       mrbinplot$plotRegion[3]<<-max(topMax,
                                 mrbinplot$plotRegion[3]+(mrbinplot$plotRegion[3]-mrbinplot$plotRegion[4])*.1)
       mrbinplot$plotRegion[4]<<-min(bottomMax,
                                 mrbinplot$plotRegion[4]+(mrbinplot$plotRegion[3]-mrbinplot$plotRegion[4])*.1)
       plotNMR()
   }
}
