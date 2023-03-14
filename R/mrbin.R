#mrbin - Collection of R functions for processing and analyzing metabolomics data.
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
#Description: This package is a collection of functions for processing and
#    analyzing metabolomics data. The mrbin function converts 1D or 2D nuclear
#    magnetic resonance data into a matrix of values
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


#' A function executed when attaching this package
#'
#' This function displays a welcome message.
#' @param  libname Library name
#' @param  pkgname Package name
#' @return {None}
#' @keywords internal
#' @noRd
#' @examples
#' \donttest{ .onAttach() }

.onAttach <- function(libname, pkgname){
    packageStartupMessage("mrbin 1.7.3\nFor instructions and examples, please type: vignette('mrbin')")
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


#' A function returning predicted values for use with the fia function.
#'
#' This function predicts group membership and returns a numeric vector with results.
#' @param model A predictive model. Make sure to have loaded all required packages before starting this function
#' @param dataSet A matrix or dataframe containing data, depending on what your predict function requires. Columns=features, rows=samples
#' @param functionNamePredict The name of the prediction function. This only needs to be changed if the prediction function is not called predict
#' @param parameterNameObject The name of the parameter for passing the model to the prediction function
#' @param parameterNameData The name of the parameter for passing the data to the prediction function
#' @param firstLevel Numeric value of first level or group. Usually 1 but for glm such as in the example this needs to be 0.
#' @param dataFrameFlag Logical value indicating whether the data object is a data frame rather than a matrix.
#' @param ... Optional, additional parameters that will be passed to the prediction function.
#' @return A numeric (integer) vector of predicted group memberships.
#' @export
#' @examples
#'  #First, define group membership and create the example feature data
#'  group<-factor(c(rep("Group1",4),rep("Group2",5)))
#'  names(group)<-paste("Sample",1:9,sep="")
#'  dataset<-data.frame(
#'    Feature1=c(5.1,5.0,6.0,2.9,4.8,4.6,4.9,3.8,5.1),
#'    Feature2=c(2.6,4.0,3.2,1.2,3.1,2.1,4.5,6.1,1.3),
#'    Feature3=c(3.1,6.1,5.8,5.1,3.8,6.1,3.4,4.0,4.4),
#'    Feature4=c(5.3,5.2,3.1,2.7,3.2,2.8,5.9,5.8,3.1),
#'    Feature5=c(3.2,4.4,4.8,4.9,6.0,3.6,6.1,3.9,3.5)
#'    )
#'  rownames(dataset)<-names(group)
#'  #train a model - here we use a logit model instead of ANN as a demonstration
#'  mod<-glm(group~Feature1+Feature2+Feature3+Feature4+Feature5,
#'    data=data.frame(group=group,dataset),family="binomial")
#'  predictWrapper(model=mod,dataSet=dataset,firstLevel=0,type="response")
predictWrapper<-function(model,dataSet,functionNamePredict="predict",firstLevel=1,
  parameterNameObject="object",parameterNameData="x",dataFrameFlag=FALSE,...){
  if(dataFrameFlag) dataSet<-data.frame(dataSet)
  predParam<-c(list(model,dataSet),...)
  names(predParam)[1]<-parameterNameObject
  names(predParam)[2]<-parameterNameData
  predictionsTMP<-do.call(functionNamePredict,predParam)
  if(is.matrix(predictionsTMP)){#tensorflow ANN output is matrix
     predictionsFinal<-apply(predictionsTMP,1,which.max)
  } else {
     predictionsFinal<-round(predictionsTMP)
  }
  predictionsFinal<-predictionsFinal+1-firstLevel
  return(predictionsFinal)
}

#' A function identifying features of importance.
#'
#' This function finds features that can change the outcomes of a model's prediction.
#' Example: fia=1.00 means single compound found in all but 0 percent of samples.
#' fia=2.45 indicates this compound is found in pairs in all but 45 percent of tested samples
#' A function named predict needs to be present for this to work. If the function name
#' of the prediction function is different, the function name has to be provided in
#' the parameter functionNamePredict.
#' @param model A predictive model. Make sure to have loaded all required packages before starting this function
#' @param dataSet An object containing data, columns=features, rows=samples. This should be either a matrix or a dataframe, depending on which of these two the specific prediction function requires
#' @param factors A factor vector with group membership of each sample in the data set. Order of levels must correspond to the number predicted by the model
#' @param nSeed Number of times that the test will be repeated, selecting different random features
#' @param numberOfSamples Number of samples that will be randomly chosen from each group
#' @param maxFeatures Maximum number of features that will be tested. Larger numbers will be split into child nodes without testing to increase speed
#' @param innerLoop Number of repeated loops to test additional child nodes
#' @param verbose A logical vector to turn messages on or off
#' @param maxNumberAllTests Combinations of features of this length or shorter will not be split in half to create two children, but into multiple children with one feature left out each. This is done make sure no combination is missed.
#' @param firstLevel Numeric value of first level or group. Usually 1 but for glm such as in the example this needs to be 0.
#' @param saveMemory Save memory by performing only two predictions per step, which will be much slower. If you are using keras, use parameter kerasClearMemory=2 instead to free more memory and be a lot faster. FALSE to turn off.
#' @param kerasClearMemory Save memory by clearing model from memory and reloading the model between chunks of predictions. Will only work when using package keras. 0=off, 1=medium (reload between repeat with different seeds), 2=maximum memory savings (reload after each run for a single sample). This will write a model file to the working directory.
#' @param functionNamePredict The name of the prediction function. This only needs to be changed if the prediction function is not called predict
#' @param parameterNameObject The name of the parameter for passing the model to the prediction function
#' @param parameterNameData The name of the parameter for passing the data to the prediction function
#' @param ... Optional, additional parameters that will be passed to the prediction function.
#' @return A list of results: scoresSummary A vector of fia scores for the whole dataset; scores contains vectors of fia scores for each predicted group; scoresIndividual A list of fia scores for each individual sample; fiaListPerSample A list of important combinations of features for each predicted sample; fiaMatrix A list of fia scores for each predicted group.
#' @export
#' @examples
#'  #First, define group membership and create the example feature data
#'  group<-factor(c(rep("Group1",4),rep("Group2",5)))
#'  names(group)<-paste("Sample",1:9,sep="")
#'  dataset<-data.frame(
#'    Feature1=c(5.1,5.0,6.0,2.9,4.8,4.6,4.9,3.8,5.1),
#'    Feature2=c(2.6,4.0,3.2,1.2,3.1,2.1,4.5,6.1,1.3),
#'    Feature3=c(3.1,6.1,5.8,5.1,3.8,6.1,3.4,4.0,4.4),
#'    Feature4=c(5.3,5.2,3.1,2.7,3.2,2.8,5.9,5.8,3.1),
#'    Feature5=c(3.2,4.4,4.8,4.9,6.0,3.6,6.1,3.9,3.5),
#'    Feature6=c(6.8,6.7,7.2,7.0,7.3,7.1,7.2,6.9,6.8)
#'    )
#'  rownames(dataset)<-names(group)
#'  #train a model - here we use a logit model instead of ANN as a demonstration
#'  mod<-glm(group~Feature1+Feature2+Feature3+Feature4+Feature5+Feature6,
#'    data=data.frame(group=group,dataset),family="binomial")
#'  fiaresults<-fia(model=mod,dataSet=dataset,factors=group,parameterNameData="newdata",
#'    firstLevel=0,type="response")
#'  fiaresults$scores
fia<-function(model,dataSet,factors,nSeed=6,numberOfSamples=100,
  maxFeatures=10000,innerLoop=300,verbose=TRUE,maxNumberAllTests=5,firstLevel=1,
  saveMemory=FALSE,kerasClearMemory=0,functionNamePredict="predict",
  parameterNameObject="object",parameterNameData="x",...){
  digitDict<-c(1:9,0,letters,LETTERS)#replace number>9 with single digit characters
  seedList=(1:nSeed-1)*100
  if(is.factor(factors)){
    #factors<-droplevels(factors)#this could change the order compared to the ANN
  } else {
    factors<-factor(factors)
  }
  dataFrameFlag<-FALSE
  if(is.data.frame(dataSet)){
    dataSet<-as.matrix(dataSet)
    dataFrameFlag<-TRUE
  }
  if(kerasClearMemory>0){
    keras::save_model_tf(object=model,filepath="fiaTMP",overwrite = TRUE,
      include_optimizer = TRUE,signatures = NULL, options = NULL)
  }
  factorsDict<-1:nlevels(factors)-1#-1 is necessary for tensorflow, otherwise usually not
  names(factorsDict)<-levels(factors)
  lVector<-apply(dataSet,2,quantile,.01)
  hVector<-apply(dataSet,2,quantile,.99)
  predTMP<-predictWrapper(model=model,dataSet=dataSet,firstLevel=firstLevel,
    functionNamePredict=functionNamePredict,parameterNameObject=parameterNameObject,
    parameterNameData=parameterNameData,dataFrameFlag=dataFrameFlag,verbose=0
    ,...
    )
  predAll<-names(factorsDict)[predTMP]
  #Create a sample list
  sampleList<-matrix(NA,nrow=nlevels(factors),ncol=min(c(numberOfSamples,
    nrow(dataSet))))
  sampleList2<-sampleList
  for(j in 1:nlevels(factors)){#pick a subset of all samples for testing
     set.seed(1)     #select samples that are predicted correctly
     samplesTMP<-factors==levels(factors)[j]&factors==predAll
     if(sum(samplesTMP)==0) message("No samples were present or correctly predicted for ",levels(factors)[j])
     repSampleList<-sample(which(samplesTMP),min(c(numberOfSamples,sum(samplesTMP))))
     sampleList[j,1:length(repSampleList)]<-repSampleList
     sampleList2[j,1:length(repSampleList)]<- rep(levels(factors)[j],length(repSampleList))
  }
  sampleVector<-as.vector(sampleList2)
  names(sampleVector)<-as.vector(sampleList)
  sampleVector<-sampleVector[!is.na(sampleVector)]
  samplesPositive<-vector("list",length(sampleVector))
  names(samplesPositive) <- names(sampleVector)
  for(i in 1:length(samplesPositive)){
    samplesPositive[[i]]<-list()
  }
  #for each sample, save single important features
  if(verbose){
    message("Testing single features\n0% ",appendLF = FALSE)
    flush.console()
  }
  stepSizePercent<-20
  steps<-1
  #save memory by doing this for each starting seed:
  #create matrix prefilled with Nfeatures (first digit of FIA). save here length 
  #of each saved pair, if lower than previous value
  fiaMatrixNew<-matrix(ncol(dataSet),ncol=ncol(dataSet),nrow=length(sampleVector))#nrow=nlevels(factors))
  colnames(fiaMatrixNew)<-colnames(dataSet)
  rownames(fiaMatrixNew)<-names(sampleVector)#levels(factors)
  
  irepSample<-1
  for (irepSample in 1:length(sampleVector)){ #loop over all selected samples
    repSample<-names(sampleVector)[irepSample]
    dataTMP<-dataSet[rep(as.numeric(repSample),2*ncol(dataSet)),,drop=FALSE]
    dataTMP[cbind(1:ncol(dataSet),1:ncol(dataSet))]<-
      lVector[1:ncol(dataSet)] #replace value by low value
    dataTMP[cbind(ncol(dataSet)+1:ncol(dataSet),1:ncol(dataSet))]<-
      hVector[1:ncol(dataSet)] #replace value by high value
    predTMP<-predictWrapper(model=model,dataSet=dataTMP,firstLevel=firstLevel,
      functionNamePredict=functionNamePredict,parameterNameObject=parameterNameObject,
      parameterNameData=parameterNameData,dataFrameFlag=dataFrameFlag
      ,verbose=0,...
      )
    pred<-names(factorsDict)[predTMP]
    iSingle<-1
    for(iSingle in 1:ncol(dataSet)){
      if(sum(!pred[c(iSingle,iSingle+ncol(dataSet))]==
        sampleVector[repSample])>0){#check low and high levels
          samplesPositive[[repSample]]<-c(
             samplesPositive[[repSample]],colnames(dataSet)[iSingle])
		  fiaMatrixNew[repSample,iSingle]<-1
      }
    }
    if(irepSample/length(sampleVector)>=steps*stepSizePercent/100){
      if(verbose){
        message(steps*stepSizePercent,"% ",appendLF = FALSE)
        flush.console()
        if(irepSample==length(sampleVector) ) message("\n",appendLF = FALSE)
      }
      steps<-steps+1
    }
  }
  if(verbose){
     message("Testing combinations of features\n0% ",appendLF = FALSE)
     flush.console()
  }
  stepSizePercent<-20
  steps<-1
  steps2<-1
  for(iSeed2 in 1:length(seedList)){ #repeat for different starting seeds
    iSeed<-seedList[iSeed2]
    fiaListL<-list()
    irepSample<-1#debug
    for (irepSample in 1:length(sampleVector)){ #loop over all selected samples
      repSample<-names(sampleVector)[irepSample]
      i2<-1
      i3<-"1."#. means do not test
      fiaResults <- vector("list", ceiling(log2(ncol(dataSet)))+10)
      for(iTMP in 1:length(fiaResults)){
        fiaResults[[iTMP]] <- list()
      }
      #remove positive single hits
      positiveSingleTMP<-setdiff(colnames(dataSet),
          samplesPositive[[repSample]])
      if(length(positiveSingleTMP)>0) fiaResults[[i2]][[i3]]<-positiveSingleTMP
      for(iInnerLoop in 1:innerLoop){
        seedInnerLoop<-(iInnerLoop-1)*1000
        if(length(fiaResults)>0){
          for(i2 in 1:length(fiaResults)){
           if(length(fiaResults[[i2]])>0){
             i3TMP<-names(fiaResults[[i2]])#to avoid skipping after deleting entries
             if(!saveMemory){
               if(i2==1|length(fiaResults[[i2]][[1]])>=maxFeatures){
                 pred2=rep("",2*length(fiaResults[[i2]]))
               } else {
                 dataTMP<-dataSet[rep(as.numeric(repSample),2*length(fiaResults[[i2]])),,drop=FALSE]
                 for(iPredTMP in 1:length(fiaResults[[i2]])){
                   dataTMP[iPredTMP,fiaResults[[i2]][[iPredTMP]]]<-
                     lVector[fiaResults[[i2]][[iPredTMP]]] #replace value by low value
                   dataTMP[length(fiaResults[[i2]])+iPredTMP,fiaResults[[i2]][[iPredTMP]]]<-
                     hVector[fiaResults[[i2]][[iPredTMP]]] #replace value by high value
                 }
                 predTMP<-predictWrapper(model=model,dataSet=dataTMP,firstLevel=firstLevel,
                   functionNamePredict=functionNamePredict,
                   parameterNameObject=parameterNameObject,
                   parameterNameData=parameterNameData,dataFrameFlag=dataFrameFlag,verbose=0
				   ,...
				   )
                 pred2<-names(factorsDict)[predTMP]
				 rm(predTMP)
               }
             }
             #i3b<-1
             for(i3b in 1:length(i3TMP)){
              i3<-i3TMP[i3b]
              if(saveMemory){
                if(i2==1|length(fiaResults[[i2]][[i3]])>=maxFeatures){
                  pred<-"" #first step(s) assumed to be positive
                } else {
                   dataTMP<-dataSet[c(as.numeric(repSample),as.numeric(repSample)),,drop=FALSE]
                   dataTMP[1,fiaResults[[i2]][[i3]]]<-lVector[fiaResults[[i2]][[i3]]] #replace value i by low value
                   dataTMP[2,fiaResults[[i2]][[i3]]]<-hVector[fiaResults[[i2]][[i3]]] #replace value i by high value
                   predTMP<-predictWrapper(model=model,dataSet=dataTMP,firstLevel=firstLevel,
                     functionNamePredict=functionNamePredict,
                     parameterNameObject=parameterNameObject,
                     parameterNameData=parameterNameData,dataFrameFlag=dataFrameFlag
					 ,verbose=0,...
					 )
                   pred<-names(factorsDict)[predTMP]
				   rm(predTMP)
                }
              } else {
                pred<-pred2[c(i3b,length(i3TMP)+i3b)]
              }
              parentNameTMP<-substr(i3,1,nchar(i3)-2)
              if(sum(!pred==sampleVector[repSample])==0){#correct prediction
                 #if all children were checked and the parent is still there,
                 #this means the parent is positive while all children are negative.
                 #in this case the parent will be saved to the list and then deleted
                 if(which(digitDict==substr(i3,nchar(i3)-1,nchar(i3)-1))==
                   (length(fiaResults[[i2]][[i3]])+1)){#this indicates that all children have been tested
                   if(i2>1){#save parent
                     if(!is.null(fiaResults[[i2-1]][[parentNameTMP]])){
					   for(ifiaMatrixNew in 1:length(fiaResults[[i2-1]][[parentNameTMP]])){#check all features
					     if(length(fiaResults[[i2-1]][[parentNameTMP]])<fiaMatrixNew[repSample,fiaResults[[i2-1]][[parentNameTMP]][ifiaMatrixNew]]){
					       fiaMatrixNew[repSample,fiaResults[[i2-1]][[parentNameTMP]][[ifiaMatrixNew]]]<-length(fiaResults[[i2-1]][[parentNameTMP]])
					     }
					   }
                       fiaResults[[i2-1]][[parentNameTMP]]<-NULL
                     }
                   }
                 }
                 fiaResults[[i2]][[i3]]<-NULL
              } else {#incorrect prediction, i.e. positive result
                if(i2>1){#delete parent
                  fiaResults[[i2-1]][[parentNameTMP]]<-NULL
                }
                #split into two children by selecting half of features
                if(length(fiaResults[[i2]][[i3]])>maxNumberAllTests){
                   set.seed(i2+iSeed+seedInnerLoop)
                   fiaResults[[i2+1]][[paste(i3,"1-",sep="",collapse="")]]<-
                     sort(sample(fiaResults[[i2]][[i3]],ceiling(length(
                     fiaResults[[i2]][[i3]])/2)),decreasing=TRUE)
                   fiaResults[[i2+1]][[paste(i3,"2-",sep="",collapse="")]]<-
                     (setdiff(fiaResults[[i2]][[i3]],fiaResults[[i2+1]][[
                     paste(i3,"1-",sep="",collapse="")]]))
                } else {#split low numbers manually by leaving one out each
                   iLowNumbers<-length(fiaResults[[i2]][[i3]])
                   for(iLowNumbers2 in 1:iLowNumbers){
                      fiaResults[[i2+1]][[paste(i3,digitDict[iLowNumbers2]
                        ,"-",sep="",collapse="")]]<-
                        fiaResults[[i2]][[i3]][setdiff(
                        1:iLowNumbers,iLowNumbers2)]
                   }
                }
              }
             }
           }
          }
        }
     }
     samplesTestedAdditionalTMP<-unlist(
       fiaResults[which(lapply(fiaResults,length)>0)],recursive=FALSE)
     if(length(samplesTestedAdditionalTMP)>0){
       for(ifiaList in 1:length(samplesTestedAdditionalTMP)){
		 for(ifiaMatrixNew in 1:length(samplesTestedAdditionalTMP[[ifiaList]])){#check all features
			 if(length(samplesTestedAdditionalTMP[ifiaList])<fiaMatrixNew[repSample,samplesTestedAdditionalTMP[[ifiaList]][[ifiaMatrixNew]]]){
			   fiaMatrixNew[repSample,samplesTestedAdditionalTMP[[ifiaList]][[ifiaMatrixNew]]]<-length(samplesTestedAdditionalTMP[[ifiaList]])
			 }
		 }
       }
	   rm(samplesTestedAdditionalTMP)
     }
	 
     if(steps2/(length(seedList)*length(sampleVector))>=steps*stepSizePercent/100){
       if(verbose){
         message(steps*stepSizePercent,"% ",appendLF = FALSE)
         flush.console()
         if(steps2>=(length(seedList)*length(sampleVector))){
           if(steps*stepSizePercent<100) message("100% ",appendLF = FALSE)
           message("\n",appendLF=FALSE)
           flush.console()
         } 
       }
       steps<-steps+1
     }
     steps2<-steps2+1
	 gc()
     if(kerasClearMemory==2){
	   keras::k_clear_session()
	   keras::load_model_tf(filepath="fiaTMP",custom_objects=NULL, compile=TRUE)
     }
   }
   if(kerasClearMemory>0){
	   keras::k_clear_session()
	   keras::load_model_tf(filepath="fiaTMP",custom_objects=NULL, compile=TRUE)
   }
  }
  #calculate fia scores for the whole data set
  fiaMatrix<-fiaMatrixNew
  fiaMatrixTMP<-fiaMatrix
  fiaAllTMP<-apply(fiaMatrixTMP,2,min,na.rm=TRUE)
  fiaAll<-sort(fiaAllTMP+(1-apply(fiaMatrixTMP==matrix(rep(fiaAllTMP,nrow(fiaMatrixTMP)),
     nrow=nrow(fiaMatrixTMP),byrow=TRUE),2,sum,na.rm=TRUE)/nrow(fiaMatrixTMP)))
  #fia scores per group
  scores<-vector("list",nlevels(factors))
  names(scores)<-levels(factors)
  for(iScores in 1:length(scores)){
    if(sum(sampleVector[rownames(fiaMatrix)]==names(scores)[iScores])>0){
      fiaMatrixTMP<-fiaMatrix[sampleVector[rownames(fiaMatrix)]==names(scores)[
        iScores],,drop=FALSE]
      fiaAllTMP<-apply(fiaMatrixTMP,2,min,na.rm=TRUE)
      scores[[iScores]]<-sort(fiaAllTMP+(1-apply(fiaMatrixTMP==matrix(rep(
         fiaAllTMP,nrow(fiaMatrixTMP)),
         nrow=nrow(fiaMatrixTMP),byrow=TRUE),2,sum,na.rm=TRUE)/nrow(fiaMatrixTMP)))
    }
  }
  #fia scores per sample
  scoresIndividual<-vector("list",nrow(fiaMatrix))
  names(scoresIndividual)<-rownames(fiaMatrix)
  for(iScores in 1:length(scoresIndividual)){
      fiaMatrixTMP<-fiaMatrix[iScores,,drop=FALSE]
      fiaAllTMP<-apply(fiaMatrixTMP,2,min,na.rm=TRUE)
      scoresIndividual[[iScores]]<-sort(fiaAllTMP+(1-apply(fiaMatrixTMP==matrix(rep(
         fiaAllTMP,nrow(fiaMatrixTMP)),
         nrow=nrow(fiaMatrixTMP),byrow=TRUE),2,sum,na.rm=TRUE)/nrow(fiaMatrixTMP)))
  }
  return(list(
    scores=scores,
    scoresSummary=fiaAll,
    scoresIndividual=scoresIndividual,
    fiaMatrix=fiaMatrix))
}


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
#' @param NMRdata A matrix or mrbin object containing NMR data. For matrix: columns=frequencies,rows=samples
#' @param noiseLevels A vector (can be omitted if NMRdata is an mrbin object)
#' @param verbose Should a summary be displayed if NMRdata is an mrbin object
#' @param errorsAsWarnings If TRUE, errors will be turned into warnings. Should be used with care, as errors indicate undocumented changes to the data.
#' @return An invisible matrix or mrbin object containing NMR data without negative values.
#' @export
#' @examples
#'  resetEnv()
#'  Example<-mrbin(silent=TRUE,
#'                    parameters=list(verbose=TRUE,dimension="1D",PQNScaling="No",
#'                    binwidth1D=0.005,signal_to_noise1D=1,PCA="No",binRegion=c(9.5,7.5,10,156),
#'                    saveFiles="No",referenceScaling="No",noiseRemoval="No",
#'                    fixNegatives="No",logTrafo="No",noiseThreshold=.05,
#'                    NMRfolders=c(system.file("extdata/2/10/pdata/10",package="mrbin"),
#'                               system.file("extdata/3/10/pdata/10",package="mrbin"))
#'                    ))
#'  sum(Example$bins<=0)
#'  exampleNMRpositive<-atnv(NMRdata=Example$bins, noiseLevels=Example$parameters$noise_level_adjusted)
#'  sum(exampleNMRpositive<=0)

atnv<-function(NMRdata,noiseLevels=NULL,verbose=TRUE,errorsAsWarnings=FALSE){
     NMRdata2<-NMRdata
     if(methods::is(NMRdata,"mrbin")){
        transformations="Atnv transformed"
        if(transformations %in% NMRdata$transformations){
          if(!errorsAsWarnings) stop("Data has been atnv transformed previously, this could corrupt the data.")
          warning("Data has been atnv transformed previously, this could corrupt the data.")
        }
        if("Log transformed" %in% NMRdata$transformations){
          if(!errorsAsWarnings) stop("Data has been log transformed previously, this would corrupt the data.")
          warning("Data has been log transformed previously, this would corrupt the data.")
        }
         NMRdataTMP <- NMRdata$bins
         noiseLevels <- apply(NMRdata$parameters$noise_level,1,median)
     } else {
       NMRdataTMP<-NMRdata
     }
     if(is.null(noiseLevels)){
         noiseTMP<-sort(NMRdataTMP[NMRdataTMP>0])[ceiling(.01*length(NMRdataTMP[
                        NMRdataTMP>0]))]
     }
     for(i in 1:ncol(NMRdataTMP)){
        negatives<-NMRdataTMP[,i]<=0
        if(sum(negatives)>0){
          if(!is.null(noiseLevels)){#If noise levels are available, restrict range to below noise
               noiseTMP<-stats::median(noiseLevels[negatives])
          }
            minTMP<-min(NMRdataTMP[negatives,i])#select lowest bin
            if(sum(!negatives)>0){
              maxTMP<-min(noiseTMP,min(NMRdataTMP[!negatives,i]))#select lowest bin above 0
            } else {
              maxTMP<-noiseTMP
            }
            NMRdataTMP[negatives,i]<-(NMRdataTMP[negatives,i]+(maxTMP-minTMP))/
                                          (maxTMP-minTMP)*(maxTMP*.99)+
                                           maxTMP*.01
        }
     }
     if(methods::is(NMRdata,"mrbin")){
          if(nrow(NMRdata$bins)==1){
            NMRdata$bins<-matrix(NMRdataTMP,nrow=1)
            rownames(NMRdata$bins)<-rownames(NMRdataTMP)
            colnames(NMRdata$bins)<-colnames(NMRdataTMP)
          } else {
            NMRdata$bins<-NMRdataTMP
          }
     } else {
        NMRdata<-NMRdataTMP
     }
     if(methods::is(NMRdata,"mrbin")){
       NMRdata2<-editmrbin(mrbinObject=NMRdata2,functionName="mrbin::atnv",
           versionNumber=as.character(utils::packageVersion("mrbin")),
           bins=NMRdata$bins, parameters=NMRdata$parameters,transformations=transformations,
           verbose=verbose)
       NMRdata<-NMRdata2
     }
     invisible(NMRdata)
 #}
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
    #assign("bins",NULL,mrbin.env)
    #assign("paramChangeFlag",FALSE,mrbin.env)
    assign("mrbinTMP",list(
               #mrbinversion=as.character(utils::packageVersion("mrbin")),
               #binsRaw=NULL,
               medians=NULL,
               noise_level_TMP=NULL,
               noise_level_Raw_TMP=NULL,
               meanNumberOfPointsPerBin=NULL,
               meanNumberOfPointsPerBin_TMP=NULL,
               #PCA=NULL,
               binTMP=NULL,
               binNames=NULL,
               nbins1=NULL,
               nbins2=NULL,
               nbins=NULL,
               currentFolder=NULL,
               currentSpectrum=NULL,
               currentSpectrumOriginal=NULL,
               currentSpectrumName=NULL,
               currentSpectrumFolderName=NULL,
               currentSpectrumEXPNO=NULL,
               currentSpectrumFolderName_EXPNO=NULL,
               currentSpectrumTitle=NULL,
               i=1,
               additionalPlots1D=NULL,
               additionalPlots1DMetadata=NULL,
               additionalPlots2D=NULL,
               additionalPlots2DMetadata=NULL,
               spectrumListPlotTMP=NULL,
               timeEstimate=0,
               scaleFactorTMP1=NULL,
               scaleFactorTMP2=NULL,
               scaleFactorTMP3=NULL
    ),mrbin.env)
    assign("requiredParam",c(
               "dimension","binMethod","specialBinList","binRegion","referenceScaling",
               "removeSolvent","solventRegion","removeAreas","removeAreaList",
               "sumBins","sumBinList",
               "noiseRemoval","noiseThreshold","dilutionCorrection","PQNScaling",
               "fixNegatives","logTrafo","unitVarianceScaling","PQNminimumFeatures",
               "PQNIgnoreSugarArea","PQNsugarArea","saveFiles","useAsNames","outputFileName",
               "PCAtitlelength","PCA","tryParallel","NMRfolders"
               ),mrbin.env)
    assign("requiredParam1D",c(
               "binwidth1D","reference1D","signal_to_noise1D","noiseRange1d",
               mrbin.env$requiredParam
               ),mrbin.env)
    assign("requiredParam2D",c(
               "binwidth2D","binheight","reference2D","signal_to_noise2D","cropHSQC",
               "noiseRange2d","croptopRight","croptopLeft","cropbottomRight","cropbottomLeft",
               mrbin.env$requiredParam
               ),mrbin.env)
    assign("requiredMetadata",c(
               "projectTitle","projectIdentifier","projectDescription","projectAuthors",
               "projectApprovals","sampleType","organism","solvent","pH","dilutionFactors",
               "factors","metaData","metaboliteIdentities"
               ),mrbin.env)
    assign("mrbin", createmrbin(),mrbin.env)
    assign("parameters_copy",mrbin.env$mrbin$parameters,mrbin.env)
    assign("mrbinplot",list(
               lowestContour=.01,
               plotRegion=NULL,
               intensityScale=1,
               intensityScale2D=1,
               intensityOffset=0,
               nContours=30,#60
               heatmap=FALSE),mrbin.env)
}

#' A function setting the parameters and performing binning and data processing
#'
#' This function guides the user through the set-up of parameters, starts binning
#' and performs the chosen data processing steps.
#' If a list of parameters is provided and silent is set to TRUE, no user input
#' is requested and binning and data processing are performed silently.
#' @param parameters Optional: A list of parameters, see examples for details. If omitted, the user will be asked through a series of question to set the parameters.
#' @param metadata Optional: A list of metadata. If omitted, the user can add metadata after generating bin data.
#' @param silent If TRUE, the user will be asked no questions and binning and data analysis will run according to the current parameters. Defaults to FALSE.
#' @param setDefault If TRUE, all current parameters will be replaced by the default parameters (before loading any provided parameters sets). Defaults to FALSE.
#' @return An invisible object of type "mrbin" containing bins (data after processing), parameters, and factors
#' @export
#' @examples
#' # Set parameters in command line.
#' mrbinResults<-mrbin(silent=TRUE,setDefault=TRUE,parameters=list(
#'                 dimension="1D",binwidth1D=0.01,tryParallel=FALSE,
#'                 signal_to_noise1D=25,noiseThreshold=0.75,
#'                 NMRfolders=c(
#'                 system.file("extdata/1/10/pdata/10",package="mrbin"),
#'                 system.file("extdata/2/10/pdata/10",package="mrbin"),
#'                 system.file("extdata/3/10/pdata/10",package="mrbin"))
#'                 ))

mrbin<-function(silent=FALSE,setDefault=FALSE,parameters=NULL,metadata=NULL){
 if(!exists("mrbin.env", mode="environment")) .onLoad()
 if(setDefault) resetEnv()
 if(!is.null(parameters)) setParam(parameters=parameters)
 if(!is.null(metadata)) setParam(metadata=metadata)
 if(!is.null(mrbin.env$mrbin$parameters$Factors)){
   mrbin.env$mrbin$metadata$factors<-mrbin.env$mrbin$parameters$Factors  #for backward compatibility
   mrbin.env$mrbin$parameters$Factors<-NULL
 }
 stopTMP<-FALSE
 selectionRepeat<-""
 #if(silent){
   startmrbin<-"Start binning now"
 #}
 restart<-TRUE
 while(restart){#for restarting during data review
  restart<-FALSE
  mrbin.env$mrbin$parameters$warningMessages<-NULL
  if(!silent){
   if(Sys.info()['sysname']=='Darwin'){#On Apple or Mac computer, display a hint for installing Quartz
     if(mrbin.env$mrbin$parameters$verbose){
       message("Hint: If you see text lists instead of dialog boxes, please install xquartz \nfrom https://www.xquartz.org")
       utils::flush.console()
     }
   }
   selectStep<--3
   lastStepDone<-FALSE
   while(!lastStepDone&!stopTMP){
     if(selectStep==-3){#Show hints
        if(!mrbin.env$mrbin$parameters$verbose){
           selectTMPNo<-"Do not show verbose hints and results"
           selectTMPYes<-"Show verbose hints and results (recommended)"
           selection<-utils::select.list(c(selectTMPYes,selectTMPNo),preselect=selectTMPYes,
                      title="Show verbose hints?",graphics=TRUE)
           if(length(selection)==0|selection=="") stopTMP<-TRUE
           if(!stopTMP){
            if(selection==selectTMPYes) mrbin.env$mrbin$parameters$verbose<-TRUE
            if(selection==selectTMPNo) mrbin.env$mrbin$parameters$verbose<-FALSE
           }
        }
        if(!stopTMP) selectStep<-selectStep+1
     }
     if(selectStep==-2){#Use parallel
      if(!mrbin.env$mrbin$parameters$tryParallel){
       selectTMPNo<-"Do not use parallel computing"
       selectTMPYes<-"Use parallel package for speed"
       selection<-utils::select.list(c(selectTMPYes,selectTMPNo),preselect=selectTMPYes,
                  title="Try parallel computing for speed?",graphics=TRUE)
       if(length(selection)==0|selection=="") stopTMP<-TRUE
       if(!stopTMP){
        if(selection==selectTMPYes) mrbin.env$mrbin$parameters$tryParallel<-TRUE
        if(selection==selectTMPNo) mrbin.env$mrbin$parameters$tryParallel<-FALSE
       }
      }
      if(!stopTMP) selectStep<-selectStep+1
     }
     if(selectStep==-1){#Show previews or not
      if(mrbin.env$mrbin$parameters$showSpectrumPreview=="No"){
       selectTMPNo<-"Do not show previews (e.g. for slow hardware)"
       selectTMPYes<-"Show spectrum previews (recommended)"
       selection<-utils::select.list(c(selectTMPYes,selectTMPNo#,"Go back"
                  ),preselect=selectTMPYes,
                  title="Show spectrum previews?",graphics=TRUE)
       if(length(selection)==0|selection=="") stopTMP<-TRUE
       if(!stopTMP){
          if(selection==selectTMPYes) mrbin.env$mrbin$parameters$showSpectrumPreview<-"Yes"
          if(selection==selectTMPNo) mrbin.env$mrbin$parameters$showSpectrumPreview<-"No"
       }
      }
      if(!stopTMP) selectStep<-selectStep+1
     }
     if(selectStep==0){#Set parameters
       selectionNewTMP<-NULL
       if(!is.null(mrbin.env$mrbin$parameters$NMRfolders)) selectionNewTMP<-c(
         selectionNewTMP,"Use current parameters without changes")
       selectionNewTMP<-c(selectionNewTMP,"Review parameters","Reload from file")
       selectionRepeat<-utils::select.list(c(selectionNewTMP),preselect="Review parameters",
                                          title="Edit parameters or use existing?",graphics=TRUE)
       if(length(selectionRepeat)==0|selectionRepeat=="") stopTMP<-TRUE
       if(selectionRepeat=="Reload from file"&!stopTMP){
         recreatemrbin()
         selectionRepeat2<-utils::select.list(c("Edit parameters","Use parameters from file without changes",
                                          "Go back"),
                                          preselect="Edit parameters",
                                          title="Edit parameters or use as is?",graphics=TRUE)
         if(length(selectionRepeat2)==0|selectionRepeat=="") stopTMP<-TRUE
         if(selectionRepeat2=="Review parameters"&!stopTMP) selectionRepeat<-"Review parameters"
         if(selectionRepeat2=="Use parameters from file without changes"&!stopTMP) selectionRepeat<-"Use current parameters"
       }
       if(!stopTMP&selectionRepeat=="Go back") selectStep<-selectStep-2
       if(!stopTMP) selectStep<-selectStep+1
     }
     if(selectionRepeat=="Use current parameters without changes"&!stopTMP){
       selectStep<-14
       lastStepDone<-TRUE
     }
     if(selectionRepeat=="Review parameters"&!stopTMP){
       if(selectStep==1){#1D or 2D data?
           dimension<-utils::select.list(c("1D","2D","Go back"),
                                   preselect=mrbin.env$mrbin$parameters$dimension,
                                   title="1D or 2D spectra?",graphics=TRUE)
           if(length(dimension)==0|dimension=="") stopTMP<-TRUE
           if(!stopTMP&!dimension=="Go back"){
             if(dimension%in%c("1D","2D")){
               mrbin.env$mrbin$parameters$dimension<-dimension
               if(dimension=="1D") dimlength<-2
               if(dimension=="2D") dimlength<-4
             }
           }
           if(!stopTMP&dimension=="Go back") selectStep<-selectStep-2
           if(!stopTMP) selectStep<-selectStep+1
       }
       if(selectStep==2){#Select folders
         if(!stopTMP){
           addFoldersTMP<-""
           if(length(mrbin.env$mrbin$parameters$NMRfolders)>0){
                selectionFoldersYes<-paste("Keep current list (",length(mrbin.env$mrbin$parameters$NMRfolders),
                                     " spectra)",sep="")
                selectionFolders<-utils::select.list(c("Create new spectra list",selectionFoldersYes,
                                  "Add or remove spectra from current list","Go back"),
                                  preselect=selectionFoldersYes,
                                  title="Use current spectra list?",graphics=TRUE)
                if(length(selectionFolders)==0|selectionFolders=="") stopTMP<-TRUE
                if(!stopTMP){
                  if(selectionFolders=="Create new spectra list"){
                    selectionFolders<-selectFolders()
                    if(selectionFolders=="stop")  stopTMP<-TRUE
                  }
                }
                if(!stopTMP){
                  if(selectionFolders=="Add or remove spectra from current list"){
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
                #if(selectionFolders=="")  stopTMP<-TRUE
           }
         }
         if(!stopTMP){
          if((selectionFolders=="Go back"|addFoldersTMP=="Go back")){
           selectStep<-selectStep-2
          } else {
           if(length(mrbin.env$mrbin$parameters$NMRfolders)<1) selectStep<-selectStep-1
          }
         }
         if(!stopTMP) selectStep<-selectStep+1
       }
       if(selectStep==3){
          if(!stopTMP){#load spectra for previews
            if(mrbin.env$mrbin$parameters$verbose){
                message("Hint: Review spectra to spot quality issues")
                utils::flush.console()
            }
            previewTMP<-utils::select.list(c("Review spectra","Do not review spectra","Go back"),
                       preselect="Review spectra",
                       ,title ="Review spectra?",graphics=TRUE)
            if(length(previewTMP)==0|previewTMP=="") stopTMP<-TRUE
            if(!stopTMP){
              if(!previewTMP=="Go back"){
                if(previewTMP=="Review spectra"){
                  ipreviewTMP<-1
                  while(ipreviewTMP <=length(mrbin.env$mrbin$parameters$NMRfolders)){
                    mrbin.env$mrbinTMP$currentFolder<-mrbin.env$mrbin$parameters$NMRfolders[ipreviewTMP]
                    readNMR2()
                    plotTitleTMP<-mrbin.env$mrbinTMP$currentFolder
                    if(nchar(plotTitleTMP)>56) plotTitleTMP<-paste("...",substr(
                      plotTitleTMP,nchar(plotTitleTMP)-52,nchar(plotTitleTMP)),sep="")
                    plotNMR(plotTitle=plotTitleTMP,region="all",manualScale=FALSE)
                    previewTMP2List<-NULL
                    if(ipreviewTMP<length(mrbin.env$mrbin$parameters$NMRfolders)){
                      nextTMP<-"Show next spectrum"
                    } else {
                      nextTMP<-"Continue"
                    }
                    previewTMP2List<-c(previewTMP2List,nextTMP)
                    if(ipreviewTMP>1) previewTMP2List<-c(previewTMP2List,"Show previous spectrum")
                    previewTMP2List<-c(previewTMP2List,"Skip review")
                    previewTMP2<-utils::select.list(previewTMP2List,
                       preselect=nextTMP,
                       ,title =paste("Spectrum ",ipreviewTMP," quality okay?",sep=""),graphics=TRUE)
                    if(length(previewTMP2)==0|previewTMP2==""){
                      ipreviewTMP<-length(mrbin.env$mrbin$parameters$NMRfolders)+1
                    } else {
                      if(previewTMP2==nextTMP) ipreviewTMP<-ipreviewTMP+1
                      if(previewTMP2=="Show previous spectrum") ipreviewTMP<-ipreviewTMP-1
                      if(previewTMP2=="Skip review") ipreviewTMP<-length(mrbin.env$mrbin$parameters$NMRfolders)+1
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
       if(selectStep==4){
          if(!stopTMP){#load spectra for previews
              if(mrbin.env$mrbin$parameters$verbose){
                message("Loading spectrum preview...")
                utils::flush.console()
              }
              mrbin.env$mrbinTMP$currentFolder<-mrbin.env$mrbin$parameters$NMRfolders[1]
              mrbin.env$mrbinTMP$timeEstimate<-max(.001,system.time(readNMR2())[1])
              if(mrbin.env$mrbin$parameters$showSpectrumPreview=="Yes"|mrbin.env$mrbin$parameters$PCA=="Yes"){
                #if(mrbin.env$mrbin$parameters$dimension=="2D"){
                #      #Set values below noise to zero, noise estimated from 90% percentile
                #      if(sum(mrbin.env$mrbinTMP$currentSpectrumOriginal>(3.5*sort(
                #      mrbin.env$mrbinTMP$currentSpectrumOriginal)[floor(length(
                #      mrbin.env$mrbinTMP$currentSpectrumOriginal)*.9)]))>300){
                #        mrbin.env$mrbinTMP$currentSpectrumOriginal[
                #        mrbin.env$mrbinTMP$currentSpectrumOriginal<=(3.5*sort(
                #        mrbin.env$mrbinTMP$currentSpectrumOriginal)[floor(length(
                #        mrbin.env$mrbinTMP$currentSpectrumOriginal)*.9)])]<-1e-8
                #      } else {
                #        mrbin.env$mrbinTMP$currentSpectrumOriginal[
                #        mrbin.env$mrbinTMP$currentSpectrumOriginal<(sort(
                #        mrbin.env$mrbinTMP$currentSpectrumOriginal,decreasing=TRUE)[300])]<-1e-8
                #      }
                #}
                mrbin.env$mrbinTMP$additionalPlots1D<-NULL
                mrbin.env$mrbinTMP$additionalPlots1DMetadata<-NULL
                mrbin.env$mrbinTMP$additionalPlots2D<-NULL
                mrbin.env$mrbinTMP$additionalPlots2DMetadata<-NULL
                #Find 3 more spectra: 33 percentile, 66 percentile, last spectrum
                mrbin.env$mrbinTMP$spectrumListPlotTMP<-setdiff(unique(c(
                  ceiling(length(mrbin.env$mrbin$parameters$NMRfolders)*.33),
                  ceiling(length(mrbin.env$mrbin$parameters$NMRfolders)*.66),
                  length(mrbin.env$mrbin$parameters$NMRfolders))),1)
                if(length(mrbin.env$mrbinTMP$spectrumListPlotTMP)>0){
                  #for 2D load and plot less spectra to save time
                  if(mrbin.env$mrbin$parameters$dimension=="2D"){
                    mrbin.env$mrbinTMP$spectrumListPlotTMP<-mrbin.env$mrbinTMP$spectrumListPlotTMP[1:
                      min(mrbin.env$mrbin$parameters$maxPreviewPlots2D-1,length(mrbin.env$mrbinTMP$spectrumListPlotTMP))]
                  }
                  for(ispectrumListPlotTMP in 1:length(mrbin.env$mrbinTMP$spectrumListPlotTMP)){
                    addToPlot(folder=mrbin.env$mrbin$parameters$NMRfolders[
                      mrbin.env$mrbinTMP$spectrumListPlotTMP[ispectrumListPlotTMP]],
                       dimension=mrbin.env$mrbin$parameters$dimension,
                       NMRvendor=mrbin.env$mrbin$parameters$NMRvendor,
                       useAsNames=mrbin.env$mrbin$parameters$useAsNames)
                  }
                  #if(mrbin.env$mrbin$parameters$dimension=="2D"){
                  #   for(iReducePlotsTMP in 1:length(mrbin.env$mrbinTMP$additionalPlots2D)){
                  #      #Set values below noise to zero, noise estimated from 90% percentile
                  #      if(sum(mrbin.env$mrbinTMP$additionalPlots2D[iReducePlotsTMP]>(3.5*sort(
                  #      mrbin.env$mrbinTMP$additionalPlots2D[iReducePlotsTMP])[floor(length(
                  #      mrbin.env$mrbinTMP$additionalPlots2D[iReducePlotsTMP])*.9)]))>300){
                  #        mrbin.env$mrbinTMP$additionalPlots2D[iReducePlotsTMP][
                  #        mrbin.env$mrbinTMP$additionalPlots2D[iReducePlotsTMP]<=(3.5*sort(
                  #        mrbin.env$mrbinTMP$additionalPlots2D[iReducePlotsTMP])[floor(length(
                  #        mrbin.env$mrbinTMP$additionalPlots2D[iReducePlotsTMP])*.9)])]<-1e-8
                  #      } else {
                  ##        mrbin.env$mrbinTMP$additionalPlots2D[iReducePlotsTMP][
                  #        mrbin.env$mrbinTMP$additionalPlots2D[iReducePlotsTMP]<(sort(
                  #        mrbin.env$mrbinTMP$additionalPlots2D[iReducePlotsTMP],decreasing=TRUE)[300])]<-1e-8
                  #      }
                  #   }
                  #}
                }
              }
          }
          if(!stopTMP){#Use rectangular bins or use special bin list, e.g. for lipids
            binMethodpreSelect<-mrbin.env$mrbin$parameters$binMethod
            userDefBinTMP<-"User defined bin list, e.g. for lipid analysis"
            if(binMethodpreSelect=="Custom bin list") binMethodpreSelect<-userDefBinTMP
            binMethod<-utils::select.list(c("Rectangular bins",userDefBinTMP,"Go back"),
                       preselect=binMethodpreSelect,
                       ,title ="Binning method: ",graphics=TRUE)
            if(length(binMethod)==0|binMethod=="") stopTMP<-TRUE
            if(!stopTMP){
              if(!binMethod=="Go back"){
                if(binMethod==userDefBinTMP) binMethod<-"Custom bin list"
                mrbin.env$mrbin$parameters$binMethod<-binMethod
                #Bin region
                adjRegion<-""
                if(mrbin.env$mrbin$parameters$binMethod=="Rectangular bins"){
                  if(mrbin.env$mrbin$parameters$verbose){
                    message("Hint: Include all visible peaks, excluding reference")
                    utils::flush.console()
                  }
                  accept<-FALSE
                  while(!accept&!stopTMP){
                    binRegionText<-paste(paste(c("left=","ppm, right=","ppm, top=","ppm, bottom=")[1:dimlength],
                                      mrbin.env$mrbin$parameters$binRegion[1:dimlength],collapse="",sep=""),"ppm",sep="")
                    if(mrbin.env$mrbin$parameters$showSpectrumPreview=="Yes"){
                       par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1))
                       plotMultiNMR(region="all",
                            rectangleRegions=matrix(mrbin.env$mrbin$parameters$binRegion,ncol=4),
                            color="black",manualScale=FALSE,maxPlots=2,
                            plotTitle=paste("Bin region\n",binRegionText,
                            sep=""))
                    }
                    adjRegion<-utils::select.list(c(paste("Keep: ",binRegionText,collapse=""),
                               "Change..."),preselect=paste("Use ",binRegionText,collapse=""),
                               title ="Bin region [ppm]: ",graphics=TRUE)
                    if(length(adjRegion)==0|adjRegion=="") stopTMP<-TRUE
                    if(!stopTMP){
                      if(adjRegion=="Change..."&!stopTMP){
                        regionTMP<-readline(prompt=paste("New left border, press enter to keep ",
                                  mrbin.env$mrbin$parameters$binRegion[1],": ",sep=""))
                        if(!regionTMP==""){
                            mrbin.env$mrbin$parameters$binRegion[1]<-as.numeric(regionTMP)
                        }
                        regionTMP<-readline(prompt=paste("New right border, press enter to keep ",
                                  mrbin.env$mrbin$parameters$binRegion[2],": ",sep=""))
                        if(!regionTMP==""){
                            mrbin.env$mrbin$parameters$binRegion[2]<-as.numeric(regionTMP)
                        }
                        if(mrbin.env$mrbin$parameters$dimension=="2D"){
                            regionTMP<-readline(prompt=paste("New top border, press enter to keep ",
                                      mrbin.env$mrbin$parameters$binRegion[3],": ",sep=""))
                            if(!regionTMP=="") {
                                mrbin.env$mrbin$parameters$binRegion[3]<-as.numeric(regionTMP)
                            }
                            regionTMP<-readline(prompt=paste("New bottom border, press enter to keep ",
                                      mrbin.env$mrbin$parameters$binRegion[4],": ",sep=""))
                            if(!regionTMP=="") {
                                mrbin.env$mrbin$parameters$binRegion[4]<-as.numeric(regionTMP)
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
        if(selectStep==5){ #Define bin width and height
          adjbinRegion<-""
          addbinRegion<-""
          if(mrbin.env$mrbin$parameters$dimension=="1D"&!stopTMP&mrbin.env$mrbin$parameters$binMethod=="Rectangular bins"){
              if(mrbin.env$mrbin$parameters$verbose){
                message("Hint: Should be broader than a single peak and include some margin. Gray\ncrosses indicate data point locations")
                utils::flush.console()
              }
              accept<-FALSE
              widthAdjust<-""
              while(!accept&!stopTMP&!widthAdjust=="Go back"){
                if(mrbin.env$mrbin$parameters$showSpectrumPreview=="Yes"){
                  par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1))
                  plotMultiNMR(region=mrbin.env$mrbin$parameters$previewRegion1D,
                        rectangleRegions=matrix(c((mrbin.env$mrbin$parameters$previewRegion1D[1]+
                          mrbin.env$mrbin$parameters$previewRegion1D[2])/2+as.numeric(mrbin.env$mrbin$parameters$binwidth1D)/2,
                          (mrbin.env$mrbin$parameters$previewRegion1D[1]+
                          mrbin.env$mrbin$parameters$previewRegion1D[2])/2-as.numeric(mrbin.env$mrbin$parameters$binwidth1D)/2,
                          21,21+1),ncol=4),
                        color="black", showGrid=TRUE,maxPlots=2,
                        manualScale=FALSE,
                        plotTitle=paste("Bin size\nwidth=",mrbin.env$mrbin$parameters$binwidth1D,"ppm",
                        sep=""),restrictToRange=TRUE)
                }
                binWidthTitle<-paste("Keep: width=",mrbin.env$mrbin$parameters$binwidth1D,"ppm",sep="")
                widthAdjust<-utils::select.list(c(binWidthTitle,"Change...",
                   "Show different part of spectrum...","Go back"),
                             preselect=binWidthTitle,
                             title ="Bin width [ppm]: ",graphics=TRUE)
                if(length(widthAdjust)==0|widthAdjust=="") stopTMP<-TRUE
                if(widthAdjust=="Show different part of spectrum..."){
                     widthTMP<-readline(prompt=paste("New left border, press enter to keep ",
                               mrbin.env$mrbin$parameters$previewRegion1D[1],": ",sep=""))
                     if(!widthTMP==""&!is.na(as.numeric(widthTMP))) {
                         mrbin.env$mrbin$parameters$previewRegion1D[1]<-as.numeric(widthTMP)
                     }
                     widthTMP<-readline(prompt=paste("New right border, press enter to keep ",
                               mrbin.env$mrbin$parameters$previewRegion1D[2],": ",sep=""))
                     if(!widthTMP==""&!is.na(as.numeric(widthTMP))) {
                         mrbin.env$mrbin$parameters$previewRegion1D[2]<-as.numeric(widthTMP)
                     }
                }
                if(widthAdjust=="Change..."){
                     widthTMP<-readline(prompt=paste("New 1D bin width, press enter to keep ",
                               mrbin.env$mrbin$parameters$binwidth1D,": ",sep=""))
                     if(!widthTMP==""&!is.na(as.numeric(widthTMP))) {
                         mrbin.env$mrbin$parameters$binwidth1D<-as.numeric(widthTMP)
                     }
                }
                if(widthAdjust==binWidthTitle) accept<-TRUE
              }
              if(widthAdjust=="Go back"){
                 selectStep<-selectStep-2
              }
              if(!stopTMP) selectStep<-selectStep+1
          }
          if(mrbin.env$mrbin$parameters$dimension=="2D"&!stopTMP&mrbin.env$mrbin$parameters$binMethod=="Rectangular bins"){
              if(mrbin.env$mrbin$parameters$verbose){
                message("Hint: Should be broader than a single peak and include some margin. Gray\ncrosses indicate data point locations")
                utils::flush.console()
              }
              accept<-FALSE
              widthAdjust<-""
              while(!accept&!stopTMP&!widthAdjust=="Go back"){
                if(mrbin.env$mrbin$parameters$showSpectrumPreview=="Yes") plotMultiNMR(region=
                        mrbin.env$mrbin$parameters$previewRegion2D,
                        rectangleRegions=matrix(c((mrbin.env$mrbin$parameters$previewRegion2D[1]+
                          mrbin.env$mrbin$parameters$previewRegion2D[2])/2+as.numeric(mrbin.env$mrbin$parameters$binwidth2D)/2,
                          (mrbin.env$mrbin$parameters$previewRegion2D[1]+
                          mrbin.env$mrbin$parameters$previewRegion2D[2])/2-as.numeric(mrbin.env$mrbin$parameters$binwidth2D)/2,
                          (mrbin.env$mrbin$parameters$previewRegion2D[3]+
                          mrbin.env$mrbin$parameters$previewRegion2D[4])/2+as.numeric(mrbin.env$mrbin$parameters$binheight)/2,
                          (mrbin.env$mrbin$parameters$previewRegion2D[3]+
                          mrbin.env$mrbin$parameters$previewRegion2D[4])/2-as.numeric(mrbin.env$mrbin$parameters$binheight)/2
                          ),ncol=4),
                        color="black",manualScale=FALSE, showGrid=TRUE,maxPlots=2,
                        plotTitle=paste("Bin size\nwidth=",mrbin.env$mrbin$parameters$binwidth2D,
                                  "ppm, height=",mrbin.env$mrbin$parameters$binheight,"ppm",sep=""),
                        restrictToRange=TRUE)
                currentBinSize<-paste("Keep: width=",mrbin.env$mrbin$parameters$binwidth2D,
                                  "ppm, height=",mrbin.env$mrbin$parameters$binheight,"ppm",sep="")
                widthAdjust<-utils::select.list(c(currentBinSize,"Change...",
                  "Show different part of spectrum...","Go back"),
                             preselect=currentBinSize,
                             title ="Bin size [ppm]: ",graphics=TRUE)
                if(length(widthAdjust)==0|widthAdjust=="") stopTMP<-TRUE
                if(widthAdjust==currentBinSize) accept<-TRUE
                if(widthAdjust=="Show different part of spectrum..."){
                     widthTMP<-readline(prompt=paste("New left border, press enter to keep ",
                               mrbin.env$mrbin$parameters$previewRegion2D[1],": ",sep=""))
                     if(!widthTMP==""&!is.na(as.numeric(widthTMP))) {
                         mrbin.env$mrbin$parameters$previewRegion2D[1]<-as.numeric(widthTMP)
                     }
                     widthTMP<-readline(prompt=paste("New right border, press enter to keep ",
                               mrbin.env$mrbin$parameters$previewRegion2D[2],": ",sep=""))
                     if(!widthTMP==""&!is.na(as.numeric(widthTMP))) {
                         mrbin.env$mrbin$parameters$previewRegion2D[2]<-as.numeric(widthTMP)
                     }
                     widthTMP<-readline(prompt=paste("New top border, press enter to keep ",
                               mrbin.env$mrbin$parameters$previewRegion2D[3],": ",sep=""))
                     if(!widthTMP==""&!is.na(as.numeric(widthTMP))) {
                         mrbin.env$mrbin$parameters$previewRegion2D[3]<-as.numeric(widthTMP)
                     }
                     widthTMP<-readline(prompt=paste("New bottom border, press enter to keep ",
                               mrbin.env$mrbin$parameters$previewRegion2D[4],": ",sep=""))
                     if(!widthTMP==""&!is.na(as.numeric(widthTMP))) {
                         mrbin.env$mrbin$parameters$previewRegion2D[4]<-as.numeric(widthTMP)
                     }
                }
                if(widthAdjust=="Change..."){
                     widthTMP<-readline(prompt=paste("New 2D bin width, press enter to keep ",
                               mrbin.env$mrbin$parameters$binwidth2D,": ",sep=""))
                     if(!widthTMP=="") {
                         #mrbin.env$paramChangeFlag<-TRUE
                         mrbin.env$mrbin$parameters$binwidth2D<-as.numeric(widthTMP)
                     }
                     heightTMP<-readline(prompt=paste("New 2D bin height, press enter to keep ",
                                mrbin.env$mrbin$parameters$binheight,": ",sep=""))
                     if(!heightTMP=="") {
                       #mrbin.env$paramChangeFlag<-TRUE
                       mrbin.env$mrbin$parameters$binheight<-as.numeric(heightTMP)
                     }
                }
              }
              if(widthAdjust=="Go back"){
                 selectStep<-selectStep-2
              }
              if(!stopTMP) selectStep<-selectStep+1
          }#Set custom bin list
          if(!stopTMP&mrbin.env$mrbin$parameters$binMethod=="Custom bin list"){
            adjbinRegion<-""
            if(!is.null(mrbin.env$mrbin$parameters$specialBinList)){
              if(nrow(mrbin.env$mrbin$parameters$specialBinList)==0) mrbin.env$mrbin$parameters$specialBinList<-NULL
            }
            adjbinRegionSelect<-""
            adjbinRegionAccept<-""
            if(!is.null(mrbin.env$mrbin$parameters$specialBinList)){
              if(nrow(mrbin.env$mrbin$parameters$specialBinList)==1){
                specialBinList_s<-""
                specialBinList_dots<-""
              } else {
                specialBinList_s<-"s"
                specialBinList_dots<-", ..."
                }
              keepbinRegionText<-paste(paste(c("left=","ppm, right=","ppm, top=","ppm, bottom=")[1:dimlength],
                          mrbin.env$mrbin$parameters$specialBinList[1:dimlength],collapse="",sep=""),"ppm",sep="")
              keepbinRegionYes<-paste("Keep previous bin list (",nrow(mrbin.env$mrbin$parameters$specialBinList),
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
            if(mrbin.env$mrbin$parameters$showSpectrumPreview=="Yes") plotMultiNMR(region="all",
                    rectangleRegions=mrbin.env$mrbin$parameters$specialBinList,color="black",
                    manualScale=FALSE,rectangleColors="green",maxPlots=2,
                    plotTitle=paste("Bin regions\n",sep=""),restrictToRange=TRUE)
            adjbinRegionSelect<-utils::select.list(c("Create new bin list",keepbinRegionYes,
                              editbinRegionYes,"Go back")[keepbinRegionIndex],
                              preselect=preselectbinRegion,
                              title="Create new bin list?",graphics=TRUE)
            if(length(adjbinRegionSelect)==0|adjbinRegionSelect=="") stopTMP<-TRUE
            if(!stopTMP){
              if(adjbinRegionSelect=="Create new bin list"){
                mrbin.env$mrbin$parameters$specialBinList<-NULL
              }
            }
            if(!stopTMP){
              if(adjbinRegionSelect=="Create new bin list"|adjbinRegionSelect==editbinRegionYes){
                if(mrbin.env$mrbin$parameters$verbose){
                  message("Hint: Should be broader than a single peak and include some margin")
                  utils::flush.console()
                }
                if(is.null(mrbin.env$mrbin$parameters$specialBinList)){
                        mrbin.env$mrbin$parameters$specialBinList<-matrix(ncol=4,nrow=0,dimnames=list(NULL,c("left","right","top","bottom")))
                }
                ibinRegions <- 1
                adjbinRegion<-""
                addbinRegion<-""
                while(ibinRegions <= (nrow(mrbin.env$mrbin$parameters$specialBinList)+1)&!stopTMP&
                      !adjbinRegion=="Go back"&!addbinRegion=="Go back"&!addbinRegion=="No"){
                  if(!stopTMP&!adjbinRegion=="Go back"){
                    if(ibinRegions>nrow(mrbin.env$mrbin$parameters$specialBinList)){
                      addbinRegion<-utils::select.list(c("Yes","No","Go back"),preselect="Yes",
                                 title ="Add a new bin?",graphics=TRUE)
                      if(length(addbinRegion)==0|addbinRegion==""){
                        stopTMP<-TRUE
                        addbinRegion<-""
                      }
                      if(!stopTMP){
                        if(addbinRegion=="Yes"){
                          mrbin.env$mrbin$parameters$specialBinList<-rbind(mrbin.env$mrbin$parameters$specialBinList,c(0,0,0,0))
                          if(nrow(mrbin.env$mrbin$parameters$specialBinList)==1){
                            rownames(mrbin.env$mrbin$parameters$specialBinList)<-""
                          } else {
                            rownames(mrbin.env$mrbin$parameters$specialBinList)[ibinRegions]<-""
                          }

                        }
                      }
                    }
                    if(!stopTMP&!adjbinRegion=="Go back"&!addbinRegion=="Go back"&!addbinRegion=="No"){
                      mean1<-mean(mrbin.env$mrbin$parameters$specialBinList[ibinRegions,1:2])
                      range1<-mrbin.env$mrbin$parameters$specialBinList[ibinRegions,1]-mrbin.env$mrbin$parameters$specialBinList[ibinRegions,2]
                      mean2<-mean(mrbin.env$mrbin$parameters$specialBinList[ibinRegions,3:4])
                      range2<-mrbin.env$mrbin$parameters$specialBinList[ibinRegions,4]-mrbin.env$mrbin$parameters$specialBinList[ibinRegions,3]
                      regionTMP<-c(mean1+3.5*range1,mean1-3.5*range1,mean2-3.5*range2,mean2+3.5*range2)
                      showGridTMP<-TRUE
                      if(sum(mrbin.env$mrbin$parameters$specialBinList[ibinRegions,]==0)==4){
                        regionTMP<-"all"
                        showGridTMP<-FALSE
                      }
                      if(mrbin.env$mrbin$parameters$dimension=="1D"){
                        if(!range1==0){
                          if(mrbin.env$mrbin$parameters$showSpectrumPreview=="Yes") plotMultiNMR(region=
                                  regionTMP,
                                  rectangleRegions=matrix(c(mrbin.env$mrbin$parameters$specialBinList[ibinRegions,1],
                                      mrbin.env$mrbin$parameters$specialBinList[ibinRegions,2],0,2),ncol=4),
                                  color="black", showGrid=showGridTMP,manualScale=FALSE,maxPlots=2,
                                  plotTitle=paste("Custom bins\nleft=",mrbin.env$mrbin$parameters$specialBinList[ibinRegions,1],
                                    "ppm, right=",mrbin.env$mrbin$parameters$specialBinList[ibinRegions,2],"ppm",sep=""),
                                  restrictToRange=TRUE)
                        }
                      }
                      if(mrbin.env$mrbin$parameters$dimension=="2D"){
                        if(!range1==0&!range2==0){
                          if(mrbin.env$mrbin$parameters$showSpectrumPreview=="Yes") plotMultiNMR(region=regionTMP,
                                  rectangleRegions=matrix(mrbin.env$mrbin$parameters$specialBinList[ibinRegions,1:4],ncol=4),
                                  color="black",manualScale=FALSE, showGrid=showGridTMP,maxPlots=2,
                                  plotTitle=paste("Custom bins\nleft=",mrbin.env$mrbin$parameters$specialBinList[ibinRegions,1],
                                    "ppm, right=",mrbin.env$mrbin$parameters$specialBinList[ibinRegions,2],"ppm",sep=""),
                                  restrictToRange=TRUE)
                        }
                      }
                      adjbinRegionAccept<-paste(paste(c("Keep left=","ppm, right=","ppm, top=","ppm, bottom=")[1:dimlength],
                                        mrbin.env$mrbin$parameters$specialBinList[ibinRegions,1:dimlength],collapse="",sep=""),"ppm",sep="")
                      if(sum(mrbin.env$mrbin$parameters$specialBinList[ibinRegions,]==0)==4){
                        adjbinRegion<-"Change..."
                      } else {
                        if(rownames(mrbin.env$mrbin$parameters$specialBinList)[ibinRegions]==""){
                          binTitleTMP<-""
                        } else {
                          binTitleTMP<-paste(" (\"",rownames(mrbin.env$mrbin$parameters$specialBinList)[ibinRegions],"\")",sep="")
                        }
                        adjbinRegion<-utils::select.list(c(adjbinRegionAccept,"Change...","Remove bin","Go back"),
                                 preselect=adjbinRegionAccept,
                                 title =paste("Edit bin ",ibinRegions,binTitleTMP,"?",sep=""),graphics=TRUE)
                      }
                      if(length(adjbinRegion)==0|adjbinRegion=="") stopTMP<-TRUE
                    }
                  }
                  if(adjbinRegion=="Change..."&!stopTMP&!adjbinRegion=="Go back"){

                    if(rownames(mrbin.env$mrbin$parameters$specialBinList)[ibinRegions]==""){
                      promptTMP<-paste("New bin name, press enter for no name: ",sep="")
                    } else {
                      promptTMP<-paste("New bin name, press enter to keep ",
                                rownames(mrbin.env$mrbin$parameters$specialBinList)[ibinRegions],": ",sep="")
                    }
                    nameTMP<-readline(prompt=promptTMP)
                    if(!nameTMP=="") {
                           rownames(mrbin.env$mrbin$parameters$specialBinList)[ibinRegions]<-nameTMP
                    }
                    regionTMP<-readline(prompt=paste("Bin ",ibinRegions,": left border, press enter to keep ",
                              mrbin.env$mrbin$parameters$specialBinList[ibinRegions,1],": ",sep=""))
                    if(!regionTMP=="") {
                           mrbin.env$mrbin$parameters$specialBinList[ibinRegions,1]<-as.numeric(regionTMP)
                    }
                    regionTMP<-readline(prompt=paste("Bin ",ibinRegions,": right border, press enter to keep ",
                              mrbin.env$mrbin$parameters$specialBinList[ibinRegions,2],": ",sep=""))
                    if(!regionTMP=="") {
                           mrbin.env$mrbin$parameters$specialBinList[ibinRegions,2]<-as.numeric(regionTMP)
                    }
                    if(mrbin.env$mrbin$parameters$specialBinList[ibinRegions,1]<mrbin.env$mrbin$parameters$specialBinList[ibinRegions,2]){
                      TMP<-mrbin.env$mrbin$parameters$specialBinList[ibinRegions,1]
                      mrbin.env$mrbin$parameters$specialBinList[ibinRegions,1]<-mrbin.env$mrbin$parameters$specialBinList[ibinRegions,2]
                      mrbin.env$mrbin$parameters$specialBinList[ibinRegions,2]<-TMP
                    }
                    if(mrbin.env$mrbin$parameters$dimension=="2D"&!stopTMP){
                      regionTMP<-readline(prompt=paste("Bin ",ibinRegions,": top border, press enter to keep ",
                                mrbin.env$mrbin$parameters$specialBinList[ibinRegions,3],": ",sep=""))
                      if(!regionTMP=="") {
                             mrbin.env$mrbin$parameters$specialBinList[ibinRegions,3]<-as.numeric(regionTMP)
                      }
                      regionTMP<-readline(prompt=paste("Bin ",ibinRegions,": bottom border, press enter to keep ",
                                mrbin.env$mrbin$parameters$specialBinList[ibinRegions,4],": ",sep=""))
                      if(!regionTMP=="") {
                             mrbin.env$mrbin$parameters$specialBinList[ibinRegions,4]<-as.numeric(regionTMP)
                      }
                      if(mrbin.env$mrbin$parameters$specialBinList[ibinRegions,4]<mrbin.env$mrbin$parameters$specialBinList[ibinRegions,3]){
                        TMP<-mrbin.env$mrbin$parameters$specialBinList[ibinRegions,3]
                        mrbin.env$mrbin$parameters$specialBinList[ibinRegions,3]<-mrbin.env$mrbin$parameters$specialBinList[ibinRegions,4]
                        mrbin.env$mrbin$parameters$specialBinList[ibinRegions,4]<-TMP
                      }
                    }
                  }
                  if(adjbinRegion=="Remove bin"&!stopTMP&!adjbinRegion=="Go back"){
                    mrbin.env$mrbin$parameters$specialBinList<-mrbin.env$mrbin$parameters$specialBinList[-ibinRegions,,drop=FALSE]
                  }
                  if(adjbinRegion==adjbinRegionAccept&!stopTMP&!adjbinRegion=="Go back"){
                    ibinRegions <- ibinRegions+1
                  }
                }
                if(nrow(mrbin.env$mrbin$parameters$specialBinList)==0){
                  mrbin.env$mrbin$parameters$specialBinList<-NULL
                  adjbinRegion<-"Go back"
                }
              }
            }
            if(!stopTMP&!addbinRegion=="Go back"&!adjbinRegion=="Go back"){
              if(is.null(mrbin.env$mrbin$parameters$specialBinList)){
                addbinRegion<-"Go back"
              }
            }
            if(adjbinRegion=="Go back"|addbinRegion=="Go back"|adjbinRegionSelect=="Go back"){
               selectStep<-selectStep-2
            }
            if(!stopTMP) selectStep<-selectStep+1
          }
        }
        if(selectStep==6){#Scale to reference
          if(!stopTMP){
            adjRegion<-""
            referenceScaling<-utils::select.list(c("Yes","No","Go back"),
                                     preselect=mrbin.env$mrbin$parameters$referenceScaling,
                                          title = "Scale to reference signal?",graphics=TRUE)
            if(length(referenceScaling)==0|referenceScaling=="") stopTMP<-TRUE
            if(!stopTMP&!referenceScaling=="Go back"){
              mrbin.env$mrbin$parameters$referenceScaling<-referenceScaling
              if(mrbin.env$mrbin$parameters$referenceScaling=="Yes"){
                if(mrbin.env$mrbin$parameters$verbose){
                  message("Hint: Include reference signal with some margin")
                  utils::flush.console()
                }
                if(mrbin.env$mrbin$parameters$dimension=="1D"){
                accept<-FALSE
                adjRegion<-""
                while(!accept&!stopTMP&!adjRegion=="Go back"){
                  mean1<-mean(mrbin.env$mrbin$parameters$reference1D)
                  range1<-max(mrbin.env$mrbin$parameters$reference1D)-min(mrbin.env$mrbin$parameters$reference1D)
                  if(mrbin.env$mrbin$parameters$showSpectrumPreview=="Yes") plotMultiNMR(
                          region=c(mean1+4*range1,mean1-4*range1,-10,10),
                          rectangleRegions=matrix(c(mrbin.env$mrbin$parameters$reference1D[1],
                                                  mrbin.env$mrbin$parameters$reference1D[2],0,2),ncol=4),
                          color="black",manualScale=FALSE,restrictToRange=TRUE,maxPlots=2,
                          plotTitle=paste("Reference region\nleft=",mrbin.env$mrbin$parameters$reference1D[1],
                                    "ppm, right=",mrbin.env$mrbin$parameters$reference1D[2],"ppm",sep=""))
                  RefRegionTitle<-paste("Keep: left=",mrbin.env$mrbin$parameters$reference1D[1],
                                        "ppm, right=",mrbin.env$mrbin$parameters$reference1D[2],"ppm",sep="")
                  adjRegion<-utils::select.list(c(RefRegionTitle,
                             "Change...","Go back"),
                             preselect=RefRegionTitle,title ="Reference region [ppm]: ",graphics=TRUE)
                  if(length(adjRegion)==0|adjRegion=="") stopTMP<-TRUE
                  if(adjRegion==RefRegionTitle) accept<-TRUE
                  if(adjRegion=="Change..."){
                    regionTMP<-readline(prompt=paste("New left border, press enter to keep ",
                              mrbin.env$mrbin$parameters$reference1D[1],": ",sep=""))
                    if(!regionTMP=="") {
                           mrbin.env$mrbin$parameters$reference1D[1]<-as.numeric(regionTMP)
                    }
                    regionTMP<-readline(prompt=paste("New right border, press enter to keep ",
                              mrbin.env$mrbin$parameters$reference1D[2],": ",sep=""))
                    if(!regionTMP=="") {
                           mrbin.env$mrbin$parameters$reference1D[2]<-as.numeric(regionTMP)
                    }
                  }
                  }
                }
                if(mrbin.env$mrbin$parameters$dimension=="2D"){
                  accept<-FALSE
                  adjRegion<-""
                  while(!accept&!stopTMP&!adjRegion=="Go back"){
                    mean1<-mean(mrbin.env$mrbin$parameters$reference2D[1:2])
                    range1<-max(mrbin.env$mrbin$parameters$reference2D[1:2])-min(mrbin.env$mrbin$parameters$reference2D[1:2])
                    mean2<-mean(mrbin.env$mrbin$parameters$reference2D[3:4])
                    range2<-max(mrbin.env$mrbin$parameters$reference2D[3:4])-min(mrbin.env$mrbin$parameters$reference2D[3:4])
                    if(mrbin.env$mrbin$parameters$showSpectrumPreview=="Yes") plotMultiNMR(
                            region=c(mean1+4*range1,mean1-4*range1,
                                     mean2-4*range2,mean2+4*range2),
                            rectangleRegions=matrix(mrbin.env$mrbin$parameters$reference2D,ncol=4),
                            color="black",manualScale=FALSE,restrictToRange=TRUE,maxPlots=2,
                            plotTitle=paste("Reference region\nleft=",mrbin.env$mrbin$parameters$reference2D[1],
                                      "ppm, right=",mrbin.env$mrbin$parameters$reference2D[2],
                                      "ppm, top=",mrbin.env$mrbin$parameters$reference2D[3],
                                      "ppm, bottom=",mrbin.env$mrbin$parameters$reference2D[4],"ppm",sep=""))
                    RefRegionTitle<-paste("Keep: left=",mrbin.env$mrbin$parameters$reference2D[1],
                                          "ppm, right=",mrbin.env$mrbin$parameters$reference2D[2],
                                          "ppm, top=",mrbin.env$mrbin$parameters$reference2D[3],
                                          "ppm, bottom=",mrbin.env$mrbin$parameters$reference2D[4],"ppm",sep="")
                    adjRegion<-utils::select.list(c(RefRegionTitle,
                               "Change...","Go back"),
                               preselect=RefRegionTitle,title ="Reference region [ppm]: ",graphics=TRUE)
                    if(length(adjRegion)==0|adjRegion=="") stopTMP<-TRUE
                    if(adjRegion==RefRegionTitle) accept<-TRUE
                    if(adjRegion=="Change..."){
                      regionTMP<-readline(prompt=paste("New left border, press enter to keep ",
                                mrbin.env$mrbin$parameters$reference2D[1],": ",sep=""))
                      if(!regionTMP=="") {
                             mrbin.env$mrbin$parameters$reference2D[1]<-as.numeric(regionTMP)
                      }
                      regionTMP<-readline(prompt=paste("New right border, press enter to keep ",
                                mrbin.env$mrbin$parameters$reference2D[2],": ",sep=""))
                      if(!regionTMP=="") {
                             mrbin.env$mrbin$parameters$reference2D[2]<-as.numeric(regionTMP)
                      }
                      regionTMP<-readline(prompt=paste("New top border, press enter to keep ",
                                mrbin.env$mrbin$parameters$reference2D[3],": ",sep=""))
                      if(!regionTMP=="") {
                             mrbin.env$mrbin$parameters$reference2D[3]<-as.numeric(regionTMP)
                      }
                      regionTMP<-readline(prompt=paste("New bottom border, press enter to keep ",
                                mrbin.env$mrbin$parameters$reference2D[4],": ",sep=""))
                      if(!regionTMP=="") {
                             mrbin.env$mrbin$parameters$reference2D[4]<-as.numeric(regionTMP)
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
        if(selectStep==7){#Remove solvent
          if(!stopTMP){
            adjRegion<-""
            removeSolvent<-utils::select.list(c("Yes","No","Go back"),
                                     preselect=mrbin.env$mrbin$parameters$removeSolvent,
                                     title = "Remove solvent area?",graphics=TRUE)
            if(length(removeSolvent)==0|removeSolvent=="") stopTMP<-TRUE
            if(!stopTMP&!removeSolvent=="Go back"){
              mrbin.env$mrbin$parameters$removeSolvent<-removeSolvent
              if(mrbin.env$mrbin$parameters$removeSolvent=="Yes"){
                if(mrbin.env$mrbin$parameters$verbose){
                  message("Hint: Include solvent signal and some margin")
                  utils::flush.console()
                }
                accept<-FALSE
                adjRegion<-""
                while(!accept&!stopTMP&!adjRegion=="Go back"){
                  mean1<-mean(mrbin.env$mrbin$parameters$solventRegion[1:2])
                  range1<-max(mrbin.env$mrbin$parameters$solventRegion[1:2])-min(mrbin.env$mrbin$parameters$solventRegion[1:2])
                  if(mrbin.env$mrbin$parameters$showSpectrumPreview=="Yes") plotMultiNMR(
                          region=c(mean1+6*range1,mean1-6*range1,-10,160),rectangleColors="orange",
                          rectangleRegions=matrix(c(mrbin.env$mrbin$parameters$solventRegion[1],
                                                  mrbin.env$mrbin$parameters$solventRegion[2],-1000,1000),ncol=4),
                          color="black",manualScale=FALSE,restrictToRange=TRUE,maxPlots=2,
                          plotTitle=paste("Solvent region\nleft=",mrbin.env$mrbin$parameters$solventRegion[1],
                                    "ppm, right=",mrbin.env$mrbin$parameters$solventRegion[2],"ppm",sep=""))
                  SolventRegionTitle<-paste("Keep: left=",mrbin.env$mrbin$parameters$solventRegion[1],
                                        "ppm, right=",mrbin.env$mrbin$parameters$solventRegion[2],"ppm",sep="")
                  adjRegion<-utils::select.list(c(SolventRegionTitle,
                             "Change...","Go back"),
                             preselect=SolventRegionTitle,title ="Solvent region to be removed: ",graphics=TRUE)
                  if(length(adjRegion)==0|adjRegion=="") stopTMP<-TRUE
                  if(adjRegion==SolventRegionTitle) accept=TRUE
                  if(adjRegion=="Change..."&!stopTMP){
                    regionTMP<-readline(prompt=paste("New left border, press enter to keep ",
                              mrbin.env$mrbin$parameters$solventRegion[1],": ",sep=""))
                    if(!regionTMP=="") {
                           mrbin.env$mrbin$parameters$solventRegion[1]<-as.numeric(regionTMP)
                    }
                    regionTMP<-readline(prompt=paste("New right border, press enter to keep ",
                              mrbin.env$mrbin$parameters$solventRegion[2],": ",sep=""))
                    if(!regionTMP=="") {
                           mrbin.env$mrbin$parameters$solventRegion[2]<-as.numeric(regionTMP)
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
        if(selectStep==8){#Remove additional areas
          if(!stopTMP){
            removeAreaListTMP<-""
            adjbinRegion<-""
            addbinRegion<-""
            if(mrbin.env$mrbin$parameters$verbose){
              message("Hint: Remove spectral artifacts and solvent and contaminant signals")
              utils::flush.console()
            }
            removeAreas<-utils::select.list(c("Yes","No","Go back"),preselect=mrbin.env$mrbin$parameters$removeAreas,
                                     title = "Remove additional areas?",graphics=TRUE)
            if(length(removeAreas)==0|removeAreas=="") stopTMP<-TRUE
            if(!stopTMP&!removeAreas=="Go back"){
              mrbin.env$mrbin$parameters$removeAreas<-removeAreas
              if(mrbin.env$mrbin$parameters$removeAreas=="Yes"){
                addAreasFlag<-TRUE
                if(!is.null( mrbin.env$mrbin$parameters$removeAreaList)){
                  if(nrow(mrbin.env$mrbin$parameters$removeAreaList)==0){
                    mrbin.env$mrbin$parameters$removeAreaList<-NULL
                  }
                }
                if(!is.null(mrbin.env$mrbin$parameters$removeAreaList)){
                  if(nrow(mrbin.env$mrbin$parameters$removeAreaList)>0){
                      addAreasFlag<-FALSE
                      if(nrow(mrbin.env$mrbin$parameters$removeAreaList)==1){
                        regions_s<-""
                        regions_dots<-""
                      }
                      if(nrow(mrbin.env$mrbin$parameters$removeAreaList)>1){
                        regions_s<-"s"
                        regions_dots<-", ..."
                      }
                      preselectKeepTMP<-paste("Keep current list (",nrow(mrbin.env$mrbin$parameters$removeAreaList)," region",regions_s,", ",
                                        paste(c("left=","ppm, right=","ppm, top=","ppm ,bottom=")[1:dimlength],
                                          mrbin.env$mrbin$parameters$removeAreaList[1,1:dimlength],
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
                if(mrbin.env$mrbin$parameters$showSpectrumPreview=="Yes") plotMultiNMR(region="all",
                  rectangleRegions=mrbin.env$mrbin$parameters$removeAreaList,color="black",
                        manualScale=FALSE,rectangleColors="orange",maxPlots=2,
                        plotTitle=paste("Removed areas\n",
                        sep=""))
                removeAreaListTMP<-utils::select.list(c("Create new list",preselectKeepTMP,"Edit current list","Go back")[keepbinRegionIndex],
                                   preselect=preselectKeepTMPYes,
                                   title = "Use previous area list or define new?",graphics=TRUE)
                if(length(removeAreaListTMP)==0|removeAreaListTMP=="") stopTMP<-TRUE
                if(!removeAreaListTMP==preselectKeepTMP&!stopTMP&!removeAreaListTMP=="Go back"){
                  addAreasFlag<-TRUE
                  if(removeAreaListTMP=="Create new list"&!stopTMP){
                    mrbin.env$mrbin$parameters$removeAreaList<-matrix(ncol=4,nrow=0,dimnames=list(NULL,c("left","right","top","bottom")))
                  }
                }
                if(!stopTMP){
                  if(removeAreaListTMP=="Create new list"|removeAreaListTMP=="Edit current list"){
                    if(is.null(mrbin.env$mrbin$parameters$removeAreaList)){
                            mrbin.env$mrbin$parameters$removeAreaList<-matrix(ncol=4,nrow=0,dimnames=list(NULL,c("left","right","top","bottom")))
                    }
                    ibinRegions <- 1
                    adjbinRegion<-""
                    addbinRegion<-""
                    adjbinRegionAccept<-""
                    while(ibinRegions <= (nrow(mrbin.env$mrbin$parameters$removeAreaList)+1)&!stopTMP&
                          !adjbinRegion=="Go back"&!addbinRegion=="Go back"&!addbinRegion=="No"){
                      if(!stopTMP&!adjbinRegion=="Go back"){
                        if(ibinRegions>nrow(mrbin.env$mrbin$parameters$removeAreaList)){
                          addbinRegion<-utils::select.list(c("Yes","No","Go back"),preselect="No",
                                     title ="Add a new region?",graphics=TRUE)
                          if(length(addbinRegion)==0|addbinRegion==""){
                            stopTMP<-TRUE
                            addbinRegion<-""
                          }
                          if(!stopTMP){
                            if(addbinRegion=="Yes"){
                              mrbin.env$mrbin$parameters$removeAreaList<-rbind(mrbin.env$mrbin$parameters$removeAreaList,c(0,0,0,0))
                            }
                          }
                        }
                        if(!stopTMP&!adjbinRegion=="Go back"&!addbinRegion=="Go back"&!addbinRegion=="No"){
                          mean1<-mean(mrbin.env$mrbin$parameters$removeAreaList[ibinRegions,1:2])
                          range1<-mrbin.env$mrbin$parameters$removeAreaList[ibinRegions,1]-mrbin.env$mrbin$parameters$removeAreaList[ibinRegions,2]
                          mean2<-mean(mrbin.env$mrbin$parameters$removeAreaList[ibinRegions,3:4])
                          range2<-mrbin.env$mrbin$parameters$removeAreaList[ibinRegions,4]-mrbin.env$mrbin$parameters$removeAreaList[ibinRegions,3]
                          regionTMP<-c(mean1+4*range1,mean1-4*range1,mean2-4*range2,mean2+4*range2)
                          if(sum(mrbin.env$mrbin$parameters$removeAreaList[ibinRegions,]==0)==4) regionTMP<-"all"
                          if(mrbin.env$mrbin$parameters$dimension=="1D"){
                            if(mrbin.env$mrbin$parameters$showSpectrumPreview=="Yes") plotMultiNMR(
                                    region=regionTMP,rectangleColors="orange",
                                    rectangleRegions=matrix(c(mrbin.env$mrbin$parameters$removeAreaList[ibinRegions,1],
                                                            mrbin.env$mrbin$parameters$removeAreaList[ibinRegions,2],0,2),ncol=4),
                                    color="black",manualScale=FALSE,maxPlots=2,
                                    plotTitle=paste("Remove area\nleft=",mrbin.env$mrbin$parameters$removeAreaList[ibinRegions,1],
                                      "ppm, right=",mrbin.env$mrbin$parameters$removeAreaList[ibinRegions,2],"ppm",sep=""),
                                    restrictToRange=TRUE)
                          }
                          if(mrbin.env$mrbin$parameters$dimension=="2D"){
                            if(mrbin.env$mrbin$parameters$showSpectrumPreview=="Yes") plotMultiNMR(
                                    region=regionTMP,rectangleColors="orange",
                                    rectangleRegions=matrix(mrbin.env$mrbin$parameters$removeAreaList[ibinRegions,1:4],ncol=4),
                                    color="black",manualScale=FALSE,maxPlots=2,
                                    plotTitle=paste("Remove area\nleft=",mrbin.env$mrbin$parameters$removeAreaList[ibinRegions,1],
                                    "ppm, right=",mrbin.env$mrbin$parameters$removeAreaList[ibinRegions,2],"ppm",sep=""),
                                    restrictToRange=TRUE)
                          }
                          adjbinRegionAccept<-paste(paste(c("Keep left=","ppm, right=","ppm, top=","ppm, bottom=")[1:dimlength],
                                            mrbin.env$mrbin$parameters$removeAreaList[ibinRegions,1:dimlength],collapse="",sep=""),"ppm",sep="")
                          if(sum(mrbin.env$mrbin$parameters$removeAreaList[ibinRegions,]==0)==4){
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
                                  mrbin.env$mrbin$parameters$removeAreaList[ibinRegions,1],": ",sep=""))
                        if(!regionTMP=="") {
                               mrbin.env$mrbin$parameters$removeAreaList[ibinRegions,1]<-as.numeric(regionTMP)
                        }
                        regionTMP<-readline(prompt=paste("Region ",ibinRegions,": right border, press enter to keep ",
                                  mrbin.env$mrbin$parameters$removeAreaList[ibinRegions,2],": ",sep=""))
                        if(!regionTMP=="") {
                               mrbin.env$mrbin$parameters$removeAreaList[ibinRegions,2]<-as.numeric(regionTMP)
                        }
                        if(mrbin.env$mrbin$parameters$removeAreaList[ibinRegions,1]<mrbin.env$mrbin$parameters$removeAreaList[ibinRegions,2]){
                          TMP<-mrbin.env$mrbin$parameters$removeAreaList[ibinRegions,1]
                          mrbin.env$mrbin$parameters$removeAreaList[ibinRegions,1]<-mrbin.env$mrbin$parameters$removeAreaList[ibinRegions,2]
                          mrbin.env$mrbin$parameters$removeAreaList[ibinRegions,2]<-TMP
                        }
                      if(mrbin.env$mrbin$parameters$dimension=="2D"&!stopTMP){
                        regionTMP<-readline(prompt=paste("Region ",ibinRegions,": top border, press enter to keep ",
                                  mrbin.env$mrbin$parameters$removeAreaList[ibinRegions,3],": ",sep=""))
                        if(!regionTMP=="") {
                               mrbin.env$mrbin$parameters$removeAreaList[ibinRegions,3]<-as.numeric(regionTMP)
                        }
                        regionTMP<-readline(prompt=paste("Region ",ibinRegions,": bottom border, press enter to keep ",
                                  mrbin.env$mrbin$parameters$removeAreaList[ibinRegions,4],": ",sep=""))
                        if(!regionTMP=="") {
                               mrbin.env$mrbin$parameters$removeAreaList[ibinRegions,4]<-as.numeric(regionTMP)
                        }
                      if(mrbin.env$mrbin$parameters$removeAreaList[ibinRegions,4]<mrbin.env$mrbin$parameters$removeAreaList[ibinRegions,3]){
                        TMP<-mrbin.env$mrbin$parameters$removeAreaList[ibinRegions,3]
                        mrbin.env$mrbin$parameters$removeAreaList[ibinRegions,3]<-mrbin.env$mrbin$parameters$removeAreaList[ibinRegions,4]
                        mrbin.env$mrbin$parameters$removeAreaList[ibinRegions,4]<-TMP
                      }
                     }
                    }
                    if(adjbinRegion=="Remove entry"&!stopTMP&!adjbinRegion=="Go back"){
                      mrbin.env$mrbin$parameters$removeAreaList<-mrbin.env$mrbin$parameters$removeAreaList[-ibinRegions,,drop=FALSE]
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
          if(!is.null(mrbin.env$mrbin$parameters$removeAreaList)){
            if(nrow(mrbin.env$mrbin$parameters$removeAreaList)==0){
              mrbin.env$mrbin$parameters$removeAreaList<-NULL
              adjbinRegion<-"Go back"
            }
          }
          if(adjbinRegion=="Go back"|addbinRegion=="Go back"|removeAreaListTMP=="Go back"|removeAreas=="Go back"){
             selectStep<-selectStep-2
          }
          if(!stopTMP) selectStep<-selectStep+1
        }
        if(selectStep==9){#Merge bins containing unstable peaks
          if(!stopTMP){
            if(mrbin.env$mrbin$parameters$verbose){
              message("Hint: Signals that differ in chemical shift from sample to sample")
              utils::flush.console()
            }
            sumBinListTMP<-""
            removeAreaListTMP<-""
            adjbinRegion<-""
            addbinRegion<-""
            sumBins<-utils::select.list(c("Merge bins of unstable peaks (e.g. citrate)","No","Go back"),
                                 preselect=mrbin.env$mrbin$parameters$sumBins,
                                 title = "Merge bins of unstable peaks?",graphics=TRUE)
            if(length(sumBins)==0|sumBins=="") stopTMP<-TRUE
            if(!stopTMP&!sumBins=="Go back"){
              if(sumBins=="Merge bins of unstable peaks (e.g. citrate)"){
                mrbin.env$mrbin$parameters$sumBins<-"Yes"
              } else {
                mrbin.env$mrbin$parameters$sumBins<-sumBins
              }
              if(mrbin.env$mrbin$parameters$sumBins=="Yes"&!stopTMP){
                  addAreasFlag<-TRUE
                  if(!is.null( mrbin.env$mrbin$parameters$sumBinList)){
                    if(nrow(mrbin.env$mrbin$parameters$sumBinList)==0){
                      mrbin.env$mrbin$parameters$sumBinList<-NULL
                    }
                  }
                  if(!is.null( mrbin.env$mrbin$parameters$sumBinList)){
                    if(nrow(mrbin.env$mrbin$parameters$sumBinList)>0&!stopTMP){
                        addAreasFlag<-FALSE
                        if(nrow(mrbin.env$mrbin$parameters$sumBinList)==1){
                          regions_s<-""
                          regions_dots<-""
                        }
                        if(nrow(mrbin.env$mrbin$parameters$sumBinList)>1){
                          regions_s<-"s"
                          regions_dots<-", ..."
                        }
                        preselectKeepTMP<-paste("Keep current list (",nrow(mrbin.env$mrbin$parameters$sumBinList)," region",regions_s,", ",
                                          paste(c("left=","ppm, right=","ppm, top=","ppm ,bottom=")[1:dimlength],
                                            mrbin.env$mrbin$parameters$sumBinList[1,1:dimlength],
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
                  if(mrbin.env$mrbin$parameters$showSpectrumPreview=="Yes") plotMultiNMR(region="all",
                          rectangleRegions=mrbin.env$mrbin$parameters$sumBinList,color="black",
                          manualScale=FALSE,rectangleColors="green",maxPlots=2,
                          plotTitle=paste("Summed areas\n",sep=""),restrictToRange=TRUE)
                  sumBinListTMP<-utils::select.list(c("Create new list",preselectKeepTMP,"Edit current list","Go back")[keepbinRegionIndex],
                                 preselect=preselectKeepTMPYes,
                                 title = "Use previous area list or define new?",graphics=TRUE)
                  if(length(sumBinListTMP)==0|sumBinListTMP=="") stopTMP<-TRUE
                  if(!sumBinListTMP==preselectKeepTMP&!stopTMP&!sumBinListTMP=="Go back"){
                    addAreasFlag<-TRUE
                    if(sumBinListTMP=="Create new list"&!stopTMP){
                        mrbin.env$mrbin$parameters$sumBinList<-matrix(ncol=4,nrow=0,dimnames=list(NULL,c("left","right","top","bottom")))
                    }
                  }
                if(!stopTMP){
                  if(sumBinListTMP=="Create new list"|sumBinListTMP=="Edit current list"){
                    if(is.null(mrbin.env$mrbin$parameters$sumBinList)){
                            mrbin.env$mrbin$parameters$sumBinList<-matrix(ncol=4,nrow=0,dimnames=list(NULL,c("left","right","top","bottom")))
                    }
                    ibinRegions <- 1
                    adjbinRegion<-""
                    addbinRegion<-""
                    adjbinRegionAccept<-""
                    while(ibinRegions <= (nrow(mrbin.env$mrbin$parameters$sumBinList)+1)&!stopTMP&
                          !adjbinRegion=="Go back"&!addbinRegion=="Go back"&!addbinRegion=="No"){
                      if(!stopTMP&!adjbinRegion=="Go back"){
                        if(ibinRegions>nrow(mrbin.env$mrbin$parameters$sumBinList)){
                          addbinRegion<-utils::select.list(c("Yes","No","Go back"),preselect="No",
                                     title ="Add a new region?",graphics=TRUE)
                          if(length(addbinRegion)==0|addbinRegion==""){
                            stopTMP<-TRUE
                            addbinRegion<-""
                          }
                          if(!stopTMP){
                            if(addbinRegion=="Yes"){
                              mrbin.env$mrbin$parameters$sumBinList<-rbind(mrbin.env$mrbin$parameters$sumBinList,c(0,0,0,0))
                            }
                          }
                        }
                        if(!stopTMP&!adjbinRegion=="Go back"&!addbinRegion=="Go back"&!addbinRegion=="No"){
                          mean1<-mean(mrbin.env$mrbin$parameters$sumBinList[ibinRegions,1:2])
                          range1<-mrbin.env$mrbin$parameters$sumBinList[ibinRegions,1]-mrbin.env$mrbin$parameters$sumBinList[ibinRegions,2]
                          mean2<-mean(mrbin.env$mrbin$parameters$sumBinList[ibinRegions,3:4])
                          range2<-mrbin.env$mrbin$parameters$sumBinList[ibinRegions,4]-mrbin.env$mrbin$parameters$sumBinList[ibinRegions,3]
                          regionTMP<-c(mean1+4*range1,mean1-4*range1,mean2-4*range2,mean2+4*range2)
                          if(sum(mrbin.env$mrbin$parameters$sumBinList[ibinRegions,]==0)==4) regionTMP<-"all"
                          if(mrbin.env$mrbin$parameters$dimension=="1D"){
                            if(mrbin.env$mrbin$parameters$showSpectrumPreview=="Yes") plotMultiNMR(
                                    region=regionTMP,rectangleRegions=matrix(c(
                                    mrbin.env$mrbin$parameters$sumBinList[ibinRegions,1],
                                    mrbin.env$mrbin$parameters$sumBinList[ibinRegions,2],0,2),ncol=4),
                                    color="black",manualScale=FALSE,maxPlots=2,
                                    plotTitle=paste("Sum area\nleft=",mrbin.env$mrbin$parameters$sumBinList[ibinRegions,1],
                                    "ppm, right=",mrbin.env$mrbin$parameters$sumBinList[ibinRegions,2],"ppm",sep=""),
                                    restrictToRange=TRUE)
                          }
                          if(mrbin.env$mrbin$parameters$dimension=="2D"){
                            if(mrbin.env$mrbin$parameters$showSpectrumPreview=="Yes") plotMultiNMR(
                                    region=regionTMP,
                                    rectangleRegions=matrix(mrbin.env$mrbin$parameters$sumBinList[ibinRegions,1:4],ncol=4),
                                    color="black", manualScale=FALSE,maxPlots=2,
                                    plotTitle=paste("Sum area\nleft=",mrbin.env$mrbin$parameters$sumBinList[ibinRegions,1],
                                    "ppm, right=",mrbin.env$mrbin$parameters$sumBinList[ibinRegions,2],"ppm",sep=""),
                                    restrictToRange=TRUE)
                          }
                          adjbinRegionAccept<-paste(paste(c("Keep left=","ppm, right=","ppm, top=","ppm, bottom=")[1:dimlength],
                                            mrbin.env$mrbin$parameters$sumBinList[ibinRegions,1:dimlength],collapse="",sep=""),"ppm",sep="")
                          if(sum(mrbin.env$mrbin$parameters$sumBinList[ibinRegions,]==0)==4){
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
                                  mrbin.env$mrbin$parameters$sumBinList[ibinRegions,1],": ",sep=""))
                        if(!regionTMP=="") {
                               mrbin.env$mrbin$parameters$sumBinList[ibinRegions,1]<-as.numeric(regionTMP)
                        }
                        regionTMP<-readline(prompt=paste("Region ",ibinRegions,": right border, press enter to keep ",
                                  mrbin.env$mrbin$parameters$sumBinList[ibinRegions,2],": ",sep=""))
                        if(!regionTMP=="") {
                               mrbin.env$mrbin$parameters$sumBinList[ibinRegions,2]<-as.numeric(regionTMP)
                        }
                        if(mrbin.env$mrbin$parameters$sumBinList[ibinRegions,1]<mrbin.env$mrbin$parameters$sumBinList[ibinRegions,2]){
                          TMP<-mrbin.env$mrbin$parameters$sumBinList[ibinRegions,1]
                          mrbin.env$mrbin$parameters$sumBinList[ibinRegions,1]<-mrbin.env$mrbin$parameters$sumBinList[ibinRegions,2]
                          mrbin.env$mrbin$parameters$sumBinList[ibinRegions,2]<-TMP
                        }
                      if(mrbin.env$mrbin$parameters$dimension=="2D"&!stopTMP){
                        regionTMP<-readline(prompt=paste("Region ",ibinRegions,": top border, press enter to keep ",
                                  mrbin.env$mrbin$parameters$sumBinList[ibinRegions,3],": ",sep=""))
                        if(!regionTMP=="") {
                               mrbin.env$mrbin$parameters$sumBinList[ibinRegions,3]<-as.numeric(regionTMP)
                        }
                        regionTMP<-readline(prompt=paste("Region ",ibinRegions,": bottom border, press enter to keep ",
                                  mrbin.env$mrbin$parameters$sumBinList[ibinRegions,4],": ",sep=""))
                        if(!regionTMP=="") {
                               mrbin.env$mrbin$parameters$sumBinList[ibinRegions,4]<-as.numeric(regionTMP)
                        }
                      if(mrbin.env$mrbin$parameters$sumBinList[ibinRegions,4]<mrbin.env$mrbin$parameters$sumBinList[ibinRegions,3]){
                        TMP<-mrbin.env$mrbin$parameters$sumBinList[ibinRegions,3]
                        mrbin.env$mrbin$parameters$sumBinList[ibinRegions,3]<-mrbin.env$mrbin$parameters$sumBinList[ibinRegions,4]
                        mrbin.env$mrbin$parameters$sumBinList[ibinRegions,4]<-TMP
                      }
                     }
                    }
                    if(adjbinRegion=="Remove"&!stopTMP&!adjbinRegion=="Go back"){
                      mrbin.env$mrbin$parameters$sumBinList<-mrbin.env$mrbin$parameters$sumBinList[-ibinRegions,,drop=FALSE]
                    }
                    if(adjbinRegion==adjbinRegionAccept&!stopTMP){
                      ibinRegions <- ibinRegions+1
                    }
                  }
                }
              }
            if(is.null(mrbin.env$mrbin$parameters$sumBinList)){
              adjbinRegion<-"Go back"
            } else {
              if(nrow(mrbin.env$mrbin$parameters$sumBinList)==0){
                mrbin.env$mrbin$parameters$sumBinList<-NULL
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
        if(selectStep==10){#Crop HSQCs
          if(mrbin.env$mrbin$parameters$dimension=="2D"&!stopTMP){
             if(!stopTMP&mrbin.env$mrbin$parameters$verbose){
               message("Hint: Cropping may remove noisy areas, optimized for HSQCs")
               utils::flush.console()
             }
             if(mrbin.env$mrbin$parameters$showSpectrumPreview=="Yes") plotMultiNMR(region="all",
                    polygonRegion=matrix(c(mrbin.env$mrbin$parameters$croptopRight,
                                  mrbin.env$mrbin$parameters$croptopLeft,
                                  mrbin.env$mrbin$parameters$cropbottomLeft,
                                  mrbin.env$mrbin$parameters$cropbottomRight),
                                  ncol=2,byrow=TRUE),
                    color="black",manualScale=FALSE,maxPlots=2,
                    plotTitle=paste("Crop spectrum to diagonal",sep=""))
             cropHSQC<-utils::select.list(c("Yes","No","Go back"),
                                   preselect=mrbin.env$mrbin$parameters$cropHSQC,
                                   title="Crop spectra?",graphics=TRUE)
             if(length(cropHSQC)==0|cropHSQC=="") stopTMP<-TRUE
             if(!stopTMP&!cropHSQC=="Go back"){
               mrbin.env$mrbin$parameters$cropHSQC<-cropHSQC
             }
            if(cropHSQC=="Go back"){
               selectStep<-selectStep-2
            }
          }
          if(!stopTMP) selectStep<-selectStep+1
        }
        if(selectStep==11){#Define sample names
          if(!stopTMP){
            if(!stopTMP&mrbin.env$mrbin$parameters$verbose){
              message("Hint: If only EXPNO differs choose Folder names and EXPNO")
              utils::flush.console()
            }
            NamesDictTMP<-c("Folder names","Spectrum titles","Folder names and EXPNO")
            names(NamesDictTMP)<-paste(NamesDictTMP," (\"",c(mrbin.env$mrbinTMP$currentSpectrumFolderName,
                                 mrbin.env$mrbinTMP$currentSpectrumTitle,
                                 mrbin.env$mrbinTMP$currentSpectrumFolderName_EXPNO),
                                 "\", ...)",sep="")
            NamesDictTMP2<-names(NamesDictTMP)
            names(NamesDictTMP2)<-NamesDictTMP
            useAsNames<-utils::select.list(c(names(NamesDictTMP),"Go back"),
                                      preselect=NamesDictTMP2[mrbin.env$mrbin$parameters$useAsNames],
                                      title = "Create sample names from",graphics=TRUE)
            if(length(useAsNames)==0|useAsNames=="") stopTMP<-TRUE
            if(!stopTMP&!useAsNames=="Go back"){
              mrbin.env$mrbin$parameters$useAsNames<-NamesDictTMP[useAsNames]
            }
            if(!stopTMP&useAsNames=="Go back"){
               selectStep<-selectStep-2
            }
            if(!stopTMP) selectStep<-selectStep+1
          }
        }
        if(selectStep==12){#Plot results
          if(!stopTMP){
           PCA<-"Yes"
           if(mrbin.env$mrbin$parameters$PCA=="No"){
            if(!stopTMP&mrbin.env$mrbin$parameters$verbose){
              message("Hint: Recommended for quality control")
              utils::flush.console()
            }
            PCAtitlelength<-""
            PCA<-utils::select.list(c("Yes","No","Go back"),
                              preselect = "Yes",#mrbin.env$mrbin$parameters$PCA,
                              title = "Create result plot?",graphics=TRUE)
            if(length(PCA)==0|PCA=="") stopTMP<-TRUE
            if(!stopTMP&!PCA=="Go back"){
              mrbin.env$mrbin$parameters$PCA<-PCA
            }
           }
           if(!stopTMP&!PCA=="Go back"){
              if(!stopTMP&mrbin.env$mrbin$parameters$PCA=="Yes"){
                currentPCAtitlelength<-as.character(mrbin.env$mrbin$parameters$PCAtitlelength)
                if(mrbin.env$mrbin$parameters$useAsNames=="Spectrum titles") Title<-mrbin.env$mrbinTMP$currentSpectrumTitle
                if(mrbin.env$mrbin$parameters$useAsNames=="Folder names") Title<-mrbin.env$mrbinTMP$currentSpectrumFolderName
                if(mrbin.env$mrbin$parameters$useAsNames=="Folder names and EXPNO") Title<-mrbin.env$mrbinTMP$currentSpectrumFolderName_EXPNO
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
                if(!stopTMP&mrbin.env$mrbin$parameters$verbose){
                  message("Hint: Recommended for a nicer plot, make sure names are unique")
                  utils::flush.console()
                }
                PCAtitlelength<-utils::select.list(names(TitleListTMP2),
                                  preselect = names(TitleListTMP2)[TitleListTMP4==as.character(mrbin.env$mrbin$parameters$PCAtitlelength)],
                                  title = "Crop titles for plot?",graphics=TRUE)
                if(length(PCAtitlelength)==0|PCAtitlelength=="") stopTMP<-TRUE
                if(!stopTMP&!PCAtitlelength=="Go back"){
                  if(TitleListTMP2[PCAtitlelength]=="All"){
                    mrbin.env$mrbin$parameters$PCAtitlelength<-500
                  } else {
                    if(PCAtitlelength=="Custom..."){
                        PCAtitlelengthTMP<-readline(prompt=paste("New title length, press enter to keep ",
                            mrbin.env$mrbin$parameters$PCAtitlelength,": ",sep=""))
                        if(!PCAtitlelengthTMP=="") {
                            mrbin.env$mrbin$parameters$PCAtitlelength<-as.numeric(PCAtitlelengthTMP)
                        }
                    } else {
                       mrbin.env$mrbin$parameters$PCAtitlelength<-as.numeric(TitleListTMP2[PCAtitlelength])
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
       if(selectStep==13){#Save output files to hard drive?
         if(!stopTMP){
           saveFilesTMP2<-"Select new folder and file name"
           saveFilesTMP<-utils::select.list(c("Yes","No","Go back"),
                       preselect=mrbin.env$mrbin$parameters$saveFiles,
                       title ="Save output to disk?",graphics=TRUE)
           if(length(saveFilesTMP)==0|saveFilesTMP=="") stopTMP<-TRUE
           if(!stopTMP&!saveFilesTMP=="Go back"){
            mrbin.env$mrbin$parameters$saveFiles<-saveFilesTMP
            if(mrbin.env$mrbin$parameters$saveFiles=="Yes"&!stopTMP){
              if(!is.null(mrbin.env$mrbin$parameters$outputFileName)){
               keepFileTMP<-paste("Keep ",mrbin.env$mrbin$parameters$outputFileName,sep="")
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
                   filenameTMP<-utils::select.list(c(paste("mrbin_",gsub(":","-",gsub(" ","_",Sys.Date())),
                               sep=""),"Change..."),
                               title ="Output file name: ",graphics=TRUE)
                   if(length(filenameTMP)==0|filenameTMP=="") stopTMP<-TRUE
                   if(!stopTMP){
                    if(filenameTMP=="Change..."&!stopTMP){
                      filenameTMP<-readline(prompt=paste("New file name, press enter to use ",
                                paste("mrbin_",gsub(":","-",gsub(" ","_",Sys.Date())),sep=""),": \n",sep=""))
                      if(filenameTMP=="") filenameTMP<-paste("mrbin_",gsub(":","-",gsub(" ","_",Sys.Date())),
                                 sep="")
                     }
                     mrbin.env$mrbin$parameters$outputFileName<-gsub("//","/",paste(
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
       if(selectStep==14){
         if(!stopTMP){#make some example calculations to estimate speed of binning
          createBinNumbers()#necessary here for time estimate
          createBinRegions()#necessary here for time estimate
          if(mrbin.env$mrbin$parameters$dimension=="1D") coverageRatioTMP<-nrow(
            mrbin.env$mrbin$parameters$binRegions)/((max(as.numeric(
            names(mrbin.env$mrbinTMP$currentSpectrum)))-min(as.numeric(
            names(mrbin.env$mrbinTMP$currentSpectrum))))/mrbin.env$mrbin$parameters$binwidth1D)
          if(mrbin.env$mrbin$parameters$dimension=="2D") coverageRatioTMP<-nrow(
            mrbin.env$mrbin$parameters$binRegions)/((max(as.numeric(
            colnames(mrbin.env$mrbinTMP$currentSpectrum)))-min(as.numeric(
            colnames(mrbin.env$mrbinTMP$currentSpectrum))))*
            (max(as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrum)))-min(as.numeric(
            rownames(mrbin.env$mrbinTMP$currentSpectrum))))/
            (mrbin.env$mrbin$parameters$binwidth2D*mrbin.env$mrbin$parameters$binheight))
          dimScaleTMP<-1
          if(mrbin.env$mrbin$parameters$dimension=="2D") dimScaleTMP<-.5#2D spectra take less time empirically
          NrowTMP<-400#mock number of bins
          NrowTMP2<-nrow(mrbin.env$mrbin$parameters$binRegions)/NrowTMP#ratio of real number of bins to mock number of bins
          NpointsTMP<-100#mock number of data points per bin
          NpointsTMP2<-ceiling(1+(coverageRatioTMP*length(mrbin.env$mrbinTMP$currentSpectrum)/
                                    nrow(mrbin.env$mrbin$parameters$binRegions)))/NpointsTMP#ratio of estimated number of data points per bin to mock data point number
          #calculate some sums in a loop as a mock binning example for time estimation
          mrbin.env$mrbinTMP$timeEstimate<-mrbin.env$mrbinTMP$timeEstimate+max(.001,system.time(
            for(i in 1:NrowTMP) sum((1:(NpointsTMP*NrowTMP))[(1:(NpointsTMP*NrowTMP))<(i*NpointsTMP)&
                                                               (1:(NpointsTMP*NrowTMP))>((i-1)*NpointsTMP)]))[1])*
            (dimScaleTMP*20*NrowTMP2*NpointsTMP2^.5)
          if(mrbin.env$mrbin$parameters$tryParallel){
           mrbin.env$mrbinTMP$timeEstimate<-dimScaleTMP*mrbin.env$mrbinTMP$timeEstimate*ceiling(length(
             mrbin.env$mrbin$parameters$NMRfolders)/max(1,(parallel::detectCores()-1)))
          } else {
           mrbin.env$mrbinTMP$timeEstimate<-dimScaleTMP*mrbin.env$mrbinTMP$timeEstimate*length(
             mrbin.env$mrbin$parameters$NMRfolders)
          }
           if(mrbin.env$mrbin$parameters$verbose){
              message(paste("Hint: Estimated binning time: ",round(
              mrbin.env$mrbinTMP$timeEstimate/60,0)," minutes or more",sep=""))
              utils::flush.console()
            }
           startmrbin<-utils::select.list(c("Start binning now","I'll do it later","Go back"),
                      preselect="Start binning now",
                      title =paste("Start now? Estimate: >",round(
                      mrbin.env$mrbinTMP$timeEstimate/60,0)," min",sep=""),graphics=TRUE)
           if(length(startmrbin)==0|startmrbin==""|startmrbin=="I'll do it later") stopTMP<-TRUE
           if(!stopTMP&startmrbin=="Go back"){
             selectStep<-selectStep-2
           }
           if(!stopTMP) selectStep<-selectStep+1
           if(startmrbin=="Start binning now") lastStepDone<-TRUE
         }
       }
     }
    }
   }#end if(!silent)
   if(!stopTMP){
    if(startmrbin=="Start binning now"){
     #Check if files or folders exist first to avoid long waiting due to binning failure
     for(iCheckFiles in 1:length(mrbin.env$mrbin$parameters$NMRfolders)){
       readNMR(folder=mrbin.env$mrbin$parameters$NMRfolders[iCheckFiles],
         dimension=mrbin.env$mrbin$parameters$dimension,checkFiles=TRUE)
     }
     mrbin.env$mrbin<-mrbinrun(createbins=TRUE,process=FALSE,silent=silent)
     if(mrbin.env$mrbin$parameters$verbose){
      if(!is.null(mrbin.env$mrbin$parameters$warningMessages)){
       for(iWarningTMP in 1:length(mrbin.env$mrbin$parameters$warningMessages)){
          message("Warning: ",mrbin.env$mrbin$parameters$warningMessages[iWarningTMP])
       }
       utils::flush.console()
      }
     }
    }
   }
   if(!silent){
    lastStepDone<-FALSE
    if(selectionRepeat=="Use current parameters without changes"&!stopTMP){
       selectStep<-25
       lastStepDone<-TRUE
    }
    while(!lastStepDone&!stopTMP&!restart){
     if(selectionRepeat=="Review parameters"&!stopTMP){
           if(selectStep==15){
              if(!stopTMP){
                if(mrbin.env$mrbin$parameters$verbose){
                 message("Hint: Review plots for data quality. If issues are present, such as phasing\nor baseline issues, fix the spectra, e.g. in Topspin, then run mrbin again")
                 utils::flush.console()
                }
                plotReview<-utils::select.list(c("I will fix issues, if any, and then run mrbin again",
                                                 "I have fixed all spectrum issues and wish to continue",
                                                 "I would like to restart to adjust parameters"),
                                         preselect="I will fix issues, if any, and then run mrbin again",
                                         title="Please review data for quality issues",graphics=TRUE)
                if(length(plotReview)==0|plotReview==""|plotReview=="I will fix issues, if any, and then run mrbin again") stopTMP<-TRUE
                if(!stopTMP&plotReview=="I would like to restart to adjust parameters"){
                  restart<-TRUE
                  selectStep<--4#start from beginning
                }
                if(!stopTMP) selectStep<-selectStep+1
              }
           }
           if(selectStep==16){#data quality plots
              if(!stopTMP){
               if(!is.null(mrbin.env$mrbin$parameters$warningMessages)){
                 if(mrbin.env$mrbin$parameters$verbose){
                   message("Hint: Run warnings(), fix issues if possible, then run mrbin again\n")
                   utils::flush.console()
                 }
                 listWarningTMP<-c("I will check warnings() and then run mrbin again",
                                   "I have fixed all issues and wish to continue",
                                   "I would like to restart to adjust parameters")
                 plotReview<-utils::select.list(listWarningTMP,
                                           preselect="I will check warnings() and then run mrbin again",
                                           title="There were warning messages",graphics=TRUE)
                 if(length(plotReview)==0|plotReview==""|plotReview=="I will check warnings() and then run mrbin again") stopTMP<-TRUE
                }
                if(!stopTMP&plotReview=="I would like to restart to adjust parameters"){
                  selectStep<--4#start from beginning
                  restart<-TRUE
                }
                if(!stopTMP) selectStep<-selectStep+1
              }
           }
           if(selectStep==17){#Remove noise
              if(!stopTMP){
                if(mrbin.env$mrbin$parameters$verbose){
                  message("Hint: All processing steps can be performed later if you skip them now")
                  message("Hint: Remove noise to increase statistical power")
                  utils::flush.console()
                }
                noiseRemoval<-utils::select.list(c("Yes","No"),
                                         preselect=mrbin.env$mrbin$parameters$noiseRemoval,
                                         title="Remove noise?",graphics=TRUE)
                if(length(noiseRemoval)==0|noiseRemoval=="") stopTMP<-TRUE
                if(!stopTMP){
                  if(!noiseRemoval==mrbin.env$mrbin$parameters$noiseRemoval){
                    mrbin.env$mrbin<-editmrbin(mrbinObject=mrbin.env$mrbin,functionName="mrbin::mrbin",
                      versionNumber=as.character(utils::packageVersion("mrbin")),
                      parameters=list(noiseRemoval=noiseRemoval),verbose=FALSE)
                  }
                  if(mrbin.env$mrbin$parameters$noiseRemoval=="Yes"){
                    mrbin.env$mrbin<-setNoiseLevels(mrbin.env$mrbin,plotOnly=FALSE)
                  }
                  selectStep<-selectStep+1
                }
              }
            }
            if(selectStep==18){#Dilution correction scaling
              if(!stopTMP){
                if(!stopTMP&mrbin.env$mrbin$parameters$verbose){
                  message("Hint: Use if volumes or weights differ between samples")
                  utils::flush.console()
                }
                dilutionCorrection<-utils::select.list(c("Yes","No","Go back"),
                       preselect=mrbin.env$mrbin$parameters$dilutionCorrection,
                       title = "Dilution correction?",graphics=TRUE)
                if(length(dilutionCorrection)==0|dilutionCorrection=="") stopTMP<-TRUE
                if(!stopTMP&!dilutionCorrection=="Go back"){
                  if(!dilutionCorrection==mrbin.env$mrbin$parameters$dilutionCorrection){
                   mrbin.env$mrbin<-editmrbin(mrbinObject=mrbin.env$mrbin,
                     functionName="mrbin::mrbin",
                     versionNumber=as.character(utils::packageVersion("mrbin")),
                     parameters=list(dilutionCorrection=dilutionCorrection),verbose=FALSE)
                  }
                  if(mrbin.env$mrbin$parameters$dilutionCorrection=="Yes"){
                    mrbin.env$mrbin<-setDilutionFactors(mrbin.env$mrbin)
                  }
                }
                if(!stopTMP&(dilutionCorrection=="Go back")){
                   selectStep<-selectStep-2
                }
                if(!stopTMP) selectStep<-selectStep+1
              }
            }
            if(selectStep==19){#PQN scaling
              if(!stopTMP){
                PQNScalingIgnoreSugar<-""
                if(!stopTMP&mrbin.env$mrbin$parameters$verbose){
                  message("Hint: Recommended for urine and tissue extracts")
                  utils::flush.console()
                }
                PQNScaling<-utils::select.list(c("Yes","No","Go back"),
                                        preselect=mrbin.env$mrbin$parameters$PQNScaling,
                                        title = "PQN normalization?",graphics=TRUE)
                if(length(PQNScaling)==0|PQNScaling=="") stopTMP<-TRUE
                if(!stopTMP&!PQNScaling=="Go back"){
                  if(!PQNScaling==mrbin.env$mrbin$parameters$PQNScaling){
                   mrbin.env$mrbin<-editmrbin(mrbinObject=mrbin.env$mrbin,functionName="mrbin::mrbin",
                     versionNumber=as.character(utils::packageVersion("mrbin")),
                     parameters=list(PQNScaling=PQNScaling),verbose=FALSE)
                  }
                  if(mrbin.env$mrbin$parameters$PQNScaling=="Yes"){
                    if(mrbin.env$mrbin$parameters$verbose){
                      message("Hint: Improves PQN but works only for 1H and 1H-13C spectra")
                      utils::flush.console()
                    }
                    PQNScalingIgnoreSugar<-utils::select.list(c("Yes","No","Go back"),
                                            preselect=mrbin.env$mrbin$parameters$PQNIgnoreSugarArea,
                                            title = "Ignore glucose for PQN?",graphics=TRUE)
                    if(length(PQNScalingIgnoreSugar)==0|PQNScalingIgnoreSugar=="") stopTMP<-TRUE
                    if(!stopTMP&!PQNScalingIgnoreSugar=="Go back"&
                     !PQNScalingIgnoreSugar==mrbin.env$mrbin$parameters$PQNIgnoreSugarArea){
                     mrbin.env$mrbin<-editmrbin(mrbinObject=mrbin.env$mrbin,
                       functionName="mrbin::mrbin",
                       versionNumber=as.character(utils::packageVersion("mrbin")),
                       parameters=list(PQNIgnoreSugarArea=PQNScalingIgnoreSugar),
                       verbose=FALSE)
                    }
                  }
                }
                if(!stopTMP&(PQNScaling=="Go back"|PQNScalingIgnoreSugar=="Go back")){
                   selectStep<-selectStep-2
                }
                if(!stopTMP) selectStep<-selectStep+1
              }
            }
            if(selectStep==20){#Replace negative values
              if(!stopTMP&mrbin.env$mrbin$parameters$verbose){
                message("Hint: Replace negative values if you plan to do log transform")
                utils::flush.console()
              }
              if(!stopTMP){
                fixNegatives<-utils::select.list(c("Yes","No","Go back"),
                                          preselect=mrbin.env$mrbin$parameters$fixNegatives,
                                          title="Fix negative values (atnv)",graphics=TRUE)
                if(length(fixNegatives)==0|fixNegatives=="") stopTMP<-TRUE
                if(!stopTMP&!fixNegatives=="Go back"&!fixNegatives==mrbin.env$mrbin$parameters$fixNegatives){
                 mrbin.env$mrbin<-editmrbin(mrbinObject=mrbin.env$mrbin,
                   functionName="mrbin::mrbin",
                   versionNumber=as.character(utils::packageVersion("mrbin")),
                   parameters=list(fixNegatives=fixNegatives),
                   verbose=FALSE)
                }
                if(!stopTMP&fixNegatives=="Go back"){
                   selectStep<-selectStep-2
                }
                if(!stopTMP) selectStep<-selectStep+1
              }
            }
            if(selectStep==21){#Log scaling
              if(!stopTMP){
               #if(mrbin.env$mrbin$parameters$fixNegatives=="Yes"|mrbin.env$mrbin$parameters$logTrafo=="Yes"){
                if(!stopTMP&mrbin.env$mrbin$parameters$verbose){
                  message("Hint: Makes data more normal but breaks linearity. Requires positive data.")
                  utils::flush.console()
                }
                preselectTMP<-mrbin.env$mrbin$parameters$logTrafo
                #if(!mrbin.env$mrbin$parameters$fixNegatives=="Yes"){
                #  preselectTMP<-"No"
                #}
                logTrafo<-utils::select.list(c("Yes","No","Go back"),preselect=preselectTMP,
                                      title="Log transformation?",graphics=TRUE)
                if(length(logTrafo)==0|logTrafo=="") stopTMP<-TRUE
                if(!stopTMP&!logTrafo=="Go back"&!logTrafo==mrbin.env$mrbin$parameters$logTrafo){
                     mrbin.env$mrbin<-editmrbin(mrbinObject=mrbin.env$mrbin,
                       functionName="mrbin::mrbin",
                       versionNumber=as.character(utils::packageVersion("mrbin")),
                       parameters=list(logTrafo=logTrafo),
                       verbose=FALSE)
                }
                if(!stopTMP&logTrafo=="Go back"){
                   selectStep<-selectStep-2
                }
               #}
               if(!stopTMP) selectStep<-selectStep+1
              }
            }
            if(selectStep==22){#Unit variance scaling
              if(!stopTMP){
                if(!stopTMP&mrbin.env$mrbin$parameters$verbose){
                  message("Hint: Usually not required. Breaks linearity.")
                  utils::flush.console()
                }
                preselectTMP<-mrbin.env$mrbin$parameters$unitVarianceScaling
                unitVarianceScaling<-utils::select.list(c("Yes","No","Go back"),preselect=preselectTMP,
                                      title="Unit variance scaling?",graphics=TRUE)
                if(length(unitVarianceScaling)==0|unitVarianceScaling=="") stopTMP<-TRUE
                if(!stopTMP&!unitVarianceScaling=="Go back"&!unitVarianceScaling==mrbin.env$mrbin$parameters$unitVarianceScaling){
                     mrbin.env$mrbin<-editmrbin(mrbinObject=mrbin.env$mrbin,
                       functionName="mrbin::mrbin",
                       versionNumber=as.character(utils::packageVersion("mrbin")),
                       parameters=list(unitVarianceScaling=unitVarianceScaling),
                       verbose=FALSE)
                }
                if(!stopTMP&unitVarianceScaling=="Go back"){
                   selectStep<-selectStep-2
                }
               #}
               if(!stopTMP) selectStep<-selectStep+1
              }
            }
       if(selectStep>=23){
         if(!stopTMP&!restart)  lastStepDone<-TRUE
       }
     }
    }
   }#end if(!silent)
 }
 if(!stopTMP){
   mrbin.env$mrbin<-mrbinrun(createbins=FALSE,process=TRUE,mrbinResults=mrbin.env$mrbin,
     silent=silent)
   if(mrbin.env$mrbin$parameters$verbose){
       mrbin.env$mrbin<-editmrbin(mrbinObject=mrbin.env$mrbin,
         functionName="mrbin::mrbin",
         versionNumber=as.character(utils::packageVersion("mrbin")),
         parameters=list(createCode=printParameters()),verbose=FALSE)
   }
   invisible(mrbin.env$mrbin)
 }
}

#' A function performing all data read and processing steps.
#'
#' This function reads parameters from the global variable mrbin.env$mrbin$parameters and
#' performs the following operations:
#' Reading NMR files, creating bins, removing solvent area, removing additional
#' user-defined areas, summing up bins that contain unstable peaks such as
#' citric acid, removes noise bins, crops HSQC spectra to the diagonal area,
#' performs PQN scaling, replaces negative values, log transforms and displays a
#' PCA plot. Parameters are then saved in a text file. These can be recreated
#' using recreatemrbin().
#' @param createbins If TRUE, new bin data is generated
#' @param process If TRUE, bin data is processed, e.g. by noise removal, atnv, etc.
#' @param mrbinResults An mrbin object. Needs to be provided only if createbins is FALSE
#' @param silent If set to FALSE, no new time calculation is performed
#' @return An invisible mrbin object
#' @export
#' @examples
#' resetEnv()
#' setParam(parameters=list(dimension="2D",binwidth2D=0.1,binheight=5,
#'    binRegion=c(8,1,15,140),PQNScaling="No",tryParallel=FALSE,
#'    fixNegatives="No",logTrafo="No",signal_to_noise2D=10,solventRegion=c(5.5,4.2),
#'    NMRfolders=c(system.file("extdata/1/12/pdata/10",package="mrbin"),
#'                 system.file("extdata/2/12/pdata/10",package="mrbin"))))
#' mrbinrun()

mrbinrun<-function(createbins=TRUE,process=TRUE,mrbinResults=NULL,silent=TRUE){
  defineGroups<-FALSE
  if(!is.null(mrbin.env$mrbin$parameters$NMRfolders)){
    if(createbins){
      mrbin.env$mrbinTMP$scaleFactorTMP1<-NULL
      mrbin.env$mrbinTMP$scaleFactorTMP2<-NULL
      mrbin.env$mrbinTMP$scaleFactorTMP3<-NULL
      if(mrbin.env$mrbin$parameters$verbose){
        message("\nPreparing parameters... ", appendLF = FALSE)
        utils::flush.console()
      }
      mrbin.env$mrbinTMP$currentFolder<-mrbin.env$mrbin$parameters$NMRfolders[1]
      readNMR2()
      createBinNumbers()
      createBinRegions()
      mrbin.env$mrbin$parameters$numberOfFeaturesRaw<-nrow(mrbin.env$mrbin$parameters$binRegions)
      if(mrbin.env$mrbin$parameters$removeSolvent=="Yes") removeSolvent()
      if(mrbin.env$mrbin$parameters$removeAreas=="Yes") removeAreas()
      if(mrbin.env$mrbin$parameters$sumBins=="Yes") sumBins()
      if(mrbin.env$mrbin$parameters$cropHSQC=="Yes"&mrbin.env$mrbin$parameters$dimension=="2D") cropNMR()
      #Sort bins
      if(nrow(mrbin.env$mrbin$parameters$binRegions)>1){
        binRegionsTMP<-mrbin.env$mrbin$parameters$binRegions
        binRegionsTMP<-binRegionsTMP[rev(order(binRegionsTMP[,3])),,drop=FALSE]
      }
      if(mrbin.env$mrbin$parameters$verbose){
        message("done. \n", appendLF = FALSE)
        utils::flush.console()
      }
      if(mrbin.env$mrbin$parameters$verbose){
        message("Binning spectra... ", appendLF = FALSE)
        utils::flush.console()
      }
      if(silent|mrbin.env$mrbinTMP$timeEstimate==0){#time needs to be estimated now if mrbin was not run in interactive mode
          mrbin.env$mrbinTMP$currentFolder<-mrbin.env$mrbin$parameters$NMRfolders[1]
          readNMR2()#necessary to fill the variables, otherwise the next line is "too fast" and estimate will be off
          mrbin.env$mrbinTMP$timeEstimate<-max(.001,system.time(readNMR2())[1])
          if(mrbin.env$mrbin$parameters$dimension=="1D") coverageRatioTMP<-nrow(
            mrbin.env$mrbin$parameters$binRegions)/((max(as.numeric(
            names(mrbin.env$mrbinTMP$currentSpectrum)))-min(as.numeric(
            names(mrbin.env$mrbinTMP$currentSpectrum))))/mrbin.env$mrbin$parameters$binwidth1D)
          if(mrbin.env$mrbin$parameters$dimension=="2D") coverageRatioTMP<-nrow(
            mrbin.env$mrbin$parameters$binRegions)/((max(as.numeric(
            colnames(mrbin.env$mrbinTMP$currentSpectrum)))-min(as.numeric(
            colnames(mrbin.env$mrbinTMP$currentSpectrum))))*
            (max(as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrum)))-min(as.numeric(
            rownames(mrbin.env$mrbinTMP$currentSpectrum))))/
            (mrbin.env$mrbin$parameters$binwidth2D*mrbin.env$mrbin$parameters$binheight))
          dimScaleTMP<-1
          if(mrbin.env$mrbin$parameters$dimension=="2D") dimScaleTMP<-.5
          NrowTMP<-400#mock number of bins
          NrowTMP2<-nrow(mrbin.env$mrbin$parameters$binRegions)/NrowTMP#ratio of real number of bins to mock number of bins
          NpointsTMP<-100#mock number of data points per bin
          NpointsTMP2<-ceiling(1+(coverageRatioTMP*length(mrbin.env$mrbinTMP$currentSpectrum)/
                                    nrow(mrbin.env$mrbin$parameters$binRegions)))/NpointsTMP#ratio of estimated number of data points per bin to mock data point number
          #calculate some sums in a loop as a mock binning example for time estimation
          mrbin.env$mrbinTMP$timeEstimate<-mrbin.env$mrbinTMP$timeEstimate+max(.001,system.time(
            for(i in 1:NrowTMP) sum((1:(NpointsTMP*NrowTMP))[(1:(NpointsTMP*NrowTMP))<(i*NpointsTMP)&
                                                               (1:(NpointsTMP*NrowTMP))>((i-1)*NpointsTMP)]))[1])*
            (dimScaleTMP*20*NrowTMP2*NpointsTMP2^.5)
          if(mrbin.env$mrbin$parameters$tryParallel){
            mrbin.env$mrbinTMP$timeEstimate<-dimScaleTMP*mrbin.env$mrbinTMP$timeEstimate*ceiling(length(
              mrbin.env$mrbin$parameters$NMRfolders)/max(1,(parallel::detectCores()-1)))
          } else {
            mrbin.env$mrbinTMP$timeEstimate<-dimScaleTMP*mrbin.env$mrbinTMP$timeEstimate*length(
              mrbin.env$mrbin$parameters$NMRfolders)
          }
      }
      if(mrbin.env$mrbin$parameters$verbose){
        message(paste("(estimated time: ",round(
        mrbin.env$mrbinTMP$timeEstimate/60,0)," min or more) ",sep=""), appendLF = FALSE)
        utils::flush.console()
      }
      #This creates the actual bins:
      mrbinResults<-binMultiNMR()
      if(mrbinResults$parameters$trimZeros=="Yes") mrbinResults<-trimZeros(mrbinResults)
      if(mrbinResults$parameters$verbose) message("done.\n", appendLF = FALSE)
      utils::flush.console()
    }
    if(process){
      if(mrbinResults$parameters$verbose){
        message("Processing data... ", appendLF = FALSE)
        utils::flush.console()
      }
        mrbin.env$mrbinTMP$additionalPlots1D<-NULL
        mrbin.env$mrbinTMP$additionalPlots1DMetadata<-NULL
        mrbin.env$mrbinTMP$additionalPlots2D<-NULL
        mrbin.env$mrbinTMP$additionalPlots2DMetadata<-NULL
        #Find 3 more spectra: 33 percentile, 66 percentile, last spectrum
        mrbin.env$mrbinTMP$spectrumListPlotTMP<-setdiff(unique(c(
          ceiling(length(mrbinResults$parameters$NMRfolders)*.33),
          ceiling(length(mrbinResults$parameters$NMRfolders)*.66),
          length(mrbinResults$parameters$NMRfolders))),1)
        if(length(mrbin.env$mrbinTMP$spectrumListPlotTMP)>0){
          #for 2D spectra, only plot one additional spectrum to increase speed
          for(ispectrumListPlotTMP in 1:length(mrbin.env$mrbinTMP$spectrumListPlotTMP)){
            addToPlot(folder=mrbinResults$parameters$NMRfolders[
              mrbin.env$mrbinTMP$spectrumListPlotTMP[ispectrumListPlotTMP]],
               dimension=mrbinResults$parameters$dimension,
               NMRvendor=mrbinResults$parameters$NMRvendor,
               useAsNames=mrbinResults$parameters$useAsNames)
          }
        }
      #create and save noise plots
      if(silent&mrbinResults$parameters$saveFiles=="Yes") setNoiseLevels(mrbinResults,plotOnly=TRUE,silent=silent)
      if(mrbinResults$parameters$noiseRemoval=="Yes") mrbinResults<-removeNoise(mrbinResults,verbose=FALSE)
      if(mrbinResults$parameters$dilutionCorrection=="Yes") mrbinResults<-dilutionCorrection(mrbinResults)
      if(mrbinResults$parameters$fixNegatives=="Yes") mrbinResults<-atnv(mrbinResults,verbose=FALSE)
      if(mrbinResults$parameters$PQNScaling=="Yes") mrbinResults<-PQNScaling(mrbinResults,verbose=FALSE)
      if(mrbinResults$parameters$logTrafo=="Yes") mrbinResults<-logTrafo(mrbinResults,verbose=FALSE)
      if(mrbinResults$parameters$unitVarianceScaling=="Yes") mrbinResults<-unitVarianceScaling(mrbinResults,verbose=FALSE)
      if(mrbinResults$parameters$verbose){
         message("done.\n", appendLF = FALSE)
         utils::flush.console()
      }
      if(mrbinResults$parameters$saveFiles=="Yes"){
        save(mrbinResults,file=paste(mrbinResults$parameters$outputFileName,".Rdata",sep=""))
        dput(mrbinResults$parameters, file = paste(mrbinResults$parameters$outputFileName,".txt",sep=""))
        #utils::write.csv(mrbinResults$bins,file=paste(mrbinResults$parameters$outputFileName,"bins.csv",sep=""))
      }
       resultOutputTMP<-c("\nNumber of spectra: ",nrow(mrbinResults$bins),"\n",
           "Number of bins at start: ",mrbinResults$parameters$numberOfFeaturesRaw,"\n")
       if(!is.null(mrbinResults$parameters$numberOfFeaturesAfterRemovingSolvent)&
         mrbinResults$parameters$removeSolvent=="Yes"){
            resultOutputTMP<-c(resultOutputTMP,"Number of bins after removing solvent: ",
               mrbinResults$parameters$numberOfFeaturesAfterRemovingSolvent,"\n")
       }
       if(!is.null(mrbinResults$parameters$numberOfFeaturesAfterRemovingAreas)&
         mrbinResults$parameters$removeAreas=="Yes"){
            resultOutputTMP<-c(resultOutputTMP,"Number of bins after removing areas: ",
            mrbinResults$parameters$numberOfFeaturesAfterRemovingAreas,"\n")
       }
       if(!is.null(mrbinResults$parameters$numberOfFeaturesAfterSummingBins)&mrbinResults$parameters$sumBins=="Yes"){
            resultOutputTMP<-c(resultOutputTMP,
               "Number of bins after summing bins: ",mrbinResults$parameters$numberOfFeaturesAfterSummingBins,"\n")
       }
       if(!is.null(mrbinResults$parameters$numberOfFeaturesAfterCropping)&
         mrbinResults$parameters$cropHSQC=="Yes"&mrbinResults$parameters$dimension=="2D"){
            resultOutputTMP<-c(resultOutputTMP,
               "Number of bins after cropping: ",mrbinResults$parameters$numberOfFeaturesAfterCropping,"\n")
       }
       if(!is.null(mrbinResults$parameters$numberOfFeaturesAfterTrimmingZeros)&
         mrbinResults$parameters$trimZeros=="Yes"){
            resultOutputTMP<-c(resultOutputTMP,"Number of bins after trimming zero-value bins: ",
               mrbinResults$parameters$numberOfFeaturesAfterTrimmingZeros,"\n")
       }
       if(!is.null(mrbinResults$parameters$numberOfFeaturesAfterNoiseRemoval)&
         mrbinResults$parameters$noiseRemoval=="Yes"){
            resultOutputTMP<-c(resultOutputTMP,"Number of bins after noise removal: ",
               mrbinResults$parameters$numberOfFeaturesAfterNoiseRemoval,"\n")
       }
       resultOutputTMP<-paste(resultOutputTMP,sep="")
       if(mrbinResults$parameters$verbose){
         message(resultOutputTMP, appendLF = FALSE)
         utils::flush.console()
       }
    }
    if(mrbinResults$parameters$PCA=="Yes"){
      if(mrbin.env$mrbin$parameters$saveFiles=="Yes"|process|!silent){
        plotResults(mrbinResults,defineGroups=defineGroups,process=process,silent=silent)
      }
    }
    invisible(mrbinResults)
  }
}

#' A function for setting dilution factors.
#'
#' This function edits the dilution factors of an mrbin object but does not change the bin data.
#' @param mrbinObject An mrbin object
#' @param dilutionFactors An optional vector of dilution factors. If provided, no user input is requested
#' @param errorsAsWarnings If TRUE, errors will be turned into warnings. Should be used with care, as errors indicate undocumented changes to the data.
#' @param alwaysShowOptionKeep If TRUE, you will be asked to keep current values even if they do not match the current dataset.
#' @return An invisible mrbin object
#' @export
#' @examples
#'  mrbinObject<-mrbin(silent=TRUE,
#'                    parameters=list(verbose=TRUE,dimension="1D",PQNScaling="No",
#'                    binwidth1D=0.04,signal_to_noise1D=1,PCA="No",binRegion=c(9.5,0.5,10,156),
#'                    saveFiles="No",referenceScaling="No",noiseRemoval="No",
#'                    fixNegatives="No",logTrafo="No",noiseThreshold=.05,tryParallel=FALSE,
#'                    NMRfolders=c(system.file("extdata/2/10/pdata/10",package="mrbin"),
#'                               system.file("extdata/3/10/pdata/10",package="mrbin"))
#'                    ))
#'  mrbinObject<-setDilutionFactors(mrbinObject,dilutionFactors=c(1.5,2))

setDilutionFactors<-function(mrbinObject,dilutionFactors=NULL,
  errorsAsWarnings=FALSE,alwaysShowOptionKeep=FALSE){
  stopTMP<-FALSE
  if(is.null(dilutionFactors)){
     optionTMP<-"Create new list"
	 yesOptionTMP<-optionTMP
	 if(alwaysShowOptionKeep) yesOptionTMP<-c("Keep current dilution factors",
	   yesOptionTMP)
     if(length(mrbinObject$metadata$dilutionFactors)==
        nrow(mrbinObject$bins)){
         yesOptionTMP<-unique(c(yesOptionTMP,"Keep current dilution factors",
		 "Edit current dilution factors"))
     }
     optionTMP<-yesOptionTMP[1]
     selectionFactors<-utils::select.list(yesOptionTMP,
                     preselect=optionTMP,
                     title = "Use previous dilution factors?",graphics=TRUE)
     if(length(selectionFactors)==0|selectionFactors=="") stopTMP<-TRUE
     if(!stopTMP&!selectionFactors=="Go back"){
       #if(selectionFactors=="Edit current dilution factors"){
         dilutionFactors<-mrbinObject$metadata$dilutionFactors
       #}
       if(selectionFactors=="Create new list"){
         dilutionFactors<-rep(1,nrow(mrbinObject$bins))
       }
       if(!selectionFactors=="Keep current dilution factors"){
        for(idilutionFactors in 1:length(dilutionFactors)){
           TMP<-readline(prompt=paste("Enter dilution factor for sample ",
             rownames(mrbinObject$bins)[idilutionFactors],
             ", enter to keep ",dilutionFactors[idilutionFactors],
             ": ",sep=""))
           if(!TMP==""&!is.na(as.numeric(TMP))) dilutionFactors[idilutionFactors]<-as.numeric(TMP)
        }
       }
     }
  } else {
	  if(!length(dilutionFactors)==nrow(mrbinObject$bins)){
		if(!errorsAsWarnings) stop("Data dimension does not match the number of provided dilution factors.")
		warning("Data dimension does not match the number of provided dilution factors.")
	  }
  }
  if(!stopTMP&!identical(dilutionFactors,mrbinObject$metadata$dilutionFactors)){
    mrbinObject<-editmrbin(mrbinObject=mrbinObject,
         functionName="mrbin::setDilutionFactors",
         versionNumber=as.character(utils::packageVersion("mrbin")),
         metadata=list(dilutionFactors=dilutionFactors),
         verbose=FALSE)
  }
  invisible(mrbinObject)
}

#' A function for setting and plotting noise levels.
#'
#' This function reads parameters from the global variable mrbin.env$mrbin$parameters and
#' plots exemplary spectra and respective noise levels. Plots will be saved if saveFiles is set to "Yes".
#' @param mrbinObject An mrbin object
#' @param plotOnly Should only noise plots be generated (TRUE), or should noise levels be adjusted interactively (FALSE)
#' @param showSpectrumPreview Should plots be shown? If not provided, this value will be taken from the mrbin object parameters
#' @param silent If set to TRUE, plots will not be shown but might still be saved
#' @return An invisible mrbin object
#' @export
#' @examples
#'  mrbinObject<-mrbin(silent=TRUE,
#'                    parameters=list(verbose=TRUE,dimension="1D",PQNScaling="No",
#'                    binwidth1D=0.04,signal_to_noise1D=1,PCA="No",binRegion=c(9.5,0.5,10,156),
#'                    saveFiles="No",referenceScaling="No",noiseRemoval="No",
#'                    fixNegatives="No",logTrafo="No",noiseThreshold=.05,tryParallel=FALSE,
#'                    NMRfolders=c(system.file("extdata/2/10/pdata/10",package="mrbin"),
#'                               system.file("extdata/3/10/pdata/10",package="mrbin"))
#'                    ))
#'  mrbinObject<-setNoiseLevels(mrbinObject,plotOnly=TRUE)

setNoiseLevels<-function(mrbinObject,plotOnly=FALSE,
  showSpectrumPreview=NULL,silent=FALSE){
    if(is.null(showSpectrumPreview)) showSpectrumPreview<-mrbinObject$parameters$showSpectrumPreview
    stopTMP<-FALSE

    if(mrbinObject$parameters$dimension=="1D") dimlength<-2
    if(mrbinObject$parameters$dimension=="2D") dimlength<-4
    if(showSpectrumPreview=="Yes"){
       if(!identical(mrbin.env$mrbin$parameters$dimension,mrbinObject$parameters$dimension)|
        (mrbinObject$parameters$dimension=="1D"&is.null(mrbin.env$mrbinTMP$additionalPlots1D))|
        (mrbinObject$parameters$dimension=="2D"&is.null(mrbin.env$mrbinTMP$additionalPlots2D))|
        !identical(mrbin.env$mrbinTMP$spectrumListPlotTMP[1],
        ceiling(length(mrbinObject$parameters$NMRfolders)*.33))|#check if current environemtn data is from this dataset
        !identical(mrbinObject$parameters$NMRfolders,
        mrbin.env$mrbin$parameters$NMRfolders)
        ){
          mrbin.env$mrbinTMP$currentFolder<-mrbinObject$parameters$NMRfolders[1]
          readNMR2(dimension=mrbinObject$parameters$dimension,
                   NMRvendor=mrbinObject$parameters$NMRvendor,
                   useAsNames=mrbinObject$parameters$useAsNames)
          mrbin.env$mrbinTMP$additionalPlots1D<-NULL
          mrbin.env$mrbinTMP$additionalPlots1DMetadata<-NULL
          mrbin.env$mrbinTMP$additionalPlots2D<-NULL
          mrbin.env$mrbinTMP$additionalPlots2DMetadata<-NULL
          #Find 3 more spectra: 33 percentile, 66 percentile, last spectrum
          mrbin.env$mrbinTMP$spectrumListPlotTMP<-setdiff(unique(c(
            ceiling(length(mrbinObject$parameters$NMRfolders)*.33),
            ceiling(length(mrbinObject$parameters$NMRfolders)*.66),
            length(mrbinObject$parameters$NMRfolders))),1)
          if(length(mrbin.env$mrbinTMP$spectrumListPlotTMP)>0){
            for(ispectrumListPlotTMP in 1:length(mrbin.env$mrbinTMP$spectrumListPlotTMP)){
              addToPlot(folder=mrbinObject$parameters$NMRfolders[
                mrbin.env$mrbinTMP$spectrumListPlotTMP[ispectrumListPlotTMP]],
                 dimension=mrbinObject$parameters$dimension,
                 NMRvendor=mrbinObject$parameters$NMRvendor,
                 useAsNames=mrbinObject$parameters$useAsNames)
            }
          }
        }
   }
   if(!plotOnly){
     adjNoiseRegionAcceptFlag<-FALSE
     adjNoiseRegion<-""
     if(mrbinObject$parameters$verbose){
        message("Hint: Region should have no real signals")
        utils::flush.console()
     }
     while(!adjNoiseRegionAcceptFlag&!stopTMP){
      if(mrbinObject$parameters$showSpectrumPreview=="Yes") par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1))
      if(mrbinObject$parameters$dimension=="1D"){
          if(mrbinObject$parameters$showSpectrumPreview=="Yes") plotMultiNMR(#region=regionTMP,
                  rectangleRegions=matrix(c(mrbinObject$parameters$noiseRange1d,-1000,1000),ncol=4),
                  color="black",manualScale=FALSE,maxPlots=2,
                  plotTitle=paste("Noise region\nleft=",mrbinObject$parameters$noiseRange1d[1],
                            "ppm, right=",mrbinObject$parameters$noiseRange1d[2],"ppm",sep=""),
                  restrictToRange=TRUE,enableSplit=FALSE,dimension=mrbinObject$parameters$dimension)
          adjNoiseRegionAccept<-paste(paste(c("Keep left=","ppm, right=","ppm, top=","ppm, bottom=")[1:dimlength],
                            mrbinObject$parameters$noiseRange1d,collapse="",sep=""),"ppm",sep="")
        } else {
          if(mrbinObject$parameters$showSpectrumPreview=="Yes") plotMultiNMR(#region=regionTMP,
                  rectangleRegions=matrix(mrbinObject$parameters$noiseRange2d,ncol=4),
                  color="black",manualScale=FALSE,maxPlots=2,
                  plotTitle=paste("Noise region\nleft=",mrbinObject$parameters$noiseRange2d[1],
                            "ppm, right=",mrbinObject$parameters$noiseRange2d[2],"ppm, top=",
                            mrbinObject$parameters$noiseRange2d[3],
                            "ppm, bottom=",mrbinObject$parameters$noiseRange2d[4],"ppm",sep=""),
                  restrictToRange=TRUE,enableSplit=FALSE,dimension=mrbinObject$parameters$dimension)

          adjNoiseRegionAccept<-paste(paste(c("Keep left=","ppm, right=","ppm, top=","ppm, bottom=")[1:dimlength],
                            mrbinObject$parameters$noiseRange2d,collapse="",sep=""),"ppm",sep="")
        }
      adjNoiseRegion<-utils::select.list(c(adjNoiseRegionAccept,"Change..."#,"Go back"
               ),
               preselect=adjNoiseRegionAccept,
               title ="Keep noise region?",graphics=TRUE)
      if(length(adjNoiseRegion)==0|adjNoiseRegion=="") stopTMP<-TRUE
       if(!stopTMP){
            if(adjNoiseRegion==adjNoiseRegionAccept) adjNoiseRegionAcceptFlag<-TRUE
            if(adjNoiseRegion=="Change..."&!adjNoiseRegion=="Go back"){
             if(mrbinObject$parameters$dimension=="1D"){
              regionTMP<-readline(prompt=paste("Noise: left border, press enter to keep ",
                        mrbinObject$parameters$noiseRange1d[1],": ",sep=""))
              if(!regionTMP=="") {
                     noiseRange1d<-mrbinObject$parameters$noiseRange1d
                     noiseRange1d[1]<-as.numeric(regionTMP)
                     mrbinObject<-editmrbin(mrbinObject=mrbinObject,functionName="mrbin::mrbin",
                       versionNumber=as.character(utils::packageVersion("mrbin")),
                       parameters=list(noiseRange1d=noiseRange1d),verbose=FALSE)
              }
              regionTMP<-readline(prompt=paste("Noise: right border, press enter to keep ",
                        mrbinObject$parameters$noiseRange1d[2],": ",sep=""))
              if(!regionTMP=="") {
                     noiseRange1d<-mrbinObject$parameters$noiseRange1d
                     noiseRange1d[2]<-as.numeric(regionTMP)
                     mrbinObject<-editmrbin(mrbinObject=mrbinObject,functionName="mrbin::mrbin",
                       versionNumber=as.character(utils::packageVersion("mrbin")),
                       parameters=list(noiseRange1d=noiseRange1d),verbose=FALSE)
              }
              if(mrbinObject$parameters$noiseRange1d[1]<mrbinObject$parameters$noiseRange1d[2]){
                TMP<-mrbinObject$parameters$noiseRange1d[1]
                noiseRange1d<-mrbinObject$parameters$noiseRange1d
                noiseRange1d[1]<-mrbinObject$parameters$noiseRange1d[2]
                noiseRange1d[2]<-TMP
                mrbinObject<-editmrbin(mrbinObject=mrbinObject,functionName="mrbin::mrbin",
                       versionNumber=as.character(utils::packageVersion("mrbin")),
                       parameters=list(noiseRange1d=noiseRange1d),verbose=FALSE)
              }
             }
            if(mrbinObject$parameters$dimension=="2D"&!stopTMP){
              regionTMP<-readline(prompt=paste("Noise: left border, press enter to keep ",
                        mrbinObject$parameters$noiseRange2d[1],": ",sep=""))
              if(!regionTMP=="") {
                     noiseRange2d<-mrbinObject$parameters$noiseRange2d
                     noiseRange2d[1]<-as.numeric(regionTMP)
                     mrbinObject<-editmrbin(mrbinObject=mrbinObject,functionName="mrbin::mrbin",
                       versionNumber=as.character(utils::packageVersion("mrbin")),
                       parameters=list(noiseRange2d=noiseRange2d),verbose=FALSE)
              }
              regionTMP<-readline(prompt=paste("Noise: right border, press enter to keep ",
                        mrbinObject$parameters$noiseRange2d[2],": ",sep=""))
              if(!regionTMP=="") {
                     noiseRange2d<-mrbinObject$parameters$noiseRange2d
                     noiseRange2d[2]<-as.numeric(regionTMP)
                     mrbinObject<-editmrbin(mrbinObject=mrbinObject,functionName="mrbin::mrbin",
                       versionNumber=as.character(utils::packageVersion("mrbin")),
                       parameters=list(noiseRange2d=noiseRange2d),verbose=FALSE)
              }
              if(mrbinObject$parameters$noiseRange2d[1]<mrbinObject$parameters$noiseRange2d[2]){
                TMP<-mrbinObject$parameters$noiseRange2d[1]
                noiseRange2d<-mrbinObject$parameters$noiseRange2d
                noiseRange2d[1]<-mrbinObject$parameters$noiseRange2d[2]
                noiseRange2d[2]<-TMP
               mrbinObject<-editmrbin(mrbinObject=mrbinObject,functionName="mrbin::mrbin",
                 versionNumber=as.character(utils::packageVersion("mrbin")),
                 parameters=list(noiseRange2d=noiseRange2d),verbose=FALSE)
              }
              regionTMP<-readline(prompt=paste("Noise: top border, press enter to keep ",
                        mrbinObject$parameters$noiseRange2d[3],": ",sep=""))
              if(!regionTMP=="") {
                     noiseRange2d<-mrbinObject$parameters$noiseRange2d
                     noiseRange2d[3]<-as.numeric(regionTMP)
                     mrbinObject<-editmrbin(mrbinObject=mrbinObject,functionName="mrbin::mrbin",
                       versionNumber=as.character(utils::packageVersion("mrbin")),
                       parameters=list(noiseRange2d=noiseRange2d),verbose=FALSE)
              }
              regionTMP<-readline(prompt=paste("Noise: bottom border, press enter to keep ",
                        mrbinObject$parameters$noiseRange2d[4],": ",sep=""))
              if(!regionTMP=="") {
                     noiseRange2d<-mrbinObject$parameters$noiseRange2d
                     noiseRange2d[4]<-as.numeric(regionTMP)
                     mrbinObject<-editmrbin(mrbinObject=mrbinObject,functionName="mrbin::mrbin",
                       versionNumber=as.character(utils::packageVersion("mrbin")),
                       parameters=list(noiseRange2d=noiseRange2d),verbose=FALSE)
              }
            if(mrbinObject$parameters$noiseRange2d[4]<mrbinObject$parameters$noiseRange2d[3]){
              TMP<-mrbinObject$parameters$noiseRange2d[3]
              noiseRange2d<-mrbinObject$parameters$noiseRange2d
              noiseRange2d[3]<-mrbinObject$parameters$noiseRange2d[4]
              noiseRange2d[4]<-TMP
              mrbinObject<-editmrbin(mrbinObject=mrbinObject,functionName="mrbin::mrbin",
                 versionNumber=as.character(utils::packageVersion("mrbin")),
                 parameters=list(noiseRange2d=noiseRange2d),verbose=FALSE)
            }
           }
          }
       }
      }
      if(mrbinObject$parameters$verbose){
        message("Hint: Higher SNR removes more noise, lower SNR allows low intensity signals.\nShown: Example regions and approximate SNR (red)")
        utils::flush.console()
      }
    }
    nextTMP<-FALSE
    stopTMP<-FALSE
    signal_to_noise1D<-mrbinObject$parameters$signal_to_noise1D
    signal_to_noise2D<-mrbinObject$parameters$signal_to_noise2D
    #load spectra to package environment for plotting if they are missing
    while(!nextTMP&!stopTMP){
      if(mrbinObject$parameters$saveFiles=="Yes"){
        grDevices::pdf(paste(mrbinObject$parameters$outputFileName,"NoisePreviews.pdf",sep=""))
      } else {
        grDevices::pdf(NULL)
      }
      grDevices::dev.control(displaylist="enable")
      oldpar<-graphics::par("mar","mfrow","mgp")
      on.exit(graphics::par(oldpar))
      devAskNewPage(ask = FALSE)
      if(mrbinObject$parameters$dimension=="1D"){
       #first, scale spectrum to bin intensities to make noise level comparable
       #find max value for each bin region, divide by bin intentisty, then take median
       if(is.null(mrbin.env$mrbinTMP$scaleFactorTMP1)){
        maxTMP1<-rep(0,nrow(mrbinObject$parameters$binRegions))
        maxTMP2<-rep(0,nrow(mrbinObject$parameters$binRegions))
        maxTMP3<-rep(0,nrow(mrbinObject$parameters$binRegions))
        NMRspectrumNames1<-as.numeric(names(mrbin.env$mrbinTMP$currentSpectrumOriginal))
        NMRspectrumNames2<-as.numeric(names(mrbin.env$mrbinTMP$additionalPlots1D[[1]]))
        if(length(mrbin.env$mrbinTMP$additionalPlots1D)>1){
         NMRspectrumNames3<-as.numeric(names(mrbin.env$mrbinTMP$additionalPlots1D[[2]]))
        }
        for(ibinTMP in 1:nrow(mrbinObject$parameters$binRegions)){
            indexTMP1<-NMRspectrumNames1<=mrbinObject$parameters$binRegions[ibinTMP,1]&
                      NMRspectrumNames1>mrbinObject$parameters$binRegions[ibinTMP,2]
            if(sum(indexTMP1)>0){
               maxTMP1[ibinTMP]<-max(mrbin.env$mrbinTMP$currentSpectrumOriginal[indexTMP1],na.rm=TRUE)
            }
            indexTMP2<-NMRspectrumNames2<=mrbinObject$parameters$binRegions[ibinTMP,1]&
                      NMRspectrumNames2>mrbinObject$parameters$binRegions[ibinTMP,2]
            if(sum(indexTMP2)>0){
               maxTMP2[ibinTMP]<-max(mrbin.env$mrbinTMP$additionalPlots1D[[1]][indexTMP2],na.rm=TRUE)
            }
            if(length(mrbin.env$mrbinTMP$additionalPlots1D)>1){
              indexTMP3<-NMRspectrumNames3<=mrbinObject$parameters$binRegions[ibinTMP,1]&
                      NMRspectrumNames3>mrbinObject$parameters$binRegions[ibinTMP,2]
              if(sum(indexTMP3)>0){
               maxTMP3[ibinTMP]<-max(mrbin.env$mrbinTMP$additionalPlots1D[[2]][indexTMP3],na.rm=TRUE)
              }
            }
       }
       mrbin.env$mrbinTMP$scaleFactorTMP1<-median((maxTMP1/mrbinObject$bins[1,])[
         mrbinObject$bins[1,]>0],na.rm=TRUE)
       mrbin.env$mrbinTMP$scaleFactorTMP2<-median((maxTMP2/
         mrbinObject$bins[mrbin.env$mrbinTMP$spectrumListPlotTMP[1],])[
         mrbinObject$bins[mrbin.env$mrbinTMP$spectrumListPlotTMP[1],]>0],na.rm=TRUE)
       if(length(mrbin.env$mrbinTMP$additionalPlots1D)>1){
         mrbin.env$mrbinTMP$scaleFactorTMP3<-median((maxTMP3/
           mrbinObject$bins[mrbin.env$mrbinTMP$spectrumListPlotTMP[2],])[
         mrbinObject$bins[mrbin.env$mrbinTMP$spectrumListPlotTMP[2],]>0],na.rm=TRUE)
       }
      }
      par(mfrow=c(3,min(3,(1+length(mrbin.env$mrbinTMP$additionalPlots1D)))),mar=c(2,2,1.1,.2))
      for(inoiseRegionTMP in 1:nrow(mrbinObject$parameters$noisePreviewRegion1D)){
        #plot 3 different regions of 3 spectra with current S-to-N level to help with deciding
          plotTitleTMP1<-""
          plotTitleTMP2<-""
          plotTitleTMP3<-""
          if(inoiseRegionTMP==1){
            plotTitleTMP1<-substr(rownames(mrbinObject$bins)[1],1,
              mrbinObject$parameters$PCAtitlelength)
            plotTitleTMP2<-substr(
              rownames(mrbinObject$bins)[mrbin.env$mrbinTMP$spectrumListPlotTMP[1]],1,
              mrbinObject$parameters$PCAtitlelength)
            if(length(mrbin.env$mrbinTMP$additionalPlots1D)>1) plotTitleTMP3<-
              substr(rownames(mrbinObject$bins)[mrbin.env$mrbinTMP$spectrumListPlotTMP[2]],1,
              mrbinObject$parameters$PCAtitlelength)

          }
          plotNMR(noise=1.0*mrbinObject$parameters$noise_level_adjusted[1]*signal_to_noise1D,
                perspective=TRUE,region=mrbinObject$parameters$noisePreviewRegion1D[inoiseRegionTMP,],
                color="green4",manualScale=FALSE,plotTitle=plotTitleTMP1,
                currentSpectrumOriginal=mrbin.env$mrbinTMP$currentSpectrumOriginal/
                mrbin.env$mrbinTMP$scaleFactorTMP1)
          plotNMR(noise=1.0*mrbinObject$parameters$noise_level_adjusted[
                mrbin.env$mrbinTMP$spectrumListPlotTMP[1]]*signal_to_noise1D,
                perspective=TRUE,region=mrbinObject$parameters$noisePreviewRegion1D[inoiseRegionTMP,],
                color="green4",manualScale=FALSE,plotTitle=plotTitleTMP2,
                currentSpectrumOriginal=mrbin.env$mrbinTMP$additionalPlots1D[[1]]/
                mrbin.env$mrbinTMP$scaleFactorTMP2)
          if(length(mrbin.env$mrbinTMP$additionalPlots1D)>1){
            plotNMR(noise=1.0*mrbinObject$parameters$noise_level_adjusted[
               mrbin.env$mrbinTMP$spectrumListPlotTMP[2]]*signal_to_noise1D,
                perspective=TRUE,region=mrbinObject$parameters$noisePreviewRegion1D[inoiseRegionTMP,],
                color="green4",manualScale=FALSE,plotTitle=plotTitleTMP3,
                currentSpectrumOriginal=mrbin.env$mrbinTMP$additionalPlots1D[[2]]/
                mrbin.env$mrbinTMP$scaleFactorTMP3)
          }
        }
     } else {#2D
      #first, scale spectrum to bin intensities to make noise level comparable
      #find max value for each bin region, divide by bin intentisty, then take median
      if(is.null(mrbin.env$mrbinTMP$scaleFactorTMP1)){
       maxTMP1<-rep(0,nrow(mrbinObject$parameters$binRegions))
       maxTMP2<-rep(0,nrow(mrbinObject$parameters$binRegions))
       maxTMP3<-rep(0,nrow(mrbinObject$parameters$binRegions))
       NMRspectrumRownames1<-as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrumOriginal))
       NMRspectrumColnames1<-as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrumOriginal))
       NMRspectrumRownames2<-as.numeric(rownames(mrbin.env$mrbinTMP$additionalPlots2D[[1]]))
       NMRspectrumColnames2<-as.numeric(colnames(mrbin.env$mrbinTMP$additionalPlots2D[[1]]))
       if(length(mrbin.env$mrbinTMP$additionalPlots2D)>1){
         NMRspectrumRownames3<-as.numeric(rownames(mrbin.env$mrbinTMP$additionalPlots2D[[2]]))
         NMRspectrumColnames3<-as.numeric(colnames(mrbin.env$mrbinTMP$additionalPlots2D[[2]]))
       }
       for(ibinTMP in 1:nrow(mrbinObject$parameters$binRegions)){
            rowsTMP1<-NMRspectrumRownames1<=mrbinObject$parameters$binRegions[ibinTMP,4]&
                    NMRspectrumRownames1>mrbinObject$parameters$binRegions[ibinTMP,3]
            colsTMP1<-NMRspectrumColnames1<=mrbinObject$parameters$binRegions[ibinTMP,1]&
                    NMRspectrumColnames1>mrbinObject$parameters$binRegions[ibinTMP,2]
            numberOfPointsPerBinTMP1<-(sum(rowsTMP1)*sum(colsTMP1))-sum(is.na(
                        mrbin.env$mrbinTMP$currentSpectrumOriginal[rowsTMP1,colsTMP1]))
            if(numberOfPointsPerBinTMP1>0){
               maxTMP1[ibinTMP]<-max(mrbin.env$mrbinTMP$currentSpectrumOriginal[rowsTMP1,colsTMP1],na.rm=TRUE)
            }
            rowsTMP2<-NMRspectrumRownames2<=mrbinObject$parameters$binRegions[ibinTMP,4]&
                    NMRspectrumRownames2>mrbinObject$parameters$binRegions[ibinTMP,3]
            colsTMP2<-NMRspectrumColnames2<=mrbinObject$parameters$binRegions[ibinTMP,1]&
                    NMRspectrumColnames2>mrbinObject$parameters$binRegions[ibinTMP,2]
            numberOfPointsPerBinTMP2<-(sum(rowsTMP2)*sum(colsTMP2))-sum(is.na(
                        mrbin.env$mrbinTMP$additionalPlots2D[[1]][rowsTMP2,colsTMP2]))
            if(numberOfPointsPerBinTMP2>0){
               maxTMP2[ibinTMP]<-max(mrbin.env$mrbinTMP$additionalPlots2D[[1]][rowsTMP2,colsTMP2],na.rm=TRUE)
            }
            if(length(mrbin.env$mrbinTMP$additionalPlots2D)>1){
              rowsTMP3<-NMRspectrumRownames3<=mrbinObject$parameters$binRegions[ibinTMP,4]&
                    NMRspectrumRownames3>mrbinObject$parameters$binRegions[ibinTMP,3]
              colsTMP3<-NMRspectrumColnames3<=mrbinObject$parameters$binRegions[ibinTMP,1]&
                    NMRspectrumColnames3>mrbinObject$parameters$binRegions[ibinTMP,2]
              numberOfPointsPerBinTMP3<-(sum(rowsTMP3)*sum(colsTMP3))-sum(is.na(
                        mrbin.env$mrbinTMP$additionalPlots2D[[2]][rowsTMP3,colsTMP3]))
              if(numberOfPointsPerBinTMP3>0){
                maxTMP3[ibinTMP]<-max(mrbin.env$mrbinTMP$additionalPlots2D[[2]][rowsTMP3,colsTMP3],na.rm=TRUE)
              }
            }
       }
       mrbin.env$mrbinTMP$scaleFactorTMP1<-median((maxTMP1/mrbinObject$bins[1,])[
         mrbinObject$bins[1,]>0],na.rm=TRUE)
       mrbin.env$mrbinTMP$scaleFactorTMP2<-median((maxTMP2/
         mrbinObject$bins[mrbin.env$mrbinTMP$spectrumListPlotTMP[1],])[
         mrbinObject$bins[mrbin.env$mrbinTMP$spectrumListPlotTMP[1],]>0],na.rm=TRUE)
       if(length(mrbin.env$mrbinTMP$additionalPlots2D)>1){
         mrbin.env$mrbinTMP$scaleFactorTMP3<-median((maxTMP3/
         mrbinObject$bins[mrbin.env$mrbinTMP$spectrumListPlotTMP[2],])[
         mrbinObject$bins[mrbin.env$mrbinTMP$spectrumListPlotTMP[2],]>0],na.rm=TRUE)
       }
      }
      par(mfrow=c(3,min(3,(1+length(mrbin.env$mrbinTMP$additionalPlots2D)))),mar=c(2,2,1.1,.2))
      for(inoiseRegionTMP in 1:nrow(mrbinObject$parameters$noisePreviewRegion2D)){
          plotTitleTMP1<-""
          plotTitleTMP2<-""
          plotTitleTMP3<-""
          if(inoiseRegionTMP==1){
            plotTitleTMP1<-substr(rownames(mrbinObject$bins)[1],1,
              mrbinObject$parameters$PCAtitlelength)
            plotTitleTMP2<-substr(rownames(mrbinObject$bins)[mrbin.env$mrbinTMP$spectrumListPlotTMP[1]],1,
              mrbinObject$parameters$PCAtitlelength)
            if(length(mrbin.env$mrbinTMP$additionalPlots2D)>1) plotTitleTMP3<-substr(
              rownames(mrbinObject$bins)[mrbin.env$mrbinTMP$spectrumListPlotTMP[2]],1,
              mrbinObject$parameters$PCAtitlelength)
          }
          plotNMR(noise=1.0*mrbinObject$parameters$noise_level_adjusted[1]*signal_to_noise2D,
                perspective=TRUE,region=mrbinObject$parameters$noisePreviewRegion2D[inoiseRegionTMP,],
                color="black",manualScale=FALSE,plotTitle=plotTitleTMP1,
                currentSpectrumOriginal=mrbin.env$mrbinTMP$currentSpectrumOriginal/
                mrbin.env$mrbinTMP$scaleFactorTMP1)
          plotNMR(noise=1.0*mrbinObject$parameters$noise_level_adjusted[
                mrbin.env$mrbinTMP$spectrumListPlotTMP[1]]*signal_to_noise2D,
                perspective=TRUE,region=mrbinObject$parameters$noisePreviewRegion2D[inoiseRegionTMP,],
                color="black",manualScale=FALSE,plotTitle=plotTitleTMP2,
                currentSpectrumOriginal=mrbin.env$mrbinTMP$additionalPlots2D[[1]]/
                mrbin.env$mrbinTMP$scaleFactorTMP2)
          if(length(mrbin.env$mrbinTMP$additionalPlots2D)>1){
            plotNMR(noise=1.0*mrbinObject$parameters$noise_level_adjusted[
                mrbin.env$mrbinTMP$spectrumListPlotTMP[2]]*signal_to_noise2D,
                perspective=TRUE,region=mrbinObject$parameters$noisePreviewRegion2D[inoiseRegionTMP,],
                color="black",manualScale=FALSE,plotTitle=plotTitleTMP3,
                currentSpectrumOriginal=mrbin.env$mrbinTMP$additionalPlots2D[[2]]/
                mrbin.env$mrbinTMP$scaleFactorTMP3)
          }
      }
    }
    plotRecordingTMP<-grDevices::recordPlot()
    grDevices::dev.off()
    if(showSpectrumPreview=="Yes"&!silent){
      grDevices::replayPlot(plotRecordingTMP)
    }
    if(plotOnly){
        nextTMP<-TRUE
    } else {
     if(mrbinObject$parameters$dimension=="1D"&!stopTMP){
      SNRTMP<-utils::select.list(c(paste("Keep ",signal_to_noise1D,sep=""),
              "Increase by 2","Reduce by 2","Change..."),preselect=paste("Keep ",signal_to_noise1D,sep=""),
              title="Signal-to-noise ratio (SNR):",graphics=TRUE)
      if(length(SNRTMP)==0|SNRTMP=="") stopTMP<-TRUE
      if(!stopTMP){
        if(SNRTMP==paste("Keep ",signal_to_noise1D,sep="")){
          nextTMP<-TRUE
          if(!signal_to_noise1D==mrbinObject$parameters$signal_to_noise1D){
           mrbinObject<-editmrbin(mrbinObject=mrbinObject,
             functionName="mrbin::setNoiseLevels",
             versionNumber=as.character(utils::packageVersion("mrbin")),
             parameters=list(signal_to_noise1D=signal_to_noise1D),
             verbose=FALSE)
          }
        }
        if(SNRTMP=="Increase by 2") signal_to_noise1D<-signal_to_noise1D+2
        if(SNRTMP=="Reduce by 2") signal_to_noise1D<-signal_to_noise1D-2
        if(SNRTMP=="Change..."){
          SNRTMP<-readline(prompt=paste("New signal to noise ratio, press enter to keep ",
                             signal_to_noise1D,": ",sep=""))
          if(!SNRTMP==""&!is.null(as.numeric(SNRTMP))){
           signal_to_noise1D<-as.numeric(SNRTMP)
          }
        }
      }
     }
     if(mrbinObject$parameters$dimension=="2D"&!stopTMP){
      SNRTMP<-utils::select.list(c(paste("Keep ",signal_to_noise2D,sep=""),
              "Increase by 2","Reduce by 2","Change..."),preselect=paste("Keep ",signal_to_noise2D,sep=""),
              title="Signal-to-noise ratio (SNR):",graphics=TRUE)
      if(length(SNRTMP)==0|SNRTMP=="") stopTMP<-TRUE
      if(!stopTMP){
        if(SNRTMP==paste("Keep ",signal_to_noise2D,sep="")){
          nextTMP<-TRUE
          if(!signal_to_noise2D==mrbinObject$parameters$signal_to_noise2D){
           mrbinObject<-editmrbin(mrbinObject=mrbinObject,
             functionName="mrbin::setNoiseLevels",
             versionNumber=as.character(utils::packageVersion("mrbin")),
             parameters=list(signal_to_noise2D=signal_to_noise2D),
             verbose=FALSE)
          }
        }
        if(SNRTMP=="Increase by 2") signal_to_noise2D<-signal_to_noise2D+2
        if(SNRTMP=="Reduce by 2") signal_to_noise2D<-signal_to_noise2D-2
        if(SNRTMP=="Change..."){
          SNRTMP<-readline(prompt=paste("New signal to noise ratio, press enter to keep ",
                             signal_to_noise2D,": ",sep=""))
          if(!SNRTMP==""&!is.null(as.numeric(SNRTMP))){
            signal_to_noise2D<-as.numeric(SNRTMP)
          }
        }
      }
     }
    }
   }
    if(!plotOnly){
      if(!stopTMP&mrbinObject$parameters$verbose){
        message("Hint: 0.75 if metabolites are always present (serum), 0.2 if metabolites are \nabsent in some samples (urine, cell culture)")
        utils::flush.console()
      }
      if(showSpectrumPreview=="Yes"){#show a plot of number of bins left at different thresholds 0.05-1
        thresholdValuesTMP<-sort(unique(c(as.numeric(mrbinObject$parameters$noiseThreshold),
          0.05*(1:20))))
        nBins<-rep(NA,length(thresholdValuesTMP))
        for(i in 1:length(thresholdValuesTMP)){
          mrbinObjectTMP<-editmrbin(mrbinObject=mrbinObject,
             functionName="mrbin::setNoiseLevels",
             versionNumber="0",
             parameters=list(noiseThreshold=thresholdValuesTMP[i]),
             verbose=FALSE)
          nBinsTMP<-ncol(removeNoise(mrbinObjectTMP,verbose=FALSE,errorsAsWarnings=TRUE)$bins)
          if(i==1){
            nBins[i]<-nBinsTMP
          } else {
            if(nBinsTMP==ncol(mrbinObject$bins)&
             nBinsTMP>nBins[i-1]) nBinsTMP<-0#this usually means no bins left after noise removal
            nBins[i]<-nBinsTMP
          }
        }
        dev.off()
        graphics::plot(thresholdValuesTMP,nBins,type="l",main="Noise threshold preview",
          xlab="(More noise)                                                                      (Less noise)\nNoise threshold",
          ylab="Number of bins after noise removal",
          xlim=c(0,1),ylim=c(.96*min(nBins),1.04*max(nBins)))
        graphics::abline(h=ncol(mrbinObject$bins),col="blue",lty=2)
        graphics::abline(v=#unique(c(
          as.numeric(mrbinObject$parameters$noiseThreshold)#,0.2,0.75,0.05))
          ,col="red",lty=2)
      }
      if(!stopTMP){
        noiseTMP<-utils::select.list(unique(c(as.character(mrbinObject$parameters$noiseThreshold),
                               "0.2","0.75","0.05","Change...")),
                               preselect=as.character(mrbinObject$parameters$noiseThreshold),
                               title="Minimum percentage > noise",graphics=TRUE)
        if(length(noiseTMP)==0|noiseTMP=="") stopTMP<-TRUE
        if(!stopTMP){
          if(noiseTMP=="Change..."&!stopTMP){
                noiseTMP<-readline(prompt=paste("New noise threshold, press enter to keep ",
                          mrbinObject$parameters$noiseThreshold,": ",sep=""))
          }
          if(!noiseTMP==""&!is.na(as.numeric(noiseTMP))&!as.numeric(noiseTMP)==
             mrbinObject$parameters$noiseThreshold)  {
             mrbinObject<-editmrbin(mrbinObject=mrbinObject,
               functionName="mrbin::setNoiseLevels",
               versionNumber=as.character(utils::packageVersion("mrbin")),
               parameters=list(noiseThreshold=as.numeric(noiseTMP)),
               verbose=FALSE)
          }
        }
      }
    }
    invisible(mrbinObject)
}


#' A function for interactively editing metadata of mrbin objects.
#'
#' This function edits interactively or non-interactively the metadata filed of the provided mrbin object.
#' @param mrbinResults An mrbin object
#' @param metadata An optional list of objects to be changed. If provided, interactive mode is deactivated
#' @return An invisible mrbin object
#' @export
#' @examples
#'  mrbinObject<-mrbin(silent=TRUE,
#'                    parameters=list(verbose=TRUE,dimension="1D",PQNScaling="No",
#'                    binwidth1D=0.04,signal_to_noise1D=1,PCA="No",binRegion=c(9.5,0.5,10,156),
#'                    saveFiles="No",referenceScaling="No",noiseRemoval="No",
#'                    fixNegatives="No",logTrafo="No",noiseThreshold=.05,tryParallel=FALSE,
#'                    NMRfolders=c(system.file("extdata/2/10/pdata/10",package="mrbin"),
#'                               system.file("extdata/3/10/pdata/10",package="mrbin"))
#'                    ))
#'  mrbinObject<-metadatamrbin(mrbinObject,metadata=list(projectTitle="Test project"))

metadatamrbin<-function(mrbinResults,metadata=NULL){
  if(!is.null(metadata)){
     mrbinResults<-editmrbin(mrbinResults,functionName="mrbin::metadatamrbin",
              versionNumber=as.character(utils::packageVersion("mrbin")),
              metadata=metadata,verbose=FALSE)
  } else {
    i<-1
    stopTMP<-FALSE
    iValues<-setdiff(names(mrbinResults$metadata),c("annotations","metaData"))
    while(i <= length(iValues)){
     if(!stopTMP){
      if(i<1) i<-1
       message(paste("Preview of ",iValues[i],": ",
         substr(paste(mrbinResults$metadata[[iValues[i]]],sep=" ",collapse=" "),1,45),sep=""))
       utils::flush.console()
       if(iValues[i]=="dilutionFactors"){#Define dilution
          mrbinResults<-setDilutionFactors(mrbinResults,alwaysShowOptionKeep=TRUE)
       } else {
         if(iValues[i]=="factors"){#Define groups
          selectionFactors<-""
          defineGroups<-utils::select.list(c("Edit group memberships",#"Use a column from metadata data.frame",
                                    "Keep","Go back"),
                                    preselect="Keep",
                                    title = "Edit classes/groups for PCA",graphics=TRUE)
          if(length(defineGroups)==0|defineGroups=="") stopTMP<-TRUE
          if(!stopTMP){
           if(defineGroups=="Go back") i<-i-2
           if(defineGroups=="Edit group memberships"){
               if(!is.null(mrbinResults$metadata$factors)){
                 if(length(mrbinResults$metadata$factors)==length(mrbinResults$parameters$NMRfolders)){
                   yesOptionTMP<-paste("Use previous factor list (",
                                       paste(mrbinResults$metadata$factors[1:min(3,length(mrbinResults$metadata$factors))],
                                             sep=", ",collapse=", "),
                                       ", ...)",sep="")
                   selectionFactors<-utils::select.list(c(yesOptionTMP,"No"),
                                     preselect=yesOptionTMP,
                                     title = "Use previous factor list?",graphics=TRUE)
                   if(length(selectionFactors)==0|selectionFactors=="") stopTMP<-TRUE
                   if(!stopTMP&selectionFactors=="No")  mrbinResults<-setFactors(mrbinResults)
                 } else {
                   mrbinResults<-setFactors(mrbinResults)
                 }
              } else {
                mrbinResults<-setFactors(mrbinResults)
              }
            }
          }
       } else {
          defineGroups<-utils::select.list(c(paste("Edit ",iValues[i],sep=""),"Keep","Replace","Go back"),
                                preselect="Keep",title=paste("Edit ",iValues[i],sep=""),graphics=TRUE)
          if(length(defineGroups)==0|defineGroups=="") stopTMP<-TRUE
          if(!stopTMP){
           if(defineGroups=="Go back") i<-i-2
           if(!defineGroups=="Keep"){
            TMP<-NULL
            if(defineGroups==paste("Edit ",iValues[i],sep="")){
             TMP<-utils::edit(mrbinResults$metadata[[iValues[i]]])
            }
            if(defineGroups=="Replace"){
             TMP<-readline(prompt=paste("Type new content, press enter to keep current content: ",sep=""))
            }
            if(!is.null(TMP)){
               TMPlist<-list(TMP)
               names(TMPlist)<-iValues[i]
               mrbinResults<-editmrbin(mrbinResults,functionName="mrbin::metadatamrbin",
                        versionNumber=as.character(utils::packageVersion("mrbin")),
                        metadata=TMPlist,verbose=FALSE)

            }
           }
          }
         }
       }
     }
     i<-i+1
    }
  }
  mrbinResults<-annotatemrbin(mrbinResults)
  invisible(mrbinResults)
}

#' A function for printing parameters to the screen.
#'
#' This function reads parameters from the global variable mrbin.env$mrbin$parameters and
#' prints the required R code for creating a data set to the screen.
#' @return {None}
#' @export
#' @examples
#' printParameters()

printParameters<-function(){
  if(!exists("mrbin.env", mode="environment")) .onLoad()
  if(mrbin.env$mrbin$parameters$verbose){
    if(mrbin.env$mrbin$parameters$dimension=="1D") paramNamesTMP0<-mrbin.env$requiredParam1D
    if(mrbin.env$mrbin$parameters$dimension=="2D") paramNamesTMP0<-mrbin.env$requiredParam2D
    paramNamesTMP<-NULL
    for(irequiredMetadata in paramNamesTMP0){
      if(!is.null(mrbin.env$mrbin$parameters[[irequiredMetadata]])){
        if(!identical(mrbin.env$mrbin$parameters[[irequiredMetadata]],"")){
          if(is.matrix(mrbin.env$mrbin$parameters[[irequiredMetadata]])){
            if(nrow(mrbin.env$mrbin$parameters[[irequiredMetadata]])>0){
              paramNamesTMP<-c(paramNamesTMP,irequiredMetadata)
            }
          } else {
            paramNamesTMP<-c(paramNamesTMP,irequiredMetadata)
          }
        }
      }
    }
    printTMPline<-paste("results<-mrbin(silent=FALSE,parameters=list(")
    printTMP<-NULL
    counter<-0
    for(i in paramNamesTMP){
      counter<-counter+1
      if(is.character(mrbin.env$mrbin$parameters[[i]])){
        sepSymbol<-"\""
      } else {
        sepSymbol<-""
      }
      returnSymbol<-"\n  "
      if(!i=="NMRfolders"){
        returnSymbol<-""
      }
      vectorSymbol1<-""
      vectorSymbol2<-""
      if(is.matrix(mrbin.env$mrbin$parameters[[i]])){
        if(nrow(mrbin.env$mrbin$parameters[[i]])==0){
          vectorSymbol1<-"matrix(nrow=0"
          vectorSymbol2<-paste(",ncol=",ncol(mrbin.env$mrbin$parameters[[i]]),")",sep="")
        } else {
          vectorSymbol1<-"matrix(c(\n  "
          if(is.null(rownames(mrbin.env$mrbin$parameters[[i]]))){
            vectorSymbol2<-paste("\n  ),ncol=",ncol(mrbin.env$mrbin$parameters[[i]]),",byrow=TRUE)",sep="")
          } else {
            vectorSymbol2<-paste("\n  ),ncol=",ncol(mrbin.env$mrbin$parameters[[i]]),
                           ",dimnames=list(c(\n  \"",
                           paste(rownames(mrbin.env$mrbin$parameters[[i]]),sep="\",\"",collapse="\",\""),
                           "\"\n  ),NULL),byrow=TRUE",
                           ")",sep="")
          }
        }
      }
      if(length(mrbin.env$mrbin$parameters[[i]])>1){
        if(is.vector(mrbin.env$mrbin$parameters[[i]])){
          vectorSymbol1<-"c("
          vectorSymbol2<-")"
          if(i=="NMRfolders"){
            vectorSymbol1<-"c(\n  "
            vectorSymbol2<-"\n )"
          }
        }
      }
      if(is.factor(mrbin.env$mrbin$parameters[[i]])){
        vectorSymbol1<-"factor(c("
        vectorSymbol2<-"))"
        if(nchar(paste(mrbin.env$mrbin$parameters[[i]],sep=",",collapse=","))>60){
          vectorSymbol1<-"factor(c(\n  "
          vectorSymbol2<-"\n  ))"
        }
        sepSymbol<-"\""
      }
      sepTMP<-paste(sepSymbol,",",returnSymbol,sepSymbol,sep="")
      sepTMPReturn<-paste(sepSymbol,",","\n  ",sepSymbol,sep="")
      if(is.null(mrbin.env$mrbin$parameters[[i]])){
          valueTMP<-"NULL"
      } else {
        if(is.factor(mrbin.env$mrbin$parameters[[i]])){
          lengthCutOff2<-8
          if(length(mrbin.env$mrbin$parameters[[i]])<=lengthCutOff2){
            valueTMP<-paste(as.character(mrbin.env$mrbin$parameters[[i]]),sep=sepTMP,collapse=sepTMP)
          } else {
            valueTMP<-NULL
            for(iFactorLengthTMP in 1:ceiling(length(mrbin.env$mrbin$parameters[[i]])/lengthCutOff2)){
              sepFactorTMPTMP<-sepTMPReturn
              if(iFactorLengthTMP==ceiling(length(mrbin.env$mrbin$parameters[[i]])/lengthCutOff2)){
                sepFactorTMPTMP<-""
              }
              valueTMP<-paste(valueTMP,
                        paste(as.character(mrbin.env$mrbin$parameters[[i]])[((iFactorLengthTMP-1)*lengthCutOff2+1):
                               min(length(mrbin.env$mrbin$parameters[[i]]),(iFactorLengthTMP-1)*lengthCutOff2+lengthCutOff2)],
                          sep=sepTMP,collapse=sepTMP),sepFactorTMPTMP,
                        sep="",collapse="")
            }
          }
        } else {
          lengthCutOff<-12
          if(length(mrbin.env$mrbin$parameters[[i]])<=lengthCutOff){
            valueTMP<-paste(t(mrbin.env$mrbin$parameters[[i]]),sep=sepTMP,collapse=sepTMP)
          } else {
            valueTMP<-NULL
            for(iFactorLengthTMP in 1:ceiling(length(mrbin.env$mrbin$parameters[[i]])/lengthCutOff)){
              sepFactorTMPTMP<-sepTMPReturn
              if(iFactorLengthTMP==ceiling(length(mrbin.env$mrbin$parameters[[i]])/lengthCutOff)){
                sepFactorTMPTMP<-""
              }
              valueTMP<-paste(valueTMP,
                        paste(t(mrbin.env$mrbin$parameters[[i]])[((iFactorLengthTMP-1)*lengthCutOff+1):
                               min(length(mrbin.env$mrbin$parameters[[i]]),(iFactorLengthTMP-1)*lengthCutOff+lengthCutOff)],
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
        if((nchar(printTMPline)+nchar(TMPline))<=61){
          printTMP<-paste(printTMP,printTMPline,TMPline,")",sep="")
        } else {
           printTMP<-paste(printTMP,printTMPline,"\n ",sep="")
           printTMP<-paste(printTMP,TMPline,")",sep="")
        }
      }
    }
    #now metadata:
    paramNamesTMP<-NULL
    for(irequiredMetadata in mrbin.env$requiredMetadata){
      if(!is.null(mrbin.env$mrbin$metadata[[irequiredMetadata]])){
        if(!identical(mrbin.env$mrbin$metadata[[irequiredMetadata]],"")){
          if(is.matrix(mrbin.env$mrbin$metadata[[irequiredMetadata]])){
            if(nrow(mrbin.env$mrbin$metadata[[irequiredMetadata]])>0){
              paramNamesTMP<-c(paramNamesTMP,irequiredMetadata)
            }
          } else {
            paramNamesTMP<-c(paramNamesTMP,irequiredMetadata)
          }
        }
      }
    }
    if(!is.null(paramNamesTMP)){
        printTMPline<-",\n metadata=list("
        counter<-0
        for(i in paramNamesTMP){
          counter<-counter+1
          if(is.character(mrbin.env$mrbin$metadata[[i]])){
            sepSymbol<-"\""
          } else {
            sepSymbol<-""
          }
          returnSymbol<-""
          vectorSymbol1<-""
          vectorSymbol2<-""
          if(is.matrix(mrbin.env$mrbin$metadata[[i]])){
            if(nrow(mrbin.env$mrbin$metadata[[i]])==0){
              vectorSymbol1<-"matrix(nrow=0"
              vectorSymbol2<-paste(",ncol=",ncol(mrbin.env$mrbin$metadata[[i]]),")",sep="")
            } else {
              vectorSymbol1<-"matrix(c(\n  "
              if(is.null(rownames(mrbin.env$mrbin$metadata[[i]]))){
                vectorSymbol2<-paste("\n  ),ncol=",ncol(mrbin.env$mrbin$metadata[[i]]),
                ",byrow=TRUE)",sep="")
              } else {
                vectorSymbol2<-paste("\n  ),ncol=",ncol(mrbin.env$mrbin$metadata[[i]]),
                               ",dimnames=list(c(\n  \"",
                               paste(rownames(mrbin.env$mrbin$metadata[[i]]),sep="\",\"",
                               collapse="\",\""),
                               "\"\n  ), ",
                               "NULL",
                               #"c(\n  \"",paste(colnames(mrbin.env$mrbin$metadata[[i]]),sep="\",\"",
                               #collapse="\",\""),
                               #"\")",
                               "),byrow=TRUE",
                               ")",sep="")
              }
            }
          }
          if(length(mrbin.env$mrbin$metadata[[i]])>1){
            if(is.vector(mrbin.env$mrbin$metadata[[i]])){
              vectorSymbol1<-"c("
              vectorSymbol2<-")"
            }
          }
          if(is.factor(mrbin.env$mrbin$metadata[[i]])){
            vectorSymbol1<-"factor(c("
            vectorSymbol2<-"))"
            if(nchar(paste(mrbin.env$mrbin$metadata[[i]],sep=",",collapse=","))>60){
              vectorSymbol1<-"factor(c(\n  "
              vectorSymbol2<-"\n  ))"
            }
            sepSymbol<-"\""
          }
          sepTMP<-paste(sepSymbol,",",returnSymbol,sepSymbol,sep="")
          sepTMPReturn<-paste(sepSymbol,",","\n  ",sepSymbol,sep="")
          if(is.null(mrbin.env$mrbin$metadata[[i]])){
              valueTMP<-"NULL"
          } else {
            if(is.factor(mrbin.env$mrbin$metadata[[i]])){
              lengthCutOff2<-8
              if(length(mrbin.env$mrbin$metadata[[i]])<=lengthCutOff2){
                valueTMP<-paste(as.character(mrbin.env$mrbin$metadata[[i]]),sep=sepTMP,
                collapse=sepTMP)
              } else {
                valueTMP<-NULL
                for(iFactorLengthTMP in 1:ceiling(length(mrbin.env$mrbin$metadata[[i]])/
                  lengthCutOff2)){
                  sepFactorTMPTMP<-sepTMPReturn
                  if(iFactorLengthTMP==ceiling(length(mrbin.env$mrbin$metadata[[i]])/
                    lengthCutOff2)){
                    sepFactorTMPTMP<-""
                  }
                  valueTMP<-paste(valueTMP,
                            paste(as.character(mrbin.env$mrbin$metadata[[i]])[
                                   ((iFactorLengthTMP-1)*lengthCutOff2+1):
                                   min(length(mrbin.env$mrbin$metadata[[i]]),
                                   (iFactorLengthTMP-1)*lengthCutOff2+lengthCutOff2)],
                              sep=sepTMP,collapse=sepTMP),sepFactorTMPTMP,
                            sep="",collapse="")
                }
              }
            } else {
              lengthCutOff<-12
              if(length(mrbin.env$mrbin$metadata[[i]])<=lengthCutOff){
                valueTMP<-paste(t(mrbin.env$mrbin$metadata[[i]]),sep=sepTMP,collapse=sepTMP)
              } else {
                valueTMP<-NULL
                for(iFactorLengthTMP in 1:ceiling(length(mrbin.env$mrbin$metadata[[i]])/
                  lengthCutOff)){
                  sepFactorTMPTMP<-sepTMPReturn
                  if(iFactorLengthTMP==ceiling(length(mrbin.env$mrbin$metadata[[i]])/
                    lengthCutOff)){
                    sepFactorTMPTMP<-""
                  }
                  valueTMP<-paste(valueTMP,
                            paste(t(mrbin.env$mrbin$metadata[[i]])[((iFactorLengthTMP-1)*
                                   lengthCutOff+1):
                                   min(length(mrbin.env$mrbin$metadata[[i]]),
                                   (iFactorLengthTMP-1)*lengthCutOff+lengthCutOff)],
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
              printTMP<-paste(printTMP,printTMPline,TMPline,")",sep="")
            } else {
               printTMP<-paste(printTMP,printTMPline,"\n ",sep="")
               printTMP<-paste(printTMP,TMPline,")",sep="")
            }
          }
        }
      }
      printTMP<-paste(
        "\n\nTo recreate this data set, use the following code:\n\n##################################################\n\n",
        printTMP,sep="")
      printTMP<-paste(printTMP,
        ")\n\n##################################################\n\nTo recreate this data set, use the code above this line.\n",
        sep="")
      if(mrbin.env$mrbin$parameters$saveFiles=="Yes"){
        printTMP<-paste(printTMP,"To load the saved data set, use the following code:\n\n",
          #" data<-read.csv(\n  \"",mrbin.env$mrbin$parameters$outputFileName,
          #      "bins.csv\",\n",
          #"  check.names = FALSE,row.names=1)\n\n",
          " load(\n  \"",mrbin.env$mrbin$parameters$outputFileName,
                ".Rdata\"\n )\n\n",
          sep="")
      }
      #printTMP<-paste(printTMP,"))\n",sep="")
      #message("\nTo recreate this data set, use the following code:\n\n##################################################\n\n", appendLF = FALSE)
      #Necessary as message() cuts off strings after around 8187 characters (Windows 10)
      messageMaxLength<-2000
      for(iLengthMessage in 1:ceiling(nchar(printTMP)/messageMaxLength)){
        message(substr(printTMP,
                      (iLengthMessage-1)*messageMaxLength+1,
                      (iLengthMessage-1)*messageMaxLength+messageMaxLength),
                appendLF = FALSE)
      }
      #message("\n##################################################\n\nTo recreate this data set, use the code above this line.\n", appendLF = FALSE)
      #if(mrbin.env$mrbin$parameters$saveFiles=="Yes"){
      #  message("To load the saved data set, use the following code:\n\n", appendLF = FALSE)
      #  message(" data<-read.csv(\n  \"",mrbin.env$mrbin$parameters$outputFileName,
      #          "bins.csv\",\n", appendLF = FALSE)
      #  message("  check.names = FALSE,row.names=1)\n\n", appendLF = FALSE)
      #}
      utils::flush.console()
      invisible(printTMP)
  }
}


#' A function recreating parameters from previous runs.
#'
#' This function reads parameters from a text file that was created during a
#' previous run or mrbin(). After reading, the data can be recreated using
#' mrbin(). File names in $parameters might need to be updated.
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
#' @param metadata List of metadata to be set
#' @return {None}
#' @export
#' @examples
#' setParam(parameters=list(dimension="1D"))

setParam<-function(parameters=NULL,metadata=NULL){
  if(!is.null(parameters)){
    mrbin.env$mrbin$parameters_copy<-mrbin.env$mrbin$parameters
    if("dimension"%in%names(parameters)){
      dimTMP<-parameters$dimension
    } else {
      dimTMP<-mrbin.env$mrbin$parameters$dimension
    }
    if(dimTMP=="1D") diffSet3<-setdiff(mrbin.env$requiredParam1D,names(parameters))
    if(dimTMP=="2D") diffSet3<-setdiff(mrbin.env$requiredParam2D,names(parameters))
    diffSet2<-setdiff(names(parameters)[!sapply(parameters, is.null)],names(mrbin.env$mrbin$parameters_copy))
    intersectSet<-intersect(names(parameters),names(mrbin.env$mrbin$parameters_copy))
    if(length(diffSet2)>0){
       warning(paste("Unexpected parameters: ",
           paste(diffSet2,sep=", ", collapse=", "),"\n",
           "These parameters are not used. Potentially they were created in a different mrbin version.",
           sep=""))
    }
    if(length(intersectSet)>0){
       for(iintersectSet in intersectSet){
            mrbin.env$mrbin$parameters[[iintersectSet]]<-parameters[[iintersectSet]]
       }
    }
    if(!as.character(utils::packageVersion("mrbin"))==mrbin.env$mrbin$parameters$mrbinversionTracker){
       warning(paste("Imported file was created using another mrbin version: ",mrbin.env$mrbin$parameters$mrbinversionTracker,
           ".\n For exact reproduction of results, please get the old version at: www.kleinomicslab.com",
           sep=""))
    }
  }
   if(!is.null(metadata)){
    if(length(metadata)>0){
       for(iintersectSet in names(metadata)){
            mrbin.env$mrbin$metadata[[iintersectSet]]<-metadata[[iintersectSet]]
       }
    }
  }
}

#' A function for log transforming data.
#'
#' This function performs logarithm transformation. Will not work with negative data.
#' @param mrbinResults An mrbin object
#' @param verbose Should a summary be printed?
#' @param errorsAsWarnings If TRUE, errors will be turned into warnings. Should be used with care, as errors indicate undocumented changes to the data.
#' @return An invisible mrbin object
#' @export
#' @examples
#' resetEnv()
#' mrbinResults<-mrbin(silent=TRUE,setDefault=TRUE,parameters=list(dimension="1D", logTrafo="No",
#'                     binwidth1D=0.05,signal_to_noise1D=50,verbose=TRUE,PCA="No",tryParallel=FALSE,
#'                     NMRfolders=c(system.file("extdata/1/10/pdata/10",package="mrbin"),
#'                                 system.file("extdata/2/10/pdata/10",package="mrbin"))))
#' mrbinResults<-logTrafo(mrbinResults)

logTrafo<-function(mrbinResults,verbose=TRUE,errorsAsWarnings=FALSE){
  mrbinResults2<-mrbinResults
    transformations="Log transformed"
    if(transformations %in% mrbinResults$transformations){
      if(!errorsAsWarnings) stop("Data has been log transformed previously, this could corrupt the data.")
      warning("Data has been log transformed previously, this could corrupt the data.")
    }
   if(sum(mrbinResults$bins<=0)>0){
       stop("Log transform does not work with negative values.")
   } else {
      if(nrow(mrbinResults$bins)==1){
        rownamesTMP<-rownames(mrbinResults$bins)
        colnamesTMP<-colnames(mrbinResults$bins)
        mrbinResults$bins<-matrix(log(mrbinResults$bins),nrow=1)
        rownames(mrbinResults$bins)<-rownamesTMP
        colnames(mrbinResults$bins)<-colnamesTMP
      } else {
       mrbinResults$bins<-log(mrbinResults$bins)
      }
   }
   mrbinResults2<-editmrbin(mrbinObject=mrbinResults2,functionName="mrbin::logTrafo",
       versionNumber=as.character(utils::packageVersion("mrbin")),
       bins=mrbinResults$bins, parameters=mrbinResults$parameters,transformations=transformations,
       verbose=verbose)
  invisible(mrbinResults2)
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
       newCurrent<-utils::select.list(mrbin.env$mrbin$parameters$NMRfolders,preselect=mrbin.env$mrbinTMP$currentFolder,
                                          title="Select new current spectrum.",graphics=TRUE)
     } else {
       if(length(mrbin.env$mrbin$parameters$NMRfolders)>=spectrumNumber){
         newCurrent<- mrbin.env$mrbin$parameters$NMRfolders[spectrumNumber]
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
#' analysis.
#' @param mrbinResults An mrbin object. If not provided, the function works on the package environment
#' @param spectra Character vector with NMR folder names to be excluded. If provided, no interactive selection will be shown
#' @param verbose Should a summary be printed?
#' @param errorsAsWarnings If TRUE, errors will be turned into warnings. Should be used with care, as errors indicate undocumented changes to the data.
#' @return An invisible mrbin object (only if an mrbin object was provided)
#' @export
#' @examples
#' mrbinResults<-mrbin(silent=TRUE,setDefault=TRUE,parameters=list(dimension="1D",
#'          binwidth1D=0.05,PQNScaling="No",PCA="No",tryParallel=FALSE,logTrafo="No",
#'          noiseRemoval="No",
#'          NMRfolders=c(system.file("extdata/1/10/pdata/10",package="mrbin"),
#'                       system.file("extdata/2/10/pdata/10",package="mrbin"),
#'                       system.file("extdata/3/10/pdata/10",package="mrbin"))))
#' mrbinResults<-removeSpectrum(mrbinResults,
#'  spectra=c(system.file("extdata/2/10/pdata/10",package="mrbin")))

removeSpectrum<-function(mrbinResults=NULL,spectra=NULL,verbose=TRUE,errorsAsWarnings=FALSE){
 if(!is.null(mrbinResults)){
    transformations="Spectra removed"
    if("Noise removed" %in% mrbinResults$transformations){
      if(!errorsAsWarnings) stop("Data has been noise removed previously, this would corrupt the data.")
      warning("Data has been noise removed previously, this would corrupt the data.")
    }
    if("Unit variance scaled" %in% mrbinResults$transformations){
      if(!errorsAsWarnings) stop("Data has been unit variance scaled previously, this would corrupt the data.")
      warning("Data has been unit variance scaled previously, this would corrupt the data.")
    }
    if("Atnv transformed"%in% mrbinResults$transformations){
      if(!errorsAsWarnings) stop("Data has been atnv transformed previously, this would corrupt the data.")
      warning("Data has been atnv transformed previously, this would corrupt the data.")
    }
    if("PQN scaled" %in% mrbinResults$transformations){
      if(!errorsAsWarnings) stop("Data has been PQN transformed previously, this would corrupt the data.")
      warning("Data has been PQN transformed previously, this would corrupt the data.")
    }
    NMRfolders<-mrbinResults$parameters$NMRfolders
    factors<-mrbinResults$metadata$factors
 } else {
   NMRfolders<-mrbin.env$mrbin$parameters$NMRfolders
   factors<-mrbin.env$mrbin$metadata$factors
 }
   listTMP<-spectra
   if(is.null(spectra)){
    listTMP<-utils::select.list(NMRfolders,preselect=NULL,
      multiple=TRUE,title ="Remove spectra, cancel to keep all",graphics=TRUE)
   }
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
         if(length(factors)==length(NMRfolders)){
           factors<-factors[-which(NMRfolders%in%listTMP),drop=FALSE]
         }
         NMRfolders<-NMRfolders[-which(NMRfolders%in%listTMP),drop=FALSE]
         if(!is.null(mrbinResults)){
          mrbinResults<-editmrbin(mrbinObject=mrbinResults,
           functionName="mrbin::removeSpectrum",
           versionNumber=as.character(utils::packageVersion("mrbin")),
           bins=mrbinResults$bins[-which(mrbinResults$parameters$NMRfolders%in%listTMP),,drop=FALSE],
           parameters=list(NMRfolders=NMRfolders),
           metadata=list(factors=factors),
           transformations=transformations,verbose=verbose)
         } else {
             mrbin.env$mrbin$parameters$NMRfolders<-NMRfolders
             mrbin.env$mrbin$metadata$factors<-factors
         }
        }
    }
   if(!is.null(mrbinResults)) invisible(mrbinResults)
}

#' A function for setting group members.
#'
#' This function lets the user pick samples from a list to assign them to
#' groups. This function is meant only for use within the mrbin function.
#' @return An invisible mrbin object
#' @keywords internal
#' @noRd
#' @examples
#' \donttest{ setFactors() }

setFactors<-function(mrbinResults){
   if(length(mrbinResults$parameters$NMRfolders)>0){
      Factors<-rep("Group 0",length(mrbinResults$parameters$NMRfolders))
      names(Factors)<-mrbinResults$parameters$NMRfolders
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
      mrbinResults<-editmrbin(mrbinResults,functionName="mrbin::setFactors",
           versionNumber=as.character(utils::packageVersion("mrbin")),metadata=list(factors=Factors),
           verbose=FALSE)

   }
   invisible(mrbinResults)
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
      if(mrbin.env$mrbin$parameters$NMRvendor=="Bruker"){
          selectionFolders<-selectBrukerFolders(keep)
      }  else {
          stop(paste("No folder selection function defined for vendor ",mrbin.env$mrbin$parameters$NMRvendor,".\n",sep=""))
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
    mrbin.env$mrbin$parameters$NMRfolders<-NULL
  }
  NMRfoldersTMP<-NULL
  datanameDict<-c("1r","2rr")
  names(datanameDict)<-c("1D","2D")
  datanameTmp<-datanameDict[mrbin.env$mrbin$parameters$dimension]
  singleFolderFlag<-FALSE
  enterFolders<-utils::select.list(c("Browse...",#"Enter parent folder path manually",
                 "Go back"),
                 preselect="Browse...",title="Set NMR parent folder:",graphics=TRUE)
  if(enterFolders==""|length(enterFolders)==0){
       if(is.null(mrbin.env$mrbin$parameters$NMRfolders)) selectionFolders<-"stop"
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
             if(mrbin.env$mrbin$parameters$verbose) message("Looking for NMR data in folder...")
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
               if(!is.null(NMRfoldersTMP)) mrbin.env$mrbin$parameters$NMRfolders<-unique(c(mrbin.env$mrbin$parameters$NMRfolders,NMRfoldersTMP))
               #if(mrbin.env$mrbin$parameters$verbose) message(paste("Adding to list:\n",paste(NMRfoldersTMP,"\n",
               #                                     sep="",collapse="")), appendLF = FALSE)
               addSpectrumTMP<-TRUE
               while(addSpectrumTMP){
                 yesornoPreSelect<-paste("Keep current spectra list (",length(mrbin.env$mrbin$parameters$NMRfolders)," spectra)",sep="")
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
      #if(!is.null(NMRfolders)) mrbin.env$mrbin$parameters$NMRfolders<-unique(c(mrbin.env$mrbin$parameters$NMRfolders,NMRfolders))
    }
  }
  invisible(selectionFolders)
}

#' A function for binning multiple NMR spectra.
#'
#' This function creates bins for each spectrum in mrbin.env$mrbin$parameters$NMRfolders and
#' returns the bins. Should only be run within the mrbin function.
#' @return An invisible mrbin object
#' @keywords internal
#' @noRd
#' @examples
#' \donttest{ binMultiNMR() }

binMultiNMR<-function(){
 if(!is.null(mrbin.env$mrbin$parameters$NMRfolders)){
    mrbin.env$mrbinTMP$binNames<-NULL
    #Open and bin all spectra
    #Before binning first spectrum
    mrbin.env$mrbin$parameters$AcquPars<-list(
      NS=rep(0,length(mrbin.env$mrbin$parameters$NMRfolders)),
      BF1=rep(0,length(mrbin.env$mrbin$parameters$NMRfolders)),
      P1=rep(0,length(mrbin.env$mrbin$parameters$NMRfolders)),
      RG=rep(0,length(mrbin.env$mrbin$parameters$NMRfolders)),
      PULPROG=rep("",length(mrbin.env$mrbin$parameters$NMRfolders)),
      SOLVENT=rep("",length(mrbin.env$mrbin$parameters$NMRfolders)))
    mrbin.env$mrbin$parameters$noise_level_Raw<-rep(NA,length(mrbin.env$mrbin$parameters$NMRfolders))
    mrbin.env$mrbin$parameters$noise_level<-matrix(rep(NA,length(mrbin.env$mrbin$parameters$NMRfolders)*
      nrow(mrbin.env$mrbin$parameters$binRegions)),ncol=nrow(mrbin.env$mrbin$parameters$binRegions))
    mrbin.env$mrbinTMP$meanNumberOfPointsPerBin<-matrix(rep(NA,length(mrbin.env$mrbin$parameters$NMRfolders)*
      nrow(mrbin.env$mrbin$parameters$binRegions)),ncol=nrow(mrbin.env$mrbin$parameters$binRegions))
    mrbin.env$mrbinTMP$binTMP<-NULL
    mrbin.env$mrbinTMP$currentFolder<-mrbin.env$mrbin$parameters$NMRfolders[1]
    transformations<-NULL
    if(mrbin.env$mrbin$parameters$referenceScaling=="Yes"){
      transformations<-"Reference scaled"
    } else {
      transformations<-"Not scaled"
    }
    useParallel<-FALSE
    if(mrbin.env$mrbin$parameters$tryParallel){
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
      try(
        parallel::clusterExport(cluster, c(
          "readNMR","readBruker","referenceScaling","removeSolvent2",
          "removeAreas2","binSingleNMR","calculateNoise",
          "checkBaseline"))
      ,silent=TRUE)
      try(
        binData<-parallel::parLapply(cluster,
          mrbin.env$mrbin$parameters$NMRfolders,binMultiNMR2,
          dimension=mrbin.env$mrbin$parameters$dimension,
          binRegions=mrbin.env$mrbin$parameters$binRegions,
          referenceScaling=mrbin.env$mrbin$parameters$referenceScaling,
          removeSolvent=mrbin.env$mrbin$parameters$removeSolvent,
          removeAreas=mrbin.env$mrbin$parameters$removeAreas,
          reference1D=mrbin.env$mrbin$parameters$reference1D,
          reference2D=mrbin.env$mrbin$parameters$reference2D,
          solventRegion=mrbin.env$mrbin$parameters$solventRegion,
          removeAreaList=mrbin.env$mrbin$parameters$removeAreaList,
          NMRvendor=mrbin.env$mrbin$parameters$NMRvendor,
          noiseRange1d=mrbin.env$mrbin$parameters$noiseRange1d,
          noiseRange2d=mrbin.env$mrbin$parameters$noiseRange2d,
          binMethod=mrbin.env$mrbin$parameters$binMethod,
          useAsNames=mrbin.env$mrbin$parameters$useAsNames,
          useMeanIntensityForBins=mrbin.env$mrbin$parameters$useMeanIntensityForBins
          #,readAcqus=TRUE
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
        mrbin.env$mrbin$parameters$NMRfolders,binMultiNMR2,
        dimension=mrbin.env$mrbin$parameters$dimension,
        binRegions=mrbin.env$mrbin$parameters$binRegions,
        referenceScaling=mrbin.env$mrbin$parameters$referenceScaling,
        removeSolvent=mrbin.env$mrbin$parameters$removeSolvent,
        removeAreas=mrbin.env$mrbin$parameters$removeAreas,
        reference1D=mrbin.env$mrbin$parameters$reference1D,
        reference2D=mrbin.env$mrbin$parameters$reference2D,
        solventRegion=mrbin.env$mrbin$parameters$solventRegion,
        removeAreaList=mrbin.env$mrbin$parameters$removeAreaList,
        NMRvendor=mrbin.env$mrbin$parameters$NMRvendor,
        noiseRange1d=mrbin.env$mrbin$parameters$noiseRange1d,
        noiseRange2d=mrbin.env$mrbin$parameters$noiseRange2d,
        binMethod=mrbin.env$mrbin$parameters$binMethod,
        useAsNames=mrbin.env$mrbin$parameters$useAsNames,
        useMeanIntensityForBins=mrbin.env$mrbin$parameters$useMeanIntensityForBins
        #,readAcqus=TRUE
        )
    }
    binsRaw<-matrix(rep(0,nrow(mrbin.env$mrbin$parameters$binRegions)*
                                      length(mrbin.env$mrbin$parameters$NMRfolders)),
                                      nrow=length(mrbin.env$mrbin$parameters$NMRfolders))
    currentSpectrumNameTMP<-paste("TemporaryRowName_",1:length(mrbin.env$mrbin$parameters$NMRfolders),sep="")
    for(ibinData in 1:length(binData)){
      if(!is.null(binData[[ibinData]]$warningMessage)){
        warning(binData[[ibinData]]$warningMessage)
        #save warning messages to parameters
         mrbin.env$mrbin$parameters$warningMessages<-
            c(mrbin.env$mrbin$parameters$warningMessages,
              binData[[ibinData]]$warningMessage)
      }
      #save acquisition parameters
      mrbin.env$mrbin$parameters$AcquPars$RG[ibinData]<-binData[[ibinData]]$AcquPars$RG
      mrbin.env$mrbin$parameters$AcquPars$NS[ibinData]<-binData[[ibinData]]$AcquPars$NS
      mrbin.env$mrbin$parameters$AcquPars$BF1[ibinData]<-binData[[ibinData]]$AcquPars$BF1
      mrbin.env$mrbin$parameters$AcquPars$PULPROG[ibinData]<-binData[[ibinData]]$AcquPars$PULPROG
      mrbin.env$mrbin$parameters$AcquPars$SOLVENT[ibinData]<-binData[[ibinData]]$AcquPars$SOLVENT
      mrbin.env$mrbin$parameters$AcquPars$P1[ibinData]<-binData[[ibinData]]$AcquPars$P1

      binsRaw[ibinData,]<-binData[[ibinData]]$binTMP
      mrbin.env$mrbinTMP$meanNumberOfPointsPerBin[ibinData,]<-
        binData[[ibinData]]$meanNumberOfPointsPerBin_TMP
      mrbin.env$mrbin$parameters$noise_level_Raw[ibinData]<-
        binData[[ibinData]]$noise_level_Raw_TMP#noise/sample before reference scaling
      mrbin.env$mrbin$parameters$noise_level[ibinData,]<-
        binData[[ibinData]]$noise_level_TMP
      mrbin.env$mrbin$parameters$noise_level_adjusted[ibinData]<-
        median(binData[[ibinData]]$noise_level_TMP)
      currentSpectrumNameTMP[ibinData]<-binData[[ibinData]]$currentSpectrumName
    }
    if(length(unique(mrbin.env$mrbin$parameters$AcquPars$NS))>1){
      warningTMP<-paste("NS mismatch: ",
        paste(unique(mrbin.env$mrbin$parameters$AcquPars$NS),sep=", ",collapse=", ")
        ,". Number of scans differs between spectra. Data could be inconsistent.",sep="")
      warning(warningTMP)#save warning message to parameters
      mrbin.env$mrbin$parameters$warningMessages<-
            c(mrbin.env$mrbin$parameters$warningMessages,warningTMP)
    }
    if(length(unique(mrbin.env$mrbin$parameters$AcquPars$BF1))>1){
      warningTMP<-paste("BF1 mismatch: ",
        paste(unique(mrbin.env$mrbin$parameters$AcquPars$BF1),sep=", ",collapse=", ")
        ,". Magnets differ between spectra. Data could be inconsistent.",sep="")
      warning(warningTMP)#save warning message to parameters
      mrbin.env$mrbin$parameters$warningMessages<-
            c(mrbin.env$mrbin$parameters$warningMessages,warningTMP)
    }
    if(length(unique(mrbin.env$mrbin$parameters$AcquPars$PULPROG))>1){
      warningTMP<-paste("PULPROG mismatch: ",
        paste(unique(mrbin.env$mrbin$parameters$AcquPars$PULPROG),sep=", ",collapse=", ")
        ,". Pulse programs differ between spectra. Data could be inconsistent.",sep="")
      warning(warningTMP)#save warning message to parameters
      mrbin.env$mrbin$parameters$warningMessages<-
            c(mrbin.env$mrbin$parameters$warningMessages,warningTMP)
    }
    if(length(unique(mrbin.env$mrbin$parameters$AcquPars$SOLVENT))>1){
      warningTMP<-paste("SOLVENT mismatch: ",
        paste(unique(mrbin.env$mrbin$parameters$AcquPars$SOLVENT),sep=", ",collapse=", ")
        ,". Solvents differ between spectra. Data could be inconsistent.",sep="")
      warning(warningTMP)#save warning message to parameters
      mrbin.env$mrbin$parameters$warningMessages<-
            c(mrbin.env$mrbin$parameters$warningMessages,warningTMP)
    }
    #i_currentSpectrumNameTMP<-1
    duplicatedFlag<-FALSE
    while(sum(duplicated(currentSpectrumNameTMP))>0){
      duplicatedFlag<-TRUE
      duplicateListTMP<-which(currentSpectrumNameTMP==currentSpectrumNameTMP[which(duplicated(currentSpectrumNameTMP))[1]])
      for(i_currentSpectrumNameTMP in 1:length(duplicateListTMP)){
       currentSpectrumNameTMP[duplicateListTMP[i_currentSpectrumNameTMP]]<-
            paste(currentSpectrumNameTMP[duplicateListTMP[i_currentSpectrumNameTMP]],".",
            i_currentSpectrumNameTMP,sep="")
       #i_currentSpectrumNameTMP<-i_currentSpectrumNameTMP+1
      }
    }
    if(duplicatedFlag){
      warningDuplicateTMP<-paste(
        "Renamed duplicate spectrum titles. Please use a different naming method.",
        sep="")
      warning(warningDuplicateTMP)
      #save warning message to parameters
      mrbin.env$mrbin$parameters$warningMessages<-
            c(mrbin.env$mrbin$parameters$warningMessages,
              warningDuplicateTMP)
    }
    rownames(binsRaw)<-currentSpectrumNameTMP
    colnames(binsRaw)<-names(binData[[1]]$binTMP)
    binsRaw<-createBinNames(binsRaw)
    #Store bins in package environment for noise previews etc.
    mrbin.env$mrbin$bins<-binsRaw
    #create object of type "mrbin"
    mrbinResults<-createmrbin()
    mrbinResults<-editmrbin(mrbinObject=mrbinResults,
      functionName="mrbin::binMultiNMR",
      versionNumber=as.character(utils::packageVersion("mrbin")),
      bins=binsRaw, parameters=mrbin.env$mrbin$parameters,
      metadata=mrbin.env$mrbin$metadata,
      transformations=transformations,
      verbose=FALSE)
    invisible(mrbinResults)
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
   if(mrbin.env$mrbin$parameters$binMethod=="Rectangular bins"){
       mrbin.env$mrbin$parameters$binRegions<-matrix(ncol=4,
                                        nrow=mrbin.env$mrbinTMP$nbins,
                                        dimnames=list(NULL,c("left","right","top","bottom")))
       if(mrbin.env$mrbin$parameters$dimension=="1D"){
          decimalDigits<-max(nchar(strsplit(as.character(mrbin.env$mrbin$parameters$binRegion[1]),"[.]")[[1]][2]),
                             nchar(strsplit(as.character(mrbin.env$mrbin$parameters$binwidth1D),"[.]")[[1]][2]),
                             0,na.rm=TRUE)+4
          mrbin.env$mrbin$parameters$binRegions[,1]<-round(mrbin.env$mrbin$parameters$binRegion[1]-((1:mrbin.env$mrbinTMP$nbins)-1)*mrbin.env$mrbin$parameters$binwidth1D,decimalDigits)
          mrbin.env$mrbin$parameters$binRegions[,2]<-round(mrbin.env$mrbin$parameters$binRegion[1]-((1:mrbin.env$mrbinTMP$nbins))*mrbin.env$mrbin$parameters$binwidth1D,decimalDigits)
       }
       if(mrbin.env$mrbin$parameters$dimension=="2D"){
          decimalDigits1<-max(nchar(strsplit(as.character(mrbin.env$mrbin$parameters$binRegion[1]),"[.]")[[1]][2]),
                              nchar(strsplit(as.character(mrbin.env$mrbin$parameters$binwidth2D),"[.]")[[1]][2]),
                              0,na.rm=TRUE)+2
          decimalDigits2<-max(nchar(strsplit(as.character(mrbin.env$mrbin$parameters$binRegion[4]),"[.]")[[1]][2]),
                              nchar(strsplit(as.character(mrbin.env$mrbin$parameters$binheight),"[.]")[[1]][2]),
                              0,na.rm=TRUE)+2
          mrbin.env$mrbin$parameters$binRegions[,1]<-round(rep(mrbin.env$mrbin$parameters$binRegion[1]-((1:mrbin.env$mrbinTMP$nbins2)-1)*mrbin.env$mrbin$parameters$binwidth2D,
                                             mrbin.env$mrbinTMP$nbins1),decimalDigits1)
          mrbin.env$mrbin$parameters$binRegions[,2]<-round(rep(mrbin.env$mrbin$parameters$binRegion[1]-((1:mrbin.env$mrbinTMP$nbins2))*mrbin.env$mrbin$parameters$binwidth2D,
                                             mrbin.env$mrbinTMP$nbins1),decimalDigits1)
          mrbin.env$mrbin$parameters$binRegions[,3]<-round(sort(rep(mrbin.env$mrbin$parameters$binRegion[4]-((1:mrbin.env$mrbinTMP$nbins1))*mrbin.env$mrbin$parameters$binheight,
                                             mrbin.env$mrbinTMP$nbins2),decreasing=TRUE),decimalDigits2)
          mrbin.env$mrbin$parameters$binRegions[,4]<-round(sort(rep(mrbin.env$mrbin$parameters$binRegion[4]-((1:mrbin.env$mrbinTMP$nbins1)-1)*mrbin.env$mrbin$parameters$binheight,
                                             mrbin.env$mrbinTMP$nbins2),decreasing=TRUE),decimalDigits2)
       }
  }
  if(mrbin.env$mrbin$parameters$binMethod=="Custom bin list"){
    mrbin.env$mrbin$parameters$binRegions<-mrbin.env$mrbin$parameters$specialBinList
    if(!is.matrix(mrbin.env$mrbin$parameters$binRegions)) mrbin.env$mrbin$parameters$binRegions<-matrix(
      mrbin.env$mrbin$parameters$binRegions,ncol=4)
  }
}

#' A function for creating bin titles.
#'
#' This function creates titles for the bins to represent their ppm range. This function is
#' meant only for use within the mrbin function.
#' @param binsRaw A matrix of binned data
#' @return An invisible matrix of binned data
#' @keywords internal
#' @noRd
#' @examples
#' \donttest{ createBinNames() }

createBinNames<-function(binsRaw){
   rownamesSepSignTMP<-""
   if(!is.null(rownames(mrbin.env$mrbin$parameters$binRegions))){
     if(sum(rownames(mrbin.env$mrbin$parameters$binRegions)=="")<nrow(mrbin.env$mrbin$parameters$binRegions)){#For custom bin lists
        rownamesSepSignTMP<-";"
     }
   }
   if(mrbin.env$mrbin$parameters$dimension=="1D"){
     if(nrow(mrbin.env$mrbin$parameters$binRegions)>1){
          namesTMP<-apply(mrbin.env$mrbin$parameters$binRegions[,1:2],1,paste,collapse=",")
          rownames(mrbin.env$mrbin$parameters$binRegions)<-paste(rownames(mrbin.env$mrbin$parameters$binRegions),namesTMP,
                                                         sep=rownamesSepSignTMP)
     } else {
          namesTMP<-paste(mrbin.env$mrbin$parameters$binRegions[,1:2],collapse=",")
          rownames(mrbin.env$mrbin$parameters$binRegions)<-paste(rownames(mrbin.env$mrbin$parameters$binRegions),
                                                         namesTMP,
                                                         sep=rownamesSepSignTMP)
     }
   }
   if(mrbin.env$mrbin$parameters$dimension=="2D"){
     if(nrow(mrbin.env$mrbin$parameters$binRegions)>1){
          namesTMP<-apply(mrbin.env$mrbin$parameters$binRegions,1,paste,collapse=",")
          rownames(mrbin.env$mrbin$parameters$binRegions)<-paste(rownames(mrbin.env$mrbin$parameters$binRegions),
                                                         namesTMP,
                                                         sep=rownamesSepSignTMP)
     } else {
          namesTMP<-paste(mrbin.env$mrbin$parameters$binRegions,collapse=",")
          rownames(mrbin.env$mrbin$parameters$binRegions)<-paste(rownames(mrbin.env$mrbin$parameters$binRegions),
                                                         namesTMP,
                                                         sep=rownamesSepSignTMP)
     }
   }
   mrbin.env$mrbinTMP$binNames<-rownames(mrbin.env$mrbin$parameters$binRegions)
   colnames(binsRaw)<-rownames(mrbin.env$mrbin$parameters$binRegions)
   invisible(binsRaw)
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
   if(mrbin.env$mrbin$parameters$binMethod=="Custom bin list"){
          mrbin.env$mrbinTMP$nbins<-nrow(mrbin.env$mrbin$parameters$specialBinList)
   }
   if(mrbin.env$mrbin$parameters$binMethod=="Rectangular bins"){
       if(mrbin.env$mrbin$parameters$dimension=="2D"){
          mrbin.env$mrbinTMP$nbins2<-ceiling((mrbin.env$mrbin$parameters$binRegion[1]-mrbin.env$mrbin$parameters$binRegion[2])/mrbin.env$mrbin$parameters$binwidth2D)
          mrbin.env$mrbinTMP$nbins1<-ceiling((mrbin.env$mrbin$parameters$binRegion[4]-mrbin.env$mrbin$parameters$binRegion[3])/mrbin.env$mrbin$parameters$binheight)
          mrbin.env$mrbinTMP$nbins<-mrbin.env$mrbinTMP$nbins2*mrbin.env$mrbinTMP$nbins1
      }
      if(mrbin.env$mrbin$parameters$dimension=="1D"){
          mrbin.env$mrbinTMP$nbins<-ceiling((mrbin.env$mrbin$parameters$binRegion[1]-mrbin.env$mrbin$parameters$binRegion[2])/mrbin.env$mrbin$parameters$binwidth1D)
      }
   }
}



#' A function for reading NMR spectra.
#'
#' This function picks the correct NMR reading function, based on vendor. This function is
#' meant only for use within the mrbin function.
#' @param onlyTitles Read only spectrum titles, but no data. Defaults to FALSE
#' @param dimension If not provided, this is taken from the package environment
#' @param NMRvendor If not provided, this is taken from the package environment
#' @param useAsNames If not provided, this is taken from the package environment
#' @param folder If not provided, this is taken from the package environment
#' @return {none}
#' @keywords internal
#' @noRd
#' @examples
#' \donttest{ readNMR2() }

readNMR2<-function(onlyTitles=FALSE,dimension=NULL,NMRvendor=NULL,
  useAsNames=NULL,folder=NULL#,readAcqus=FALSE
  ){#Read NMR spectral data
   if(is.null(dimension)) dimension<-mrbin.env$mrbin$parameters$dimension
   if(is.null(NMRvendor)) NMRvendor<-mrbin.env$mrbin$parameters$NMRvendor
   if(is.null(useAsNames)) useAsNames<-mrbin.env$mrbin$parameters$useAsNames
   if(is.null(folder)) folder<-mrbin.env$mrbinTMP$currentFolder

   NMRdataList<-readNMR(onlyTitles=onlyTitles,folder=folder,
           dimension=dimension,NMRvendor=NMRvendor,
           useAsNames=useAsNames#,readAcqus=readAcqus
           )
   if(!onlyTitles){
     mrbin.env$mrbinTMP$currentSpectrum<-NMRdataList$currentSpectrum
     mrbin.env$mrbinTMP$currentSpectrumOriginal<-NMRdataList$currentSpectrum
     #<-NMRdataList$AcquPars
   }
   mrbin.env$mrbinTMP$currentSpectrumTitle<-NMRdataList$currentSpectrumTitle
   mrbin.env$mrbinTMP$currentSpectrumFolderName<-NMRdataList$currentSpectrumFolderName
   mrbin.env$mrbinTMP$currentSpectrumEXPNO<-NMRdataList$currentSpectrumEXPNO
   mrbin.env$mrbinTMP$currentSpectrumFolderName_EXPNO<-NMRdataList$currentSpectrumFolderName_EXPNO
   mrbin.env$mrbinTMP$currentSpectrumName<-NMRdataList$currentSpectrumName
}

#' A function for adding NMR spectra to the plot list.
#'
#' This function adds a spectrum to the plot list.
#' @param folder Defines the exact NMR data folder. If NULL, mrbin parameter set is used
#' @param dimension Defines the data dimension, "1D" or "2D". Only used if not NULL
#' @param NMRvendor Defines the NMR manufacturer, default is "Bruker"
#' @param useAsNames How should sample names be generated
#' @param add Add spectra to existing list, or replace existing spectra. Default is TRUE
#' @return {none}
#' @export
#' @examples
#' \donttest{ addToPlot() }

addToPlot<-function(folder=NULL,
           dimension="1D",
           NMRvendor="Bruker",
           useAsNames="Folder names",
           add=TRUE){#Read NMR spectral data
  if(is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){#if no spectrum has ever been loaded
    readNMR2(folder=folder,
           dimension=dimension,
           NMRvendor=NMRvendor,
           useAsNames=useAsNames)
  } else {
    NMRdataList<-readNMR(folder=folder,
           dimension=dimension,
           NMRvendor=NMRvendor,
           useAsNames=useAsNames)
    if(dimension=="1D"){
     if(!add){
       mrbin.env$mrbinTMP$additionalPlots1D<-NULL
       mrbin.env$mrbinTMP$additionalPlots1DMetadata<-NULL
     }
     if( is.null(mrbin.env$mrbinTMP$additionalPlots1D)){
       mrbin.env$mrbinTMP$additionalPlots1D<-list(NMRdataList$currentSpectrum)
     } else {
       mrbin.env$mrbinTMP$additionalPlots1D[[
         length(mrbin.env$mrbinTMP$additionalPlots1D)+1]]<-NMRdataList$currentSpectrum
     }
     mrbin.env$mrbinTMP$additionalPlots1DMetadata<-rbind(
       mrbin.env$mrbinTMP$additionalPlots1DMetadata,
       c(folder,
         NMRdataList$currentSpectrumTitle,
         NMRdataList$currentSpectrumFolderName,
         NMRdataList$currentSpectrumEXPNO,
         NMRdataList$currentSpectrumFolderName_EXPNO,
         NMRdataList$currentSpectrumName,
         1,#This is a scaling factor
         setdiff(unique(c(grDevices::colors()[c(552,26,547,81,498)],
           grDevices::colors())),c("white","black",
           mrbin.env$mrbinTMP$additionalPlots1DMetadata[,8]))[1],#new color
         0,#offset
         FALSE#hide this spectrum
         ))
   }
   if(dimension=="2D"){
     if(!add){
       mrbin.env$mrbinTMP$additionalPlots2D<-NULL
       mrbin.env$mrbinTMP$additionalPlots2DMetadata<-NULL
     }
     if( is.null(mrbin.env$mrbinTMP$additionalPlots2D)){
       mrbin.env$mrbinTMP$additionalPlots2D<-list(NMRdataList$currentSpectrum)
     } else {
       mrbin.env$mrbinTMP$additionalPlots2D[[
         length(mrbin.env$mrbinTMP$additionalPlots2D)+1]]<-NMRdataList$currentSpectrum
     }
     mrbin.env$mrbinTMP$additionalPlots2DMetadata<-rbind(
       mrbin.env$mrbinTMP$additionalPlots2DMetadata,
       c(folder,
         NMRdataList$currentSpectrumTitle,
         NMRdataList$currentSpectrumFolderName,
         NMRdataList$currentSpectrumEXPNO,
         NMRdataList$currentSpectrumFolderName_EXPNO,
         NMRdataList$currentSpectrumName,
         1,#This is a scaling factor
         setdiff(unique(c(grDevices::colors()[c(552,26,547,81,498)],
           grDevices::colors())),c("white","black",
           mrbin.env$mrbinTMP$additionalPlots2DMetadata[,8]))[1],#new color
         0,#offset
         FALSE#hide this spectrum
         ))
   }
  }
}


#' A function for removing NMR spectra from the plot list.
#'
#' This function removes a spectrum from the plot list.
#' @param folder Defines the exact NMR data folder.
#' @param dimension Defines the data dimension, "1D" or "2D".
#' @return {none}
#' @export
#' @examples
#' \donttest{ removeFromPlot() }

removeFromPlot<-function(folder=NULL,
           dimension="1D"){#
   if(dimension=="1D"){
     if(!is.null(mrbin.env$mrbinTMP$additionalPlots1D)){
       if(folder%in%mrbin.env$mrbinTMP$additionalPlots1DMetadata[,5]){
         if(nrow(mrbin.env$mrbinTMP$additionalPlots1DMetadata)==1){
            mrbin.env$mrbinTMP$additionalPlots1DMetadata<-NULL
            mrbin.env$mrbinTMP$additionalPlots1D<-NULL
         } else {
           deleteTMP<-which(mrbin.env$mrbinTMP$additionalPlots1DMetadata[,5]==folder)
           mrbin.env$mrbinTMP$additionalPlots1D[[deleteTMP]]<-NULL
           mrbin.env$mrbinTMP$additionalPlots1DMetadata<-
             mrbin.env$mrbinTMP$additionalPlots1DMetadata[-deleteTMP,,drop=FALSE]
         }
       } else {message("folder name not found. ")}
     }
   }
   if(dimension=="2D"){
     if(!is.null(mrbin.env$mrbinTMP$additionalPlots2D)){
       if(folder%in%mrbin.env$mrbinTMP$additionalPlots2DMetadata[,5]){
         if(nrow(mrbin.env$mrbinTMP$additionalPlots2DMetadata)==1){
            mrbin.env$mrbinTMP$additionalPlots2DMetadata<-NULL
            mrbin.env$mrbinTMP$additionalPlots2D<-NULL
         } else {
           deleteTMP<-which(mrbin.env$mrbinTMP$additionalPlots2DMetadata[,5]==folder)
           mrbin.env$mrbinTMP$additionalPlots2D[[deleteTMP]]<-NULL
           mrbin.env$mrbinTMP$additionalPlots2DMetadata<-
             mrbin.env$mrbinTMP$additionalPlots2DMetadata[-deleteTMP,,drop=FALSE]
         }
       }
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
#' @keywords internal
#' @noRd
#' @examples
#' \donttest{ sumBins() }

sumBins<-function(){#sum up regions with shifting peaks and remove remaining bins
 #if(!is.null(mrbin.env$mrbin$bins)){
  if(nrow(mrbin.env$mrbin$parameters$sumBinList)>0&nrow(mrbin.env$mrbin$parameters$binRegions)>1){
    for(i in 1:nrow(mrbin.env$mrbin$parameters$sumBinList)){
       limits<-mrbin.env$mrbin$parameters$sumBinList[i,]
       if(mrbin.env$mrbin$parameters$dimension=="1D"){
          #Find partially overlapping bins
          TMP_left<-mrbin.env$mrbin$parameters$binRegions[,2]<limits[1]&mrbin.env$mrbin$parameters$binRegions[,1]>=limits[1]
          mrbin.env$mrbin$parameters$binRegions[TMP_left,2]<-limits[1]
          TMP_right<-mrbin.env$mrbin$parameters$binRegions[,1]>limits[2]&mrbin.env$mrbin$parameters$binRegions[,2]<=limits[2]
          mrbin.env$mrbin$parameters$binRegions[TMP_right,1]<-limits[2]
          #Find completely overlapping bins
          TMP<-mrbin.env$mrbin$parameters$binRegions[,2]>limits[2]&mrbin.env$mrbin$parameters$binRegions[,1]<limits[1]
          if(sum(TMP)>0){
              #i_TMP<-quantile(x=1:sum(TMP), probs = .5,type=3)#define "middle" bin. This one will be kept
              #i_TMP2<-which(TMP)[i_TMP]
              #mrbin.env$mrbin$parameters$binRegions[i_TMP2,1]<-limits[1]
              #mrbin.env$mrbin$parameters$binRegions[i_TMP2,2]<-limits[2]
              #if(sum(TMP)>1)
              mrbin.env$mrbin$parameters$binRegions<-mrbin.env$mrbin$parameters$binRegions[!TMP,,drop=FALSE]
          } #else {
              #mrbin.env$mrbin$parameters$binRegions<-rbind(c(0,0,0,0),mrbin.env$mrbin$parameters$binRegions)
              #i_TMP2<-1#nrow(mrbin.env$mrbin$parameters$binRegions)
              #mrbin.env$mrbin$parameters$binRegions[i_TMP2,1]<-limits[1]
              #mrbin.env$mrbin$parameters$binRegions[i_TMP2,2]<-limits[2]
          #}
          #if(!is.matrix(mrbin.env$mrbin$parameters$binRegions)) mrbin.env$mrbin$parameters$binRegions<-matrix(mrbin.env$mrbin$parameters$binRegions,ncol=4)

       } else {#2D limits=c(4.04,4.08,58,60)
           #NMRdataNames<-cbind(apply(mrbin.env$mrbin$parameters$binRegions[,3:4],1,mean),apply(mrbin.env$mrbin$parameters$binRegions[,1:2],1,mean))
           #Find partially overlapping bins ("corners" are an issue here)
           TMP_left<-mrbin.env$mrbin$parameters$binRegions[,2]<limits[1]&mrbin.env$mrbin$parameters$binRegions[,1]>limits[1]&
                     mrbin.env$mrbin$parameters$binRegions[,3]>=limits[3]&mrbin.env$mrbin$parameters$binRegions[,4]<=limits[4]
           mrbin.env$mrbin$parameters$binRegions[TMP_left,2]<-limits[1]
           TMP_right<-mrbin.env$mrbin$parameters$binRegions[,1]>limits[2]&mrbin.env$mrbin$parameters$binRegions[,2]<limits[2]&
                     mrbin.env$mrbin$parameters$binRegions[,3]>=limits[3]&mrbin.env$mrbin$parameters$binRegions[,4]<=limits[4]
           mrbin.env$mrbin$parameters$binRegions[TMP_right,1]<-limits[2]
           TMP_top<-mrbin.env$mrbin$parameters$binRegions[,1]<=limits[1]&mrbin.env$mrbin$parameters$binRegions[,2]>=limits[2]&
                     mrbin.env$mrbin$parameters$binRegions[,3]<limits[3]&mrbin.env$mrbin$parameters$binRegions[,4]>limits[3]
           mrbin.env$mrbin$parameters$binRegions[TMP_top,3]<-limits[3]
           TMP_bottom<-mrbin.env$mrbin$parameters$binRegions[,1]<=limits[1]&mrbin.env$mrbin$parameters$binRegions[,2]>=limits[2]&
                     mrbin.env$mrbin$parameters$binRegions[,4]>limits[4]&mrbin.env$mrbin$parameters$binRegions[,3]<limits[4]
           mrbin.env$mrbin$parameters$binRegions[TMP_bottom,3]<-limits[4]
           #Find completely overlapping bins
           TMP<-mrbin.env$mrbin$parameters$binRegions[,2]>limits[2]&mrbin.env$mrbin$parameters$binRegions[,1]<limits[1]&
             mrbin.env$mrbin$parameters$binRegions[,3]>limits[3]& mrbin.env$mrbin$parameters$binRegions[,4]<limits[4]
           if(sum(TMP)>0){
              #i_TMP<-quantile(x=1:sum(TMP), probs = .5,type=3)#define "middle" bin. This one will be kept
              #mrbin.env$mrbin$bins[,which(TMP)[i_TMP]]<-apply(mrbin.env$mrbin$bins[,which(TMP)],1,sum)
              #if(nrow(mrbin.env$mrbin$bins)==1){
              #  rownamesTMP<-rownames(mrbin.env$mrbin$bins)
              #  colnamesTMP<-colnames(mrbin.env$mrbin$bins)[-which(TMP)[-i_TMP]]
              #  mrbin.env$mrbin$bins<-matrix(mrbin.env$mrbin$bins[,-which(TMP)[-i_TMP]],nrow=1)
              #  rownames(mrbin.env$mrbin$bins)<-rownamesTMP
              #  colnames(mrbin.env$mrbin$bins)<-colnamesTMP
              #} else {
              #  mrbin.env$mrbin$bins<-mrbin.env$mrbin$bins[,-which(TMP)[-i_TMP]]
              #}
              #mrbin.env$mrbin$parameters$binRegions[which(TMP)[i_TMP],]<-limits
              #if(sum(TMP)>1)
              mrbin.env$mrbin$parameters$binRegions<-mrbin.env$mrbin$parameters$binRegions[!TMP,,drop=FALSE]
           } #else {
           #   mrbin.env$mrbin$parameters$binRegions<-rbind(mrbin.env$mrbin$parameters$binRegions,c(0,0,0,0))
           #   i_TMP2<-nrow(mrbin.env$mrbin$parameters$binRegions)
           #   mrbin.env$mrbin$parameters$binRegions[i_TMP2,]<-limits
           #}
          # if(!is.matrix(mrbin.env$mrbin$parameters$binRegions)) mrbin.env$mrbin$parameters$binRegions<-matrix(mrbin.env$mrbin$parameters$binRegions,ncol=4)
       }
       #New bins are added on top of list to ensure they see all data points before they are set to NA
       mrbin.env$mrbin$parameters$binRegions<-rbind(c(0,0,0,0),mrbin.env$mrbin$parameters$binRegions)
       mrbin.env$mrbin$parameters$binRegions[1,]<-limits
    }
  }
  mrbin.env$mrbin$parameters$numberOfFeaturesAfterSummingBins<-nrow(mrbin.env$mrbin$parameters$binRegions)
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
   solventTMP<-mrbin.env$mrbin$parameters$binRegions[,2]>mrbin.env$mrbin$parameters$solventRegion[2]&
               mrbin.env$mrbin$parameters$binRegions[,1]<mrbin.env$mrbin$parameters$solventRegion[1]
   if(sum(solventTMP)>0){
      mrbin.env$mrbin$parameters$binRegions<-mrbin.env$mrbin$parameters$binRegions[-which(solventTMP),,drop=FALSE]
      if(!is.matrix(mrbin.env$mrbin$parameters$binRegions)) mrbin.env$mrbin$parameters$binRegions<-matrix(mrbin.env$mrbin$parameters$binRegions,ncol=4)
   }
   mrbin.env$mrbin$parameters$numberOfFeaturesAfterRemovingSolvent<-nrow(mrbin.env$mrbin$parameters$binRegions)
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
  if(nrow(mrbin.env$mrbin$parameters$removeAreaList)>0){
     removeTMP<-NULL
     for(i in 1:nrow(mrbin.env$mrbin$parameters$removeAreaList)){
         limits<-mrbin.env$mrbin$parameters$removeAreaList[i,]
         if(mrbin.env$mrbin$parameters$dimension=="1D"){
             removeTMP2<-mrbin.env$mrbin$parameters$binRegions[,2]>limits[2]&mrbin.env$mrbin$parameters$binRegions[,1]<limits[1]
         }
         if(mrbin.env$mrbin$parameters$dimension=="2D"){
             removeTMP2<-mrbin.env$mrbin$parameters$binRegions[,2]>limits[2]&mrbin.env$mrbin$parameters$binRegions[,1]<limits[1]&
                mrbin.env$mrbin$parameters$binRegions[,3]>limits[3]& mrbin.env$mrbin$parameters$binRegions[,4]<limits[4]
         }
         if(sum(removeTMP2)>0)        removeTMP<-c(removeTMP,which(removeTMP2))
     }
     if(!is.null(removeTMP)){
         removeTMP<-unique(removeTMP)
         mrbin.env$mrbin$parameters$binRegions<-mrbin.env$mrbin$parameters$binRegions[-removeTMP,,drop=FALSE]
         if(!is.matrix(mrbin.env$mrbin$parameters$binRegions)) mrbin.env$mrbin$parameters$binRegions<-matrix(mrbin.env$mrbin$parameters$binRegions,ncol=4)
     }
  }
  mrbin.env$mrbin$parameters$numberOfFeaturesAfterRemovingAreas<-nrow(mrbin.env$mrbin$parameters$binRegions)
 #}
}


#' A function for trimming zero-values bins.
#'
#' This function removes zero-values bins. These might be created during removal of
#' solvent and additional areas, or at the edges of the spectrum.
#' @param mrbinResults An mrbin object
#' @return An invisible mrbin object
#' @export
#' @examples
#' mrbinResults<-mrbin(silent=TRUE,setDefault=TRUE,parameters=list(dimension="1D", logTrafo="No",
#'                     binwidth1D=0.05,signal_to_noise1D=50, verbose=TRUE, PCA="No",
#'                       trimZeros="No",tryParallel=FALSE,
#'                     NMRfolders=c(system.file("extdata/1/10/pdata/10",package="mrbin"),
#'                                 system.file("extdata/2/10/pdata/10",package="mrbin"))))
#' mrbinResults<-trimZeros(mrbinResults)

trimZeros<-function(mrbinResults){
  #mrbin.env$mrbin$binsTEST<-mrbin.env$mrbin$bins
  mrbinResultsTMP<-mrbinResults
  if(nrow(mrbinResultsTMP$bins)>1){
    TMP<-apply(mrbinResultsTMP$bins==0,2,sum)/nrow(mrbinResultsTMP$bins)<.75
    if(sum(TMP)>0&sum(TMP)<ncol(mrbinResultsTMP$bins)){
      mrbinResultsTMP$bins<-mrbinResultsTMP$bins[,which(TMP),drop=FALSE]
      mrbinResultsTMP$parameters$binRegions<-mrbinResultsTMP$parameters$binRegions[which(TMP),,drop=FALSE]
    }
  } else {
    #if(sum(sum(mrbin.env$mrbin$bins==0)/nrow(mrbin.env$mrbin$bins)<.75)>0){
    #  mrbin.env$mrbin$bins<-mrbin.env$mrbin$bins[,which(sum(mrbin.env$mrbin$bins==0)/nrow(mrbin.env$mrbin$bins)<.75),drop=FALSE]
    #}
    warning("Too few samples for zero trimming.")
  }
   #mrbin.env$mrbin$parameters$numberOfFeaturesAfterTrimmingZeros<-ncol(mrbinResults$bins)
   mrbinResultsTMP$parameters$numberOfFeaturesAfterTrimmingZeros<-ncol(mrbinResults$bins)
   mrbinResults<-editmrbin(mrbinObject=mrbinResults,functionName="mrbin::trimZeros",
       versionNumber=as.character(utils::packageVersion("mrbin")),
       bins=mrbinResultsTMP$bins, parameters=mrbinResultsTMP$parameters,verbose=FALSE,
       transformations="Zero trimmed")
   invisible(mrbinResults)
}


#' A function for removing bins below noise level.
#'
#' This function checks for each bin (column) whether its level is below the
#' individual noise level times the signal-to-noise ratio. If less than the
#' defined threshold level are above noise*SNR, the whole bin is removed.
#' @param mrbinResults An mrbin object
#' @param verbose Should a summary be printed?
#' @param errorsAsWarnings If TRUE, errors will be turned into warnings. Should be used with care, as errors indicate undocumented changes to the data.
#' @return An invisible mrbin object
#' @export
#' @examples
#' mrbinResults<-mrbin(silent=TRUE,setDefault=TRUE,parameters=list(dimension="1D",
#'                     binwidth1D=0.05,noiseRemoval="No",PQNScaling="No",tryParallel=FALSE,
#'                     fixNegatives="No",logTrafo="No",PCA="No",verbose=TRUE,
#'                     NMRfolders=c(system.file("extdata/1/10/pdata/10",package="mrbin"),
#'                                 system.file("extdata/2/10/pdata/10",package="mrbin"),
#'                                 system.file("extdata/3/10/pdata/10",package="mrbin"))))
#' mrbinResults<-removeNoise(mrbinResults)

removeNoise<-function(mrbinResults,verbose=TRUE,errorsAsWarnings=FALSE){#remove noise peaks
  transformations="Noise removed"
  if(transformations %in% mrbinResults$transformations){
   if(!errorsAsWarnings) stop("Data has been noise filtered previously, this could corrupt the data.")
   warning("Data has been noise filtered previously, this could corrupt the data.")
  }
  if("Log transformed" %in% mrbinResults$transformations){
    if(!errorsAsWarnings) stop("Data has been log transformed previously, this would corrupt the data.")
    warning("Data has been log transformed previously, this would corrupt the data.")
  }
  if("PQN scaled" %in% mrbinResults$transformations){
    if(!errorsAsWarnings) stop("Data has been PQN transformed previously, this could corrupt the data.")
    warning("Data has been PQN transformed previously, this could corrupt the data.")
  }
  mrbinResults2<-mrbinResults
  if(ncol(mrbinResults$bins)>1){
    minimumNumber<-max(1,floor(mrbinResults$parameters$noiseThreshold*nrow(mrbinResults$bins)))
    colnames_NMRdata_no_noise<-NULL
    if(mrbinResults$parameters$dimension=="2D"){
         NMRdataNames<-cbind(apply(mrbinResults$parameters$binRegions[,3:4],1,mean),
           apply(mrbinResults$parameters$binRegions[,1:2],1,mean))
         SNR<-mrbinResults$parameters$signal_to_noise2D
    } else {
         SNR<-mrbinResults$parameters$signal_to_noise1D
    }
    for(i in 1:ncol(mrbinResults$bins)){#Keep only bins where at least X spectra are > SNR
          if(sum(mrbinResults$bins[,i]>(mrbinResults$parameters$noise_level[,i]*SNR))>=minimumNumber){
              colnames_NMRdata_no_noise<-c(colnames_NMRdata_no_noise,i)
          }
    }
    if(!is.null(colnames_NMRdata_no_noise)){
        mrbinResults$bins<-mrbinResults$bins[,colnames_NMRdata_no_noise,drop=FALSE]
        mrbinResults$parameters$binRegions<-
          mrbinResults$parameters$binRegions[colnames_NMRdata_no_noise,,drop=FALSE]
        if(!is.matrix(mrbinResults$parameters$binRegions)) mrbinResults$parameters$binRegions<-
           matrix(mrbinResults$parameters$binRegions,ncol=4)
    } else {
        if(!errorsAsWarnings) warning("No bins above noise level. Noise removal stopped.")
    }
  mrbinResults$parameters$numberOfFeaturesAfterNoiseRemoval<-ncol(mrbinResults$bins)
 } else {
   if(!errorsAsWarnings) warning("Too few bins for noise removal. Noise removal stopped.")
 }
 mrbinResults2<-editmrbin(mrbinObject=mrbinResults2,functionName="mrbin::removeNoise",
       versionNumber=as.character(utils::packageVersion("mrbin")),
       bins=mrbinResults$bins, parameters=mrbinResults$parameters,
        transformations=transformations,verbose=verbose)
 invisible(mrbinResults2)
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
      warningMessage<-paste("Baseline distorted (noise region): ",currentSpectrumName,
                  ". Please check phase and baseline. Results may be corrupted.",
                  sep="")
    }
    invisible(warningMessage)
 }
}



#' A function for cropping HSQC spectra.
#'
#' This function crops HSQC spectra to the region along the diagonal to remove
#' uninformative signals. Will work only for 1H-13C HSQC spectra.
#' @return {None}
#' @export
#' @examples
#' resetEnv()
#' Example<-mrbin(silent=TRUE,
#'          parameters=list(dimension="2D",binwidth2D=1,binheight=4,cropHSQC="No",PCA="No",
#'          PQNScaling="No",noiseRemoval="No",removeSolvent="No",verbose=TRUE,tryParallel=FALSE,
#'          NMRfolders=c(system.file("extdata/1/12/pdata/10",package="mrbin"))))
#' cropNMR()

cropNMR<-function(){
if(nrow(mrbin.env$mrbin$parameters$binRegions)>1){
 if(mrbin.env$mrbin$parameters$dimension=="2D"){
   selectedCols<-NULL
   coordTmp2<-cbind(apply(mrbin.env$mrbin$parameters$binRegions[,3:4],1,mean),apply(mrbin.env$mrbin$parameters$binRegions[,1:2],1,mean))
   for(j in 1:nrow(mrbin.env$mrbin$parameters$binRegions)){
      coordTmp<-coordTmp2[j,]#1=C,2=H
      if(((coordTmp[2]-mrbin.env$mrbin$parameters$croptopLeft[2])*(mrbin.env$mrbin$parameters$cropbottomLeft[1]-mrbin.env$mrbin$parameters$croptopLeft[1])-
          (coordTmp[1]-mrbin.env$mrbin$parameters$croptopLeft[1])*(mrbin.env$mrbin$parameters$cropbottomLeft[2]-mrbin.env$mrbin$parameters$croptopLeft[2]))<0&
         ((coordTmp[2]-mrbin.env$mrbin$parameters$croptopRight[2])*(mrbin.env$mrbin$parameters$cropbottomRight[1]-mrbin.env$mrbin$parameters$croptopRight[1])-
          (coordTmp[1]-mrbin.env$mrbin$parameters$croptopRight[1])*(mrbin.env$mrbin$parameters$cropbottomRight[2]-mrbin.env$mrbin$parameters$croptopRight[2]))>0
        ){#along diagonal? outer product
          selectedCols<-c(selectedCols,j)
      }
  }
  mrbin.env$mrbin$parameters$binRegions<-mrbin.env$mrbin$parameters$binRegions[selectedCols,,drop=FALSE]
  if(!is.matrix(mrbin.env$mrbin$parameters$binRegions)) mrbin.env$mrbin$parameters$binRegions<-matrix(mrbin.env$mrbin$parameters$binRegions,ncol=4)
  mrbin.env$mrbin$parameters$numberOfFeaturesAfterCropping<-nrow(mrbin.env$mrbin$parameters$binRegions)
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
#' @param NMRdata A matrix containing NMR data or an mrbin object. Columns=frequencies,rows=samples
#' @param ignoreGlucose A character value ("Yes" or "No")
#' @param dimension A character value ("1D" or "2D")
#' @param ppmNames A character value ("borders" or "mean")
#' @param sugarArea A numeric vector defining the the borders of glucose area
#' @param minimumFeatures A numeric value defining minimum feature number used
#' @param showHist A logical value, default is FALSE
#' @param verbose Should a summary be printed?
#' @param errorsAsWarnings If TRUE, errors will be turned into warnings. Should be used with care, as errors indicate undocumented changes to the data.
#' @return An invisible matrix or mrbin object containing scaled NMR data.
#' @export
#' @examples
#' mrbinResults<-mrbin(silent=TRUE,setDefault=TRUE,parameters=list(dimension="1D",
#'                     binwidth1D=0.05,PQNScaling="No",PCA="No",tryParallel=FALSE,logTrafo="No",
#'                     NMRfolders=c(system.file("extdata/1/10/pdata/10",package="mrbin"),
#'                                 system.file("extdata/2/10/pdata/10",package="mrbin"),
#'                                 system.file("extdata/3/10/pdata/10",package="mrbin"))))
#' mrbinResults<-PQNScaling(mrbinResults)

PQNScaling<-function(NMRdata,ignoreGlucose="Yes",dimension="1D",
                     ppmNames="borders",sugarArea=c(5.4,3.35,72,100),
                     minimumFeatures=40,showHist=FALSE,verbose=TRUE,errorsAsWarnings=FALSE){
  dataProvidedtoFunction<-TRUE
  NMRdata2<-NMRdata
  if(methods::is(NMRdata,"mrbin")){
    transformations="PQN scaled"
    if(transformations %in% NMRdata$transformations){
      if(!errorsAsWarnings) stop("Data has been PQN scaled previously, this could corrupt the data.")
      warning("Data has been PQN scaled previously, this could corrupt the data.")
    }
    if("Log transformed" %in% NMRdata$transformations){
      if(!errorsAsWarnings) stop("Data has been log transformed previously, this would corrupt the data.")
      warning("Data has been log transformed previously, this would corrupt the data.")
    }
    dataProvidedtoFunction<-FALSE
    NMRdataTMP<-NMRdata$bins
    ignoreGlucose<-NMRdata$parameters$PQNIgnoreSugarArea
    dimension<-NMRdata$parameters$dimension
    ppmNames<-"borders"
    sugarArea<-NMRdata$parameters$PQNsugarArea
    minimumFeatures<-NMRdata$parameters$PQNminimumFeatures
    showHist<-NMRdata$parameters$PQNshowHist
  } else {
    NMRdataTMP<-as.matrix(NMRdata)
  }
  if(ncol(NMRdataTMP)>1){
    if(nrow(NMRdataTMP)>1){#Create synthetic median spectrum by averaging all spectra
      NMRdataTmp<-rbind(NMRdataTMP,apply(NMRdataTMP,2,mean))
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
              coordTmpAll<-cbind(apply(mrbin.env$mrbin$parameters$binRegions[,3:4],1,mean),
                  apply(mrbin.env$mrbin$parameters$binRegions[,1:2],1,mean))
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
              coordTmpAll<-apply(mrbin.env$mrbin$parameters$binRegions[,1:2],1,mean)
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
    }  else {
      NMRdataTmp2<-NMRdataTmp
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
      NMRdata$bins<-NMRdataTmp_scaledMedian
      NMRdata$parameters$medians<-medianFoldChanges
    } else {
      NMRdata<-NMRdataTmp_scaledMedian
    }
   } else {
      warning("Too few samples to perform PQN normalization.")
   }
 } else {
   warning("Too few bins to perform PQN scaling.")
 }
 if(methods::is(NMRdata,"mrbin")){
   NMRdata2<-editmrbin(mrbinObject=NMRdata2,functionName="mrbin::PQNScaling",
       versionNumber=as.character(utils::packageVersion("mrbin")),
       bins=NMRdata$bins, parameters=NMRdata$parameters, transformations=transformations,verbose=verbose)
   NMRdata<-NMRdata2
 }
 invisible(NMRdata)
}

#' A function for scaling to individual dilution factors.
#'
#' This function performs sample-wise scaling of binned data to correct for dilution
#' through different sample volumes used, or for different sample weights. All bin
#' values of one sample are multiplied by the corresponding dilution factor.
#' @param mrbinResults An mrbin object
#' @param errorsAsWarnings If TRUE, errors will be turned into warnings. Should be used with care, as errors indicate undocumented changes to the data.
#' @param verbose Should a summary be printed?
#' @return An invisible mrbin object containing scaled NMR data.
#' @export
#' @examples
#' mrbinResults<-mrbin(silent=TRUE,setDefault=TRUE,parameters=list(dimension="1D",
#'          binwidth1D=0.05,PQNScaling="No",PCA="No",tryParallel=FALSE,logTrafo="No",
#'          NMRfolders=c(system.file("extdata/1/10/pdata/10",package="mrbin"),
#'                       system.file("extdata/2/10/pdata/10",package="mrbin"),
#'                       system.file("extdata/3/10/pdata/10",package="mrbin"))),
#'          metadata=list(dilutionFactors=c(.75,1,.5)))
#' mrbinResults<-dilutionCorrection(mrbinResults)

dilutionCorrection<-function(mrbinResults,verbose=TRUE,errorsAsWarnings=FALSE){
    transformations="Dilution corrected"
    if(transformations %in% mrbinResults$transformations){
      if(!errorsAsWarnings) stop("Data has been dilution corrected previously, this would corrupt the data.")
      warning("Data has been dilution corrected previously, this would corrupt the data.")
    }
    if("Log transformed" %in% mrbinResults$transformations){
      if(!errorsAsWarnings) stop("Data has been log transformed previously, this would corrupt the data.")
      warning("Data has been log transformed previously, this would corrupt the data.")
    }
    if("PQN scaled" %in% mrbinResults$transformations){
      if(!errorsAsWarnings) stop("Data has been PQN transformed previously, this would corrupt the data.")
      warning("Data has been PQN transformed previously, this would corrupt the data.")
    }
    NMRdataTMP<-mrbinResults$bins
 if(length(mrbinResults$metadata$dilutionFactors)==nrow(NMRdataTMP)){
   for(i in 1:length(mrbinResults$metadata$dilutionFactors)){
     NMRdataTMP[i,]<-NMRdataTMP[i,]*mrbinResults$metadata$dilutionFactors[i]
   }
   mrbinResults<-editmrbin(mrbinObject=mrbinResults,
       functionName="mrbin::dilutionCorrection",
       versionNumber=as.character(utils::packageVersion("mrbin")),
       bins=NMRdataTMP,
       transformations=transformations,verbose=verbose)
 } else {
    if(!errorsAsWarnings) stop("Data dimension does not match the number of provided dilution factors.")
    warning("Data dimension does not match the number of provided dilution factors.")
 }
 invisible(mrbinResults)
}

#' A function for scaling to unit variance.
#'
#' This function performs scaling of binned data to unit variance so that each
#' bin has variance 1 and mean 0. This is
#' rarely necessary, but might be advantageous, e.g. in artificial neural networks.
#' @param mrbinResults An mrbin object
#' @param errorsAsWarnings If TRUE, errors will be turned into warnings. Should be used with care, as errors indicate undocumented changes to the data.
#' @param verbose Should a summary be printed?
#' @return An invisible mrbin object containing scaled NMR data.
#' @export
#' @examples
#' mrbinResults<-mrbin(silent=TRUE,setDefault=TRUE,parameters=list(dimension="1D",
#'          binwidth1D=0.05,PQNScaling="No",PCA="No",tryParallel=FALSE,logTrafo="No",
#'          NMRfolders=c(system.file("extdata/1/10/pdata/10",package="mrbin"),
#'                       system.file("extdata/2/10/pdata/10",package="mrbin"),
#'                       system.file("extdata/3/10/pdata/10",package="mrbin"))))
#' mrbinResults<-unitVarianceScaling(mrbinResults)

unitVarianceScaling<-function(mrbinResults,verbose=TRUE,errorsAsWarnings=FALSE){
    transformations="Unit variance scaled"
    if(transformations %in% mrbinResults$transformations){
      if(!errorsAsWarnings) stop("Data has been unit variance scaled previously, this would corrupt the data.")
      warning("Data has been unit variance scaled previously, this would corrupt the data.")
    }
    #if("Log transformed" %in% mrbinResults$transformations){
    #  if(!errorsAsWarnings) stop("Data has been log transformed previously, this would corrupt the data.")
    #  warning("Data has been log transformed previously, this would corrupt the data.")
    #}
    #if("PQN scaled" %in% mrbinResults$transformations){
    #  if(!errorsAsWarnings) stop("Data has been PQN transformed previously, this would corrupt the data.")
    #  warning("Data has been PQN transformed previously, this would corrupt the data.")
    #}
    NMRdataTMP<-mrbinResults$bins
   for(i in 1:ncol(mrbinResults$bins)){
     if(stats::var(NMRdataTMP[,i])==0) stop(paste("Variance was 0 for feature ",colnames(NMRdataTMP)[i],sep=""))
     NMRdataTMP[,i]<-(NMRdataTMP[,i]-mean(NMRdataTMP[,i]))/stats::var(NMRdataTMP[,i])^.5
   }
   mrbinResults<-editmrbin(mrbinObject=mrbinResults,
       functionName="mrbin::unitVarianceScaling",
       versionNumber=as.character(utils::packageVersion("mrbin")),
       bins=NMRdataTMP,
       transformations=transformations,verbose=verbose)
 invisible(mrbinResults)
}


#' A function for plotting PCA plots.
#'
#' This function performs PCA, then plots PC1 and PC2.
#' @param mrbinResults An mrbin object
#' @param defineGroups Should groups be colored differently?
#' @param loadings Should loadings be plotted instead of scores?
#' @param legendPosition Where should the legend be plotted, Defaults to "left", other options include "top", "topright", etc.
#' @param annotate Should loadings be annotated with metabolite identities, if available in $metadata?
#' @param verbose Should a summary be displayed?
#' @return An invisible prcomp result object
#' @export
#' @examples
#' mrbinResults<-mrbin(silent=TRUE,setDefault=FALSE,parameters=list(dimension="2D",
#'     binRegion=c(8,1,15,140),binwidth2D=0.1,binheight=4,solventRegion=c(5.5,4.2),
#'     PQNScaling="No",noiseRemoval="Yes",trimZeros="Yes",tryParallel=FALSE,
#'     fixNegatives="No",logTrafo="No",PCA="No",signal_to_noise2D=10,
#'     NMRfolders=c(system.file("extdata/1/12/pdata/10",package="mrbin"),
#'                  system.file("extdata/2/12/pdata/10",package="mrbin"))))
#' plotPCA(mrbinResults)
plotPCA<-function(mrbinResults,defineGroups=TRUE,loadings=FALSE,legendPosition="bottomleft",
 annotate=TRUE,verbose=TRUE){
 PCA<-NULL
 if(!is.null(mrbinResults$bins)){
    checkmrbin(mrbinResults,verbose=verbose)
    if(is.null(mrbinResults$metadata$factors)){
      FactorsTMP<-factor(rep("Group 0",nrow(mrbinResults$bins)))
    } else {
      FactorsTMP<-mrbinResults$metadata$factors
    }
    colorPalette<-c("black",grDevices::rainbow(length(levels(FactorsTMP))-1))
    if(annotate){
      mrbinResults<-annotatemrbin(mrbinResults)
      if(length(mrbinResults$metadata$annotations)==ncol(mrbinResults$bins)){
        colnames(mrbinResults$bins)<-mrbinResults$metadata$annotations
      }
    }
    if(nrow(mrbinResults$bins)>1){
      graphics::par(mgp=c(1,1,0))
      numlevels<-NULL
      if(defineGroups){
        for(i in 1:nlevels(FactorsTMP)) numlevels<-c(numlevels,as.numeric(
                         FactorsTMP[which(FactorsTMP==levels(FactorsTMP)[i])][1]))
      }
      numlevels2<-numlevels
      #move PCA to a new function in future
      PCA<-stats::prcomp(mrbinResults$bins)
      addTMP1<-.15*(max(PCA$x[,1])-min(PCA$x[,1]))
      xlimTMP<-c(min(PCA$x[,1])-addTMP1,max(PCA$x[,1])+addTMP1)
      addTMP1rot<-.25*(max(PCA$rotation[,1])-min(PCA$rotation[,1]))
      xlimTMProt<-c(min(PCA$rotation[,1])-addTMP1rot,max(PCA$rotation[,1])+addTMP1rot)
      xlabTMP<-paste("PC1 (",round(100*(PCA$sdev[1]^2)/sum(PCA$sdev^2),1),"%)",sep="")
      if(ncol(PCA$x)>1){
        addTMP2<-.15*(max(PCA$x[,2])-min(PCA$x[,2]))
        ylimTMP<-c(min(PCA$x[,2])-addTMP2,max(PCA$x[,2])+addTMP2)
        ylabTMP<-paste("PC2 (",round(100*(PCA$sdev[2]^2)/sum(PCA$sdev^2),1),"%)",sep="")
        addTMP2rot<-.25*(max(PCA$rotation[,2])-min(PCA$rotation[,2]))
        ylimTMProt<-c(min(PCA$rotation[,2])-addTMP2rot,max(PCA$rotation[,2])+addTMP2rot)
      } else {
        ylimTMP<-NULL
        ylabTMP<-""
        ylimTMProt<-NULL
      }
      if(defineGroups){
        PCAFactors<-as.numeric(FactorsTMP)
        PCAFactors2<-PCAFactors+14
        PCAFactors2[PCAFactors2>25]<-PCAFactors2[PCAFactors2>25]+7#pch 26 through 32 are missing
        numlevels2<-numlevels2+14
        numlevels2[numlevels2>25]<-numlevels2[numlevels2>25]+7#pch 26 through 32 are missing
      } else {
        PCAFactors<-rep(1,nrow(mrbinResults$bins))
        PCAFactors2<-PCAFactors+14
        numlevels2<-numlevels2+14
      }
      if(loadings){#plot loadings: do only if parameter set
        graphics::plot(PCA$rotation,xlim=xlimTMProt,ylim=ylimTMProt,xlab="PC1",ylab="PC2",
                       pch=16,cex=.75,main="PCA Loadings Plot", ask=FALSE, xaxt='n', yaxt='n')
        #display only top feature names
        topFeatures<-order(PCA$rotation[,1]^2+PCA$rotation[,2]^2,
          decreasing=TRUE)[1:min(nrow(PCA$rotation),50)]
        graphics::text(PCA$rotation[topFeatures,],labels=rownames(PCA$rotation)[topFeatures]
                       ,pos=4,cex=1)
      } else {
        graphics::plot(PCA$x
             ,xlim=xlimTMP,
             ylim=ylimTMP,
             pch=PCAFactors2,
             xaxt='n', yaxt='n',
             col=colorPalette[PCAFactors],
             main="PCA",
             xlab=xlabTMP,
             ylab=ylabTMP,
             cex=.75, ask=FALSE)
        graphics::text(PCA$x,labels=paste(substr(rownames(PCA$x),1,
          mrbinResults$parameters$PCAtitlelength)),pos=3,cex=1,
               col=colorPalette[PCAFactors])
          if(defineGroups) graphics::legend(legendPosition,
                  legend=levels(FactorsTMP),
                  col=colorPalette[numlevels],pch=numlevels2,cex=1,bg=NULL)#"white"
      }
      utils::flush.console()
    } else {
       warning("Too few samples to perform PCA.")
    }
 }
 invisible(PCA)
}

#' A function for plotting quality indicators, including PCA plots.
#'
#' This function plots boxplots (bin-wise and sample-wise) as visual quality indicators.
#' It also performs PCA, then plots PC1 and PC2 and loading plots.
#' @param mrbinResults An mrbin object
#' @param defineGroups Should group membership be highlighted in PCA?
#' @param process If set to FALSE, the file name will be extended by "Raw" to indicate that data has not been processed yet
#' @param silent If set to TRUE, plots will be saved but not shown for the binning step for speed purposes
#' @return {None}
#' @export
#' @examples
#' mrbinResults<-mrbin(silent=TRUE,setDefault=FALSE,parameters=list(dimension="2D",
#'     binRegion=c(8,1,15,140),binwidth2D=0.2,binheight=4,solventRegion=c(5.5,4.2),
#'     PQNScaling="No",noiseRemoval="Yes",trimZeros="Yes",tryParallel=FALSE,
#'     fixNegatives="No",logTrafo="No",PCA="No",signal_to_noise2D=10,
#'     NMRfolders=c(system.file("extdata/1/12/pdata/10",package="mrbin"),
#'                  system.file("extdata/2/12/pdata/10",package="mrbin"))))
#' plotResults(mrbinResults)

plotResults<-function(mrbinResults,defineGroups=TRUE,process=TRUE,silent=FALSE){
 if(!is.null(mrbinResults$bins)){
    #first write plot to object, then plot object (speed!)
    outputFileName<-mrbinResults$parameters$outputFileName
    mainTitle<-"Processed data"
    if(!process){
      outputFileName<-paste(outputFileName,"Raw",sep="")
      mainTitle<-"Raw data"
    }
    if(mrbinResults$parameters$saveFiles=="Yes"){
      grDevices::pdf(paste(outputFileName,"plot.pdf",sep=""))
    } else {
      grDevices::pdf(NULL)
    }
    grDevices::dev.control(displaylist="enable")
    oldpar<-graphics::par("mar","mfrow","mgp")
    on.exit(graphics::par(oldpar))
    devAskNewPage(ask = FALSE)
    graphics::par(mfrow=c(2,2),mar=c(3.1,2,2,0.5))
    #To save memory: avoid plotting thousand of small rectangles. Saving as pdf might take a long time.
    #This might happen when no or little noise removal was used
    #simplify rectangles by merging rectangles rowwise and then columnwise(2D)
    binRegionsTMP<-mrbinResults$parameters$binRegions[1,,drop=FALSE]
    if(nrow(mrbinResults$parameters$binRegions)>1){
      for(i in 2:nrow(mrbinResults$parameters$binRegions)){
         if(identical(mrbinResults$parameters$binRegions[i,3:4],binRegionsTMP[nrow(binRegionsTMP),3:4])&#if top and bottom are same
            mrbinResults$parameters$binRegions[i,1]==binRegionsTMP[nrow(binRegionsTMP),2]){#expand to right
              binRegionsTMP[nrow(binRegionsTMP),2]<-mrbinResults$parameters$binRegions[i,2]
         } else {#expand to left
           if(identical(mrbinResults$parameters$binRegions[i,3:4],binRegionsTMP[nrow(binRegionsTMP),3:4])&#if top and bottom are same
              mrbinResults$parameters$binRegions[i,2]==binRegionsTMP[nrow(binRegionsTMP),1]){
                binRegionsTMP[nrow(binRegionsTMP),1]<-mrbinResults$parameters$binRegions[i,1]
           } else {
              binRegionsTMP<-rbind(binRegionsTMP, mrbinResults$parameters$binRegions[i,])
           }
         }
      }
    }
    binRegionsTMP2<-binRegionsTMP[1,,drop=FALSE]
    if(nrow(binRegionsTMP)>1){
      for(i in 2:nrow(binRegionsTMP)){#expand down
         if(identical(binRegionsTMP[i,1:2],binRegionsTMP2[nrow(binRegionsTMP2),1:2])&#if l,r are same
            binRegionsTMP[i,3]==binRegionsTMP2[nrow(binRegionsTMP2),4]){#expand to bottom
              binRegionsTMP2[nrow(binRegionsTMP2),4]<-binRegionsTMP[i,4]
         } else {#expand to top
           if(identical(binRegionsTMP[i,1:2],binRegionsTMP2[nrow(binRegionsTMP2),1:2])&#if lr are same
              binRegionsTMP[i,4]==binRegionsTMP2[nrow(binRegionsTMP2),3]){
                binRegionsTMP2[nrow(binRegionsTMP2),3]<-binRegionsTMP[i,3]
           } else {
              binRegionsTMP2<-rbind(binRegionsTMP2, binRegionsTMP[i,])
           }
         }
      }
    }
    rectangleRegionsTMP<-NULL
    if(mrbinResults$parameters$dimension=="1D"){
      region<-c(max(mrbinResults$parameters$binRegions[,1:2],na.rm=TRUE),
        min(mrbinResults$parameters$binRegions[,1:2],na.rm=TRUE),0,100)
    } else {#2D
      region<-c(max(mrbinResults$parameters$binRegions[,1:2],na.rm=TRUE),
        min(mrbinResults$parameters$binRegions[,1:2],na.rm=TRUE),
        min(mrbinResults$parameters$binRegions[,3:4],na.rm=TRUE),
        max(mrbinResults$parameters$binRegions[,3:4],na.rm=TRUE))
    }
    if(nrow(binRegionsTMP2)<9000){
      rectangleRegionsTMP<-binRegionsTMP
    }
    if(mrbinResults$parameters$dimension=="1D"){
      plotMultiNMR(region=region,rectangleRegions=rectangleRegionsTMP,
          color="black",rectangleColors="green", rectangleFront=FALSE,
          manualScale=FALSE,maxPlots=2,plotTitle=mainTitle,restrictToRange=TRUE,
          dimension=mrbinResults$parameters$dimension,enableSplit=FALSE)
    }
    if(mrbinResults$parameters$dimension=="2D"){
      plotMultiNMR(rectangleRegions=rectangleRegionsTMP,
          color="darkgray",rectangleColors="blue", rectangleFront=FALSE,#TRUE,
          manualScale=FALSE,maxPlots=2,plotTitle=mainTitle,restrictToRange=TRUE,
          dimension=mrbinResults$parameters$dimension,enableSplit=FALSE)
    }
    axisCex1<-.1
    if(ncol(mrbinResults$bins)<100) axisCex1<-.25
    if(ncol(mrbinResults$bins)<25) axisCex1<-.5
    if(ncol(mrbinResults$bins)<15) axisCex1<-1
    axisCex2<-.1
    if(nrow(mrbinResults$bins)<100) axisCex2<-.25
    if(nrow(mrbinResults$bins)<25) axisCex2<-.5
    if(nrow(mrbinResults$bins)<15) axisCex2<-1
    graphics::boxplot(mrbinResults$bins[,1:min(ncol(mrbinResults$bins),5000)],main="Features",
      xlab="Bins",ylab="",boxwex=1,ask=FALSE,xaxt="n")
    graphics::axis(1,las=2,cex.axis=axisCex1,at=1:ncol(mrbinResults$bins),labels=colnames(mrbinResults$bins))
    graphics::boxplot(t(mrbinResults$bins),main="Samples",
      xlab="Samples",ylab="",boxwex=1,ask=FALSE,xaxt="n")
    graphics::axis(1,las=2,cex.axis=axisCex2,at=1:nrow(mrbinResults$bins),
      labels=substr(rownames(mrbinResults$bins),1,mrbinResults$parameters$PCAtitlelength))
    utils::flush.console()
    plotPCA(mrbinResults,defineGroups=defineGroups,loadings=FALSE,verbose=FALSE)
    plotRecordingTMP<-grDevices::recordPlot()
    invisible(grDevices::dev.off())
    if(process|(!process&!silent)) grDevices::replayPlot(plotRecordingTMP)
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
#' @param restrictToRange Restrict plot area to range of available data points. Defaults to FALSE
#' @param currentSpectrumOriginal Optional spectral data. If omitted, data from the environment variables is used
#' @param perspective If TRUE, a perspective plot will be displayed for 2D data instead of the regular topographic view
#' @param noise If provided, a line or plane at this level will be added to the plot to indicate noise level
#' @param dimension "1D' or "2D". If not provided, this will be deduced from the data
#' @return {None}
#' @export
#' @examples
#' mrbin(silent=TRUE,setDefault=TRUE,parameters=list(dimension="1D",binwidth1D=.1,
#'          PQNScaling="No",noiseRemoval="No",trimZeros="No",tryParallel=FALSE,
#'          fixNegatives="No",logTrafo="No",PCA="No",verbose=TRUE,
#'          NMRfolders=system.file("extdata/1/10/pdata/10",package="mrbin")))
#' plotNMR()

plotNMR<-function(region=NULL,rectangleRegions=NULL,
                   rectangleColors=c("green","orange","blue","red","yellow","gray","purple"),
                   rectangleFront=FALSE,
                   polygonRegion=NULL,
                   color=NULL,add=FALSE,showGrid=FALSE,
                   manualScale=TRUE,plotTitle="",
                   restrictToRange=FALSE,currentSpectrumOriginal=NULL,
                   perspective=FALSE,noise=NULL,dimension=NULL){
 #if(!is.null(noise)) noise<-1.1*noise#increase noise level, as whole bin must exceed this level, not just top of a peak
 if(is.null(currentSpectrumOriginal)){
   currentSpectrumOriginal<-mrbin.env$mrbinTMP$currentSpectrumOriginal
 }
 if(is.null(dimension)){
   if(is.matrix(currentSpectrumOriginal)){
     dimension<-"2D"
   } else {
     dimension<-"1D"
   }
 } else {
   if((dimension=="1D"&is.matrix(currentSpectrumOriginal))|
      (dimension=="2D"&!is.matrix(currentSpectrumOriginal))){#mismatch between dimension and current loaded data
      currentSpectrumOriginal<-NULL
   }
 }
 #if(!is.null(currentSpectrumOriginal)){
   devAskNewPage(ask = FALSE)
   if(is.null(region)){
       if(is.null(mrbin.env$mrbinplot$plotRegion)){
          if(!is.null(currentSpectrumOriginal)){
            if(dimension=="2D"){
                mrbin.env$mrbinplot$plotRegion<-
                          c(max(as.numeric(colnames(currentSpectrumOriginal))),
                            min(as.numeric(colnames(currentSpectrumOriginal))),
                            min(as.numeric(rownames(currentSpectrumOriginal))),
                            max(as.numeric(rownames(currentSpectrumOriginal))))
            } else {
                mrbin.env$mrbinplot$plotRegion<-
                         c(max(as.numeric(names(currentSpectrumOriginal))),
                         min(as.numeric(names(currentSpectrumOriginal))),
                         -10,160)
            }
          } else {
            mrbin.env$mrbinplot$plotRegion<-c(12,-1,-10,160)
          }
       } #else {
       region<-mrbin.env$mrbinplot$plotRegion
       #}
   }
   if(length(region)==1){
     if(region=="all"){
       if(!is.null(currentSpectrumOriginal)){
          if(dimension=="2D"){
                region<-c(max(as.numeric(colnames(currentSpectrumOriginal))),
                            min(as.numeric(colnames(currentSpectrumOriginal))),
                            min(as.numeric(rownames(currentSpectrumOriginal))),
                            max(as.numeric(rownames(currentSpectrumOriginal))))
          } else {
                region<-c(max(as.numeric(names(currentSpectrumOriginal))),
                         min(as.numeric(names(currentSpectrumOriginal))),
                         -10,160)
          }
       } else {
          region<-c(12,-1,-10,160)
       }
     }
   }
  if(restrictToRange){
   if(!is.null(currentSpectrumOriginal)){
      if(dimension=="2D"){#2d spectra
          region[1]<-min(region[1],max(as.numeric(colnames(currentSpectrumOriginal))))
          region[2]<-max(region[2],min(as.numeric(colnames(currentSpectrumOriginal))))
          region[3]<-max(region[3],min(as.numeric(rownames(currentSpectrumOriginal))))
          region[4]<-min(region[4],max(as.numeric(rownames(currentSpectrumOriginal))))
      } else {
         region[1]<-min(region[1],max(as.numeric(names(currentSpectrumOriginal))))
         region[2]<-max(region[2],min(as.numeric(names(currentSpectrumOriginal))))
      }
   }
  }
   if(dimension=="2D"#is.matrix(currentSpectrumOriginal)
   ){#2D spectra
    if(is.null(color)) color<-rev(grDevices::rainbow(mrbin.env$mrbinplot$nContours))
    spectrumTMP<-NULL
    if(!is.null(currentSpectrumOriginal)){
      TMP_1<-which(as.numeric(rownames(currentSpectrumOriginal))<region[4]&
                   as.numeric(rownames(currentSpectrumOriginal))>=region[3])
      if(length(TMP_1)==0)  TMP_1<-1:nrow(currentSpectrumOriginal)
      TMP_2<-which(as.numeric(colnames(currentSpectrumOriginal))>=region[2]&
                         as.numeric(colnames(currentSpectrumOriginal))<region[1])
      if(length(TMP_2)==0)  TMP_2<-1:ncol(currentSpectrumOriginal)
      spectrumTMP<-currentSpectrumOriginal[TMP_1,TMP_2]
    }
    if(perspective){
       zlimTMP<-c(min(spectrumTMP,na.rm=TRUE),max(spectrumTMP,na.rm=TRUE))
       if(!is.null(noise)){
         zlimTMP[1]<-sort(spectrumTMP)[ceiling(.25*length(spectrumTMP))]#quantile(spectrumTMP,.15)
         zlimTMP<-c(min(0,zlimTMP[1]),noise)
       }
       zlimTMP[2]<-zlimTMP[2]+0.8*abs(zlimTMP[2])
       if(length(spectrumTMP)>0){
         spectrumTMP2<-spectrumTMP
         spectrumTMP2[spectrumTMP2<zlimTMP[1]]<-zlimTMP[1]+0.01*abs(zlimTMP[1])#avoid warning messages
         spectrumTMP2[spectrumTMP2>zlimTMP[2]]<-1.0*zlimTMP[2]#avoid warning messages
         pmat<-graphics::persp(#x=-as.numeric(rownames(spectrumTMP)),
                #y=-as.numeric(colnames(spectrumTMP)),
                z=t(spectrumTMP2),zlim=zlimTMP,zlab="",shade=.5,#.7,
                border=NA,col="green4",
                phi=5,ltheta = -65, lphi = 55,main=plotTitle,#theta=-45,
                #ticktype="detailed",
                axes=FALSE)
         axis(1,at=0,line=-1.35,lwd=0,
           labels=paste(region[1],"-",region[2],"ppm, ",region[3],"-",region[4],"ppm",sep="",collapse=""))
         if(!is.null(noise)){
           graphics::lines(grDevices::trans3d(x=c(0,0,1,1),y=c(1,0,0,1),z=noise,pmat=pmat),col="red",lwd=2)
           graphics::lines(grDevices::trans3d(x=c(0,1),y=c(1,1),z=noise,pmat=pmat),col="red",lty=3)
         }
         utils::flush.console()
       }
    } else {
      if(manualScale){
       if(!is.null(currentSpectrumOriginal)){
        if(sum(spectrumTMP<(mrbin.env$mrbinplot$lowestContour*max(
             currentSpectrumOriginal)))>0){
             spectrumTMP[spectrumTMP<(mrbin.env$mrbinplot$lowestContour*max(
             currentSpectrumOriginal))]<-0
        }
       }
      } else {
       if(!is.null(currentSpectrumOriginal)){
        #Set values below noise to zero, noise estimated from 90% percentile
        if(sum(spectrumTMP>(3.5*sort(spectrumTMP)[floor(length(spectrumTMP)*.9)]))>10){
          spectrumTMP[spectrumTMP<=(3.5*sort(spectrumTMP)[floor(length(spectrumTMP)*.9)])]<-1e-8
        } else {
          spectrumTMP[spectrumTMP<=(sort(spectrumTMP,decreasing=TRUE)[10])]<-1e-8
        }
        #spectrumTMP[spectrumTMP<=0]<-1e-8
        spectrumTMP<-log(100000*spectrumTMP/max(spectrumTMP))
       }
      }
      displaySize<-256#512#Reduce resolution for faster plotting of high-res spectra
      if(!is.null(currentSpectrumOriginal)){
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
      }
      if(!manualScale){ #mrbin.env$mrbinTMP$currentSpectrumOriginal mrbin.env$mrbinTMP$additionalPlots2D
        #display only a limited number of peaks to avoid crowded, noisy, slow plots
        if(!is.null(currentSpectrumOriginal)){
          while(sum(spectrumTMP>0)>max(c(.0001*length(spectrumTMP),min(c(150,.005*length(spectrumTMP)))))){
            spectrumTMP[spectrumTMP>0][spectrumTMP[spectrumTMP>0]<sort(
                       spectrumTMP[spectrumTMP>0])[floor(.5*length(spectrumTMP[spectrumTMP>0]))]]<-0
          }
        }
      }
      options(max.contour.segments=40)
      if(!add){
          graphics::plot(NULL,NULL,
            type="l",xlim=c(-region[1],-region[2]),
            ylim=c(-region[4],-region[3]),main=plotTitle,xaxt="n",yaxt="n",
            xlab="Chemical shift [ppm]",ylab="Chemical shift [ppm]")
          if(!is.null(currentSpectrumOriginal)){
              magnitude2<-10^round(log(max(as.numeric(colnames(spectrumTMP)))-
                          min(as.numeric(colnames(spectrumTMP))),base=10))/10
              magnitude1<-10^round(log(max(as.numeric(rownames(spectrumTMP)))-
                           min(as.numeric(rownames(spectrumTMP))),base=10))/10
              at2<--(0:100*magnitude1+floor(min(as.numeric(rownames(spectrumTMP)))/magnitude1)*magnitude1)
              labels2<-(0:100*magnitude1+floor(min(as.numeric(rownames(spectrumTMP)))/magnitude1)*magnitude1)
              at1<-(0:100*magnitude2+floor(min(-as.numeric(colnames(spectrumTMP)))/magnitude2)*magnitude2)
              labels1<--(0:100*magnitude2+floor(min(-as.numeric(colnames(spectrumTMP)))/magnitude2)*magnitude2)
          } else {
              magnitude2<-10^round(log(max(region[1:2])-min(region[1:2]),base=10))/10
              magnitude1<-10^round(log(max(region[3:4])-min(region[3:4]),base=10))/10
              at2<--(0:100*magnitude1+floor(min(region[3:4])/magnitude1)*magnitude1)
              labels2<-(0:100*magnitude1+floor(min(region[3:4])/magnitude1)*magnitude1)
              at1<-(0:100*magnitude2+floor(min(-region[1:2])/magnitude2)*magnitude2)
              labels1<--(0:100*magnitude2+floor(min(-region[1:2])/magnitude2)*magnitude2)
          }
          graphics::axis(2,at=at2,labels=labels2)
          graphics::axis(1,at=at1,labels=labels1)
      }
      if(!add&!is.null(rectangleRegions)&!rectangleFront){
          graphics::rect(xleft=-rectangleRegions[,1], ybottom=-rectangleRegions[,4],
                         xright=-rectangleRegions[,2], ytop=-rectangleRegions[,3],
                         col = rectangleColors, border = rectangleColors,lwd=1)
          graphics::box()
      }
      if(showGrid){
        if(!is.null(currentSpectrumOriginal)){
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
      }
      if(!add&!is.null(polygonRegion)){
          graphics::polygon(-1*polygonRegion[,2],-1*polygonRegion[,1],
                    #c(0,10,11,1),c(6e+5,0,1000,6e+5)
                    #col=NA, border = rectangleColors
                    col=rectangleColors, border = NA
                    )
          graphics::box()
      }
      if(!is.null(currentSpectrumOriginal)){
        suppressWarnings(graphics::contour(x = -as.numeric(colnames(spectrumTMP)),
          y = -as.numeric(rownames(spectrumTMP)),
          z = t(spectrumTMP)*mrbin.env$mrbinplot$intensityScale2D,
          levels =  (mrbin.env$mrbinplot$lowestContour+
                        .8*(1:mrbin.env$mrbinplot$nContours-1)/
                        (mrbin.env$mrbinplot$nContours)*(1-
                        mrbin.env$mrbinplot$lowestContour))* max(spectrumTMP)
          ,drawlabels=FALSE,col=color,lwd=1,add=TRUE))
      }
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
    }
   } else {  #1D
      if(is.null(color)) color<-"black"
      if(!is.null(currentSpectrumOriginal)){
        if(manualScale){
          ymin<-min(currentSpectrumOriginal)/mrbin.env$mrbinplot$intensityScale+mrbin.env$mrbinplot$intensityOffset
          ymax<-max(currentSpectrumOriginal)/mrbin.env$mrbinplot$intensityScale+mrbin.env$mrbinplot$intensityOffset
        } else {
          ymin<-0
          spectrumTMP<-currentSpectrumOriginal[as.numeric(names(currentSpectrumOriginal))<
                     region[1]&
                     as.numeric(names(currentSpectrumOriginal))>region[2]]
          ymax<-sort(spectrumTMP,decreasing=TRUE)[ceiling(.01*length(spectrumTMP))]
          if(!is.null(noise)){
           ymax<-max(min(spectrumTMP)*6,7*noise)
           ymin<-sort(spectrumTMP)[ceiling(.25*length(spectrumTMP))]#quantile(spectrumTMP,.15)
           ymin<-min(0,ymin)
          }
        }
      } else {
        ymin<-0
        ymax<-1
      }
      if(!add){
          graphics::plot(NULL,NULL,
            type="l",xlim=c(region[1],region[2]),
            ylim=c(ymin,ymax),main=plotTitle,
            xlab="Chemical shift [ppm]",ylab="Intensity")
      }
      if(!add&!is.null(rectangleRegions)){
          graphics::rect(xleft=rectangleRegions[,1], ybottom=0,
                    xright=rectangleRegions[,2], ytop=ymax*2,
                    col = rectangleColors, border = NA)
          graphics::box()
      }
      if(!is.null(currentSpectrumOriginal)){
        graphics::lines(as.numeric(names(currentSpectrumOriginal)),
                currentSpectrumOriginal,
                col=color)
        if(showGrid){
          GridxCoord<-as.numeric(names(currentSpectrumOriginal))
          Gridx<-currentSpectrumOriginal[GridxCoord<region[1]&GridxCoord>region[2]]
          GridxCoord<-GridxCoord[GridxCoord<region[1]&GridxCoord>region[2]]
          graphics::points(GridxCoord,Gridx,
                         col = "darkgray",pch="+",cex=.85)
        }
      }
      if(!is.null(noise)){
         graphics::abline(h=noise,col="red",lwd=2)
      }
      utils::flush.console()
    }
 #}
 Sys.sleep(0)#This forces RStudio to update the plot
}

#' A function for changing plotNMR plots.
#'
#' This function increases the intensity of the current NMR spectrum plot.
#' @param refreshPlot Refresh plot automatically. Defaults to TRUE
#' @param dimension Dimension to use. Defaults to "1D"
#' @return {None}
#' @export
#' @examples
#' mrbin(silent=TRUE,setDefault=TRUE,parameters=list(dimension="1D",binwidth1D=.1,
#'          PQNScaling="No",noiseRemoval="No",trimZeros="No",tryParallel=FALSE,
#'          fixNegatives="No",logTrafo="No",PCA="No",verbose=TRUE,
#'          NMRfolders=system.file("extdata/1/10/pdata/10",package="mrbin")))
#' plotNMR()
#' intPlus()

intPlus<-function(dimension="1D",refreshPlot=TRUE){#increase plot intensity
 if(dimension=="1D"){
   mrbin.env$mrbinplot$intensityScale<-mrbin.env$mrbinplot$intensityScale*2
 } else {
   mrbin.env$mrbinplot$intensityScale2D<-mrbin.env$mrbinplot$intensityScale2D*2
 }
 if(refreshPlot)  plotNMR()
 #}
}

#' A function for changing plotNMR plots.
#'
#' This function decreases the intensity of the current NMR spectrum plot.
#' @param dimension Dimension to use. Defaults to "1D"
#' @param refreshPlot Refresh plot automatically. Defaults to TRUE
#' @param value Set exact value. Defaults to NULL
#' @return {None}
#' @export
#' @examples
#' mrbin(silent=TRUE,setDefault=TRUE,parameters=list(dimension="1D",binwidth1D=.1,
#'          PQNScaling="No",noiseRemoval="No",trimZeros="No",tryParallel=FALSE,
#'          fixNegatives="No",logTrafo="No",PCA="No",verbose=TRUE,
#'          NMRfolders=system.file("extdata/1/10/pdata/10",package="mrbin")))
#' plotNMR()
#' intMin()

intMin<-function(dimension="1D",refreshPlot=TRUE,value=NULL){#decrease plot intensity
 if(dimension=="1D"){
  if(is.null(value)){
   mrbin.env$mrbinplot$intensityScale<-mrbin.env$mrbinplot$intensityScale*0.5
  } else {
   mrbin.env$mrbinplot$intensityScale<-as.numeric(value)
  }
 } else {
  if(is.null(value)){
   mrbin.env$mrbinplot$intensityScale2D<-mrbin.env$mrbinplot$intensityScale2D*0.5
  } else {
   mrbin.env$mrbinplot$intensityScale2D<-as.numeric(value)
  }

 }
 if(refreshPlot)  plotNMR()
 #}
}

#' A function for changing plotNMR plots.
#'
#' This function increases the minimum contour level of the current 2D NMR
#' spectrum plot.
#' @param refreshPlot Refresh plot automatically. Defaults to TRUE
#' @return {None}
#' @export
#' @examples
#' resetEnv()
#' addToPlot(folder=system.file("extdata/1/12/pdata/10",package="mrbin"),dimension="2D")
#' plotNMR()
#' contPlus()

contPlus<-function(refreshPlot=TRUE){#decrease plot intensity
 #if(!is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
   mrbin.env$mrbinplot$lowestContour<-mrbin.env$mrbinplot$lowestContour*1.5
   if(refreshPlot)  plotNMR()
 #}
}

#' A function for changing plotNMR plots.
#'
#' This function decreases the minimum contour level of the current 2D NMR
#' spectrum plot.
#' @param refreshPlot Refresh plot automatically. Defaults to TRUE
#' @return {None}
#' @export
#' @examples
#' resetEnv()
#' mrbin(silent=TRUE,parameters=list(dimension="2D",binwidth2D=0.5,
#'          binheight=3,PQNScaling="No",referenceScaling="No",binRegion=c(4,3,60,65),
#'          noiseRemoval="No",trimZeros="No",cropHSQC="No",tryParallel=FALSE,
#'          fixNegatives="No",logTrafo="No",PCA="No",verbose=TRUE,saveFiles="No",
#'          NMRfolders=c(system.file("extdata/1/12/pdata/10",package="mrbin"))))
#' plotNMR()
#' contMin()

contMin<-function(refreshPlot=TRUE){#decrease plot intensity
 #if(!is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
   mrbin.env$mrbinplot$lowestContour<-mrbin.env$mrbinplot$lowestContour*0.75
   if(refreshPlot)  plotNMR()
 #}
}

#' A function for changing plotNMR plots.
#'
#' This function changes the plot region of the current NMR plot. Can be called with
#' no arguments: zoom(). In this case the user will be asked for manual input.
#' @param left New left boundary
#' @param right New right boundary
#' @param top New top boundary
#' @param bottom New bottom boundary
#' @param refreshPlot Refresh plot automatically. Defaults to TRUE
#' @return {None}
#' @export
#' @examples
#' mrbin(silent=TRUE,setDefault=TRUE,parameters=list(dimension="1D",binwidth1D=.1,
#'          PQNScaling="No",noiseRemoval="No",trimZeros="No",tryParallel=FALSE,
#'          fixNegatives="No",logTrafo="No",PCA="No",verbose=TRUE,
#'          NMRfolders=system.file("extdata/1/10/pdata/10",package="mrbin")))
#' plotNMR()
#' zoom(left=4.6,right=2,top=10,bottom=150)

zoom<-function(left=NULL,right=NULL,top=NULL,bottom=NULL,refreshPlot=TRUE){
 #if(!is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
   if(is.null(left)){
      left<-readline(prompt=paste("Please set left limit, press enter to keep ",
                              mrbin.env$mrbinplot$plotRegion[1],": ",sep=""))
      if(left=="") left<-mrbin.env$mrbinplot$plotRegion[1]
   }
   if(is.null(right)){
      right<-readline(prompt=paste("Please set right limit, press enter to keep ",
                              mrbin.env$mrbinplot$plotRegion[2],": ",sep=""))
      if(right=="") right<-mrbin.env$mrbinplot$plotRegion[2]
   }
   #if(is.matrix(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
   if(is.null(top)){
      top<-readline(prompt=paste("Please set top limit, press enter to keep ",
                              mrbin.env$mrbinplot$plotRegion[3],": ",sep=""))
      if(top=="") top<-mrbin.env$mrbinplot$plotRegion[3]
   }
   if(is.null(bottom)){
      bottom<-readline(prompt=paste("Please set bottom limit, press enter to keep ",
                              mrbin.env$mrbinplot$plotRegion[4],": ",sep=""))
      if(bottom=="") bottom<-mrbin.env$mrbinplot$plotRegion[4]
   }
   #}
   mrbin.env$mrbinplot$plotRegion[1]<-as.numeric(left)
   mrbin.env$mrbinplot$plotRegion[2]<-as.numeric(right)
   mrbin.env$mrbinplot$plotRegion[3]<-as.numeric(top)
   mrbin.env$mrbinplot$plotRegion[4]<-as.numeric(bottom)
   if(refreshPlot)  plotNMR()
 #}
}

#' A function for changing plotNMR plots.
#'
#' This function zooms into the plot region of the current NMR plot.
#' @param refreshPlot Refresh plot automatically. Defaults to TRUE
#' @param x Change x axis? Defaults to TRUE
#' @param y Change y axis? Defaults to TRUE
#' @return {None}
#' @export
#' @examples
#' mrbin(silent=TRUE,setDefault=TRUE,parameters=list(dimension="1D",binwidth1D=.1,
#'          PQNScaling="No",noiseRemoval="No",trimZeros="No",tryParallel=FALSE,
#'          fixNegatives="No",logTrafo="No",PCA="No",verbose=TRUE,
#'          NMRfolders=system.file("extdata/1/10/pdata/10",package="mrbin")))
#' plotNMR()
#' zoomIn()

zoomIn<-function(refreshPlot=TRUE,x=TRUE,y=TRUE){#Zoom into NMR spectrum plot
 #if(!is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
   leftMax<-Inf
   rightMax<--Inf
   bottomMax<-Inf
   topMax<--Inf
   if(mrbin.env$mrbin$parameters$dimension=="1D"){
       if(!is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
         currentSpectrumOriginal<-mrbin.env$mrbinTMP$currentSpectrumOriginal
       } else {
         currentSpectrumOriginal<-mrbin.env$mrbinTMP$additionalPlots1D[[1]]
       }
       leftMax<-max(as.numeric(names(currentSpectrumOriginal)))
       rightMax<-min(as.numeric(names(currentSpectrumOriginal)))
       topMax<--10
       bottomMax<-160
   }
   if(mrbin.env$mrbin$parameters$dimension=="2D"){
       if(!is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
         currentSpectrumOriginal<-mrbin.env$mrbinTMP$currentSpectrumOriginal
       } else {
         currentSpectrumOriginal<-mrbin.env$mrbinTMP$additionalPlots2D[[1]]
       }
       leftMax<-max(as.numeric(colnames(currentSpectrumOriginal)))
       rightMax<-min(as.numeric(colnames(currentSpectrumOriginal)))
       topMax<-min(as.numeric(rownames(currentSpectrumOriginal)))
       bottomMax<-max(as.numeric(rownames(currentSpectrumOriginal)))
   }
   if(x){
     mrbin.env$mrbinplot$plotRegion[1]<-min(leftMax,
                               mrbin.env$mrbinplot$plotRegion[1]-(mrbin.env$mrbinplot$plotRegion[1]-mrbin.env$mrbinplot$plotRegion[2])/4)
     mrbin.env$mrbinplot$plotRegion[2]<-max(rightMax,
                               mrbin.env$mrbinplot$plotRegion[1]-(mrbin.env$mrbinplot$plotRegion[1]-mrbin.env$mrbinplot$plotRegion[2])*3/4)
   }
   if(y){
     mrbin.env$mrbinplot$plotRegion[3]<-max(topMax,
                               mrbin.env$mrbinplot$plotRegion[4]-(mrbin.env$mrbinplot$plotRegion[4]-mrbin.env$mrbinplot$plotRegion[3])*3/4)
     mrbin.env$mrbinplot$plotRegion[4]<-min(bottomMax,
                               mrbin.env$mrbinplot$plotRegion[4]-(mrbin.env$mrbinplot$plotRegion[4]-mrbin.env$mrbinplot$plotRegion[3])/4)
   }
   if(refreshPlot)  plotNMR()
 #}
}

#' A function for changing plotNMR plots.
#'
#' This function zooms out from the plot region of the current NMR plot.
#' @param refreshPlot Refresh plot automatically. Defaults to TRUE
#' @param x Change x axis? Defaults to TRUE
#' @param y Change y axis? Defaults to TRUE
#' @return {None}
#' @export
#' @examples
#' mrbin(silent=TRUE,setDefault=TRUE,parameters=list(dimension="1D",binwidth1D=.1,
#'          PQNScaling="No",noiseRemoval="No",trimZeros="No",tryParallel=FALSE,
#'          fixNegatives="No",logTrafo="No",PCA="No",verbose=TRUE,
#'          NMRfolders=system.file("extdata/1/10/pdata/10",package="mrbin")))
#' plotNMR()
#' zoomIn()
#' zoomOut()

zoomOut<-function(refreshPlot=TRUE,x=TRUE,y=TRUE){#Zoom out from NMR spectrum plot
 #if(!is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
   leftMax<-Inf
   rightMax<--Inf
   bottomMax<-Inf
   topMax<--Inf
   if(mrbin.env$mrbin$parameters$dimension=="1D"){
       if(!is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
         currentSpectrumOriginal<-mrbin.env$mrbinTMP$currentSpectrumOriginal
       } else {
         currentSpectrumOriginal<-mrbin.env$mrbinTMP$additionalPlots1D[[1]]
       }
       leftMax<-max(as.numeric(names(currentSpectrumOriginal)))
       rightMax<-min(as.numeric(names(currentSpectrumOriginal)))
       topMax<--10
       bottomMax<-160
   }
   if(mrbin.env$mrbin$parameters$dimension=="2D"){
       if(!is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
         currentSpectrumOriginal<-mrbin.env$mrbinTMP$currentSpectrumOriginal
       } else {
         currentSpectrumOriginal<-mrbin.env$mrbinTMP$additionalPlots2D[[1]]
       }
       leftMax<-max(as.numeric(colnames(currentSpectrumOriginal)))
       rightMax<-min(as.numeric(colnames(currentSpectrumOriginal)))
       topMax<-min(as.numeric(rownames(currentSpectrumOriginal)))
       bottomMax<-max(as.numeric(rownames(currentSpectrumOriginal)))
   }
   if(x){
     mrbin.env$mrbinplot$plotRegion[1]<-min(leftMax,
                               mrbin.env$mrbinplot$plotRegion[1]+(mrbin.env$mrbinplot$plotRegion[1]-mrbin.env$mrbinplot$plotRegion[2])/2)
     mrbin.env$mrbinplot$plotRegion[2]<-max(rightMax,
                               mrbin.env$mrbinplot$plotRegion[1]-(mrbin.env$mrbinplot$plotRegion[1]-mrbin.env$mrbinplot$plotRegion[2])*3/2)
   }
   if(y){
     mrbin.env$mrbinplot$plotRegion[3]<-max(topMax,
                               mrbin.env$mrbinplot$plotRegion[4]-(mrbin.env$mrbinplot$plotRegion[4]-mrbin.env$mrbinplot$plotRegion[3])*3/2)
     mrbin.env$mrbinplot$plotRegion[4]<-min(bottomMax,
                               mrbin.env$mrbinplot$plotRegion[4]+(mrbin.env$mrbinplot$plotRegion[4]-mrbin.env$mrbinplot$plotRegion[3])/2)
   }
   if(refreshPlot)  plotNMR()
 #}
}


#' A function for changing plotNMR plots.
#'
#' This function moves up or down the 1D plot region of the current NMR plot.
#' @param offsetValue The new offset value. Defaults to NULL
#' @return {None}
#' @export
#' @examples
#' setOffset(0)
setOffset<-function(offsetValue=NULL){
  if(!is.null(offsetValue)){
    mrbin.env$mrbinplot$intensityOffset<-offsetValue
  }
}


#' A function for changing plotNMR plots.
#'
#' This function moves left the plot region of the current NMR plot.
#' @param refreshPlot Refresh plot automatically. Defaults to TRUE
#' @return {None}
#' @export
#' @examples
#' resetEnv()
#' mrbin(silent=TRUE,parameters=list(dimension="1D",binwidth1D=.5,
#'          noiseRemoval="No",trimZeros="No",tryParallel=FALSE,
#'          PQNScaling="No",saveFiles="No",referenceScaling="No",
#'          fixNegatives="No",logTrafo="No",PCA="No",verbose=TRUE,
#'          NMRfolders=system.file("extdata/1/10/pdata/10",package="mrbin")))
#' plotNMR()
#' zoomIn()
#' left()

left<-function(refreshPlot=TRUE){
 #if(!is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
   leftMax<-Inf
   rightMax<--Inf
   bottomMax<-Inf
   topMax<--Inf
   if(mrbin.env$mrbin$parameters$dimension=="1D"&
     !is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
       leftMax<-max(as.numeric(names(mrbin.env$mrbinTMP$currentSpectrumOriginal)))
       rightMax<-min(as.numeric(names(mrbin.env$mrbinTMP$currentSpectrumOriginal)))
   }
   if(mrbin.env$mrbin$parameters$dimension=="2D"&
     !is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
       leftMax<-max(as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrumOriginal)))
       rightMax<-min(as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrumOriginal)))
   }
   mrbin.env$mrbinplot$plotRegion[1]<-min(leftMax,
                             mrbin.env$mrbinplot$plotRegion[1]+(mrbin.env$mrbinplot$plotRegion[1]-mrbin.env$mrbinplot$plotRegion[2])*.1)
   mrbin.env$mrbinplot$plotRegion[2]<-max(rightMax,
                             mrbin.env$mrbinplot$plotRegion[2]+(mrbin.env$mrbinplot$plotRegion[1]-mrbin.env$mrbinplot$plotRegion[2])*.1)
   if(refreshPlot)  plotNMR()
 #}
}

#' A function for changing plotNMR plots.
#'
#' This function moves right the plot region of the current NMR plot.
#' @param refreshPlot Refresh plot automatically. Defaults to TRUE
#' @return {None}
#' @export
#' @examples
#' resetEnv()
#' mrbin(silent=TRUE,parameters=list(dimension="1D",binwidth1D=.5,
#'          noiseRemoval="No",trimZeros="No",tryParallel=FALSE,
#'          PQNScaling="No",saveFiles="No",referenceScaling="No",
#'          fixNegatives="No",logTrafo="No",PCA="No",verbose=TRUE,
#'          NMRfolders=system.file("extdata/1/10/pdata/10",package="mrbin")))
#' plotNMR()
#' zoomIn()
#' right()

right<-function(refreshPlot=TRUE){
 #if(!is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
   leftMax<-Inf
   rightMax<--Inf
   bottomMax<-Inf
   topMax<--Inf
   if(mrbin.env$mrbin$parameters$dimension=="1D"&
     !is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
       leftMax<-max(as.numeric(names(mrbin.env$mrbinTMP$currentSpectrumOriginal)))
       rightMax<-min(as.numeric(names(mrbin.env$mrbinTMP$currentSpectrumOriginal)))
   }
   if(mrbin.env$mrbin$parameters$dimension=="2D"&
     !is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
       leftMax<-max(as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrumOriginal)))
       rightMax<-min(as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrumOriginal)))
   }
   mrbin.env$mrbinplot$plotRegion[1]<-min(leftMax,
                             mrbin.env$mrbinplot$plotRegion[1]-(mrbin.env$mrbinplot$plotRegion[1]-mrbin.env$mrbinplot$plotRegion[2])*.1)
   mrbin.env$mrbinplot$plotRegion[2]<-max(rightMax,
                             mrbin.env$mrbinplot$plotRegion[2]-(mrbin.env$mrbinplot$plotRegion[1]-mrbin.env$mrbinplot$plotRegion[2])*.1)
   if(refreshPlot)  plotNMR()
 #}
}

#' A function for changing plotNMR plots.
#'
#' This function moves down the plot region of the current NMR plot (only 2D).
#' @param refreshPlot Refresh plot automatically. Defaults to TRUE
#' @return {None}
#' @export
#' @examples
#' resetEnv()
#' mrbin(silent=TRUE,parameters=list(dimension="2D",binwidth2D=0.5,
#'          binheight=3,PQNScaling="No",referenceScaling="No",binRegion=c(4,3,60,65),
#'          noiseRemoval="No",trimZeros="No",cropHSQC="No",tryParallel=FALSE,
#'          fixNegatives="No",logTrafo="No",PCA="No",verbose=TRUE,saveFiles="No",
#'          NMRfolders=c(system.file("extdata/1/12/pdata/10",package="mrbin"))))
#' plotNMR()
#' zoomIn()
#' down()

down<-function(refreshPlot=TRUE){
 #if(!is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
   bottomMax<-Inf
   topMax<--Inf
   if(mrbin.env$mrbin$parameters$dimension=="2D"&
     !is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
       topMax<-min(as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrumOriginal)))
       bottomMax<-max(as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrumOriginal)))
   }
   mrbin.env$mrbinplot$plotRegion[3]<-max(topMax,
                             mrbin.env$mrbinplot$plotRegion[3]-(mrbin.env$mrbinplot$plotRegion[3]-mrbin.env$mrbinplot$plotRegion[4])*.1)
   mrbin.env$mrbinplot$plotRegion[4]<-min(bottomMax,
                             mrbin.env$mrbinplot$plotRegion[4]-(mrbin.env$mrbinplot$plotRegion[3]-mrbin.env$mrbinplot$plotRegion[4])*.1)
   if(refreshPlot)  plotNMR()
 #}
}

#' A function for changing plotNMR plots.
#'
#' This function moves up the plot region of the current NMR plot (only 2D).
#' @param refreshPlot Refresh plot automatically. Defaults to TRUE
#' @return {None}
#' @export
#' @examples
#' resetEnv()
#' mrbin(silent=TRUE,parameters=list(dimension="2D",binwidth2D=0.5,
#'          binheight=3,PQNScaling="No",referenceScaling="No",binRegion=c(4,3,60,65),
#'          noiseRemoval="No",trimZeros="No",cropHSQC="No",tryParallel=FALSE,
#'          fixNegatives="No",logTrafo="No",PCA="No",verbose=TRUE,saveFiles="No",
#'          NMRfolders=c(system.file("extdata/1/12/pdata/10",package="mrbin"))))
#' plotNMR()
#' zoomIn()
#' up()

up<-function(refreshPlot=TRUE){
 #if(!is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
   bottomMax<-Inf
   topMax<--Inf
   if(mrbin.env$mrbin$parameters$dimension=="2D"&
     !is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
       topMax<-min(as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrumOriginal)))
       bottomMax<-max(as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrumOriginal)))
   }
   mrbin.env$mrbinplot$plotRegion[3]<-max(topMax,
                             mrbin.env$mrbinplot$plotRegion[3]+(mrbin.env$mrbinplot$plotRegion[3]-mrbin.env$mrbinplot$plotRegion[4])*.1)
   mrbin.env$mrbinplot$plotRegion[4]<-min(bottomMax,
                             mrbin.env$mrbinplot$plotRegion[4]+(mrbin.env$mrbinplot$plotRegion[3]-mrbin.env$mrbinplot$plotRegion[4])*.1)
   if(refreshPlot)  plotNMR()
 #}
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
   useAsNames="Spectrum titles",useMeanIntensityForBins=FALSE
   ){#Bin NMR spectral data
  if(!is.null(folder)){
    warningMessage<-NULL
    NMRdataList<-readNMR(folder=folder,dimension=dimension,
                  NMRvendor=NMRvendor,useAsNames=useAsNames
                  )
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
          paste("Reference signal is negative: ",NMRdataList$currentSpectrumName,
          ". Please check phase and baseline. Results may be corrupted.",
          sep="")
      }
    }
    if(removeSolvent=="Yes") NMRdata<-removeSolvent2(NMRdata=NMRdata,
               dimension=dimension,solventRegion=solventRegion)
    if(removeAreas=="Yes") NMRdata<-removeAreas2(NMRdata=NMRdata,
               dimension=dimension,removeAreaList=removeAreaList)
    binData<-binSingleNMR(currentSpectrum=NMRdata,dimension=dimension,
              binRegions=binRegions,binMethod=binMethod,
              useMeanIntensityForBins=useMeanIntensityForBins,
              spectrumTitle=NMRdataList$currentSpectrumName)
    noiseData<-calculateNoise(NMRdata=NMRdataOriginal,
               pointsPerBin=binData$pointsPerBin,dimension=dimension,
               noiseRange1d=noiseRange1d,noiseRange2d=noiseRange2d,
               binRegions=binRegions,#must be original bin regions??????????????????????????? before removing solvent and areas
               useMeanIntensityForBins=useMeanIntensityForBins)#list
    if(referenceScaling=="Yes"){
      if(scalingFactor<(.2*noiseData$noise_level)){#3
        warningMessage<-c(warningMessage,paste(
          "Reference signal very low: ",NMRdataList$currentSpectrumName,
          ". Please check if reference peak is at 0ppm. Results may be corrupted.",
          sep=""))
      }
    }
    warningMessageTMP<-
       checkBaseline(NMRdata=NMRdataOriginal,dimension=dimension,
               currentSpectrumName=NMRdataList$currentSpectrumName,
               noiseRange1d=noiseRange1d,noiseRange2d=noiseRange2d)
    if(!is.null(warningMessageTMP)){
      warningMessage<-c(warningMessage,warningMessageTMP)
    }
    invisible(list(binTMP=binData$binTMP,
       meanNumberOfPointsPerBin_TMP=binData$pointsPerBin,#mrbin.env$mrbinTMP$meanNumberOfPointsPerBin_TMP
       noise_level_Raw_TMP=noiseData$noise_level,#noise level per sample before reference scaling
       noise_level_TMP=noiseData$noise_level_TMP/scalingFactor,#noise level per bin and sample after reference scaling #noiseDataScaled$noise_level_TMP,
       currentSpectrumName=NMRdataList$currentSpectrumName,
       AcquPars=NMRdataList$AcquPars,
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
#' @param checkFiles Only check if the folder exists or contains NMR data. Defaults to FALSE
#' @return An (invisible) list containing spectral data and the spectrum name
#' @keywords internal
#' @noRd
#' @examples
#' \donttest{ readNMR() }

readNMR<-function(folder=NULL,dimension=NULL,onlyTitles=FALSE,
          NMRvendor="Bruker",useAsNames="Spectrum titles",
          checkFiles=FALSE#,readAcqus=FALSE
          ){#Read NMR spectral data
 #if(!is.null(mrbin.env$mrbinTMP$currentFolder)){
  if(NMRvendor=="Bruker"){
      currentSpectrum<-readBruker(folder=folder,dimension=dimension,
                      onlyTitles=onlyTitles,useAsNames=useAsNames,
                      checkFiles=checkFiles#,readAcqus=readAcqus
                      )
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
#' @param checkFiles Only check if the folder exists or contains NMR data. Defaults to FALSE
#' @return An (invisible) list containing spectral data and the spectrum name
#' @export
#' @examples
#' exampleData<-readBruker(folder=system.file("extdata/1/10/pdata/10",package="mrbin"),
#'                         dimension="1D")

readBruker<-function(folder=NULL,dimension=NULL,onlyTitles=FALSE,
  useAsNames="Spectrum titles",checkFiles=FALSE#,readAcqus=FALSE
  ){#Read Bruker NMR spectral data
 AcquPars<-list(NS=0,BF1=0,P1=0,RG=0,PULPROG="",SOLVENT="")
 datanameDict<-c("1r","2rr")
 names(datanameDict)<-c("1D","2D")
 spectrum_proc_path<-folder
 datanameTmp<-datanameDict[dimension]
 if(!is.null(spectrum_proc_path)){
   BYTORDP_Dict<-c("little","big")
   names(BYTORDP_Dict)<-c(0,1)
   TITLE<-""
   try(TITLE<-scan(file=paste(spectrum_proc_path,"/title",sep=""),what="character",sep="\n",quiet=TRUE)[1],
     silent=TRUE)
   if(checkFiles){
     currentSpectrum <- NULL
     titleFinal <- NULL
     currentSpectrumTitle <- NULL
     currentSpectrumFolderName <- NULL
     currentSpectrumEXPNO <- NULL
     currentSpectrumFolderName_EXPNO <- NULL
     list.filesTMP<-list.files(spectrum_proc_path)
     if(!"procs"%in%list.filesTMP){
        stop(paste("Could not open spectrum",spectrum_proc_path))
     }
   }  else {
     if(!onlyTitles){
       #if(readAcqus){
       #try reading acquisition parameters for reporting
         spectrum_acqu_path<-NULL
         spectrum_acqu_pathTMP<-strsplit(spectrum_proc_path,"/")[[1]]
         spectrum_acqu_path<-paste(spectrum_acqu_pathTMP[1:
           (length(spectrum_acqu_pathTMP)-2)],sep="/",collapse="/")
         acqusTMP<-NULL
         try(acqusTMP<-scan(file=paste(spectrum_acqu_path,"/acqu",
           sep=""),what="character",sep="\n",quiet=TRUE),
           silent=TRUE)
         if(!(is.null(acqusTMP)|(length(acqusTMP)==0))){
           acqusTMP<-gsub("#","",acqusTMP)
           acqusTMP<-gsub("\\$"," ",acqusTMP)
           RGpos<-grep(" RG=",acqusTMP)
           if(length(RGpos)>0) AcquPars$RG<-as.numeric(strsplit(
             acqusTMP[RGpos],split="= ")[[1]][2])#receiver gain
           NSpos<-grep(" NS=",acqusTMP)
           if(length(NSpos)>0) AcquPars$NS<-as.numeric(strsplit(
             acqusTMP[NSpos],split="= ")[[1]][2])#NS number of scans
           BF1pos<-grep(" BF1=",acqusTMP)
           if(length(BF1pos)>0) AcquPars$BF1<-as.numeric(strsplit(
             acqusTMP[BF1pos],split="= ")[[1]][2])#BF1 Specrometer frequency
           PULPROGpos<-grep(" PULPROG=",acqusTMP)
           if(length(PULPROGpos)>0) AcquPars$PULPROG<-strsplit(
             acqusTMP[PULPROGpos],split="= ")[[1]][2]#pulse program name
           SOLVENTpos<-grep(" SOLVENT=",acqusTMP)
           if(length(SOLVENTpos)>0) AcquPars$SOLVENT<-strsplit(
             acqusTMP[SOLVENTpos],split="= ")[[1]][2]#solvent name
           P1pos<-grep(" P=",acqusTMP)
           if(length(P1pos)>0) AcquPars$P1<-as.numeric(#substr(
             strsplit(
             acqusTMP[P1pos+1],
             ,split=" ")[[1]][2]#entry 2 is P1
             #1,5)
             )#90 degree pulse length #$P= (0..63)\n 13.33
         }
       #}
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
   }
   invisible(list(currentSpectrum=currentSpectrum,currentSpectrumName=titleFinal,
           currentSpectrumTitle=currentSpectrumTitle,
           currentSpectrumFolderName=currentSpectrumFolderName,
           currentSpectrumEXPNO=currentSpectrumEXPNO,
           currentSpectrumFolderName_EXPNO=currentSpectrumFolderName_EXPNO,
           AcquPars=AcquPars))
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
#' mrbin(silent=TRUE,setDefault=TRUE,parameters=list(dimension="1D",tryParallel=FALSE,
#'          referenceScaling="No",binwidth1D=0.05,PQNScaling="No",PCA="No",
#'          NMRfolders=system.file("extdata/1/10/pdata/10",package="mrbin")))
#' referenceScaling()
#' }

referenceScaling<-function(NMRdata=NULL,reference1D=NULL,reference2D=NULL,dimension="1D"){
  dataProvided<-TRUE
  if(dimension=="1D"){
       scalingFactor<-mean(NMRdata[which(as.numeric(names(NMRdata))<=
                           reference1D[1]&as.numeric(names(NMRdata))>
                                       reference1D[2])],na.rm=TRUE)*
                                       abs(reference1D[1]-reference1D[2])
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
       scalingFactor<-mean(NMRdata[selectedTMP1,selectedTMP2],na.rm=TRUE)*
                             abs(reference2D[1]-reference2D[2])*
                             abs(reference2D[3]-reference2D[4])
       scaledSpectrum<-NMRdata/scalingFactor
  }
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
           #               as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrum))>mrbin.env$mrbin$parameters$removeAreaList[i,4]|
           #               as.numeric(rownames(mrbin.env$mrbinTMP$currentSpectrum))<mrbin.env$mrbin$parameters$removeAreaList[i,3]
           #               ),which(
           #               as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrum))>mrbin.env$mrbin$parameters$removeAreaList[i,1]|
           #               as.numeric(colnames(mrbin.env$mrbinTMP$currentSpectrum))<mrbin.env$mrbin$parameters$removeAreaList[i,2]
           #               )]
       }
     }
    invisible(newNMRdata)
  }
}

#' A function for creating mrbin objects.
#'
#' This function creates an mrbin object and returns it.
#' @return An (invisible) mrbin object
#' @export
#' @examples
#' mrbinObject<-createmrbin()

createmrbin<-function(){
  mrbinObject<-structure(list(
    bins=NULL,#matrix of bin data, rows=samples, columns=features
    parameters=list(
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
               dilutionCorrection="No",
               PQNScaling="No",
               fixNegatives="No",
               logTrafo="No",
               unitVarianceScaling="No",
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
               useMeanIntensityForBins=FALSE,
               mrbinversionTracker=as.character(utils::packageVersion("mrbin")),
               previewRegion1D=c(1.5,1.15,16.5,27),
               previewRegion2D=c(2,1,16.5,27),
               noisePreviewRegion1D=rbind(c(8.1,7.5,125,135),c(3.3,2.7,38,50),c(1.1,.5,10,27)),
               noisePreviewRegion2D=rbind(c(8.0,7.2,123,137),c(3.7,2.9,38,50),c(1.6,.8,9,28)),
               maxPreviewPlots2D=2,
               saveFiles="No",
               verbose=TRUE,
               outputFileName=NULL,
               NMRfolders=NULL,
               medians=NULL,
               noise_level_Raw=NULL,
               noise_level_adjusted=NULL,
               noise_level=NULL,#noise levels for each bin of each spectrum
               binRegions=matrix(ncol=4,nrow=0,dimnames=list(NULL,c("left","right","top","bottom"))),
               AcquPars=list(NS=0,BF1=0,P1=0,RG=0,PULPROG="",SOLVENT=""),
               numberOfFeaturesRaw=NULL,
               numberOfFeaturesAfterRemovingSolvent=NULL,
               numberOfFeaturesAfterRemovingAreas=NULL,
               numberOfFeaturesAfterSummingBins=NULL,
               numberOfFeaturesAfterTrimmingZeros=NULL,
               numberOfFeaturesAfterNoiseRemoval=NULL,
               numberOfFeaturesAfterCropping=NULL,
               tryParallel=TRUE,
               createCode="",#Code for recreating this dataset can be stored here
               warningMessages=NULL,
               errorsAsWarnings=FALSE,
               Factors=NULL,#for backward conpatibility
               defineGroups=NULL#to be removed in future versions
    ),
    metadata=list(
               projectTitle="",#Project title
               projectIdentifier="",#Project identifier
               projectDescription="",#Project description
               projectAuthors="",#Project author list
               projectApprovals="",#Project approval numbers, e.g. IRB protocol number
               sampleType="",#e.g. serum, plasma, urine, blood, CSF, tissue extract, etc.
               organism="",#e.g. human, mouse, rat, E.coli, etc.
               solvent="",#e.g. H2O, D2O, CDCl3, etc.
               pH="",#e.g. 7.4
               dilutionFactors=NULL,#vector of dilution correction values
               annotations="",#Vector of features with potential metabolite identities, corresponds to columns in $bins
               factors=NULL,
               metaData=NULL,#for storing matrix of metadata
               metaboliteIdentities=matrix(nrow=0,ncol=4,dimnames=list(NULL,c('left','right','top','bottom')))#Potential metabolite identities and ppm values. data.frame with metaboliteName, left, right, top, bottom
    ),
    transformations=NULL,#data transformation or scaling is documented here for later review, including: reference scaling, log trafo, PQN, atnv, etc.
    changeLog=NULL,#changes to data, parameters, or metadata are documented here, including: time, function name, version, names of changed parameters
    changeValues=list(NULL)
  ), class = "mrbin")

  mrbinObject<-timeStampMrbin(mrbinObject,steps=0)
  mrbinObject<-timeStampMrbin(mrbinObject,
    functionName="mrbin::createmrbin",
    versionNumber=as.character(utils::packageVersion("mrbin")),
    changeDetails="Created mrbin object",
    steps=1
  )

  invisible(mrbinObject)
}

#' A function for time stamping mrbin objects.
#'
#' This function adds time stamps to an mrbin object and returns it. Is used only within functions making changes to mrbin objects.
#' @param mrbinObject An mrbin object
#' @param functionName Name of the package and function calling this command
#' @param versionNumber Version number of the package calling this command
#' @param changeDetails Details of changes made to the mrbin object
#' @param steps Indicates which step to perform: 0 (only pre-change), 1 (only post-change)
#' @param comment An optional character vector describing the change
#' @return An (invisible) mrbin object
#' @export
#' @examples
#' mrbinObject<-createmrbin()
#' mrbinObject<-timeStampMrbin(mrbinObject)

timeStampMrbin<-function(mrbinObject,functionName="InProgress...",versionNumber="0",
  changeDetails="InProgress...",steps=0,comment=""){
  rowTMP<-1
  changeValues<-vector("list",0)
  if(length(nrow(mrbinObject$changeLog))>0){
    rowTMP<-nrow(mrbinObject$changeLog)+1
  }
  if(steps==1){
    rowTMP<-nrow(mrbinObject$changeLog)
    changeValues <- mrbinObject$changeValues[[rowTMP]]
  }
  if(!is.null(mrbinObject$bins)) {
    checkValue1<-(sum(mrbinObject$bins*(1:length(mrbinObject$bins)),na.rm=TRUE)-
        sum(mrbinObject$bins,na.rm=TRUE)+
        sum(nchar(rownames(mrbinObject$bins))*(1:nrow(mrbinObject$bins)),na.rm=TRUE)+
        sum(nchar(colnames(mrbinObject$bins))*(1:ncol(mrbinObject$bins)),na.rm=TRUE)-
        sum(nchar(rownames(mrbinObject$bins)),na.rm=TRUE)-
        sum(nchar(colnames(mrbinObject$bins)),na.rm=TRUE))*cos(rowTMP)
  } else {
    checkValue1<-0
  }
  changeValues[["bins"]][1]<-checkValue1
  if(steps==0) changeValues[["bins"]][2]<-checkValue1  #parameters old
  checkValue3<-length(mrbinObject$parameters)*cos(rowTMP)
  changeValues[["numberOfParameters"]][1]<-checkValue3
  if(steps==0) changeValues[["numberOfParameters"]][2]<-checkValue3  #parameters old
  if(!is.null(mrbinObject$parameters)) {
    for(i in 1:length(mrbinObject$parameters)){
        checkValue3<-length(mrbinObject$parameters[[i]])
        if(is.character(mrbinObject$parameters[[i]])){
          checkValue3<-checkValue3+sum(nchar(mrbinObject$parameters[[i]],keepNA=FALSE),na.rm=TRUE)
        }
        if(is.numeric(mrbinObject$parameters[[i]])){
          checkValue3<-checkValue3+sum(mrbinObject$parameters[[i]],na.rm=TRUE)
        }
        if(!is.null(names(mrbinObject$parameters[[i]]))){
          checkValue3<-checkValue3+sum(nchar(names(mrbinObject$parameters[[i]]),keepNA=FALSE),na.rm=TRUE)
        }
        if(!is.null(colnames(mrbinObject$parameters[[i]]))){
          checkValue3<-checkValue3+sum(nchar(colnames(mrbinObject$parameters[[i]]),keepNA=FALSE),na.rm=TRUE)
        }
        if(!is.null(rownames(mrbinObject$parameters[[i]]))){
          checkValue3<-checkValue3+sum(nchar(rownames(mrbinObject$parameters[[i]]),keepNA=FALSE),na.rm=TRUE)
        }
      checkValue3<-checkValue3*cos(rowTMP)
      changeValues[[names(mrbinObject$parameters)[i]]][1]<-checkValue3
      if(steps==0) changeValues[[names(mrbinObject$parameters)[i]]][2]<-checkValue3 #parameters old
    }
  }
  if(!is.null(mrbinObject$metadata)) {
    for(i in 1:length(mrbinObject$metadata)){
        checkValue3<-length(mrbinObject$metadata[[i]])
        if(is.character(mrbinObject$metadata[[i]])){
          checkValue3<-checkValue3+sum(nchar(mrbinObject$metadata[[i]],keepNA=FALSE),na.rm=TRUE)
        }
        if(is.numeric(mrbinObject$metadata[[i]])){
          checkValue3<-checkValue3+sum(mrbinObject$metadata[[i]],na.rm=TRUE)
        }
        if(!is.null(names(mrbinObject$metadata[[i]]))){
          checkValue3<-checkValue3+sum(nchar(names(mrbinObject$metadata[[i]]),keepNA=FALSE),na.rm=TRUE)
        }
        if(!is.null(colnames(mrbinObject$metadata[[i]]))){
          checkValue3<-checkValue3+sum(nchar(colnames(mrbinObject$metadata[[i]]),keepNA=FALSE),na.rm=TRUE)
        }
        if(!is.null(rownames(mrbinObject$metadata[[i]]))){
          checkValue3<-checkValue3+sum(nchar(rownames(mrbinObject$metadata[[i]]),keepNA=FALSE),na.rm=TRUE)
        }
      checkValue3<-checkValue3*cos(rowTMP)
      changeValues[[names(mrbinObject$metadata)[i]]][1]<-checkValue3
      if(steps==0) changeValues[[names(mrbinObject$metadata)[i]]][2]<-checkValue3 #parameters old
    }
  }
  #}
  #message(paste(checkValue1,checkValue2,checkValue3,checkValue4))
  if(steps==0){
    mrbinObject$changeLog<-rbind(mrbinObject$changeLog,
      data.frame(
        Time1=Sys.time(),
        Time2=Sys.time(),
        Comment=comment,
        Function=functionName,#"mrbin::createmrbin",
        Version=versionNumber,
        System=paste0(R.Version()$language,R.Version()$major,".",R.Version()$minor,",",R.Version()$platform),
        Change=changeDetails,
        stringsAsFactors=FALSE
      )
    )
    mrbinObject$changeValues[[rowTMP]]<-changeValues
  }
  if(steps==1){
     mrbinObject$changeLog[rowTMP,"Time2"]<-Sys.time()
     mrbinObject$changeLog[rowTMP,"Comment"]<-comment
     mrbinObject$changeLog[rowTMP,"Function"]<-functionName
     mrbinObject$changeLog[rowTMP,"Version"]<-versionNumber
     mrbinObject$changeLog[rowTMP,"System"]<-paste0(R.Version()$language,R.Version()$major,".",
       R.Version()$minor,",",R.Version()$platform)
     mrbinObject$changeLog[rowTMP,"Change"]<-changeDetails
     mrbinObject$changeValues[[rowTMP]]<-changeValues
  }
  invisible(mrbinObject)
}

#' A function for checking mrbin objects.
#'
#' This function checks an mrbin object and returns warning if changes were not documented
#' @param mrbinObject An mrbin object
#' @param verbose Should a summary be displayed? (Warnings will be displayed even when setting verbose to FALSE)
#' @param errorsAsWarnings If TRUE, errors will be turned into warnings. Should be used with care, as errors indicate undocumented changes to the data. If not provided, this will be taken from the mrbinObject.
#' @return An (invisible) character vector of warnings
#' @export
#' @examples
#' mrbinObject<-createmrbin()
#' mrbinObject<-checkmrbin(mrbinObject)

checkmrbin<-function(mrbinObject,verbose=TRUE,errorsAsWarnings=NULL){
 if(methods::is(mrbinObject,"mrbin")){
  if(is.null(errorsAsWarnings)) errorsAsWarnings<-mrbinObject$parameters$errorsAsWarnings
  mrbinObject<-timeStampMrbin(mrbinObject,steps=0)
  mrbinObject<-timeStampMrbin(mrbinObject,steps=1)
  warningMessages<-NULL
  if(nrow(mrbinObject$changeLog)>1){
    if(!identical(nrow(mrbinObject$changeLog),length(mrbinObject$changeValues))){
      warningMessages<-c(warningMessages,
        "Change logs have been altered. Data and parameters might not be accurate.")
    } else {
      for(i in 2:nrow(mrbinObject$changeLog)){
          namesTMP<-NULL
          for(j in names(mrbinObject$changeValues[[i]])){
           if(j %in% names(mrbinObject$changeValues[[i-1]])){
             if(!identical(#round values to avoid numeric instabilities
                 signif(mrbinObject$changeValues[[i]][[j]][2],6),
                 signif(mrbinObject$changeValues[[i-1]][[j]][1]*cos(i)/cos(i-1),6))
                 ){
                 namesTMP<-c(namesTMP,j)
             }
           } else {
             namesTMP<-c(namesTMP,j)
           }
         }
         if(!is.null(namesTMP)){
           warningMessages<-c(warningMessages,paste("Unstated changes to: ",
               paste(namesTMP,sep=", ",collapse=", ")," between ",
               mrbinObject$changeLog[i-1,"Time2"]," and ",mrbinObject$changeLog[i,"Time1"],
               sep="",collapse=""))
         }
      }
    }
  }
  if(verbose){
    message("Data was processed as follows:\n",paste(mrbinObject$transformations,
               sep=", ",collapse=", "))
  }
  if(!is.null(warningMessages)){
    #for(i in 1:length(warningMessages)){
     if(!errorsAsWarnings){
       stop(paste(warningMessages,sep="\n",collapse="\n"))#[i])
     } else {
       warning(paste(warningMessages,sep="\n",collapse="\n"))#[i])
     }
    #}
  } else {
   if(verbose){
     message("Data seems to be valid. All changes have been documented.")
   }
  }
  invisible(warningMessages)
 } else {
  stop("Not an mrbin object.")
 }
}

#' A function for editing metabolite identities.
#'
#' This function edits the metabolite list within an mrbin object and returns it
#' @param mrbinObject An mrbin object
#' @param borders A matrix of signal borders. 1D: two columns: left, right 2D: four columns: left, right, top, bottom
#' @param metabolitenames A character vector of metabolite identities
#' @param add Should the new metabolite list be added to an existing list, or replace the current list?
#' @return An (invisible) mrbin object
#' @export
#' @examples
#'  mrbinObject<-mrbin(silent=TRUE,
#'                    parameters=list(verbose=TRUE,dimension="1D",PQNScaling="No",
#'                    binwidth1D=0.04,signal_to_noise1D=1,PCA="No",binRegion=c(9.5,0.5,10,156),
#'                    saveFiles="No",referenceScaling="No",noiseRemoval="No",
#'                    fixNegatives="No",logTrafo="No",noiseThreshold=.05,tryParallel=FALSE,
#'                    NMRfolders=c(system.file("extdata/2/10/pdata/10",package="mrbin"),
#'                               system.file("extdata/3/10/pdata/10",package="mrbin"))
#'                    ))
#'  mrbinObject<-editmetabolitesmrbin(mrbinObject,borders=matrix(c(
#'      1.346,1.324,
#'      3.052,3.043,
#'      4.066,4.059
#'    ),ncol=2,byrow=TRUE),metabolitenames=c(
#'    "Lactate",
#'    "Creatinine",
#'    "Creatinine"
#'    ))
#' mrbinObject$parameters$metaboliteIdentities

editmetabolitesmrbin<-function(mrbinObject,borders,metabolitenames,add=TRUE){ #Define metabolite names
  if(!nrow(borders)==length(metabolitenames)){
    stop("Metabolite names do not match number of rows in borders matrix.")
  }
  if(ncol(borders)==2){
     borders<-cbind(borders,rep(-10,nrow(borders)),rep(200,nrow(borders)))
  }
  colnames(borders)<-c('left','right','top','bottom')
  metaboliteIdentities<-mrbinObject$metadata$metaboliteIdentities
  if(!add){
    metaboliteIdentities<-NULL
  }
  rownames(borders)<-metabolitenames
  metaboliteIdentities=rbind(metaboliteIdentities,borders)
  mrbinObject<-editmrbin(mrbinObject,metadata=list(
    metaboliteIdentities=metaboliteIdentities),
    functionName="mrbin::editmetabolitesmrbin",
    versionNumber=as.character(utils::packageVersion("mrbin")),
    verbose=FALSE)
  mrbinObject<-annotatemrbin(mrbinObject)
  invisible(mrbinObject)
}


#' A function for annotating mrbin objects.
#'
#' This function annotates an mrbin object and returns it with updated $annotations vector
#' @param mrbinObject An mrbin object
#' @return An (invisible) mrbin object
#' @export
#' @examples
#'  mrbinObject<-mrbin(silent=TRUE,
#'                    parameters=list(verbose=TRUE,dimension="1D",PQNScaling="No",
#'                    binwidth1D=0.04,signal_to_noise1D=1,PCA="No",binRegion=c(9.5,0.5,10,156),
#'                    saveFiles="No",referenceScaling="No",noiseRemoval="No",
#'                    fixNegatives="No",logTrafo="No",noiseThreshold=.05,tryParallel=FALSE,
#'                    NMRfolders=c(system.file("extdata/2/10/pdata/10",package="mrbin"),
#'                               system.file("extdata/3/10/pdata/10",package="mrbin"))),
#'                    metadata=list(metaboliteIdentities=matrix(c(
#'                    1.346,1.324,21,23,
#'                    3.052,3.043,30.5,33.5,
#'                    4.066,4.059,57,59.5
#'                    ),ncol=4,byrow=TRUE,dimnames=list(
#'                    c("Lactate","Creatinine","Creatinine"),NULL))))
#' mrbinObject<-annotatemrbin(mrbinObject)
#' mrbinObject$metadata$annotations[125:135]
#' plotPCA(mrbinObject,loadings=TRUE)

annotatemrbin<-function(mrbinObject){ #Define metabolite names
  if(nrow(mrbinObject$metadata$metaboliteIdentities)>0){
    annotations<-rep("",nrow(mrbinObject$parameters$binRegions))
    names(annotations)<-colnames(mrbinObject$bins)
    if(mrbinObject$parameters$dimension=="1D"){
      #if(nrow(mrbinObject$metadata$metaboliteIdentities)>0){
        for(i in 1:nrow(mrbinObject$metadata$metaboliteIdentities)){
          testTMP<-(mrbinObject$parameters$binRegions[,1]<=mrbinObject$metadata$metaboliteIdentities[i,1]&
                    mrbinObject$parameters$binRegions[,1]>=mrbinObject$metadata$metaboliteIdentities[i,2])|
                   (mrbinObject$parameters$binRegions[,2]<=mrbinObject$metadata$metaboliteIdentities[i,1]&
                    mrbinObject$parameters$binRegions[,2]>=mrbinObject$metadata$metaboliteIdentities[i,2])|
                   (mrbinObject$parameters$binRegions[,1]>=mrbinObject$metadata$metaboliteIdentities[i,1]&
                    mrbinObject$parameters$binRegions[,2]<=mrbinObject$metadata$metaboliteIdentities[i,2])
          if(sum(testTMP)>0){
           metNameTMP<-rownames(mrbinObject$metadata$metaboliteIdentities)[i]
           for(j in which(testTMP)){
            if(annotations[j]==""){
              annotations[j]<-metNameTMP
            } else {
              if(!(annotations[j]==metNameTMP|
                 substr(annotations[j],1,nchar(metNameTMP)+2)==paste(metNameTMP,", ",sep="")|
                 substr(annotations[j],max(1,nchar(annotations[j])-nchar(metNameTMP)-2),nchar(annotations[j]))==
                 paste(", ",metNameTMP,sep="")|
                 grepl(paste(", ",metNameTMP,", ",sep=""),annotations[j]))
                 ){
                 annotations[j]<-paste(annotations[j],metNameTMP,sep=", ")
              }
            }
           }
          }
        }
      #}
    }
    if(mrbinObject$parameters$dimension=="2D"){
      #if(nrow(mrbinObject$metadata$metaboliteIdentities)>0){
        for(i in 1:nrow(mrbinObject$metadata$metaboliteIdentities)){
          testTMP<-((mrbinObject$parameters$binRegions[,1]<=mrbinObject$metadata$metaboliteIdentities[i,1]&
                    mrbinObject$parameters$binRegions[,1]>=mrbinObject$metadata$metaboliteIdentities[i,2])|
                   (mrbinObject$parameters$binRegions[,2]<=mrbinObject$metadata$metaboliteIdentities[i,1]&
                    mrbinObject$parameters$binRegions[,2]>=mrbinObject$metadata$metaboliteIdentities[i,2])|
                   (mrbinObject$parameters$binRegions[,1]>=mrbinObject$metadata$metaboliteIdentities[i,1]&
                    mrbinObject$parameters$binRegions[,2]<=mrbinObject$metadata$metaboliteIdentities[i,2]))&#2D
                   ((mrbinObject$parameters$binRegions[,4]<=mrbinObject$metadata$metaboliteIdentities[i,4]&
                    mrbinObject$parameters$binRegions[,4]>=mrbinObject$metadata$metaboliteIdentities[i,3])|
                   (mrbinObject$parameters$binRegions[,3]<=mrbinObject$metadata$metaboliteIdentities[i,4]&
                    mrbinObject$parameters$binRegions[,3]>=mrbinObject$metadata$metaboliteIdentities[i,3])|
                   (mrbinObject$parameters$binRegions[,4]>=mrbinObject$metadata$metaboliteIdentities[i,4]&
                    mrbinObject$parameters$binRegions[,3]<=mrbinObject$metadata$metaboliteIdentities[i,3]))
          if(sum(testTMP)>0){
           metNameTMP<-rownames(mrbinObject$metadata$metaboliteIdentities)[i]
           for(j in which(testTMP)){
            if(annotations[j]==""){
              annotations[j]<-metNameTMP
            } else {
              if(!(annotations[j]==metNameTMP|
                 substr(annotations[j],1,nchar(metNameTMP)+2)==paste(metNameTMP,", ",sep="")|
                 substr(annotations[j],max(1,nchar(annotations[j])-nchar(metNameTMP)-2),nchar(annotations[j]))==
                 paste(", ",metNameTMP,sep="")|grepl(paste(", ",metNameTMP,", ",sep=""),annotations[j]))
                 ){
                 annotations[j]<-paste(annotations[j],metNameTMP,sep=", ")
              }
            }
           }
          }
        }
      #}
    }
    if(sum(annotations=="")>0) annotations[annotations==""]<-colnames(mrbinObject$bins)[annotations==""]
    mrbinObject<-editmrbin(mrbinObject,functionName="mrbin::annotatemrbin",
             versionNumber=as.character(utils::packageVersion("mrbin")),
      metadata=list(annotations=trimws(annotations)),verbose=FALSE)
  }
  invisible(mrbinObject)
}


#' A function for editing mrbin objects.
#'
#' This function edits an mrbin object and returns it. This is the only documented way to edit mrbin objects, all other ways of editing such object might cause warning message
#' @param mrbinObject An mrbin object
#' @param functionName Name of the package and function calling this command
#' @param versionNumber Version number of the package calling this command
#' @param bins A matrix containing values to be written to the mrbin object
#' @param parameters A list containing values to be written to the mrbin object parameters, names must be names of the mrbin object, e.g. dimension
#' @param metadata A list containing values to be written to the mrbin object parameters, names must be names of the mrbin object
#' @param transformations An optional character vector describing any used data transformations or scaling such as reference scaling, PQN, log, atnv, etc.
#' @param comment An optional character vector describing the change
#' @param verbose Should a summary be displayed?
#' @return An (invisible) mrbin object
#' @export
#' @examples
#' mrbinObject<-createmrbin()
#' mrbinObject<-editmrbin(mrbinObject)

editmrbin<-function(mrbinObject,functionName="mrbin::editmrbin",
  versionNumber=as.character(utils::packageVersion("mrbin")),
  bins=NULL,parameters=NULL,metadata=NULL,transformations=NULL,comment="",
  verbose=TRUE){

  mrbinObject<-timeStampMrbin(mrbinObject,steps=0)
  #changeDetails<-""
  parametersTMP<-NULL
  if(!is.null(bins)){
    if(!identical(bins,mrbinObject$bins)){
      if(is.null(transformations)){
        message("Changes to bin data should be accompanied by a brief explanation in the parameter transformations")
      }
      mrbinObject$bins<-bins
      parametersTMP<-c(parametersTMP,"bins")
      if(is.null(transformations)){
        commentTMP<-comment#if bin data is changed, add this to $transformations
        if(comment=="") commentTMP<-paste("Unknown change in ",functionName,sep="")
        transformations<-c(transformations,commentTMP)
      }
    } else {
       transformations<-NULL#If bins did not change, do not save the transformation info
    }
  } else {
     transformations<-NULL#If bins did not change, do not save the transformation info
  }

  if(!is.null(transformations)){
    mrbinObject$transformations<-c(mrbinObject$transformations,transformations)
    parametersTMP<-c(parametersTMP,"transformations")
  }
  if(!is.null(parameters)){
    for(i in 1:length(parameters)){
      if(!identical(parameters[[i]],mrbinObject$parameters[[names(parameters)[i]]])){
        mrbinObject$parameters[[names(parameters)[i]]]<-parameters[[i]]
        parametersTMP<-c(parametersTMP,names(parameters)[i])
      }
    }
  }
  if(!is.null(metadata)){
    for(i in 1:length(metadata)){
      if(!identical(metadata[[i]],mrbinObject$metadata[[names(metadata)[i]]])){
        mrbinObject$metadata[[names(metadata)[i]]]<-metadata[[i]]
        parametersTMP<-c(parametersTMP,names(metadata)[i])
      }
    }
  }
  if(!is.null(parametersTMP)){
    changeDetails=paste("Changes: ",paste(parametersTMP,sep=", ",collapse=", "),
      sep="",collapse="")
  } else {
    changeDetails<-"No changes."
  }
  mrbinObject2<-timeStampMrbin(mrbinObject,
    functionName=functionName,
    versionNumber=versionNumber,
    changeDetails=changeDetails,
    steps=1,
    comment=comment
  )
  mrbin.env$mrbinObject2<-mrbinObject2
  checkmrbin(mrbinObject2,verbose=verbose)
  invisible(mrbinObject2)
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
  binRegions=NULL,binMethod="Rectangular bins",
  useMeanIntensityForBins=FALSE,spectrumTitle=NULL){
  numberOfPointsPerBin<-NULL
  binTMP<-rep(0,nrow(binRegions))
  names(binTMP)<-rownames(binRegions)
  if(dimension=="2D"){#2d spectra
    #Create index of signals in each bin
    NMRspectrumRownames<-as.numeric(rownames(currentSpectrum))
    NMRspectrumColnames<-as.numeric(colnames(currentSpectrum))
    for(ibinTMP in 1:nrow(binRegions)){
      rowsTMP<-NMRspectrumRownames<=binRegions[ibinTMP,4]&
                            NMRspectrumRownames>binRegions[ibinTMP,3]
      colsTMP<-NMRspectrumColnames<=binRegions[ibinTMP,1]&
                    NMRspectrumColnames>binRegions[ibinTMP,2]
      numberOfPointsPerBinTMP<-(sum(rowsTMP)*sum(colsTMP))-sum(is.na(
        currentSpectrum[rowsTMP,colsTMP]))
      numberOfPointsPerBin<-c(numberOfPointsPerBin,numberOfPointsPerBinTMP)
      if(numberOfPointsPerBinTMP>0){
         adjustTMP<-1
         if(!useMeanIntensityForBins){
           adjustTMP<-abs(binRegions[ibinTMP,4]-binRegions[ibinTMP,3])*
                             abs(binRegions[ibinTMP,2]-binRegions[ibinTMP,1])
         }
        binTMP[ibinTMP]<-(sum(currentSpectrum[rowsTMP,colsTMP],na.rm=TRUE)/
                         numberOfPointsPerBinTMP)*adjustTMP
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
           adjustTMP<-1
           if(!useMeanIntensityForBins){
             adjustTMP<-abs(binRegions[ibinTMP,2]-binRegions[ibinTMP,1])
           }
           binTMP[ibinTMP]<-(sum(currentSpectrum[
                            indexTMP])/numberOfPointsPerBinTMP)*adjustTMP
           #Set to NA: Make sure each point is counted only once for rectangular bins. For custom bins lists, double counting may be on purpose
           #if(binMethod=="Rectangular bins") currentSpectrum[NMRspectrumNames<=binRegions[ibinTMP,1]&
           #                        NMRspectrumNames>binRegions[ibinTMP,2]]<-NA
        }
     }
   }
   if(FALSE){
     #this code would avoid a time-consuming loop, but it uses too much memory
     #in parallel mode. Without parallel, it might save only little time
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
  }
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
   noiseRange1d=NULL,noiseRange2d=NULL,binRegions=NULL,
   useMeanIntensityForBins=FALSE){
  if(!is.null(NMRdata)){
    if(dimension=="1D"){
         adjustTMP<-1
         if(!useMeanIntensityForBins){
           adjustTMP<-abs(binRegions[,2]-binRegions[,1])
           #message(adjustTMP)
         }
         noise_level<-stats::sd(NMRdata[
                 which(as.numeric(names(NMRdata))<=max(noiseRange1d[1:2])&
                      as.numeric(names(NMRdata))>=min(noiseRange1d[1:2]))])
         noise_level_TMP<-noise_level*(pointsPerBin^(-.5))*adjustTMP
    }
    if(dimension=="2D"){
         adjustTMP<-1
         if(!useMeanIntensityForBins){
           adjustTMP<-abs(binRegions[,2]-binRegions[,1])*abs(binRegions[,4]-binRegions[,3])
         }
         noise_level<-stats::sd(NMRdata[
               which(as.numeric(rownames(NMRdata))>=min(noiseRange2d[3:4])&
                 as.numeric(rownames(NMRdata))<=max(noiseRange2d[3:4])),
               which(as.numeric(colnames(NMRdata))<=max(noiseRange2d[1:2])&
                 as.numeric(colnames(NMRdata))>=min(noiseRange2d[1:2]))])
         #correct noise for size of noise range:
         #noise_level<-
         noise_level_TMP<-noise_level*(pointsPerBin^(-.5))*adjustTMP
    }
    if(sum(is.infinite(noise_level_TMP))>0){
       noise_level_TMP[is.infinite(noise_level_TMP)]<-NaN
    }
    if(sum(is.nan(noise_level_TMP))>0){
       noise_level_TMP[is.nan(noise_level_TMP)]<-median(noise_level_TMP,na.rm=TRUE)#Inf
    }
    invisible(list(noise_level=noise_level,#noise level per sample before reference scaling
      noise_level_TMP=noise_level_TMP))#noise level per bin and sample before reference scaling
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
#' @param maxPlots The maximum number of 2D plots to be overlaid
#' @param color Defines the color of the spectrum plot. If NULL, a rainbow theme is used for 2D NMR
#' @param add If TRUE, additional spectrum plots are overlaid with the current plot
#' @param manualScale If TRUE, scaling factor is taken from environment variables
#' @param plotTitle Defines the main title of the plot
#' @param showGrid Shows a grid of data points. Defaults to FALSE
#' @param restrictToRange Restrict plot area to range of available data points. Defaults to FALSE
#' @param enableSplit Allow split plots for showing 1D and 2D spectra simultaneously
#' @param dimension If not provided, this will be taken from package environment
#' @return {None}
#' @export
#' @examples
#' mrbin(silent=TRUE,setDefault=TRUE,parameters=list(dimension="1D",binwidth1D=.1,
#'          PQNScaling="No",noiseRemoval="No",trimZeros="No",tryParallel=FALSE,
#'          fixNegatives="No",logTrafo="No",PCA="No",verbose=TRUE,
#'          NMRfolders=system.file("extdata/1/10/pdata/10",package="mrbin")))
#' plotMultiNMR()

plotMultiNMR<-function(region=NULL,rectangleRegions=NULL,
                   rectangleColors=c("green","orange","blue","red","yellow",
                     "gray","purple"),
                   rectangleFront=FALSE,
                   polygonRegion=NULL,maxPlots=Inf,
                   color=NULL,add=FALSE,showGrid=FALSE,
                   manualScale=TRUE,plotTitle="",
                   restrictToRange=FALSE,enableSplit=TRUE,dimension=NULL){
  if(is.null(dimension)) dimension<-mrbin.env$mrbin$parameters$dimension
  #if 1D and 2D spectra are present, create a plot with two rows
  TMP1D<-FALSE
  TMP2D<-FALSE
  if(dimension=="1D"#&!is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)
  ){
    TMP1D<-TRUE
  }
  if(dimension=="2D"#&!is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)
  ){
    TMP2D<-TRUE
  }
  if(enableSplit){
    if(!is.null(mrbin.env$mrbinTMP$additionalPlots1D)) TMP1D<-TRUE
    if(!is.null(mrbin.env$mrbinTMP$additionalPlots2D)) TMP2D<-TRUE
    if(TMP1D&TMP2D){
      par(mfrow=c(2,1),mar=c(2.3, 2.1, 0.5, 0.8))#bottom, left, top, and right
    }
  }
  if(TMP2D){
    counterTMP<-0
    #if(!is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
      plotNMR(region=region,rectangleRegions=rectangleRegions,
                   rectangleColors=rectangleColors,
                   rectangleFront=rectangleFront,
                   polygonRegion=polygonRegion,
                   color=color,add=add,showGrid=showGrid,
                   manualScale=manualScale,plotTitle=plotTitle,
                   restrictToRange=restrictToRange,dimension=dimension)
    if(!is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
      counterTMP<-counterTMP+1
    }
    if(!is.null(mrbin.env$mrbinTMP$additionalPlots2D)){
      for(iTMP in 1:min(maxPlots-1,length(mrbin.env$mrbinTMP$additionalPlots2D))){
        if(counterTMP==0){
          plotNMR(region=region,rectangleRegions=NULL,
                       rectangleColors=rectangleColors,
                       rectangleFront=FALSE,
                       polygonRegion=polygonRegion,
                       color=color,add=add,showGrid=FALSE,
                       manualScale=manualScale,plotTitle=plotTitle,
                       restrictToRange=restrictToRange,
                       currentSpectrumOriginal=
                       mrbin.env$mrbinTMP$additionalPlots2D[[iTMP]]*as.numeric(
                       mrbin.env$mrbinTMP$additionalPlots2DMetadata[iTMP,7])+as.numeric(
                       mrbin.env$mrbinTMP$additionalPlots2DMetadata[iTMP,9]),dimension=dimension)
          counterTMP<-counterTMP+1
        } else {
          plotNMR(region=region,rectangleRegions=NULL,
                       rectangleColors=rectangleColors,
                       rectangleFront=FALSE,
                       polygonRegion=polygonRegion,
                       color=mrbin.env$mrbinTMP$additionalPlots2DMetadata[iTMP,8],
                       add=TRUE,showGrid=FALSE,
                       manualScale=manualScale,plotTitle=plotTitle,
                       restrictToRange=restrictToRange,
                       currentSpectrumOriginal=
                       mrbin.env$mrbinTMP$additionalPlots2D[[iTMP]]*as.numeric(
                       mrbin.env$mrbinTMP$additionalPlots2DMetadata[iTMP,7])+as.numeric(
                       mrbin.env$mrbinTMP$additionalPlots2DMetadata[iTMP,9]),dimension=dimension)
          counterTMP<-counterTMP+1
        }
      }
    }
  }
  if(TMP1D){
    counterTMP<-0
    #if(!is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
      plotNMR(region=region,rectangleRegions=rectangleRegions,
                   rectangleColors=rectangleColors,
                   rectangleFront=rectangleFront,
                   polygonRegion=polygonRegion,
                   color=color,add=add,showGrid=showGrid,
                   manualScale=manualScale,plotTitle=plotTitle,
                   restrictToRange=restrictToRange,dimension=dimension)
    if(!is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
      counterTMP<-counterTMP+1
    }
    if(!is.null(mrbin.env$mrbinTMP$additionalPlots1D)){
      for(iTMP in 1:length(mrbin.env$mrbinTMP$additionalPlots1D)){
        if(counterTMP==0){
          plotNMR(region=region,rectangleRegions=NULL,
                       rectangleColors=rectangleColors,
                       rectangleFront=FALSE,
                       polygonRegion=polygonRegion,
                       color=color,add=add,showGrid=FALSE,
                       manualScale=manualScale,plotTitle=plotTitle,
                       restrictToRange=restrictToRange,
                       currentSpectrumOriginal=
                       mrbin.env$mrbinTMP$additionalPlots1D[[iTMP]]*as.numeric(
                       mrbin.env$mrbinTMP$additionalPlots1DMetadata[iTMP,7])+as.numeric(
                       mrbin.env$mrbinTMP$additionalPlots1DMetadata[iTMP,9]),dimension=dimension)
          counterTMP<-counterTMP+1
        } else {
          plotNMR(region=region,rectangleRegions=NULL,
                       rectangleColors=rectangleColors,
                       rectangleFront=FALSE,
                       polygonRegion=polygonRegion,
                       color=mrbin.env$mrbinTMP$additionalPlots1DMetadata[iTMP,8],
                       add=TRUE,showGrid=FALSE,
                       manualScale=manualScale,plotTitle=plotTitle,
                       restrictToRange=restrictToRange,
                       currentSpectrumOriginal=
                       mrbin.env$mrbinTMP$additionalPlots1D[[iTMP]]*as.numeric(
                       mrbin.env$mrbinTMP$additionalPlots1DMetadata[iTMP,7])+as.numeric(
                       mrbin.env$mrbinTMP$additionalPlots1DMetadata[iTMP,9]),dimension=dimension)
          counterTMP<-counterTMP+1
        }
      }
    }
  }
  gc()
}


#' A function for plotting NMR spectra.
#'
#' This function plots NMR spectra. A menu of commands is displayed to edit the
#' plot view and add spectra. Multiple spectra will be overlaid, and if both
#' 1D and 2D spectra are selected, they are shown in two plots with matched ranges.
#' @param hideMenu Do not show the menu. Defaults to FALSE
#' @param folders Optional vector of folder names of spectra to load. Defaults to NULL
#' @param dimensions Optional vector dimensions of spectra to load. Defaults to NULL
#' @param zoom Optional vector of initial zoom area. Defaults to NULL
#' @param intensity1D Optional value of initial 1D intensity. Defaults to NULL
#' @return {None}
#' @export
#' @examples
#' resetEnv()
#' mrplot(folders=c(system.file("extdata/1/12/pdata/10",package="mrbin"),
#'                  system.file("extdata/1/10/pdata/10",package="mrbin"),
#'                  system.file("extdata/2/10/pdata/10",package="mrbin"),
#'                  system.file("extdata/3/10/pdata/10",package="mrbin")),
#'        dimensions=c("2D","1D","1D","1D"),zoom=c(2.8,2.3,20,55),
#'        intensity1D=32,hideMenu=TRUE)

mrplot<-function(hideMenu=FALSE,folders=NULL,dimensions=NULL,intensity1D=NULL,
  zoom=NULL){
  mrbin.env$mrbinTMP$currentSpectrumOriginal<-NULL
  if(!is.null(folders)){
    if(is.null(dimensions)) dimensions<-rep("1D",length(folders))
    addToPlot(folder=folders[1],dimension=dimensions[1],add=FALSE)
    if(length(folders)>1){
      for(iTMP in 2:length(folders)){
         addToPlot(folder=folders[iTMP],dimension=dimensions[iTMP],add=TRUE)
      }
    }
  }
  if(!is.null(zoom)) zoom(zoom[1],zoom[2],zoom[3],zoom[4],refreshPlot=FALSE)
  if(!is.null(intensity1D)) intMin(value=intensity1D,refreshPlot=FALSE)
  plotMultiNMR()
  if(!hideMenu){
    stopTMP<-FALSE
    menuItemAdd<-"Add or remove spectra..."
    menuItemZoomInX<-"Zoom in (X)"
    menuItemZoomOutX<-"Zoom out (X)"
    menuItemZoomInY<-"Zoom in (Y)"
    menuItemZoomOutY<-"Zoom out (Y)"
    menuItemSetZoom<-"Set zoom..."
    currentItem<-menuItemAdd
    while(!stopTMP){
      menuList<-c(menuItemAdd,menuItemZoomInX,menuItemZoomOutX,menuItemZoomInY,
          menuItemZoomOutY,
          "Move left",
          "Move right",
          "Move up (1D)",
          "Move down (1D)",
          "Increase intensity (1D)",
          "Decrease intensity (1D)",
          "Move up (2D)",
          "Move down (2D)",
          "Increase contour level (2D)",
          "Decrease contour level (2D)",
          "Increase intensity (2D)",
          "Decrease intensity (2D)",
          menuItemSetZoom,
          "Set offset and scale of selected spectra..."#,
          #"Change colors..."
          )
      #if not 2D: remove up and down
      TMP1D<-FALSE
      TMP2D<-FALSE
      if(mrbin.env$mrbin$parameters$dimension=="1D"&
           !is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
        TMP1D<-TRUE
      }
      if(mrbin.env$mrbin$parameters$dimension=="2D"&
           !is.null(mrbin.env$mrbinTMP$currentSpectrumOriginal)){
        TMP2D<-TRUE
      }
      if(!is.null(mrbin.env$mrbinTMP$additionalPlots1D)) TMP1D<-TRUE
      if(!is.null(mrbin.env$mrbinTMP$additionalPlots2D)) TMP2D<-TRUE
      if(!TMP2D){
        menuList<-menuList[!menuList%in%c(menuItemZoomInY,menuItemZoomOutY,
                            "Move up (2D)","Move down (2D)",
                            "Increase contour level (2D)",
                            "Decrease contour level (2D)",
                            "Increase intensity (2D)",
                            "Decrease intensity (2D)"
                            )]
      }
      if(!TMP1D){
        menuList<-menuList[!menuList%in%c("Move up (1D)",
          "Move down (1D)")]
      }
      menuItem<-utils::select.list(menuList,
                   preselect=currentItem,
                   title ="Bin size [ppm]: ",graphics=TRUE)
      if(length(menuItem)==0|menuItem==""){
        stopTMP<-TRUE
      } else {
        currentItem<-menuItem
      }
      if(menuItem=="Move left") left(refreshPlot=FALSE)
      if(menuItem=="Move right") right(refreshPlot=FALSE)
      if(menuItem=="Move up (2D)") up(refreshPlot=FALSE)
      if(menuItem=="Move down (2D)") down(refreshPlot=FALSE)
      if(menuItem=="Move up (1D)"){
        #mrbin.env$mrbinplot$intensityOffset<-
        setOffset(mrbin.env$mrbinplot$intensityOffset-
            sort(mrbin.env$mrbinTMP$additionalPlots1D[[1]])[ceiling(length(mrbin.env$mrbinTMP$additionalPlots1D[[1]])*.99)]/5)*mrbin.env$mrbinplot$intensityScale
      }
      if(menuItem=="Move down (1D)"){
        #mrbin.env$mrbinplot$intensityOffset<-
        setOffset(mrbin.env$mrbinplot$intensityOffset+
             sort(mrbin.env$mrbinTMP$additionalPlots1D[[1]])[ceiling(length(mrbin.env$mrbinTMP$additionalPlots1D[[1]])*.99)]/5)*mrbin.env$mrbinplot$intensityScale
      }

      if(menuItem==menuItemZoomOutX) zoomOut(refreshPlot=FALSE,y=FALSE)
      if(menuItem==menuItemZoomOutY) zoomOut(refreshPlot=FALSE,x=FALSE)
      if(menuItem==menuItemZoomInX) zoomIn(refreshPlot=FALSE,y=FALSE)
      if(menuItem==menuItemZoomInY) zoomIn(refreshPlot=FALSE,x=FALSE)
      if(menuItem==menuItemSetZoom) zoom(refreshPlot=FALSE)
      if(menuItem=="Increase contour level (2D)") contPlus(refreshPlot=FALSE)
      if(menuItem=="Decrease contour level (2D)") contMin(refreshPlot=FALSE)
      if(menuItem=="Decrease intensity (1D)") intMin(refreshPlot=FALSE)
      if(menuItem=="Increase intensity (1D)") intPlus(refreshPlot=FALSE)
      if(menuItem=="Decrease intensity (2D)") intMin(refreshPlot=FALSE,dimension="2D")
      if(menuItem=="Increase intensity (2D)") intPlus(refreshPlot=FALSE,dimension="2D")
      if(menuItem=="Set offset and scale of selected spectra..."){
          if(TMP1D&TMP2D){
            menuItem3<-utils::select.list(c("Scale or offset 1D spectra",
                       "Scale or offset 2D spectra"),
                       preselect="Scale or offset 1D spectra",
                       title ="Please select",graphics=TRUE)
          } else {
            if(TMP1D) menuItem3<-"Scale or offset 1D spectra"
            if(TMP2D) menuItem3<-"Scale or offset 2D spectra"
          }
          if(menuItem3=="Scale or offset 1D spectra"&!
            is.null(mrbin.env$mrbinTMP$additionalPlots1D)){
            menuItem4<-utils::select.list(
              mrbin.env$mrbinTMP$additionalPlots1DMetadata[,5],
                     #preselect="Remove 1D spectra",
                     title ="Select spectrum",graphics=TRUE)
            if(!(length(menuItem4)==0|menuItem4=="")){
              scaleTMP<-readline(prompt=paste("New scaling factor, press enter to keep ",
                                      mrbin.env$mrbinTMP$additionalPlots1DMetadata[
                                        which(mrbin.env$mrbinTMP$additionalPlots1DMetadata[,5]==menuItem4),7]#scaling
                                      ,": ",sep=""))
              if(!scaleTMP==""){
                  mrbin.env$mrbinTMP$additionalPlots1DMetadata[
                                        which(mrbin.env$mrbinTMP$additionalPlots1DMetadata[,5]==menuItem4),7]<-scaleTMP
              }
              offsetTMP<-readline(prompt=paste("New offset factor, press enter to keep ",
                                   mrbin.env$mrbinTMP$additionalPlots1DMetadata[
                                        which(mrbin.env$mrbinTMP$additionalPlots1DMetadata[,5]==menuItem4),9]#offset
                                      ,": ",sep=""))
              if(!offsetTMP==""){
                  mrbin.env$mrbinTMP$additionalPlots1DMetadata[
                                        which(mrbin.env$mrbinTMP$additionalPlots1DMetadata[,5]==menuItem4),9]<-offsetTMP
              }
            }
          }
          if(menuItem3=="Scale or offset 2D spectra"&!
            is.null(mrbin.env$mrbinTMP$additionalPlots2D)){
            menuItem4<-utils::select.list(
              mrbin.env$mrbinTMP$additionalPlots2DMetadata[,5],
                     title ="Select spectrum",graphics=TRUE)
            if(!(length(menuItem4)==0|menuItem4=="")){
              scaleTMP<-readline(prompt=paste("New scaling factor, press enter to keep ",
                                      mrbin.env$mrbinTMP$additionalPlots2DMetadata[
                                        which(mrbin.env$mrbinTMP$additionalPlots2DMetadata[,5]==menuItem4),7]#scaling
                                      ,": ",sep=""))
              if(!scaleTMP==""){
                  mrbin.env$mrbinTMP$additionalPlots2DMetadata[
                                        which(mrbin.env$mrbinTMP$additionalPlots2DMetadata[,5]==menuItem4),7]<-scaleTMP
              }
              offsetTMP<-readline(prompt=paste("New offset factor, press enter to keep ",
                                   mrbin.env$mrbinTMP$additionalPlots2DMetadata[
                                        which(mrbin.env$mrbinTMP$additionalPlots2DMetadata[,5]==menuItem4),9]#offset
                                      ,": ",sep=""))
              if(!offsetTMP==""){
                  mrbin.env$mrbinTMP$additionalPlots2DMetadata[
                                        which(mrbin.env$mrbinTMP$additionalPlots2DMetadata[,5]==menuItem4),9]<-offsetTMP
              }
            }
          }

      }
      if(menuItem==menuItemAdd){
        menuItem2<-utils::select.list(c("Add spectra",
          "Remove spectra"#,"Hide or unhide spectra"
          ),
                     preselect="Add spectra",
                     title ="Please select",graphics=TRUE)
        if(menuItem2=="Add spectra"){
          #Ask for dimension
          menuItem3<-utils::select.list(c("Add 1D spectra",
          "Add 2D spectra"),
                     preselect="Add 1D spectra",
                     title ="Please select",graphics=TRUE)
          if(menuItem3=="Add 1D spectra") mrbin.env$mrbin$parameters$dimension<-"1D"
          if(menuItem3=="Add 2D spectra") mrbin.env$mrbin$parameters$dimension<-"2D"
          if(!(length(menuItem3)==0|menuItem3=="")){
            foldersTMP<-selectFolders()
            if(!foldersTMP=="stop"&length(mrbin.env$mrbin$parameters$NMRfolders)>0){
              for(iFolders in 1:length(mrbin.env$mrbin$parameters$NMRfolders)){
                addToPlot(folder=mrbin.env$mrbin$parameters$NMRfolders[iFolders]
                  #,dimension=dimension
                  )
              }
            }
          }
        }
        if(menuItem2=="Remove spectra"){
          if(TMP1D&TMP2D){
          menuItem3<-utils::select.list(c("Remove 1D spectra",
                     "Remove 2D spectra"),
                     preselect="Remove 1D spectra",
                     title ="Please select",graphics=TRUE)
          } else {
            if(TMP1D) menuItem3<-"Remove 1D spectra"
            if(TMP2D) menuItem3<-"Remove 2D spectra"
          }
          if(menuItem3=="Remove 1D spectra"&!
            is.null(mrbin.env$mrbinTMP$additionalPlots1D)){
            menuItem4<-utils::select.list(
              mrbin.env$mrbinTMP$additionalPlots1DMetadata[,5],
                     #preselect="Remove 1D spectra",
                     title ="Remove spectrum",graphics=TRUE)
            if(!(length(menuItem4)==0|menuItem4=="")){
              removeFromPlot(folder=menuItem4,
                dimension="1D")
            }
          }
          if(menuItem3=="Remove 2D spectra"&!
            is.null(mrbin.env$mrbinTMP$additionalPlots2D)){
            menuItem4<-utils::select.list(
              mrbin.env$mrbinTMP$additionalPlots2DMetadata[,5],
                     #preselect="Remove 1D spectra",
                     title ="Remove spectrum",graphics=TRUE)
            if(!(length(menuItem4)==0|menuItem4=="")){
              removeFromPlot(folder=menuItem4,
                dimension="2D")
            }
          }
        }
        #if(menuItem2=="Hide or unhide spectra"){
        #}
      }
      #if(menuItem=="Remove spectra...")

      plotMultiNMR()
      #utils::flush.console()
    }
  }

}
