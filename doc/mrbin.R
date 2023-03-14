### R code from vignette source 'mrbin.Rnw'

###################################################
### code chunk number 1: mrbin.Rnw:53-54 (eval = FALSE)
###################################################
## install.packages("mrbin")


###################################################
### code chunk number 2: mrbin.Rnw:59-60
###################################################
library(mrbin)


###################################################
### code chunk number 3: mrbin.Rnw:100-101 (eval = FALSE)
###################################################
## mrbinResults<-mrbin()


###################################################
### code chunk number 4: mrbin.Rnw:167-182 (eval = FALSE)
###################################################
## mrbinResults<-mrbin(parameters=list(dimension="1D",
##              binwidth1D=.02,
##              referenceScaling="Yes",
##              removeSolvent="Yes",
##              solventRegion=c(4.95,4.65),
##              noiseRemoval="Yes",
##              signal_to_noise1D=35,
##              noiseThreshold=0.75,
##              PQNScaling="No",
##              fixNegatives="No",
##              logTrafo="No",
##              NMRfolders=c(system.file("extdata/1/10/pdata/10",package="mrbin"),
##                           system.file("extdata/2/10/pdata/10",package="mrbin"),
##                           system.file("extdata/3/10/pdata/10",package="mrbin"))
##      ))


###################################################
### code chunk number 5: mrbin.Rnw:185-203
###################################################
mrbinResults<-mrbin(silent=TRUE,#this suppresses all interactive prompts
     setDefault=FALSE,
     parameters=list(verbose=TRUE,
             dimension="1D",
             binwidth1D=.02,
             referenceScaling="Yes",
             removeSolvent="Yes",
             solventRegion=c(4.95,4.65),
             noiseRemoval="Yes",
             signal_to_noise1D=35,
             noiseThreshold=0.75,
             PQNScaling="No",
             fixNegatives="No",
             logTrafo="No",
             NMRfolders=c(system.file("extdata/1/10/pdata/10",package="mrbin"),
                          system.file("extdata/2/10/pdata/10",package="mrbin"),
                          system.file("extdata/3/10/pdata/10",package="mrbin"))
     ))


###################################################
### code chunk number 6: mrbin.Rnw:211-212
###################################################
setNoiseLevels(mrbinResults,plotOnly=TRUE)


###################################################
### code chunk number 7: mrbin.Rnw:222-227 (eval = FALSE)
###################################################
## mrbinResults<-mrbin(parameters=list(
##    NMRfolders=c("C:/Bruker/TopSpin3.6.1/data/guest/nmr/sample_1/10/pdata/10",
##              "C:/Bruker/TopSpin3.6.1/data/guest/nmr/sample_2/10/pdata/10",
##              "C:/Bruker/TopSpin3.6.1/data/guest/nmr/sample_3/10/pdata/10")
##    ))


###################################################
### code chunk number 8: mrbin.Rnw:246-247 (eval = FALSE)
###################################################
## mrbinResults<-setNoiseLevels(mrbinResults)


###################################################
### code chunk number 9: mrbin.Rnw:259-260 (eval = FALSE)
###################################################
## mrbinResults<-removeNoise(mrbinResults)


###################################################
### code chunk number 10: mrbin.Rnw:265-266 (eval = FALSE)
###################################################
## mrbinResults<-atnv(mrbinResults)


###################################################
### code chunk number 11: mrbin.Rnw:273-275 (eval = FALSE)
###################################################
## mrbinResults<-setDilutionFactors(mrbinResults)
## mrbinResults<-dilutionCorrection(mrbinResults)


###################################################
### code chunk number 12: mrbin.Rnw:281-282 (eval = FALSE)
###################################################
## mrbinResults<-PQNScaling(mrbinResults)


###################################################
### code chunk number 13: mrbin.Rnw:289-290 (eval = FALSE)
###################################################
## mrbinResults<-logTrafo(mrbinResults)


###################################################
### code chunk number 14: mrbin.Rnw:297-298 (eval = FALSE)
###################################################
## mrbinResults<-unitVarianceScaling(mrbinResults)


###################################################
### code chunk number 15: mrbin.Rnw:303-304 (eval = FALSE)
###################################################
## plotResults(mrbinResults)


###################################################
### code chunk number 16: mrbin.Rnw:312-313 (eval = FALSE)
###################################################
## mrbinResults<-removeSpectrum(mrbinResults)


###################################################
### code chunk number 17: mrbin.Rnw:321-322 (eval = FALSE)
###################################################
## mrbinResults<-metadatamrbin(mrbinResults)


###################################################
### code chunk number 18: mrbin.Rnw:331-334
###################################################
mrbinResults<-metadatamrbin(mrbinResults,metadata=list(
  projectTitle="Test project",
  factors=factor(c("Control","Control","Treatment"))))


###################################################
### code chunk number 19: mrbin.Rnw:340-341
###################################################
plotPCA(mrbinResults)


###################################################
### code chunk number 20: mrbin.Rnw:351-360
###################################################
mrbinResults<-editmetabolitesmrbin(mrbinResults,borders=matrix(c(
      1.346,1.324,21,23,
      3.052,3.043,30.5,33.5,
      4.066,4.059,57,59.5
    ),ncol=4,byrow=TRUE),metabolitenames=c(
    "Lactate",
    "Creatinine",
    "Creatinine"
    ))


###################################################
### code chunk number 21: mrbin.Rnw:365-366
###################################################
plotPCA(mrbinResults,loadings=TRUE,annotate=TRUE)


###################################################
### code chunk number 22: mrbin.Rnw:377-385 (eval = FALSE)
###################################################
## mrbinResults<-editmrbin(
##   mrbinResults,
##   bins=mrbinResults$bins,#omit this line if no changes are made to bins
##   parameters=list(noiseThreshold=0.75),
##   metadata=list(projectTitle="Test project"),
##   comment="Changed title and noise parameters",
##   transformations="Scaling"#If you change bins, provide a short explanation
##   )


###################################################
### code chunk number 23: mrbin.Rnw:390-391 (eval = FALSE)
###################################################
## mrbinResults$changeLog


###################################################
### code chunk number 24: mrbin.Rnw:397-398
###################################################
checkmrbin(mrbinResults)


###################################################
### code chunk number 25: mrbin.Rnw:406-421 (eval = FALSE)
###################################################
## mrbinResults2D<-mrbin(setDefault=TRUE,parameters=list(dimension="2D",
##                binwidth2D=0.3,
##                binheight=4,
##                removeSolvent="Yes",
##                solventRegion=c(4.95,4.65),
##                noiseRemoval="Yes",
##                signal_to_noise2D=5,
##                noiseThreshold=0.75,
##                PQNScaling="No",
##                fixNegatives="No",
##                logTrafo="No",
##                NMRfolders=c(system.file("extdata/1/12/pdata/10",package="mrbin"),
##                        system.file("extdata/2/12/pdata/10",package="mrbin"),
##                        system.file("extdata/3/12/pdata/10",package="mrbin"))
##                ))


###################################################
### code chunk number 26: mrbin.Rnw:423-441
###################################################
mrbinResults2D<-mrbin(silent=TRUE,#this suppresses all interactive prompts
     setDefault=TRUE,
     parameters=list(verbose=TRUE,
               dimension="2D",
               binwidth2D=0.3,
               binheight=4,
               removeSolvent="Yes",
               solventRegion=c(4.95,4.65),
               noiseRemoval="Yes",
               signal_to_noise2D=5,
               noiseThreshold=0.75,
               PQNScaling="No",
               fixNegatives="No",
               logTrafo="No",
               NMRfolders=c(system.file("extdata/1/12/pdata/10",package="mrbin"),
                       system.file("extdata/2/12/pdata/10",package="mrbin"),
                       system.file("extdata/3/12/pdata/10",package="mrbin"))
               ))


###################################################
### code chunk number 27: mrbin.Rnw:448-449
###################################################
setNoiseLevels(mrbinResults2D,plotOnly=TRUE)


###################################################
### code chunk number 28: mrbin.Rnw:461-462 (eval = FALSE)
###################################################
## load("C:/Users/User/Documents/mrbin.Rdata")


###################################################
### code chunk number 29: mrbin.Rnw:470-473
###################################################
mrbinResults2D<-metadatamrbin(mrbinResults,metadata=list(
  projectTitle="Test project",
  factors=factor(c("Control","Control","Treatment"))))


###################################################
### code chunk number 30: mrbin.Rnw:479-480
###################################################
plotPCA(mrbinResults2D)


###################################################
### code chunk number 31: mrbin.Rnw:486-495
###################################################
mrbinResults2D<-editmetabolitesmrbin(mrbinResults2D,borders=matrix(c(
      1.346,1.324,21,23,
      3.052,3.043,30.5,33.5,
      4.066,4.059,57,59.5
    ),ncol=4,byrow=TRUE),metabolitenames=c(
    "Lactate",
    "Creatinine",
    "Creatinine"
    ))


###################################################
### code chunk number 32: mrbin.Rnw:500-501
###################################################
plotPCA(mrbinResults2D,loadings=TRUE,annotate=TRUE)


###################################################
### code chunk number 33: mrbin.Rnw:510-533 (eval = FALSE)
###################################################
## results <- mrbin(parameters=list(dimension="1D",binMethod="Custom bin list",
##  specialBinList=matrix(c(5.45,5.2,0,160,
##                          2.9,2.74,0,160,
##                          2.14,1.93,0,160,
##                          1.41,1.2,0,160,
##                          0.94,0.8,0,160,
##                          2.44,2.2,0,160,
##                          4.325,4.26,0,160
##                          ),ncol=4,byrow=TRUE,dimnames=list(c(
##                          "-CH=CH- Methene",
##                          "=CH-CH2-CH= Diallylic",
##                          "-CH2-CH=CH- Allylic",
##                          "-CH2- Methylene",
##                          "-CH3 Methyl",
##                          "COO-CH2-CH2- Methylene_to_carboxyl",
##                          "Glycerol"
##                          ),NULL)),
##  referenceScaling="Yes",reference1D=c(0.03,-0.03),removeSolvent="No",
##  noiseRemoval="No",PQNScaling="No",fixNegatives="Yes",logTrafo="No",
##  NMRfolders=c(system.file("extdata/1/10/pdata/10",package="mrbin"),
##               system.file("extdata/2/10/pdata/10",package="mrbin"),
##               system.file("extdata/3/10/pdata/10",package="mrbin"))
##  ))


###################################################
### code chunk number 34: mrbin.Rnw:553-554 (eval = FALSE)
###################################################
## mrbinResults<-mrbin()


###################################################
### code chunk number 35: mrbin.Rnw:613-614 (eval = FALSE)
###################################################
## atnv(NMRdataMatrix,noiseLevelVector)


###################################################
### code chunk number 36: mrbin.Rnw:629-630 (eval = FALSE)
###################################################
## mrplot()


###################################################
### code chunk number 37: mrbin.Rnw:639-642 (eval = FALSE)
###################################################
## readBruker(dimension="1D",
##   folder=system.file("extdata/1/10/pdata/10",package="mrbin"))
## plotNMR()


###################################################
### code chunk number 38: mrbin.Rnw:648-651 (eval = FALSE)
###################################################
## addToPlot(dimension="1D",
##   folder="C:/Bruker/TopSpin3.6.1/data/guest/nmr/sample_1/12/pdata/10")
## plotNMR()


###################################################
### code chunk number 39: mrbin.Rnw:657-660 (eval = FALSE)
###################################################
## addToPlot(dimension="2D",
##   folder=system.file("extdata/1/12/pdata/10",package="mrbin"))
## plotNMR()


###################################################
### code chunk number 40: mrbin.Rnw:665-672 (eval = FALSE)
###################################################
## zoom(left=4.6, right=2, top=10, bottom=150) #Exact zoom
## zoomIn() #Zoom in
## zoomOut() #Zoom out
## intPlus() #Increase intensity
## intMin() #Decrease intensity
## left() #Move spectrum to the left
## right() #Move spectrum to the right


###################################################
### code chunk number 41: mrbin.Rnw:677-681 (eval = FALSE)
###################################################
## contMin() #Decrease minimum contour level (show more small peaks)
## contPlus() #Increase minimum contour level (remove small peaks)
## up() #Move spectrum up
## down() #Move spectrum down


###################################################
### code chunk number 42: mrbin.Rnw:701-719 (eval = FALSE)
###################################################
## #First, define group membership and create the example feature data
## group<-factor(c(rep("Group1",4),rep("Group2",5)))
## names(group)<-paste("Sample",1:9,sep="")
## dataset<-data.frame(
##     Feature1=c(5.1,5.0,6.0,2.9,4.8,4.6,4.9,3.8,5.1),
##     Feature2=c(2.6,4.0,3.2,1.2,3.1,2.1,4.5,6.1,1.3),
##     Feature3=c(3.1,6.1,5.8,5.1,3.8,6.1,3.4,4.0,4.4),
##     Feature4=c(5.3,5.2,3.1,2.7,3.2,2.8,5.9,5.8,3.1),
##     Feature5=c(3.2,4.4,4.8,4.9,6.0,3.6,6.1,3.9,3.5),
##     Feature6=c(6.8,6.7,7.2,7.0,7.3,7.1,7.2,6.9,6.8)
## )
## rownames(dataset)<-names(group)
## #train the logit model
## mod<-glm(group~Feature1+Feature2+Feature3+Feature4+Feature5+Feature6,
##    data=data.frame(group=group,dataset),family="binomial")
## fiaresults<-fia(model=mod,dataSet=dataset,factors=group,
##   parameterNameData="newdata",firstLevel=0,type="response")
## fiaresults$scores


