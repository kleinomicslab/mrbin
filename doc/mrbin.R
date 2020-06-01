### R code from vignette source 'mrbin.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: mrbin.Rnw:65-66
###################################################
library(mrbin)


###################################################
### code chunk number 2: mrbin.Rnw:156-175
###################################################
mrbinResults<-mrbin(silent=TRUE,
     setDefault=FALSE,
     parameters=list(verbose=TRUE,
             dimension="1D",
             binMethod="Rectangular bins",
             binwidth1D=.01,
             referenceScaling="Yes",
             removeSolvent="Yes",
             removeAreas="No",
             sumBins="No",
             noiseRemoval="Yes",
             PQNScaling="Yes",
             fixNegatives="Yes",
             logTrafo="Yes",
             saveFiles="Yes",
             NMRfolders=c(system.file("extdata/1/10/pdata/10",package="mrbin"),
                          system.file("extdata/2/10/pdata/10",package="mrbin"),
                          system.file("extdata/3/10/pdata/10",package="mrbin"))
     ))


###################################################
### code chunk number 3: mrbin.Rnw:192-206
###################################################
mrbinResults<-mrbin(silent=TRUE,
     setDefault=FALSE,
     parameters=list(verbose=TRUE,
               dimension="2D",
               binwidth2D=0.1,
               binheight=3,
               PQNScaling="No",
               fixNegatives="No",
               logTrafo="No",
               signal_to_noise2D=20,
               NMRfolders=c(system.file("extdata/1/12/pdata/10",package="mrbin"),
                       system.file("extdata/2/12/pdata/10",package="mrbin"),
                       system.file("extdata/3/12/pdata/10",package="mrbin"))
               ))


###################################################
### code chunk number 4: mrbin.Rnw:240-266
###################################################
results <- mrbin(silent=TRUE,parameters=list(binMethod="Custom bin list",
 dimension="1D",specialBinList=matrix(c(
                               5.45,5.2,0,160,
                               2.9,2.74,0,160,
                               2.14,1.93,0,160,
                               1.41,1.2,0,160,
                               0.94,0.8,0,160,
                               2.44,2.2,0,160,
                               4.325,4.26,0,160
                               ),ncol=4,byrow=TRUE,dimnames=list(c(
                               "-CH=CH- Methene",
                               "=CH-CH2-CH= Diallylic",
                               "-CH2-CH=CH- Allylic",
                               "-CH2- Methylene",
                               "-CH3 Methyl",
                               "COO-CH2-CH2- Methylene_to_carboxyl",
                               "Glycerol"
                               ),NULL)),
 referenceScaling="Yes",reference1D=c(0.03,-0.03),removeSolvent="No",
 removeAreas="No",sumBins="No",trimZeros="Yes",noiseRemoval="No",
 PQNScaling="No",fixNegatives="Yes",logTrafo="No",defineGroups="No",PCA="Yes",
 createBins="Yes",useAsNames="Folder names",saveFiles="No",verbose=TRUE,
 NMRfolders=c(system.file("extdata/1/10/pdata/10",package="mrbin"),
              system.file("extdata/2/10/pdata/10",package="mrbin"),
              system.file("extdata/3/10/pdata/10",package="mrbin"))
 ))


