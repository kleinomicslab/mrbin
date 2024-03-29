% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrbin.R
\name{fia}
\alias{fia}
\title{A function identifying features of importance.}
\usage{
fia(
  model,
  dataSet,
  factors,
  nSeed = 6,
  numberOfSamples = 100,
  maxFeatures = 10000,
  innerLoop = 300,
  verbose = TRUE,
  maxNumberAllTests = 5,
  firstLevel = 1,
  saveMemory = FALSE,
  kerasClearMemory = 0,
  functionNamePredict = "predict",
  parameterNameObject = "object",
  parameterNameData = "x",
  ...
)
}
\arguments{
\item{model}{A predictive model. Make sure to have loaded all required packages before starting this function}

\item{dataSet}{An object containing data, columns=features, rows=samples. This should be either a matrix or a dataframe, depending on which of these two the specific prediction function requires}

\item{factors}{A factor vector with group membership of each sample in the data set. Order of levels must correspond to the number predicted by the model}

\item{nSeed}{Number of times that the test will be repeated, selecting different random features}

\item{numberOfSamples}{Number of samples that will be randomly chosen from each group}

\item{maxFeatures}{Maximum number of features that will be tested. Larger numbers will be split into child nodes without testing to increase speed}

\item{innerLoop}{Number of repeated loops to test additional child nodes}

\item{verbose}{A logical vector to turn messages on or off}

\item{maxNumberAllTests}{Combinations of features of this length or shorter will not be split in half to create two children, but into multiple children with one feature left out each. This is done make sure no combination is missed.}

\item{firstLevel}{Numeric value of first level or group. Usually 1 but for glm such as in the example this needs to be 0.}

\item{saveMemory}{Save memory by performing only two predictions per step, which will be much slower. If you are using keras, use parameter kerasClearMemory=2 instead to free more memory and be a lot faster. FALSE to turn off.}

\item{kerasClearMemory}{Save memory by clearing model from memory and reloading the model between chunks of predictions. Will only work when using package keras. 0=off, 1=medium (reload between repeat with different seeds), 2=maximum memory savings (reload after each run for a single sample). This will write a model file to the working directory.}

\item{functionNamePredict}{The name of the prediction function. This only needs to be changed if the prediction function is not called predict}

\item{parameterNameObject}{The name of the parameter for passing the model to the prediction function}

\item{parameterNameData}{The name of the parameter for passing the data to the prediction function}

\item{...}{Optional, additional parameters that will be passed to the prediction function.}
}
\value{
A list of results: scoresSummary A vector of fia scores for the whole dataset; scores contains vectors of fia scores for each predicted group; scoresIndividual A list of fia scores for each individual sample; fiaListPerSample A list of important combinations of features for each predicted sample; fiaMatrix A list of fia scores for each predicted group.
}
\description{
This function finds features that can change the outcomes of a model's prediction.
Example: fia=1.00 means single compound found in all but 0 percent of samples.
fia=2.45 indicates this compound is found in pairs in all but 45 percent of tested samples
A function named predict needs to be present for this to work. If the function name
of the prediction function is different, the function name has to be provided in
the parameter functionNamePredict.
}
\examples{
 #First, define group membership and create the example feature data
 group<-factor(c(rep("Group1",4),rep("Group2",5)))
 names(group)<-paste("Sample",1:9,sep="")
 dataset<-data.frame(
   Feature1=c(5.1,5.0,6.0,2.9,4.8,4.6,4.9,3.8,5.1),
   Feature2=c(2.6,4.0,3.2,1.2,3.1,2.1,4.5,6.1,1.3),
   Feature3=c(3.1,6.1,5.8,5.1,3.8,6.1,3.4,4.0,4.4),
   Feature4=c(5.3,5.2,3.1,2.7,3.2,2.8,5.9,5.8,3.1),
   Feature5=c(3.2,4.4,4.8,4.9,6.0,3.6,6.1,3.9,3.5),
   Feature6=c(6.8,6.7,7.2,7.0,7.3,7.1,7.2,6.9,6.8)
   )
 rownames(dataset)<-names(group)
 #train a model - here we use a logit model instead of ANN as a demonstration
 mod<-glm(group~Feature1+Feature2+Feature3+Feature4+Feature5+Feature6,
   data=data.frame(group=group,dataset),family="binomial")
 fiaresults<-fia(model=mod,dataSet=dataset,factors=group,parameterNameData="newdata",
   firstLevel=0,type="response")
 fiaresults$scores
}
