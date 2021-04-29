# mrbin - Magnetic Resonance Binning, Integration and Normalization

Nuclear Magnetic Resonance (NMR) is widely used for metabolomics research. This 
package uses spectral binning to convert 1D or 2D NMR data into a matrix of values 
suitable for further data analysis and performs basic processing steps in a 
reproducible way. Negative values, a common issue in NMR data, are replaced by 
positive values. All used parameters are stored in a readable text file and can 
be restored from that file to enable exact reproduction of the data at a later 
time.

## Citation
If you are using mrbin in a publication, please cite the following manuscript:
Klein, M.S. (2021): Affine Transformation of Negative Values for NMR Metabolomics 
Using the mrbin R Package. J. Proteome Res. 20(2):1397-1404, 
DOI: 10.1021/acs.jproteome.0c00684

## Getting Started

The main functions of this package are controlled via the mrbin() function. Most 
other functions in the package will not usually be ever called by the user, but 
serve internal purposes. Results returned include the final bin list and a set 
of used parameters.

### Installation

To install  mrbin, please install the latest version of R first. Then install 
the package as follows:

To install the latest stable version from CRAN:

```
install.packages("mrbin")
library("mrbin")
```

You can get more help in the vignette file:

```
vignette("mrbin")
```

#### Development Version

To install the latest development version of mrbin from Github, use this code:

```
library(devtools)
install_github("kleinomicslab/mrbin")
```

To be able to run devtools, you may need to install additional software.

### Running

To use this package, you will need your NMR data in the Bruker file format 
accessible on your computer. Please make sure your data is Fourier transformed, 
phase corrected, baseline corrected, and correctly referenced. The data has to 
be stored in folders according to standard Bruker folders, that means 
foldername/1/pdata/1 etc. Experiment numbers and processing numbers can be 
freely chosen.

This package has been tested for 1D NOESY and 2D 1H-13C HSQC spectra.

Before starting mrbin, take a look at your NMR data, for example in Bruker 
Topspin, and decide on the following:
* Bin area: Area where signals are observed in your data set
* Bin width: Should match roughly the width of a singlet peak in your data set. Given in ppm.
* Bin height (only 2D): Should match roughly the height of a singlet peak in your data set. Given in ppm.
* Solvent area: Area to exclude to remove solvent artifacts
* Additional areas to be removed: Any other area containing artifacts, such as streaks surrounding strong peaks.

mrbin will also show you preview plots for these parameters during the run.

You can start mrbin using the following code:

```
mrbinResults<-mrbin()
```

This will start a series of questions that will guide you through the parameters to be used. 

The sequence of data processing is as follows:

* Gathering all parameters from user
* Creating a set with coordinates of each bin 
* Removing solvent region
* Removing additional regions
* Cropping of HSQC spectra to the region along the diagonal
* Summing or merging regions containing peaks with unstable positions such as citric acid
* Reading Bruker NMR data
* Scaling to reference region
* Binning 
* Removal of bins containing mostly noise
* PQN transformation
* Replacement of negative values
* Log transform
* Plotting a quality control plot, including a PCA plot
* Saving bins, parameters and the plot to the hard drive

mrbin() also returns an (invisible) list containing three variables: 

* bins: A matrix containing bin data for all samples, Depending on the option you chose, the data will be cleaned up and scaled.
* parameters: A list containing all parameters used to create the bin matrix.
* factors: A vector containing group names for all samples.

Up to three files may be written to the chosen directory:
* A .txt file containing all parameters. This file can be reloaded to R using recreatemrbin("filename"). This will enable reusing parameters used in a previous run and can help increase reproducibility.
* A .csv file containing the bin data for use in other software tools.
* A .pdf file containing quality control plots, including a PCA plot 

### Recreating Data and Parameters
In order to create reproducible results, mrbin will save the used parameters to a text file. Please keep this file. You may want to share this file in a data repository when publishing your findings.

While it is fine to view the parameter text file in a text editor, please do not 
change its contents, as this may break its formatting.


In order to recreate a previous data set, or to reload previously used parameters, 
use:

```
mrbin()
```

and select ""Reload from file" when asked "Set parameters or use existing 
parameters?". This will restore all parameters that were previously used. If the 
file was created using an older version of mrbin, this may cause inconsistencies. 
Missing parameters will be added using standard parameters. Ideally, download the 
older mrbin version at kleinomicslab.com and use the old version to recreated the 
data in an exact way.

Please be aware that bins will have to be recalculated, so the original NMR 
spectra will have to be present to do this.

### Submitting Parameters at the Command Line
Parameters can be submitted at the command line, using the following syntax:

```
mrbin(silent=TRUE,
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
             NMRfolders=c("C:/NMR/Sample01/1/pdata/1",
                          "C:/NMR/Sample19/1/pdata/1",
                          "C:/NMR/Sample61/1/pdata/1")
     ))
```

This will set up all parameters and run all steps without asking for user input.

When setting silent=FALSE, the user will be guided through the user input 
questionnaire to make adjustments to the parameters.

### Affine Transformation of Negative Values

The function atnv replaces (column-wise) negative values by a small positive
number. The number is calculated as an affine transformation to the range of
the lowest positive number to 0,01*the lowest positive number (of this
column). Ranks stay unchanged. Positive numbers are not altered.
If sample-wise noise levels are available, the median noise level of samples
with negative values is calculated and replaces the lowest positive number in
case it is smaller. If no noise data is available, the 1% percentile of all
positive values in the data set is used as an estimate.
It is recommended to use this function AFTER noise removal and other data
clean-up methods, as it may alter (reduce) the noise level of the binned data.
If no NMR data and noise levels are provided as arguments, the function will
use NMR data and noise levels from the global variables mrbin.env$bins and
mrbin.env$mrbinTMP.

To use own data:
```
atnv(NMRdataMatrix,noiseLevelVector)
```
To use current mrbin data, use the following syntax. This requires data loaded 
using mrbin(). This is usually not necessary as it is included in the mrbin work 
flow.
```
atnv()
```

## Known Issues

### Pop-Up Windows
mrbin is set up to ask for user input through pop-up windows. This requires 
graphics support, otherwise the use input will be asked through command line 
menus, which is less user friendly but still offers the full functionality.

### Apple/Mac Computers And RStudio
In some cases, running mrbin from within RStudio on Apple computers will not 
generate pop-up windows. To enable pop-up windows, it might be helpful to install 
the newest version of xquartz from https://www.xquartz.org.

## Built With

* [roxygen2]
* [devtools]


## Authors

* **Matthias S. Klein** - [KleinOmicsLab](https://github.com/kleinomicslab/)


## License

This project is licensed under GPL-3.0.
