# mrbin - Metabolomics Data Analysis Functions

This package is a collection of functions for processing and analyzing
metabolomics data.

The namesake function mrbin() uses spectral binning to convert 1D or 2D Nuclear
Magnetic Resonance (NMR) data into a matrix of values
suitable for further data analysis and performs basic processing steps in a
reproducible way. Negative values, a common issue in NMR data, are replaced by
positive values. All used parameters are stored in a readable text file and can
be restored from that file to enable exact reproduction of the data at a later
time.

The atnv algorithm for replacing negative values in NMR data sets can be
employed using atnv().

NMR plotting functions are found in mrplot().

Artificial Neural Network features can be analyzed using Feature Impact
Assessment (FIA) using the function fia().


## Installation

To install  mrbin, please install the latest version of R first. Then install 
the package as follows:

To install the latest stable version from CRAN:

```
install.packages("mrbin")
library("mrbin")
```

You can find more details and examples in the vignette file:

```
vignette("mrbin")
```


## mrbin: Magnetic Resonance Binning, Integration and Normalization

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

mrbin() returns an (invisible) mrbin object containing the following variables: 
* bins: A matrix containing bin data for all samples, Depending on the option you chose, the data will be cleaned up and scaled.
* parameters: A list containing all parameters used to create the bin matrix.
* metadata: A list containing metadata, if provided.
* transformations: A character vector containing information on the data transformation and scaling that has been performed, for example reference scaling, PQN, atnv, log transform, etc.
* changeLog: A data.frame containing information on documented changes that were made to the data, including time stamps.
* changeValues: A list containing control values, enabling verifying changes by checkmrbin(mrbinResults)

Several files may be written to the chosen directory:
* A .txt file containing all parameters and potential warning messages from the mrbin run. This file can be reloaded to R using recreatemrbin("filename"). This will enable reusing parameters used in a previous run and can help increase reproducibility.
* A .Rdata file containing the generated mrbin data object.
* A .pdf file containing quality control plots of the raw binned data
* A .pdf file containing quality control plots of the processed binned data
* A .pdf file containing preview plots of the chosen signal-to-noise levels of selected spectra

### Submitting Parameters at the Command Line
Parameters can be submitted at the command line, using the following syntax:

```
mrbin(parameters=list(dimension="1D",
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

## atnv: Affine Transformation of Negative Values

Please find details on the atnv algorithm in vignette("mrbin").


## fia: Feature Impact Assessment

Please find details on the fia algorithm in vignette("mrbin").


## Known Issues


### Firewall Warnings
If parallel computing is turned on and the package parallel is installed,
mrbin will try to use the socket approach for computing. This requires
establishing network connections to the local cluster, which might
trigger the firewall. It is safe to unblock these connections.

### Pop-Up Windows
mrbin is set up to ask for user input through pop-up windows. This requires
graphics support, otherwise the user input will be asked through command line
menus, which is less user friendly but still offers the full functionality.

### Apple/Mac Computers And RStudio
In some cases, running mrbin from within RStudio on Apple computers will not
generate pop-up windows. To enable pop-up windows, it might be helpful to install
the newest version of xquartz from https://www.xquartz.org.

### Spectra are Missing
If a Bruker spectrum is not shown during browsing, please make sure a file
with filename title is present in the PROCNO folder of that spectrum. You
can create a title file by opening the spectrum in Bruker Topspin, selecting
the Title tab, entering a title and clicking the disk symbol for saving.


## Built With

* [roxygen2]
* [devtools]


## Authors

* **Matthias S. Klein** - [KleinOmicsLab](https://github.com/kleinomicslab/)


## Citation
If you are using mrbin in a publication, please cite the following manuscript:
Klein, M.S. (2021): Affine Transformation of Negative Values for NMR Metabolomics 
Using the mrbin R Package. J. Proteome Res. 20(2):1397-1404, 
DOI: 10.1021/acs.jproteome.0c00684


## License

This project is licensed under GPL-3.0.
