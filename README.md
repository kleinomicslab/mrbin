# mrbin

Nuclear Magnetic Resonance (NMR) is widely used for metabolomics research. This package converts 1D or 2D NMR data into a matrix of values suitable for further data analysis and performs basic processing steps in a reproducible way. Negative values, a common issue in NMR data, are replaced by positive values. All used parameters are stored in a readable text file and can be restored from that file to enable exact reproduction of the data at a later time.

## Getting Started

The main functions of this package are controlled via the mrbin() function. Most other functions in the package will not usually be ever called by the user, but serve internal purposes. Results returned include the final bin list and a set of used parameters.

### Installation

To install  mrbin, please install the latest version of R first. Then install the package as follows:

```
install.packages("mrbin")
```

### Running

To use this package, you will need your NMR data in the Bruker file format accessible on your computer. Please make sure your data is Fourier transformed, phase corrected, baseline corrected, and correctly referenced. The data has to be stored in folders according to standard Bruker folders, that means foldername/1/pdata/1 etc. Experiment numbers and processing numbers can be freely chosen.

This package has been tested for 1D NOESY and 2D 1H-13C HSQC spectra.

Before starting mrbin, take a look at your NMR data, for example in Bruker Topspin, and decide on the following:
* Bin area: Area where signals are observed in your data set
* Bin width: Should match roughly the width of a singlet peak in your data set. Given in ppm.
* Bin height (only 2D): Should match roughly the height of a singlet peak in your data set. Given in ppm.
* Solvent area: Area to exclude to remove solvent artifacts
* Additional areas to be removed: Any other area containing artifacts, such as streaks surrounding strong peaks.

Once you have decided on these values, you can start mrbin using the following code:

```
library("mrbin")
mrbinResults<-mrbin()
```

This will start a series of questions that will guide you through the parameters to be used. 

The sequence of data processing is as follows:

* Gathering all parameters from user
* Reading Bruker NMR data
* Scaling to reference region
* Binning 
* Removing solvent region
* Removing additional regions
* Summing up region containing peaks with unstable positions such as citric acid
* Removal of bins containing mostly noise
* Cropping of HSQC spectra to the region along the diagonal
* PQN transformation
* Removal of negative values
* Log transform
* Plotting a PCA
* Saving bins and parameters to the hard drive

After finishing, a PCA can be displayed. mrbin() also returns a list containing two variables: 

* bins: A matrix containing bin data for all samples, Depending on the option you chose, the data will be cleaned up and scaled.
* parameters: A list containing all parameters used to create the bin matrix.

Two text files may be written to the working directory:
* A .txt file containing all parameters. This file can be reloaded to R using recreatemrbin("filename"). This will enable reusing parameters used in a previous run and can help increase reproducibility.
* A .csv file containing the bin data for use in other software tools.

### Recreating Data and Parameters
In order to create reproducible results, mrbin will save the used parameters to a text file. Please keep this file. You may want to share this file in a data repository when publishing your findings.

While it is fine to view the parameter text file in a text editor, please do not change its contents, as this may break its formatting.


In order to recreate a previous data set, or to reload previously used parameters, use:

```
mrbin()
```

and select ""Reload from file" when asked "Set parameters or use existing parameters?". This will restore all parameters that were previously used. If the file was created using an older version of mrbin, this may cause inconsistencies. Missing parameters will be added using standard parameters. Ideally, download the older mrbin version at kleinomicslab.com and use the old version to recreated the data in an exact way.

Please be aware that bins will have to be recalculated, so the original NMR spectra will have to be present to do this.


## Built With

* [roxygen2]
* [devtools]


## Authors

* **Matthias S. Klein** - [KleinOmicsLab](https://github.com/kleinomicslab)


## License

This project is licensed under GPL-3.0.
