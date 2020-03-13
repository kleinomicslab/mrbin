# mrbin News

## Version 1.3.0

* Noise is now estimated from raw data points (not bin values) and mean number of data points per bin
* New noise calculations tend to be lower than the previous estimates. SNR levels may need to be increased by a factor of 2.5. Default SNR levels have been raised to reflect this.
* Solvent region change: Any bins that contain any part of these regions are removed (previously only if mean of bin was within region)
* Added new way of creating sample names from EXPNO and folder name
* Some minor and major fixes to increase user friendliness
* Some minor and major bug fixes to increase stability
* Speed improvements

## Version 1.2.0

* Default 1D noise range was changed to 10 - 9.5 ppm (was 10-9.4). Default 1D binning area changed to 10-0.5 (was 9.5-0.5)
* Some minor and major fixes to increase user friendliness
* Some minor and major bug fixes to increase stability


## Version 1.1.0

* All parameters can now be conveniently set from the command line by passing a list of parameters to mrbin()
* Special bin regions can now be defined instead of rectangular bins, e.g. for lipid analysis
* mrbin() now returns the group factors as an additional element "factors" in the returned list
* New function setCurrentSpectrum() to set current spectrum, e.g. for plotting 
* Reloading a parameter set from a file now checks whether all required parameters have been set and adds missing ones
* Results plots are nicer now. A new plot shows which regions were left after all processing steps.
* Some minor and major fixes to increase user friendliness
* Some minor and major bug fixes to increase stability
