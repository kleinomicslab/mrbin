# mrbin News

## Version 1.5.2

* Speed improvements for 1D spectra
* New quality check function checkBaseline: For severe baseline distortions in noise region, a warning is displayed
* New quality check in function binMultiNMR2: For severe baseline distortions in the reference region, the reference peak integral might be negative. In this case,a warning is displayed and the absolute value is used instead
* New quality check function binMultiNMR2: If reference scaling is used and the reference signal intensity is unexpectedly low, a warning is displayed
* If quality check warnings are generated during a run, the warnings are saved in the parameter text file for later viewing


## Version 1.5.1

* Small improvement in spectrum browsing. Experiments with no title file are no longer displayed during browsing, as these are usually uninformative pulse calibration data
* Added a message recommending installing xquartz on Apple computers
* Citation information was added


## Version 1.5.0

* An error in reference scaling was fixed, this prevented reference scaling in some cases
* PQN scaling is now performed after fixing negative values by atnv() (previously: before atnv). This might slightly change PQN scaling results in rare cases
* Speed improvements by adding optional support for parallel computing using the parallel package


## Version 1.4.4

* Function PQNScaling now works with externally provided data
* Minor improvements for user friendliness


## Version 1.4.3

* Noise levels per bin are no longer saved in the .txt file to keep this file clear and tidy. Instead, raw noise levels and median noise levels (adjusted for bin size) are now saved in the output .txt file


## Version 1.4.2

* Summed bins are now calculated with the exact defined borders. Surrounding bins are truncated
* Noise levels are now calculated separately for each bin in each spectrum, correcting for differences in data points per bin
* PCA plot labels now show percentage of variance instead of percentage of standard deviation
* Some minor and major fixes to increase user friendliness
* Some minor and major bug fixes to increase stability, especially when analyzing single samples and/or single bins


## Version 1.4.1

* Vignette was expanded
* Bin borders are now rounded to avoid floating point inconsistencies
* PQN normalization now also ignores glucose signals in 1D spectra. This behavior can be turned off now as well.
* Some minor and major fixes to increase user friendliness
* Some minor and major bug fixes to increase stability, especially when analyzing single samples and/or single bins


## Version 1.4.0

* Solvent region change: Only bins that are completely within the solvent region are removed
* Excluded regions: No NMR data from excluded regions or solvent regions is now used to calculate bins, even if the bin overlaps with the excluded region
* New trimming function to remove bins that have mostly values of zero; These are created at the edges of the spectrum and at edges of removed regions
* A vignette was added
* Some minor and major fixes to increase user friendliness
* Some minor and major bug fixes to increase stability

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
