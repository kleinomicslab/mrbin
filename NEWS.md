# mrbin News

## Version 1.9.0

* mrbin: Nicer and speedier plots
* mrbin: Improved bin annotation system
* mrplot: A vector of mrbin bin names can now be provided to mrplot and the bin positions will be highlighted
* mrplot: Small bugfixes
* mrbin: Small bugfixes
* mrbin: Clarifications of prompts and reply options
* mrbin: Clarifications in the vignette
* mrheatmap: New function for heatmaps
* examples: A new data set was included for better example results


## Version 1.8.0

* mrbin: Speed improvements
* mrbin: Pop-up windows for list selection can now be turned off using graphics=FALSE 
* plotNMR: Improvements for quicker plot display - For full spectra, only a subset of the data points is plotted
* mrplot: Improved plot manipulation, e.g. zooming in/out, moving left/right
* mrplot: Added options for adjusting the plot, e.g. setting a background color
* mrplot: Small bug fixes
* mrbin: Small bug fixes and clarifications of prompts and reply options


## Version 1.7.5

* mrbin: Bug fixes especially for plots in (Mac) RStudio
* mrbin: The number of bins left after trimming zero-value bins was displayed incorrectly during an mrbin run, this is now fixed
* mrbin: All function examples now use parallel computing by default for consistency
* mrbin: Small bug fixes and clarifications of prompts and reply options
* mrbin: Bin data is now additionally saved as a csv file to disk (again)


## Version 1.7.4

* mrbin: Bugfixes especially for plots in (Mac) RStudio


## Version 1.7.3

* mrbin: Warnings are now displayed during binning
* mrbin: New warnings if spectra differ in field strength, pulse program, solvent, or number of scans. This could cause inconsistent data
* mrbin: Some acquisition parameters such as p1 pulse length are now saved in parameters$AcquPars
* mrbin: Small potential bugs fixed in spectrum selection step
* mrbin: Bugfix: In rare cases the preview of features after noise removal was incorrect, this is now fixed
* fia: Predictions now show less output (using verbose=0 within the predict function)
* fia: Memory usage improvements. If using package keras, consider using kerasClearMemory=2 for massive memory and speed improvements
* fia: Parameter innerLoop is now default at 300 (was 100) to increase quality of results. This increases calculation times. 


## Version 1.7.2

* In setNoiseLevels() there is now a preview of the number of bins left for noise thresholds between 0.05 and 1
* Renaming of duplicated spectrum names was improved, this also fixed a rare bug that could cause mrbin to get stuck while binning
* Previews of all spectra can now be reviewed while running mrbin() to identify quality issues
* Speed improvements in 2D plots by reducing the number of displayed data points
* plotResults() will now show a preview of spectral region evens if no spectral data has been loaded before. This is useful when using mrbin objects that have been previously created
* Small improvements to documentations and vignette


## Version 1.7.1

* Bug fix: annotatemrbin() sometimes added extra commas in metabolite annotations, this has been fixed
* Bug fix: checkmrbin(), in rare cases, erroneously flagged unstated changes to an unchanged mrbin object. This was caused by numeric instabilities and has now been fixed
* Bug fix: mrbin() did not properly add metadata to the result object, if added at the command line. This is now fixed
* New functions setDilutionFactors() and dilutionCorrection() for individual dilution correction
* Variable metaboliteIdentities in metadata was changed to matrix (was data.frame)
* printParameters() now omits empty (NULL, "", nrow==0) variables, but prints non-empty metadata now
* New function unitVarianceScaling()


## Version 1.7.0

* mrbin: results are now returned as an object of the new class "mrbin", containing bin data, metadata, and change logs. These objects should only be changed using functions from this package. Individual edits can be performed using editmrbin()
* mrbin: Noise levels are now visually shown for a few spectra to help pick reasonable signal-to-noise levels
* mrbin: Warnings during binning now cause a pop-up message to give users the chance to review and fix data issues such as phase, reference, and baseline errors
* mrbin: PCA plots do no longer show group memberships during the data processing phase. To display color-coded groups, add metadata using mrbinResults<-metadatamrbin(mrbinResults) and then use plotPCA(mrbinResults)
* mrbin: Changes to mrbin objects are now logged in $changeLog to ensure workflow reproducibility. Undocumented changes can be identified using checkmrbin(mrbinResults)
* mrbin: New functionality to annotate mrbin data with molecule identities was added
* mrbin: Improvements in memory usage to better handle large data sets. If more than 5000 bins are created, these will not be plotted in the output plot.
* mrbin: The variable noise_level was moved to parameters list in mrbin environment
* mrbin: Bug fix: Excluding the sugar area during PQN could not be turned off within mrbin, this is now fixed
* mrbin: Users will no longer be asked whether to use parallel, previews and hints, these will be always used now unless turned off by manually setting the respective parameters
* mrbin: Time estimates are now being displayed before and during binning


## Version 1.6.5

* New output of fia function: scoresIndividual shows fia score for each tested sample
* Speed improvements in the fia function



## Version 1.6.4

* Minor changes in calculations in the fia function. Features of lowest impact now are assigned a (large) FIA score instead of an NA value


## Version 1.6.3

* Small bug fixes in the fia function
* Small bug fixes in the mrplot function


## Version 1.6.2

* New function fia for analyzing predictive metabolomics models
* Small bug fixes in the mrplot function


## Version 1.6.1

* Small bug fixes in the mrplot function


## Version 1.6.0

* Bin intensities are now pseudo-integrals (mean of data points times range (1D), or times area (2D)). So far, intensities were the mean of data points in the bin range. For simple rectangular bins, results are identical apart from a scaling factor, so most data analysis methods will give identical results. To enforce the old behavior, use mrbin(parameters=list(useMeanIntensityForBins=TRUE))
* Multiple spectra are now shown in previews to visualize sample-to-sample variability 
* PCA plots: When using more groups than available symbols (from pch), letter are now used to mark the groups. This also removes the warning messages generated in these cases.
* Reading spectra without title does no longer create an error in readBruker(), however, such spectra are still not shown in the spectrum browser as they are usually not of interest
* Before binning, folder list is checked to find missing folders. This is to avoid long calculation that end with a "could not open file" error. This uses the new readBruker() parameter checkFiles.
* Speed improvements for displaying and saving the final results plot 
* New plotting function mrplot()
* RStudio plotting issue solved: Plots are now refreshed so the current plot is displayed timely


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
