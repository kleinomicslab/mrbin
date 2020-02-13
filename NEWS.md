# mrbin News


## Version 1.0.0.9000

* All parameters can now be conveniently set from the command line by passing a list of parameters to mrbin()
* Special bin regions can now be defined instead of rectangular bins, e.g. for lipid analysis
* mrbin() now returns the group factors as an additional element "factors" in the returned list
* New function setCurrentSpectrum() to set current spectrum, e.g. for plotting 
* Region borders for each bin are now saved in a separate matrix 
* Reloading a parameter set from a file now checks whether all required parameters have been set and adds missing ones
* Some minor fixes to increase user friendliness
* Some minor bug fixes to increase stability
