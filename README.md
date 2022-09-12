# Antiproton_Analysis_Scripts
A series of scripts used in antiproton analysis resulting from MDRange simulations

This is a collection of files used in antiproton analysis in conjunction with MDRange.

Due to licensing, it is not possible to provide example MDRange scripts nor the code. Please contact the author of MDRange for the program,
I may be able provide examples privately if contacted and given permission by the program author.

Elstop_Fit.py is crucial to the running of MDRange, it produces an electronic stopping file suitable for use within MDRange, by using previously measured
stopping powers of antiprotons within given materials. Used in conjunction with the DFT routine (also available on github.com/stevepadden/Antiproton_DFT_Setup)
this provides both a nuclear and electronic stopping power to MDRange, required for an accurate description.

Range_profile plots the initial range profile of a 100 keV antiproton beam returned from MDRange into an infintiley thick foil.
It then does the same for a 5 keV beam, subtracts the modal centers and returns a suggested start width used for finding the optimised thickness.

By running multiple MDRange simulations of foil thickness' around this "suggested" thickness, a collection of measurements
can be made which shows how each foil degrades a 100 kev beam. Width_Find uses this collection of measurements to produce a gaussian
fitting routine, of which the centroid is the optimised thickness.

Armed with the optimised thickness, an MDRange simulation can be made at this value. Once it has been completed, Analyse.py 
will return a number of useful graphs, it also uses bootstrap resampling to provide error's resulting from multiple runs.

Pbar_Proton_Heatmap can be used to return a spatial heatmap of how differently protons and antiprotons behave in optimised foils.

Finally Pkl2Track.py converts the particles which traverse a foil of a given thickness, into an input beam usable by G4Beamline.
