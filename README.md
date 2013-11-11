hzz4l-bayesian
==============

This class provides a statistical tool for treating statistical uncertainties in binned data.  If we represent a model as one or more histograms which are created from simulated data, then each bin in each histogram is associated with a statistical uncertainty.   In order to properly account for the effect of these uncertainties, the RooHistPoissonGamma tool performs a bayesian integration over statistical uncertainties. 

This class is written in ROOT using the ROOFIT libraries.  The class is written in C++, and testing scripts are written in python.  It is necessary to have ROOT with RooFit and pyROOT installed. 

Weighted histograms: It is crucial that the event weights be set so that the statistical uncertainties are properly represented in the histograms. Thus, the option Sumw2 must always be called for each template histogram. 
