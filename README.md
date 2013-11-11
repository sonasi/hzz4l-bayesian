hzz4l-bayesian
==============

We wish to represent our real and simulated events by one or more histograms in variables which are sensitive to our parameters of interest. I refer to these histograms as templates. To fill the histograms for the various signal, background, and data samples, there is a (rather ugly) python script:

The RooHistPoissonGamma tool performs a bayesian integration over statistical uncertainties. It is therefore crucial that the event weights be set so that the statistical uncertainties are properly represented in the templates. Thus, the option Sumw2 must always be called for each template histogram. 
