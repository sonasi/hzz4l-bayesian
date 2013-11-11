#!/usr/bin/env python

#------------------------------------------------------------------------------
# File: test_pgamma_hist.txt
# Created: 11-Nov-2013 Joe Bochenek
# $Revision$
# Description: This provides a simple example of the RooHistPoissonGamma class
#			   in the context of a mass search.  We generate fake data at a 
# 			    particular mass and then perform a mass scan in one dimension.
#------------------------------------------------------------------------------


import os, sys
from ROOT import *
from time import sleep

from array import array

# Load Poisson-Gamma RooThing
gROOT.ProcessLine(".L ../src/RooHistPoissonGamma.cxx+")

