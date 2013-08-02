#!/usr/bin/env python
#-----------------------------------------------------------------------------
# Example of plotting systematics as confidence bands
# Created: sometime this century HBP
# Updated: 20 Apr. 2013 HBP implement polygon clipping in mkpline
#-----------------------------------------------------------------------------
import os, sys, re
from ROOT import *
from string import *
from histutil import *
from array import array
from time import sleep
#-----------------------------------------------------------------------------
def main():
	setStyle()

	# select boundary of histogram in order to
	# test polygon clipping in mkpline
	
	xmin = 105     # Signal
	xmax = 155
	ymin = 0
	ymax = 0.6
#	xmin = 80      # Background
#	xmax = 200
#	ymin = 0
#	ymax = 0.12
	boundary = (xmin, xmax, ymin, ymax)

	# Use previously-saved systematics plots
	filename = "plots.root"
	nhists = 100        # Number of histograms in ROOT file
	rfile = TFile(filename)
	
	if not rfile.IsOpen():
		print "can't open %s" % filename
		sys.exit(0)

	plot0 = rfile.Get("plot0")
	nbins= plot0.GetNbinsX()
	
	count = plot0.Integral()
	scale = 1/count

	print "count:", count
	print "scale:", scale

	plot0.Scale(scale)
	
#	count= data.Integral()
#	xsect= c000.Integral()
#	scale= xsect/count
#
#	# scale data to QCD prediction
#	data.Scale(scale)
#
#	print 
#	print "count: %7d" % count
#	print "xsect: %10.2f pb" % xsect
#	print "Number of bins: %d" % nbins

	# get bin centers
	pt = array('d')
	hpt= plot0.Clone("hpt") # create histogram of bin widths
	for ii in xrange(nbins):
		pT = plot0.GetBinLowEdge(ii+1) + 0.5*plot0.GetBinWidth(ii+1)
		pt.append(pT)
		hpt.SetBinContent(ii+1, plot0.GetBinWidth(ii+1))

#	# divide by bin width
#	data.Divide(hpt)
	
	# ----------------------------------
	# Add histograms to percentive curve
	# ----------------------------------
	pc = PercentileCurve(nbins)
	
	# loop over histograms
	for ii in xrange(nhists):
		nn = ii+1
		hname = "plot%d" % nn
		h = rfile.Get(hname)
		if h == None:
			print "*** can't find %s" % hname
			sys.exit(0)
		h.Scale(scale)
		
#		# divide by bin width to get dsigma/dp_T
#		h.Divide(hpt)
		pc.add(h)
#		if nn % 100 == 0:
#			print hname
			
	# for each percentile, create an array of points
	# representing the y values of the percentile curve
	curve = []
	for p in PERCENT:
		c = array('d')
		c.fromlist( pc(p) )
		curve.append(c)
		
	# ----------------------------
	# Now plot
	# ----------------------------
	cp = TCanvas("fig_spectrum", "spectrum", 10, 10, 500, 500)
	cp.cd()
#	cp.SetLogy()
	
	# create 95% C.L. band
	pl95 = mkpline(pt, curve[0], curve[-1], boundary, color=kYellow)

	# create 68% C.L. band
	pl68 = mkpline(pt, curve[1], curve[-2], boundary, color=kGreen)

	# set plot style and boundary
	htmp = mkhist1("htmp", "m_{4L} (GeV)",
		       "d#sigma/dm_{4L} (norm.)", 100,
		       xmin, xmax, ymin=ymin, ymax=ymax)
	htmp.Draw()

	# plot 95% band first
	pl95.Draw("f")

	# superimpose 68% band
	pl68.Draw("f")

	# superimpose QCD curve
	plot0.SetLineColor(kRed+1)
	plot0.SetLineWidth(2)
	plot0.Draw("hist c same")

	# superimpose data
#	data.Draw("e same")
#	# add horizontal "error" bars to indicate bin widths
#	setex = TExec("setex","gStyle->SetErrorX(0.5)")
#	setex.Draw()

	scribe = addTitle("CMS Simulation                 #sqrt{s} = 8 TeV")
	scribe.vspace(0.8)
#	scribe.write('|#eta| < 0.5', 0.2)
	
	# add legend
	xx = 0.60
	yy = 0.60
	xw = 0.28
	yw = 0.30
	lg = mklegend(xx, yy, xw, yw)
#	lg.AddEntry(data, "Data", "ep")
	lg.AddEntry(plot0, "No Systematics", "l")
	lg.AddEntry(pl68, "68% C.L.", "f")
	lg.AddEntry(pl95, "95% C.L.", "f")
	lg.Draw()
	
	cp.Update()
	cp.SaveAs(".png")
#	cp.SaveAs(".pdf")

	# NB: Use Canvas menu to exit cleanly!
	gApplication.Run()
#-----------------------------------------------------------------------------
try:
	main()
except KeyboardInterrupt:
	print
	print "\tBye!"

