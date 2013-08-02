#!/usr/bin/env python

#------------------------------------------------------------------------------
# File: likelihoodscan.py
# Created: 24-April-2013 Joe Bochenek
# $Revision$
# Description: Set of H->ZZ->4l likelihood and do a 2D likelihood likelihood
#              scan in m_4l and mu_gg 
#------------------------------------------------------------------------------

import os, sys
from ROOT import *
from time import sleep
import plot_util
import histutil

histutil.setStyle()

# Load Poisson-Gamma RooThing
gROOT.ProcessLine(".L ../src/RooHistPoissonGamma.cxx+")

# S e t u p
# ---------------------------------------------

# Signal Strength Modifiers (define for each signal and background)
mu1 = RooRealVar("mu1","mu1",1.,0,20) 
mu2 = RooRealVar("mu2","mu2",1.,0,20) 
mu3 = RooRealVar("mu3","mu3",1.,0,20) 
mu4 = RooRealVar("mu4","mu4",1.,0,20) 
mu5 = RooRealVar("mu5","mu5",1.,0,20) 
mu6 = RooRealVar("mu6","mu6",1.,0,20) 

# Input parameters
histVar = RooRealVar("histVar","histVar",0., 1000.) 

# Prepare input histograms
fin = TFile( "/home/jbochenek/data/HZZ4l_2013_paper/templates/hzz_templates_0001_00000000.root")

# Set constants
channel = "2e2mu"
era = "8TeV"
mva = "cat1_2dvbf"

masses = ["126"]
stepsize = 0.3
xmin = 0.1
ymin = 0.1
channels = ["2e2mu", "4mu", "4e" ]

mvas = [ "3dto1d_cat1", "3dto1d_cat2_200" ]
#mvas = [ "3dto1d_mass" ]

nullvaltot = 0

c1 = TCanvas()
likelihood_ave = []



for mass in masses:
    likelihoods = []
    
    nullave = 0
    for mva in mvas:
        for channel in channels:

            likelihoods.append(TGraph2D())

            print "mass: {}, plot: {}, channel: {}".format(mass, mva, channel)

            fin.cd()

            # Get qq->ZZ bkg template
            histname = "qqzz_{0}_{1}_{2}".format(channel, era, mva)
            bkg_qqzz = gDirectory.Get(histname)
            print histname
        
            histname = "ggzz_{0}_{1}_{2}".format(channel, era, mva)
            bkg_ggzz = gDirectory.Get(histname)
            print histname


            histname = "ggH_{0}_{1}_{2}_{3}".format(mass, channel, era, mva)
            sig_ggH = gDirectory.Get(histname)    
            print histname
            
            histname = "qqH_{0}_{1}_{2}_{3}".format(mass, channel, era, mva)
            sig_qqH = gDirectory.Get(histname)
            print histname

            # Use VBF for WH and VH 
            sig_VH = sig_qqH.Clone()
            sig_VH.Scale(4.8/7)

            
            # ZX
            histname = "zx_{0}_{1}_{2}".format(channel, era, mva)
            bkg_zx = gDirectory.Get(histname)
            print histname

            histname = "data_{0}_{1}_{2}".format(channel, era, mva)
            data = gDirectory.Get(histname)
            print histname

            print "sig_ggH: ", sig_ggH.Integral()
            print "sig_qqH: ", sig_qqH.Integral()

            print "bkg_qqzz: ", bkg_qqzz.Integral()
            print "bkg_ggzz: ", bkg_ggzz.Integral()
            print "bkg_zjets: ", bkg_zx.Integral()

            print "data: ", data.Integral()
    

            # Convert TH2s to RooDataHists and RooHistPdfs
            hbkg_qqzz =  RooDataHist("bkg_qqzz","bkg_qqzz",  RooArgList(histVar),   bkg_qqzz)
            rbkg_qqzz =  RooHistPdf("rbkg_qqzz", "rbkg_qqzz",   RooArgSet(histVar), hbkg_qqzz)

            hbkg_ggzz =  RooDataHist("bkg_ggzz","bkg_ggzz",  RooArgList(histVar),   bkg_ggzz)
            rbkg_ggzz =  RooHistPdf("rbkg_ggzz", "rbkg_ggzz",   RooArgSet(histVar), hbkg_ggzz)

            rdata  =  RooDataHist("data","data",    RooArgList(histVar), data)

            hsig_qqH = RooDataHist("hsig_qqH","hsig_qqH",RooArgList(histVar),sig_qqH)
            rsig_qqH = RooHistPdf( "rsig_qqH","rsig_qqH",RooArgSet(histVar) ,hsig_qqH)

            hsig_ggH = RooDataHist( "hsig_ggH","hsig_ggH",RooArgList(histVar),sig_ggH)
            rsig_ggH = RooHistPdf(  "rsig_ggH","rsig_ggH",RooArgSet(histVar) ,hsig_ggH)

            hbkg_zx= RooDataHist( "hsig_zx","hsig_zx",RooArgList(histVar),bkg_zx)
            rbkg_zx = RooHistPdf( "rsig_zx","rsig_zx",RooArgSet(histVar), hbkg_zx)

            hsig_VH = RooDataHist("hsig_VH","hsig_VH",RooArgList(histVar),sig_VH)
            rsig_VH = RooHistPdf( "rsig_VH","rsig_VH",RooArgSet(histVar) ,hsig_VH)



            bkgfactor = 1.
            mu3.setVal(bkgfactor)
            mu4.setVal(bkgfactor)
            mu5.setVal(bkgfactor)
        

            null_roohistpg = RooHistPoissonGamma("null_pghistpdf","null_pghistpdf",RooArgSet(histVar), rdata, RooArgList(rbkg_qqzz, rbkg_ggzz, rbkg_zx), RooArgList(mu3       , mu4    , mu5  ), 0)
            nullzval = null_roohistpg.getVal()
            print "Null Likelihood: {}".format(nullzval)
            nullave += 2*nullzval


            samples = RooArgList(rsig_qqH, rsig_ggH,  rbkg_qqzz, rbkg_ggzz, rbkg_zx,  rsig_VH)
            scales  = RooArgList(mu1,      mu2      , mu3      , mu4    , mu5 , mu6 )



            # Do a likelihood scan over mH and mu_gg
            for binmu_vbf in xrange(80):
                for binmu_gg in xrange(20):
                    mu_gg_point  = xmin + binmu_gg * 0.1
                    mu_vbf_point  = ymin + binmu_vbf * 0.1

                    
                    mu2.setVal(mu_gg_point) 
                    mu1.setVal(mu_vbf_point)        
                    mu6.setVal(mu_vbf_point)        

                    #, RooArgList(mu1, mu2)
                    roohistpg = RooHistPoissonGamma("pghistpdf","pghistpdf",RooArgSet(histVar), rdata, samples, scales, 0)
                    zval = 2*roohistpg.getVal() # nullzval - 
                    if zval > 0:
                        likelihoods[len(likelihoods)-1].SetPoint( likelihoods[len(likelihoods)-1].GetN(),  mu_gg_point,  mu_vbf_point, zval);
                    else:
                        likelihoods[len(likelihoods)-1].SetPoint( likelihoods[len(likelihoods)-1].GetN(),  mu_gg_point,  mu_vbf_point, zval);

            likelihoods[len(likelihoods)-1].Draw("colz")
            c1.Update()

    #dt = plot_util.multiply_likelihood(likelihoods)
    dt = plot_util.add_loglikelihood(likelihoods, 1)
    likelihood_ave.append(dt)
    nullvaltot += nullave/len(likelihood_ave)

dt = plot_util.ave_ln_loglikelihood(likelihood_ave)
#reset to zero 


print "Total Null Likelihood {}".format(nullvaltot)

signif = plot_util.significance(dt, nullvaltot)

signif.SetTitle("")
signif.GetXaxis().SetTitle("#mu_{F}")
signif.GetYaxis().SetTitle("#mu_{V}")
signif.GetZaxis().SetTitle("Significance")
signif.GetZaxis().SetTitleOffset(1.3)

# Print 2d posteror
c = TCanvas()
c.cd()
signif.Draw("colz")
c.SaveAs("../output/significance_muv_muf.png")
c.SaveAs("../output/significance_muv_muf.pdf")


dt = plot_util.setzero_likelihood(dt)


c = TCanvas()
c.Divide(2,2)
c.cd(1)
likelihoods[0].Draw("colz")
c.cd(2)
likelihoods[1].Draw("colz")
c.cd(3)
likelihoods[2].Draw("colz")
c.cd(4)
dt.Draw("colz")
c.SaveAs("../output/nllscan_muvbf_mugg_ind.png")
    
    
c = TCanvas()
c.cd()
dt.Draw("cont2")
c.SaveAs("../output/nllscan_muvbf_mugg_contz.png")

# Write stuff out
fout = TFile("../output/lh_mu_output.root", "recreate")
fout.cd()
dt.SetName("nll_mu_gg_vs_mH")
dt.SetTitle("")
dt.GetXaxis().SetTitle("#mu_{gg}")
dt.GetYaxis().SetTitle("#mu_{VBF}")
dt.GetXaxis().SetTitleOffset(0.5)
dt.GetYaxis().SetTitleOffset(0.5)
dt.GetXaxis().SetTitleSize(0.06)
dt.GetYaxis().SetTitleSize(0.06)
dt.GetHistogram().SetContour(11)
dt.Write()


c = TCanvas()
c.cd()
dt.Draw("colz")
dt.Draw("cont2 same")
c.SaveAs("../output/nllscan_muvbf_mugg.png")

fout.Close()
fin.Close()