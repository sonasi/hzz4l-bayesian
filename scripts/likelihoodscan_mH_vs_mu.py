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

from array import array

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

# Input parameters
m4l = RooRealVar("m4l","m4l",115, 180) 
D = RooRealVar("bnn","bnn",0, 1) 



# Prepare input histograms
fin = TFile( "/home/jbochenek/data/HZZ4l_2013_paper/templates/hzz_templates_0000_00000000.root")

# Set constants
channel = "2e2mu"
era = "8TeV"
mva = "2d"

masses = ["118", "120", "122", "123", "124", "125", "126", "127", "128", "130", "135"]
massint = map(int, masses)
stepsize = 0.1
xmin = 0.1
ymin = 0.1
points = 30
channels = [ "4e", "2e2mu", "4mu" ]

likelihoods = []
likelihoods1d = []

fout = TFile("../output/lh_mh_output.root", "recreate")

c1 = TCanvas()
c1.Divide(2,2)

nullvaltot = 0


c_1d = TCanvas()
c_1d.cd()


# Sanity Check: print 1d histograms with data
hb = THStack("hb","")
hd = THStack("hd","")
hs = THStack("hs","")


for channel in channels:
    fin.cd()
    # ZX
    histname = "zx_{0}_{1}_{2}".format(channel, era, mva)
    bkg_zx = gDirectory.Get(histname)
    print histname
    bkg_zx.SetFillColor(kGreen)
    bkg_zx.SetLineColor(1)   
    hb.Add(bkg_zx.ProjectionX())

for channel in channels:
    fin.cd()
    # Get qq->ZZ bkg template
    histname = "qqzz_{0}_{1}_{2}".format(channel, era, mva)
    bkg_qqzz = gDirectory.Get(histname)
    histname = "ggzz_{0}_{1}_{2}".format(channel, era, mva)
    bkg_ggzz = gDirectory.Get(histname)
    bkg_qqzz.SetLineColor(1)
    bkg_qqzz.SetFillColor(kBlue)
    bkg_ggzz.SetFillColor(kBlue)
    hb.Add(bkg_qqzz.ProjectionX())
    hb.Add(bkg_ggzz.ProjectionX())

for channel in channels:
    fin.cd()
    mass = "126"
    sig = gDirectory.Get( "ggH_" + str(mass) + "_" + channel + "_" + era + "_" + mva   )
    sig.SetLineColor(kRed)
    print sig.Integral()
    hb.Add(sig.ProjectionX())
    sig2 = gDirectory.Get( "qqH_" + str(mass) + "_" + channel + "_" + era + "_" + mva   )
    sig2.SetLineColor(kRed)
    print sig2.Integral()

for channel in channels:
    fin.cd()
    # Get the data    
    histname = "data_{0}_{1}_{2}".format(channel, era, mva)
    data = gDirectory.Get(histname)
    data.ProjectionX().SetMarkerStyle(20)
    data.ProjectionX().SetMarkerSize(1)
    data.SetMarkerStyle(20)
    data.SetMarkerSize(1)
    hd.Add(data.ProjectionX(), "p")
    
hd.Draw("hist e0p")
hb.Draw("hist same")
hd.Draw("hist e0p same")

print hd.GetHistogram().Integral()

c_1d.Update()
c_1d.SaveAs("../output/masscheck.png")

print "PRINTING 1D m4l"
sleep(10)



for channel in channels:
    fin.cd()
    likelihoods.append(TGraph2D())

    # Get qq->ZZ bkg template
    histname = "qqzz_{0}_{1}_{2}".format(channel, era, mva)
    bkg_qqzz = gDirectory.Get(histname)
    print histname

    histname = "ggzz_{0}_{1}_{2}".format(channel, era, mva)
    bkg_ggzz = gDirectory.Get(histname)
    print histname

    # ZX
    histname = "zx_{0}_{1}_{2}".format(channel, era, mva)
    bkg_zx = gDirectory.Get(histname)
    print histname

    # Get the data    
    histname = "data_{0}_{1}_{2}".format(channel, era, mva)
    print histname
    data = gDirectory.Get(histname)

    

    # Convert TH2s to RooDataHists and RooHistPdfs
    hbkg_qqzz =  RooDataHist("bkg_qqzz","bkg_qqzz",  RooArgList(m4l,D),   bkg_qqzz)
    rbkg_qqzz =  RooHistPdf("rbkg_qqzz", "rbkg_qqzz",   RooArgSet(m4l,D), hbkg_qqzz)

    hbkg_ggzz =  RooDataHist("bkg_ggzz","bkg_ggzz",  RooArgList(m4l,D),   bkg_ggzz)
    rbkg_ggzz =  RooHistPdf("rbkg_ggzz", "rbkg_ggzz",   RooArgSet(m4l,D), hbkg_ggzz)

    rdata  =  RooDataHist("data","data",    RooArgList(m4l,D), data)

    hbkg_zx =  RooDataHist("bkg_zx","bkg_zx",  RooArgList(m4l,D),   bkg_zx)
    rbkg_zx =  RooHistPdf("rbkg_zx", "rbkg_zx",   RooArgSet(m4l,D), hbkg_zx)


    # Compute likelihood for 'null' hypothesis
    null_roohistpg = RooHistPoissonGamma("null_pghistpdf","null_pghistpdf",RooArgSet(m4l, D), rdata, RooArgList(rbkg_qqzz, rbkg_ggzz, rbkg_zx), RooArgList(mu3       , mu4    , mu5  ), 0)
    nullzval = null_roohistpg.getVal()
    print "Null Likelihood: {}".format(nullzval)
    nullvaltot += 2*nullzval

    scanmu1d = [0]*len(masses)

    for i, mass in enumerate(masses):
            print "Channel: {}, Mass: {}".format(channel, mass)
    
            # Get the signal
            sig_ggH = gDirectory.Get(    "ggH_" + str(mass) + "_" + channel + "_" + era + "_" + mva   )    
            sig_qqH = gDirectory.Get(    "qqH_" + str(mass) + "_" + channel + "_" + era + "_" + mva   )
    
            # Convert TH2s to RooDataHists and RooHistPdfs
            hsig_qqH = RooDataHist("hsig_qqH","hsig_qqH",RooArgList(m4l,D),sig_qqH)
            rsig_qqH = RooHistPdf( "rsig_qqH","rsig_qqH",RooArgSet(m4l,D) ,hsig_qqH)

            hsig_ggH = RooDataHist("hsig_ggH","hsig_ggH",RooArgList(m4l,D),sig_ggH)
            rsig_ggH = RooHistPdf( "rsig_ggH","rsig_ggH",RooArgSet(m4l,D) ,hsig_ggH)

            samples = RooArgList(rsig_qqH, rsig_ggH,  rbkg_qqzz, rbkg_ggzz, rbkg_zx)
            scales  = RooArgList(mu1,      mu2      , mu3       , mu4    , mu5  )

            # Do a likelihood scan over mH and mu_gg
            for binmu_gg in xrange(points):
                mu_gg_point  = xmin + binmu_gg * stepsize
                mu2.setVal(mu_gg_point)
                mu1.setVal((1.7)*mu_gg_point)        

                #, RooArgList(mu1, mu2)
                roohistpg = RooHistPoissonGamma("pghistpdf","pghistpdf",RooArgSet(m4l, D), rdata, samples, scales, 0)
                zval = 2*roohistpg.getVal() # 2*(nullzval - roohistpg.getVal())
                if zval > 0:
                    likelihoods[len(likelihoods)-1].SetPoint( likelihoods[len(likelihoods)-1].GetN(), float(mass), mu_gg_point, zval);
                else:
                    likelihoods[len(likelihoods)-1].SetPoint( likelihoods[len(likelihoods)-1].GetN(), float(mass), mu_gg_point, zval);

            mu2.setVal(1.)        
            mu1.setVal(1.7)        
            zval = 2*roohistpg.getVal() # 2*(nullzval - roohistpg.getVal())
            scanmu1d[i] += zval

    for i, mass in enumerate(masses):
        scanmu1d[i] = scanmu1d[i]/points
    print len(masses)
    print len(massint)
    print massint
    print scanmu1d    
    gr = TGraph(len(masses), array("d", massint) , array("d",scanmu1d) )
    gr.SetName("nll_mu_mH_{}".format(channel))
    c1.cd(1)
    likelihoods[len(likelihoods)-1].Draw("colz")
    c1.cd(2)
    gr.Draw("AC*")
    c1.cd(3)
    hd.Draw("lep")
    hb.Draw("hist same")
    c1.Update()
    sleep(4)
    fout.cd()
    gr.Write()
    likelihoods1d.append(gr)


#dt = plot_util.multiply_likelihood(likelihoods)
dt = plot_util.add_loglikelihood(likelihoods, 1)

print "Total Null Likelihood: {}".format(nullvaltot)

signif = plot_util.significance(dt, nullvaltot)

signif.SetTitle("")
signif.GetXaxis().SetTitle("m_{H} (GeV)")
signif.GetYaxis().SetTitle("#mu")
signif.GetZaxis().SetTitle("Significance")
signif.GetZaxis().SetTitleOffset(1.3)

# Print 2d posteror
c = TCanvas()
c.cd()
signif.Draw("colz")
c.SaveAs("../output/significance_m4l_mugg.png")
c.SaveAs("../output/significance_m4l_mugg.pdf")


#reset to zero 
dt = plot_util.setzero_likelihood(dt)

dt1d = plot_util.add_tgraph(likelihoods1d)
dt1d.SetName("nll_mH_1d")
dt1d.Write()

dt1d.Draw()
sleep(4)
    
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
c.SaveAs("../output/nllscan_m4l_mugg_ind.png")


    

dt.SetName("nll_mu_gg_vs_mH")
dt.SetTitle("")
dt.GetXaxis().SetTitle("m_{4l}")
dt.GetYaxis().SetTitle("#mu_{gg}")
dt.GetHistogram().SetContour(  25  )
dt.Write()

    
c = TCanvas()
c.cd()
dt.Draw("colz")
dt.Draw("cont2 same")
c.SaveAs("../output/nllscan_m4l_mugg_contz.png")
sleep(2)


c = TCanvas()
c.cd()
dt.Draw("colz")
c.SaveAs("../output/nllscan_m4l_mugg.png")


fout.Close()
fin.Close()



