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

import glob

histutil.setStyle()

from array import array

# Load Poisson-Gamma RooThing
gROOT.ProcessLine(".include /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.02/x86_64-slc5-gcc43-opt/root/include/")
gROOT.ProcessLine(".L ../src/RooHistPoissonGamma.cxx+")

# Read seed from command line
# If the seed number is 0 then no systematics are applied
argv =  sys.argv[1:]
if len(argv)==0:
	TOYMASS= 126
	JOBNUM= 1234
elif len(argv)==1:
	TOYMASS = atoi(argv[0])
	JOBNUM= 1234
else:
	TOYMASS = atoi(argv[0])
	JOBNUM = atoi(argv[1])

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
histVar = RooRealVar("histVar","histVar",0., 1.) 


f = open("/lustre/cms/store/user/jpb/output/lh_mh_pull_{0}.txt".format(TOYMASS), "w")

# Prepare input histograms

# Set constants
channel = "4mu"
era = "8TeV"

masses = ["118", "120", "122", "123", "124", "125", "126", "127", "128", "130", "135"]
massint = map(int, masses)
stepsize = 0.1
xmin = 0.1
ymin = 0.1
points = 40
channels = [ "2e2mu", "4mu", "4e"]

likelihoods = []
likelihoods1d = []

fout = TFile("/lustre/cms/store/user/jpb/output/lh_mh_pull_{0}.root".format(TOYMASS), "recreate")

c1 = TCanvas()
c1.Divide(2,2)

nullvaltot = 0





mvas = [ "3dto1d_cat1", "3dto1d_cat2" ]
mvas = [ "3dto1d_mass" ]


systfiles = glob.glob('/lustre/cms/store/user/jpb/templates/syst/hzz_templates_*.root')
likelihoods_iter = []

signif_ave = 0

for iter2 in range(5):
    for iter, systfile in enumerate(systfiles):
        nullvaltot = 0

        print "Iteration: {0}/{1}, file {2}".format( iter, systfile, len(systfiles))
        fin = TFile( systfile )
        print systfile 
    
        fin.cd()

        likelihoods = []

        for mva in mvas:
            for channel in channels:

                fin.cd()
                likelihoods.append(TGraph2D())

                # Get qq->ZZ bkg template
                histname = "qqzz_{0}_{1}_{2}".format(channel, era, mva)
                bkg_qqzz = gDirectory.Get(histname)
                print "{0} Integral: {1}".format(histname, bkg_qqzz.Integral())

                histname = "ggzz_{0}_{1}_{2}".format(channel, era, mva)
                bkg_ggzz = gDirectory.Get(histname)
                print "{0} Integral: {1}".format(histname, bkg_ggzz.Integral())
    
    

                # ZX
                histname = "zx_{0}_{1}_{2}".format(channel, era, mva)
                bkg_zx = gDirectory.Get(histname)
                print "{0} Integral: {1}".format(histname, bkg_zx.Integral())

                # Get the data    
                histname = "data_{0}_{1}_{2}".format(channel, era, mva)
                data = gDirectory.Get(histname)
                print "{0} Integral: {1}".format(histname, data.Integral())

    

                # Convert TH2s to RooDataHists and RooHistPdfs
                hbkg_qqzz =  RooDataHist("bkg_qqzz","bkg_qqzz",  RooArgList(histVar),   bkg_qqzz)
                rbkg_qqzz =  RooHistPdf("rbkg_qqzz", "rbkg_qqzz",   RooArgSet(histVar), hbkg_qqzz)

                hbkg_ggzz =  RooDataHist("bkg_ggzz","bkg_ggzz",  RooArgList(histVar),   bkg_ggzz)
                rbkg_ggzz =  RooHistPdf("rbkg_ggzz", "rbkg_ggzz",   RooArgSet(histVar), hbkg_ggzz)

                rdata  =  RooDataHist("data","data",    RooArgList(histVar), data)

                hbkg_zx =  RooDataHist("bkg_zx","bkg_zx",  RooArgList(histVar),   bkg_zx)
                rbkg_zx =  RooHistPdf("rbkg_zx", "rbkg_zx",   RooArgSet(histVar), hbkg_zx)

                scanmu1d = [0]*len(masses)
            
                # Get the signal
                sig_ggH = gDirectory.Get(    "ggH_" + str(TOYMASS) + "_" + channel + "_" + era + "_" + mva   )    
                sig_qqH = gDirectory.Get(    "qqH_" + str(TOYMASS) + "_" + channel + "_" + era + "_" + mva   )

                # Use VBF for WH and VH 
                sig_VH = sig_qqH.Clone()
                sig_VH.Scale(4.8/7)

                # Convert TH2s to RooDataHists and RooHistPdfs
                hsig_qqH = RooDataHist("hsig_qqH","hsig_qqH",RooArgList(histVar),sig_qqH)
                rsig_qqH = RooHistPdf( "rsig_qqH","rsig_qqH",RooArgSet(histVar) ,hsig_qqH)

                hsig_ggH = RooDataHist("hsig_ggH","hsig_ggH",RooArgList(histVar),sig_ggH)
                rsig_ggH = RooHistPdf( "rsig_ggH","rsig_ggH",RooArgSet(histVar) ,hsig_ggH)

                hsig_VH = RooDataHist("hsig_VH","hsig_VH",RooArgList(histVar),sig_VH)
                rsig_VH = RooHistPdf( "rsig_VH","rsig_VH",RooArgSet(histVar) ,hsig_VH)

                samples = RooArgList(rsig_qqH, rsig_ggH,  rbkg_qqzz, rbkg_ggzz, rbkg_zx,  rsig_VH)
                scales  = RooArgList(mu1,      mu2      , mu3      , mu4      , mu5    ,  mu6 )


                #, RooArgList(mu1, mu2)
                roohistpg = RooHistPoissonGamma("pghistpdf","pghistpdf",RooArgSet(histVar), rdata, samples, scales, 0)

                zval = roohistpg.genHist(JOBNUM*iter2)
                rdata  =  roohistpg.dataHist()

                hdata = rdata.createHistogram("fakedata_"+channel, histVar, RooFit.Binning(500) )
                hdata.SetName("fakedata_{0}_{1}".format(TOYMASS,JOBNUM))

                fout.cd()
                hdata.Write()

                fin.cd()
            
            
                mu2.setVal(0.)
                mu1.setVal(0.)        
                mu6.setVal(0.)  

                nullzval = roohistpg.getVal()
                print "Null Likelihood: {0}".format(nullzval)
                nullvaltot += 2*nullzval


                for i, mass in enumerate(masses):
                        print "Channel: {0}, Mass: {1}".format(channel, mass)
    
                        # Get the signal
                        sig_ggH = gDirectory.Get( "ggH_" + str(mass) + "_" + channel + "_" + era + "_" + mva   )    
                        sig_qqH = gDirectory.Get( "qqH_" + str(mass) + "_" + channel + "_" + era + "_" + mva   )
    
                        # Use VBF for WH and VH 
                        sig_VH = sig_qqH.Clone()
                        sig_VH.Scale(4.8/7)

    
                        print "ggH_" + str(mass) + "_" + channel + "_" + era + "_" + mva   
                        print "Integral: {0}".format(sig_ggH.Integral())
                        print "qqH_" + str(mass) + "_" + channel + "_" + era + "_" + mva   
                        print "Integral: {0}".format(sig_qqH.Integral())

                        # Convert TH2s to RooDataHists and RooHistPdfs
                        hsig_qqH = RooDataHist("hsig_qqH","hsig_qqH",RooArgList(histVar),sig_qqH)
                        rsig_qqH = RooHistPdf( "rsig_qqH","rsig_qqH",RooArgSet(histVar) ,hsig_qqH)

                        hsig_ggH = RooDataHist("hsig_ggH","hsig_ggH",RooArgList(histVar),sig_ggH)
                        rsig_ggH = RooHistPdf( "rsig_ggH","rsig_ggH",RooArgSet(histVar) ,hsig_ggH)

                        hsig_VH = RooDataHist("hsig_VH","hsig_VH",RooArgList(histVar),sig_VH)
                        rsig_VH = RooHistPdf( "rsig_VH","rsig_VH",RooArgSet(histVar) ,hsig_VH)

                        samples = RooArgList(rsig_qqH, rsig_ggH,  rbkg_qqzz, rbkg_ggzz, rbkg_zx,  rsig_VH)
                        scales  = RooArgList(mu1,      mu2      , mu3      , mu4    , mu5 , mu6 )



                        # Do a likelihood scan over mH and mu_gg
                        for binmu_gg in xrange(points):
                            mu_gg_point  = xmin + binmu_gg * stepsize
                            mu2.setVal(mu_gg_point)
                            mu1.setVal(mu_gg_point)        
                            mu6.setVal(mu_gg_point)        

                            #, RooArgList(mu1, mu2)
                            roohistpg = RooHistPoissonGamma("pghistpdf","pghistpdf",RooArgSet(histVar), rdata, samples, scales, 0)
                            zval = 2*roohistpg.getVal() # 2*(nullzval - roohistpg.getVal())
                            if zval > 0:
                                likelihoods[len(likelihoods)-1].SetPoint( likelihoods[len(likelihoods)-1].GetN(), float(mass), mu_gg_point, zval);
                            else:
                                likelihoods[len(likelihoods)-1].SetPoint( likelihoods[len(likelihoods)-1].GetN(), float(mass), mu_gg_point, zval);

                        mu2.setVal(1.)        
                        mu1.setVal(1.)        
                        mu6.setVal(1.)        

                        zval = 2*roohistpg.getVal() # 2*(nullzval - roohistpg.getVal())
                        scanmu1d[i] += zval

                for i, mass in enumerate(masses):
                    scanmu1d[i] = scanmu1d[i]/points
                print len(masses)
                print len(massint)
                print massint
                print scanmu1d    
                gr = TGraph(len(masses), array("d", massint) , array("d",scanmu1d) )
                gr.SetName("nll_mu_mH_{0}".format(channel))
                c1.cd(1)
                likelihoods[len(likelihoods)-1].Draw("colz")
                c1.cd(2)
                gr.Draw("AC*")
                c1.cd(3)
                c1.Update()
                fout.cd()
                gr.Write()
                likelihoods1d.append(gr)


        #dt = plot_util.multiply_likelihood(likelihoods)
        dt = plot_util.add_loglikelihood(likelihoods, 1)

        print "Total Null Likelihood: {0}".format(nullvaltot)

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

        signif.SetName("signif_mu_gg_vs_mH_{0}".format(iter+iter2*len(systfiles)))

        #reset to zero 
        dt = plot_util.setzero_likelihood(dt)

        dt.SetName("nll_mu_gg_vs_mH_{0}".format(iter+iter2*len(systfiles)))
        fout.cd()
        dt.Write()
        signif.Write()
        maxsignif = plot_util.maxsignificance(signif)
        print "Maximum Significance: {0}".format(maxsignif )
        f.write("{0}\n".format(maxsignif ))
        signif_ave += maxsignif

f.write("Ave Significance: {0}".format(signif_ave / len(systfiles) ))
f.close()
