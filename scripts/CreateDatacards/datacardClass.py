#! /usr/bin/env python
import sys
import os
import re
import math
from scipy.special import erf
from ROOT import *
import ROOT
from array import array
from systematicsClass import *
from inputReader import *

## ------------------------------------
##  card and workspace class
## ------------------------------------

class datacardClass:

    def __init__(self):
    
        self.ID_4mu = 1
        self.ID_4e  = 2
        self.ID_2e2mu = 3    
        self.isFSR = True

    def loadIncludes(self):
        
        ROOT.gSystem.AddIncludePath("-I$ROOFITSYS/include/")
        ROOT.gSystem.AddIncludePath("-Iinclude/")
        ROOT.gROOT.ProcessLine(".L include/tdrstyle.cc")
        ROOT.gSystem.Load("libRooFit")
        ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
        ROOT.gSystem.Load("include/HiggsCSandWidth_cc.so")
        ROOT.gSystem.Load("include/HiggsCSandWidthSM4_cc.so")

    # cross section filter for 7 TeV efficiency
    def csFilter(self,hmass):

        a = 80.85
        b = 50.42
        
        f = 0.5 + 0.5*erf( (hmass - a)/b )
        
        return f

    # cs x br function 
    def makeXsBrFunction(self,signalProc,rrvMH):
            
        procName = "ggH"
        if(signalProc == 0): procName = "ggH" #dummy, when you sum up all the 5 chans
        if(signalProc == 1): procName = "ggH"
        if(signalProc == 2): procName = "qqH"
        if(signalProc == 3): procName = "WH"
        if(signalProc == 4): procName = "ZH"
        if(signalProc == 5): procName = "ttH"

        
        
        
        channelName = ""
        if (self.channel == self.ID_4mu): channelName = "4mu"
        elif (self.channel == self.ID_4e): channelName = "4e"
        elif (self.channel == self.ID_2e2mu): channelName = "2e2mu"
        else: print "Input Error: Unknown channel! (4mu = 1, 4e = 2, 2e2mu = 3)" 

     
        
        myCSWrhf = HiggsCSandWidth()
        
        histXsBr = ROOT.TH1F("hsmxsbr_{0}_{1}".format(procName,channelName),"", 1781, 109.75, 1000)
        
        for i in range(1,1782):
            
            mHVal = histXsBr.GetBinCenter(i)
            BR = 0.0 
            if (self.channel == self.ID_2e2mu):
                BR = myCSWrhf.HiggsBR(13,mHVal)
            else:
                BR = myCSWrhf.HiggsBR(12,mHVal)

            if (signalProc==0):
                totXs=0
                for ch in range(1,6):
                    totXs+=myCSWrhf.HiggsCS(ch, mHVal, self.sqrts)
                histXsBr.SetBinContent(i, totXs * BR)
            else:
                histXsBr.SetBinContent(i, myCSWrhf.HiggsCS(signalProc, mHVal, self.sqrts) * BR)

            #print '\nmakeXsBrFunction : procName=',procName,'   signalProc=',signalProc,'  mH (input)=',rrvMH.getVal(),
            #print '   CS=',myCSWrhf.HiggsCS(signalProc, mHVal, self.sqrts),'   BR=',BR
            
        rdhname = "rdhXsBr_{0}_{1}_{2}".format(procName,self.channel,self.sqrts)
        rdhXsBr = RooDataHist(rdhname,rdhname, ROOT.RooArgList(rrvMH), histXsBr)  
        
        return rdhXsBr
    
    # return trueVar if testStatement else return falseVar
    def getVariable(self,trueVar,falseVar,testStatement):

        if (testStatement): 
            return trueVar
        else:
            return falseVar
    
    # main datacard and workspace function
    def makeCardsWorkspaces(self, theMH, theis2D, theOutputDir, theInputs,theTemplateDir="templates2D", theIncludingError=False):

        ## --------------- SETTINGS AND DECLARATIONS --------------- ##
        DEBUG = False
        self.mH = theMH
        self.lumi = theInputs['lumi']
        self.sqrts = theInputs['sqrts']
        self.channel = theInputs['decayChannel']
        self.is2D = theis2D
        self.outputDir = theOutputDir
        self.sigMorph = theInputs['useCMS_zz4l_sigMELA']
        self.bkgMorph = theInputs['useCMS_zz4l_bkgMELA']
        self.templateDir = theTemplateDir
	self.bIncludingError=theIncludingError
        
        FactorizedShapes = False

        self.all_chan = theInputs['all']
        self.ggH_chan = theInputs['ggH']
        self.qqH_chan = theInputs['qqH']
        self.WH_chan = theInputs['WH']
        self.ZH_chan = theInputs['ZH']
        self.ttH_chan = theInputs['ttH']
        self.qqZZ_chan = theInputs['qqZZ']
        self.ggZZ_chan = theInputs['ggZZ']
        self.zjets_chan = theInputs['zjets']
        self.ttbar_chan = theInputs['ttbar']
        self.zbb_chan = theInputs['zbb']
        
        ## ---------------- SET PLOTTING STYLE ---------------- ## 
        ROOT.setTDRStyle(True)
        ROOT.gStyle.SetPalette(1)
        ROOT.gStyle.SetPadLeftMargin(0.16)        

        ## ---------------- VARIABLES FOR LATER --------------- ##
        self.bUseCBnoConvolution = False
        ForXSxBR = False

        myCSW = HiggsCSandWidth()
                
        ## ----------------- WIDTH AND RANGES ----------------- ##
        self.widthHVal =  myCSW.HiggsWidth(0,self.mH)
        if(self.widthHVal < 0.12):
            self.bUseCBnoConvolution = True
        self.isHighMass = False
        if self.mH >= 400:
            if theInputs['useHighMassReweightedShapes']:
                self.isHighMass = True
            else: print "useHighMassReweightedShapes set to FALSE, using non-reweighted shapes!"

            
        if(DEBUG): print "width: ",self.widthHVal
        
        self.windowVal = max( self.widthHVal, 1.0)
        lowside = 100.0
        if (self.mH >= 275):
            lowside = 180.0
        else:
            lowside = 100.0
        
        self.low_M = max( (self.mH - 20.*self.windowVal), lowside)
        self.high_M = min( (self.mH + 15.*self.windowVal), 1000)

        #self.low_M = 100.0
        #self.high_M = 800.0
       
        if (self.channel == self.ID_4mu): self.appendName = '4mu'
        elif (self.channel == self.ID_4e): self.appendName = '4e'
        elif (self.channel == self.ID_2e2mu): self.appendName = '2e2mu'
        else: print "Input Error: Unknown channel! (4mu = 1, 4e = 2, 2e2mu = 3)"

        self.isAltSig = False
        if (theInputs['doHypTest']):
            self.isAltSig = True
            
        if self.isAltSig and not self.all_chan :
            raise RuntimeError, "You asked to prepare DC and WS for Hyp Test but you did not want to sum over all signal channels. This is forbidden. Check inputs ! (it should have already send you this error message, strange that  you are here...)"

        if (self.isAltSig and not (self.is2D==1)):
            raise RunTimeError, "Cannot perform hypothesis testing without a 2D analysis, feature not supported yet. Exiting."
        

        self.appendHypType = theInputs['altHypLabel']
        if self.isAltSig and self.appendHypType=="" :
            self.appendHypType = "_ALT"
            
        
        ## ------------------------- SYSTEMATICS CLASSES ----------------------------- ##
    
        systematics = systematicsClass( self.mH, False, self.isFSR, theInputs)
        systematics_forXSxBR = systematicsClass( self.mH, True, self.isFSR,theInputs)

        ## -------------------------- SIGNAL SHAPE ----------------------------------- ##
    
        bins = 1000
        if(self.bUseCBnoConvolution): bins = 200
        if(self.bIncludingError): bins = 200

        CMS_zz4l_mass = ROOT.RooRealVar("CMS_zz4l_mass","CMS_zz4l_mass",self.low_M,self.high_M)
        CMS_zz4l_mass.setBins(bins,"fft") 

        self.LUMI = ROOT.RooRealVar("LUMI_{0:.0f}".format(self.sqrts),"LUMI_{0:.0f}".format(self.sqrts),self.lumi)
        self.LUMI.setConstant(True)
    
        self.MH = ROOT.RooRealVar("MH","MH",self.mH)
        self.MH.setConstant(True)
    
	# bIncludingError 
	CMS_zz4l_massErr = ROOT.RooRealVar("CMS_zz4l_massErr", "CMS_zz4l_massErr", 0.01*self.low_M/3., 0.01*self.high_M*5 );
	CMS_zz4l_massErr.setBins(100);

	# n2, alpha2 are right side parameters of DoubleCB
	# n, alpha are left side parameters of DoubleCB

        n_CB_d = 0.0
        alpha_CB_d = 0.0
        n2_CB_d = 0.0
        alpha2_CB_d = 0.0
        mean_CB_d = 0.0
        sigma_CB_d = 0.0
        mean_BW_d = self.mH
        gamma_BW_d = 0.0
        
        if(self.all_chan):
            rdhXsBrFuncV_1 = self.makeXsBrFunction(0,self.MH)
        else:
            rdhXsBrFuncV_1 = self.makeXsBrFunction(1,self.MH)
        rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}".format("ggH",self.channel,self.sqrts)
        rhfXsBrFuncV_1 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), rdhXsBrFuncV_1, 1)
        
        rdhXsBrFuncV_2 = self.makeXsBrFunction(2,self.MH)
        rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}".format("VBF",self.channel,self.sqrts)
        rhfXsBrFuncV_2 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), rdhXsBrFuncV_2, 1)
        
        rdhXsBrFuncV_3 = self.makeXsBrFunction(3,self.MH)
        rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}".format("WH",self.channel,self.sqrts)
        rhfXsBrFuncV_3 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), rdhXsBrFuncV_3, 1)
        
        rdhXsBrFuncV_4 = self.makeXsBrFunction(4,self.MH)
        rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}".format("ZH",self.channel,self.sqrts)
        rhfXsBrFuncV_4 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), rdhXsBrFuncV_4, 1)
        
        rdhXsBrFuncV_5 = self.makeXsBrFunction(5,self.MH)
        rhfname = "rhfXsBr_{0}_{1:.0f}_{2:.0f}".format("ttH",self.channel,self.sqrts)
        rhfXsBrFuncV_5 = ROOT.RooHistFunc(rhfname,rhfname, ROOT.RooArgSet(self.MH), rdhXsBrFuncV_5, 1)

    
        ## -------- Variable Definitions -------- ##
    
        name = "CMS_zz4l_mean_e_sig"
        CMS_zz4l_mean_e_sig = ROOT.RooRealVar(name,"CMS_zz4l_mean_e_sig",0.0,-10.0,10.0)
        name = "CMS_zz4l_sigma_e_sig"
        CMS_zz4l_sigma_e_sig = ROOT.RooRealVar(name,"CMS_zz4l_sigma_sig",3.0,0.0,30.0)
        
        name = "CMS_zz4l_mean_m_sig"
        CMS_zz4l_mean_m_sig = ROOT.RooRealVar(name,"CMS_zz4l_mean_sig",0.0,-10.0,10.0)
        name = "CMS_zz4l_sigma_m_sig"
        CMS_zz4l_sigma_m_sig = ROOT.RooRealVar(name,"CMS_zz4l_sigma_sig",3.0,0.0,30.0)
        
        name = "CMS_zz4l_alpha2_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_alpha2 = ROOT.RooRealVar(name,"CMS_zz4l_alpha2",1.,-10.,10.)
        name = "CMS_zz4l_n2_sig_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_n2 = ROOT.RooRealVar(name,"CMS_zz4l_n2",2.,-10.,10.)
        name = "CMS_zz4l_alpha_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_alpha = ROOT.RooRealVar(name,"CMS_zz4l_alpha",1.,-10.,10.)
        name = "CMS_zz4l_n_sig_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_n = ROOT.RooRealVar(name,"CMS_zz4l_n",2.,-10.,10.)
        name = "CMS_zz4l_mean_BW_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_mean_BW = ROOT.RooRealVar(name,"CMS_zz4l_mean_BW",self.mH,self.low_M,self.high_M)
        name = "CMS_zz4l_gamma_sig_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_gamma = ROOT.RooRealVar(name,"CMS_zz4l_gamma",10.,0.001,1000.)
        name = "CMS_zz4l_widthScale_{0}_{1:.0f}".format(self.channel,self.sqrts)
        CMS_zz4l_widthScale = ROOT.RooRealVar(name,"CMS_zz4l_widthScale",1.0)
        one = ROOT.RooRealVar("one","one",1.0)
        one.setConstant(True)
    
        CMS_zz4l_mean_BW.setVal( mean_BW_d )
        CMS_zz4l_gamma.setVal(0)
        CMS_zz4l_mean_e_sig.setVal(0)
        CMS_zz4l_sigma_e_sig.setVal(0)
        CMS_zz4l_mean_m_sig.setVal(0)
        CMS_zz4l_sigma_m_sig.setVal(0)
        CMS_zz4l_alpha.setVal(0)
        CMS_zz4l_n.setVal(0)
        CMS_zz4l_alpha2.setVal(0)
        CMS_zz4l_n2.setVal(0)
    
        CMS_zz4l_widthScale.setConstant(True)
        #CMS_zz4l_alpha.setConstant(True)  # also read from input file
        CMS_zz4l_mean_BW.setConstant(True)
        #CMS_zz4l_gamma_BW.setConstant(True)

        print "HEEERRRRRRRRRRRRRRRRREEEEEEE"

        print "mean_BW ", CMS_zz4l_mean_BW.getVal()
        print "gamma_BW ", CMS_zz4l_gamma.getVal()
        print "mean_e_sig ", CMS_zz4l_mean_e_sig.getVal()
        print "sigma_e ", CMS_zz4l_sigma_e_sig.getVal()
        print "mean_m_sig ",CMS_zz4l_mean_m_sig.getVal()
        print "sigma_m ", CMS_zz4l_sigma_m_sig.getVal()
        print "alpha ", CMS_zz4l_alpha.getVal()
        print "n ", CMS_zz4l_n.getVal()
        print "alpha2 ", CMS_zz4l_alpha2.getVal()
        print "n2 ", CMS_zz4l_n2.getVal()

                                                                


        ## -------------------- RooFormulaVar's -------------------- ##
        rfv_n_CB = ROOT.RooFormulaVar()
        rfv_alpha_CB = ROOT.RooFormulaVar()
        rfv_n2_CB = ROOT.RooFormulaVar()
        rfv_alpha2_CB = ROOT.RooFormulaVar()
        rfv_mean_CB = ROOT.RooFormulaVar()
        rfv_sigma_CB = ROOT.RooFormulaVar()
        
        name = "CMS_zz4l_n_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
        if self.isHighMass : rfv_n_CB = ROOT.RooFormulaVar(name,"("+theInputs['n_CB_shape_HM']+")"+"*(1+@1)",ROOT.RooArgList(self.MH,CMS_zz4l_n))
        else : rfv_n_CB = ROOT.RooFormulaVar(name,"("+theInputs['n_CB_shape']+")"+"*(1+@1)",ROOT.RooArgList(self.MH,CMS_zz4l_n))

        name = "CMS_zz4l_alpha_{0:.0f}_centralValue".format(self.channel)
        if self.isHighMass : rfv_alpha_CB = ROOT.RooFormulaVar(name,theInputs['alpha_CB_shape_HM'], ROOT.RooArgList(self.MH))
        else : rfv_alpha_CB = ROOT.RooFormulaVar(name,theInputs['alpha_CB_shape'], ROOT.RooArgList(self.MH))

        name = "CMS_zz4l_n2_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
        #if self.isHighMass : rfv_n2_CB = ROOT.RooFormulaVar(name,"("+theInputs['n2_CB_shape_HM']+")"+"*(1+@1)",ROOT.RooArgList(self.MH,CMS_zz4l_n2))
        #else : rfv_n2_CB = ROOT.RooFormulaVar(name,"("+theInputs['n2_CB_shape']+")"+"*(1+@1)",ROOT.RooArgList(self.MH,CMS_zz4l_n2))
        if self.isHighMass : rfv_n2_CB = ROOT.RooFormulaVar(name,"("+theInputs['n2_CB_shape_HM']+")",ROOT.RooArgList(self.MH))
        else : rfv_n2_CB = ROOT.RooFormulaVar(name,"("+theInputs['n2_CB_shape']+")",ROOT.RooArgList(self.MH))

        name = "CMS_zz4l_alpha2_{0:.0f}_centralValue".format(self.channel)
        if self.isHighMass : rfv_alpha2_CB = ROOT.RooFormulaVar(name,theInputs['alpha2_CB_shape_HM'], ROOT.RooArgList(self.MH))
        else : rfv_alpha2_CB = ROOT.RooFormulaVar(name,theInputs['alpha2_CB_shape'], ROOT.RooArgList(self.MH))

        name = "CMS_zz4l_mean_sig_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
        if (self.channel == self.ID_4mu) :
            if self.isHighMass : rfv_mean_CB = ROOT.RooFormulaVar(name,"("+theInputs['mean_CB_shape_HM']+")"+"+@0*@1", ROOT.RooArgList(self.MH, CMS_zz4l_mean_m_sig))
            else : rfv_mean_CB = ROOT.RooFormulaVar(name,"("+theInputs['mean_CB_shape']+")"+"+@0*@1", ROOT.RooArgList(self.MH, CMS_zz4l_mean_m_sig))
        elif (self.channel == self.ID_4e) :
            if self.isHighMass : rfv_mean_CB = ROOT.RooFormulaVar(name,"("+theInputs['mean_CB_shape_HM']+")"+"+@0*@1", ROOT.RooArgList(self.MH, CMS_zz4l_mean_e_sig))
            else : rfv_mean_CB = ROOT.RooFormulaVar(name,"("+theInputs['mean_CB_shape']+")"+"+@0*@1", ROOT.RooArgList(self.MH, CMS_zz4l_mean_e_sig))
        elif (self.channel == self.ID_2e2mu) :
            if self.isHighMass : rfv_mean_CB = ROOT.RooFormulaVar(name,"("+theInputs['mean_CB_shape_HM']+")"+"+ @0*@1 + @0*@2", ROOT.RooArgList(self.MH, CMS_zz4l_mean_m_sig,CMS_zz4l_mean_e_sig))
            else : rfv_mean_CB = ROOT.RooFormulaVar(name,"("+theInputs['mean_CB_shape']+")"+"+ @0*@1 + @0*@2", ROOT.RooArgList(self.MH, CMS_zz4l_mean_m_sig,CMS_zz4l_mean_e_sig))
        

        name = "CMS_zz4l_sigma_sig_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
        if (self.channel == self.ID_4mu) :
            if self.isHighMass : rfv_sigma_CB = ROOT.RooFormulaVar(name,"("+theInputs['sigma_CB_shape_HM']+")"+"*(1+@1)", ROOT.RooArgList(self.MH, CMS_zz4l_sigma_m_sig))
            else : rfv_sigma_CB = ROOT.RooFormulaVar(name,"("+theInputs['sigma_CB_shape']+")"+"*(1+@1)", ROOT.RooArgList(self.MH, CMS_zz4l_sigma_m_sig))
        elif (self.channel == self.ID_4e) :
            if self.isHighMass : rfv_sigma_CB = ROOT.RooFormulaVar(name,"("+theInputs['sigma_CB_shape_HM']+")"+"*(1+@1)", ROOT.RooArgList(self.MH, CMS_zz4l_sigma_e_sig))
            else : rfv_sigma_CB = ROOT.RooFormulaVar(name,"("+theInputs['sigma_CB_shape']+")"+"*(1+@1)", ROOT.RooArgList(self.MH, CMS_zz4l_sigma_e_sig))
        elif (self.channel == self.ID_2e2mu) :
            if self.isHighMass : rfv_sigma_CB = ROOT.RooFormulaVar(name,"("+theInputs['sigma_CB_shape_HM']+")"+"*(1+@1)*(1+@2)", ROOT.RooArgList(self.MH, CMS_zz4l_sigma_m_sig,CMS_zz4l_sigma_e_sig))
            else : rfv_sigma_CB = ROOT.RooFormulaVar(name,"("+theInputs['sigma_CB_shape']+")"+"*(1+@1)*(1+@2)", ROOT.RooArgList(self.MH, CMS_zz4l_sigma_m_sig,CMS_zz4l_sigma_e_sig))


        name = "CMS_zz4l_gamma_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
        rfv_gamma_BW = ROOT.RooFormulaVar(name,"("+theInputs['gamma_BW_shape_HM']+")"+"*(1+@1)",ROOT.RooArgList(self.MH,CMS_zz4l_gamma))

        if (DEBUG): print " DEBUG *********  ", theInputs['sigma_CB_shape'] 

        print "n_CB ", rfv_n_CB.getVal()
        print "alpha_CB ", rfv_alpha_CB.getVal()
        print "n2_CB ", rfv_n2_CB.getVal()
        print "alpha2_CB ", rfv_alpha2_CB.getVal()
        print "mean_CB ", rfv_mean_CB.getVal()
        print "sigma_CB ", rfv_sigma_CB.getVal()
        print "gamma_BW ", rfv_gamma_BW.getVal()    

        
        CMS_zz4l_mean_sig_NoConv = ROOT.RooFormulaVar("CMS_zz4l_mean_sig_NoConv_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts),"@0+@1", ROOT.RooArgList(rfv_mean_CB, self.MH))

        print "mean_sig_NoConv ", CMS_zz4l_mean_sig_NoConv.getVal()
        
        ## --------------------- SHAPE FUNCTIONS ---------------------- ##
    
        signalCB_ggH = ROOT.RooDoubleCB("signalCB_ggH","signalCB_ggH",CMS_zz4l_mass, self.getVariable(CMS_zz4l_mean_sig_NoConv,rfv_mean_CB, self.bUseCBnoConvolution) , self.getVariable(CMS_zz4l_massErr,rfv_sigma_CB, self.bIncludingError),rfv_alpha_CB,rfv_n_CB, rfv_alpha2_CB, rfv_n2_CB)
        #Low mass pdf
        signalBW_ggH = ROOT.RooRelBWUFParam("signalBW_ggH", "signalBW_ggH",CMS_zz4l_mass,CMS_zz4l_mean_BW,CMS_zz4l_widthScale)
        sig_ggH =  ROOT.RooFFTConvPdf("sig_ggH","BW (X) CB",CMS_zz4l_mass,signalBW_ggH,signalCB_ggH, 2)
        #High mass pdf
        signalBW_ggH_HM = ROOT.RooRelBWHighMass("signalBW_ggH", "signalBW_ggH",CMS_zz4l_mass,CMS_zz4l_mean_BW,rfv_gamma_BW)
        sig_ggH_HM =  ROOT.RooFFTConvPdf("sig_ggH","BW (X) CB",CMS_zz4l_mass,signalBW_ggH_HM,signalCB_ggH, 2)
  
        
        signalCB_VBF = ROOT.RooDoubleCB("signalCB_VBF","signalCB_VBF",CMS_zz4l_mass,self.getVariable(CMS_zz4l_mean_sig_NoConv,rfv_mean_CB,self.bUseCBnoConvolution),self.getVariable(CMS_zz4l_massErr,rfv_sigma_CB, self.bIncludingError),rfv_alpha_CB,rfv_n_CB, rfv_alpha2_CB, rfv_n2_CB)
        #Low mass pdf
        signalBW_VBF = ROOT.RooRelBWUFParam("signalBW_VBF", "signalBW_VBF",CMS_zz4l_mass,CMS_zz4l_mean_BW,CMS_zz4l_widthScale)
        sig_VBF = ROOT.RooFFTConvPdf("sig_VBF","BW (X) CB",CMS_zz4l_mass,signalBW_VBF,signalCB_VBF, 2)
        #High mass pdf
        signalBW_VBF_HM = ROOT.RooRelBWHighMass("signalBW_VBF", "signalBW_VBF",CMS_zz4l_mass,CMS_zz4l_mean_BW,rfv_gamma_BW)
        sig_VBF_HM = ROOT.RooFFTConvPdf("sig_VBF","BW (X) CB",CMS_zz4l_mass,signalBW_VBF_HM,signalCB_VBF, 2)
                       
        
        signalCB_WH = ROOT.RooDoubleCB("signalCB_WH","signalCB_WH",CMS_zz4l_mass,self.getVariable(CMS_zz4l_mean_sig_NoConv,rfv_mean_CB,self.bUseCBnoConvolution),self.getVariable(CMS_zz4l_massErr,rfv_sigma_CB, self.bIncludingError),rfv_alpha_CB,rfv_n_CB, rfv_alpha2_CB, rfv_n2_CB)
        #Low mass pdf
        signalBW_WH = ROOT.RooRelBWUFParam("signalBW_WH", "signalBW_WH",CMS_zz4l_mass,CMS_zz4l_mean_BW,CMS_zz4l_widthScale)
        sig_WH = ROOT.RooFFTConvPdf("sig_WH","BW (X) CB",CMS_zz4l_mass,signalBW_WH,signalCB_WH, 2)
        #High mass pdf
        signalBW_WH_HM = ROOT.RooRelBWHighMass("signalBW_WH", "signalBW_WH",CMS_zz4l_mass,CMS_zz4l_mean_BW,rfv_gamma_BW)
        sig_WH_HM = ROOT.RooFFTConvPdf("sig_WH","BW (X) CB",CMS_zz4l_mass,signalBW_WH_HM,signalCB_WH, 2)

        
        signalCB_ZH = ROOT.RooDoubleCB("signalCB_ZH","signalCB_ZH",CMS_zz4l_mass,self.getVariable(CMS_zz4l_mean_sig_NoConv,rfv_mean_CB,self.bUseCBnoConvolution),self.getVariable(CMS_zz4l_massErr,rfv_sigma_CB, self.bIncludingError),rfv_alpha_CB,rfv_n_CB, rfv_alpha2_CB, rfv_n2_CB)
        #Low mass pdf
        signalBW_ZH = ROOT.RooRelBWUFParam("signalBW_ZH", "signalBW_ZH",CMS_zz4l_mass,CMS_zz4l_mean_BW,CMS_zz4l_widthScale)
        sig_ZH = ROOT.RooFFTConvPdf("sig_ZH","BW (X) CB",CMS_zz4l_mass,signalBW_ZH,signalCB_ZH, 2)
        #High mass pdf
        signalBW_ZH_HM = ROOT.RooRelBWHighMass("signalBW_ZH", "signalBW_ZH",CMS_zz4l_mass,CMS_zz4l_mean_BW,rfv_gamma_BW)
        sig_ZH_HM = ROOT.RooFFTConvPdf("sig_ZH","BW (X) CB",CMS_zz4l_mass,signalBW_ZH_HM,signalCB_ZH, 2)

        
        signalCB_ttH = ROOT.RooDoubleCB("signalCB_ttH","signalCB_ttH",CMS_zz4l_mass,self.getVariable(CMS_zz4l_mean_sig_NoConv,rfv_mean_CB,self.bUseCBnoConvolution),self.getVariable(CMS_zz4l_massErr,rfv_sigma_CB, self.bIncludingError),rfv_alpha_CB,rfv_n_CB, rfv_alpha2_CB, rfv_n2_CB)
        #Low mass pdf
        signalBW_ttH = ROOT.RooRelBWUFParam("signalBW_ttH", "signalBW_ttH",CMS_zz4l_mass,CMS_zz4l_mean_BW,CMS_zz4l_widthScale)
        sig_ttH = ROOT.RooFFTConvPdf("sig_ttH","BW (X) CB",CMS_zz4l_mass,signalBW_ttH,signalCB_ttH, 2) 
        #High mass pdf
        signalBW_ttH_HM = ROOT.RooRelBWHighMass("signalBW_ttH", "signalBW_ttH",CMS_zz4l_mass,CMS_zz4l_mean_BW,rfv_gamma_BW)
        sig_ttH_HM = ROOT.RooFFTConvPdf("sig_ttH","BW (X) CB",CMS_zz4l_mass,signalBW_ttH_HM,signalCB_ttH, 2)
        
        
        ## Buffer fraction for cyclical behavior
        sig_ggH.setBufferFraction(0.2)
        sig_VBF.setBufferFraction(0.2)
        sig_WH.setBufferFraction(0.2)
        sig_ZH.setBufferFraction(0.2)
        sig_ttH.setBufferFraction(0.2)
        
        sig_ggH_HM.setBufferFraction(0.2)
        sig_VBF_HM.setBufferFraction(0.2)
        sig_WH_HM.setBufferFraction(0.2)
        sig_ZH_HM.setBufferFraction(0.2)
        sig_ttH_HM.setBufferFraction(0.2)


	#------------------------------------------------begin  bIncludingError 
	#(g_channel == ID_2e2mu):
	# sprintf( name, "CMS_zz4l_sigma_sig_%i_centralValue", g_channel );
	# rfv_sigma_CB = new RooFormulaVar(name, "(4.2e-10*@0*@0*@0*@0 - 2.5e-07*@0*@0*@0 + 3.2e-05*@0*@0 + 0.013*@0)*(1+@1)", RooArgList(mH, CMS_zz4l_sigma_sig));
	#	sprintf( name, "CMS_zz4l_sigmaB_mean_%i_centralValue", g_channel ); // for zz bkg, the average error is larger than signal by 10%, so multiply by a factor of 1.1
	#	rfv_sigmaB_mean =  new RooFormulaVar(name, "(4.2e-10*@0*@0*@0*@0 - 2.5e-07*@0*@0*@0 + 3.2e-05*@0*@0 + 0.013*@0)*1.1*(1+@1)", RooArgList(CMS_zz4l_mass, CMS_zz4l_sigma_sig));

        name = "CMS_zz4l_sigmaB_sig_{0:.0f}_{1:.0f}_centralValue".format(self.channel,self.sqrts)
	rfv_sigmaB_mean = ROOT.RooFormulaVar(name,"("+theInputs['sigma_CB_shape']+")"+"*(1+@1)", ROOT.RooArgList(CMS_zz4l_mass, CMS_zz4l_sigma_m_sig))

	name = "CMS_zz4l_massErrS_kappa_{0:.0f}".format(self.channel);
	rfv_sigma_kappa = ROOT.RooFormulaVar(name, "@0*0 + 1.4", ROOT.RooArgList(self.MH)); #the kappa should be parametrized as a function of MH  --> TBD
	pdfErrS = ROOT.RooLognormal("pdfErrS", "pdfErrS", CMS_zz4l_massErr,  rfv_sigma_CB, rfv_sigma_kappa);
	name = "CMS_zz4l_massErrB_kappa_{0:.0f}".format(self.channel);
	rfv_sigmaB_kappa = ROOT.RooFormulaVar(name, "@0*0 + 1.4", ROOT.RooArgList(CMS_zz4l_mass));  #for bkg,  the kappa should be parametrized as a function of m4l --> TBD
	pdfErrB = ROOT.RooLognormal("pdfErrB", "pdfErrB", CMS_zz4l_massErr,  rfv_sigmaB_mean, rfv_sigmaB_kappa);

	#sig_ggHErr = new RooProdPdf("sig_ggHErr","BW (X) CB * pdfErr", pdfErrS, Conditional(bUseCBnoConvolution?(signalCB):(*sig_ggH), CMS_zz4l_mass));
	sig_ggHErr = ROOT.RooProdPdf("sig_ggHErr","BW (X) CB * pdfErr", ROOT.RooArgSet(pdfErrS), ROOT.RooFit.Conditional(ROOT.RooArgSet(signalCB_ggH), ROOT.RooArgSet(CMS_zz4l_mass)));
	sig_VBFErr = ROOT.RooProdPdf("sig_VBFErr","BW (X) CB * pdfErr", ROOT.RooArgSet(pdfErrS), ROOT.RooFit.Conditional(ROOT.RooArgSet(signalCB_VBF), ROOT.RooArgSet(CMS_zz4l_mass)));
	sig_WHErr = ROOT.RooProdPdf("sig_WHErr","BW (X) CB * pdfErr", ROOT.RooArgSet(pdfErrS), ROOT.RooFit.Conditional(ROOT.RooArgSet(signalCB_WH), ROOT.RooArgSet(CMS_zz4l_mass)));
	sig_ZHErr = ROOT.RooProdPdf("sig_ZHErr","BW (X) CB * pdfErr", ROOT.RooArgSet(pdfErrS), ROOT.RooFit.Conditional(ROOT.RooArgSet(signalCB_ZH), ROOT.RooArgSet(CMS_zz4l_mass)));
	sig_ttHErr = ROOT.RooProdPdf("sig_ttHErr","BW (X) CB * pdfErr", ROOT.RooArgSet(pdfErrS), ROOT.RooFit.Conditional(ROOT.RooArgSet(signalCB_ttH), ROOT.RooArgSet(CMS_zz4l_mass)));

	#------------------------------------------------end bIncludingError 

        ## --------------------------- MELA 2D PDFS ------------------------- ##

        discVarName = "melaLD"
        D = ROOT.RooRealVar(discVarName,discVarName,0,1)
    
        templateSigName = "{0}/Dsignal_{1}.root".format(self.templateDir ,self.appendName)
        
        sigTempFile = ROOT.TFile(templateSigName)
        sigTemplate = sigTempFile.Get("h_mzzD")
        sigTemplate_Up = sigTempFile.Get("h_mzzD_up")
        sigTemplate_Down = sigTempFile.Get("h_mzzD_dn")
        
        TemplateName = "sigTempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,D),sigTemplate)
        TemplateName = "sigTempDataHist_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTempDataHist_Up = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,D),sigTemplate_Up)
        TemplateName = "sigTempDataHist_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTempDataHist_Down = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,D),sigTemplate_Down)

        
        TemplateName = "sigTemplatePdf_ggH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplatePdf_ggH = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist)
        TemplateName = "sigTemplatePdf_ggH_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplatePdf_ggH_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Up)
        TemplateName = "sigTemplatePdf_ggH_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplatePdf_ggH_Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Down)
        
        TemplateName = "sigTemplatePdf_VBF_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplatePdf_VBF = ROOT.RooHistPdf(TemplateName,TemplateName,RooArgSet(CMS_zz4l_mass,D),sigTempDataHist)
        TemplateName = "sigTemplatePdf_VBF_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplatePdf_VBF_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Up)
        TemplateName = "sigTemplatePdf_VBF_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplatePdf_VBF_Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Down)
        
        TemplateName = "sigTemplatePdf_WH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplatePdf_WH = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist)
        TemplateName = "sigTemplatePdf_WH_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplatePdf_WH_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Up)
        TemplateName = "sigTemplatePdf_WH_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplatePdf_WH_Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Down)
        
        TemplateName = "sigTemplatePdf_ZH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplatePdf_ZH = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist)
        TemplateName = "sigTemplatePdf_ZH_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplatePdf_ZH_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Up)
        TemplateName = "sigTemplatePdf_ZH_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplatePdf_ZH_Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Down)
        
        TemplateName = "sigTemplatePdf_ZH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplatePdf_ttH = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist)
        TemplateName = "sigTemplatePdf_ZH_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplatePdf_ttH_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Up)
        TemplateName = "sigTemplatePdf_ZH_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplatePdf_ttH_Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_Down)

        funcList_ggH = ROOT.RooArgList()  
        funcList_VBF = ROOT.RooArgList()
        funcList_WH  = ROOT.RooArgList()
        funcList_ZH  = ROOT.RooArgList()
        funcList_ttH = ROOT.RooArgList()


        if(self.isAltSig):
            #only ggH because if we do hypothesis testing we sum up over the channels in any case
              templateSigName = "{0}/Dsignal{2}_{1}.root".format(self.templateDir,self.appendName, self.appendHypType)
              print 'Taking 2D template for ALT signal from ',templateSigName
              sigTempFile = ROOT.TFile(templateSigName)
              sigTemplate = sigTempFile.Get("h_mzzD")
              sigTemplate_Up = sigTempFile.Get("h_mzzD_up")
              sigTemplate_Down = sigTempFile.Get("h_mzzD_dn")
              
              TemplateName = "sigTempDataHist_{0:.0f}_{1:.0f}{2}".format(self.channel,self.sqrts, self.appendHypType)
              sigTempDataHist_ALT = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,D),sigTemplate)
              TemplateName = "sigTempDataHist_Up_{0:.0f}_{1:.0f}{2}".format(self.channel,self.sqrts, self.appendHypType)
              sigTempDataHist_ALT_Up = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,D),sigTemplate_Up)
              TemplateName = "sigTempDataHist_Down_{0:.0f}_{1:.0f}{2}".format(self.channel,self.sqrts, self.appendHypType)
              sigTempDataHist_ALT_Down = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,D),sigTemplate_Down)
              TemplateName = "sigTemplatePdf_ggH{2}_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts, self.appendHypType)
              sigTemplatePdf_ggH_ALT = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_ALT)
              TemplateName = "sigTemplatePdf_ggH{2}_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts, self.appendHypType)
              sigTemplatePdf_ggH_ALT_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_ALT_Up)
              TemplateName = "sigTemplatePdf_ggH{2}_Down_{0:.0f}_{1:.0f}{2}".format(self.channel,self.sqrts, self.appendHypType)
              sigTemplatePdf_ggH_ALT_Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),sigTempDataHist_ALT_Down)

        funcList_ggH_ALT = ROOT.RooArgList() 

        if(self.sigMorph):
            
            funcList_ggH.add(sigTemplatePdf_ggH)
            funcList_ggH.add(sigTemplatePdf_ggH_Up)
            funcList_ggH.add(sigTemplatePdf_ggH_Down)  
            
            funcList_VBF.add(sigTemplatePdf_VBF)
            funcList_VBF.add(sigTemplatePdf_VBF_Up)
            funcList_VBF.add(sigTemplatePdf_VBF_Down)  
            
            funcList_WH.add(sigTemplatePdf_WH)
            funcList_WH.add(sigTemplatePdf_WH_Up)
            funcList_WH.add(sigTemplatePdf_WH_Down)  
            
            funcList_ZH.add(sigTemplatePdf_ZH)
            funcList_ZH.add(sigTemplatePdf_ZH_Up)
            funcList_ZH.add(sigTemplatePdf_ZH_Down)  
            
            funcList_ttH.add(sigTemplatePdf_ttH)
            funcList_ttH.add(sigTemplatePdf_ttH_Up)
            funcList_ttH.add(sigTemplatePdf_ttH_Down)  
            if(self.isAltSig):
                funcList_ggH_ALT.add(sigTemplatePdf_ggH_ALT)
                funcList_ggH_ALT.add(sigTemplatePdf_ggH_ALT_Up)
                funcList_ggH_ALT.add(sigTemplatePdf_ggH_ALT_Down)
        else:
            
            funcList_ggH.add(sigTemplatePdf_ggH)
            funcList_VBF.add(sigTemplatePdf_VBF)
            funcList_WH.add(sigTemplatePdf_WH)
            funcList_ZH.add(sigTemplatePdf_ZH)
            funcList_ttH.add(sigTemplatePdf_ttH)
            if(self.isAltSig):
                funcList_ggH_ALT.add(sigTemplatePdf_ggH_ALT)
    
        morphSigVarName = "CMS_zz4l_sigMELA_{0:.0f}".format(self.channel)
        alphaMorphSig = ROOT.RooRealVar(morphSigVarName,morphSigVarName,0,-20,20)
        if(self.sigMorph): alphaMorphSig.setConstant(False)
        else: alphaMorphSig.setConstant(True)
        
        morphVarListSig = ROOT.RooArgList()
    
        if(self.sigMorph): morphVarListSig.add(alphaMorphSig)  ## just one morphing for all signal processes (fully correlated systematics)
        
        TemplateName = "sigTemplateMorphPdf_ggH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplateMorphPdf_ggH = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,D,true,funcList_ggH,morphVarListSig,1.0,1)
        
        TemplateName = "sigTemplateMorphPdf_VBF_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplateMorphPdf_VBF = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,D,true,funcList_VBF,morphVarListSig,1.0,1)
        
        TemplateName = "sigTemplateMorphPdf_WH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplateMorphPdf_WH = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,D,true,funcList_WH,morphVarListSig,1.0,1)
        
        TemplateName = "sigTemplateMorphPdf_ZH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplateMorphPdf_ZH = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,D,true,funcList_ZH,morphVarListSig,1.0,1)
        
        TemplateName = "sigTemplateMorphPdf_ttH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplateMorphPdf_ttH = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,D,true,funcList_ttH,morphVarListSig,1.0,1)
        if(self.isAltSig):
            TemplateName = "sigTemplateMorphPdf_ggH{2}_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts,self.appendHypType)
            sigTemplateMorphPdf_ggH_ALT = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,D,true,funcList_ggH_ALT,morphVarListSig,1.0,1)
    
        sig2d_ggH = ROOT.RooProdPdf("sig2d_ggH","sig2d_ggH",ROOT.RooArgSet(self.getVariable(sig_ggH_HM,sig_ggH,self.isHighMass)),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_ggH),ROOT.RooArgSet(D)))
        sig2d_VBF = ROOT.RooProdPdf("sig2d_VBF","sig2d_VBF",ROOT.RooArgSet(self.getVariable(sig_VBF_HM,sig_VBF,self.isHighMass)),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_VBF),ROOT.RooArgSet(D)))
        sig2d_WH = ROOT.RooProdPdf("sig2d_WH","sig2d_WH",ROOT.RooArgSet(self.getVariable(sig_WH_HM,sig_WH,self.isHighMass)),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_WH),ROOT.RooArgSet(D)))
        sig2d_ZH = ROOT.RooProdPdf("sig2d_ZH","sig2d_ZH",ROOT.RooArgSet(self.getVariable(sig_ZH_HM,sig_ZH,self.isHighMass)),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_ZH),ROOT.RooArgSet(D)))
        sig2d_ttH = ROOT.RooProdPdf("sig2d_ttH","sig2d_ttH",ROOT.RooArgSet(self.getVariable(sig_ttH_HM,sig_ttH,self.isHighMass)),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_ttH),ROOT.RooArgSet(D)))
                
        sigCB2d_ggH = ROOT.RooProdPdf("sigCB2d_ggH","sigCB2d_ggH",ROOT.RooArgSet(self.getVariable(sig_ggHErr,signalCB_ggH,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_ggH),ROOT.RooArgSet(D)))
        sigCB2d_VBF = ROOT.RooProdPdf("sigCB2d_VBF","sigCB2d_VBF",ROOT.RooArgSet(self.getVariable(sig_VBFErr,signalCB_VBF,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_VBF),ROOT.RooArgSet(D)))
        sigCB2d_WH = ROOT.RooProdPdf("sigCB2d_WH","sigCB2d_WH",ROOT.RooArgSet(self.getVariable(sig_WHErr,signalCB_WH,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_WH),ROOT.RooArgSet(D)))
        sigCB2d_ZH = ROOT.RooProdPdf("sigCB2d_ZH","sigCB2d_ZH",ROOT.RooArgSet(self.getVariable(sig_ZHErr,signalCB_ZH,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_ZH),ROOT.RooArgSet(D)))
        sigCB2d_ttH = ROOT.RooProdPdf("sigCB2d_ttH","sigCB2d_ttH",ROOT.RooArgSet(self.getVariable(sig_ttHErr,signalCB_ttH,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_ttH),ROOT.RooArgSet(D)))
        
        if(self.isAltSig):
            sig2d_ggH_ALT = ROOT.RooProdPdf("sig2d_ggH{0}".format(self.appendHypType),"sig2d_ggH{0}".format(self.appendHypType),ROOT.RooArgSet(sig_ggH),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_ggH_ALT),ROOT.RooArgSet(D)))
            sigCB2d_ggH_ALT = ROOT.RooProdPdf("sigCB2d_ggH{0}".format(self.appendHypType),"sigCB2d_ggH{0}".format(self.appendHypType),ROOT.RooArgSet(signalCB_ggH),ROOT.RooFit.Conditional(ROOT.RooArgSet(sigTemplateMorphPdf_ggH_ALT),ROOT.RooArgSet(D)))            



        ## --------------------------- superMELA 1D PDFS ------------------------- ##

        superDiscVarName = "supermelaLD"
        SD = ROOT.RooRealVar(superDiscVarName,superDiscVarName,0,1)
    
        templateSDSigName = "{0}/Dsignal_superMELA_{1}.root".format(self.templateDir ,self.appendName)
        sigTempSDFile = ROOT.TFile(templateSDSigName)
        sigTemplateSD = sigTempSDFile.Get("hSuperD_sig") 
        
        TemplateSDName = "sigTempSDDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTempSDDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(SD),sigTemplateSD)
        
        TemplateSDName = "sigTemplateSDPdf_ggH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplateSDPdf_ggH = ROOT.RooHistPdf(TemplateSDName,TemplateSDName,ROOT.RooArgSet(SD),sigTempSDDataHist)
        
        TemplateSDName = "sigTemplateSDPdf_VBF_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplateSDPdf_VBF = ROOT.RooHistPdf(TemplateSDName,TemplateSDName,RooArgSet(SD),sigTempSDDataHist)
        
        TemplateSDName = "sigTemplateSDPdf_WH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplateSDPdf_WH = ROOT.RooHistPdf(TemplateSDName,TemplateSDName,ROOT.RooArgSet(SD),sigTempSDDataHist)
        
        TemplateSDName = "sigTemplateSDPdf_ZH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplateSDPdf_ZH = ROOT.RooHistPdf(TemplateSDName,TemplateSDName,ROOT.RooArgSet(SD),sigTempSDDataHist)
        
        TemplateSDName = "sigTemplateSDPdf_ZH_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        sigTemplateSDPdf_ttH = ROOT.RooHistPdf(TemplateSDName,TemplateSDName,ROOT.RooArgSet(SD),sigTempSDDataHist)
        print sigTemplateSDPdf_ttH

        ##--------------##

        ## -------------------------- BACKGROUND SHAPES ---------------------------------- ##
    
        ## qqZZ contribution
        name = "CMS_qqzzbkg_a0_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a0 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a0",115.3,0.,200.)
        name = "CMS_qqzzbkg_a1_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a1 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a1",21.96,0.,200.)
        name = "CMS_qqzzbkg_a2_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a2 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a2",122.8,0.,200.)
        name = "CMS_qqzzbkg_a3_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a3 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a3",0.03479,0.,1.)
        name = "CMS_qqzzbkg_a4_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a4 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a4",185.5,0.,200.)
        name = "CMS_qqzzbkg_a5_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a5 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a5",12.67,0.,200.)
        name = "CMS_qqzzbkg_a6_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a6 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a6",34.81,0.,100.)
        name = "CMS_qqzzbkg_a7_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a7 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a7",0.1393,0.,1.)
        name = "CMS_qqzzbkg_a8_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a8 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a8",66.,0.,200.)
        name = "CMS_qqzzbkg_a9_{0:.0f}_{1:.0f}".format( self.channel,self.sqrts )
        CMS_qqzzbkg_a9 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a9",0.07191,0.,1.)
        name = "CMS_qqzzbkg_a10_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts )
        CMS_qqzzbkg_a10 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a10",94.11,0.,200.)
        name = "CMS_qqzzbkg_a11_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts )
        CMS_qqzzbkg_a11 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a11",-5.111,-100.,100.)
        name = "CMS_qqzzbkg_a12_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts )
        CMS_qqzzbkg_a12 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a12",4834,0.,10000.)
        name = "CMS_qqzzbkg_a13_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts )
        CMS_qqzzbkg_a13 = ROOT.RooRealVar(name,"CMS_qqzzbkg_a13",0.2543,0.,1.)
        

        if (DEBUG) :
            print "qqZZshape_a0 = ",theInputs['qqZZshape_a0']
            print "qqZZshape_a1 = ",theInputs['qqZZshape_a1']
            print "qqZZshape_a2 = ",theInputs['qqZZshape_a2']
            print "qqZZshape_a3 = ",theInputs['qqZZshape_a3']
            print "qqZZshape_a4 = ",theInputs['qqZZshape_a4']
            print "qqZZshape_a5 = ",theInputs['qqZZshape_a5']
            print "qqZZshape_a6 = ",theInputs['qqZZshape_a6']
            print "qqZZshape_a7 = ",theInputs['qqZZshape_a7']
            print "qqZZshape_a8 = ",theInputs['qqZZshape_a8']
            print "qqZZshape_a9 = ",theInputs['qqZZshape_a9']
            print "qqZZshape_a10 = ",theInputs['qqZZshape_a10']
            print "qqZZshape_a11 = ",theInputs['qqZZshape_a11']
            print "qqZZshape_a12 = ",theInputs['qqZZshape_a12']
            print "qqZZshape_a13 = ",theInputs['qqZZshape_a13']

        
        CMS_qqzzbkg_a0.setVal(theInputs['qqZZshape_a0'])
        CMS_qqzzbkg_a1.setVal(theInputs['qqZZshape_a1'])
        CMS_qqzzbkg_a2.setVal(theInputs['qqZZshape_a2'])
        CMS_qqzzbkg_a3.setVal(theInputs['qqZZshape_a3'])
        CMS_qqzzbkg_a4.setVal(theInputs['qqZZshape_a4'])
        CMS_qqzzbkg_a5.setVal(theInputs['qqZZshape_a5'])
        CMS_qqzzbkg_a6.setVal(theInputs['qqZZshape_a6'])
        CMS_qqzzbkg_a7.setVal(theInputs['qqZZshape_a7'])
        CMS_qqzzbkg_a8.setVal(theInputs['qqZZshape_a8'])
        CMS_qqzzbkg_a9.setVal(theInputs['qqZZshape_a9'])
        CMS_qqzzbkg_a10.setVal(theInputs['qqZZshape_a10'])
        CMS_qqzzbkg_a11.setVal(theInputs['qqZZshape_a11'])
        CMS_qqzzbkg_a12.setVal(theInputs['qqZZshape_a12'])
        CMS_qqzzbkg_a13.setVal(theInputs['qqZZshape_a13'])
        
        CMS_qqzzbkg_a0.setConstant(True)
        CMS_qqzzbkg_a1.setConstant(True)
        CMS_qqzzbkg_a2.setConstant(True)
        CMS_qqzzbkg_a3.setConstant(True)
        CMS_qqzzbkg_a4.setConstant(True)
        CMS_qqzzbkg_a5.setConstant(True)
        CMS_qqzzbkg_a6.setConstant(True)
        CMS_qqzzbkg_a7.setConstant(True)
        CMS_qqzzbkg_a8.setConstant(True)
        CMS_qqzzbkg_a9.setConstant(True)
        CMS_qqzzbkg_a10.setConstant(True)
        CMS_qqzzbkg_a11.setConstant(True)
        CMS_qqzzbkg_a12.setConstant(True)
        CMS_qqzzbkg_a13.setConstant(True)
        
        bkg_qqzz = ROOT.RooqqZZPdf_v2("bkg_qqzzTmp","bkg_qqzzTmp",CMS_zz4l_mass,CMS_qqzzbkg_a0,CMS_qqzzbkg_a1,CMS_qqzzbkg_a2,CMS_qqzzbkg_a3,CMS_qqzzbkg_a4,CMS_qqzzbkg_a5,CMS_qqzzbkg_a6,CMS_qqzzbkg_a7,CMS_qqzzbkg_a8,CMS_qqzzbkg_a9,CMS_qqzzbkg_a10,CMS_qqzzbkg_a11,CMS_qqzzbkg_a12,CMS_qqzzbkg_a13)
        
        ## ggZZ contribution
        name = "CMS_ggzzbkg_a0_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        CMS_ggzzbkg_a0 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a0",115.3,0.,200.)
        name = "CMS_ggzzbkg_a1_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        CMS_ggzzbkg_a1 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a1",21.96,0.,200.)
        name = "CMS_ggzzbkg_a2_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        CMS_ggzzbkg_a2 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a2",122.8,0.,200.)
        name = "CMS_ggzzbkg_a3_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        CMS_ggzzbkg_a3 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a3",0.03479,0.,1.)
        name = "CMS_ggzzbkg_a4_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts )
        CMS_ggzzbkg_a4 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a4",185.5,0.,200.)
        name = "CMS_ggzzbkg_a5_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        CMS_ggzzbkg_a5 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a5",12.67,0.,200.)
        name = "CMS_ggzzbkg_a6_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        CMS_ggzzbkg_a6 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a6",34.81,0.,100.)
        name = "CMS_ggzzbkg_a7_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        CMS_ggzzbkg_a7 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a7",0.1393,0.,1.)
        name = "CMS_ggzzbkg_a8_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts ) 
        CMS_ggzzbkg_a8 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a8",66.,0.,200.)
        name = "CMS_ggzzbkg_a9_{0:.0f}_{1:.0f}".format( self.channel, self.sqrts )
        CMS_ggzzbkg_a9 = ROOT.RooRealVar(name,"CMS_ggzzbkg_a9",0.07191,0.,1.)
        
        CMS_ggzzbkg_a0.setVal(theInputs['ggZZshape_a0'])
        CMS_ggzzbkg_a1.setVal(theInputs['ggZZshape_a1'])
        CMS_ggzzbkg_a2.setVal(theInputs['ggZZshape_a2'])
        CMS_ggzzbkg_a3.setVal(theInputs['ggZZshape_a3'])
        CMS_ggzzbkg_a4.setVal(theInputs['ggZZshape_a4'])
        CMS_ggzzbkg_a5.setVal(theInputs['ggZZshape_a5'])
        CMS_ggzzbkg_a6.setVal(theInputs['ggZZshape_a6'])
        CMS_ggzzbkg_a7.setVal(theInputs['ggZZshape_a7'])
        CMS_ggzzbkg_a8.setVal(theInputs['ggZZshape_a8'])
        CMS_ggzzbkg_a9.setVal(theInputs['ggZZshape_a9'])
        
        CMS_ggzzbkg_a0.setConstant(True)
        CMS_ggzzbkg_a1.setConstant(True)
        CMS_ggzzbkg_a2.setConstant(True)
        CMS_ggzzbkg_a3.setConstant(True)
        CMS_ggzzbkg_a4.setConstant(True)
        CMS_ggzzbkg_a5.setConstant(True)
        CMS_ggzzbkg_a6.setConstant(True)
        CMS_ggzzbkg_a7.setConstant(True)
        CMS_ggzzbkg_a8.setConstant(True)
        CMS_ggzzbkg_a9.setConstant(True)

        if (DEBUG) :
            print "ggZZshape_a0 = ",theInputs['ggZZshape_a0']
            print "ggZZshape_a1 = ",theInputs['ggZZshape_a1']
            print "ggZZshape_a2 = ",theInputs['ggZZshape_a2']
            print "ggZZshape_a3 = ",theInputs['ggZZshape_a3']
            print "ggZZshape_a4 = ",theInputs['ggZZshape_a4']
            print "ggZZshape_a5 = ",theInputs['ggZZshape_a5']
            print "ggZZshape_a6 = ",theInputs['ggZZshape_a6']
            print "ggZZshape_a7 = ",theInputs['ggZZshape_a7']
            print "ggZZshape_a8 = ",theInputs['ggZZshape_a8']
            print "ggZZshape_a9 = ",theInputs['ggZZshape_a9']
                   
        
        bkg_ggzz = ROOT.RooggZZPdf_v2("bkg_ggzzTmp","bkg_ggzzTmp",CMS_zz4l_mass,CMS_ggzzbkg_a0,CMS_ggzzbkg_a1,CMS_ggzzbkg_a2,CMS_ggzzbkg_a3,CMS_ggzzbkg_a4,CMS_ggzzbkg_a5,CMS_ggzzbkg_a6,CMS_ggzzbkg_a7,CMS_ggzzbkg_a8,CMS_ggzzbkg_a9)
    
        ## Reducible backgrounds
        val_meanL = float(theInputs['zjetsShape_mean'])
        val_sigmaL = float(theInputs['zjetsShape_sigma'])
 
        name = "mlZjet_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        mlZjet = ROOT.RooRealVar(name,"mean landau Zjet",val_meanL)
        name = "slZjet_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        slZjet = ROOT.RooRealVar(name,"sigma landau Zjet",val_sigmaL)
        bkg_zjets = ROOT.RooLandau("bkg_zjetsTmp","bkg_zjetsTmp",CMS_zz4l_mass,mlZjet,slZjet) 


 
	bkg_qqzzErr = ROOT.RooProdPdf("bkg_qqzzErr","bkg_qqzzErr", ROOT.RooArgSet(bkg_qqzz), ROOT.RooFit.Conditional(ROOT.RooArgSet(pdfErrB), ROOT.RooArgSet(CMS_zz4l_massErr)));
	bkg_ggzzErr = ROOT.RooProdPdf("bkg_ggzzErr","bkg_ggzzErr", ROOT.RooArgSet(bkg_ggzz), ROOT.RooFit.Conditional(ROOT.RooArgSet(pdfErrB), ROOT.RooArgSet(CMS_zz4l_massErr)));
	bkg_zjetsErr = ROOT.RooProdPdf("bkg_zjetsErr","bkg_zjetsErr", ROOT.RooArgSet(bkg_zjets), ROOT.RooFit.Conditional(ROOT.RooArgSet(pdfErrB), ROOT.RooArgSet(CMS_zz4l_massErr)));

      ## ----------------- 2D BACKGROUND SHAPES --------------- ##
        
        templateBkgName = "{0}/Dbackground_qqZZ_{1}.root".format(self.templateDir ,self.appendName)
        bkgTempFile = ROOT.TFile(templateBkgName)
        bkgTemplate = bkgTempFile.Get("h_mzzD")
        bkgTemplate_Up = bkgTempFile.Get("h_mzzD_up")
        bkgTemplate_Down = bkgTempFile.Get("h_mzzD_dn")
        
        TemplateName = "bkgTempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        bkgTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,D),bkgTemplate)
        TemplateName = "bkgTempDataHist_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        bkgTempDataHist_Up = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,D),bkgTemplate_Up)
        TemplateName = "bkgTempDataHist_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        bkgTempDataHist_Down = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,D),bkgTemplate_Down)
        
        templateggBkgName = "{0}/Dbackground_ggZZ_{1}.root".format(self.templateDir ,self.appendName)
        ggbkgTempFile = ROOT.TFile(templateggBkgName)
        ggbkgTemplate = ggbkgTempFile.Get("h_mzzD")
        TemplateName = "ggbkgTempDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        ggbkgTempDataHist = ROOT.RooDataHist(TemplateName,TemplateName,ROOT.RooArgList(CMS_zz4l_mass,D),ggbkgTemplate)
        
        TemplateName = "bkgTemplatePdf_qqzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        bkgTemplatePdf_qqzz = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),bkgTempDataHist)
        TemplateName = "bkgTemplatePdf_ggzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        bkgTemplatePdf_ggzz = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),ggbkgTempDataHist)
        TemplateName = "bkgTemplatePdf_zjets_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        bkgTemplatePdf_zjets = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),bkgTempDataHist)
        TemplateName,"bkgTemplatePdf_zjets_Up_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        bkgTemplatePdf_zjets_Up = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),bkgTempDataHist_Up)
        TemplateName = "bkgTemplatePdf_zjets_Down_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        bkgTemplatePdf_zjets_Down = ROOT.RooHistPdf(TemplateName,TemplateName,ROOT.RooArgSet(CMS_zz4l_mass,D),bkgTempDataHist_Down)
        
        funcList_zjets = ROOT.RooArgList()  
        morphBkgVarName = "CMS_zz4l_bkgMELA"    
        alphaMorphBkg = ROOT.RooRealVar(morphBkgVarName,morphBkgVarName,0,-20,20)
        morphVarListBkg = ROOT.RooArgList()
        
        if(self.bkgMorph):
            funcList_zjets.add(bkgTemplatePdf_zjets)
            funcList_zjets.add(bkgTemplatePdf_zjets_Up)
            funcList_zjets.add(bkgTemplatePdf_zjets_Down)  
            alphaMorphBkg.setConstant(False)
            morphVarListBkg.add(alphaMorphBkg)  
        else:
            funcList_zjets.add(bkgTemplatePdf_zjets)
            alphaMorphBkg.setConstant(True)


        TemplateName = "bkgTemplateMorphPdf_qqzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        bkgTemplateMorphPdf_qqzz = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,D,true,ROOT.RooArgList(bkgTemplatePdf_qqzz),ROOT.RooArgList(),1.0,1)
        TemplateName = "bkgTemplateMorphPdf_ggzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        bkgTemplateMorphPdf_ggzz = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,D,true,ROOT.RooArgList(bkgTemplatePdf_ggzz),ROOT.RooArgList(),1.0,1)
        TemplateName = "bkgTemplateMorphPdf_zjets_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        bkgTemplateMorphPdf_zjets = ROOT.FastVerticalInterpHistPdf2D(TemplateName,TemplateName,CMS_zz4l_mass,D,true,funcList_zjets,morphVarListBkg,1.0,1)

        bkg2d_qqzz = ROOT.RooProdPdf("bkg2d_qqzz","bkg2d_qqzz",ROOT.RooArgSet(self.getVariable(bkg_qqzzErr,bkg_qqzz,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(bkgTemplateMorphPdf_qqzz),ROOT.RooArgSet(D)))
        bkg2d_ggzz = ROOT.RooProdPdf("bkg2d_ggzz","bkg2d_ggzz",ROOT.RooArgSet(self.getVariable(bkg_ggzzErr,bkg_ggzz,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(bkgTemplateMorphPdf_ggzz),ROOT.RooArgSet(D)))
        bkg2d_zjets = ROOT.RooProdPdf("bkg2d_zjets","bkg2d_zjets",ROOT.RooArgSet(self.getVariable(bkg_zjetsErr,bkg_zjets,self.bIncludingError)),ROOT.RooFit.Conditional(ROOT.RooArgSet(bkgTemplateMorphPdf_zjets),ROOT.RooArgSet(D)))

        ## ----------------- SUPERMELA BACKGROUND SHAPES --------------- ##
        
        templateSDBkgName = "{0}/Dbackground_qqZZ_superMELA_{1}.root".format(self.templateDir ,self.appendName) 
        bkgTempSDFile = ROOT.TFile(templateSDBkgName)
        bkgTemplateSD = bkgTempSDFile.Get("hSuperD_bkg") 
        
        TemplateSDName = "bkgTempSDDataHist_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        bkgTempSDDataHist = ROOT.RooDataHist(TemplateSDName,TemplateSDName,ROOT.RooArgList(SD),bkgTemplateSD)
        
        TemplateSDName = "bkgTemplateSDPdf_qqzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        bkgTemplateSDPdf_qqzz = ROOT.RooHistPdf(TemplateSDName,TemplateSDName,ROOT.RooArgSet(SD),bkgTempSDDataHist)

        TemplateSDName = "bkgTemplateSDPdf_ggzz_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        bkgTemplateSDPdf_ggzz = ROOT.RooHistPdf(TemplateSDName,TemplateSDName,ROOT.RooArgSet(SD),bkgTempSDDataHist)
        TemplateSDName = "bkgTemplateSDPdf_zjets_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)    
        bkgTemplateSDPdf_zjets = ROOT.RooHistPdf(TemplateSDName,TemplateSDName,ROOT.RooArgSet(SD),bkgTempSDDataHist)
        

        ## ----------------------- PLOTS FOR SANITY CHECKS -------------------------- ##
        
        czz = ROOT.TCanvas( "czz", "czz", 750, 700 )
        czz.cd()
        zzframe_s = CMS_zz4l_mass.frame(45)
        if self.bUseCBnoConvolution: super(RooDoubleCB,signalCB_ggH).plotOn(zzframe_s, ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineColor(1) )
        elif self.isHighMass : super(ROOT.RooFFTConvPdf,sig_ggH_HM).plotOn(zzframe_s, ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineColor(1) )
        else : super(ROOT.RooFFTConvPdf,sig_ggH).plotOn(zzframe_s, ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineColor(1) )
        super(ROOT.RooqqZZPdf_v2,bkg_qqzz).plotOn(zzframe_s, ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineColor(4) )
        super(ROOT.RooggZZPdf_v2,bkg_ggzz).plotOn(zzframe_s, ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineColor(6) )
        super(ROOT.RooLandau,bkg_zjets).plotOn(zzframe_s, ROOT.RooFit.LineStyle(2), ROOT.RooFit.LineColor(6) )
        zzframe_s.Draw()
        figName = "{0}/figs/mzz_{1}_{2}.png".format(self.outputDir, self.mH, self.appendName)
        czz.SaveAs(figName)
        del czz
        
        ## ------------------- LUMI -------------------- ##
        
        rrvLumi = ROOT.RooRealVar("cmshzz4l_lumi","cmshzz4l_lumi",self.lumi)  
        
        ## ----------------------- SIGNAL RATES ----------------------- ##
        
        CMS_zz4l_mass.setRange("shape",self.low_M,self.high_M)
        
        fr_low_M = self.low_M
        fr_high_M = self.high_M        
        if (self.mH >= 450): 
            fr_low_M = 100
            fr_high_M = 1000
        if (self.mH >= 750):
            fr_low_M = 100
            fr_high_M = 1400
            

        CMS_zz4l_mass.setRange("fullrangesignal",fr_low_M,fr_high_M)
        CMS_zz4l_mass.setRange("fullrange",100,1400)
        
        rfvCsFilter = RooFormulaVar()
        filterName = "cmshzz4l_csFilter_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        if(self.sqrts == 7): 
            rfvCsFilter = ROOT.RooFormulaVar(filterName,"0.5+0.5*TMath::Erf((@0 - 80.85)/50.42)", ROOT.RooArgList(self.MH) )
        else:
            rfvCsFilter = ROOT.RooFormulaVar(filterName,"@0",ROOT.RooArgList(one))

        if(DEBUG):
            print "@@@@@@@ rfvCsFilter = ",rfvCsFilter.getVal()

        sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_a1".format(self.channel,self.sqrts)
        rrva1 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_a1'])
        sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_a2".format(self.channel,self.sqrts)
        rrva2 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_a2'])
        sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_a3".format(self.channel,self.sqrts)
        rrva3 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_a3'])
        sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_a4".format(self.channel,self.sqrts)
        rrva4 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_a4'])
        sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_b1".format(self.channel,self.sqrts)
        rrvb1 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_b1'])
        sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_b2".format(self.channel,self.sqrts)
        rrvb2 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_b2'])
        sigEffName = "hzz4lsigeff_{0:.0f}_{1:.0f}_b3".format(self.channel,self.sqrts)
        rrvb3 = ROOT.RooRealVar(sigEffName,sigEffName, theInputs['sigEff_b3'])

        if(DEBUG):
            print "sigEff_a1 = ",theInputs['sigEff_a1']
            print "sigEff_a2 = ",theInputs['sigEff_a2']
            print "sigEff_a3 = ",theInputs['sigEff_a3']
            print "sigEff_a4 = ",theInputs['sigEff_a4']
            print "sigEff_b1 = ",theInputs['sigEff_b1']
            print "sigEff_b2 = ",theInputs['sigEff_b2']
            print "sigEff_b3 = ",theInputs['sigEff_b3']
           
    
        sigEffName = "hzz4lsignaleff_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)

        rfvSigEff = ROOT.RooFormulaVar(sigEffName,"(@0+@1*TMath::Erf((@7-@2)/@3))*(@4+@5*@7+@6*@7*@7)", ROOT.RooArgList(rrva1,rrva2,rrva3,rrva4,rrvb1,rrvb2,rrvb3,self.MH))
        ## following printout is needed ,  dont remove it
        print " @@@@@@@@ sigeff ",rfvSigEff.getVal()
    
        CS_ggH = myCSW.HiggsCS(1,self.mH,self.sqrts)
        CS_VBF = myCSW.HiggsCS(2,self.mH,self.sqrts)
        CS_WH = myCSW.HiggsCS(3,self.mH,self.sqrts)
        CS_ZH = myCSW.HiggsCS(4,self.mH,self.sqrts)
        CS_ttH = myCSW.HiggsCS(5,self.mH,self.sqrts)
    
        BRH2e2mu = myCSW.HiggsBR(13,self.mH)
        BRH4mu = myCSW.HiggsBR(12,self.mH)
        BRH4e = myCSW.HiggsBR(12,self.mH)
        BR = 0.0
        if( self.channel == self.ID_4mu ): BR = BRH4mu
        if( self.channel == self.ID_4e ): BR = BRH4e
        if( self.channel == self.ID_2e2mu ): BR = BRH2e2mu
    
        sigEfficiency = rfvSigEff.getVal()

        if(DEBUG):
            print "CS_ggH: ",CS_ggH,", CS_VBF: ",CS_VBF,", CS_WH: ",CS_WH,", CS_ZH: ",CS_ZH
            print ", CS_ttH: ",CS_ttH,", BRH2e2mu: ",BRH2e2mu,", BRH4mu: ",BRH4mu,", BRH4e: ",BRH4e

        csCorrection = 1.0
        if(self.sqrts == 7): csCorrection = self.csFilter(self.mH)

        ## SIG YIELDS
        sigRate_ggH = csCorrection*CS_ggH*BR*sigEfficiency*1000.*self.lumi
        sigRate_VBF = csCorrection*CS_VBF*BR*sigEfficiency*1000.*self.lumi
        sigRate_WH = csCorrection*CS_WH*BR*sigEfficiency*1000.*self.lumi
        sigRate_ZH = csCorrection*CS_ZH*BR*sigEfficiency*1000.*self.lumi
        sigRate_ttH = csCorrection*CS_ttH*BR*sigEfficiency*1000.*self.lumi
       
        tmpNormSigNoConv = signalCB_ggH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("fullrangesignal") ).getVal()
        tmpNormSigConv = sig_ggH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("fullrangesignal") ).getVal()
        tmpNormSigHM   = sig_ggH_HM.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("fullrangesignal") ).getVal()
      
        normalizationSignal = 0.0
        if self.isHighMass : normalizationSignal = tmpNormSigHM
        else : normalizationSignal = self.getVariable(tmpNormSigNoConv,tmpNormSigConv,self.bUseCBnoConvolution)
            
        print "#################### ",signalCB_ggH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("fullrangesignal") ).getVal()
        print "#################### ",signalCB_ggH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        print "#################### ",sig_ggH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("fullrangesignal") ).getVal()
        
        sclFactorSig_ggH = sigRate_ggH/normalizationSignal
        sclFactorSig_VBF = sigRate_VBF/normalizationSignal
        sclFactorSig_WH = sigRate_WH/normalizationSignal
        sclFactorSig_ZH = sigRate_ZH/normalizationSignal
        sclFactorSig_ttH = sigRate_ttH/normalizationSignal

        integral_ggH = 0.0
        integral_VBF = 0.0
        integral_WH  = 0.0
        integral_ZH  = 0.0
        integral_ttH = 0.0

        if self.isHighMass : integral_ggH = sig_ggH_HM.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        else : integral_ggH = self.getVariable(signalCB_ggH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),sig_ggH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),self.bUseCBnoConvolution)

        if self.isHighMass : integral_VBF = sig_VBF_HM.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        else : integral_VBF = self.getVariable(signalCB_VBF.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),sig_VBF.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),self.bUseCBnoConvolution)

        if self.isHighMass : integral_WH = sig_WH_HM.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        else : integral_WH = self.getVariable(signalCB_WH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),sig_WH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),self.bUseCBnoConvolution)

        if self.isHighMass : integral_ZH = sig_ZH_HM.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        else : integral_ZH = self.getVariable(signalCB_ZH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),sig_ZH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),self.bUseCBnoConvolution)

        if self.isHighMass : integral_ttH = sig_ttH_HM.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        else : integral_ttH = self.getVariable(signalCB_ttH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),sig_ttH.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal(),self.bUseCBnoConvolution)
        
        sigRate_ggH_Shape = sclFactorSig_ggH*integral_ggH
        sigRate_VBF_Shape = sclFactorSig_VBF*integral_VBF
        sigRate_WH_Shape = sclFactorSig_WH*integral_WH
        sigRate_ZH_Shape = sclFactorSig_ZH*integral_ZH
        sigRate_ttH_Shape = sclFactorSig_ttH*integral_ttH
        
        
        normSigName = "cmshzz4l_normalizationSignal_{0:.0f}_{1:.0f}".format(self.channel,self.sqrts)
        rrvNormSig = ROOT.RooRealVar()

        if self.isHighMass :
            rrvNormSig = ROOT.RooRealVar(normSigName,normSigName, sig_ggH_HM.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass)).getVal())
        else :
            rrvNormSig = ROOT.RooRealVar(normSigName,normSigName, self.getVariable(signalCB_ggH.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass)).getVal(),sig_ggH.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass)).getVal(),self.bUseCBnoConvolution))
        rrvNormSig.setConstant(True)

      #  rfvSigRate_ggH = ROOT.RooFormulaVar("ggH_norm","@0*@1*@2*1000*{0}*{2}/{1}".format(self.lumi,rrvNormSig.getVal(),self.getVariable(signalCB_ggH.createIntegral(RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),sig_ggH.createIntegral(RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),self.bUseCBnoConvolution)),ROOT.RooArgList(rfvCsFilter,rfvSigEff, rhfXsBrFuncV_1))

        rfvSigRate_ggH = ROOT.RooFormulaVar("ggH_norm","@0*@1*@2*1000*{0}*{2}/{1}".format(self.lumi,rrvNormSig.getVal(),integral_ggH),ROOT.RooArgList(rfvCsFilter,rfvSigEff, rhfXsBrFuncV_1))

        print "Compare integrals: integral_ggH=",integral_ggH,"  ; calculated=",self.getVariable(signalCB_ggH.createIntegral(RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),sig_ggH.createIntegral(RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),self.bUseCBnoConvolution)
        
        rfvSigRate_VBF = ROOT.RooFormulaVar("qqH_norm","@0*@1*@2*1000*{0}*{2}/{1}".format(self.lumi,rrvNormSig.getVal(),integral_VBF),ROOT.RooArgList(rfvCsFilter,rfvSigEff, rhfXsBrFuncV_2))
                         

        rfvSigRate_WH = ROOT.RooFormulaVar("WH_norm","@0*@1*@2*1000*{0}*{2}/{1}".format(self.lumi,rrvNormSig.getVal(),integral_WH),ROOT.RooArgList(rfvCsFilter,rfvSigEff, rhfXsBrFuncV_3))
                         

        rfvSigRate_ZH = ROOT.RooFormulaVar("ZH_norm","@0*@1*@2*1000*{0}*{2}/{1}".format(self.lumi,rrvNormSig.getVal(),integral_ZH),ROOT.RooArgList(rfvCsFilter,rfvSigEff, rhfXsBrFuncV_4))
                         

        rfvSigRate_ttH = ROOT.RooFormulaVar("ttH_norm","@0*@1*@2*1000*{0}*{2}/{1}".format(self.lumi,rrvNormSig.getVal(),integral_ttH),ROOT.RooArgList(rfvCsFilter,rfvSigEff, rhfXsBrFuncV_5))
                         


        print signalCB_ggH.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass)).getVal(),"   ",sig_ggH.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass)).getVal()
        print signalCB_ggH.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal(),"   ",sig_ggH.createIntegral(ROOT.RooArgSet(CMS_zz4l_mass),ROOT.RooFit.Range("shape")).getVal()
        if (self.all_chan):
            print "Requested to sum up over the 5 chans: the norm in rfvSigRate_ggH should be the sum of the values of sigRate_XYZ_Shape variables:"
        print " @@@@@@@ norm sig = ",rrvNormSig.getVal()
        print " @@@@@@@ rfvSigRate_ggH = ",rfvSigRate_ggH.getVal()
        print " sigRate_ggH_Shape=",sigRate_ggH_Shape
        print " @@@@@@@ rfvSigRate_VBF = ",rfvSigRate_VBF.getVal()
        print " sigRate_VBF_Shape=",sigRate_VBF_Shape
        print " @@@@@@@ rfvSigRate_WH = ",rfvSigRate_WH.getVal()
        print " sigRate_WH_Shape=",sigRate_WH_Shape
        print " @@@@@@@ rfvSigRate_ZH = ",rfvSigRate_ZH.getVal()
        print " sigRate_ZH_Shape=",sigRate_ZH_Shape
        print " @@@@@@@ rfvSigRate_ttH = ",rfvSigRate_ttH.getVal()
        print " sigRate_ttH_Shape=",sigRate_ttH_Shape
        print "Sum of sigRate_XYZ_Shape=",sigRate_ggH_Shape+sigRate_VBF_Shape+sigRate_WH_Shape+sigRate_ZH_Shape+sigRate_ttH_Shape
        ## SET RATES TO 1 
        ## DC RATES WILL BE MULTIPLIED
        ## BY RATES IMPORTED TO WS
        sigRate_ggH_Shape = 1
        sigRate_VBF_Shape = 1
        sigRate_WH_Shape = 1
        sigRate_ZH_Shape = 1
        sigRate_ttH_Shape = 1

             
        ## ----------------------- BACKGROUND RATES ----------------------- ##

        ## rates per lumi for scaling
        bkgRate_qqzz = theInputs['qqZZ_rate']/theInputs['qqZZ_lumi']
        bkgRate_ggzz = theInputs['ggZZ_rate']/theInputs['qqZZ_lumi']
        bkgRate_zjets = theInputs['zjets_rate']/theInputs['zjets_lumi']
        
        ## Get Normalizations
        normalizationBackground_qqzz = bkg_qqzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("fullrange") ).getVal()
        normalizationBackground_ggzz = bkg_ggzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("fullrange") ).getVal()
        normalizationBackground_zjets = bkg_zjets.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("fullrange") ).getVal()
        
        sclFactorBkg_qqzz = self.lumi*bkgRate_qqzz/normalizationBackground_qqzz
        sclFactorBkg_ggzz = self.lumi*bkgRate_ggzz/normalizationBackground_ggzz
        sclFactorBkg_zjets = self.lumi*bkgRate_zjets/normalizationBackground_zjets
               
        bkgRate_qqzz_Shape = sclFactorBkg_qqzz * bkg_qqzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        bkgRate_ggzz_Shape = sclFactorBkg_ggzz * bkg_ggzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        bkgRate_zjets_Shape = sclFactorBkg_zjets * bkg_zjets.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("shape") ).getVal()
        
        if(DEBUG):
            print "Shape signal rate: ",sigRate_ggH_Shape,", background rate: ",bkgRate_qqzz_Shape,", ",bkgRate_zjets_Shape," in ",low_M," - ",high_M
            CMS_zz4l_mass.setRange("lowmassregion",100.,160.)
            bkgRate_qqzz_lowmassregion = sclFactorBkg_qqzz * bkg_qqzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("lowmassregion") ).getVal()
            bkgRate_ggzz_lowmassregion = sclFactorBkg_ggzz * bkg_ggzz.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("lowmassregion") ).getVal()
            bkgRate_zjets_lowmassregion = sclFactorBkg_zjets * bkg_zjets.createIntegral( ROOT.RooArgSet(CMS_zz4l_mass), ROOT.RooFit.Range("lowmassregion") ).getVal()
            lowmassyield = bkgRate_qqzz_lowmassregion + bkgRate_ggzz_lowmassregion + bkgRate_zjets_lowmassregion
            print "low mass yield: ",lowmassyield
        
        ## --------------------------- DATASET --------------------------- ##

        dataFileDir = "CMSdata"
        dataTreeName = "data_obs" 
        dataFileName = "{0}/hzz{1}_{2}.root".format(dataFileDir,self.appendName,self.lumi)
        if (DEBUG): print dataFileName," ",dataTreeName 
        data_obs_file = ROOT.TFile(dataFileName)

        print data_obs_file.Get(dataTreeName)
        
        if not (data_obs_file.Get(dataTreeName)):
            print "File, \"",dataFileName,"\", or tree, \"",dataTreeName,"\", not found" 
            print "Exiting..."
            sys.exit()
        
        data_obs_tree = data_obs_file.Get(dataTreeName)
        data_obs = ROOT.RooDataSet()
        datasetName = "data_obs_{0}".format(self.appendName)
        
        if (self.is2D == 0):
            if(self.bIncludingError): data_obs = ROOT.RooDataSet(datasetName,datasetName,data_obs_tree,ROOT.RooArgSet(CMS_zz4l_mass, CMS_zz4l_massErr))
            else: data_obs = ROOT.RooDataSet(datasetName,datasetName,data_obs_tree,ROOT.RooArgSet(CMS_zz4l_mass))
		
        if (self.is2D == 1):
            if(self.bIncludingError): data_obs = ROOT.RooDataSet(datasetName,datasetName,data_obs_tree,ROOT.RooArgSet(CMS_zz4l_mass,D, CMS_zz4l_massErr))
            else: data_obs = ROOT.RooDataSet(datasetName,datasetName,data_obs_tree,ROOT.RooArgSet(CMS_zz4l_mass,D))

        if (self.is2D == 2):
            data_obs = ROOT.RooDataSet(datasetName,datasetName,data_obs_tree,ROOT.RooArgSet(SD))
            
        ## --------------------------- WORKSPACE -------------------------- ##

        endsInP5 = False
        tmpMH = self.mH
        if ( math.fabs(math.floor(tmpMH)-self.mH) == 0.5): endsInP5 = True
        if (DEBUG): print "ENDS IN P5  ",endsInP5

        name_Shape = ""
        name_ShapeWS = ""
        name_ShapeWS2 = ""
        name_ShapeWSXSBR = ""
        
        if (endsInP5): name_Shape = "{0}/HCG/{1}/hzz4l_{2}S_{3:.0f}TeV.txt".format(self.outputDir,self.mH,self.appendName,self.sqrts)
        else: name_Shape = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV.txt".format(self.outputDir,self.mH,self.appendName,self.sqrts)

        if (endsInP5): name_ShapeWS = "{0}/HCG/{1}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts)
        else: name_ShapeWS = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts)

        if (endsInP5): name_ShapeWSXSBR = "{0}/HCG_XSxBR/{1}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts)
        else: name_ShapeWSXSBR = "{0}/HCG_XSxBR/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV.input.root".format(self.outputDir,self.mH,self.appendName,self.sqrts)

        name_ShapeWS2 = "hzz4l_{0}S_{1:.0f}TeV.input.root".format(self.appendName,self.sqrts)

        if(DEBUG): print name_Shape,"  ",name_ShapeWS2
        
        w = ROOT.RooWorkspace("w","w")
        
        w.importClassCode(RooqqZZPdf_v2.Class(),True)
        w.importClassCode(RooggZZPdf_v2.Class(),True)
        w.importClassCode(RooRelBWUFParam.Class(),True)
        w.importClassCode(RooDoubleCB.Class(),True)
        w.importClassCode(RooFormulaVar.Class(),True)
        if self.isHighMass :
            w.importClassCode(RooRelBWHighMass.Class(),True)

        if( FactorizedShapes ):
            if( self.channel == self.ID_4mu ):
                w.importClassCode(RooFourMuMassShapePdf2.Class(),True)
                w.importClassCode(RooFourMuMassRes.Class(),True)
            elif( self.channel == self.ID_4e ):
                w.importClassCode(RooFourEMassShapePdf2.Class(),True)
                w.importClassCode(RooFourEMassRes.Class(),True)
            elif( self.channel == self.ID_2e2mu ):
                w.importClassCode(RooTwoETwoMuMassShapePdf2.Class(),True)
                w.importClassCode(RooTwoETwoMuMassRes.Class(),True)
            
                
                
        getattr(w,'import')(data_obs,ROOT.RooFit.Rename("data_obs")) ### Should this be renamed?
    
        if(self.bUseCBnoConvolution) :
            if (self.is2D == 0):
		if not self.bIncludingError:
                	signalCB_ggH.SetNameTitle("ggH","ggH")
                	signalCB_VBF.SetNameTitle("qqH","qqH")
                	signalCB_WH.SetNameTitle("WH","WH")
                	signalCB_ZH.SetNameTitle("ZH","ZH")
                	signalCB_ttH.SetNameTitle("ttH","ttH")
                
                	getattr(w,'import')(signalCB_ggH, ROOT.RooFit.RecycleConflictNodes())
                	getattr(w,'import')(signalCB_VBF, ROOT.RooFit.RecycleConflictNodes())
               		getattr(w,'import')(signalCB_WH, ROOT.RooFit.RecycleConflictNodes())
                	getattr(w,'import')(signalCB_ZH, ROOT.RooFit.RecycleConflictNodes())
                	getattr(w,'import')(signalCB_ttH, ROOT.RooFit.RecycleConflictNodes())
		else:
                	sig_ggHErr.SetNameTitle("ggH","ggH")
                	sig_VBFErr.SetNameTitle("qqH","qqH")
                	sig_WHErr.SetNameTitle("WH","WH")
                	sig_ZHErr.SetNameTitle("ZH","ZH")
                	sig_ttHErr.SetNameTitle("ttH","ttH")
                
                	getattr(w,'import')(sig_ggHErr, ROOT.RooFit.RecycleConflictNodes())
                	getattr(w,'import')(sig_VBFErr, ROOT.RooFit.RecycleConflictNodes())
               		getattr(w,'import')(sig_WHErr, ROOT.RooFit.RecycleConflictNodes())
                	getattr(w,'import')(sig_ZHErr, ROOT.RooFit.RecycleConflictNodes())
                	getattr(w,'import')(sig_ttHErr, ROOT.RooFit.RecycleConflictNodes())
                
            if (self.is2D == 1):
                sigCB2d_ggH.SetNameTitle("ggH","ggH")
                sigCB2d_VBF.SetNameTitle("qqH","qqH")
                sigCB2d_WH.SetNameTitle("WH","WH")
                sigCB2d_ZH.SetNameTitle("ZH","ZH")
                sigCB2d_ttH.SetNameTitle("ttH","ttH")
                
                getattr(w,'import')(sigCB2d_ggH, ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(sigCB2d_VBF, ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(sigCB2d_WH, ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(sigCB2d_ZH, ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(sigCB2d_ttH, ROOT.RooFit.RecycleConflictNodes())
                if(self.isAltSig):
                    sigCB2d_ggH_ALT.SetNameTitle("ggH{0}".format(self.appendHypType),"ggH{0}".format(self.appendHypType))
                    getattr(w,'import')(sigCB2d_ggH_ALT, ROOT.RooFit.RecycleConflictNodes())

            if (self.is2D == 2):
                sigTemplateSDPdf_ggH.SetNameTitle("ggH","ggH")
                sigTemplateSDPdf_VBF.SetNameTitle("qqH","qqH")
                sigTemplateSDPdf_WH.SetNameTitle("WH","WH")
                sigTemplateSDPdf_ZH.SetNameTitle("ZH","ZH")
                sigTemplateSDPdf_ttH.SetNameTitle("ttH","ttH")
                
                getattr(w,'import')(sigTemplateSDPdf_ggH, ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(sigTemplateSDPdf_VBF, ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(sigTemplateSDPdf_WH, ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(sigTemplateSDPdf_ZH, ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(sigTemplateSDPdf_ttH, ROOT.RooFit.RecycleConflictNodes())

        else:
                
            if (self.is2D == 0):

                if self.isHighMass:
                    sig_ggH_HM.SetNameTitle("ggH","ggH")
                    sig_VBF_HM.SetNameTitle("qqH","qqH")
                    sig_WH_HM.SetNameTitle("WH","WH")
                    sig_ZH_HM.SetNameTitle("ZH","ZH")
                    sig_ttH_HM.SetNameTitle("ttH","ttH")
                    
                    getattr(w,'import')(sig_ggH_HM, ROOT.RooFit.RecycleConflictNodes())
                    getattr(w,'import')(sig_VBF_HM, ROOT.RooFit.RecycleConflictNodes())
                    getattr(w,'import')(sig_WH_HM, ROOT.RooFit.RecycleConflictNodes())
                    getattr(w,'import')(sig_ZH_HM, ROOT.RooFit.RecycleConflictNodes())
                    getattr(w,'import')(sig_ttH_HM, ROOT.RooFit.RecycleConflictNodes())

                else :
                    sig_ggH.SetNameTitle("ggH","ggH")
                    sig_VBF.SetNameTitle("qqH","qqH")
                    sig_WH.SetNameTitle("WH","WH")
                    sig_ZH.SetNameTitle("ZH","ZH")
                    sig_ttH.SetNameTitle("ttH","ttH")
                    
                    getattr(w,'import')(sig_ggH, ROOT.RooFit.RecycleConflictNodes())
                    getattr(w,'import')(sig_VBF, ROOT.RooFit.RecycleConflictNodes())
                    getattr(w,'import')(sig_WH, ROOT.RooFit.RecycleConflictNodes())
                    getattr(w,'import')(sig_ZH, ROOT.RooFit.RecycleConflictNodes())
                    getattr(w,'import')(sig_ttH, ROOT.RooFit.RecycleConflictNodes())
                    
            if (self.is2D == 1):
                sig2d_ggH.SetNameTitle("ggH","ggH")
                sig2d_VBF.SetNameTitle("qqH","qqH")
                sig2d_WH.SetNameTitle("WH","WH")
                sig2d_ZH.SetNameTitle("ZH","ZH")
                sig2d_ttH.SetNameTitle("ttH","ttH")     
                
                getattr(w,'import')(sig2d_ggH,ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(sig2d_VBF,ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(sig2d_WH,ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(sig2d_ZH,ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(sig2d_ttH,ROOT.RooFit.RecycleConflictNodes())            
        
                if(self.isAltSig):
                    sig2d_ggH_ALT.SetNameTitle("ggH{0}".format(self.appendHypType),"ggH{0}".format(self.appendHypType))
                    getattr(w,'import')(sig2d_ggH_ALT, ROOT.RooFit.RecycleConflictNodes())

            if (self.is2D == 2): 
                sigTemplateSDPdf_ggH.SetNameTitle("ggH","ggH")
                sigTemplateSDPdf_VBF.SetNameTitle("qqH","qqH")
                sigTemplateSDPdf_WH.SetNameTitle("WH","WH")
                sigTemplateSDPdf_ZH.SetNameTitle("ZH","ZH")
                sigTemplateSDPdf_ttH.SetNameTitle("ttH","ttH")
                
                getattr(w,'import')(sigTemplateSDPdf_ggH, ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(sigTemplateSDPdf_VBF, ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(sigTemplateSDPdf_WH, ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(sigTemplateSDPdf_ZH, ROOT.RooFit.RecycleConflictNodes())
                getattr(w,'import')(sigTemplateSDPdf_ttH, ROOT.RooFit.RecycleConflictNodes())

 
        if (self.is2D == 0):
		if not self.bIncludingError:
			bkg_qqzz.SetNameTitle("bkg_qqzz","bkg_qqzz")
			bkg_ggzz.SetNameTitle("bkg_ggzz","bkg_ggzz")
			bkg_zjets.SetNameTitle("bkg_zjets","bkg_zjets")
            		getattr(w,'import')(bkg_qqzz, ROOT.RooFit.RecycleConflictNodes())
            		getattr(w,'import')(bkg_ggzz, ROOT.RooFit.RecycleConflictNodes())
            		getattr(w,'import')(bkg_zjets, ROOT.RooFit.RecycleConflictNodes())
		else:
			bkg_qqzzErr.SetNameTitle("bkg_qqzz","bkg_qqzz")
			bkg_ggzzErr.SetNameTitle("bkg_ggzz","bkg_ggzz")
			bkg_zjetsErr.SetNameTitle("bkg_zjets","bkg_zjets")
            		getattr(w,'import')(bkg_qqzzErr, ROOT.RooFit.RecycleConflictNodes())
            		getattr(w,'import')(bkg_ggzzErr, ROOT.RooFit.RecycleConflictNodes())
            		getattr(w,'import')(bkg_zjetsErr, ROOT.RooFit.RecycleConflictNodes())
            
        if (self.is2D == 1):
            getattr(w,'import')(bkg2d_qqzz,ROOT.RooFit.RecycleConflictNodes())
            getattr(w,'import')(bkg2d_ggzz,ROOT.RooFit.RecycleConflictNodes())
            getattr(w,'import')(bkg2d_zjets,ROOT.RooFit.RecycleConflictNodes())

        if (self.is2D == 2): 
            bkgTemplateSDPdf_qqzz.SetNameTitle("bkg_qqzz","bkg_qqzz")
            bkgTemplateSDPdf_ggzz.SetNameTitle("bkg_ggzz","bkg_ggzz")
            bkgTemplateSDPdf_zjets.SetNameTitle("bkg_zjets","bkg_zjets")
            getattr(w,'import')(bkgTemplateSDPdf_ggzz, ROOT.RooFit.RecycleConflictNodes())
            getattr(w,'import')(bkgTemplateSDPdf_qqzz, ROOT.RooFit.RecycleConflictNodes())
            getattr(w,'import')(bkgTemplateSDPdf_zjets, ROOT.RooFit.RecycleConflictNodes())

        
        getattr(w,'import')(rfvSigRate_ggH, ROOT.RooFit.RecycleConflictNodes())
        getattr(w,'import')(rfvSigRate_VBF, ROOT.RooFit.RecycleConflictNodes())
        getattr(w,'import')(rfvSigRate_WH, ROOT.RooFit.RecycleConflictNodes())
        getattr(w,'import')(rfvSigRate_ZH, ROOT.RooFit.RecycleConflictNodes())
        getattr(w,'import')(rfvSigRate_ttH, ROOT.RooFit.RecycleConflictNodes())
        if(self.isAltSig):
            rfvSigRate_ggH_ALT=ROOT.RooFormulaVar(rfvSigRate_ggH,"ggH{0}_norm".format(self.appendHypType))
            print 'Compare signal rates: STD=',rfvSigRate_ggH.getVal(),"   ALT=",rfvSigRate_ggH_ALT.getVal()
            getattr(w,'import')(rfvSigRate_ggH_ALT, ROOT.RooFit.RecycleConflictNodes())
            
        w.writeToFile(name_ShapeWS)
        w.writeToFile(name_ShapeWSXSBR)
        
        ## --------------------------- DATACARDS -------------------------- ##

        systematics.setSystematics(bkgRate_qqzz_Shape, bkgRate_ggzz_Shape, bkgRate_zjets_Shape)
        systematics_forXSxBR.setSystematics(bkgRate_qqzz_Shape, bkgRate_ggzz_Shape,bkgRate_zjets_Shape)

        ## If the channel is not declared in inputs, set rate = 0
        if not self.ggH_chan and not self.all_chan :  sigRate_ggH_Shape = 0
        if not self.qqH_chan:  sigRate_VBF_Shape = 0
        if not self.WH_chan:   sigRate_WH_Shape = 0
        if not self.ZH_chan:   sigRate_ZH_Shape = 0
        if not self.ttH_chan:  sigRate_ttH_Shape = 0

        if not self.qqZZ_chan:  bkgRate_qqzz_Shape = 0
        if not self.ggZZ_chan:  bkgRate_ggzz_Shape = 0
        if not self.zjets_chan: bkgRate_zjets_Shape = 0

        rates = {}
        rates['ggH'] = sigRate_ggH_Shape
        rates['qqH'] = sigRate_VBF_Shape
        rates['WH']  = sigRate_WH_Shape
        rates['ZH']  = sigRate_ZH_Shape
        rates['ttH'] = sigRate_ttH_Shape

        rates['qqZZ']  = bkgRate_qqzz_Shape
        rates['ggZZ']  = bkgRate_ggzz_Shape
        rates['zjets'] = bkgRate_zjets_Shape
        rates['ttbar'] = 0
        rates['zbb']   = 0
        

        ## Write Datacards
        fo = open( name_Shape, "wb")
        self.WriteDatacard(fo,theInputs, name_ShapeWS2, rates, data_obs.numEntries(), self.is2D )
        systematics.WriteSystematics(fo, theInputs)
        systematics.WriteShapeSystematics(fo,theInputs)
        fo.close()

        if(self.isAltSig):
            if (endsInP5): name_Shape = "{0}/HCG/{1}/hzz4l_{2}S_{3:.0f}TeV{4}.txt".format(self.outputDir,self.mH,self.appendName,self.sqrts,self.appendHypType)
            else: name_Shape = "{0}/HCG/{1:.0f}/hzz4l_{2}S_{3:.0f}TeV{4}.txt".format(self.outputDir,self.mH,self.appendName,self.sqrts,self.appendHypType)
            fo = open( name_Shape, "wb")
            self.WriteDatacard(fo,theInputs, name_ShapeWS2, rates, data_obs.numEntries(), self.is2D,True,self.appendHypType )
            systematics.WriteSystematics(fo, theInputs)
            systematics.WriteShapeSystematics(fo,theInputs)
            fo.close()

        ## forXSxBR
        if (endsInP5): name_Shape = "{0}/HCG_XSxBR/{2}/hzz4l_{1}S_{3:.0f}TeV.txt".format(self.outputDir,self.appendName,self.mH,self.sqrts)	
        else: name_Shape = "{0}/HCG_XSxBR/{2:.0f}/hzz4l_{1}S_{3:.0f}TeV.txt".format(self.outputDir,self.appendName,self.mH,self.sqrts)
            
        fo = open( name_Shape, "wb" )
        self.WriteDatacard(fo, theInputs,name_ShapeWS2, rates, data_obs.numEntries(), self.is2D )
        systematics_forXSxBR.WriteSystematics(fo, theInputs)
        systematics_forXSxBR.WriteShapeSystematics(fo,theInputs)
        fo.close()

        if(self.isAltSig):
            if (endsInP5): name_Shape = "{0}/HCG_XSxBR/{2}/hzz4l_{1}S_{3:.0f}TeV{4}.txt".format(self.outputDir,self.appendName,self.mH,self.sqrts,self.appendHypType)	
            else: name_Shape = "{0}/HCG_XSxBR/{2:.0f}/hzz4l_{1}S_{3:.0f}TeV{4}.txt".format(self.outputDir,self.appendName,self.mH,self.sqrts,self.appendHypType)
            fo = open( name_Shape, "wb")
            self.WriteDatacard(fo,theInputs,name_ShapeWS2,rates,data_obs.numEntries(),self.is2D,True,self.appendHypType )
            systematics.WriteSystematics(fo, theInputs)
            systematics.WriteShapeSystematics(fo,theInputs)
            fo.close()
        


    def WriteDatacard(self,file,theInputs,nameWS,theRates,obsEvents,is2D,isAltCard=False,AltLabel=""):

        numberSig = self.numberOfSigChan(theInputs)
        numberBg  = self.numberOfBgChan(theInputs)

        file.write("imax 1\n")
        file.write("jmax {0}\n".format(numberSig+numberBg-1))
        file.write("kmax *\n")
        
        file.write("------------\n")
        file.write("shapes * * {0} w:$PROCESS \n".format(nameWS))
        file.write("------------\n")
        
        file.write("bin a{0} \n".format(self.channel))
        file.write("observation {0} \n".format(obsEvents))
        
        file.write("------------\n")
        file.write("## mass window [{0},{1}] \n".format(self.low_M,self.high_M))
        file.write("bin ")        

        channelList=['ggH','qqH','WH','ZH','ttH','qqZZ','ggZZ','zjets','ttbar','zbb']

        channelName1D=['ggH','qqH','WH','ZH','ttH','bkg_qqzz','bkg_ggzz','bkg_zjets','bkg_ttbar','bkg_zbb']
        channelName2D=['ggH','qqH','WH','ZH','ttH','bkg2d_qqzz','bkg2d_ggzz','bkg2d_zjets','bkg2d_ttbar','bkg2d_zbb']

#            channelList=['ggH{0}'.format(AltLabel),'qqZZ','ggZZ','zjets','ttbar','zbb']

        if theInputs["all"]:
            channelList=['ggH','qqZZ','ggZZ','zjets','ttbar','zbb']
            if isAltCard :
                channelName2D=['ggH{0}'.format(AltLabel),'bkg2d_qqzz','bkg2d_ggzz','bkg2d_zjets','bkg2d_ttbar','bkg2d_zbb']
            else:
                channelName2D=['ggH','bkg2d_qqzz','bkg2d_ggzz','bkg2d_zjets','bkg2d_ttbar','bkg2d_zbb']
          
         
        for chan in channelList:
            if theInputs[chan]:
                file.write("a{0} ".format(self.channel))
            else:
                if chan.startswith("ggH") and theInputs["all"] :
                    file.write("a{0} ".format(self.channel))
        file.write("\n")
                                        
        file.write("process ")

        i=0
        if not (self.is2D == 1):
            for chan in channelList:
                if theInputs[chan]:
                    file.write("{0} ".format(channelName1D[i]))
                i+=1
        else:
            for chan in channelList:
#                print 'checking if ',chan,' is in the list of to-do'
                if theInputs[chan]:
                    file.write("{0} ".format(channelName2D[i]))
#                    print 'writing in card index=',i,'  chan=',chan
                    i+=1
                else:
                    if chan.startswith("ggH") and theInputs["all"] :
                        file.write("{0} ".format(channelName2D[i]))
#                        print 'writing in card TOTAL SUM, index=',i,'  chan=',chan,'  ',channelName2D[i]
                        i+=1
        
        file.write("\n")
            
        processLine = "process "

        for x in range(-numberSig+1,1):
            processLine += "{0} ".format(x)

        for y in range(1,numberBg+1):
            processLine += "{0} ".format(y)

        file.write(processLine)
        file.write("\n")
            
        file.write("rate ")
        for chan in channelList:
            if theInputs[chan] or (chan.startswith("ggH") and theInputs["all"]):
                file.write("{0:.4f} ".format(theRates[chan]))
        file.write("\n")
        file.write("------------\n")


        
    def numberOfSigChan(self,inputs):

        counter=0

        if inputs['ggH']: counter+=1
        if inputs['qqH']: counter+=1
        if inputs['WH']:  counter+=1
        if inputs['ZH']:  counter+=1
        if inputs['ttH']: counter+=1
        if inputs['all']: counter+=1
        
        return counter

    def numberOfBgChan(self,inputs):

        counter=0

        if inputs['qqZZ']:  counter+=1
        if inputs['ggZZ']:  counter+=1
        if inputs['zjets']: counter+=1
        if inputs['ttbar']: counter+=1
        if inputs['zbb']:   counter+=1
        
        return counter

