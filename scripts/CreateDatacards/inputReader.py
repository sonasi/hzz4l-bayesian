#!/usr/bin/python
import os
import re
import math
import collections
from ROOT import *
from array import array

## ---------------------------------------------------------------
## card reader class
## ---------------------------------------------------------------

class inputReader:

    def __init__(self, inputTextFile):

        if not os.path.exists(inputTextFile):
            raise RuntimeError, "File {0} does not exist!!!".format(inputTextFile)
        
        # input file
        self.theInput = inputTextFile
        # model
        self.model = ""
        # decay channel
        self.decayChan = -999.9
        # lumi
        self.lumi = -999.9
        # sqrts
        self.sqrts = -999.9
        # channels
        self.all_chan = False
        self.ggH_chan = False
        self.qqH_chan = False
        self.WH_chan = False
        self.ZH_chan = False
        self.ttH_chan = False
        self.qqZZ_chan = False
        self.ggZZ_chan = False
        self.zjets_chan = False
        self.ttbar_chan = False
        self.zbb_chan = False
        # rates
        self.qqZZ_rate = -999.9
        self.ggZZ_rate = -999.9
        self.zjets_rate = -999.9
        self.ttbar_rate = -999.9
        self.zbb_rate = -999.9
        self.qqZZ_lumi = -999.9
        self.ggZZ_lumi = -999.9
        self.zjets_lumi = -999.9
        self.ttbar_lumi = -999.9
        self.zbb_lumi = -999.9
        # signal shapes
        self.useHighMassReweightedShapes = True
        self.n_CB_shape = -999.9
        self.alpha_CB_shape = -999.9
        self.n2_CB_shape = -999.9
        self.alpha2_CB_shape = -999.9
        self.mean_CB_shape = -999.9
        self.sigma_CB_shape = -999.9
        self.n_CB_shape_HM = -999.9
        self.alpha_CB_shape_HM = -999.9
        self.n2_CB_shape_HM = -999.9
        self.alpha2_CB_shape_HM = -999.9
        self.mean_CB_shape_HM = -999.9
        self.sigma_CB_shape_HM = -999.9
        self.gamma_BW_shape_HM = -999.9
        # signal efficiency params
        self.sigeff_a1 = -999.9
        self.sigeff_a2 = -999.9
        self.sigeff_a3 = -999.9
        self.sigeff_a4 = -999.9
        self.sigeff_b1 = -999.9
        self.sigeff_b2 = -999.9
        self.sigeff_b3 = -999.9
        # qqZZ shape
        self.qqZZshape_a0 = -999.9
        self.qqZZshape_a1 = -999.9
        self.qqZZshape_a2 = -999.9
        self.qqZZshape_a3 = -999.9
        self.qqZZshape_a4 = -999.9
        self.qqZZshape_a5 = -999.9
        self.qqZZshape_a6 = -999.9
        self.qqZZshape_a7 = -999.9
        self.qqZZshape_a8 = -999.9
        self.qqZZshape_a9 = -999.9
        self.qqZZshape_a10 = -999.9
        self.qqZZshape_a11 = -999.9
        self.qqZZshape_a12 = -999.9
        self.qqZZshape_a13 = -999.9
        # ggZZ shape
        self.ggZZshape_a0 = -999.9
        self.ggZZshape_a1 = -999.9
        self.ggZZshape_a2 = -999.9
        self.ggZZshape_a3 = -999.9
        self.ggZZshape_a4 = -999.9
        self.ggZZshape_a5 = -999.9
        self.ggZZshape_a6 = -999.9
        self.ggZZshape_a7 = -999.9
        self.ggZZshape_a8 = -999.9
        self.ggZZshape_a9 = -999.9
        # zjets shape
        self.zjetsShape_mean = -999.9
        self.zjetsShape_sigma = -999.9
        # systematics 
        self.zjetsKappaLow = -999.9
        self.zjetsKappaHigh = -999.9
        self.lumiUnc = -999.9
        self.muonFullUnc = -999.9
        self.muonFullUnc_HM = -999.9
        self.muonFullCutoff = -999.9
        self.muonTrigUnc = -999.9
        self.muonTrigUnc_HM = -999.9
        self.muonTrigCutoff = -999.9
        self.elecFullUnc = -999.9
        self.elecFullUnc_HM = -999.9
        self.elecFullCutoff = -999.9
        self.elecTrigUnc = -999.9
        self.elecTrigUnc_HM = -999.9
        self.elecTrigCutoff = -999.9

        self.CMS_zz4l_mean_m_sig = -999.9
        self.CMS_zz4l_sigma_m_sig = -999.9
        self.CMS_zz4l_mean_e_sig = -999.9
        self.CMS_zz4l_sigma_e_sig = -999.9
        self.CMS_zz4l_n_sig = -999.9
        self.CMSSW_zz4l_gamma_sig = -999.9

        self.useLumiUnc = False
        self.usePdf_gg = False
        self.usePdf_qqbar = False
        self.usePdf_hzz4l_accept = False
        self.useQCDscale_ggH = False
        self.useQCDscale_qqH = False
        self.useQCDscale_VH = False
        self.useQCDscale_ttH = False
        self.useTheoryUncXS_HighMH = False
        self.useQCDscale_ggVV = False
        self.useQCDscale_VV = False
        self.useBRhiggs_hzz4l = False
        self.useCMS_eff = False
        self.useCMS_hzz4l_Zjets = False 
        self.useCMS_zz4l_bkgMELA = False
        self.useCMS_zz4l_sigMELA = False
        self.useCMS_zz4l_mean = False
        self.useCMS_zz4l_sigma = False
        self.useCMS_zz4l_n = False
        self.useCMS_zz4l_gamma = False
        self.doHypTest = False
        self.altHypLabel = ""
        
	# ---  mekd stuffs
	self.mekd_sig_a0_shape = -999.
	self.mekd_sig_a1_shape = -999.
	self.mekd_sig_a2_shape = -999.
	self.mekd_sig_a3_shape = -999.
	self.mekd_sig_a4_shape = -999.
	self.mekd_qqZZ_a0_shape = -999.
	self.mekd_qqZZ_a1_shape = -999.
	self.mekd_qqZZ_a2_shape = -999.
	self.mekd_qqZZ_a3_shape = -999.
	self.mekd_qqZZ_a4_shape = -999.
	# --- end mekd

    def goodEntry(self,variable):
        if variable == -999.9:
            return False
        else:
            return True


    def parseBoolString(self,theString):

        return theString[0].upper()=='T'



    def readInputs(self):
        
        for line in open(self.theInput,'r'):
            f = line.split()
            if len(f) < 1: continue

            if f[0].startswith("#"): continue
            
            if f[0].lower().startswith("model"):
                
                if f[1].upper() == "SM": self.model = "SM"
                elif f[1].upper() == "SM4": self.model = "SM4"
                elif f[1].upper() == "FF" or f[1].upper() == "FP": self.model = "FF"
                else : raise RuntimeError, "Unknow model {0}, choices are SM, SM4, FF".format(f[1].upper()) 

            if f[0].lower().startswith("decay"):

                if f[1] == "4mu": self.decayChan = 1
                elif f[1] == "4e": self.decayChan = 2
                elif f[1] == "2e2mu": self.decayChan = 3
                elif f[1] == "2mu2e": self.decayChan = 3
                else : raise RuntimeError, "Unknown decay channel {0}, choices are 4mu, 4e, or 2e2mu".format(f[1])
                
            if f[0].lower().startswith("channels"):
                for chan in f:
                    if chan == f[0]: continue
                    if chan.lower().startswith("ggh"):     self.ggH_chan = True
                    elif chan.lower().startswith("qqh"):   self.qqH_chan = True
                    elif chan.lower().startswith("wh"):    self.WH_chan = True
                    elif chan.lower().startswith("zh"):    self.ZH_chan = True
                    elif chan.lower().startswith("tth"):   self.ttH_chan = True
                    elif chan.lower().startswith("qqzz"):  self.qqZZ_chan = True
                    elif chan.lower().startswith("ggzz"):  self.ggZZ_chan = True
                    elif chan.lower().startswith("zjets"): self.zjets_chan = True
                    elif chan.lower().startswith("ttbar"): self.ttbar_chan = True
                    elif chan.lower().startswith("zbb"):   self.zbb_chan = True
                    elif chan.lower().startswith("all"):   self.all_chan = True
                    else : raise RuntimeError, "Unknown channel {0}, choices are ggH, qqH, WH, ZH, ttH, qqZZ, ggZZ, zjets".format(chan)
          
            
            if f[0].lower().startswith("rate"):
                
                if f[1].lower().startswith("qqzz"):
                    self.qqZZ_rate = float(f[2])
                    if len(f) == 4: self.qqZZ_lumi = float(f[3])
                if f[1].lower().startswith("ggzz"):
                    self.ggZZ_rate = float(f[2])
                    if len(f) == 4: self.ggZZ_lumi = float(f[3])
                if f[1].lower().startswith("zjets"):
                    self.zjets_rate = float(f[2])
                    if len(f) == 4: self.zjets_lumi = float(f[3])
                if f[1].lower().startswith("ttbar"):
                    self.ttbar_rate = float(f[2])
                    if len(f) == 4: self.ttbar_lumi = float(f[3])
                if f[1].lower().startswith("zbb"):
                    self.zbb_rate = float(f[2])
                    if len(f) == 4: self.zbb_lumi = float(f[3])

            if f[0].lower().startswith("usehighmassreweightedshapes"):
                self.useHighMassReweightedShapes = True
                    
            if f[0].lower().startswith("signalshape"):

                if f[1].lower().startswith("n_cb"): 
                    if len(f) > 3 : raise RuntimeError, "{0} has a space in the formula!  Please check!".format(f[1])
                    else: self.n_CB_shape = f[2]
                if f[1].lower().startswith("alpha_cb"): 
                    if len(f) > 3 : raise RuntimeError, "{0} has a space in the formula!  Please check!".format(f[1])
                    else: self.alpha_CB_shape = f[2]
                if f[1].lower().startswith("n2_cb"): 
                    if len(f) > 3 : raise RuntimeError, "{0} has a space in the formula!  Please check!".format(f[1])
                    else: self.n2_CB_shape = f[2]
                if f[1].lower().startswith("alpha2_cb"): 
                    if len(f) > 3 : raise RuntimeError, "{0} has a space in the formula!  Please check!".format(f[1])
                    else: self.alpha2_CB_shape = f[2]
                if f[1].lower().startswith("mean_cb"): 
                    if len(f) > 3 : raise RuntimeError, "{0} has a space in the formula!  Please check!".format(f[1])
                    else: self.mean_CB_shape = f[2]
                if f[1].lower().startswith("sigma_cb"): 
                    if len(f) > 3 : raise RuntimeError, "{0} has a space in the formula!  Please check!".format(f[1])
                    else: self.sigma_CB_shape = f[2]
                if f[1].lower().startswith("mekd_sig_a0"): 
                    if len(f) > 3 : raise RuntimeError, "{0} has a space in the formula!  Please check!".format(f[1])
                    else: self.mekd_sig_a0_shape = f[2]
                if f[1].lower().startswith("mekd_sig_a1"): 
                    if len(f) > 3 : raise RuntimeError, "{0} has a space in the formula!  Please check!".format(f[1])
                    else: self.mekd_sig_a1_shape = f[2]
                if f[1].lower().startswith("mekd_sig_a2"): 
                    if len(f) > 3 : raise RuntimeError, "{0} has a space in the formula!  Please check!".format(f[1])
                    else: self.mekd_sig_a2_shape = f[2]
                if f[1].lower().startswith("mekd_sig_a3"): 
                    if len(f) > 3 : raise RuntimeError, "{0} has a space in the formula!  Please check!".format(f[1])
                    else: self.mekd_sig_a3_shape = f[2]
                if f[1].lower().startswith("mekd_sig_a4"): 
                    if len(f) > 3 : raise RuntimeError, "{0} has a space in the formula!  Please check!".format(f[1])
                    else: self.mekd_sig_a4_shape = f[2]

            if f[0].lower().startswith("highmasssignalshape") or f[0].lower().startswith("hmsignalshape"):

                if f[1].lower().startswith("n_cb"): 
                    if len(f) > 3 : raise RuntimeError, "{0} has a space in the formula!  Please check!".format(f[1])
                    else: self.n_CB_shape_HM = f[2]
                if f[1].lower().startswith("alpha_cb"): 
                    if len(f) > 3 : raise RuntimeError, "{0} has a space in the formula!  Please check!".format(f[1])
                    else: self.alpha_CB_shape_HM = f[2]
                if f[1].lower().startswith("n2_cb"): 
                    if len(f) > 3 : raise RuntimeError, "{0} has a space in the formula!  Please check!".format(f[1])
                    else: self.n2_CB_shape_HM = f[2]
                if f[1].lower().startswith("alpha2_cb"): 
                    if len(f) > 3 : raise RuntimeError, "{0} has a space in the formula!  Please check!".format(f[1])
                    else: self.alpha2_CB_shape_HM = f[2]
                if f[1].lower().startswith("mean_cb"): 
                    if len(f) > 3 : raise RuntimeError, "{0} has a space in the formula!  Please check!".format(f[1])
                    else: self.mean_CB_shape_HM = f[2]
                if f[1].lower().startswith("sigma_cb"): 
                    if len(f) > 3 : raise RuntimeError, "{0} has a space in the formula!  Please check!".format(f[1])
                    else: self.sigma_CB_shape_HM = f[2]
                if f[1].lower().startswith("gamma_bw"): 
                    if len(f) > 3 : raise RuntimeError, "{0} has a space in the formula!  Please check!".format(f[1])
                    else: self.gamma_BW_shape_HM = f[2]
                    
            if f[0].lower().startswith("signaleff"):

                if f[1].lower().startswith("a1"): self.sigeff_a1 = float(f[2])
                if f[1].lower().startswith("a2"): self.sigeff_a2 = float(f[2])
                if f[1].lower().startswith("a3"): self.sigeff_a3 = float(f[2])
                if f[1].lower().startswith("a4"): self.sigeff_a4 = float(f[2])
                if f[1].lower().startswith("b1"): self.sigeff_b1 = float(f[2])
                if f[1].lower().startswith("b2"): self.sigeff_b2 = float(f[2])
                if f[1].lower().startswith("b3"): self.sigeff_b3 = float(f[2])

            if f[0].lower().startswith("qqzzshape"):

                if f[1].lower().startswith("a0"): self.qqZZshape_a0 = float(f[2])
                if f[1].lower().startswith("a1_") or f[1].lower().startswith("a1 "): self.qqZZshape_a1 = float(f[2])
                if f[1].lower().startswith("a2"): self.qqZZshape_a2 = float(f[2])
                if f[1].lower().startswith("a3"): self.qqZZshape_a3 = float(f[2])
                if f[1].lower().startswith("a4"): self.qqZZshape_a4 = float(f[2])
                if f[1].lower().startswith("a5"): self.qqZZshape_a5 = float(f[2])
                if f[1].lower().startswith("a6"): self.qqZZshape_a6 = float(f[2])
                if f[1].lower().startswith("a7"): self.qqZZshape_a7 = float(f[2])
                if f[1].lower().startswith("a8"): self.qqZZshape_a8 = float(f[2])
                if f[1].lower().startswith("a9"): self.qqZZshape_a9 = float(f[2])
                if f[1].lower().startswith("a10"): self.qqZZshape_a10 = float(f[2])
                if f[1].lower().startswith("a11"): self.qqZZshape_a11 = float(f[2])
                if f[1].lower().startswith("a12"): self.qqZZshape_a12 = float(f[2])
                if f[1].lower().startswith("a13"): self.qqZZshape_a13 = float(f[2])
                
                if f[1].lower().startswith("mekd_qqzz_a0"): 
                    if len(f) > 3 : raise RuntimeError, "{0} has a space in the formula!  Please check!".format(f[1])
                    else: self.mekd_qqZZ_a0_shape = f[2]; print f[2]; print self.mekd_qqZZ_a0_shape
                if f[1].lower().startswith("mekd_qqzz_a1"): 
                    if len(f) > 3 : raise RuntimeError, "{0} has a space in the formula!  Please check!".format(f[1])
                    else: self.mekd_qqZZ_a1_shape = f[2]
                if f[1].lower().startswith("mekd_qqzz_a2"): 
                    if len(f) > 3 : raise RuntimeError, "{0} has a space in the formula!  Please check!".format(f[1])
                    else: self.mekd_qqZZ_a2_shape = f[2]
                if f[1].lower().startswith("mekd_qqzz_a3"): 
                    if len(f) > 3 : raise RuntimeError, "{0} has a space in the formula!  Please check!".format(f[1])
                    else: self.mekd_qqZZ_a3_shape = f[2]
                if f[1].lower().startswith("mekd_qqzz_a4"): 
                    if len(f) > 3 : raise RuntimeError, "{0} has a space in the formula!  Please check!".format(f[1])
                    else: self.mekd_qqZZ_a4_shape = f[2]

            if f[0].lower().startswith("ggzzshape"):

                if f[1].lower().startswith("a0"): self.ggZZshape_a0 = float(f[2])
                if f[1].lower().startswith("a1"): self.ggZZshape_a1 = float(f[2])
                if f[1].lower().startswith("a2"): self.ggZZshape_a2 = float(f[2])
                if f[1].lower().startswith("a3"): self.ggZZshape_a3 = float(f[2])
                if f[1].lower().startswith("a4"): self.ggZZshape_a4 = float(f[2])
                if f[1].lower().startswith("a5"): self.ggZZshape_a5 = float(f[2])
                if f[1].lower().startswith("a6"): self.ggZZshape_a6 = float(f[2])
                if f[1].lower().startswith("a7"): self.ggZZshape_a7 = float(f[2])
                if f[1].lower().startswith("a8"): self.ggZZshape_a8 = float(f[2])
                if f[1].lower().startswith("a9"): self.ggZZshape_a9 = float(f[2])
               
            if f[0].lower().startswith("zjetsshape"):

                if f[1].lower().startswith("mean"):  self.zjetsShape_mean = f[2]
                if f[1].lower().startswith("sigma"): self.zjetsShape_sigma = f[2]
                

            if f[0].lower().startswith("systematic"):
                
                if f[1].lower().startswith("zjet") and f[1].lower().find("kappalow") >= 0 :
                    self.zjetsKappaLow = f[2]
                if f[1].lower().startswith("zjet") and f[1].lower().find("kappahigh") >= 0 :
                    self.zjetsKappaHigh = f[2]
                if f[1].lower().startswith("lumiunc"):
                    self.lumiUnc = f[2]
                if f[1].lower().startswith("muon_full") or f[1].lower().startswith("muonfull"):
                    self.muonFullUnc = f[2]
                    if len(f) > 3:
                        self.muonFullUnc_HM = f[3]
                        self.muonFullCutoff = f[4]
                if f[1].lower().startswith("muon_trig") or f[1].lower().startswith("muontrig"):
                    self.muonTrigUnc = f[2]
                    if len(f) > 3:
                        self.muonTrigUnc_HM = f[3]
                        self.muonTrigCutoff = f[4]
                if f[1].lower().startswith("elec_full") or f[1].lower().startswith("elecfull"):
                    self.elecFullUnc = f[2]
                    if len(f) > 3:
                        self.elecFullUnc_HM = f[3]
                        self.elecFullCutoff = f[4]
                if f[1].lower().startswith("elec_trig") or f[1].lower().startswith("electrig"):
                    self.elecTrigUnc = f[2]
                    if len(f) > 3:
                        self.elecTrigUnc_HM = f[3]
                        self.elecTrigCutoff = f[4]

                if f[1].lower().startswith("param"):
                    if f[2].lower().startswith("cms_zz4l_mean_m_sig"):
                        self.CMS_zz4l_mean_m_sig = f[3]
                    if f[2].lower().startswith("cms_zz4l_sigma_m_sig"):
                        self.CMS_zz4l_sigma_m_sig = f[3]
                    if f[2].lower().startswith("cms_zz4l_mean_e_sig"):
                        self.CMS_zz4l_mean_e_sig = f[3]
                    if f[2].lower().startswith("cms_zz4l_sigma_e_sig"):
                        self.CMS_zz4l_sigma_e_sig = f[3]
                    if f[2].lower().startswith("cms_zz4l_n_sig"):
                        self.CMS_zz4l_n_sig = f[3]
                    if f[2].lower().startswith("cms_zz4l_gamma_sig"):
                        self.CMS_zz4l_gamma_sig = f[3]
                        
                if f[1].lower().startswith("luminosity"):
                    self.useLumiUnc = self.parseBoolString(f[2])
                if f[1].lower().startswith("pdf_gg"):
                    self.usePdf_gg = self.parseBoolString(f[2])
                if f[1].lower().startswith("pdf_qqbar"):
                    self.usePdf_qqbar = self.parseBoolString(f[2])
                if f[1].lower().startswith("qcdscale_ggh"):
                    self.useQCDscale_ggH = self.parseBoolString(f[2])
                if f[1].lower().startswith("qcdscale_qqh"):
                    self.useQCDscale_qqH = self.parseBoolString(f[2])
                if f[1].lower().startswith("qcdscale_vh"):
                    self.useQCDscale_VH = self.parseBoolString(f[2])
                if f[1].lower().startswith("qcdscale_tth"):
                    self.useQCDscale_ttH = self.parseBoolString(f[2])
                if f[1].lower().startswith("pdf_hzz4l_accept"):
                    self.usePdf_hzz4l_accept = self.parseBoolString(f[2])
                if f[1].lower().startswith("theoryuncxs_highmh"):
                    self.useTheoryUncXS_HighMH = self.parseBoolString(f[2])
                if f[1].lower().startswith("qcdscale_ggvv"):
                    self.useQCDscale_ggVV = self.parseBoolString(f[2])
                if f[1].lower().startswith("qcdscale_vv"):
                    self.useQCDscale_VV = self.parseBoolString(f[2])
                if f[1].lower().startswith("brhiggs_hzz4l"):
                    self.useBRhiggs_hzz4l = self.parseBoolString(f[2])
                if f[1].lower().startswith("cms_eff"):
                    self.useCMS_eff = self.parseBoolString(f[2])
                if f[1].lower().startswith("cms_hzz4l_zjets"):
                    self.useCMS_hzz4l_Zjets = self.parseBoolString(f[2])
                if f[1].lower().startswith("cms_zz4l_bkgmela"):
                    self.useCMS_zz4l_bkgMELA = self.parseBoolString(f[2])
                if f[1].lower().startswith("cms_zz4l_sigmela"):
                    self.useCMS_zz4l_sigMELA = self.parseBoolString(f[2])
                if f[1].lower().startswith("cms_zz4l_mean"):
                    self.useCMS_zz4l_mean = self.parseBoolString(f[2])
                if f[1].lower().startswith("cms_zz4l_sigma"):
                    self.useCMS_zz4l_sigma = self.parseBoolString(f[2])
                if f[1].lower().startswith("cms_zz4l_n"):
                    self.useCMS_zz4l_n = self.parseBoolString(f[2])
                if f[1].lower().startswith("cms_zz4l_gamma"):
                    self.useCMS_zz4l_gamma = self.parseBoolString(f[2])
                
                    
            if f[0].lower().startswith("lumi"):
                self.lumi = float(f[1])

            if f[0].lower().startswith("sqrts"):
                self.sqrts = float(f[1])
                
            if f[0].lower().startswith("dohyptest"):
                self.doHypTest = self.parseBoolString(f[1])
            if f[0].lower().startswith("altHypLabel"):
                self.altHypLabel = f[1]

    def getInputs(self):

        dict = {}

        ## check settings ##
        if self.all_chan and ( self.qqH_chan or self.ggH_chan or self.WH_chan or self.ZH_chan or self.ttH_chan):
            raise RuntimeError, "You cannot request to execute ALL signal channels and single channels at the same time. Check inputs!"
        if not self.goodEntry(self.sqrts): raise RuntimeError, "{0} is not set.  Check inputs!".format("sqrts")

        if self.qqZZ_chan and not self.goodEntry(self.qqZZ_rate): raise RuntimeError, "{0} is not set.  Check inputs!".format("qqZZ_rate")
        if self.ggZZ_chan and not self.goodEntry(self.ggZZ_rate): raise RuntimeError, "{0} is not set.  Check inputs!".format("ggZZ_rate")
        if self.zjets_chan and not self.goodEntry(self.zjets_rate): raise RuntimeError, "{0} is not set.  Check inputs!".format("zjets_rate")
        if self.zbb_chan and not self.goodEntry(self.zbb_rate): raise RuntimeError, "{0} is not set.  Check inputs!".format("zbb_rate")
        if self.ttbar_chan and not self.goodEntry(self.ttbar_rate): raise RuntimeError, "{0} is not set.  Check inputs!".format("ttbar_rate")

        if not self.goodEntry(self.n_CB_shape): raise RuntimeError, "{0} is not set.  Check inputs!".format("n_CB_shape")
        if not self.goodEntry(self.alpha_CB_shape): raise RuntimeError, "{0} is not set.  Check inputs!".format("alpha_CB_shape")
        if not self.goodEntry(self.n2_CB_shape): raise RuntimeError, "{0} is not set.  Check inputs!".format("n2_CB_shape")
        if not self.goodEntry(self.alpha2_CB_shape): raise RuntimeError, "{0} is not set.  Check inputs!".format("alpha2_CB_shape")
        if not self.goodEntry(self.mean_CB_shape): raise RuntimeError, "{0} is not set.  Check inputs!".format("mean_CB_shape")
        if not self.goodEntry(self.sigma_CB_shape): raise RuntimeError, "{0} is not set.  Check inputs!".format("sigma_CB_shape")

        
        if not self.goodEntry(self.n_CB_shape_HM):
            print "{0} is not set. Using {1} for {0}.".format("n_CB_shape_HM","n_CB_shape")
            self.n_CB_shape_HM = self.n_CB_shape
        if not self.goodEntry(self.alpha_CB_shape_HM):
            print "{0} is not set. Using {1} for {0}.".format("alpha_CB_shape_HM","alpha_CB_shape")
            self.alpha_CB_shape_HM = self.alpha_CB_shape
        if not self.goodEntry(self.n2_CB_shape_HM):
            print "{0} is not set. Using {1} for {0}.".format("n2_CB_shape_HM","n2_CB_shape")
            self.n2_CB_shape_HM = self.n2_CB_shape
        if not self.goodEntry(self.alpha2_CB_shape_HM):
            print "{0} is not set. Using {1} for {0}.".format("alpha2_CB_shape_HM","alpha2_CB_shape")
            self.alpha2_CB_shape_HM = self.alpha2_CB_shape
        if not self.goodEntry(self.n_CB_shape_HM):
            print "{0} is not set. Using {1} for {0}.".format("n_CB_shape_HM","n_CB_shape")
            self.n_CB_shape_HM = self.n_CB_shape
        if not self.goodEntry(self.alpha_CB_shape_HM):
            print "{0} is not set. Using {1} for {0}.".format("alpha_CB_shape_HM","alpha_CB_shape")
            self.alpha_CB_shape_HM = self.alpha_CB_shape
        if not self.goodEntry(self.mean_CB_shape_HM):
            print "{0} is not set. Using {1} for {0}.".format("mean_CB_shape_HM","mean_CB_shape")
            self.mean_CB_shape_HM = self.mean_CB_shape
        if not self.goodEntry(self.sigma_CB_shape_HM):
            print "{0} is not set. Using {1} for {0}.".format("sigma_CB_shape_HM","sigma_CB_shape")
            self.sigma_CB_shape_HM = self.sigma_CB_shape
                        
        if not self.goodEntry(self.gamma_BW_shape_HM): raise RuntimeError, "{0} is not set.  Check inputs!".format("gamma_BW_shape_HM")

        if not self.goodEntry(self.sigeff_a1): raise RuntimeError, "{0} is not set.  Check inputs!".format("sigEff_a1")
        if not self.goodEntry(self.sigeff_a2): raise RuntimeError, "{0} is not set.  Check inputs!".format("sigEff_a2")
        if not self.goodEntry(self.sigeff_a3): raise RuntimeError, "{0} is not set.  Check inputs!".format("sigEff_a3")
        if not self.goodEntry(self.sigeff_a4): raise RuntimeError, "{0} is not set.  Check inputs!".format("sigEff_a4")
        if not self.goodEntry(self.sigeff_b1): raise RuntimeError, "{0} is not set.  Check inputs!".format("sigEff_b1")
        if not self.goodEntry(self.sigeff_b2): raise RuntimeError, "{0} is not set.  Check inputs!".format("sigEff_b2")
        if not self.goodEntry(self.sigeff_b3): raise RuntimeError, "{0} is not set.  Check inputs!".format("sigEff_b3")

        if not self.goodEntry(self.qqZZshape_a0): raise RuntimeError, "{0} is not set.  Check inputs!".format("qqZZshape_a0")
        if not self.goodEntry(self.qqZZshape_a1): raise RuntimeError, "{0} is not set.  Check inputs!".format("qqZZshape_a1")
        if not self.goodEntry(self.qqZZshape_a2): raise RuntimeError, "{0} is not set.  Check inputs!".format("qqZZshape_a2")
        if not self.goodEntry(self.qqZZshape_a3): raise RuntimeError, "{0} is not set.  Check inputs!".format("qqZZshape_a3")
        if not self.goodEntry(self.qqZZshape_a4): raise RuntimeError, "{0} is not set.  Check inputs!".format("qqZZshape_a4")
        if not self.goodEntry(self.qqZZshape_a5): raise RuntimeError, "{0} is not set.  Check inputs!".format("qqZZshape_a5")
        if not self.goodEntry(self.qqZZshape_a6): raise RuntimeError, "{0} is not set.  Check inputs!".format("qqZZshape_a6")
        if not self.goodEntry(self.qqZZshape_a7): raise RuntimeError, "{0} is not set.  Check inputs!".format("qqZZshape_a7")
        if not self.goodEntry(self.qqZZshape_a8): raise RuntimeError, "{0} is not set.  Check inputs!".format("qqZZshape_a8")
        if not self.goodEntry(self.qqZZshape_a9): raise RuntimeError, "{0} is not set.  Check inputs!".format("qqZZshape_a9")
        if not self.goodEntry(self.qqZZshape_a10): raise RuntimeError, "{0} is not set.  Check inputs!".format("qqZZshape_a10")
        if not self.goodEntry(self.qqZZshape_a11): raise RuntimeError, "{0} is not set.  Check inputs!".format("qqZZshape_a11")
        if not self.goodEntry(self.qqZZshape_a12): raise RuntimeError, "{0} is not set.  Check inputs!".format("qqZZshape_a12")
        if not self.goodEntry(self.qqZZshape_a13): raise RuntimeError, "{0} is not set.  Check inputs!".format("qqZZshape_a13")

        if not self.goodEntry(self.ggZZshape_a0): raise RuntimeError, "{0} is not set.  Check inputs!".format("ggZZshape_a0")
        if not self.goodEntry(self.ggZZshape_a1): raise RuntimeError, "{0} is not set.  Check inputs!".format("ggZZshape_a1")
        if not self.goodEntry(self.ggZZshape_a2): raise RuntimeError, "{0} is not set.  Check inputs!".format("ggZZshape_a2")
        if not self.goodEntry(self.ggZZshape_a3): raise RuntimeError, "{0} is not set.  Check inputs!".format("ggZZshape_a3")
        if not self.goodEntry(self.ggZZshape_a4): raise RuntimeError, "{0} is not set.  Check inputs!".format("ggZZshape_a4")
        if not self.goodEntry(self.ggZZshape_a5): raise RuntimeError, "{0} is not set.  Check inputs!".format("ggZZshape_a5")
        if not self.goodEntry(self.ggZZshape_a6): raise RuntimeError, "{0} is not set.  Check inputs!".format("ggZZshape_a6")
        if not self.goodEntry(self.ggZZshape_a7): raise RuntimeError, "{0} is not set.  Check inputs!".format("ggZZshape_a7")
        if not self.goodEntry(self.ggZZshape_a8): raise RuntimeError, "{0} is not set.  Check inputs!".format("ggZZshape_a8")
        if not self.goodEntry(self.ggZZshape_a9): raise RuntimeError, "{0} is not set.  Check inputs!".format("ggZZshape_a9")

        if not self.goodEntry(self.zjetsShape_mean): raise RuntimeError, "{0} is not set.  Check inputs!".format("zjetsShape_mean")
        if not self.goodEntry(self.zjetsShape_sigma): raise RuntimeError, "{0} is not set.  Check inputs!".format("zjetsShape_sigma")
        
        if not self.goodEntry(self.zjetsKappaLow): raise RuntimeError, "{0} is not set.  Check inputs!".format("self.zjetsKappaLow")
        if not self.goodEntry(self.zjetsKappaHigh): raise RuntimeError, "{0} is not set.  Check inputs!".format("self.zjetsKappaHigh")

        if not self.goodEntry(self.qqZZ_lumi):  self.qqZZ_lumi = self.lumi
        if not self.goodEntry(self.ggZZ_lumi):  self.ggZZ_lumi = self.lumi
        if not self.goodEntry(self.zjets_lumi): self.zjets_lumi = self.lumi
        if not self.goodEntry(self.zbb_lumi):   self.zbb_lumi = self.lumi
        if not self.goodEntry(self.ttbar_lumi): self.ttbar_lumi = self.lumi

        if not self.goodEntry(self.lumiUnc) and self.useLumiUnc:
            raise RuntimeError, "{0} is not set. Check systematic inputs!".format("Lumi Uncertiatiny")
        if not self.goodEntry(self.muonFullUnc) and self.useCMS_eff and (self.decayChan == 1 or self.decayChan == 3):
            raise RuntimeError,"{0} is not set. Check systematic inputs!".format("muon_full")
        if not self.goodEntry(self.muonTrigUnc) and self.useCMS_eff and (self.decayChan == 1 or self.decayChan == 3):
            raise RuntimeError,"{0} is not set. Check systematic inputs!".format("muon_trig")
        if not self.goodEntry(self.elecFullUnc) and self.useCMS_eff and (self.decayChan == 2 or self.decayChan == 3):
            raise RuntimeError,"{0} is not set. Check systematic inputs!".format("elec_full")
        if not self.goodEntry(self.elecTrigUnc) and self.useCMS_eff and (self.decayChan == 2 or self.decayChan == 3):
            raise RuntimeError,"{0} is not set. Check systematic inputs!".format("elec_trig")

        if not self.goodEntry(self.CMS_zz4l_mean_m_sig) and self.useCMS_zz4l_mean and (self.decayChan == 1 or self.decayChan == 3):
            raise RuntimeError,"{0} is not set. Check systematic inputs!".format("CMS_zz4l_mean_m_sig")
        if not self.goodEntry(self.CMS_zz4l_mean_e_sig) and self.useCMS_zz4l_mean and (self.decayChan == 2 or self.decayChan == 3):
            raise RuntimeError,"{0} is not set. Check systematic inputs!".format("CMS_zz4l_mean_e_sig")
        if not self.goodEntry(self.CMS_zz4l_sigma_m_sig) and self.useCMS_zz4l_sigma and (self.decayChan == 1 or self.decayChan == 3):
            raise RuntimeError,"{0} is not set. Check systematic inputs!".format("CMS_zz4l_sigma_m_sig")
        if not self.goodEntry(self.CMS_zz4l_sigma_e_sig) and self.useCMS_zz4l_sigma and (self.decayChan == 2 or self.decayChan == 3):
            raise RuntimeError,"{0} is not set. Check systematic inputs!".format("CMS_zz4l_sigma_e_sig")
        if not self.goodEntry(self.CMS_zz4l_n_sig) and self.useCMS_zz4l_n:
            raise RuntimeError,"{0} is not set. Check systematic inputs!".format("CMS_zz4l_n_sig")
        if not self.goodEntry(self.CMS_zz4l_gamma_sig) and self.useCMS_zz4l_gamma:
            raise RuntimeError,"{0} is not set. Check systematic inputs!".format("CMS_zz4l_gamma_sig")

        if self.doHypTest:
            print "!!! HYPTOTHESIS TESTING !!!"
  
        if self.doHypTest and not self.all_chan:
            raise RuntimeError,"You asked to prepare DC and WS for Hyp Test but you did not want to sum over all signal channels. This is forbidden. Check inputs !"
      ## Set dictionary entries to be passed to datacard class ##
        
        dict['decayChannel'] = int(self.decayChan)
        dict['model'] = str(self.model)
        dict['lumi'] = float(self.lumi)
        dict['sqrts'] = float(self.sqrts)

        dict['all'] = self.all_chan
        dict['ggH'] = self.ggH_chan
        dict['qqH'] = self.qqH_chan
        dict['WH'] = self.WH_chan
        dict['ZH'] = self.ZH_chan
        dict['ttH'] = self.ttH_chan
        dict['qqZZ'] = self.qqZZ_chan
        dict['ggZZ'] = self.ggZZ_chan
        dict['zjets'] = self.zjets_chan
        dict['ttbar'] = self.ttbar_chan
        dict['zbb'] = self.zbb_chan
       
        dict['qqZZ_rate'] = self.qqZZ_rate
        dict['ggZZ_rate'] = self.ggZZ_rate
        dict['zjets_rate'] = self.zjets_rate
        dict['ttbar_rate'] = self.ttbar_rate
        dict['zbb_rate'] = self.zbb_rate

        dict['qqZZ_lumi'] = float(self.qqZZ_lumi)
        dict['ggZZ_lumi'] = float(self.ggZZ_lumi) 
        dict['zjets_lumi'] = float(self.zjets_lumi)
        dict['ttbar_lumi'] = float(self.ttbar_lumi)
        dict['zbb_lumi'] = float(self.zbb_lumi)

        dict['n_CB_shape'] = self.n_CB_shape
        dict['alpha_CB_shape'] = self.alpha_CB_shape
        dict['n2_CB_shape'] = self.n2_CB_shape
        dict['alpha2_CB_shape'] = self.alpha2_CB_shape
        dict['mean_CB_shape'] = self.mean_CB_shape
        dict['sigma_CB_shape'] = self.sigma_CB_shape

        dict['useHighMassReweightedShapes'] = self.useHighMassReweightedShapes
        dict['n_CB_shape_HM'] = self.n_CB_shape_HM
        dict['alpha_CB_shape_HM'] = self.alpha_CB_shape_HM
        dict['n2_CB_shape_HM'] = self.n2_CB_shape_HM
        dict['alpha2_CB_shape_HM'] = self.alpha2_CB_shape_HM
        dict['mean_CB_shape_HM'] = self.mean_CB_shape_HM
        dict['sigma_CB_shape_HM'] = self.sigma_CB_shape_HM
        dict['gamma_BW_shape_HM'] = self.gamma_BW_shape_HM

        dict['sigEff_a1'] = float(self.sigeff_a1)
        dict['sigEff_a2'] = float(self.sigeff_a2)
        dict['sigEff_a3'] = float(self.sigeff_a3)
        dict['sigEff_a4'] = float(self.sigeff_a4)
        dict['sigEff_b1'] = float(self.sigeff_b1)
        dict['sigEff_b2'] = float(self.sigeff_b2)
        dict['sigEff_b3'] = float(self.sigeff_b3)

        dict['qqZZshape_a0'] = float(self.qqZZshape_a0)
        dict['qqZZshape_a1'] = float(self.qqZZshape_a1)
        dict['qqZZshape_a2'] = float(self.qqZZshape_a2)
        dict['qqZZshape_a3'] = float(self.qqZZshape_a3)
        dict['qqZZshape_a4'] = float(self.qqZZshape_a4)
        dict['qqZZshape_a5'] = float(self.qqZZshape_a5)
        dict['qqZZshape_a6'] = float(self.qqZZshape_a6)
        dict['qqZZshape_a7'] = float(self.qqZZshape_a7)
        dict['qqZZshape_a8'] = float(self.qqZZshape_a8)
        dict['qqZZshape_a9'] = float(self.qqZZshape_a9)
        dict['qqZZshape_a10'] = float(self.qqZZshape_a10)
        dict['qqZZshape_a11'] = float(self.qqZZshape_a11)
        dict['qqZZshape_a12'] = float(self.qqZZshape_a12)
        dict['qqZZshape_a13'] = float(self.qqZZshape_a13)

        dict['ggZZshape_a0'] = float(self.ggZZshape_a0)
        dict['ggZZshape_a1'] = float(self.ggZZshape_a1)
        dict['ggZZshape_a2'] = float(self.ggZZshape_a2)
        dict['ggZZshape_a3'] = float(self.ggZZshape_a3)
        dict['ggZZshape_a4'] = float(self.ggZZshape_a4)
        dict['ggZZshape_a5'] = float(self.ggZZshape_a5)
        dict['ggZZshape_a6'] = float(self.ggZZshape_a6)
        dict['ggZZshape_a7'] = float(self.ggZZshape_a7)
        dict['ggZZshape_a8'] = float(self.ggZZshape_a8)
        dict['ggZZshape_a9'] = float(self.ggZZshape_a9)

        dict['zjetsShape_mean'] = float(self.zjetsShape_mean)
        dict['zjetsShape_sigma'] = float(self.zjetsShape_sigma)

        dict['zjetsKappaLow'] = float(self.zjetsKappaLow)
        dict['zjetsKappaHigh'] = float(self.zjetsKappaHigh)
        

        dict['lumiUnc'] = self.lumiUnc
        dict['muonFullUnc'] = float(self.muonFullUnc)
        dict['muonFullUnc_HM'] = float(self.muonFullUnc_HM) 
        dict['muonFullCutoff'] = float(self.muonFullCutoff)
        dict['elecFullUnc'] = float(self.elecFullUnc) 
        dict['elecFullUnc_HM'] = float(self.elecFullUnc_HM) 
        dict['elecFullCutoff'] = float(self.elecFullCutoff)
        
        dict['muonTrigUnc'] = float(self.muonTrigUnc)
        dict['muonTrigUnc_HM'] = float(self.muonTrigUnc_HM) 
        dict['muonTrigCutoff'] = float(self.muonTrigCutoff)
        dict['elecTrigUnc'] = float(self.elecTrigUnc) 
        dict['elecTrigUnc_HM'] = float(self.elecTrigUnc_HM)
        dict['elecTrigCutoff'] = float(self.elecTrigCutoff) 

        dict['useLumiUnc'] = self.useLumiUnc 
        dict['usePdf_gg'] = self.usePdf_gg
        dict['usePdf_qqbar'] = self.usePdf_qqbar 
        dict['usePdf_hzz4l_accept'] = self.usePdf_hzz4l_accept
        dict['useQCDscale_ggH'] = self.useQCDscale_ggH
        dict['useQCDscale_qqH'] = self.useQCDscale_qqH 
        dict['useQCDscale_VH'] = self.useQCDscale_VH
        dict['useQCDscale_ttH'] = self.useQCDscale_ttH
        dict['useTheoryUncXS_HighMH'] = self.useTheoryUncXS_HighMH
        dict['useQCDscale_ggVV'] = self.useQCDscale_ggVV
        dict['useQCDscale_VV'] = self.useQCDscale_VV
        dict['useBRhiggs_hzz4l'] = self.useBRhiggs_hzz4l
        dict['useCMS_eff'] = self.useCMS_eff
        dict['useCMS_hzz4l_Zjets'] = self.useCMS_hzz4l_Zjets 
        dict['useCMS_zz4l_bkgMELA'] = self.useCMS_zz4l_bkgMELA 
        dict['useCMS_zz4l_sigMELA'] = self.useCMS_zz4l_sigMELA
        dict['useCMS_zz4l_mean'] = self.useCMS_zz4l_mean
        dict['useCMS_zz4l_sigma'] = self.useCMS_zz4l_sigma 
        dict['useCMS_zz4l_n'] = self.useCMS_zz4l_n
        dict['useCMS_zz4l_gamma'] = self.useCMS_zz4l_gamma

        dict['CMS_zz4l_mean_m_sig'] = float(self.CMS_zz4l_mean_m_sig) 
        dict['CMS_zz4l_sigma_m_sig'] = float(self.CMS_zz4l_sigma_m_sig) 
        dict['CMS_zz4l_mean_e_sig'] = float(self.CMS_zz4l_mean_e_sig) 
        dict['CMS_zz4l_sigma_e_sig'] = float(self.CMS_zz4l_sigma_e_sig)
        dict['CMS_zz4l_n_sig'] = float(self.CMS_zz4l_n_sig)
        dict['CMS_zz4l_gamma_sig'] = float(self.CMS_zz4l_gamma_sig)

        dict['doHypTest'] = self.doHypTest
        dict['altHypLabel'] = str(self.altHypLabel)
        
	dict['mekd_sig_a0_shape'] = self.mekd_sig_a0_shape
	dict['mekd_sig_a1_shape'] = self.mekd_sig_a1_shape
	dict['mekd_sig_a2_shape'] = self.mekd_sig_a2_shape
	dict['mekd_sig_a3_shape'] = self.mekd_sig_a3_shape
	dict['mekd_sig_a4_shape'] = self.mekd_sig_a4_shape
	dict['mekd_qqZZ_a0_shape'] = self.mekd_qqZZ_a0_shape
	dict['mekd_qqZZ_a1_shape'] = self.mekd_qqZZ_a1_shape
	dict['mekd_qqZZ_a2_shape'] = self.mekd_qqZZ_a2_shape
	dict['mekd_qqZZ_a3_shape'] = self.mekd_qqZZ_a3_shape
	dict['mekd_qqZZ_a4_shape'] = self.mekd_qqZZ_a4_shape
        return dict
