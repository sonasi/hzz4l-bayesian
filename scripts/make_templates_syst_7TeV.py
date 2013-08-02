#!/usr/bin/env python

#------------------------------------------------------------------------------
# File: make_templates_ntuples.py
# Created: 24-April-2013 Joe Bochenek
# $Revision$
# Description: Apply all the corrections and reproduce ntuples with correct 
#              weight.  The combined weight is stored in f_weight.
#               Also I make input histograms
#------------------------------------------------------------------------------



from ROOT import *
import numpy as np
#import pyBNN
import sys
from time import sleep
import array, os, sys, re
from string import atoi
import math
import plot_util

# Read seed from command line
# If the seed number is 0 then no systematics are applied
argv =  sys.argv[1:]
if len(argv)==0:
	JOBNUM= 0
	SEED = 1234
#	BINS = 1000
else:
	JOBNUM = atoi(argv[0])
	SEED = atoi(argv[1])
#	BINS = atoi(argv[2])

print "\nJOB:  %d\n" % JOBNUM
print "\nSEED: %d\n" % SEED
#print "\nBINS: %d\n" % BINS


gROOT.ProcessLine(".L pyBNN.cc+")
gROOT.ProcessLine(".L Systematics.cc+")
gROOT.ProcessLine(".include CreateDatacards/include")
gROOT.ProcessLine(".L CreateDatacards/include/HiggsCSandWidthFermi.cc+")
gROOT.ProcessLine(".L CreateDatacards/include/HiggsCSandWidthSM4.cc+")
gROOT.ProcessLine(".L CreateDatacards/include/HiggsCSandWidth.cc+")

os.system("ls *.so; ls CreateDatacards/include/*.so")



#Data overlap list
overlap = []


# Create vector of variable names
varnames = vector("string")()

varlist = ["lept1_pt", "lept1_eta", "lept1_phi",
	   "lept2_pt", "lept2_eta", "lept2_phi",
	   "lept3_pt", "lept3_eta", "lept3_phi",
	   "lept4_pt", "lept4_eta", "lept4_phi",
	   "mbestH", "mZ", "mZstar"]

for name in varlist:
	varnames.push_back(name)

# Create a dictionary so that you can more easily obtain these
# values later on
dict = {}
for ind, thisvar in enumerate(varlist):
	dict[thisvar] = ind

# Create instance of Systematics class,
# using the index of the loop as the random seed
syst = Systematics(SEED, varnames)


ROOT.gSystem.AddIncludePath("-I $ROOFITSYS/include/")
ROOT.gSystem.AddIncludePath("-I CreateDatacards/include/")
#ROOT.gSystem.Load("libRooFit")
ROOT.gSystem.Load("CreateDatacards/include/HiggsCSandWidth_cc.so")
ROOT.gSystem.Load("CreateDatacards/include/HiggsCSandWidthSM4_cc.so")

#Declare BNN
myCSW = HiggsCSandWidth("CreateDatacards/include/txtFiles")
t = BayesianNN()

print myCSW.HiggsCS(1, 125, 8)


#Declare VBF MLP
mvavarset = [
"f_massjj",
"f_deltajj",
]
#Declare MLP
reader1 = TMVA.Reader()

var1_ = []
for i, var1 in enumerate(mvavarset):
	var1_.append(array.array('f',[0]))

	print "***** var1: %s" % var1

	reader1.AddVariable(var1,var1_[i])
reader1.BookMVA("MLP","weights/TMVAClassification6_cat2_MLP.weights.xml")



def checkoverlap(thisevent, overlaplist):
    bool_overlap = 0
    for ievent in overlaplist:
        if ievent == thisevent:
            bool_overlap = 1
    return bool_overlap


def makeplot3d(files, name, xmin, ymin, zmin, xmax, ymax, zmax, xbins, ybins, zbins, out_file, bool_reweight, global_reweight, syst, ftype, fnum):
    t.setChannel(channel)    
    syst = Systematics(SEED, varnames)


    #Store 1D Distributiosn For Displaying and Checking the effect of the systematics uncertainties 
    h_m4l = TH1F("plot0_m4l_"+name, "plot0_m4l_"+name, xbins, xmin, xmax)
    h_bnn = TH1F("plot0_bnn_"+name, "plot0_bnn_"+name, ybins, ymin, ymax)
    h_vbf1 = TH1F("vbf_h0_"+name,    "vbf_h0_"+name,    zbins, zmin, zmax) # < 2 jets (pt4l)
    h_vbf2 = TH1F("vbf_h0_"+name,    "vbf_h0_"+name,    zbins, zmin, zmax) # > 1 jet (bnn)

    bnn3d_cat1 = TH3F(name+"_cat1_3d", name+"_cat1_3d", xbins, xmin, xmax, ybins, ymin, ymax, zbins, zmin, zmax)
    bnn2d_cat1 = TH2F(name+"_cat1_2d", name+"_cat1_2d", xbins, xmin, xmax, ybins, ymin, ymax)
    bnn2dvbf_cat1 = TH2F(name+"_cat1_3dvbf", name+"_cat1_3dvbf", ybins, ymin,ymax, zbins, zmin, zmax)

    bnn3d_cat2 = TH3F(name+"cat2_3d", name+"_cat2_3d", xbins, xmin, xmax, ybins, ymin, ymax, zbins, zmin, zmax)
    bnn2d_cat2 = TH2F(name+"_cat2_2d", name+"_cat2_2d", xbins, xmin, xmax, ybins, ymin, ymax)
    bnn2dvbf_cat2 = TH2F(name+"_cat2_3dvbf", name+"_cat2_3dvbf", ybins, ymin,ymax, zbins, zmin, zmax)

    bnn2d = TH2F(name+"_2d", name+"_2d", xbins, xmin, xmax, ybins, ymin, ymax)
    bnn3d = TH3F(name+"_3d", name+"_3d", xbins, xmin, xmax, ybins, ymin, ymax, 10, 0, 0.04)

    # Ok now book the KDTree histograms
    fout = TFile("kdtree.root", "update")
    fout.cd()
    nbins1d = 500
    t1_3dto1d_cat1 = TH1F("t1_3dto1d_cat1_{0}".format(channel), "t1_3dto1d_cat1_{0}".format(channel), nbins1d, 0, 1)
    t1_3dto1d_cat2 = TH1F("t1_3dto1d_cat2_{0}".format(channel), "t1_3dto1d_cat2_{0}".format(channel), nbins1d, 0, 1)
    t1_3dto1d_mass = TH1F("t1_3dto1d_mass_{0}".format(channel), "t1_3dto1d_mass_{0}".format(channel), nbins1d, 0, 1)

    kdtree_mass = gDirectory.Get("thiskdtree_massErr_{0}".format(channel))
    kdtree_cat1 = gDirectory.Get("thiskdtree_cat1_{0}".format(channel))
    kdtree_cat2 = gDirectory.Get("thiskdtree_cat2_{0}".format(channel))



    nbinsX=21
    nbinsYps=25
    nbinsYgrav=29
    
    binsX = [0.000, 0.030, 0.060, 0.100, 0.200, 0.300, 0.400, 0.500, 0.550, 0.600, 0.633, 0.666, 0.700, 0.733, 0.766, 0.800, 0.833, 0.866, 0.900, 0.933, 0.966, 1.000]
    binsYps = [0.000, 0.100, 0.150, 0.200, 0.233, 0.266, 0.300, 0.333, 0.366, 0.400, 0.433, 0.466, 0.500, 0.533, 0.566, 0.600, 0.633, 0.666, 0.700, 0.733, 0.766, 0.800, 0.850, 0.900, 0.950, 1.000]
    binsYgrav = [0.000, 0.100, 0.150, 0.175 , 0.200, 0.225, 0.250, 0.275, 0.300, 0.325, 0.350, 0.375, 0.400, 0.425 , 0.450, 0.475, 0.500, 0.525, 0.575, 0.600, 0.633, 0.666, 0.700, 0.733 , 0.766, 0.800, 0.850, 0.900, 0.950, 1.000]

    abinsX = array.array('d',binsX)
    abinsYps = array.array('d',binsYps)
    abinsYgrav = array.array('d',binsYgrav)
    
    bnn_jcp_2PM = TH2F(name+"_2PM", name+"_2PM", nbinsX, abinsX, nbinsYgrav, abinsYgrav)
    bnn_jcp_0M  = TH2F(name+"_0M", name+"_0M", nbinsX, abinsX, nbinsYps, abinsYps)
#    bnn_jcp_1m  = TH2F(name+"_2d", name+"_2d", xbins, xmin, xmax, ybins, ymin, ymax)


    bnn3d_cat1.Sumw2()
    bnn2d_cat1.Sumw2()
    bnn2dvbf_cat1.Sumw2()
    bnn3d_cat2.Sumw2()
    bnn2d_cat2.Sumw2()
    bnn2dvbf_cat2.Sumw2()
    bnn2d.Sumw2()
    bnn3d.Sumw2()
    bnn_jcp_2PM.Sumw2()
    bnn_jcp_0M.Sumw2()
    t1_3dto1d_cat1.Sumw2()
    t1_3dto1d_cat2.Sumw2()
    t1_3dto1d_mass.Sumw2()

    
    int_weight_sum = 0
    int_weight_sum_cat1 = 0
    int_weight_sum_cat2 = 0
    count = 0
    cat1_count = 0
    cat2_count = 0
    weight_sum = 0
    weight_sum_nosyst = 0

    display = 0
    if display:
        c1 = TCanvas()
        c1.Divide(2,2)

    for file in files:
        f= TFile(file)


        if bool_reweight == 3:
            tree = f.Get("HZZ4LeptonsAnalysisReduced")
        else:
            tree = f.Get("HZZ4LeptonsAnalysis")

        nentry = tree.GetEntries()
    
        
        # Fill Histograms with tight and loose candidates
        for event in tree:


            # Check for redundant events in DATA
            if (bool_reweight == 3) and (checkoverlap("{0}:{1}:{2}".format(event.f_run, event.f_lumi, event.f_event), overlap)):
                print "OVERLAPPING EVENT (event mass {0}) SKIPPED".format(event.f_mass4l)
                continue
                

            if (bool_reweight == 3):
                overlap.append("{0}:{1}:{2}".format(event.f_run, event.f_lumi, event.f_event))

            for ii in xrange(4):
                varname = 'f_lept%d' % (ii+1)

                for j, postfix in enumerate(["_pt", "_eta", "_phi"]):
                    cmd = "%s = event.%s" % (varname + postfix,
                                 varname + postfix)
                    exec(cmd)

            f_njets_pass = event.f_njets_pass
            f_angle_costhetastar = event.f_angle_costhetastar
            f_angle_costheta1 = event.f_angle_costheta1
            f_angle_costheta2 = event.f_angle_costheta2
            f_angle_phi = event.f_angle_phi
            f_angle_phistar1 = event.f_angle_phistar1
            f_mass4l = event.f_mass4l
            f_pt4l = event.f_pt4l
            f_massjj = event.f_massjj
            f_deltajj = event.f_deltajj
            f_Z1mass = event.f_Z1mass
            f_Z2mass = event.f_Z2mass
        
            # Shift variables here

            varvalues = vector("double")()

            for ii in xrange(4):
                varname = 'event.f_lept%d' % (ii+1)

                for j in ["_pt", "_eta", "_phi"]:
                    cmd = varname + j

                    varvalues.push_back(eval(cmd))

            varvalues.push_back(event.f_mass4l)
            varvalues.push_back(event.f_Z1mass)
            varvalues.push_back(event.f_Z2mass)

            # Weight from "doline" method of Systematics
            sweight = 1.
            
            if SEED>0:
                if not(ftype == "data"):
                    sweight = syst.doline(varvalues, fnum, ftype)
#                print "using systematics!!!"
#            else:
#                    print "Using systematics!!!"
    
            # Update original variables
    #	    print "="*50
            for ii in xrange(4):
                varname = 'f_lept%d' % (ii+1)

                for j, postfix in enumerate(["_pt", "_eta", "_phi"]):
                    cmd = "%s = %s" % (varname + postfix,
                               varvalues[3*ii+j])
    #			    print  cmd
    #			    print "BEFORE %s = %e" % (varname+postfix,eval(varname+postfix))
                    exec(cmd)
    #			    print "AFTER %s = %e" % (varname+postfix,eval(varname+postfix))

            #print "BEFORE f_lept1_pt =", f_lept1_pt
            #setattr(event, 'f_lept1_pt', varvalues[0])
            #f_lept1_pt = varvalues[0]
            #print "AFTER f_lept1_pt =", f_lept1_pt
            #f_lept2_pt = varvalues[3]
            #f_lept3_pt = varvalues[6]
            #f_lept4_pt = varvalues[9]
            
            f_mass4l = varvalues[-3]
            if varvalues[-2] > 0:
                f_Z1mass = varvalues[-2]
            if varvalues[-1] > 0:
                f_Z2mass = varvalues[-1]


            count += 1
            if(event.f_mass4l) < 100:
                continue
            weight = 1.
                
            if math.isnan(sweight):
                sweight = 1.

            if bool_reweight==1:
                weight = event.f_weight * event.f_pu_weight * event.f_eff_weight * event.f_int_weight * sweight
            elif bool_reweight == 2:
                weight = event.f_weight_fakerate * sweight
            elif bool_reweight == 3:
                weight = 1.
            else: 
                weight = event.f_weight * event.f_pu_weight * event.f_eff_weight * sweight
            weight_sum += weight
            weight_sum_nosyst += weight / sweight



            if(event.f_mass4l) > 180:
                continue
                


            #  Fill the various histograms.  Slice the data.
            #      cat 1     cat 2
            #   --------------------
            #   |         |         |  int. reweight
            #   |         |         |
            #   --------------------
            #   |         |         |  no int. reweight
            #   |         |         |
            #   --------------------
            #    
            #   Global reweight = PU * scale * etc
            bnnvalps = t.getBNNvalueSPS(event.f_Z1mass, event.f_Z2mass, event.f_angle_costhetastar, event.f_angle_costheta1, event.f_angle_costheta2, event.f_angle_phi, event.f_angle_phistar1, event.f_mass4l)
            bnnvalgrav = t.getBNNvalueGrav(event.f_Z1mass, event.f_Z2mass, event.f_angle_costhetastar, event.f_angle_costheta1, event.f_angle_costheta2, event.f_angle_phi, event.f_angle_phistar1, event.f_mass4l)

            bnnval = t.getBNNvalue(event.f_Z1mass, event.f_Z2mass, event.f_angle_costhetastar, event.f_angle_costheta1, event.f_angle_costheta2, event.f_angle_phi, event.f_angle_phistar1, event.f_mass4l,  -1, -1)
            mlpval = 0
            bnn2d.Fill(event.f_mass4l, bnnval, weight)
            bnn3d.Fill(event.f_mass4l, bnnval, event.f_mass4lErr/event.f_mass4l, weight)

            h_m4l.Fill(event.f_mass4l, weight)
            h_bnn.Fill(bnnval, weight)

            #re-scale mass variable
            mmin = 100.0
            mmax = 180.0
            mrange= mmax-mmin

            m4l_scale  = (event.f_mass4l-mmin)/mrange

            datapoint = [ m4l_scale, bnnval, event.f_mass4lErr/event.f_mass4l ]
            adatapoint = array.array('d')
            adatapoint.fromlist(datapoint)
#            print "Node: {0} ".format( kdtree_mass.FindNode(adatapoint) - kdtree_mass.GetNNodes())  
            t1_3dto1d_mass.Fill( ( kdtree_mass.FindNode(adatapoint) - kdtree_mass.GetNNodes() + 0.5)/nbins1d , weight)



            if (((event.f_mass4l) > 121.5) & ((event.f_mass4l) < 131.5)):
                bnn_jcp_2PM.Fill(bnnval, bnnvalgrav)
                bnn_jcp_0M.Fill(bnnval, bnnvalps)


            if(event.f_njets_pass > 1):
                var1_[0][0] = float(event.f_massjj)
                var1_[1][0] = float(event.f_deltajj)            
                mlpval = reader1.EvaluateMVA("MLP")
                
                bnn3d_cat2.Fill(event.f_mass4l, bnnval, mlpval, weight)
                bnn2d_cat2.Fill(event.f_mass4l, bnnval, weight)

                bnn2dvbf_cat2.Fill(bnnval, mlpval, weight)
                h_vbf2.Fill(mlpval, weight)
                
                if (((event.f_mass4l) > 121.5) & ((event.f_mass4l) < 131.5)):
                    datapoint = [ m4l_scale, bnnval, mlpval ]
                    adatapoint = array.array('d')
                    adatapoint.fromlist(datapoint)
                    t1_3dto1d_cat2.Fill( ( kdtree_cat2.FindNode(adatapoint) - kdtree_cat2.GetNNodes() + 0.5)/nbins1d , weight)
#                print "{0:.2f};{1};{2};{3};{4};{5};{6}".format(event.f_mass4l, channel,event.f_run, event.f_lumi, event.f_event, mlpval, event.f_njets_pass)
#                print "Bin: {0}".format(( kdtree_cat2.FindNode(adatapoint) - kdtree_cat2.GetNNodes() + 0.5)/nbins1d)
            else:
                bnnval = t.getBNNvalue(event.f_Z1mass, event.f_Z2mass, event.f_angle_costhetastar, event.f_angle_costheta1, event.f_angle_costheta2, event.f_angle_phi, event.f_angle_phistar1, event.f_mass4l,  -1, -1)
                bnn3d_cat1.Fill(event.f_mass4l, bnnval, event.f_pt4l/event.f_mass4l, weight)
                bnn2d_cat1.Fill(event.f_mass4l, bnnval, weight)
                bnn2dvbf_cat1.Fill(bnnval, event.f_pt4l/event.f_mass4l, weight)
                h_vbf1.Fill(event.f_pt4l/event.f_mass4l, weight)

                if (((event.f_mass4l) > 121.5) & ((event.f_mass4l) < 131.5)):
    
                    datapoint = [ m4l_scale, bnnval, event.f_pt4l/event.f_mass4l ]
                    adatapoint = array.array('d')
                    adatapoint.fromlist(datapoint)

                    t1_3dto1d_cat1.Fill( ( kdtree_cat1.FindNode(adatapoint) - kdtree_cat1.GetNNodes() + 0.5)/nbins1d , weight)
#                print "{0:.2f};{1};{2};{3};{4};{5};{6}".format(event.f_mass4l, channel,event.f_run, event.f_lumi, event.f_event, event.f_pt4l/event.f_mass4l, event.f_njets_pass)



            if display:
                if (not(count%1000)):
                    print "{0} - mass: {1}, bnn: {2}, bnngrav: {3}, bnnps: {4}, mlp: {5}, weight: {6}".format(count, event.f_mass4l, bnnval, bnnvalps, bnnvalgrav, mlpval, weight)
                    c1.cd(1)
                    bnn2d_cat1.Draw("colz")
                    c1.cd(2)
                    bnn2d_cat2.Draw("colz")
                    c1.cd(3)
                    bnn_jcp_2PM.Draw("colz")
                    c1.cd(4)
                    bnn_jcp_0M.Draw("colz")
                    c1.Update()
    print "histogram integral before {0},\nweight_sum (> 100) {1}".format(bnn2d.Integral(), weight_sum)


    out_file.cd()

    h_m4l.Write("h_m4l_"+name)
    h_bnn.Write("h_bnn_"+name)
    h_vbf1.Write("h_vbf1_"+name)
    h_vbf2.Write("h_vbf2_"+name)

    bnn3d_cat1.Write(name+"_cat1_3d")
    bnn2d_cat1.Write(name+"_cat1_2d")
    bnn2dvbf_cat1.Write(name+"_cat1_2dvbf")

    bnn3d_cat2.Write(name+"_cat2_3d")
    bnn2d_cat2.Write(name+"_cat2_2d")
    bnn2dvbf_cat2.Write(name+"_cat2_2dvbf")

    bnn2d.Write(name+"_2d")
    bnn3d.Write(name+"_3d")

    bnn_jcp_2PM.Write(name+"_jcp_2PM")
    bnn_jcp_0M.Write(name+"_jcp_0M")

    t1_3dto1d_cat1.Write("{0}_3dto1d_cat1".format(name))
    t1_3dto1d_cat2.Write("{0}_3dto1d_cat2".format(name))
    t1_3dto1d_mass.Write("{0}_3dto1d_mass".format(name))


masses_7TeV_ggH = [115, 120, 125, 125, 126, 130, 140, 150, 160, 170, 180, 190, 200 ]
masses_7TeV_VBF = [115, 120, 125, 125, 126, 130, 140, 150, 160, 170, 180, 190, 200 ]


reweight = [1, 1, 0]
eras = ["7TeV", "8TeV"]

# X axis is m4l
xmin = 115
xmax = 180
xbins = 65

# Y axis is signal/background discriminant
ymin = 0.
zmin = 0.
ymax = 1.

# Z axis is VBF discriminant
zmax = 1.
ybins = 20
zbins = 20



dirout  = "/lustre/cms/store/user/jpb/templates/"

out_file = TFile(dirout + "hzz_templates_7TeV_%4.4d_%8.8d.root" % (JOBNUM,SEED),
		                  "recreate")
dirin =      "/lustre/cms/store/user/defilip/"
dirinzx = "/lustre/cms/store/user/jpb/data_zx_estimate/"
#out_file = TFile("hzz_templates_test.root", "recreate")
limits_dir = "CreateDatacards/"
#dirin = "/home/jbochenek/data/HZZ4l_2013_ntuples/mc/"
#dirout = "/home/jbochenek/data/HZZ4l_2013_ntuples/mc/templates/"

#dirindata = "/home/jbochenek/data/HZZ4l_2013_ntuples/data/"
Z2lBR = 0.03365
lumi7TeV = 5.051
lumi7TeV = 19.6



#We're only doing 8TeV for now
sqrts = 7
era = str(sqrts) + "TeV"




#Data Frist
overlap = []


# Convert string to integer when calling Systematics class
chanmap = {'4mu': 1,
	   '4e': 2,
	   '2e2mu': 3}

      


channels = ["4mu", "2e2mu",  "4e"]


for ichannel, channel in enumerate(channels):
    sqrts = 8
    era = str(sqrts) + "TeV"

    t.setChannel(channel)
    t.setEra(era)

    print ""
    name = "data_{0}_{1}".format(channel,era)
    print name  
    thisfile = '{0}data_{2}_{1}_bnn_new3.root'.format(dirinzx, channel, era)
    print thisfile
    evtcountfile = plot_util.getsumweight_ntuple([thisfile])

    makeplot3d([thisfile], name, xmin, ymin, zmin, xmax, ymax, zmax, xbins, ybins, zbins, out_file, 3, evtcountfile,          syst, 'data', chanmap[channel])


for ichannel, channel in enumerate(channels):
    # Reducible Background ( Z+X  )
    print ""
    name = "zx_{0}_{1}".format(channel,era)
    print name  
    thisfile = '{0}/output_ntuple_ZX_{2}_{1}_bnn.root'.format(dirinzx, channel, "8TeV")
    makeplot3d([thisfile], name, xmin, ymin, zmin, xmax, ymax, zmax, xbins, ybins, zbins, out_file, 2, 1., syst, 'bkg_ZX', chanmap[channel])


# MC Samples Now
for ichannel, channel in enumerate(channels):
    # Background ( qq -> ZZ )
    print ""
    name = "qqzz_{0}_{1}".format(channel,era)
    print name  
    thisfile = '{0}histos{1}_paper/output_ZZTo{1}_{2}-powheg-pythia6_bnn.root'.format(dirin, channel, era)
    datacard = "{0}SM_inputs_{1}/inputs_{2}.txt".format(limits_dir, era, channel) 
    sumweight = plot_util.getsumweight_ntuple([thisfile])
    rate = float(plot_util.getrate(datacard, "qqZZ"))
    print "rate = {0}".format(rate)
    print "sumwegiht = {0}".format(sumweight)
#    print "evt count / evt count file = {0}".format(evtcount / evtcountfile)

    makeplot3d([thisfile], name, xmin, ymin, zmin, xmax, ymax, zmax, xbins, ybins, zbins, out_file, 0, rate, syst, 'bkg_qqZZ', chanmap[channel])



for ichannel, channel in enumerate(channels):
    # signal VBF, qq -> H -> ZZ -> 4l
    for thismass in masses_7TeV_VBF:
        print ""
        name = "qqH_{0}_{1}_{2}".format(thismass, channel, era)
        print name
        dataset = 'VBF_HToZZTo4L_M-{0}_{1}-powheg-pythia6'.format(thismass, era)
        thisfile = '{0}histos{1}_paper/output_{2}_bnn.root'.format(dirin, channel,dataset)
        eff = plot_util.getefficiency(thisfile, "sig_input_Summer12_paper.txt", dataset)
        evtcountfile = plot_util.getsumweight_ntuple([thisfile])
        CS_ggH = myCSW.HiggsCS(2, float(thismass), sqrts)
        BR_HZZ = myCSW.HiggsBR(15, float(thismass))
        evtcount = CS_ggH*BR_HZZ*eff*1000*lumi7TeV

        print "CS_ggH ({0}) * BR_HZZ ({1}) * 9*Z2lBR({3})*Z2lBR({3}) * eff({2}) * (1000) * lumi7TeV({4}) = {5}".format(CS_ggH, BR_HZZ, eff, Z2lBR, lumi7TeV, evtcount)
        print "evt count / evt count file = {0}".format(evtcount / evtcountfile)
        doreweight = reweight[ichannel]
        if thismass > 170:
            doreweight = 0
        makeplot3d([thisfile], name, xmin, ymin, zmin, xmax, ymax, zmax, xbins, ybins, zbins, out_file, doreweight, evtcount, syst, "sig_VBF",  chanmap[channel])

        print ""




for ichannel, channel in enumerate(channels):
    # signal gg -> H -> ZZ -> 4l
    for thismass in masses_7TeV_ggH:
        print ""
        name = "ggH_{0}_{1}_{2}".format(thismass,channel,era)
        print name

        dataset = 'GluGluToHToZZTo4L_M-{0}_{1}-powheg-pythia6'.format(thismass, era)
        thisfile = '{0}histos{2}_paper/output_{1}_bnn.root'.format(dirin, dataset, channel)
        eff = plot_util.getefficiency(thisfile, "sig_input_Summer12_paper.txt", dataset)
        CS_ggH = myCSW.HiggsCS(1, float(thismass), sqrts)
        BR_HZZ = myCSW.HiggsBR(15, float(thismass))
        evtcountfile = plot_util.getsumweight_ntuple([thisfile])
        evtcount = 1*CS_ggH*BR_HZZ*eff*1000*lumi7TeV
        print "CS_ggH ({0}) * BR_HZZ ({1}) * eff({2}) * (1000) * lumi7TeV({4}) = {5}".format(CS_ggH, BR_HZZ, eff, Z2lBR, lumi7TeV, evtcount)
        print "evt count file = {0}".format(evtcount / evtcountfile)
        print ""

        doreweight = reweight[ichannel]

        if thismass > 170:
            doreweight = 0
        makeplot3d([thisfile], name, xmin, ymin, zmin, xmax, ymax, zmax, xbins, ybins, zbins, out_file, doreweight, evtcount, syst, "sig_ggH", chanmap[channel])




for ichannel, channel in enumerate(channels):
    thismass = "126"
    dataset = 'Graviton2PMToZZTo4L_M-{1}_{0}-JHUgen-PYTHIA6_Tauola'.format(era, thismass)
    thisfile = '{0}histos{2}_paper/output_{1}_bnn.root'.format(dirin, dataset, channel)

    name = "higgs2pm_{0}_{1}".format(channel,era)
    print name  
    eff = plot_util.getefficiency(thisfile, "sig_input_Summer12_paper.txt", dataset)
    CS_ggH = myCSW.HiggsCS(1, float(thismass), sqrts)
    BR_HZZ = myCSW.HiggsBR(11, float(thismass))
    evtcountfile = plot_util.getsumweight_ntuple([thisfile])
    evtcount = CS_ggH*BR_HZZ*(9*Z2lBR*Z2lBR)*eff*1000*lumi7TeV
    doreweight = 0
    makeplot3d([thisfile], name, xmin, ymin, zmin, xmax, ymax, zmax, xbins, ybins, zbins, out_file, doreweight, evtcount, syst, "sig_ggH", chanmap[channel])





for ichannel, channel in enumerate(channels):
    thismass = "126"
    dataset = 'Higgs0MToZZTo4L_M-{1}_{0}-JHUgen-pythia6'.format(era, thismass)
    thisfile = '{0}histos{2}_paper/output_{1}_bnn.root'.format(dirin, dataset, channel)
    name = "higgs0m_{0}_{1}".format(channel,era)
    print name  
    eff = plot_util.getefficiency(thisfile, "sig_input_Summer12_paper.txt", dataset)
    CS_ggH = myCSW.HiggsCS(1, float(thismass), sqrts)
    BR_HZZ = myCSW.HiggsBR(11, float(thismass))
    evtcountfile = plot_util.getsumweight_ntuple([thisfile])
    evtcount = 2*CS_ggH*BR_HZZ*(9*Z2lBR*Z2lBR)*eff*1000*lumi7TeV
    print "CS_ggH ({0}) * BR_HZZ ({1}) * 9*Z2lBR({3})*Z2lBR({3}) * eff({2}) * (1000) * lumi7TeV({4}) = {5}".format(CS_ggH, BR_HZZ, eff, Z2lBR, lumi7TeV, evtcount)

    doreweight = 0
    makeplot3d([thisfile], name, xmin, ymin, zmin, xmax, ymax, zmax, xbins, ybins, zbins, out_file, doreweight, evtcount, syst, "sig_ggH", chanmap[channel])



quit()

#  These files don't seem to exist for 7TeV

for ichannel, channel in enumerate(channels):
    #    for channel in channels:
    t.setChannel(channel)
    print "reweight? " + str(reweight[ichannel])

    # Background gg -> ZZ 
    files = [
    "{0}histos{2}_paper/output_GluGluToZZTo2L2L_TuneZ2star_{1}-gg2zz-pythia6_bnn.root".format(dirin, era, channel),
    "{0}histos{2}_paper/output_GluGluToZZTo4L_{1}-gg2zz-pythia6_bnn.root".format(dirin, era, channel)
    ]

    print ""
    name = "ggzz_{0}_{1}".format(channel,era)
    print name
    datacard = "{0}SM_inputs_{1}/inputs_{2}.txt".format(limits_dir, era, channel) 
    rate = float(plot_util.getrate(datacard, "ggZZ"))
    evtcountfile = plot_util.getsumweight_ntuple(files)
    print "rate = {0}".format(rate)
    doreweight = 0
    makeplot3d(files, name, xmin, ymin, zmin, xmax, ymax, zmax, xbins, ybins, zbins, out_file, doreweight, rate, syst, "bkg_ggZZ", chanmap[channel])


