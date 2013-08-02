#!/usr/bin/env python
#------------------------------------------------------------------------------
# File: kdtreebinning.py
# Description: Bin D vs m4l using KDtree binning
# Created: 21-Apr-2013 Harrison B. Prosper
#------------------------------------------------------------------------------
import os, sys, re
from ROOT import *
from histutil import *
from time import sleep
#------------------------------------------------------------------------------
getera = re.compile('[7-8]TeV')
def getEra(record):
    x = getera.findall(record)
    if len(x) > 0:
        return x[0]
    else:
        return None
    
getchannel = re.compile('(?<=8TeV_).+(?=_bnn)')
def getChannel(record):
    x = getchannel.findall(record)
    if len(x) > 0:
        return x[0]
    else:
        return None
        
gROOT.ProcessLine(".L pyBNN.cc+")




#Declare VBF MLP
mvavarset = [
"f_massjj",
"f_deltajj",
]
#Declare MLP
reader1 = TMVA.Reader()

var1_ = []
for i, var1 in enumerate(mvavarset):
	var1_.append(array('f',[0]))

	print "***** var1: %s" % var1

	reader1.AddVariable(var1,var1_[i])
reader1.BookMVA("MLP","/cmshome/bochenek/Higgs_systematics/scripts/weights/TMVAClassification6_cat2_MLP.weights.xml")




#------------------------------------------------------------------------------    
def main():
    print
    argv = sys.argv[1:]	# Get command line arguments
    argc = len(argv)	# Get number of arguments
    if argc < 1:
        print '''
	Usage:
	     ./kdtreebinning.py %s<root-file-name>%s %s...%s [nrows, default all]
        ''' % (BOLDRED, RESETCOLOR, BOLDBLUE, RESETCOLOR)
        sys.exit(0)

    # Get list of files
    filenames = argv
    try:
        nrows = atoi(argv[-1])
        filenames = filenames[:-1]
    except:
        nrows = None
    treename = "HZZ4LeptonsAnalysis"

    # Load BNN class
    bnn = BayesianNN()

    # Loop over files

    D   = []
    m4l = []
    m4lError = []
    mlpval = []
    weight = []
    
    D1   = []
    m4l1 = []
    m4lError1 = []
    mlpval1 = []
    weight1 = []

    D2   = []
    m4l2 = []
    m4lError2 = []
    mlpval2 = []
    weight2 = []

    
    for filename in filenames:    
        ntuple = Ntuple(filename, treename, nrows)

        era = getEra(filename)

        channel = getChannel(filename)
        print channel
        bnn.setEra(era)
        bnn.setChannel(channel)
        print "\tera: %s\tchannel: %s\n" % (era, channel)
        
        for rownumber, event in enumerate(ntuple):
            w = 1.0
            w *= event.f_weight
            w *= event.f_int_weight
            w *= event.f_pu_weight
            w *= event.f_eff_weight
            
            y = bnn.getBNNvalue(event.f_Z1mass,
                                event.f_Z2mass,
                                event.f_angle_costhetastar,
                                event.f_angle_costheta1,
                                event.f_angle_costheta2,
                                event.f_angle_phi,
                                event.f_angle_phistar1,
                                event.f_mass4l)

            var1_[0][0] = float(event.f_massjj)
            var1_[1][0] = float(event.f_deltajj)            
            y2 = reader1.EvaluateMVA("MLP")


            D.append(y)
            m4l.append(event.f_mass4l)
            mlpval.append(y2)
            m4lError.append(event.f_mass4lErr/event.f_mass4l)
            weight.append(w)


            if event.f_mass4l > 131.5:
                continue
            if event.f_mass4l < 121.5:
                continue
                   
            if event.f_njets_pass > 1:
                D2.append(y)
                m4l2.append(event.f_mass4l)
                mlpval2.append(y2)
                weight2.append(w)

            if event.f_njets_pass < 2:
                D1.append(y)
                m4l1.append(event.f_mass4l)
                mlpval1.append(event.f_pt4l/event.f_mass4l)
                weight1.append(w)
            

    #re-scale mass variable
    mmin = 100.0
    mmax = 180.0
    mrange= mmax-mmin
    m4l  = map(lambda x: (x-mmin)/mrange, m4l)
    
    # KDTreeBinning only accepts a dataset that has a size that
    # is an exact multiple of the number of bins
    # For now ignore event weights
    
    nbins = 500
    nvars = 3


    fout = TFile("kdtree.root", "update")
    fout.cd()
        
    # Do the 3D (m4l, D, m3lError)
    t1 = TH1F("th1", "th1", nbins, 0, 1)

    k = len(D) / nbins
    datasetsize = k*nbins
    print 'Datasetsize: %d' % datasetsize
    D   = D[:datasetsize]
    m4l = m4l[:datasetsize]
    m4lError = m4lError[:datasetsize]
    
    data = m4l + D + m4lError
    KDdata = array('d')
    KDdata.fromlist(data)

    kde = TKDTreeBinning(datasetsize, nvars, KDdata, nbins)
    kde.SortBinsByDensity()

    kde.Write( "kdtreebinning_massErr_{0}".format(channel) )
    kdtree = kde.GetTree()
    kdtree.Write("thiskdtree_massErr_{0}".format(channel) )

    # Create the KDTreeBinning object:



    # Do the 3D (m4l, D1, D2)  Cat 1
    t1_1 = TH1F("th1_1", "th1_1", nbins, 0, 1)

    mmin = 121.5
    mmax = 131.5
    mrange = mmax-mmin

    m4l1  = map(lambda x: (x-mmin)/mrange, m4l1)
    m4l2  = map(lambda x: (x-mmin)/mrange, m4l2)

    k1 = len(D1) / nbins
    datasetsize = k1*nbins
    print 'Datasetsize: %d' % datasetsize
    D1   = D1[:datasetsize]
    m4l1 = m4l1[:datasetsize]
    mlpval1 = mlpval1[:datasetsize]
    
    data1 = m4l1 + D1 + mlpval1
    KDdata1 = array('d')
    KDdata1.fromlist(data1)

    kde_1 = TKDTreeBinning(datasetsize, nvars, KDdata1, nbins)
    kde_1.SortBinsByDensity()
    kdtree_1 = kde_1.GetTree()

    kde_1.Write( "kdtreebinning_cat1_{0}".format(channel) )
    kdtree_1.Write("thiskdtree_cat1_{0}".format(channel))



    # Do the 3D (m4l, D1, D2)  Cat 1
    t1_2 = TH1F("th1_2", "th1_2", nbins, 0, 1)

    k2 = len(D2) / nbins
    datasetsize = k2*nbins
    print 'Datasetsize: %d' % datasetsize
    D2   = D2[:datasetsize]
    m4l2 = m4l2[:datasetsize]
    mlpval2 = mlpval2[:datasetsize]
    
    data2 = m4l2 + D2 + mlpval2
    KDdata2 = array('d')
    KDdata2.fromlist(data2)

    kde_2 = TKDTreeBinning(datasetsize, nvars, KDdata2, nbins)
    kde_2.SortBinsByDensity()
    kdtree_2 = kde_2.GetTree()

    kde_2.Write( "kdtreebinning_cat2_{0}".format(channel) )
    kdtree_2.Write( "thiskdtree_cat2_{0}".format(channel) )






    # Do the 2D (m4l, D1)  for comparison
    data = m4l + D
    KDdata3 = array('d')
    KDdata3.fromlist(data)

    kde_3 = TKDTreeBinning(datasetsize, nvars, KDdata3, nbins)
    kde_3.SortBinsByDensity()
    kdtree_3 = kde_3.GetTree()

    kde_3.Write( "kdtreebinning_2d_{0}".format(channel) )
    kdtree_3.Write( "thiskdtree_2d_{0}".format(channel) )





    # Do the 3D (m4l, D1, D2)  Cat 1 different binning
    t1_2 = TH1F("th1_2", "th1_2", nbins, 0, 1)
    nbins = 200
    k2 = len(D2) / nbins
    datasetsize = k2*nbins
    print 'Datasetsize: %d' % datasetsize
    D2   = D2[:datasetsize]
    m4l2 = m4l2[:datasetsize]
    mlpval2 = mlpval2[:datasetsize]
    
    data2 = m4l2 + D2 + mlpval2
    KDdata2 = array('d')
    KDdata2.fromlist(data2)

    kde_2 = TKDTreeBinning(datasetsize, nvars, KDdata2, nbins)
    kde_2.SortBinsByDensity()
    kdtree_2 = kde_2.GetTree()

    kde_2.Write( "kdtreebinning_cat2_200_{0}".format(channel) )
    kdtree_2.Write( "thiskdtree_cat2_200_{0}".format(channel) )


    # Do the 3D (m4l, D1, D2)  Cat 1
    t1_2 = TH1F("th1_2", "th1_2", nbins, 0, 1)

    k2 = len(D2) / nbins
    datasetsize = k2*nbins
    print 'Datasetsize: %d' % datasetsize
    D2   = D2[:datasetsize]
    m4l2 = m4l2[:datasetsize]
    mlpval2 = mlpval2[:datasetsize]
    
    data2 = m4l2 + D2 + mlpval2
    KDdata2 = array('d')
    KDdata2.fromlist(data2)

    kde_2 = TKDTreeBinning(datasetsize, nvars, KDdata2, nbins)
    kde_2.SortBinsByDensity()
    kdtree_2 = kde_2.GetTree()

    kde_2.Write( "kdtreebinning_cat2_200_{0}".format(channel) )
    kdtree_2.Write( "thiskdtree_cat2_200_{0}".format(channel) )


    # Do the 3D (m4l, D1, D2)  Cat 1 different binning
    t1_2 = TH1F("th1_2", "th1_2", nbins, 0, 1)
    nbins = 1000
    k2 = len(D2) / nbins
    datasetsize = k2*nbins
    print 'Datasetsize: %d' % datasetsize
    D2   = D2[:datasetsize]
    m4l2 = m4l2[:datasetsize]
    mlpval2 = mlpval2[:datasetsize]
    
    data2 = m4l2 + D2 + mlpval2
    KDdata2 = array('d')
    KDdata2.fromlist(data2)

    kde_2 = TKDTreeBinning(datasetsize, nvars, KDdata2, nbins)
    kde_2.SortBinsByDensity()
    kdtree_2 = kde_2.GetTree()

    kde_2.Write( "kdtreebinning_cat2_1000_{0}".format(channel) )
    kdtree_2.Write( "thiskdtree_cat2_1000_{0}".format(channel) )



    # Do the 3D (m4l, D1, D2)  Cat 1
    t1_2 = TH1F("th1_2", "th1_2", nbins, 0, 1)

    k2 = len(D2) / nbins
    datasetsize = k2*nbins
    print 'Datasetsize: %d' % datasetsize
    D2   = D2[:datasetsize]
    m4l2 = m4l2[:datasetsize]
    mlpval2 = mlpval2[:datasetsize]
    
    data2 = m4l2 + D2 + mlpval2
    KDdata2 = array('d')
    KDdata2.fromlist(data2)

    kde_2 = TKDTreeBinning(datasetsize, nvars, KDdata2, nbins)
    kde_2.SortBinsByDensity()
    kdtree_2 = kde_2.GetTree()

    kde_2.Write( "kdtreebinning_cat2_1000_{0}".format(channel) )
    kdtree_2.Write( "thiskdtree_cat2_1000_{0}".format(channel) )
    
    
    

    
    # Plot
    h2pol = TH2Poly("h2plot", "KDTree binning", 0, 1, 0, 1)
    fixhist2(h2pol, 'm_{4l}', 'D')
    
    binsMinEdges = kde_2.GetBinsMinEdges()
    binsMaxEdges = kde_2.GetBinsMaxEdges()
    
    for i in xrange(nbins):
        j = i*nvars
        h2pol.AddBin(binsMinEdges[j], binsMinEdges[j + 1],
                     binsMaxEdges[j], binsMaxEdges[j + 1])
        h2pol.SetBinContent(i+1, kde_2.GetBinDensity(i))
        datapoint = [ (binsMaxEdges[j] + binsMinEdges[j])/2,(binsMaxEdges[j+1] + binsMinEdges[j+1])/2,(binsMaxEdges[j+2] + binsMinEdges[j+2])/2]
        adatapoint = array('d')
        adatapoint.fromlist(datapoint)
        print "Node: {0} ".format( kdtree_2.FindNode(adatapoint) - kdtree_2.GetNNodes())  

        t1.SetBinContent(  (kdtree.FindNode(adatapoint) - kdtree.GetNNodes())*2 + 1, kde.GetBinContent(i))  
        t1_1.SetBinContent((kdtree_1.FindNode(adatapoint) - kdtree_1.GetNNodes())*2 + 1, kde_1.GetBinContent(i))  
        t1_2.SetBinContent((kdtree_2.FindNode(adatapoint) - kdtree_2.GetNNodes())*2 + 1, kde_2.GetBinContent(i))  

    setStyle()

    canvas = TCanvas('fig_kdtreebinning', 'KDTree binning', 10, 10, 500, 500)
    h2pol.Draw('COL')
    h2pol.Draw('same')
    canvas.SaveAs(".png")
    canvas.SaveAs(".pdf")


    h2pol2 = TH2Poly("h2plot", "KDTree binning", 0, 1, 0, 1.)
    fixhist2(h2pol2, 'm_{4l}', 'D')
    
    binsMinEdges = kde_2.GetBinsMinEdges()
    binsMaxEdges = kde_2.GetBinsMaxEdges()
    
    for i in xrange(nbins):
        j = i*nvars
        h2pol2.AddBin(binsMinEdges[j]+1, binsMinEdges[j + 2],
                     binsMaxEdges[j]+1, binsMaxEdges[j + 2])
        h2pol2.SetBinContent(i+1, kde_2.GetBinDensity(i))
    setStyle()
    canvas = TCanvas('fig_kdtreebinning_2', 'KDTree binning 2', 10, 10, 500, 500)
    h2pol2.Draw('COL')
    h2pol2.Draw('same')
    canvas.SaveAs(".png")
    canvas.SaveAs(".pdf")
    
    h2pol.Write("h2pol".format(channel) )
    h2pol2.Write("h2pol2".format(channel) )

    t1.Write("t1_massErr_{0}".format(channel) )
    t1_1.Write("t1_cat1_{0}".format(channel) )
    t1_2.Write("t1_cat2_{0}".format(channel) )


    fout.Close()


    print "bins: " 
    print kdtree.GetNNodes()

    sleep(10)
#------------------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print
    print "ciao!"
    
