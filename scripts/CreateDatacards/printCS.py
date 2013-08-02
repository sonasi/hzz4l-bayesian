from ROOT import *
ROOT.gSystem.AddIncludePath("-I$ROOFITSYS/include/")
ROOT.gSystem.AddIncludePath("-Iinclude/")
ROOT.gSystem.Load("libRooFit")
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
ROOT.gSystem.Load("include/HiggsCSandWidth_cc.so")
ROOT.gSystem.Load("include/HiggsCSandWidthSM4_cc.so")

eras = [7,8]

masses = [[115, 120, 125, 130, 140, 150, 160, 170],["115", "120", "125", "130", "135", "140", "150", "160", "170", "180"]]
for iera,era in enumerate(eras):
	for signal in range(6):
		print "cs_" + str(signal) + "_" + str(era) + " = ["
		masses_ = masses[iera]
		for index, mass in enumerate(masses_):
			myCSW = HiggsCSandWidth()
			CS_ggH = myCSW.HiggsCS(signal,float(mass),era)
			BR_HZZ = myCSW.HiggsBR(11, float(mass))
			print 'Mass {0}, CS {1}, BR(H->ZZ): {2}'.format(mass, CS_ggH, BR_HZZ)
		print "]"
