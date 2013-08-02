#!/usr/bin/env python
#-----------------------------------------------------------------------------
import os, sys, re
from string import *
from time import sleep
from random import randint
#------------------------------- ----------------------------------------------
SH = '''#!/bin/sh

#$ -S /bin/sh
#$ -o $HOME/Higgs_systematics/scripts/log -j y

export HOME=/cmshome/bochenek/
echo "JOBNAME: (%(jobname)s)"
echo "HOME DIRECTORY"
echo $HOME

cd $HOME
source setup.sh
#cd CMSSW_5_3_4/src
#cmsenv
cd $HOME

echo "=========================================="
time cp -r Higgs_systematics /tmp/
cd /tmp/Higgs_systematics/scripts
ls 
echo "PWD " `pwd`
echo "=========================================="
echo "STARTED/%(jobname)s/" `date` "/" `hostname`
time python make_templates_syst.py %(jobnum)s %(seed)s
echo "ENDED/%(jobname)s/" `date` "/" `hostname`
'''
#-----------------------------------------------------------------------------
def main():
	
	print "\n\t<=== subjobs.py ===>"

	njobs = 1000
	for index in xrange(njobs):

		seed = randint(1000, 100000)

		jobname = "higgs%3.3d" % index
				
		namemap = {'seed':  seed,
			   'jobnum': index,
			   'jobname' : jobname}					
			
		script = SH % namemap

		open("batch/%s.sh" %  jobname, "w").write(script)
				
		cmd = 'qsub -q local batch/%s.sh' % jobname
		os.system(cmd)
		sleep(5)
			
		if index % 10 == 0: os.system("qstat -u bochenek")

	sleep(5)
	os.system("qstat -u bochenek")
#-----------------------------------------------------------------------------
main()

