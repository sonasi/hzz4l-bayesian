#!/bin/sh

#$ -S /bin/sh
#$ -o $HOME/Higgs_systematics/scripts/log -j y

export HOME=/cmshome/bochenek/
echo "JOBNAME: (higgs000)"
echo "HOME DIRECTORY"
echo $HOME

cd $HOME
source setup.sh
#cd CMSSW_5_3_4/src
#cmsenv
cd $HOME

echo "=========================================="
time cp -r Higgs_systematics /tmp
cd /tmp/Higgs_systematics/scripts
echo "PWD " `pwd`
echo "=========================================="
echo "STARTED/higgs000/" `date` "/" `hostname`
time python make_templates_syst.py 0 0
echo "ENDED/higgs000/" `date` "/" `hostname`
