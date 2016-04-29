#!/bin/bash

export SCRAM_ARCH=slc6_amd64_gcc493
cd /afs/cern.ch/user/y/yohay/scratch0/CMSSW_7_6_3/src
eval `scramv1 runtime -sh`
cd -
cp /afs/cern.ch/user/y/yohay/scratch0/CMSSW_7_6_3/src/GGHAA2Mu2TauAnalysis/MuMuTauTauSkimmer/test/merge.py .
cp /afs/cern.ch/user/y/yohay/scratch0/CMSSW_7_6_3/src/GGHAA2Mu2TauAnalysis/MuMuTauTauSkimmer/test/NMSSM_ggH_a9_H1125_H2500_H3500_2mu2tau.txt .
cmsRun merge.py
/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select cp heavyHiggs_125_light_9_2mu2tau_STEP2.root /eos/cms/store/user/yohay/
rm merge.* heavyHiggs_125_light_9_2mu2tau_STEP2.root NMSSM_ggH_a9_H1125_H2500_H3500_2mu2tau.txt

exit 0
