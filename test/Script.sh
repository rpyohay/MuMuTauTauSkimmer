#!/bin/bash
#parse arguments
if [ $# -ne 6 ]
    then
    echo "Usage: ./Script.sh dir_name"
    exit 0
fi

dir_name=$1
#make directory on EOS
EOS_dir_query=`cmsLs /store/user/mshi/${dir_name}`
EOS_dir_query=`echo $EOS_dir_query | grep "No such file or directory"`
if [ "EOS_dir_query" != "" ]
    then
    cmsMkdir /store/user/mshi/${dir_name}
fi
sed -e "s%DIRNAME%${dir_name}%g"

cd ../../../
export SCRAM_ARCH=slc6_amd64_gcc493
cd /afs/cern.ch/user/m/mshi/CMSSW_7_6_3/src/GGHAA2Mu2TauAnalysis/MuMuTauTauSkimmer/test
eval `scramv1 runtime -sh`
cd -
source /afs/cern.ch/cms/caf/setup.sh
cp /afs/cern.ch/user/m/mshi/CMSSW_7_6_3/src/GGHAA2Mu2TauAnalysis/MumuTauTauSkimmer/test/tauSelectionSkim_MC.py .
cmsRun tauSelectionSkim_MC.py
eos cp -f DrellYan_Skim.root DrellYan_Mu45Selector.root DrellYan_CleanJets.root DrellYan_muHadTauSelector.root DrellYan_RECOAnalyzer.root DrellYan_FileService.root /store/user/mshi/DIRNAME
rm DrellYan_Skim.root DrellYan_Mu45Selector.root DrellYan_CleanJets.root DrellYan_muHadTauSelector.root DrellYan_RECOAnalyzer.root DrellYan_FileService.root tauSelectionSkim_MC.py
exit 0
