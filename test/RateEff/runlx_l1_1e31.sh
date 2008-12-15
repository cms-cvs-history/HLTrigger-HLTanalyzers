#!/bin/bash

cd /afs/cern.ch/user/c/chinhan/scratch0/TMD/configRates/CMSSW_2_1_12/src/HLTrigger/HLTanalyzers/test/newRateEff

if [ -n "${CMS_PATH:-}" ]; then
  echo "CMSSW computing environment already setup"
else
  export SCRAM_ARCH=`scramv1 arch`
fi
eval `scramv1 runtime -sh`

source setup.sh
./OHltRateEff l1menu_1E31_2008Dec04.cfg



