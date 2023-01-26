#!/bin/bash
echo 'start run HHbbWW with ',$0,$1,$2,$3,$4,$5
pushd /afs/cern.ch/work/t/tahuang/HHAnalysis/HeavyMassEstimator/test
# If workspace does not exist, create it once
python runHME_HHbbWW_HIG2105.py -i $1 -o $2 -nStart $3 -nEnd $4 -it $5
popd
echo 'finish HHbbWW_condor'
