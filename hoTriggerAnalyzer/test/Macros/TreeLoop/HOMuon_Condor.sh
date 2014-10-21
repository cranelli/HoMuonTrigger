#!/bin/bash

#
# variables from arguments string in jdl
#

CONDOR_CLUSTER=$1
CONDOR_PROCESS=$2
RUN_DIR=$3
#PREFIX_NAME=$4
#INPUT_NAME=$5

#
# header 
#

echo ""
echo "CMSSW on Condor"
echo ""

START_TIME=`/bin/date`
echo "started at $START_TIME"

echo ""
echo "parameter set:"
echo "CONDOR_CLUSTER: $CONDOR_CLUSTER"
echo "CONDOR_PROCESS: $CONDOR_PROCESS"
echo "RUN_DIR: $RUN_DIR"

#
# setup CMS software environment at UMD
#
export VO_CMS_SW_DIR=/sharesoft/cmssw
. $VO_CMS_SW_DIR/cmsset_default.sh
cd $RUN_DIR
eval `scramv1 runtime -sh`

#
# modify parameter-set
#

SCRATCH_DIR=`echo ${_CONDOR_SCRATCH_DIR}`

#
# run cmssw
#

echo $HOME
echo "./HOMuon_Struct.C"
root -l -b HOMuon_Struct.C



#root -l -b HOMuon_Plotter.C
echo "working"

exitcode=$?

#
# end run
#

echo ""
END_TIME=`/bin/date`
echo "finished at $END_TIME"
exit $exitcode
