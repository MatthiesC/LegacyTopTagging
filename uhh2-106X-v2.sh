#!/bin/bash

source ~/.bashrc

source /cvmfs/cms.cern.ch/cmsset_default.sh
#cmssw-el7

export SCRAM_ARCH=slc7_amd64_gcc700

cd /nfs/dust/cms/user/matthies/uhh2-106X-v2/CMSSW_10_6_28/src
cmsenv
cd /nfs/dust/cms/user/matthies/uhh2-106X-v2/SFrame
source setup.sh
cd /nfs/dust/cms/user/matthies/uhh2-106X-v2/CMSSW_10_6_28/src/UHH2

export PATH=$PATH:/nfs/dust/cms/user/matthies/SFrameBatch/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/nfs/dust/cms/user/matthies/SFrameBatch/hadd_tree_size_fix/

export PATH=$PATH:/nfs/dust/cms/user/matthies/lwtnn/lwtnn/converters

#Need to manually source this in order to use grid commands (e.g. voms, see email from 10.6.2020)
source /cvmfs/grid.desy.de/etc/profile.d/grid-ui-env.sh

function krebs() {
    source /cvmfs/cms.cern.ch/crab3/crab.sh
    source /cvmfs/cms.cern.ch/common/crab-setup.sh
    voms-proxy-init -rfc -voms cms -valid 192:00
    export X509_USER_PROXY=`voms-proxy-info -path`
    cd $CMSSW_BASE/src/UHH2/scripts/crab
}

export PS1="\[\033[1;31m\]<uhh2-106X-v2>\[\033[0m\] $PS1"
