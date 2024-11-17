#!/bin/bash

cd /nfs/dust/cms/user/matthies/uhh2-106X-v2/CMSSW_10_6_28/src/UHH2/HighPtSingleTop/config

export PATH=$PATH:/nfs/dust/cms/user/matthies/SFrameBatch_el9_python3/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/nfs/dust/cms/user/matthies/SFrameBatch/hadd_tree_size_fix/

export PS1="\[\033[1;31m\]<SFrameBatch-EL9-Python3>\[\033[0m\] $PS1"
