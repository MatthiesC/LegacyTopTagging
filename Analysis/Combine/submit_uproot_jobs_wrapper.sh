#!/bin/bash

source /cvmfs/cms.cern.ch/cmsset_default.sh

cmssw-el7 --bind /afs:/afs --bind /nfs:/nfs --bind /pnfs:/pnfs --bind /cvmfs:/cvmfs --bind /var/lib/condor:/var/lib/condor --command-to-run "source ~/uhh2-106X-v2.sh; cd /nfs/dust/cms/user/matthies/uhh2-106X-v2/CMSSW_10_6_28/src/UHH2/LegacyTopTagging/Analysis/Combine; create_root_files_for_datacards_uproot---hacked-for-Anna.py $1"
