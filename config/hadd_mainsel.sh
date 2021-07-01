#!/usr/bin/env bash

# arguments
# while getopts "y:s:" args: do
#   case "${args}" in
#     y)
#       y=${args}
year=$1
syst_type=$2

outputDir="${CMSSW_BASE}/src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel/${year}/${syst_type}/"
if [! -d "${outputDir}"]; then
  echo "${outputDir} does not exist! Exit."
  exit 1
fi

haddDir="${outputDir}/hadded/"
mkdir -p ${haddDir}

# hadd ${haddDir}/uhh2.AnalysisModuleRunner.MC.TTbar__AllMergeScenarios_UL17.root ${outputDir}/
