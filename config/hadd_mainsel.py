import os
import sys
import argparse
import glob
import numpy as np
import ROOT as root
import subprocess
import multiprocessing

years = [
'UL17',
'UL18',
]

possible_systs = [
'nominal',
'btagging_down_bc',
'btagging_down_udsg',
'btagging_up_bc',
'btagging_up_udsg',
'jec_down',
'jec_up',
'jer_down',
'jer_up',
'muf_down',
'muf_up',
'muonid_down',
'muonid_up',
'muontrigger_down',
'muontrigger_up',
'mur_down',
'mur_up',
'murmuf_down',
'murmuf_up',
'pileup_down',
'pileup_up',
'ps_FSRdown_2',
'ps_FSRup_2',
'ps_ISRdown_2',
'ps_ISRup_2',
'wp_down',
'wp_up',
]

dict_sourceFiles = {
    ### TTbar
    'TTbar__AllMergeScenarios': [
        'TTbar*__AllMergeScenarios',
    ],
    'TTbar__FullyMerged': [
        'TTbar*__FullyMerged',
    ],
    'TTbar__SemiMerged': [
        'TTbar*__WMerged',
        'TTbar*__QBMerged',
    ],
    'TTbar__WMerged': [
        'TTbar*__WMerged',
    ],
    'TTbar__QBMerged': [
        'TTbar*__QBMerged',
    ],
    'TTbar__NotMerged': [
        'TTbar*__NotMerged',
    ],
    'TTbar__Background': [
        'TTbar*__Background',
        'TTbarTo2L2Nu*',
        'TTbarToHadronic*',
    ],
    'TTbar__BkgOrNotMerged': [
        'TTbar*__NotMerged',
        'TTbar*__Background',
        'TTbarTo2L2Nu*',
        'TTbarToHadronic*',
    ],
    ### Single Top
    'ST__AllMergeScenarios': [
        'ST*__AllMergeScenarios',
    ],
    'ST__FullyMerged': [
        'ST*__FullyMerged',
    ],
    'ST__SemiMerged': [
        'ST*__WMerged',
        'ST*__QBMerged',
    ],
    'ST__WMerged': [
        'ST*__WMerged',
    ],
    'ST__QBMerged': [
        'ST*__QBMerged',
    ],
    'ST__NotMerged': [
        'ST*__NotMerged',
    ],
    'ST__Background': [
        'ST*__Background',
        'ST_sChannel_leptonDecays*',
    ],
    'ST__BkgOrNotMerged': [
        'ST*__NotMerged',
        'ST*__Background',
        'ST_sChannel_leptonDecays*',
    ],
    ### Other
    'WJetsToLNu': [
        'WJetsToLNu_HT*',
    ],
    'DYJetsToLL': [
        'DYJetsToLL_HT*',
    ],
    'Diboson': [
        'Diboson*',
    ],
    'QCD_Mu': [
        'QCD_Mu_Pt*',
    ],
    'DATA': [
        'DATA*',
    ]
}

fileNamePrefix = 'uhh2.AnalysisModuleRunner.'

parser = argparse.ArgumentParser()
parser.add_argument('-y', '--years', choices=years, nargs='*', default=[])
parser.add_argument('-s', '--systs', choices=possible_systs, nargs='*', default=['nominal'])
parser.add_argument('--all', action='store_true', help='Instead of defining the syst directories via --syst, you can hadd all of them in one go.')
args = parser.parse_args(sys.argv[1:])

if not args.all:
    args_systs = args.systs
else:
    args_systs = possible_systs
print args_systs

# sys.exit()

mainselOutputDir = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel')

hadd_tasks = list() # pairs of command_string, logFilePath (see below)

for year in args.years:
    for syst in args_systs:
        if syst=='nominal':
            nominal=True
        else:
            nominal=False
        outputDir = os.path.join(mainselOutputDir, year, ('' if nominal else 'syst_')+syst)
        haddDir = os.path.join(outputDir, 'hadded')
        logDir = os.path.join(haddDir, 'log')
        os.system('mkdir -p '+logDir)
        # print haddDir

        for key in dict_sourceFiles.keys():
            if key=='DATA':
                data=True
            else:
                data=False
            if not nominal and data: continue
            fileNamePrefix_ = fileNamePrefix+('DATA.' if data else 'MC.')
            targetFilePath = os.path.join(haddDir, fileNamePrefix_+key+'.root')
            # print targetFilePath
            sourceFilePaths = np.array([]) # list of all input root files to be hadded.
            for pattern in dict_sourceFiles[key]:
                sourceFilePaths = np.append(sourceFilePaths, glob.glob(os.path.join(outputDir, fileNamePrefix_+pattern+'_'+year+'.root')))
            sourceFilePaths.flatten()

            # In the following few lines, we sort the root files which are to be hadded by the number of events stored in their analysis trees
            # (there is a bug in hadd that would lead to not correctly added trees if the first file in the list has no events; in case that all trees are empty, we have nothing to fear)
            numbersOfEntries = np.array([])
            for rootFileName in sourceFilePaths:
                rootFile = root.TFile.Open(rootFileName, 'READ')
                numbersOfEntries = np.append(numbersOfEntries, rootFile.Get('AnalysisTree').GetEntries())
            numbersOfEntries.flatten()
            sourceFilePaths = sourceFilePaths[np.argsort(numbersOfEntries)]

            command_string = 'nice -n 10 hadd '+targetFilePath+' '+' '.join([x for x in sourceFilePaths])
            logFilePath = os.path.join(logDir, 'log.'+key+'.txt')
            hadd_tasks.append([command_string, logFilePath])

# print hadd_tasks

FNULL = open(os.devnull, 'w')

def hadd_task(args): # args = [command, logfile]
    process = subprocess.Popen([args[0]+' > '+args[1]], shell=True, stdout=FNULL, stderr=FNULL)
    process.wait()

p = multiprocessing.Pool(multiprocessing.cpu_count() / 2) # only use half of all CPU threads to be friendly to other users
p.map(hadd_task, hadd_tasks)
