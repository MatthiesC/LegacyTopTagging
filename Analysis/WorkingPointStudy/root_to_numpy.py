import os
import sys
import argparse
import ROOT
import numpy as np
import pandas as pd


options_year = ['UL16preVFP', 'UL16postVFP', 'UL17', 'UL18']

if not sys.argv[1:]: sys.exit('No arguments provided. Exit.')
parser = argparse.ArgumentParser()
parser.add_argument('-y', '--year', type=str, choices=options_year)
args = parser.parse_args(sys.argv[1:])

year = args.year

outputDirPath = os.environ.get("CMSSW_BASE")+'/src/UHH2/LegacyTopTagging/output/WorkingPointStudy/'+year+'/'
fileName_prefix = 'uhh2.AnalysisModuleRunner.MC.'
fileName_postfix = '.root.restructured'

processes = ['TTbarToHadronic', 'QCD_HT300toInf']
# processes = ['TTbarToHadronic']
# processes = ['QCD_HT300toInf']

for proc in processes:
    process = proc+'_'+year
    print 'Working on', process
    filePath = outputDirPath+fileName_prefix+process+fileName_postfix
    file = ROOT.TFile.Open(filePath, 'READ')
    tree = file.Get('my_tree')
    leafs = ['weight', 'pt', 'msd', 'subdeepcsv', 'subdeepjet', 'tau32']
    if 'TTbar' in process: leafs.append('dr')
    data = tree.AsMatrix(leafs)
    len_data_before = len(data)
    data = pd.DataFrame(data).dropna().values
    len_data_after = len(data)
    print 'Dropped', len_data_before-len_data_after, 'jets with NaN values'
    print 'Final number of jets:', len_data_after
    fileName_numpy = outputDirPath+fileName_prefix+process+'.npy'
    np.save(fileName_numpy, data)
