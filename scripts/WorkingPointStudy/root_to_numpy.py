import os
import ROOT
import numpy as np
import pandas as pd

outputDirPath = os.environ.get("CMSSW_BASE")+'/src/UHH2/LegacyTopTagging/output/WorkingPointStudy/UL17/'
fileName_prefix = 'uhh2.AnalysisModuleRunner.MC.'
fileName_postfix = '.root.restructured'

# processes = ['TTbarToHadronic_UL17', 'QCD_HT300toInf_UL17']
# processes = ['TTbarToHadronic_UL17']
processes = ['QCD_HT300toInf_UL17']

for process in processes:
    print 'Working on', process
    filePath = outputDirPath+fileName_prefix+process+fileName_postfix
    file = ROOT.TFile.Open(filePath, 'READ')
    tree = file.Get('my_tree')
    leafs = ['weight', 'pt', 'msd', 'subdeepcsv', 'tau32']
    if 'TTbar' in process: leafs.append('dr')
    data = tree.AsMatrix(leafs)
    len_data_before = len(data)
    data = pd.DataFrame(data).dropna().values
    len_data_after = len(data)
    print 'Dropped', len_data_before-len_data_after, 'jets with NaN values'
    fileName_numpy = outputDirPath+fileName_prefix+process+'.npy'
    np.save(fileName_numpy, data)
