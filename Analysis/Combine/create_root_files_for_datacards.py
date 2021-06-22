import os
import sys
import argparse
import numpy as np
import ROOT as root
import subprocess
import multiprocessing
from collections import OrderedDict
from copy import deepcopy

pt_bins = [
'Pt300to400',
'Pt400to480',
'Pt480to600',
'Pt600toInf',
'Pt300toInf',
# 'Pt400toInf',
]

jet_versions = [
'All',
# 'Mass',
'BTag',
# 'MassAndBTag',
]

wps = [
'BkgEff0p001',
'BkgEff0p005',
'BkgEff0p010',
'BkgEff0p025',
'BkgEff0p050',
'NoTau32Cut',
]

regions = [
'Pass',
'Fail',
]

variables = [
'mSD',
# 'tau32',
]

systs = OrderedDict([ # naming used in file path names: naming used for combine (None = syst is only used to calculate a proper systematic later on, see code further below)
('nominal', None),
('btagging_down_bc', 'btaggingbcDown'),
('btagging_up_bc', 'btaggingbcUp'),
('btagging_down_udsg', 'btaggingudsgDown'),
('btagging_up_udsg', 'btaggingudsgUp'),
('jec_down', 'jecDown'),
('jec_up', 'jecUp'),
('jer_down', 'jerDown'),
('jer_up', 'jerUp'),
('muf_down', None),
('muf_up', None),
# ('muonid_down', 'muonidDown'), # systematic effect of muon id efficiency very small, ignore it
# ('muonid_up', 'muonidUp'),
# ('muontrigger_down', 'muontriggerDown'), # systematic effect of muon trigger efficiency very small, ignore it
# ('muontrigger_up', 'muontriggerUp'),
('mur_down', None),
('mur_up', None),
('murmuf_down', None),
('murmuf_up', None),
('pileup_down', 'pileupDown'),
('pileup_up', 'pileupUp'),
('ps_FSRdown_2', 'fsrDown'),
('ps_FSRup_2', 'fsrUp'),
# ('ps_ISRdown_2', 'isrDown'), # ISR weights seem a bit broken, exclude them for now (some events get stupendiously large weights, especially when using the Up variation)
# ('ps_ISRup_2', 'isrUp'),
('wp_down', None),
('wp_up', None),
])

processes = OrderedDict([ # naming used in root file names: naming used for combine
### TTbar
('TTbar__FullyMerged', 'TTbar_FullyMerged'),
('TTbar__WMerged', 'TTbar_WMerged'),
('TTbar__QBMerged', 'TTbar_QBMerged'),
('TTbar__BkgOrNotMerged', 'TTbar_BkgOrNotMerged'),
### Single Top
('ST__FullyMerged', 'ST_FullyMerged'),
('ST__WMerged', 'ST_WMerged'),
('ST__QBMerged', 'ST_QBMerged'),
('ST__BkgOrNotMerged', 'ST_BkgOrNotMerged'),
### Other
('WJetsToLNu', 'WJetsToLNu'),
('DYJetsToLLAndDiboson', 'DYJetsToLLAndDiboson'),
('QCD_Mu', 'QCD_Mu'),
### Data is handeled manually in the code below
])

inputBaseDir = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel/UL17')
fileNamePrefix = 'uhh2.AnalysisModuleRunner.'

dict_inputFiles = OrderedDict()
for syst in systs.keys():
    dict_inputFiles[syst] = OrderedDict()
    syst_path_string = 'syst_'+syst if syst != 'nominal' else syst
    for proc in processes.keys():
        dict_inputFiles[syst][proc] = os.path.join(inputBaseDir, syst_path_string, 'hadded', fileNamePrefix+'MC.'+proc+'.root')
dict_inputFiles['nominal']['DATA'] = os.path.join(inputBaseDir, 'nominal', 'hadded', fileNamePrefix+'DATA.DATA.root')

# for i in dict_inputFiles.keys():
#     print i
#     for j in dict_inputFiles[i].keys():
#         print j
#         print dict_inputFiles[i][j]


outputBaseDir = os.path.join(inputBaseDir, 'combine')
os.system('mkdir -p '+outputBaseDir)

def hist_preprocessing(hist):
    # hist.Rebin(20)
    rebinning_scheme = np.array([0, 10, 25, 40, 55, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 215, 230, 250, 275, 300, 350, 400, 450, 500], dtype=float)
    hist_new = hist.Rebin(len(rebinning_scheme)-1, 'hnew', rebinning_scheme)
    for i in range(hist_new.GetNbinsX()+2): # set bins with bin content smaller zero to zero (thus, removing an artefact of negative event weights)
        if hist_new.GetBinContent(i) < 0:
            hist_new.SetBinContent(i, 0)
    return hist_new

for pt_bin in pt_bins:
    for jet_version in jet_versions:
        for wp in wps:
            task_name = '_'.join([pt_bin, jet_version, wp])
            print 'Working on', task_name
            outputDir = os.path.join(outputBaseDir, 'AK8', task_name)
            os.system('mkdir -p '+outputDir)
            outputRootFileName = task_name+'__input_histograms.root'
            outputRootFilePath = os.path.join(outputDir, outputRootFileName)
            outputRootFile = root.TFile.Open(outputRootFilePath, 'RECREATE')
            inputHistCollectionName = 'ProbeJetHists_AK8'
            for region in regions:
                for variable in variables:
                    channelName = '_'.join([task_name, region, variable])
                    target_folder = outputRootFile.mkdir(channelName)
                    inputHistName = '/'.join([inputHistCollectionName, channelName])
                    # data:
                    inputRootFile = root.TFile.Open(dict_inputFiles.get('nominal').get('DATA'), 'READ')
                    inputHist = inputRootFile.Get(inputHistName)
                    inputHist = hist_preprocessing(inputHist)
                    target_folder.cd()
                    inputHist.Write('data_obs')
                    inputRootFile.Close()
                    # MC:
                    for proc in processes.keys():
                        # nominal:
                        inputRootFile = root.TFile.Open(dict_inputFiles.get('nominal').get(proc), 'READ')
                        inputHist = inputRootFile.Get(inputHistName)
                        inputHist = hist_preprocessing(inputHist)
                        target_folder.cd()
                        inputHist.Write(processes.get(proc))
                        inputRootFile.Close()
                        # variations for systematics that don't need special computational treatment here:
                        for syst in systs.keys():
                            if systs.get(syst) == None: continue # None = this systematic needs special treatment to be correctly implemented here
                            inputRootFile = root.TFile.Open(dict_inputFiles.get(syst).get(proc), 'READ')
                            inputHist = inputRootFile.Get(inputHistName)
                            inputHist = hist_preprocessing(inputHist)
                            target_folder.cd()
                            inputHist.Write('_'.join([processes.get(proc), systs.get(syst)]))
                            inputRootFile.Close()
                        # variations for systematics needing special computational treatment:
                        ## QCD scale:
                        scale_variations = ['muf_down', 'muf_up', 'mur_down', 'mur_up', 'murmuf_down', 'murmuf_up']
                        inputHists_scale = dict()
                        for sv in scale_variations:
                            inputRootFile = root.TFile.Open(dict_inputFiles.get(sv).get(proc), 'READ')
                            inputHist = deepcopy(inputRootFile.Get(inputHistName))
                            inputHist = hist_preprocessing(inputHist)
                            inputHists_scale[sv] = deepcopy(inputHist)
                            inputRootFile.Close()
                        inputHist_muScale_Down = deepcopy(inputHists_scale.get(scale_variations[0])) # doesn't really matter which histogram is copied; we just need a TH1F template with the correct binning
                        inputHist_muScale_Up = deepcopy(inputHists_scale.get(scale_variations[0]))
                        for i in range(inputHist_muScale_Down.GetNbinsX()+2):
                            bin_contents = np.array([], dtype=float)
                            for hist in inputHists_scale.keys():
                                bin_contents = np.append(bin_contents, inputHists_scale.get(hist).GetBinContent(i))
                            inputHist_muScale_Down.SetBinContent(i, 0)
                            inputHist_muScale_Down.SetBinError(i, 0)
                            inputHist_muScale_Down.SetBinContent(i, bin_contents.min())
                            inputHist_muScale_Up.SetBinContent(i, 0)
                            inputHist_muScale_Up.SetBinError(i, 0)
                            inputHist_muScale_Up.SetBinContent(i, bin_contents.max())
                        target_folder.cd()
                        inputHist_muScale_Down.Write('_'.join([processes.get(proc), 'scaleDown']))
                        inputHist_muScale_Up.Write('_'.join([processes.get(proc), 'scaleUp']))
                        ## Tagging efficiency by number of prongs:
                        ### 3-prong:
                        is_3prong = 'FullyMerged' in proc
                        is_2prong = 'WMerged' in proc or 'QBMerged' in proc or 'SemiMerged' in proc
                        is_1prong = not (is_3prong or is_2prong)
                        prong_types = OrderedDict([
                            ('3prong', is_3prong),
                            ('2prong', is_2prong),
                            ('1prong', is_1prong),
                        ])
                        for prong_type in prong_types.keys():
                            if prong_types.get(prong_type) and variable != 'tau32': # to get the wp variation in sframe, the tau32 cut was varied. It does not make sense to have variation along the x axis (what would be the case for tau32 hists)
                                inputRootFile_prong_Down = root.TFile.Open(dict_inputFiles.get('wp_down').get(proc), 'READ')
                                inputRootFile_prong_Up = root.TFile.Open(dict_inputFiles.get('wp_up').get(proc), 'READ')
                            else:
                                inputRootFile_prong_Down = root.TFile.Open(dict_inputFiles.get('nominal').get(proc), 'READ')
                                inputRootFile_prong_Up = root.TFile.Open(dict_inputFiles.get('nominal').get(proc), 'READ')
                            inputHist_prong_Down = inputRootFile_prong_Down.Get(inputHistName)
                            inputHist_prong_Up = inputRootFile_prong_Up.Get(inputHistName)
                            inputHist_prong_Down = hist_preprocessing(inputHist_prong_Down)
                            inputHist_prong_Up = hist_preprocessing(inputHist_prong_Up)
                            target_folder.cd()
                            inputHist_prong_Down.Write('_'.join([processes.get(proc), 'tageff'+prong_type+'Down']))
                            inputHist_prong_Up.Write('_'.join([processes.get(proc), 'tageff'+prong_type+'Up']))
                            inputRootFile_prong_Down.Close()
                            inputRootFile_prong_Up.Close()
            outputRootFile.Close()
            print 'Wrote', outputRootFilePath
