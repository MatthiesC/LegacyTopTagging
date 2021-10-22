#!/usr/bin/env python

import os
import sys
import argparse
import numpy as np
import ROOT as root
# import subprocess
import multiprocessing as mp
from collections import OrderedDict
from copy import deepcopy
import itertools
import math
import re

from parallel_threading import run_with_pool

from constants import years


pt_bins = {
    'AK8': [
        'Pt300toInf', # needs to be first place, see function set_variables below
        'Pt300to400',
        'Pt400to480',
        'Pt480to600',
        'Pt600toInf',
        'Pt400toInf',
    ],
    'HOTVR': [
        'Pt200toInf', # needs to be first place, see function set_variables below
        'Pt200to250',
        'Pt250to300',
        'Pt300to400',
        'Pt400to480',
        'Pt480to600',
        'Pt600toInf',
        # 'Pt300toInf',
        'Pt400toInf',
    ],
}

jet_versions = {
    'AK8': [
        'All',
        # 'Mass',
        'BTag',
        # 'MassAndBTag',
    ],
    'HOTVR': [
        'All',
        # 'Mass',
        'HOTVRCuts',
        # 'HOTVRCutsAndMass',
    ],
}

wps = {
    'AK8': [
        'BkgEff0p001',
        'BkgEff0p005',
        'BkgEff0p010',
        'BkgEff0p025',
        'BkgEff0p050',
        'NoTau32Cut',
    ],
    'HOTVR': [
        'Standard',
        'NoTau32Cut',
    ],
}

regions = [
    'Pass',
    'Fail',
    # 'PassW',
    # 'FailW',
]

variables = {
    'AK8': [
        'pt',
        'drlepton',
        'eta',
        'phi',
        'mass',
        'mSD',
        'tau32',
        'tau21',
        'maxDeepCSV',
        'nsub',
    ],
    'HOTVR': [
        'pt',
        'drlepton',
        'eta',
        'phi',
        'mass',
        'mpair',
        'tau32',
        'tau21',
        'fpt1',
        'nsub',
    ],
}

probejet_collections = { # ProbeJetHists_X : Actual jet algorithm
    'AK8': 'AK8',
    'AK8_mSD10': 'AK8',
    'HOTVR': 'HOTVR',
}

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
    # ('tau21_down', 'tau21Down'),
    # ('tau21_up', 'tau21Up'),
    ('toppt_a_up', 'topptAUp'),
    ('toppt_a_down', 'topptADown'),
    ('toppt_b_up', 'topptBUp'),
    ('toppt_b_down', 'topptBDown'),
])

processes = OrderedDict([ # naming used in root file names: naming used for combine
    ### TTbar
    ('TTbar__AllMergeScenarios', 'TTbar'),
    ('TTbar__FullyMerged', 'TTbar_FullyMerged'),
    # ('TTbar__WMerged', 'TTbar_WMerged'),
    # ('TTbar__QBMerged', 'TTbar_QBMerged'),
    # ('TTbar__SemiMerged', 'TTbar_SemiMerged'),
    # ('TTbar__BkgOrNotMerged', 'TTbar_BkgOrNotMerged'),
    # ('TTbar__NotFullyOrWMerged', 'TTbar_NotFullyOrWMerged'),
    ('TTbar__YllufMerged', 'TTbar_YllufMerged'),
    # ('TTbar__NotMerged', 'TTbar_NotMerged'),
    ('TTbar__Background', 'TTbar_Background'),
    ### Single Top
    ('ST__AllMergeScenarios', 'ST'),
    ('ST__FullyMerged', 'ST_FullyMerged'),
    # ('ST__WMerged', 'ST_WMerged'),
    # ('ST__QBMerged', 'ST_QBMerged'),
    # ('ST__SemiMerged', 'ST_SemiMerged'),
    # ('ST__BkgOrNotMerged', 'ST_BkgOrNotMerged'),
    # ('ST__NotFullyOrWMerged', 'ST_NotFullyOrWMerged'),
    ('ST__YllufMerged', 'ST_YllufMerged'),
    # ('ST__NotMerged', 'ST_NotMerged'),
    ('ST__Background', 'ST_Background'),
    ### Other
    ('WJetsToLNu', 'WJetsToLNu'),
    ('DYJetsToLLAndDiboson', 'DYJetsToLLAndDiboson'),
    ('QCD_Mu', 'QCD_Mu'),
    ### Data is handeled manually in the code below
])

parser = argparse.ArgumentParser()
parser.add_argument('-y', '--year', choices=years.keys(), nargs=1)
args = parser.parse_args(sys.argv[1:])

year = years.get(args.year[0])

inputBaseDir = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel', year.get('short_name'))
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

def hist_preprocessing(histo, variable):

    new_histo = histo

    if variable == 'mSD' or variable == 'mass':
        # rebinning_scheme = np.array([0, 10, 25, 40, 55, 70, 85, 100, 115, 130, 145, 160, 175, 190, 210, 235, 270, 310, 350, 390, 450, 500], dtype=float) # Dennis' bins
        # Ensure to have bin edges 105, 140, 210, 220 available (for HOTVR/AK8 mass efficiency calculation)
        # rebinning_scheme = np.array([0, 10, 25, 40, 55, 70, 80, 90, 105, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 250, 275, 300, 350, 400, 450, 500], dtype=float) # before July 15, 2021
        # rebinning_scheme = np.array([0, 10, 25, 40, 55, 70, 85, 105, 120, 140, 155, 170, 185, 200, 210, 220, 230, 250, 275, 300, 350, 400, 450, 500], dtype=float) # before September 13, 2021
        rebinning_scheme = np.array([50, 70, 85, 105, 120, 140, 155, 170, 185, 200, 210, 220, 230, 250, 275, 300, 350, 400, 450, 500], dtype=float)
        # rebinning_scheme = np.array([0, 10, 25, 50, 70, 90, 105, 120, 140, 160, 180, 195, 210, 220, 240, 270, 300, 350, 400, 450, 500], dtype=float)
        new_histo = new_histo.Rebin(len(rebinning_scheme)-1, 'hnew', rebinning_scheme)
    # if variable == 'mSD':
    #     rebinning_scheme = np.array([0, 10, 30, 50, 70, 90, 105, 120, 135, 160, 175, 190, 200, 210, 230, 250, 275, 300, 350, 400, 450, 500], dtype=float)
    #     new_histo = new_histo.Rebin(len(rebinning_scheme)-1, 'hnew', rebinning_scheme)
    # elif variable == 'mass':
    #     rebinning_scheme = np.array([0, 30, 50, 70, 90, 105, 120, 140, 155, 170, 185, 200, 210, 220, 230, 250, 275, 300, 350, 400, 450, 500], dtype=float)
    #     new_histo = new_histo.Rebin(len(rebinning_scheme)-1, 'hnew', rebinning_scheme)
    elif variable == 'nsub':
        # new_histo.GetXaxis().SetNdivisions(11)
        pass
    else:
        new_histo.Rebin(20)

    # function to include overflow in last bin or include underflow in first bin:
    def include_overflow_underflow(hist, overflow=True):
        index_overflow = hist.GetNbinsX()+1 if overflow else 0
        index_last = hist.GetNbinsX() if overflow else 1
        overflow_binc = hist.GetBinContent(index_overflow)
        overflow_err = hist.GetBinError(index_overflow)
        last_binc = hist.GetBinContent(index_last)
        last_err = hist.GetBinError(index_last)
        newlast_binc = last_binc + overflow_binc
        newlast_err = math.sqrt(last_err*last_err + overflow_err*overflow_err) # Gaussian error propagation (addition case)
        hist.SetBinContent(index_last, newlast_binc)
        hist.SetBinError(index_last, newlast_err)
        hist.SetBinContent(index_overflow, 0)
        hist.SetBinError(index_overflow, 0)

    include_overflow_underflow(new_histo, True) # always include overflow
    if variable == 'eta' or variable == 'phi' or variable == 'maxDeepCSV': # only inlcude underflow bin for those variables
        include_overflow_underflow(new_histo, False)

    for i in range(new_histo.GetNbinsX()+2): # set bins with bin content smaller zero to zero (thus, removing an artefact of negative event weights)
        if new_histo.GetBinContent(i) < 0:
            new_histo.SetBinContent(i, 0)
            new_histo.SetBinError(i, 0)

    return new_histo


# for pt_bin in pt_bins:
#     for jet_version in jet_versions:
#         for wp in wps:
def create_input_histograms(probejet_coll, pt_bin, jet_version, wp, vars):
    jet_algo = probejet_collections.get(probejet_coll) # AK8 and AK8_mSD10 will be put in the same directory!
    task_name = '_'.join([probejet_collections.get(probejet_coll), pt_bin, jet_version, wp])
    print 'Working on', task_name, probejet_coll
    outputDir = os.path.join(outputBaseDir, jet_algo, task_name)
    os.system('mkdir -p '+outputDir)
    inputHistCollectionName = 'ProbeJetHists_'+probejet_coll
    outputRootFileName = task_name+'__'+inputHistCollectionName+'.root'
    outputRootFilePath = os.path.join(outputDir, outputRootFileName)
    outputRootFile = root.TFile.Open(outputRootFilePath, 'RECREATE')
    channelNames = list()
    for variable in vars:
        for region in regions:
            channelName = '_'.join([task_name, region, variable])
            channelNames.append(channelName)
            target_folder = outputRootFile.mkdir(channelName)
            inputHistName = '/'.join([inputHistCollectionName, channelName.replace(probejet_collections.get(probejet_coll)+'_', '')])
            # data:
            inputRootFile = root.TFile.Open(dict_inputFiles.get('nominal').get('DATA'), 'READ')
            inputHist = inputRootFile.Get(inputHistName)
            inputHist = hist_preprocessing(inputHist, variable)
            target_folder.cd()
            inputHist.Write('data_obs')
            inputRootFile.Close()
            # MC:
            for proc in processes.keys():
                # nominal:
                inputRootFile = root.TFile.Open(dict_inputFiles.get('nominal').get(proc), 'READ')
                inputHist = inputRootFile.Get(inputHistName)
                inputHist = hist_preprocessing(inputHist, variable)
                target_folder.cd()
                inputHist.Write(processes.get(proc))
                inputRootFile.Close()
                # variations for systematics that don't need special computational treatment here:
                for syst in systs.keys():
                    if systs.get(syst) == None: continue # None = this systematic needs special treatment to be correctly implemented here
                    inputRootFile = root.TFile.Open(dict_inputFiles.get(syst).get(proc), 'READ')
                    inputHist = inputRootFile.Get(inputHistName)
                    inputHist = hist_preprocessing(inputHist, variable)
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
                    inputHist = hist_preprocessing(inputHist, variable)
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
                # ## Tagging efficiency by number of prongs:
                # is_3prong = 'FullyMerged' in proc
                # is_2prong = ('WMerged' in proc and not 'NotFullyOrWMerged' in proc) or 'QBMerged' in proc or 'SemiMerged' in proc
                # is_1prong = not (is_3prong or is_2prong)
                # prong_types = OrderedDict([
                #     ('3prong', {'bool': is_3prong, 'eff_variation': 0.05}), # https://github.com/UHH2/TopTagging2018/blob/master/src/ProbeJetHists.cxx#L122
                #     ('2prong', {'bool': is_2prong, 'eff_variation': 0.1}),
                #     ('1prong', {'bool': is_1prong, 'eff_variation': 0.1}),
                # ])
                ## Tagging efficiency by merge category:
                is_Fully = 'FullyMerged' in proc
                is_Ylluf = 'YllufMerged' in proc
                is_Bkgrd = not (is_Fully or is_Ylluf)
                prong_types = OrderedDict([
                    ('Fully', {'bool': is_Fully, 'eff_variation': 0.05}), # https://github.com/UHH2/TopTagging2018/blob/master/src/ProbeJetHists.cxx#L122
                    ('Ylluf', {'bool': is_Ylluf, 'eff_variation': 0.1}),
                    ('Bkgrd', {'bool': is_Bkgrd, 'eff_variation': 0.1}),
                ])
                for prong_type in prong_types.keys():
                    use_wp_variation = True; # switch how to implement the "prong tagging efficiency uncertainty"
                    if use_wp_variation and prong_types.get(prong_type).get('bool') and variable != 'tau32': # to get the wp variation in sframe, the tau32 cut was varied. It does not make sense to have variation along the x axis (what would be the case for tau32 hists)
                        inputRootFile_prong_Down = root.TFile.Open(dict_inputFiles.get('wp_down').get(proc), 'READ')
                        inputRootFile_prong_Up = root.TFile.Open(dict_inputFiles.get('wp_up').get(proc), 'READ')
                    else:
                        inputRootFile_prong_Down = root.TFile.Open(dict_inputFiles.get('nominal').get(proc), 'READ')
                        inputRootFile_prong_Up = root.TFile.Open(dict_inputFiles.get('nominal').get(proc), 'READ')
                    inputHist_prong_Down = inputRootFile_prong_Down.Get(inputHistName)
                    inputHist_prong_Up = inputRootFile_prong_Up.Get(inputHistName)
                    inputHist_prong_Down = hist_preprocessing(inputHist_prong_Down, variable)
                    inputHist_prong_Up = hist_preprocessing(inputHist_prong_Up, variable)
                    if not use_wp_variation and prong_types.get(prong_type).get('bool') and wp != 'NoTau32Cut':
                        eff_variation = prong_types.get(prong_type).get('eff_variation')
                        if region == 'Pass':
                            inputHist_prong_Down.Scale(1.-eff_variation)
                            inputHist_prong_Up.Scale(1.+eff_variation)
                        # elif region == 'Fail':
                        elif region in ['Fail', 'PassW', 'FailW']:
                            inputHist_prong_Down.Scale(1.+eff_variation)
                            inputHist_prong_Up.Scale(1.-eff_variation)
                    target_folder.cd()
                    inputHist_prong_Down.Write('_'.join([processes.get(proc), 'tageff'+prong_type+'Down']))
                    inputHist_prong_Up.Write('_'.join([processes.get(proc), 'tageff'+prong_type+'Up']))
                    inputRootFile_prong_Down.Close()
                    inputRootFile_prong_Up.Close()
    outputRootFile.Close()
    print 'Wrote', outputRootFilePath
    return [outputRootFilePath, channelNames]


def create_input_histograms_task(task):
    return create_input_histograms(task.get('probejet_coll'), task.get('pt_bin'), task.get('jet_version'), task.get('wp'), task.get('variables'))

# all_tasks = list(itertools.product(pt_bins.get('AK8'), jet_versions.get('AK8'), wps.get('AK8')))
all_tasks = list(itertools.product(['AK8'], pt_bins.get('AK8'), jet_versions.get('AK8'), wps.get('AK8')))
all_tasks.extend(list((itertools.product(['HOTVR'], pt_bins.get('HOTVR'), jet_versions.get('HOTVR'), wps.get('HOTVR')))))

# print all_tasks
# print len(all_tasks)

def set_variables(x):
    vars = list()
    if x[3] == 'NoTau32Cut' and x[1] == pt_bins.get(x[0])[0]:
        vars = variables.get(x[0])
    else:
        if x[0] == 'HOTVR':
            vars.append('mass')
        elif 'AK8' in x[0]: # 'in' because x[0] could also read 'AK8_mSD10'
            vars.append('mSD')
    if len(vars) == 0: sys.exit('Check variables!')
    return vars

all_tasks = [OrderedDict([
    ('probejet_coll', x[0]),
    ('pt_bin', x[1]),
    ('jet_version', x[2]),
    ('wp', x[3]),
    ('variables', set_variables(x)), ### Schreibe wirklich nur die Variablen aus, die Du brauchst!!!!! (entferne den obigen for-loop in create_input_histograms)
]) for x in all_tasks]

print all_tasks

rootFileNames_and_channelNames = list()
def save_rootFileNames_and_channelNames(result):
    rootFileNames_and_channelNames.append(result)

## max_workers = mp.cpu_count()-1
max_workers = 12

pool = mp.Pool(max_workers)
for task in all_tasks:
    pool.apply_async(create_input_histograms_task, args=(task,), callback=save_rootFileNames_and_channelNames)
pool.close()
pool.join()

# sets_of_task_args = list()
# for task in all_tasks:
#     sets_of_task_args.append((task.get('probejet_coll'), task.get('pt_bin'), task.get('jet_version'), task.get('wp'), task.get('variables'), ))
# rootFileNames_and_channelNames = run_with_pool(function=create_input_histograms, sets_of_args=sets_of_task_args)

rootFileNames_and_channelNames_flattened = list()
for pair in rootFileNames_and_channelNames:
    for channelName in pair[1]:
        rootFileNames_and_channelNames_flattened.append([pair[0], channelName])

def write_file_with_plotting_commands(fileName, list_of_histograms):
    workDir = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/Analysis/Combine/workdir')
    os.system('mkdir -p '+workDir)
    filePath = os.path.join(workDir, fileName)
    with open(filePath, 'w') as file:
        for h in list_of_histograms:
            probejet_coll = h[0].split('/')[-1].split('ProbeJetHists_')[-1].replace('.root', '')
            # is_AK8_mSD10 = h[0].endswith('AK8_mSD10.root')
            # is_HOTVR = h[0].endswith('HOTVR.root')
            newline = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/Analysis/Combine/bin', 'plots')+' '
            args = list()
            args.append('Prefit')
            args.append(h[0])
            args.append(h[1])
            args.append(os.path.join('/'.join(h[0].split('/')[0:-2]), 'plots'))
            args.append('__'.join([probejet_coll, h[1], 'prefit.pdf']))
            args.append(probejet_collections.get(probejet_coll)+' PUPPI')
            # args.append('''{}, {:.1f} fb^{#minus1} (13 TeV)'''.format(year.get('short_name'), year.get('lumi_fb')))
            args.append(year.get('short_name')+''', '''+'''{:.1f}'''.format(year.get('lumi_fb'))+''' fb^{#minus1} (13 TeV)''')
            args.append('Work in Progress')
            args.append(str(year.get('lumi_unc')))
            for arg in args:
                arg = re.escape(arg)
                newline += ' '+arg
            newline += '\n'
            file.write(newline)
    print 'Wrote', filePath
# print rootFileNames_and_channelNames_flattened
print len(rootFileNames_and_channelNames_flattened)
write_file_with_plotting_commands('workfile_plotting_commands_prefit_'+year.get('short_name')+'.txt', rootFileNames_and_channelNames_flattened)
