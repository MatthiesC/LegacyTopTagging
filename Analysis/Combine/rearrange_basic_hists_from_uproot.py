#!/usr/bin/env python2

from __future__ import print_function

import os
from subprocess import call

import ROOT as root
from copy import deepcopy
from collections import OrderedDict
import numpy as np

import sys
sys.path.append(os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/Analysis'))
from constants import Systematics, _BANDS, _TAGGERS, _PT_INTERVALS_TANDP_AK8_T, _PT_INTERVALS_TANDP_AK8_W, _PT_INTERVALS_TANDP_HOTVR, _YEARS, get_variable_binning_xlabel_xunit, _NULL_WP

# systematics = Systematics(blacklist=['sfelec', 'sfmu_iso'])
# systematics = Systematics(blacklist=['sfmu_iso'])
systematics = Systematics(blacklist=['sfelec_trigger', 'sfmu_iso'])
# systematics = Systematics(include_jes_splits=False, blacklist=['sfelec_trigger', 'sfmu_iso'])
systs = systematics.base
# print(systs.keys())
# sys.exit()

sys.path.append(os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/NicePlots/python'))
from plotter import NiceStackWithRatio, Process, human_format

from parallel_threading import run_with_pool


# def deepcopy(x):
#     return x


all_years = [
'UL16preVFP',
# 'UL16postVFP',
# 'UL17',
# 'UL18',
]

taggers = [
'ak8_t__tau', #naf11 tmux5
# 'ak8_t_btagDJet__tau', #naf11 tmux6
# 'ak8_t_btagDCSV__tau',#naf11 tmux7
# 'hotvr_t__tau', #naf11 tmux8
# 'ak8_w__partnet',
# 'ak8_t__MDdeepak8',
]
taggers = {k: _TAGGERS[k] for k in taggers}

# # set working points:
# for tagger_k, tagger in taggers.items():
#     if tagger_k == 'ak8_t__tau' or tagger_k == 'ak8_t_btagDJet__tau' or tagger_k == 'ak8_t_btagDCSV__tau' or tagger_k == 'ak8_t_btagDCSV__tau' or tagger_k == 'ak8_t__MDdeepak8':
#         # tagger.wps = [
#         #     WorkingPoint(0.001, 0.38),
#         #     WorkingPoint(0.005, 0.47),
#         #     WorkingPoint(0.010, 0.52),
#         #     WorkingPoint(0.025, 0.61),
#         #     WorkingPoint(0.050, 0.69),
#         # ]
#         tagger.var_intervals = _PT_INTERVALS_TANDP_AK8_T
#         tagger.fit_variable = 'output_probejet_AK8_mSD'
#     if tagger_k == 'hotvr_t__tau':
#         # hotvr_wp = WorkingPoint(9.999, 0.56) # providing bkg_eff does not make sense here; I have not derived it on my own; should be 3% or so
#         # hotvr_wp.name = 'Standard'
#         # tagger.wps = [hotvr_wp]
#         tagger.var_intervals = _PT_INTERVALS_TANDP_HOTVR
#         tagger.fit_variable = 'output_probejet_HOTVR_mass'
#     if tagger_k == 'ak8_w__partnet':
#         # tagger.wps = [
#         #     WorkingPoint(0.030, 0.00), # cut value year dependent! Need to hack it in the code below where we define the cut rules...
#         # ]
#         tagger.var_intervals = _PT_INTERVALS_TANDP_AK8_W
#         tagger.fit_variable = 'output_probejet_AK8_mSD'

processes = [
# 'DATA', # treated as special case in for-loops below
'QCD',
# 'nonQCD',
'TTbar',
'TTbar__MSc_FullyMerged',
'TTbar__MSc_WMerged',
'TTbar__MSc_QBMerged',
'TTbar__MSc_NotMerged',
'TTbar__MSc_YllufMerged',
'TTbar__MSc_SemiMerged',
'TTbar__MSc_NotTopOrWMerged',
'TTbar__MSc_Background',
'ST',
'ST__MSc_FullyMerged',
'ST__MSc_WMerged',
'ST__MSc_QBMerged',
'ST__MSc_NotMerged',
'ST__MSc_YllufMerged',
'ST__MSc_SemiMerged',
'ST__MSc_NotTopOrWMerged',
'ST__MSc_Background',
'VJetsAndVV',
]

# do_histograms = True
do_histograms = False

### only relevant for plots:
do_plotting = True
# do_plotting = False
do_legend = True
# do_legend = False
# mscSplitting = 'mscNone'
# mscSplitting = 'mscTop2'
mscSplitting = 'mscTop3'
# mscSplitting = 'mscW3'

processes_Plotter = None

if mscSplitting == 'mscNone':

    processes_Plotter = [
    Process('VJetsAndVV', 'V+jets, VV', root.kSpring-3),
    Process('ST', 'Single t', root.kOrange),
    Process('TTbar', 't#bar{t}', root.kPink-3),
    Process('QCD', 'Multijet', root.kAzure+10),
    ]

elif mscSplitting == 'mscTop2':

    processes_Plotter = [
    Process('VJetsAndVV', 'V+jets, VV', root.kSpring-3),
    # Process('ST', 'Single top', root.kYellow),
    # Process('TTbar', 't#bar{t}', root.kRed),
    Process('ST__MSc_Background', 'Background', root.kOrange+4),
    Process('ST__MSc_YllufMerged', 'Not merged', root.kOrange-3),
    Process('ST__MSc_FullyMerged', 'Fully merged', root.kOrange),
    Process('TTbar__MSc_Background', 'Background', root.kPink-7),
    Process('TTbar__MSc_YllufMerged', 'Not merged', root.kPink+4),
    Process('TTbar__MSc_FullyMerged', 'Fully merged', root.kPink-3),
    Process('QCD', 'Multijet', root.kAzure+10),
    ]

elif mscSplitting == 'mscTop3':

    processes_Plotter = [
    Process('VJetsAndVV', 'V+jets, VV', root.kSpring-3),
    # Process('ST', 'Single top', root.kYellow),
    # Process('TTbar', 't#bar{t}', root.kRed),
    Process('ST__MSc_Background', 'Background', root.kOrange+4),
    Process('ST__MSc_NotMerged', 'Not merged', root.kOrange-3),
    # Process('ST__MSc_SemiMerged', 'Semi-merged', root.kOrange+6),
    Process('ST__MSc_SemiMerged', 'Semi-merged', root.kYellow),
    Process('ST__MSc_FullyMerged', 'Fully merged', root.kOrange),
    Process('TTbar__MSc_Background', 'Background', root.kPink-7),
    Process('TTbar__MSc_NotMerged', 'Not merged', root.kPink+4),
    Process('TTbar__MSc_SemiMerged', 'Semi-merged', root.kPink+6),
    Process('TTbar__MSc_FullyMerged', 'Fully merged', root.kPink-3),
    Process('QCD', 'Multijet', root.kAzure+10),
    ]

elif mscSplitting == 'mscW3':

    processes_Plotter = [
    Process('VJetsAndVV', 'V+jets, VV', root.kSpring-3),
    # Process('ST', 'Single top', root.kYellow),
    # Process('TTbar', 't#bar{t}', root.kRed),
    Process('ST__MSc_Background', 'Background', root.kOrange+4),
    Process('ST__MSc_NotTopOrWMerged', 'Not merged', root.kOrange-3),
    # Process('ST__MSc_WMerged', 'W merged', root.kOrange+6),
    Process('ST__MSc_WMerged', 'W merged', root.kYellow),
    Process('ST__MSc_FullyMerged', 't merged', root.kOrange),
    Process('TTbar__MSc_Background', 'Background', root.kPink-7),
    Process('TTbar__MSc_NotTopOrWMerged', 'Not merged', root.kPink+4),
    Process('TTbar__MSc_WMerged', 'W merged', root.kPink+6),
    Process('TTbar__MSc_FullyMerged', 't merged', root.kPink-3),
    Process('QCD', 'Multijet', root.kAzure+10),
    ]


# processes_Plotter = OrderedDict([
#     (p.name, p) for p in processes_Plotter
# ])
processes_Plotter_temp = OrderedDict()
for index, proc in enumerate(processes_Plotter):
    proc.index = index
    processes_Plotter_temp[proc.name] = proc
processes_Plotter = processes_Plotter_temp

systs_Plotter = []
for syst in systs.values():
    if syst.tandp:
        systs_Plotter.append(syst.combine_name)

regions = [
    'Pass',
    'Fail',
]

channels = [
'muo',
'ele',
]

bands = _BANDS


class DummyTH1():

    '''
    Saves the crucial information of a TH1 without the need to save a TH1 object directly.
    Can recreate a TH1 with the saved properties.

    This helper class is needed because PyROOT has some nice features called Segmentation Violations which will appear out of nowhere if too many TH1 instances are in memory at once...
    '''

    def __init__(self, th1, new_name=None, delete_orig=False):

        self.name = new_name or th1.GetName()
        self.title = new_name or th1.GetTitle()
        self.n_bins = th1.GetNbinsX()
        self.binning = np.array([], np.float32)
        self.bin_contents = np.array([], np.float32)
        self.bin_errors = np.array([], np.float32)

        for i_bin in range(1, self.n_bins + 1):
            self.binning = np.append(self.binning, [ th1.GetXaxis().GetBinLowEdge(i_bin) ])
        self.binning = np.append(self.binning, [ th1.GetXaxis().GetBinUpEdge(self.n_bins) ])
        for i_bin in range(0, self.n_bins + 2):
            self.bin_contents = np.append(self.bin_contents, [ th1.GetBinContent(i_bin) ])
            self.bin_errors = np.append(self.bin_errors, [ th1.GetBinError(i_bin) ])

        if delete_orig:
            th1.Delete()

    def recreate_th1(self):

        hnew = root.TH1F(self.name, self.title, self.n_bins, self.binning)
        for i_bin in range(0, self.n_bins + 2):
            hnew.SetBinContent(i_bin, self.bin_contents[i_bin])
            hnew.SetBinError(i_bin, self.bin_errors[i_bin])

        return hnew


def hist_preprocessing(hist, binning=None):
    if binning is not None:
        hist = hist.Rebin(len(binning)-1, 'hnew', binning)
    for i in range(hist.GetNbinsX()+2):
        if hist.GetBinContent(i) < 0:
            hist.SetBinContent(i, 0)
            hist.SetBinError(i, 0)
    return hist


def rebin_to_remove_empty_bins(hists_to_write, also_check_data=False):
    '''
    If one histogram in this dict has a bin with zero MC, we will rebin all hists such that this is not the case anymore.\n
    Actually, we only need to check the nominal hists and data_obs for zero bins, but of course then rebin all hists in the dict.
    param: also_check_data: Whether to adjust the binning also if data is zero in a bin (should be fine not to do this; only really important that total MC prediction is > 0 in each bin)
    '''

    # create total MC stack
    hists_nominal = OrderedDict()
    hist_data_obs = None
    for hist_name, hist in hists_to_write.items():
        if hist_name != 'data_obs' and not (hist_name.endswith('Up') or hist_name.endswith('Down') or '__MSc_' in hist_name):
            # hists_nominal[hist_name] = hist.Clone(hist_name+'_cloned')
            hists_nominal[hist_name] = hist.recreate_th1()
        if hist_name == 'data_obs':
            hist_data_obs = hist.recreate_th1()
    hist_total_stack_nominal = None
    for hist_name, hist in hists_nominal.items():
        if hist_total_stack_nominal is None:
            hist_total_stack_nominal = hist.Clone('hist_total_stack_nominal')
        else:
            hist_total_stack_nominal.Add(hist)
        hist.Delete()

    # get original binning and find indices of empty bins
    # hist_data_obs = hists_to_write['data_obs']
    binning = np.array([], dtype=np.float32)
    bin_contents_mc = np.array([], dtype=np.float32)
    bin_contents_data = np.array([], dtype=np.float32)
    empty_bins = []
    n_bins = hist_data_obs.GetNbinsX()
    for i_bin in range(1, n_bins + 1):
        binning = np.append(binning, [ hist_data_obs.GetXaxis().GetBinLowEdge(i_bin) ])
        bin_content_mc = hist_total_stack_nominal.GetBinContent(i_bin)
        bin_contents_mc = np.append(bin_contents_mc, [ bin_content_mc ])
        empty = bin_content_mc <= 0
        if also_check_data:
            bin_content_data = hist_data_obs.GetBinContent(i_bin)
            bin_contents_data = np.append(bin_contents_data, [ bin_content_data ])
            empty = empty and bin_content_data <= 0
        if empty:
            empty_bins.append(i_bin)
    binning = np.append(binning, [ hist_data_obs.GetXaxis().GetBinUpEdge(n_bins) ])
    # print('binning before:', binning)
    # print('bin contents MC:', bin_contents_mc)
    # if also_check_data:
    #     print('bin contents data:', bin_contents_data)
    # print('empty bins:', empty_bins)

    # remove bins with zero content by popping bin edges
    if len(empty_bins) and empty_bins[-1] == len(binning) - 1:
        empty_bins = empty_bins[:-1] # remove last element since we don't want to remove upper bin edge of last bin
    binning = np.delete(binning, empty_bins) # pop upper edges of each empty bin (except if last bin was empty)
    # print(binning)

    # rebin
    hist_total_stack_nominal_rebin = hist_total_stack_nominal.Rebin(len(binning)-1, hist_total_stack_nominal.GetTitle()+'_rebin', binning)
    hist_total_stack_nominal = deepcopy(hist_total_stack_nominal_rebin)
    hist_total_stack_nominal_rebin.Delete()
    hist_data_obs_rebin = hist_data_obs.Rebin(len(binning)-1, hist_data_obs.GetTitle()+'_rebin', binning)
    hist_data_obs = deepcopy(hist_data_obs_rebin)
    hist_data_obs_rebin.Delete()
    last_bin_empty = hist_total_stack_nominal.GetBinContent(len(binning)-1) <= 0
    if also_check_data:
        last_bin_empty = last_bin_empty and hist_data_obs.GetBinContent(len(binning)-1) <= 0
    if last_bin_empty and len(binning) > 2:
        binning = np.delete(binning, len(binning)-2)
    # print(binning)

    # rebin again
    # hist_total_stack_nominal = hist_total_stack_nominal.Rebin(len(binning)-1, hist_total_stack_nominal.GetTitle(), binning)
    # hist_data_obs = hist_data_obs.Rebin(len(binning)-1, hist_data_obs.GetTitle(), binning)
    hist_total_stack_nominal_rebin = hist_total_stack_nominal.Rebin(len(binning)-1, hist_total_stack_nominal.GetTitle()+'_rebin', binning)
    hist_total_stack_nominal = deepcopy(hist_total_stack_nominal_rebin)
    hist_total_stack_nominal_rebin.Delete()
    hist_data_obs_rebin = hist_data_obs.Rebin(len(binning)-1, hist_data_obs.GetTitle()+'_rebin', binning)
    hist_data_obs = deepcopy(hist_data_obs_rebin)
    hist_data_obs_rebin.Delete()
    # sanity check if all bins are non-empty
    for i_bin in range(1, hist_data_obs.GetNbinsX() + 1):
        empty_mc = hist_total_stack_nominal.GetBinContent(i_bin) <= 0
        empty_data = hist_data_obs.GetBinContent(i_bin) <= 0
        if empty_mc:
            sys.exit('Rebinning to remove empty bins did not succeed for MC')
        if empty_data and also_check_data:
            sys.exit('Rebinning to remove empty bins did not succeed for data')

    # rebin everything in the original dict
    # result = OrderedDict()
    # for hist_name, hist in hists_to_write.items():
    #     # result[hist_name] = hist.Rebin(len(binning)-1, hist.GetTitle(), binning)
    #     result[hist_name] = DummyTH1(hist.recreate_th1().Rebin(len(binning)-1, hist_name, binning))


    for hist_name, hist in hists_to_write.items():
        th1 = hist.recreate_th1()
        th1_rebinned = th1.Rebin(len(binning)-1, hist_name+'_rebinned', binning)
        th1.Delete()
        hist = DummyTH1(th1_rebinned)
        th1_rebinned.Delete()

    # for hist in hists_nominal.values():
    #     hist.Delete()
    hist_total_stack_nominal.Delete()
    hist_data_obs.Delete()

    # return result
    # print('binning after:', binning)
    return binning




def create_rearranged_hists(variable, tagger, year, wp, pt_bin, do_plot=False, do_hists=True):

    is_fit_template = not wp.null

    baseDir = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel', year, 'combine', tagger.name)
    inDir = os.path.join(baseDir, 'BasicHists')
    task_name = '-'.join([tagger.name, wp.name, pt_bin.name, year, variable])
    outDir = os.path.join(baseDir, task_name)
    os.system('mkdir -p '+outDir)
    outFileName = 'Templates-'+task_name+'.root'
    outFilePath = os.path.join(outDir, outFileName)
    if do_hists:
        outFile = root.TFile.Open(outFilePath, 'RECREATE')

    outFolderNames = []

    for band in bands.values():

        for region in regions:
            if region == 'Fail' and not is_fit_template:
                continue

            for channel in channels:

                inFolderName = '/'.join([band.name, region, channel])
                outFolderName = '_'.join([band.name, region, channel])
                outFolderNames.append(outFolderName)
                print(inFolderName)
                # sys.exit()

                if not do_hists:
                    continue
                outFolder = outFile.mkdir(outFolderName)

                # hists_to_write = OrderedDict()
                # hists_to_write = []
                # histnames_to_write = []
                hists_nominal = OrderedDict()

                # data
                inFileName_nominal = '-'.join(['BasicHists', tagger.name, wp.name, pt_bin.name, year, 'nominal', variable])+'.root'
                inFilePath_nominal = os.path.join(inDir, inFileName_nominal)
                inFile_nominal = root.TFile.Open(inFilePath_nominal, 'READ')
                inHistPath_data = os.path.join(inFolderName, 'DATA__nominal')
                inHist_data = inFile_nominal.Get(inHistPath_data)
                # inHist_data.SetDirectory(0)
                # inHist_data = hist_preprocessing(inHist_data)
                # outFolder.cd()
                # inHist_data.Write('data_obs')
                # hists_to_write['data_obs'] = deepcopy(inHist_data)
                # hists_to_write['data_obs'] = inHist_data.Clone('data_obs')
                # hists_to_write['data_obs'] = DummyTH1(inHist_data, 'data_obs')
                # histnames_to_write.append('data_obs')
                # hists_to_write.append(inHist_data.Clone("clone"))
                hists_nominal['data_obs'] = DummyTH1(inHist_data)

                # MC
                for process in processes:

                    # nominal
                    inHistPath_nominal = os.path.join(inFolderName, process+'__nominal')
                    inHist_nominal = inFile_nominal.Get(inHistPath_nominal)
                    hists_nominal[process] = DummyTH1(inHist_nominal)

                if tagger.fit_variable == variable:
                    rebinning_scheme = rebin_to_remove_empty_bins(hists_nominal) # HACK: Avoid this, only leads to bugs over bugs over bugs ....!!!!!
                    # rebinning_scheme = None
                else:
                    rebinning_scheme = None

                inHist_data = hist_preprocessing(inHist_data, binning=rebinning_scheme)
                outFolder.cd()
                inHist_data.Write('data_obs')

                for process in processes:

                    # nominal
                    inHistPath_nominal = os.path.join(inFolderName, process+'__nominal')
                    inHist_nominal = inFile_nominal.Get(inHistPath_nominal)
                    # inHist_nominal.SetDirectory(0)
                    print(inHistPath_nominal)
                    inHist_nominal = hist_preprocessing(inHist_nominal, binning=rebinning_scheme)
                    outFolder.cd()
                    inHist_nominal.Write(process)
                    # hists_to_write[process] = deepcopy(inHist_nominal)
                    # hists_to_write[process] = inHist_nominal.Clone(process)
                    # hists_to_write[process] = DummyTH1(inHist_nominal, process)
                    # histnames_to_write.append(process)
                    # hists_to_write.append(inHist_nominal.Clone("clone"))

                    # systematics which do not need special treatment
                    for syst_k in sorted(systs.keys()):
                        syst = systs[syst_k]

                        # print(syst_k) #######

                        if syst.name in ['cr', 'murmuf']: continue

                        inHists_variations = {}
                        at_least_one_integral_is_zero = inHist_nominal.Integral() <= 0

                        for variation_k in sorted(syst.variations.keys()):
                            variation = syst.variations[variation_k]

                            inFileName_variation = '-'.join(['BasicHists', tagger.name, wp.name, pt_bin.name, year, variation.name, variable])+'.root'
                            inFilePath_variation = os.path.join(inDir, inFileName_variation)
                            inFile_variation = root.TFile.Open(inFilePath_variation, 'READ')
                            inHistPath_variation = os.path.join(inFolderName, process+'__'+variation.name)
                            inHist_variation = inFile_variation.Get(inHistPath_variation)
                            # inHist_variation.SetDirectory(0)
                            inHist_variation = hist_preprocessing(inHist_variation, binning=rebinning_scheme)
                            if inHist_variation.Integral() <= 0:
                                at_least_one_integral_is_zero = True
                            if variation.short_name in ['up', 'mtop173p5']:
                                direction = 'Up'
                            elif variation.short_name in ['down', 'mtop171p5']:
                                direction = 'Down'
                            else:
                                sys.exit('In this code part, only systematics with simple up/down variations should be processed')
                            outHistName_variation = process+'_'+syst.combine_name+direction
                            inHists_variations[outHistName_variation] = deepcopy(inHist_variation)
                            inFile_variation.Close()

                        outFolder.cd()
                        for outHistName_variation, inHist_variation in inHists_variations.items():
                            if at_least_one_integral_is_zero:
                                inHist_nominal.Write(outHistName_variation)
                                # hists_to_write[outHistName_variation] = deepcopy(inHist_nominal)
                                # hists_to_write[outHistName_variation] = inHist_nominal.Clone(outHistName_variation)
                                # hists_to_write[outHistName_variation] = DummyTH1(inHist_nominal, outHistName_variation)
                                # histnames_to_write.append(outHistName_variation)
                                # hists_to_write.append(inHist_nominal.Clone("clone"))
                            else:
                                inHist_variation.Write(outHistName_variation)
                                # hists_to_write[outHistName_variation] = deepcopy(inHist_variation)
                                # hists_to_write[outHistName_variation] = inHist_variation.Clone(outHistName_variation)
                                # hists_to_write[outHistName_variation] = DummyTH1(inHist_variation, outHistName_variation)
                                # histnames_to_write.append(outHistName_variation)
                                # hists_to_write.append(inHist_variation.Clone("clone"))

                    # scale variation (muR/muF) and color reconnection
                     ## create envelope of all available variations and nominal
                    for syst_k in ['murmuf', 'cr']:
                        syst = systs[syst_k]
                        inHists_variations = {}

                        for variation_k in sorted(syst.variations.keys()):
                            variation = syst.variations[variation_k]

                            inFileName_variation = '-'.join(['BasicHists', tagger.name, wp.name, pt_bin.name, year, variation.name, variable])+'.root'
                            inFilePath_variation = os.path.join(inDir, inFileName_variation)
                            inFile_variation = root.TFile.Open(inFilePath_variation, 'READ')
                            inHistPath_variation = os.path.join(inFolderName, process+'__'+variation.name)
                            inHist_variation = inFile_variation.Get(inHistPath_variation)
                            # inHist_variation.SetDirectory(0)
                            inHist_variation = hist_preprocessing(inHist_variation, binning=rebinning_scheme)
                            inHists_variations[variation.name] = deepcopy(inHist_variation)

                        outHist_up = deepcopy(inHist_nominal)
                        outHist_down = deepcopy(inHist_nominal)
                        for i in range(inHist_nominal.GetNbinsX()+2):
                            bin_contents = [inHist_nominal.GetBinContent(i)]
                            for inHist_variation in inHists_variations.values():
                                bin_contents.append(inHist_variation.GetBinContent(i))
                            outHist_up.SetBinContent(i, max(bin_contents))
                            outHist_up.SetBinError(i, 0)
                            outHist_down.SetBinContent(i, min(bin_contents))
                            outHist_down.SetBinError(i, 0)

                        outFolder.cd()
                        outHist_up.Write(process+'_'+syst.combine_name+'Up')
                        outHist_down.Write(process+'_'+syst.combine_name+'Down')
                        # hists_to_write[process+'_'+syst.combine_name+'Up'] = deepcopy(outHist_up)
                        # hists_to_write[process+'_'+syst.combine_name+'Down'] = deepcopy(outHist_down)
                        # hists_to_write[process+'_'+syst.combine_name+'Up'] = outHist_up.Clone(process+'_'+syst.combine_name+'Up')
                        # hists_to_write[process+'_'+syst.combine_name+'Down'] = outHist_down.Clone(process+'_'+syst.combine_name+'Down')
                        # hists_to_write[process+'_'+syst.combine_name+'Up'] = DummyTH1(outHist_up, process+'_'+syst.combine_name+'Up')
                        # hists_to_write[process+'_'+syst.combine_name+'Down'] = DummyTH1(outHist_up, process+'_'+syst.combine_name+'Down')
                        # histnames_to_write.append(process+'_'+syst.combine_name+'Up')
                        # hists_to_write.append(outHist_up.Clone("clone"))
                        # histnames_to_write.append(process+'_'+syst.combine_name+'Down')
                        # hists_to_write.append(outHist_down.Clone("clone"))

                # if tagger.fit_variable == variable:
                #     # hists_to_write = rebin_to_remove_empty_bins(hists_to_write)
                #     rebin_to_remove_empty_bins(hists_to_write)
                # outFolder.cd()
                # for hist_name, hist in hists_to_write.items():
                #     # hist.Write(hist_name)
                #     th1 = hist.recreate_th1()
                #     th1.Write(hist_name)
                #     th1.Delete()
                #
                # # hists_to_write = rebin_to_remove_empty_bins(hists_to_write, histnames_to_write)
                # # outFolder.cd()
                # # for i in range(len(hists_to_write)):
                # #     hists_to_write[i].Write(histnames_to_write[i])

                inFile_nominal.Close()


    if do_hists:
        outFile.Close()
        print('Wrote', outFilePath)

    if do_plot:

        for outFolderName in outFolderNames:

            region = ''
            channel = ''
            for r in regions:
                if region == 'Fail' and not is_fit_template:
                    continue
                if r in outFolderName.split('_'):
                    region = r
                    for c in channels:
                        if c in outFolderName.split('_'):
                            channel = c

            # if variable.endswith('_mSD'):
            #     x_axis_title = 'Probe jet #it{m}_{SD}'
            #     x_axis_unit = 'GeV'
            # elif variable.endswith('_mass'):
            #     x_axis_title = 'Probe jet #it{m}_{jet}'
            #     x_axis_unit = 'GeV'


            if tagger.name.startswith('ak8_'):
                probejetalgo = 'AK8'
            elif tagger.name.startswith('hotvr_'):
                probejetalgo = 'HOTVR'
            _, x_axis_title, x_axis_unit, logy, leg_offset_x, leg_offset_y = get_variable_binning_xlabel_xunit(variable_name=variable.split('_'+probejetalgo+'_')[-1], tagger_name=tagger.name, fit_variable=is_fit_template)

            nice = NiceStackWithRatio(
                infile_path = outFilePath,
                infile_directory = outFolderName, # the directory within the ROOT file
                x_axis_title = x_axis_title,
                x_axis_unit = x_axis_unit,
                prepostfit = 'prefitRaw',
                processes = processes_Plotter.values(),
                syst_names = systs_Plotter,
                lumi_unc = _YEARS.get(year).get('lumi_unc'),
                # divide_by_bin_width = False,
                data_name = 'data_obs',
                # text_prelim = 'Private Work',
                text_prelim = 'Preliminary',
                # text_top_left = _YEARS.get(year).get('long_name'),
                text_top_left = 'T&P '+('e' if channel == 'ele' else '#mu')+'+jets, UL '+_YEARS.get(year).get('year'),
                text_top_right = _YEARS.get(year).get('lumi_fb_display')+' fb^{#minus1} (13 TeV)',
                # nostack = True,
                logy = logy,
            )

            plotName = 'Plot-prefitRaw-'+task_name+'-'+outFolderName+'-'+mscSplitting+('-withLeg' if do_legend else '-noLeg')+('-Preliminary' if nice.text_prelim == 'Preliminary' else '')+'.pdf'
            plotDir = os.path.join(outDir, 'plots')

            nice.plot()

            #________________
            # Some text labels

            nice.canvas.cd()

            # print(tagger.label)
            pt_text_string = tagger.label.replace('{WP_VALUE}', '{}'.format(wp.get_cut_value(year)))+' [#bf{'+region+'}]'
            # pt_text_string2 = 'T&P '+('e' if channel == 'ele' else '#mu')+'+jets, '
            pt_text_string2 = ''
            if pt_bin.max_set:
                pt_text_string2 += '#it{p}_{T} #in ({pt_min}, {pt_max}] GeV, |#eta| < 2.5'.replace('{pt_min}', '{}'.format(pt_bin.var_min)).replace('{pt_max}', '{}'.format(pt_bin.var_max))
            else:
                pt_text_string2 += '#it{p}_{T} > {pt_min} GeV, |#eta| < 2.5'.replace('{pt_min}', '{}'.format(pt_bin.var_min))

            # remove the line with the information regarding the tagger/wp definition and move the pt/eta line one line up
            if not is_fit_template:
                pt_text_string = pt_text_string2

            tlatex_pt = root.TLatex(nice.coord.graph_to_pad_x(0.95), nice.coord.graph_to_pad_y(0.84), pt_text_string if is_fit_template else pt_text_string2)
            tlatex_pt.SetTextAlign(31) # left top
            tlatex_pt.SetTextFont(42)
            tlatex_pt.SetTextSize(0.025)
            tlatex_pt.SetNDC()
            tlatex_pt.Draw()

            if is_fit_template:
                tlatex_pt2 = root.TLatex(nice.coord.graph_to_pad_x(0.95), nice.coord.graph_to_pad_y(0.78), pt_text_string2)
                tlatex_pt2.SetTextAlign(31) # left top
                tlatex_pt2.SetTextFont(42)
                tlatex_pt2.SetTextSize(0.025)
                tlatex_pt2.SetNDC()
                tlatex_pt2.Draw()

            tlatex_probejetalgo = root.TLatex(nice.coord.graph_to_pad_x(0.95), nice.coord.graph_to_pad_y(0.95), probejetalgo+' PUPPI')
            tlatex_probejetalgo.SetTextAlign(33) # right top
            tlatex_probejetalgo.SetTextFont(62)
            tlatex_probejetalgo.SetTextSize(0.05)
            tlatex_probejetalgo.SetTextColor(root.kGray+1)
            tlatex_probejetalgo.SetNDC()
            tlatex_probejetalgo.Draw()

            #________________
            # Legend

            if do_legend:

                if mscSplitting == 'mscNone':

                    legend = root.TLegend(nice.coord.graph_to_pad_x(0.45+leg_offset_x), nice.coord.graph_to_pad_y(0.5+leg_offset_y), nice.coord.graph_to_pad_x(0.7+leg_offset_x), nice.coord.graph_to_pad_y(0.73+leg_offset_y))
                    legend.AddEntry(nice.data_hist, 'Data', 'ep')
                    legend.AddEntry(nice.stack.GetStack().At(processes_Plotter['TTbar'].index), processes_Plotter['TTbar'].legend, 'f')
                    legend.AddEntry(nice.stack.GetStack().At(processes_Plotter['ST'].index), processes_Plotter['ST'].legend, 'f')
                    legend.SetTextSize(0.025)
                    legend.SetBorderSize(0)
                    legend.SetFillStyle(0)
                    legend.Draw()

                    legend2 = root.TLegend(nice.coord.graph_to_pad_x(0.45+0.22+leg_offset_x), nice.coord.graph_to_pad_y(0.5+leg_offset_y), nice.coord.graph_to_pad_x(0.7+0.22+leg_offset_x), nice.coord.graph_to_pad_y(0.73+leg_offset_y))
                    legend2.AddEntry(nice.stack.GetStack().At(processes_Plotter['VJetsAndVV'].index), processes_Plotter['VJetsAndVV'].legend, 'f')
                    legend2.AddEntry(nice.stack.GetStack().At(processes_Plotter['QCD'].index), processes_Plotter['QCD'].legend, 'f')
                    legend2.AddEntry(nice.stack_unc, 'Uncertainty', 'f')
                    legend2.SetTextSize(0.025)
                    legend2.SetBorderSize(0)
                    legend2.SetFillStyle(0)
                    legend2.Draw()

                else:

                    legend = root.TLegend(nice.coord.graph_to_pad_x(0.45+leg_offset_x), nice.coord.graph_to_pad_y(0.5+leg_offset_y), nice.coord.graph_to_pad_x(0.7+leg_offset_x), nice.coord.graph_to_pad_y(0.73+leg_offset_y))
                    if mscSplitting == 'mscTop3' or mscSplitting == 'mscW3':
                        legend.SetHeader('')
                    legend.AddEntry(nice.data_hist, 'Data', 'ep')
                    legend.AddEntry(nice.stack.GetStack().At(processes_Plotter['VJetsAndVV'].index), processes_Plotter['VJetsAndVV'].legend, 'f')
                    legend.AddEntry(nice.stack.GetStack().At(processes_Plotter['QCD'].index), processes_Plotter['QCD'].legend, 'f')
                    legend.AddEntry(nice.stack_unc, 'Uncertainty', 'f')
                    legend.SetTextSize(0.025)
                    legend.SetBorderSize(0)
                    legend.SetFillStyle(0)
                    legend.Draw()

                    legend2 = root.TLegend(nice.coord.graph_to_pad_x(0.45+0.22+leg_offset_x), nice.coord.graph_to_pad_y(0.5+leg_offset_y), nice.coord.graph_to_pad_x(0.7+0.22+leg_offset_x), nice.coord.graph_to_pad_y(0.73+leg_offset_y))
                    legend2.SetHeader('t#bar{t} / single t categories:')
                    if mscSplitting == 'mscTop2':
                        legend2.AddEntry(nice.stack.GetStack().At(processes_Plotter['TTbar__MSc_FullyMerged'].index), '/', 'f')
                        legend2.AddEntry(nice.stack.GetStack().At(processes_Plotter['TTbar__MSc_YllufMerged'].index), '/', 'f')
                        legend2.AddEntry(nice.stack.GetStack().At(processes_Plotter['TTbar__MSc_Background'].index), '/', 'f')
                    elif mscSplitting == 'mscTop3':
                        legend2.AddEntry(nice.stack.GetStack().At(processes_Plotter['TTbar__MSc_FullyMerged'].index), '/', 'f')
                        legend2.AddEntry(nice.stack.GetStack().At(processes_Plotter['TTbar__MSc_SemiMerged'].index), '/', 'f')
                        legend2.AddEntry(nice.stack.GetStack().At(processes_Plotter['TTbar__MSc_NotMerged'].index), '/', 'f')
                        legend2.AddEntry(nice.stack.GetStack().At(processes_Plotter['TTbar__MSc_Background'].index), '/', 'f')
                    elif mscSplitting == 'mscW3':
                        legend2.AddEntry(nice.stack.GetStack().At(processes_Plotter['TTbar__MSc_FullyMerged'].index), '/', 'f')
                        legend2.AddEntry(nice.stack.GetStack().At(processes_Plotter['TTbar__MSc_WMerged'].index), '/', 'f')
                        legend2.AddEntry(nice.stack.GetStack().At(processes_Plotter['TTbar__MSc_NotTopOrWMerged'].index), '/', 'f')
                        legend2.AddEntry(nice.stack.GetStack().At(processes_Plotter['TTbar__MSc_Background'].index), '/', 'f')
                    legend2.SetTextSize(0.025)
                    legend2.SetBorderSize(0)
                    legend2.SetFillStyle(0)
                    legend2.Draw()

                    legend3 = root.TLegend(nice.coord.graph_to_pad_x(0.45+0.22+0.067+leg_offset_x), nice.coord.graph_to_pad_y(0.5+leg_offset_y), nice.coord.graph_to_pad_x(0.7+0.22+0.067+leg_offset_x), nice.coord.graph_to_pad_y(0.73+leg_offset_y))
                    legend3.SetHeader('')
                    if mscSplitting == 'mscTop2':
                        legend3.AddEntry(nice.stack.GetStack().At(processes_Plotter['ST__MSc_FullyMerged'].index), processes_Plotter['ST__MSc_FullyMerged'].legend, 'f')
                        legend3.AddEntry(nice.stack.GetStack().At(processes_Plotter['ST__MSc_YllufMerged'].index), processes_Plotter['ST__MSc_YllufMerged'].legend, 'f')
                        legend3.AddEntry(nice.stack.GetStack().At(processes_Plotter['ST__MSc_Background'].index), processes_Plotter['ST__MSc_Background'].legend, 'f')
                    elif mscSplitting == 'mscTop3':
                        legend3.AddEntry(nice.stack.GetStack().At(processes_Plotter['ST__MSc_FullyMerged'].index), processes_Plotter['ST__MSc_FullyMerged'].legend, 'f')
                        legend3.AddEntry(nice.stack.GetStack().At(processes_Plotter['ST__MSc_SemiMerged'].index), processes_Plotter['ST__MSc_SemiMerged'].legend, 'f')
                        legend3.AddEntry(nice.stack.GetStack().At(processes_Plotter['ST__MSc_NotMerged'].index), processes_Plotter['ST__MSc_NotMerged'].legend, 'f')
                        legend3.AddEntry(nice.stack.GetStack().At(processes_Plotter['ST__MSc_Background'].index), processes_Plotter['ST__MSc_Background'].legend, 'f')
                    elif mscSplitting == 'mscW3':
                        legend3.AddEntry(nice.stack.GetStack().At(processes_Plotter['ST__MSc_FullyMerged'].index), processes_Plotter['ST__MSc_FullyMerged'].legend, 'f')
                        legend3.AddEntry(nice.stack.GetStack().At(processes_Plotter['ST__MSc_WMerged'].index), processes_Plotter['ST__MSc_WMerged'].legend, 'f')
                        legend3.AddEntry(nice.stack.GetStack().At(processes_Plotter['ST__MSc_NotTopOrWMerged'].index), processes_Plotter['ST__MSc_NotTopOrWMerged'].legend, 'f')
                        legend3.AddEntry(nice.stack.GetStack().At(processes_Plotter['ST__MSc_Background'].index), processes_Plotter['ST__MSc_Background'].legend, 'f')
                    legend3.SetTextSize(0.025)
                    legend3.SetBorderSize(0)
                    legend3.SetFillStyle(0)
                    legend3.Draw()

            pt_text_string3 = '#minus prefit #minus'
            tlatex_pt3 = root.TLatex(nice.coord.graph_to_pad_x(0.95), nice.coord.graph_to_pad_y(0.4 if is_fit_template else 0.77), pt_text_string3)
            tlatex_pt3.SetTextAlign(31) # left top
            tlatex_pt3.SetTextFont(72)
            tlatex_pt3.SetTextSize(nice.text_size)
            tlatex_pt3.SetNDC()
            tlatex_pt3.Draw()

            nice.save_plot(plotName, plotDir)


if __name__=='__main__':

    #__________________________________________________
    # Default code to create all fitting templates (mSD or jet mass)

    for the_tagger in taggers.values():
        print('Working on', the_tagger.name)
        for year in all_years:
            print('Working on', year)
            for wp in the_tagger.get_wp(year=year):
                print('Working on', wp.name)
                for pt_bin in the_tagger.var_intervals.values():
                    # if pt_bin.name != 'pt_300to400': continue # HACK
                    print('Working on', pt_bin.name)
                    create_rearranged_hists(the_tagger.fit_variable, the_tagger, year, wp, pt_bin, do_plot=do_plotting, do_hists=do_histograms)
                # break # HACK

    # #__________________________________________________
    # # Code to create plots of other variables (e.g. substructure)
    #
    # for the_tagger in taggers.values():
    #
    #     the_vars = []
    #     if the_tagger.name.startswith('ak8_t'):
    #         the_vars = [
    #             'output_probejet_AK8_tau32',
    #             'output_probejet_AK8_maxDeepJet',
    #             'output_probejet_AK8_maxDeepCSV',
    #             'output_probejet_AK8_MDDeepAK8_TvsQCD',
    #             'output_probejet_AK8_mSD',
    #             'output_probejet_AK8_mass',
    #             'output_probejet_AK8_pt',
    #         ]
    #     elif the_tagger.name.startswith('ak8_w'):
    #         the_vars = [
    #             'output_probejet_AK8_ParticleNet_WvsQCD',
    #             'output_probejet_AK8_mSD',
    #             'output_probejet_AK8_mass',
    #             'output_probejet_AK8_pt',
    #         ]
    #     elif the_tagger.name.startswith('hotvr_t'):
    #         the_vars = [
    #             'output_probejet_HOTVR_tau32',
    #             'output_probejet_HOTVR_nsub',
    #             'output_probejet_HOTVR_fpt1',
    #             'output_probejet_HOTVR_mpair',
    #             'output_probejet_HOTVR_mass',
    #             'output_probejet_HOTVR_pt',
    #         ]
    #
    #     print('Working on', the_tagger.name)
    #
    #     for var in the_vars:
    #
    #         print('Working on', var)
    #
    #         for year in all_years:
    #
    #             print('Working on', year)
    #             wp = _NULL_WP
    #
    #             for pt_bin in the_tagger.var_intervals.values():
    #
    #                 if not pt_bin.total_range:
    #                     continue
    #                 print('Working on', pt_bin.name)
    #                 create_rearranged_hists(var, the_tagger, year, wp, pt_bin, do_plot=do_plotting, do_hists=do_histograms)
