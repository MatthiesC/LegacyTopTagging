import os
import sys
import argparse
from tqdm import tqdm

import uproot
import numpy as np

from hist import Hist
import hist

years = ['UL16preVFP', 'UL16postVFP', 'UL17', 'UL18']

if not sys.argv[1:]: sys.exit('No arguments provided. Exit.')
parser = argparse.ArgumentParser()
parser.add_argument('-y', '--years', choices=years, nargs='*', default=[])
parser.add_argument('-t', '--test', action='store_true')
args = parser.parse_args(sys.argv[1:])

branches = [
'event_weight',
'hotvrjets_passes_jet_id',
'hotvrjets_dr_rec_to_gen',
'hotvrjets_reff',
'hotvrjets_pt_gen',
'hotvrjets_pt_rec_raw',
'hotvrjets_pt_rec_corr',
'hotvrjets_pt_rec_corr_jes_up',
'hotvrjets_pt_rec_corr_jes_down',
'hotvrjets_pt_rec_corr_jer_up',
'hotvrjets_pt_rec_corr_jer_down',
'hotvrjets_eta',
]

aliases = [x.replace('hotvrjets_', '').replace('event_', '') for x in branches]
alias_calculations = {}
for i, alias in enumerate(aliases):
    alias_calculations[alias] = '{0}'.format(branches[i])

binning = [190, 210] # point with center at 200
for i in range(19): binning.append(binning[-1] + 10)
for i in range(20): binning.append(binning[-1] + 15)
for i in range(4): binning.append(binning[-1] + 20)
for i in range(3): binning.append(binning[-1] + 40)
for i in range(2): binning.append(binning[-1] + 75)
for i in range(2): binning.append(binning[-1] + 100)
# for i in range(3): binning.append(binning[-1] + 250)
binning = np.array(binning, dtype=np.float32)
# print(binning)

class Bin():
    def __init__(self, low_edge, high_edge):
        self.low_edge = low_edge
        self.high_edge = high_edge
        self.bin_center = low_edge + 0.5 * (high_edge - low_edge)
        self.response_raw = np.array([], dtype=np.float32)
        self.response_corr = np.array([], dtype=np.float32)
        self.response_corr_jes_up = np.array([], dtype=np.float32)
        self.response_corr_jes_down = np.array([], dtype=np.float32)
        self.response_corr_jer_up = np.array([], dtype=np.float32)
        self.response_corr_jer_down = np.array([], dtype=np.float32)
        self.weights = np.array([], dtype=np.float32)
        # self.response_raw_mean = None
        # self.response_raw_std = None
        # self.response_raw_mc_unc = None
        # self.response_corr_mean = None
        # self.response_corr_std = None
        # self.response_corr_mc_unc = None

bins = []
for i_bin, low_edge in enumerate(binning):
    if i_bin == len(binning)-1: continue
    high_edge = binning[i_bin+1]
    bins.append(Bin(low_edge, high_edge))

base_cut = '(pt_gen > {0}) & (passes_jet_id == True) & (eta < 2.5) & (dr_rec_to_gen < reff)'.format(str(binning[0]))

hist_names = [
'response_gen',
'response_raw',
'response_corr',
'response_corr_jes_up',
'response_corr_jes_down',
'response_corr_jer_up',
'response_corr_jer_down',

'response_raw_mean',
'response_raw_mean_mc_unc',
'response_raw_std',
'response_raw_std_mc_unc',
'response_corr_mean',
'response_corr_mean_mc_unc',
'response_corr_mean_jes_up',
'response_corr_mean_jes_down',
'response_corr_mean_jer_up',
'response_corr_mean_jer_down',
'response_corr_std',
'response_corr_std_mc_unc',
'response_corr_std_jes_up',
'response_corr_std_jes_down',
'response_corr_std_jer_up',
'response_corr_std_jer_down',
]

def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return average, np.sqrt(variance)

def splitted_weighted_stds_with_common_average(values, weights, splits=10):
    split_values = np.array_split(values, splits)
    split_weights = np.array_split(weights, splits)
    average = np.average(values, weights=weights)
    split_stds = np.array([], dtype=np.float32)
    split_weight_sums = np.array([], dtype=np.float32)
    for i in range(splits):
        variance = np.average((split_values[i]-average)**2, weights=split_weights[i])
        split_stds = np.append(split_stds, [ np.sqrt(variance) ])
        split_weight_sums = np.append(split_weight_sums, [ np.sum(split_weights[i]) ])
    return split_stds, split_weight_sums

for year in args.years:

    hists = {}

    for hist_name in hist_names:

        hists[hist_name] = Hist(
            hist.axis.Variable(
                binning, label='generated pt', name='pt_gen', underflow=True, overflow=True
            ),
            storage=hist.storage.Weight(), # like ROOT's Sumw2()
        )

    workDir = '_'.join(['workdir', 'WorkingPointStudy', year])
    inputDir = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/output/WorkingPointStudy', year, 'nominal', workDir)
    inputFileName = 'uhh2.AnalysisModuleRunner.MC.TTbarToHadronic_'+year+'_{0}.root'
    inputFilePath = os.path.join(inputDir, inputFileName)

    if args.test:
        files_i = range(1)
        batches = [inputFilePath.format(str(i))+':AnalysisTree' for i in files_i]
        print(batches)
    else:
        batches = [inputFilePath.format('*')+':AnalysisTree']

    for batch in tqdm(uproot.iterate(batches, expressions=aliases, aliases=alias_calculations, cut=base_cut, library='pd')):

        weight = batch['weight'].to_numpy()
        pt_gen = batch['pt_gen'].to_numpy()
        hists['response_gen'].fill(pt_gen, weight=weight)
        hists['response_raw'].fill(pt_gen, weight=weight * batch['pt_rec_raw'].to_numpy() / pt_gen)
        hists['response_corr'].fill(pt_gen, weight=weight * batch['pt_rec_corr'].to_numpy() / pt_gen)
        hists['response_corr_jes_up'].fill(pt_gen, weight=weight * batch['pt_rec_corr_jes_up'].to_numpy() / pt_gen)
        hists['response_corr_jes_down'].fill(pt_gen, weight=weight * batch['pt_rec_corr_jes_down'].to_numpy() / pt_gen)
        hists['response_corr_jer_up'].fill(pt_gen, weight=weight * batch['pt_rec_corr_jer_up'].to_numpy() / pt_gen)
        hists['response_corr_jer_down'].fill(pt_gen, weight=weight * batch['pt_rec_corr_jer_down'].to_numpy() / pt_gen)

        for bin in bins:
            this_bin = batch[(batch['pt_gen'] > bin.low_edge) & (batch['pt_gen'] < bin.high_edge)]
            this_bin__pt_gen = this_bin['pt_gen'].to_numpy()
            this_bin__response_raw = this_bin['pt_rec_raw'].to_numpy() / this_bin__pt_gen
            this_bin__response_corr = this_bin['pt_rec_corr'].to_numpy() / this_bin__pt_gen
            this_bin__response_corr_jes_up = this_bin['pt_rec_corr_jes_up'].to_numpy() / this_bin__pt_gen
            this_bin__response_corr_jes_down = this_bin['pt_rec_corr_jes_down'].to_numpy() / this_bin__pt_gen
            this_bin__response_corr_jer_up = this_bin['pt_rec_corr_jer_up'].to_numpy() / this_bin__pt_gen
            this_bin__response_corr_jer_down = this_bin['pt_rec_corr_jer_down'].to_numpy() / this_bin__pt_gen
            this_bin__weights = this_bin['weight'].to_numpy()
            bin.response_raw = np.append(bin.response_raw, this_bin__response_raw)
            bin.response_corr = np.append(bin.response_corr, this_bin__response_corr)
            bin.response_corr_jes_up = np.append(bin.response_corr_jes_up, this_bin__response_corr_jes_up)
            bin.response_corr_jes_down = np.append(bin.response_corr_jes_down, this_bin__response_corr_jes_down)
            bin.response_corr_jer_up = np.append(bin.response_corr_jer_up, this_bin__response_corr_jer_up)
            bin.response_corr_jer_down = np.append(bin.response_corr_jer_down, this_bin__response_corr_jer_down)
            bin.weights = np.append(bin.weights, this_bin__weights)

    for bin in bins:
        bin.response_raw_mean, bin.response_raw_std = weighted_avg_and_std(bin.response_raw, bin.weights)
        bin.response_corr_mean, bin.response_corr_std = weighted_avg_and_std(bin.response_corr, bin.weights)
        bin.response_corr_mean_jes_up, bin.response_corr_std_jes_up = weighted_avg_and_std(bin.response_corr_jes_up, bin.weights)
        bin.response_corr_mean_jes_down, bin.response_corr_std_jes_down = weighted_avg_and_std(bin.response_corr_jes_down, bin.weights)
        bin.response_corr_mean_jer_up, bin.response_corr_std_jer_up = weighted_avg_and_std(bin.response_corr_jer_up, bin.weights)
        bin.response_corr_mean_jer_down, bin.response_corr_std_jer_down = weighted_avg_and_std(bin.response_corr_jer_down, bin.weights)
        hists['response_raw_mean'].fill([bin.bin_center], weight=[bin.response_raw_mean])
        hists['response_raw_std'].fill([bin.bin_center], weight=[bin.response_raw_std])
        hists['response_corr_mean'].fill([bin.bin_center], weight=[bin.response_corr_mean])
        hists['response_corr_mean_jes_up'].fill([bin.bin_center], weight=[bin.response_corr_mean_jes_up])
        hists['response_corr_mean_jes_down'].fill([bin.bin_center], weight=[bin.response_corr_mean_jes_down])
        hists['response_corr_mean_jer_up'].fill([bin.bin_center], weight=[bin.response_corr_mean_jer_up])
        hists['response_corr_mean_jer_down'].fill([bin.bin_center], weight=[bin.response_corr_mean_jer_down])
        hists['response_corr_std'].fill([bin.bin_center], weight=[bin.response_corr_std])
        hists['response_corr_std_jes_up'].fill([bin.bin_center], weight=[bin.response_corr_std_jes_up])
        hists['response_corr_std_jes_down'].fill([bin.bin_center], weight=[bin.response_corr_std_jes_down])
        hists['response_corr_std_jer_up'].fill([bin.bin_center], weight=[bin.response_corr_std_jer_up])
        hists['response_corr_std_jer_down'].fill([bin.bin_center], weight=[bin.response_corr_std_jer_down])

        # now estimate MC stat uncertainty on the mean by checking variance of mean values of 10 different sub samples:
        sub_means_raw = np.array([], dtype=np.float32)
        # sub_stds_raw = np.array([], dtype=np.float32)
        sub_means_corr = np.array([], dtype=np.float32)
        # sub_stds_corr = np.array([], dtype=np.float32)
        sub_weight_sums = np.array([], dtype=np.float32)
        split_mean_raw = np.array_split(bin.response_raw, 10)
        split_mean_corr = np.array_split(bin.response_corr, 10)
        split_weights = np.array_split(bin.weights, 10)
        # print(len(bin.weights))
        for i in range(len(split_mean_raw)):
            sub_means_raw = np.append(sub_means_raw, [ np.average(split_mean_raw[i], weights=split_weights[i]) ])
            sub_means_corr = np.append(sub_means_corr, [ np.average(split_mean_corr[i], weights=split_weights[i]) ])
            sub_weight_sums = np.append(sub_weight_sums, [ np.sum(split_weights[i]) ])
        avg, bin.response_raw_mean_mc_unc = weighted_avg_and_std(sub_means_raw, sub_weight_sums)
        avg, bin.response_corr_mean_mc_unc = weighted_avg_and_std(sub_means_corr, sub_weight_sums)
        hists['response_raw_mean_mc_unc'].fill([bin.bin_center], weight=[bin.response_raw_mean_mc_unc])
        hists['response_corr_mean_mc_unc'].fill([bin.bin_center], weight=[bin.response_corr_mean_mc_unc])

        sub_stds_raw, sub_stds_raw_weights = splitted_weighted_stds_with_common_average(bin.response_raw, bin.weights, 10)
        avg, bin.response_raw_std_mc_unc = weighted_avg_and_std(sub_stds_raw, sub_stds_raw_weights)
        hists['response_raw_std_mc_unc'].fill([bin.bin_center], weight=[bin.response_raw_std_mc_unc])

        sub_stds_corr, sub_stds_corr_weights = splitted_weighted_stds_with_common_average(bin.response_corr, bin.weights, 10)
        avg, bin.response_corr_std_mc_unc = weighted_avg_and_std(sub_stds_corr, sub_stds_corr_weights)
        hists['response_corr_std_mc_unc'].fill([bin.bin_center], weight=[bin.response_corr_std_mc_unc])



    # print('raw')
    # print(hists['response_raw'].view().value / hists['response_gen'].view().value)
    # print(np.sqrt(hists['response_raw'].view().variance / hists['response_gen'].view().value))
    #
    # print('corr')
    # print(hists['response_corr'].view().value / hists['response_gen'].view().value)
    # print(np.sqrt(hists['response_corr'].view().variance / hists['response_gen'].view().value))

    outputDir = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/output/WorkingPointStudy', year, 'nominal')
    outputFileName = 'HOTVR_JetPtResponse.root'
    outputFilePath = os.path.join(outputDir, outputFileName)
    with uproot.recreate(outputFilePath) as outfile:
        # outfile['response_gen'] = hists['response_gen']
        # outfile['response_raw'] = hists['response_raw']
        # outfile['response_corr'] = hists['response_corr']
        # outfile['response_corr_jes_up'] = hists['response_corr_jes_up']
        # outfile['response_corr_jes_down'] = hists['response_corr_jes_down']
        # outfile['response_corr_jer_up'] = hists['response_corr_jer_up']
        # outfile['response_corr_jer_down'] = hists['response_corr_jer_down']

        outfile['response_raw_mean'] = hists['response_raw_mean']
        outfile['response_raw_mean_mc_unc'] = hists['response_raw_mean_mc_unc']
        outfile['response_raw_std'] = hists['response_raw_std']
        outfile['response_raw_std_mc_unc'] = hists['response_raw_std_mc_unc']
        outfile['response_corr_mean'] = hists['response_corr_mean']
        outfile['response_corr_mean_mc_unc'] = hists['response_corr_mean_mc_unc']
        outfile['response_corr_mean_jes_up'] = hists['response_corr_mean_jes_up']
        outfile['response_corr_mean_jes_down'] = hists['response_corr_mean_jes_down']
        outfile['response_corr_mean_jer_up'] = hists['response_corr_mean_jer_up']
        outfile['response_corr_mean_jer_down'] = hists['response_corr_mean_jer_down']
        outfile['response_corr_std'] = hists['response_corr_std']
        outfile['response_corr_std_mc_unc'] = hists['response_corr_std_mc_unc']
        outfile['response_corr_std_jes_up'] = hists['response_corr_std_jes_up']
        outfile['response_corr_std_jes_down'] = hists['response_corr_std_jes_down']
        outfile['response_corr_std_jer_up'] = hists['response_corr_std_jer_up']
        outfile['response_corr_std_jer_down'] = hists['response_corr_std_jer_down']
