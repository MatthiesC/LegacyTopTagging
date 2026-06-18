from __future__ import print_function

import ROOT as root
import os
import json
from collections import OrderedDict
from copy import deepcopy
import numpy as np

import sys
sys.path.append(os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/Analysis'))
from constants import _YEARS

out_dir = os.path.join(os.environ.get("CMSSW_BASE"), "src/UHH2/LegacyTopTagging/Analysis/ControlPlots/plots")
data_dir = os.path.join(os.environ.get("CMSSW_BASE"), "src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel/data_noTree")
years = [
    'UL16preVFP',
    'UL16postVFP',
    'UL17',
    'UL18',
]
channels = [
    'ele',
    'muo',
]
data_files = OrderedDict()
data_hists = OrderedDict()
hist_name = "Presel_Main_Lumi/luminosity"
do_scale_last_bin_by_fraction = False
do_scale_last_bin_by_fraction = True

def scale_last_bin_by_fraction(hist, lumi):
    last_bin = hist.GetNbinsX()
    fraction = lumi - int(lumi)
    old_binc = hist.GetBinContent(last_bin)
    new_binc = old_binc / fraction
    old_bine = hist.GetBinError(last_bin)
    new_bine = old_bine / fraction
    hist.SetBinContent(last_bin, new_binc)
    hist.SetBinError(last_bin, new_bine)
    print("Old last binc/bine:", old_binc, old_bine)
    print("New last binc/bine:", new_binc, new_bine)

for channel in channels:
    for year in years:
        file_name = "uhh2.AnalysisModuleRunner.DATA.DATA.noTree_{}_{}.root".format(year, channel)
        file_path = os.path.join(data_dir, file_name)
        if not os.path.isfile(file_path):
            raise RuntimeError("File {} does not exist".format(file_path))
        data_files.setdefault(channel, {})[year] = file_path
        file = root.TFile.Open(file_path, "READ")
        hist = file.Get(hist_name)
        lumi = _YEARS[year]['lumi_fb']
        print("Bins:", hist.GetNbinsX(), "Lumi:", lumi)
        if do_scale_last_bin_by_fraction:
            scale_last_bin_by_fraction(hist, lumi)
        data_hists.setdefault(channel, {})[year] = deepcopy(hist)
        del hist
        file.Close()

print(json.dumps(data_files, indent=4))

# Important: Bins represent 1000 invpb (= 1 invfb)
# The last bin in each histogram should be scaled by the fraction of the incomplete invfb

for channel in channels:
    # Combine into one hist:
    n_bins = np.sum([data_hists[channel][year].GetNbinsX() for year in years])
    print("Total bins:", n_bins)
    hist = root.TH1F(channel, channel, n_bins, 0, n_bins)
    binc = []
    bine = []
    for year in years:
        data_hist = data_hists[channel][year]
        for i_bin in range(1, data_hist.GetNbinsX()+1):
            binc.append(data_hist.GetBinContent(i_bin))
            bine.append(data_hist.GetBinError(i_bin))
    # print(len(binc), binc)
    # print(len(bine), bine)
    assert n_bins == len(binc)
    assert n_bins == len(bine)
    for i_bin in range(1, n_bins+1):
        hist.SetBinContent(i_bin, binc[i_bin-1])
        hist.SetBinError(i_bin, bine[i_bin-1])

    if do_scale_last_bin_by_fraction:
        out_file_name = "lumi_hist_total_{}_TagAndProbe_scaled.root".format(channel)
    else:
        out_file_name = "lumi_hist_total_{}_TagAndProbe.root".format(channel)
    out_file_path = os.path.join(out_dir, out_file_name)
    out_file = root.TFile.Open(out_file_path, "RECREATE")
    out_file.cd()
    hist.Write()
    out_file.Close()

    ## Do the plotting

    canvas = root.TCanvas("canvas", "canvas", 2400, 600)
    canvas.cd()

    hist.Draw()

    if do_scale_last_bin_by_fraction:
        plot_name = "lumi_hist_total_{}_TagAndProbe_scaled.pdf".format(channel)
    else:
        plot_name = "lumi_hist_total_{}_TagAndProbe.pdf".format(channel)
    plot_path = os.path.join(out_dir, plot_name)
    canvas.SaveAs(plot_path)
