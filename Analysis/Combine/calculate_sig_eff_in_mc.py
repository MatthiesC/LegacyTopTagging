#!/usr/bin/env python2

from __future__ import print_function

import ROOT as root
import os
from collections import OrderedDict
import numpy as np

import sys
sys.path.append(os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/Analysis'))
from constants import Systematics, _BANDS, _TAGGERS, _PT_INTERVALS_TANDP_AK8_T, _PT_INTERVALS_TANDP_AK8_W, _PT_INTERVALS_TANDP_HOTVR, _YEARS, get_variable_binning_xlabel_xunit, _NULL_WP


all_years = [
'UL16preVFP',
'UL16postVFP',
'UL17',
'UL18',
]

taggers = [
'ak8_t__tau',
'ak8_t_btagDJet__tau',
# 'ak8_t_btagDCSV__tau',
'hotvr_t__tau',
# 'ak8_w__partnet',
# 'ak8_t__MDdeepak8',
]
taggers = {k: _TAGGERS[k] for k in taggers}



all_channels = [
'ele',
'muo',
]

all_regions = [
'Pass',
'Fail',
]

data_path = 'sig_eff_in_mc_BasicHists'


for the_tagger in taggers.values():

    print('Working on', the_tagger.name)

    is_hotvr = the_tagger.name.startswith('hotvr_t')
    the_var = 'mass' if is_hotvr else 'mSD'
    the_var_tagMin = 140. if is_hotvr else 105.
    the_var_tagMax = 220. if is_hotvr else 210.
    the_jetalgo = 'HOTVR' if is_hotvr else 'AK8'
    hist_name = 'TTbar__MSc_FullyMerged__nominal'
    if the_tagger.name.startswith('ak8_w'):
        the_var_tagMin = 65.
        the_var_tagMax = 105.
        hist_name = 'TTbar__MSc_WMerged__nominal'

    sorted_pt_bins = []
    for pt_bin in the_tagger.var_intervals.values():
        if pt_bin.fit == True:
            sorted_pt_bins.append(pt_bin)
    sorted_pt_bins = sorted(sorted_pt_bins, key=lambda x: x.var_min, reverse=False)
    pt_bins = OrderedDict()
    for pt_bin in sorted_pt_bins:
        pt_bins[pt_bin.name] = pt_bin

    # Create output ROOT file:
    outfileName = 'sig_eff_in_mc-'+the_tagger.name+'.root'
    outfilePath = os.path.join(data_path, outfileName)
    outfile = root.TFile.Open(outfilePath, 'RECREATE')
    
    for year in all_years:

        print('Working on', year)

        for wp in the_tagger.get_wp(year=year):

            print('Working on', wp.name)

            # Create TGraphs
            n_points = len(pt_bins)
            graph_name_base = '-'.join(['sig_eff', the_tagger.name, year, wp.name])

            null_array = np.zeros(n_points)

            graph_sig_eff_withWindow_muo = root.TGraph(n_points, null_array, null_array)
            graph_sig_eff_withWindow_ele = root.TGraph(n_points, null_array, null_array)
            graph_sig_eff_withWindow_both = root.TGraph(n_points, null_array, null_array)
            graph_sig_eff_noWindow_muo = root.TGraph(n_points, null_array, null_array)
            graph_sig_eff_noWindow_ele = root.TGraph(n_points, null_array, null_array)
            graph_sig_eff_noWindow_both = root.TGraph(n_points, null_array, null_array)

            graph_sig_eff_withWindow_muo.SetName('graph_sig_eff_withWindow_muo')
            graph_sig_eff_withWindow_ele.SetName('graph_sig_eff_withWindow_ele')
            graph_sig_eff_withWindow_both.SetName('graph_sig_eff_withWindow_both')
            graph_sig_eff_noWindow_muo.SetName('graph_sig_eff_noWindow_muo')
            graph_sig_eff_noWindow_ele.SetName('graph_sig_eff_noWindow_ele')
            graph_sig_eff_noWindow_both.SetName('graph_sig_eff_noWindow_both')

            pt_bin_i = -1

            for pt_bin in pt_bins.values():

                pt_bin_i += 1

                print('Working on', pt_bin.name)

                pt_bin_center = (pt_bin.var_max - pt_bin.var_min) * 0.5 + pt_bin.var_min
                if pt_bin_center > 600: # FIXME: leonidas' uppermost pt bins have fixed upper values; i'll just use randomly 800 GeV
                    pt_bin_center = 700
                
                # infileName = 'BasicHists-hotvr_t__tau-Standard-pt_200to250-UL16postVFP-nominal-output_probejet_HOTVR_mass.root'
                infileName = '-'.join(['BasicHists', the_tagger.name, wp.name, pt_bin.name, year, 'nominal', 'output_probejet_{}_{}'.format(the_jetalgo, the_var)])+'.root'
                infilePath = os.path.join(data_path, the_tagger.name, infileName)
                infile = root.TFile.Open(infilePath, 'READ')
                # infile.Close()

                hist_pass_muo = infile.Get('Main/Pass/muo/'+hist_name)
                hist_pass_ele = infile.Get('Main/Pass/ele/'+hist_name)
                hist_fail_muo = infile.Get('Main/Fail/muo/'+hist_name)
                hist_fail_ele = infile.Get('Main/Fail/ele/'+hist_name)

                nBins_pass_muo = hist_pass_muo.GetNbinsX()
                nBins_pass_ele = hist_pass_ele.GetNbinsX()
                nBins_fail_muo = hist_fail_muo.GetNbinsX()
                nBins_fail_ele = hist_fail_ele.GetNbinsX()

                n_pass_muo_inWindow = 0.
                n_pass_ele_inWindow = 0.
                n_pass_muo_outside = 0.
                n_pass_ele_outside = 0.
                n_fail_muo = 0.
                n_fail_ele = 0.

                for i_bin in range(nBins_pass_muo + 2):
                    bin_content = hist_pass_muo.GetBinContent(i_bin)
                    if i_bin == 0 or i_bin == nBins_pass_muo + 1:
                        n_pass_muo_outside += bin_content
                    else:
                        bin_center = hist_pass_muo.GetBinCenter(i_bin)
                        if bin_center > the_var_tagMin and bin_center < the_var_tagMax:
                            n_pass_muo_inWindow += bin_content
                        else:
                            n_pass_muo_outside += bin_content

                for i_bin in range(nBins_pass_ele + 2):
                    bin_content = hist_pass_ele.GetBinContent(i_bin)
                    if i_bin == 0 or i_bin == nBins_pass_ele + 1:
                        n_pass_ele_outside += bin_content
                    else:
                        bin_center = hist_pass_ele.GetBinCenter(i_bin)
                        if bin_center > the_var_tagMin and bin_center < the_var_tagMax:
                            n_pass_ele_inWindow += bin_content
                        else:
                            n_pass_ele_outside += bin_content

                for i_bin in range(nBins_fail_muo + 2):
                    n_fail_muo += hist_fail_muo.GetBinContent(i_bin)

                for i_bin in range(nBins_fail_ele + 2):
                    n_fail_ele += hist_fail_ele.GetBinContent(i_bin)

                n_pass_muo = n_pass_muo_inWindow + n_pass_muo_outside
                n_pass_ele = n_pass_ele_inWindow + n_pass_ele_outside

                n_tot_muo = n_pass_muo + n_fail_muo
                n_tot_ele = n_pass_ele + n_fail_ele

                sig_eff_withWindow_muo = n_pass_muo_inWindow / n_tot_muo
                sig_eff_withWindow_ele = n_pass_ele_inWindow / n_tot_ele
                sig_eff_withWindow_both = (n_pass_muo_inWindow + n_pass_ele_inWindow) / (n_tot_muo + n_tot_ele)

                sig_eff_noWindow_muo = n_pass_muo / n_tot_muo
                sig_eff_noWindow_ele = n_pass_ele / n_tot_ele
                sig_eff_noWindow_both = (n_pass_muo + n_pass_ele) / (n_tot_muo + n_tot_ele)

                print('sig_eff_withWindow_muo', sig_eff_withWindow_muo)
                print('sig_eff_withWindow_ele', sig_eff_withWindow_ele)
                print('sig_eff_withWindow_both', sig_eff_withWindow_both)
                print('sig_eff_noWindow_muo', sig_eff_noWindow_muo)
                print('sig_eff_noWindow_ele', sig_eff_noWindow_ele)
                print('sig_eff_noWindow_both', sig_eff_noWindow_both)

                graph_sig_eff_withWindow_muo.SetPoint(pt_bin_i, pt_bin_center, sig_eff_withWindow_muo)
                graph_sig_eff_withWindow_ele.SetPoint(pt_bin_i, pt_bin_center, sig_eff_withWindow_ele)
                graph_sig_eff_withWindow_both.SetPoint(pt_bin_i, pt_bin_center, sig_eff_withWindow_both)
                graph_sig_eff_noWindow_muo.SetPoint(pt_bin_i, pt_bin_center, sig_eff_noWindow_muo)
                graph_sig_eff_noWindow_ele.SetPoint(pt_bin_i, pt_bin_center, sig_eff_noWindow_ele)
                graph_sig_eff_noWindow_both.SetPoint(pt_bin_i, pt_bin_center, sig_eff_noWindow_both)

            # Write graphs to outfile:
            dirName = year+'-'+wp.name
            rootDir = outfile.mkdir(dirName)
            rootDir.cd()
            graph_sig_eff_withWindow_muo.Write()
            graph_sig_eff_withWindow_ele.Write()
            graph_sig_eff_withWindow_both.Write()
            graph_sig_eff_noWindow_muo.Write()
            graph_sig_eff_noWindow_ele.Write()
            graph_sig_eff_noWindow_both.Write()

    outfile.Close()
