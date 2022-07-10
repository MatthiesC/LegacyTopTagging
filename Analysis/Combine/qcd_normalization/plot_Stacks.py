from __future__ import print_function

import os
import sys
import argparse

from collections import OrderedDict

import ROOT as root
import numpy as np

sys.path.append(os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/NicePlots/python'))
from plotter import NiceStackWithRatio, Process, human_format
sys.path.append(os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/Analysis'))
from constants import _YEARS, _TCOLORS, _BANDS

all_years = ['UL16preVFP', 'UL16postVFP', 'UL17', 'UL18']
all_channels = ['muo']

tagger_bases = ['ak8_t', 'hotvr_t', 'ak8_w']

bands = _BANDS.keys()

if not sys.argv[1:]: sys.exit('No arguments provided. Exit.')
parser = argparse.ArgumentParser()
parser.add_argument('-y', '--years', nargs='+', choices=all_years, default=all_years)
parser.add_argument('-c', '--channels', nargs='+', choices=all_channels, default=all_channels)
parser.add_argument('-t', '--tagger-base', default='ak8_t__tau', choices=tagger_bases)
# parser.add_argument('-l', '--better-labels', action='store_true') # Use custom y axis labels on main pad
args = parser.parse_args(sys.argv[1:])

tagger_base = args.tagger_base

for year in args.years:

    for channel in args.channels:

        for band in bands:

            workDir = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel', year, channel)
            inputDir = os.path.join(workDir, 'qcd_normalization_study')

            if band == 'Main':
                qcd_mc_tcolor = _TCOLORS.get('pyplot_orange')
            elif band == 'QCD':
                qcd_mc_tcolor = _TCOLORS.get('pyplot_blue')

            processes = [
                Process('MC.nonQCD', 'Other processes (MC)', root.kGray),
                Process('MC.QCD', 'QCD multijet (MC)', qcd_mc_tcolor),
            ]
            # processes = [
            #     Process('MC.nonQCD', 'Other processes (MC)', root.kGreen),
            #     Process('MC.QCD', 'QCD multijet (MC)', root.kRed),
            # ]

            inputFileName = 'qcd_normalization_histograms__'+tagger_base+'.root'
            inputPath = os.path.join(inputDir, inputFileName)

            nice = NiceStackWithRatio(
                infile_path = inputPath,
                infile_directory = _BANDS.get(band).name, # the directory within the ROOT file
                x_axis_title = '#it{M}_{T, W}',
                x_axis_unit = 'GeV',
                prepostfit = 'prefitRaw',
                processes = processes,
                lumi_unc = _YEARS.get(year).get('lumi_unc'),
                divide_by_bin_width = False,
                data_name = 'DATA.DATA',
                text_prelim = 'Private Work',
                text_top_left = _YEARS.get(year).get('long_name'),
                text_top_right = _YEARS.get(year).get('lumi_fb_display')+' fb^{#minus1} (13 TeV)',
                # nostack = True,
            )

            outFileName = 'plot_QCDestimate_PrefitStacks__'+year+'_'+channel+'_'+_BANDS.get(band).name+'__'+tagger_base+'.pdf'
            outDir = os.path.join(inputDir, 'plots')

            nice.plot()

            nice.canvas.cd()

            if tagger_base == 'ak8_t':
                pt_text_string = 'T&P selection w/ AK8 probe'
                pt_text_string2 = '#it{p}_{T, probe} > 300 GeV, |#eta_{probe}| < 2.5'
            elif tagger_base == 'ak8_w':
                pt_text_string = 'T&P selection w/ AK8 probe'
                pt_text_string2 = '#it{p}_{T, probe} > 200 GeV, |#eta_{probe}| < 2.5'
            elif tagger_base == 'hotvr_t':
                pt_text_string = 'T&P selection w/ HOTVR probe'
                pt_text_string2 = '#it{p}_{T, probe} > 200 GeV, |#eta_{probe}| < 2.5'

            tlatex_pt = root.TLatex(nice.coord.graph_to_pad_x(0.95), nice.coord.graph_to_pad_y(0.9), pt_text_string)
            tlatex_pt.SetTextAlign(31) # left top
            tlatex_pt.SetTextFont(42)
            tlatex_pt.SetTextSize(0.03)
            tlatex_pt.SetNDC()
            tlatex_pt.Draw()

            tlatex_pt2 = root.TLatex(nice.coord.graph_to_pad_x(0.95), nice.coord.graph_to_pad_y(0.815), pt_text_string2)
            tlatex_pt2.SetTextAlign(31) # left top
            tlatex_pt2.SetTextFont(42)
            tlatex_pt2.SetTextSize(0.03)
            tlatex_pt2.SetNDC()
            tlatex_pt2.Draw()

            if band == 'Main':
                pt_text_string3 = 'Signal region, isolated muon, #bf{pre-fit}'
            elif band == 'QCD':
                pt_text_string3 = 'Sideband region, non-iso. muon, #bf{pre-fit}'

            tlatex_pt3 = root.TLatex(nice.coord.graph_to_pad_x(0.95), nice.coord.graph_to_pad_y(0.73), pt_text_string3)
            tlatex_pt3.SetTextAlign(31) # left top
            tlatex_pt3.SetTextFont(42)
            tlatex_pt3.SetTextSize(0.03)
            tlatex_pt3.SetNDC()
            tlatex_pt3.Draw()

            legend = root.TLegend(nice.coord.graph_to_pad_x(0.54), nice.coord.graph_to_pad_y(0.4), nice.coord.graph_to_pad_x(0.95), nice.coord.graph_to_pad_y(0.69))

            legend.AddEntry(nice.data_hist, 'Data', 'ep')
            # data_label = 'Data, '+_YEARS.get(year).get('year')
            # legend.AddEntry(nice.data_hist, data_label, 'ep')

            # legend.AddEntry(nice.totalprocs, 'Total prediction', 'l')
            # legend.AddEntry(nice.stack[1], processes[1].legend, 'l')
            # legend.AddEntry(nice.stack[0], '#lower[0.00]{'+processes[0].legend+'}', 'l')
            legend.AddEntry(nice.stack.GetStack().At(1), processes[1].legend, 'f')
            legend.AddEntry(nice.stack.GetStack().At(0), '#lower[0.00]{'+processes[0].legend+'}', 'f')
            # legend.AddEntry(nice.stack_unc, 'Lumi. #oplus MC stat. unc.', 'f')
            legend.AddEntry(nice.stack_unc, 'Total uncertainty', 'f')

            legend.SetTextSize(0.03)
            legend.SetBorderSize(0)
            legend.SetFillStyle(0)
            legend.Draw()

            nice.save_plot(outFileName, outDir)



            #########################
            # Now do postfit plot
            if band == 'QCD': continue

            inputFileName = 'prepostfitshapes_'+tagger_base+'_rebinned.root'
            inputPath = os.path.join(inputDir, inputFileName)

            processes = [
                Process('nonQCD', 'Other processes (MC)', root.kGreen),
                Process('QCD_DD', 'QCD multijet (DD)', root.kRed),
            ]

            nice2 = NiceStackWithRatio(
                infile_path = inputPath,
                infile_directory = 'postfit', # the directory within the ROOT file
                x_axis_title = '#it{M}_{T, W}',
                x_axis_unit = 'GeV',
                prepostfit = 'postfitCombine',
                processes = processes,
                lumi_unc = _YEARS.get(year).get('lumi_unc'),
                divide_by_bin_width = False,
                text_prelim = 'Private Work',
                text_top_left = _YEARS.get(year).get('long_name'),
                text_top_right = _YEARS.get(year).get('lumi_fb_display')+' fb^{#minus1} (13 TeV)',
                nostack = True,
            )

            outFileName = outFileName.replace('PrefitStacks', 'Postfit')

            nice2.plot()

            nice2.canvas.cd()
            tlatex_pt.Draw()
            tlatex_pt2.Draw()
            if band == 'Main':
                pt_text_string3 = 'Signal region, isolated muon, #bf{post-fit}'
            elif band == 'QCD':
                pt_text_string3 = 'Sideband region, non-iso. muon, #bf{post-fit}'

            tlatex_pt3 = root.TLatex(nice2.coord.graph_to_pad_x(0.95), nice2.coord.graph_to_pad_y(0.73), pt_text_string3)
            tlatex_pt3.SetTextAlign(31) # left top
            tlatex_pt3.SetTextFont(42)
            tlatex_pt3.SetTextSize(0.03)
            tlatex_pt3.SetNDC()
            tlatex_pt3.Draw()


            legend2 = root.TLegend(nice2.coord.graph_to_pad_x(0.54), nice2.coord.graph_to_pad_y(0.4), nice2.coord.graph_to_pad_x(0.95), nice2.coord.graph_to_pad_y(0.69))
            # legend.SetHeader('#bf{Pre-fit}')

            legend2.AddEntry(nice2.data_hist, 'Data', 'ep')
            # data_label = 'Data, '+_YEARS.get(year).get('year')
            # legend.AddEntry(nice.data_hist, data_label, 'ep')

            # for i_process, process in reversed(list(enumerate(nice.processes))):
            #     legend.AddEntry(nice.stack.GetStack().At(i_process), process.legend, 'f')
            legend2.AddEntry(nice2.totalprocs, 'Total fit', 'l')
            legend2.AddEntry(nice2.stack[1], processes[1].legend, 'l')
            legend2.AddEntry(nice2.stack[0], '#lower[0.00]{'+processes[0].legend+'}', 'l')
            # legend.AddEntry(nice.stack_unc, 'Lumi. #oplus MC stat. unc.', 'f')
            legend2.AddEntry(nice2.stack_unc, 'Total uncertainty', 'f')
            legend2.SetTextSize(0.03)
            legend2.SetBorderSize(0)
            legend2.SetFillStyle(0)
            legend2.Draw()

            nice2.save_plot(outFileName, outDir)
