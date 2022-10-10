#!/usr/bin/env python2

from __future__ import print_function

import os
from subprocess import call

import ROOT as root
from copy import deepcopy
from collections import OrderedDict

import sys
sys.path.append(os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/Analysis'))
from constants import Systematics, _BANDS, _TAGGERS, _PT_INTERVALS_TANDP_AK8_T, _PT_INTERVALS_TANDP_AK8_W, _PT_INTERVALS_TANDP_HOTVR, _YEARS

systematics = Systematics(blacklist=['sfelec', 'sfmu_iso'])
systs = systematics.base

sys.path.append(os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/NicePlots/python'))
from plotter import NiceStackWithRatio, Process, human_format

from parallel_threading import run_with_pool



all_years = [
'UL16preVFP',
'UL16postVFP',
'UL17',
'UL18',
]

taggers = [
'ak8_t__tau',
'ak8_t_btagDJet__tau',
'hotvr_t__tau',
# 'ak8_w__partnet',
]
taggers = {k: _TAGGERS[k] for k in taggers}

# set working points:
for tagger_k, tagger in taggers.items():
    if tagger_k == 'ak8_t__tau' or tagger_k == 'ak8_t_btagDJet__tau' or tagger_k == 'ak8_t_btagDCSV__tau':
        # tagger.wps = [
        #     WorkingPoint(0.001, 0.38),
        #     WorkingPoint(0.005, 0.47),
        #     WorkingPoint(0.010, 0.52),
        #     WorkingPoint(0.025, 0.61),
        #     WorkingPoint(0.050, 0.69),
        # ]
        tagger.var_intervals = _PT_INTERVALS_TANDP_AK8_T
        tagger.fit_variable = 'output_probejet_AK8_mSD'
    if tagger_k == 'hotvr_t__tau':
        # hotvr_wp = WorkingPoint(9.999, 0.56) # providing bkg_eff does not make sense here; I have not derived it on my own; should be 3% or so
        # hotvr_wp.name = 'Standard'
        # tagger.wps = [hotvr_wp]
        tagger.var_intervals = _PT_INTERVALS_TANDP_HOTVR
        tagger.fit_variable = 'output_probejet_HOTVR_mass'
    if tagger_k == 'ak8_w__partnet':
        # tagger.wps = [
        #     WorkingPoint(0.030, 0.00), # cut value year dependent! Need to hack it in the code below where we define the cut rules...
        # ]
        tagger.var_intervals = _PT_INTERVALS_TANDP_AK8_W
        tagger.fit_variable = 'output_probejet_AK8_mSD'

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
'TTbar__MSc_Background',
'ST',
'ST__MSc_FullyMerged',
'ST__MSc_WMerged',
'ST__MSc_QBMerged',
'ST__MSc_NotMerged',
'ST__MSc_YllufMerged',
'ST__MSc_SemiMerged',
'ST__MSc_Background',
'VJetsAndVV',
]

processes_Plotter = [
Process('QCD', 'Multijet', root.kAzure+10),
Process('VJetsAndVV', 'V+jets, VV', root.kSpring-3),
# Process('ST', 'Single top', root.kYellow),
# Process('TTbar', 't#bar{t}', root.kRed),
Process('ST__MSc_Background', 'Background', root.kOrange+4),
Process('ST__MSc_YllufMerged', 'Not merged', root.kOrange-3),
Process('ST__MSc_FullyMerged', 'Fully merged', root.kOrange),
Process('TTbar__MSc_Background', 'Background', root.kPink-7),
Process('TTbar__MSc_YllufMerged', 'Not merged', root.kPink+4),
Process('TTbar__MSc_FullyMerged', 'Fully merged', root.kPink-3),
]
# processes_Plotter = OrderedDict([
#     (p.name, p) for p in processes_Plotter
# ])
processes_Plotter_temp = OrderedDict()
for index, p in enumerate(processes_Plotter):
    p.index = index
    processes_Plotter_temp[p.name] = p
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
]

bands = _BANDS


def hist_preprocessing(hist):
    for i in range(hist.GetNbinsX()+2):
        if hist.GetBinContent(i) < 0:
            hist.SetBinContent(i, 0)
            hist.SetBinError(i, 0)
    return hist


def create_rearranged_hists(variable, tagger, year, wp, pt_bin, do_plot=False, do_hists=True):
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

            for channel in channels:


                inFolderName = '/'.join([band.name, region, channel])
                outFolderName = '_'.join([band.name, region, channel])
                outFolder = outFile.mkdir(outFolderName)
                outFolderNames.append(outFolderName)
                print(inFolderName)
                # sys.exit()

                if not do_hists:
                    continue

                # data
                inFileName_nominal = '-'.join(['BasicHists', tagger.name, wp.name, pt_bin.name, year, 'nominal', variable])+'.root'
                inFilePath_nominal = os.path.join(inDir, inFileName_nominal)
                inFile_nominal = root.TFile.Open(inFilePath_nominal, 'READ')
                inHistPath_data = os.path.join(inFolderName, 'DATA__nominal')
                inHist_data = inFile_nominal.Get(inHistPath_data)
                inHist_data = hist_preprocessing(inHist_data)
                outFolder.cd()
                inHist_data.Write('data_obs')

                # MC
                for process in processes:

                    # nominal
                    inHistPath_nominal = os.path.join(inFolderName, process+'__nominal')
                    inHist_nominal = inFile_nominal.Get(inHistPath_nominal)
                    print(inHistPath_nominal)
                    inHist_nominal = hist_preprocessing(inHist_nominal)
                    outFolder.cd()
                    inHist_nominal.Write(process)

                    # systematics which do not need special treatment
                    for syst_k in sorted(systs.keys()):
                        syst = systs[syst_k]

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
                            inHist_variation = hist_preprocessing(inHist_variation)
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
                            else:
                                inHist_variation.Write(outHistName_variation)

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
                            inHist_variation = hist_preprocessing(inHist_variation)
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

                inFile_nominal.Close()

    outFile.Close()
    print('Wrote', outFilePath)

    if do_plot:

        for outFolderName in outFolderNames:

            if variable.endswith('_mSD'):
                x_axis_title = 'Probe jet #it{m}_{SD}'
                x_axis_unit = 'GeV'
            elif variable.endswith('_mass'):
                x_axis_title = 'Probe jet #it{m}_{jet}'
                x_axis_unit = 'GeV'

            nice = NiceStackWithRatio(
                infile_path = outFilePath,
                infile_directory = outFolderName, # the directory within the ROOT file
                x_axis_title = x_axis_title or 'undefined',
                x_axis_unit = x_axis_unit or 'undefined',
                prepostfit = 'prefitRaw',
                processes = processes_Plotter.values(),
                syst_names = systs_Plotter,
                lumi_unc = _YEARS.get(year).get('lumi_unc'),
                # divide_by_bin_width = False,
                data_name = 'data_obs',
                text_prelim = 'Private Work',
                text_top_left = _YEARS.get(year).get('long_name'),
                text_top_right = _YEARS.get(year).get('lumi_fb_display')+' fb^{#minus1} (13 TeV)',
                # nostack = True,
            )

            plotName = 'Plot-prefitRaw-'+task_name+'-'+outFolderName+'.pdf'
            plotDir = os.path.join(outDir, 'plots')

            nice.plot()

            #________________
            # Some text labels

            nice.canvas.cd()

            if tagger.name.startswith('ak8_'):
                probejetalgo = 'AK8 PUPPI'
            elif tagger.name.startswith('hotvr_'):
                probejetalgo = 'HOTVR PUPPI'

            # print(tagger.label)
            pt_text_string = tagger.label.replace('{WP_VALUE}', '{}'.format(wp.cut_value))
            for region in regions:
                if region in outFolderName.split('_'):
                    pt_text_string += ' [#bf{'+region+'}]'
            if pt_bin.max_set:
                pt_text_string2 = '#it{p}_{T} #in ({pt_min}, {pt_max}] GeV, |#eta| < 2.5'.replace('{pt_min}', '{}'.format(pt_bin.var_min)).replace('{pt_max}', '{}'.format(pt_bin.var_max))
            else:
                pt_text_string2 = '#it{p}_{T} > {pt_min} GeV, |#eta| < 2.5'.replace('{pt_min}', '{}'.format(pt_bin.var_min))

            tlatex_pt = root.TLatex(nice.coord.graph_to_pad_x(0.95), nice.coord.graph_to_pad_y(0.84), pt_text_string)
            tlatex_pt.SetTextAlign(31) # left top
            tlatex_pt.SetTextFont(42)
            tlatex_pt.SetTextSize(0.025)
            tlatex_pt.SetNDC()
            tlatex_pt.Draw()

            tlatex_pt2 = root.TLatex(nice.coord.graph_to_pad_x(0.95), nice.coord.graph_to_pad_y(0.78), pt_text_string2)
            tlatex_pt2.SetTextAlign(31) # left top
            tlatex_pt2.SetTextFont(42)
            tlatex_pt2.SetTextSize(0.025)
            tlatex_pt2.SetNDC()
            tlatex_pt2.Draw()

            tlatex_probejetalgo = root.TLatex(nice.coord.graph_to_pad_x(0.95), nice.coord.graph_to_pad_y(0.95), probejetalgo)
            tlatex_probejetalgo.SetTextAlign(33) # right top
            tlatex_probejetalgo.SetTextFont(62)
            tlatex_probejetalgo.SetTextSize(0.05)
            tlatex_probejetalgo.SetTextColor(root.kGray+1)
            tlatex_probejetalgo.SetNDC()
            tlatex_probejetalgo.Draw()

            #________________
            # Legend

            legend = root.TLegend(nice.coord.graph_to_pad_x(0.45), nice.coord.graph_to_pad_y(0.5), nice.coord.graph_to_pad_x(0.7), nice.coord.graph_to_pad_y(0.73))
            legend.AddEntry(nice.data_hist, 'Data', 'ep')
            legend.AddEntry(nice.stack.GetStack().At(processes_Plotter['VJetsAndVV'].index), processes_Plotter['VJetsAndVV'].legend, 'f')
            legend.AddEntry(nice.stack.GetStack().At(processes_Plotter['QCD'].index), processes_Plotter['QCD'].legend, 'f')
            legend.AddEntry(nice.stack_unc, 'Uncertainty', 'f')
            legend.SetTextSize(0.025)
            legend.SetBorderSize(0)
            legend.SetFillStyle(0)
            legend.Draw()

            legend2 = root.TLegend(nice.coord.graph_to_pad_x(0.45+0.22), nice.coord.graph_to_pad_y(0.5), nice.coord.graph_to_pad_x(0.7+0.22), nice.coord.graph_to_pad_y(0.73))
            legend2.SetHeader('t#bar{t} / single t categories:')
            legend2.AddEntry(nice.stack.GetStack().At(processes_Plotter['TTbar__MSc_FullyMerged'].index), '/', 'f')
            legend2.AddEntry(nice.stack.GetStack().At(processes_Plotter['TTbar__MSc_YllufMerged'].index), '/', 'f')
            legend2.AddEntry(nice.stack.GetStack().At(processes_Plotter['TTbar__MSc_Background'].index), '/', 'f')
            legend2.SetTextSize(0.025)
            legend2.SetBorderSize(0)
            legend2.SetFillStyle(0)
            legend2.Draw()

            legend3 = root.TLegend(nice.coord.graph_to_pad_x(0.45+0.22+0.067), nice.coord.graph_to_pad_y(0.5), nice.coord.graph_to_pad_x(0.7+0.22+0.067), nice.coord.graph_to_pad_y(0.73))
            legend3.SetHeader('')
            legend3.AddEntry(nice.stack.GetStack().At(processes_Plotter['ST__MSc_FullyMerged'].index), processes_Plotter['ST__MSc_FullyMerged'].legend, 'f')
            legend3.AddEntry(nice.stack.GetStack().At(processes_Plotter['ST__MSc_YllufMerged'].index), processes_Plotter['ST__MSc_YllufMerged'].legend, 'f')
            legend3.AddEntry(nice.stack.GetStack().At(processes_Plotter['ST__MSc_Background'].index), processes_Plotter['ST__MSc_Background'].legend, 'f')
            legend3.SetTextSize(0.025)
            legend3.SetBorderSize(0)
            legend3.SetFillStyle(0)
            legend3.Draw()






            nice.save_plot(plotName, plotDir)


if __name__=='__main__':

    # # the_tagger = _TAGGERS['ak8_t_btagDJet__tau']
    # # the_tagger = _TAGGERS['hotvr_t__tau']
    # the_tagger = _TAGGERS['ak8_t__tau']
    # # year = 'UL16preVFP'
    # year = 'UL18'
    # the_tagger.var_intervals = _PT_INTERVALS_TANDP_AK8_T
    # # the_tagger.var_intervals = _PT_INTERVALS_TANDP_HOTVR
    # for bin in the_tagger.var_intervals.values():
    #     pt_bin = bin
    #     create_rearranged_hists(the_tagger.fit_variable, the_tagger, year, the_tagger.get_wp(4, year), pt_bin, do_plot=True)





    for the_tagger in taggers.values():
        print('Working on', the_tagger.name)
        for year in all_years:
            print('Working on', year)
            for wp in the_tagger.get_wp(year=year):
                print('Working on', wp.name)
                for pt_bin in the_tagger.var_intervals.values():
                    print('Working on', pt_bin.name)
                    create_rearranged_hists(the_tagger.fit_variable, the_tagger, year, wp, pt_bin, do_plot=True)

    # sets_of_args = []
    # for the_tagger in taggers.values():
    #     for year in all_years:
    #         for wp in the_tagger.get_wp(year=year):
    #             for pt_bin in the_tagger.var_intervals.values():
    #                 set_of_args = {}
    #                 set_of_args['tagger'] = the_tagger
    #                 set_of_args['year'] = year
    #                 set_of_args['variable'] = the_tagger.fit_variable
    #                 set_of_args['wp'] = wp
    #                 set_of_args['pt_bin'] = pt_bin
    #                 set_of_args['do_plot'] = True
    #                 sets_of_args.append(set_of_args)
    #
    # run_with_pool(function=create_rearranged_hists, sets_of_args=sets_of_args)
