#!/usr/bin/env python2

from __future__ import print_function

import os

import ROOT as root
from copy import deepcopy
from collections import OrderedDict
# import numpy as np

import sys
sys.path.append(os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/Analysis'))
from constants import _TAGGERS, _PT_INTERVALS_TANDP_AK8_T, _PT_INTERVALS_TANDP_AK8_W, _PT_INTERVALS_TANDP_HOTVR, _YEARS, _NULL_WP

sys.path.append(os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/NicePlots/python'))
from plotter import CoordinateConverter


all_years = [
'UL16preVFP',
'UL16postVFP',
'UL17',
'UL18',
]

all_channels = [
'muo',
'ele',
]

taggers = [
'ak8_t__tau',
'ak8_t_btagDJet__tau',
# 'ak8_t_btagDCSV__tau',
'hotvr_t__tau',
# 'ak8_w__partnet',
# 'ak8_t__MDdeepak8',
]
taggers = OrderedDict(
    [(k, _TAGGERS[k]) for k in taggers]
)

base_path = os.path.join(
    os.environ.get("CMSSW_BASE"),
    'src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel',
)

#/data/dust/user/matthies/uhh2-106X-v2/CMSSW_10_6_28/src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel/UL18/combine/ak8_t__tau/BasicHists-MSD/
#BasicHists-ak8_t__tau-NullWP-pt_600toInf-UL18-nominal-output_probejet_AK8_mSD.root

def add_histos(tagger, wp, pt_bin, process='TTbar__MSc_FullyMerged__nominal', max_value=300):

    histo = None
    variable_name = 'output_probejet_' + ('AK8_mSD' if tagger.name.startswith('ak8') else 'HOTVR_mass')

    for year in all_years:
        infile_name = '-'.join(['BasicHists', tagger.name, wp.name, pt_bin.name, year, 'nominal', variable_name]) + '.root'
        infile_path = os.path.join(
            base_path,
            year,
            'combine',
            tagger.name,
            'BasicHists-MSD',
            infile_name,
        )
        infile = root.TFile.Open(infile_path, 'READ')
        for channel in all_channels:
            inhist_name = 'Main/Pass/' + channel + '/' + process
            inhist = infile.Get(inhist_name)
            if histo is None:
                histo = deepcopy(inhist)
            else:
                histo.Add(inhist)
            inhist.Delete()
        infile.Close()

    #histo.Rebin(2)
    if max_value is not None:
        n_bins = histo.GetNbinsX()
        sum_above_max_value = histo.GetBinContent(n_bins + 1)  # overflow bin
        last_bin_before_max = None
        for i_bin in range(0, n_bins + 1):
            if histo.GetBinCenter(i_bin) > max_value:
                sum_above_max_value += histo.GetBinContent(i_bin)
                histo.SetBinContent(i_bin, 0)
            else:
                last_bin_before_max = i_bin
        histo.SetBinContent(last_bin_before_max, sum_above_max_value)
        # histo.SetBinContent(n_bins + 1, 0)

    return histo


def get_maximum(list_of_histos):
    maximum = 0.
    for histo in list_of_histos:
        for i_bin in range(1, histo.GetNbinsX() + 1):
            maximum = max(maximum, histo.GetBinContent(i_bin))
    return maximum


def normalize_histo(histo):
    integral = 0.
    n_bins = histo.GetNbinsX()
    for i_bin in range(0, n_bins + 2):
        integral += histo.GetBinContent(i_bin)
    if integral <= 0:
        raise ValueError("integral <= 0")
    for i_bin in range(0, n_bins + 2):
        histo.SetBinContent(i_bin, histo.GetBinContent(i_bin) / integral)
        histo.SetBinError(i_bin, histo.GetBinError(i_bin) / integral)


def create_plot_ak8(wp=None, normalized=False):

    tagger_names = [
        'ak8_t__tau',
        'ak8_t_btagDJet__tau',
    ]
    pt_bins = _PT_INTERVALS_TANDP_AK8_T

    histos = OrderedDict()
    for tagger_name in tagger_names:
        the_tagger = taggers[tagger_name]
        for pt_bin in pt_bins.values():
            histos.setdefault(tagger_name, OrderedDict())[pt_bin.name] = add_histos(
                tagger=the_tagger,
                wp=wp,
                pt_bin=pt_bin,
            )

    # print(
    #     histos['ak8_t__tau']['pt_300toInf'].GetBinContent(20)
    # )


    canvas = root.TCanvas('canvas', '', 600, 600)
    canvas.cd()

    root.gStyle.SetOptStat(0)

    max_value = 300

    text_size = 0.035

    margin_l = 0.15
    margin_r = 0.05
    margin_b = 0.12
    margin_t = 0.08
    tick_length = 0.015 # fraction of canvas width/height
    canvas.SetMargin(margin_l, margin_r, margin_b, margin_t) #lrbt
    coord = CoordinateConverter(margin_l, margin_r, margin_b, margin_t)
    canvas.SetTickx(1)
    canvas.SetTicky(1)


    histo_pt_300to400_nom = deepcopy(histos['ak8_t__tau']['pt_300to400'])
    histo_pt_300to400_nom.SetLineWidth(2)
    histo_pt_300to400_nom.SetLineColor(root.kRed)
    histo_pt_300to400_nom.Draw('hist')

    histo_pt_300to400_sjb = deepcopy(histos['ak8_t_btagDJet__tau']['pt_300to400'])
    histo_pt_300to400_sjb.SetLineWidth(2)
    histo_pt_300to400_sjb.SetLineColor(root.kRed)
    histo_pt_300to400_sjb.SetLineStyle(2)
    histo_pt_300to400_sjb.Draw('hist same')


    histo_pt_400to480_nom = deepcopy(histos['ak8_t__tau']['pt_400to480'])
    histo_pt_400to480_nom.SetLineWidth(2)
    histo_pt_400to480_nom.SetLineColor(root.kMagenta)
    histo_pt_400to480_nom.Draw('hist same')

    histo_pt_400to480_sjb = deepcopy(histos['ak8_t_btagDJet__tau']['pt_400to480'])
    histo_pt_400to480_sjb.SetLineWidth(2)
    histo_pt_400to480_sjb.SetLineColor(root.kMagenta)
    histo_pt_400to480_sjb.SetLineStyle(2)
    histo_pt_400to480_sjb.Draw('hist same')


    histo_pt_480to600_nom = deepcopy(histos['ak8_t__tau']['pt_480to600'])
    histo_pt_480to600_nom.SetLineWidth(2)
    histo_pt_480to600_nom.SetLineColor(root.kBlue)
    histo_pt_480to600_nom.Draw('hist same')

    histo_pt_480to600_sjb = deepcopy(histos['ak8_t_btagDJet__tau']['pt_480to600'])
    histo_pt_480to600_sjb.SetLineWidth(2)
    histo_pt_480to600_sjb.SetLineColor(root.kBlue)
    histo_pt_480to600_sjb.SetLineStyle(2)
    histo_pt_480to600_sjb.Draw('hist same')


    histo_pt_600toInf_nom = deepcopy(histos['ak8_t__tau']['pt_600toInf'])
    histo_pt_600toInf_nom.SetLineWidth(2)
    histo_pt_600toInf_nom.SetLineColor(root.kCyan)
    histo_pt_600toInf_nom.Draw('hist same')

    histo_pt_600toInf_sjb = deepcopy(histos['ak8_t_btagDJet__tau']['pt_600toInf'])
    histo_pt_600toInf_sjb.SetLineWidth(2)
    histo_pt_600toInf_sjb.SetLineColor(root.kCyan)
    histo_pt_600toInf_sjb.SetLineStyle(2)
    histo_pt_600toInf_sjb.Draw('hist same')

    all_histos = [
        histo_pt_300to400_nom,
        histo_pt_300to400_sjb,
        histo_pt_400to480_nom,
        histo_pt_400to480_sjb,
        histo_pt_480to600_nom,
        histo_pt_480to600_sjb,
        histo_pt_600toInf_nom,
        histo_pt_600toInf_sjb,
    ]

    if normalized:
        for histo in all_histos:
            normalize_histo(histo)

    maximum = get_maximum(all_histos)
    histo_pt_300to400_nom.SetMaximum(maximum * 1.15)
    histo_pt_300to400_nom.SetMinimum(0)

    if normalized:
        histo_pt_300to400_nom.GetYaxis().SetTitle('Fraction of events / 5 GeV')
    else:
        histo_pt_300to400_nom.GetYaxis().SetTitle('Events / 5 GeV')
    histo_pt_300to400_nom.GetYaxis().SetTitleSize(0.045)
    histo_pt_300to400_nom.GetYaxis().SetTitleOffset(1.4)
    # histo_pt_300to400_nom.GetYaxis().CenterTitle()

    histo_pt_300to400_nom.GetXaxis().SetRangeUser(0, max_value)
    # histo_pt_300to400_nom.GetXaxis().SetLimits(0, max_value)

    histo_pt_300to400_nom.GetXaxis().SetTitle('Probe jet #it{m}_{SD} [GeV]')
    histo_pt_300to400_nom.GetXaxis().SetTitleSize(0.045)
    histo_pt_300to400_nom.GetXaxis().SetTitleOffset(1.2)
    # histo_pt_300to400_nom.GetXaxis().CenterTitle()

    histo_pt_300to400_nom.GetXaxis().ChangeLabel(7,-1,-1,-1,-1,-1,'#geq 300');


    leg_x1 = 0.44 - 0.25
    leg_y1 = 0.25
    leg_x2 = 0.77 - 0.25
    leg_y2 = 0.5

    legend = root.TLegend(leg_x1, leg_y1, leg_x2, leg_y2)

    legend.AddEntry(histo_pt_300to400_nom, "300 < #it{p}_{T}/GeV < 400", "l")
    legend.AddEntry(histo_pt_400to480_nom, "400 < #it{p}_{T}/GeV < 480", "l")
    legend.AddEntry(histo_pt_480to600_nom, "480 < #it{p}_{T}/GeV < 600", "l")
    legend.AddEntry(histo_pt_600toInf_nom, "#it{p}_{T} > 600 GeV", "l")

    legend.SetTextSize(0.03)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    legend.Draw()


    # header_leg2 = '#bf{Fully merged top quark, t#bar{t} MC}'
    header_leg2 = '#bf{t#bar{t} MC, fully merged t}'
    leg2_x1 = 0.44 - 0.25
    leg2_y1 = 0.72 - 0.08
    leg2_x2 = 0.77 - 0.25
    leg2_y2 = 0.87 - 0.08
    legend2 = root.TLegend(leg2_x1, leg2_y1, leg2_x2, leg2_y2)

    histo_pt_300to400_nom_leg2 = deepcopy(histo_pt_300to400_nom)
    histo_pt_300to400_nom_leg2.SetLineColor(root.kBlack)
    histo_pt_300to400_sjb_leg2 = deepcopy(histo_pt_300to400_sjb)
    histo_pt_300to400_sjb_leg2.SetLineColor(root.kBlack)

    legend2.SetHeader(header_leg2)
    legend2.AddEntry(histo_pt_300to400_nom_leg2, "w/o subjet b tag")
    # legend2.AddEntry(histo_pt_300to400_sjb_leg2, "w/ loose DeepJet subjet b tag")
    legend2.AddEntry(histo_pt_300to400_sjb_leg2, "w/ subjet b tag")

    legend2.SetTextSize(0.03)
    legend2.SetBorderSize(0)
    legend2.SetFillStyle(0)

    legend2.Draw()


    text_pw1 = root.TLatex(coord.graph_to_pad_x(0.05), coord.graph_to_pad_y(0.95), 'Private work')
    text_pw1.SetTextAlign(13)
    text_pw1.SetTextFont(52)
    text_pw1.SetTextSize(0.03)
    text_pw1.SetNDC()
    text_pw1.Draw()

    text_pw2 = root.TLatex(coord.graph_to_pad_x(0.05), coord.graph_to_pad_y(0.90), '(CMS data/simulation)')
    text_pw2.SetTextAlign(13)
    text_pw2.SetTextFont(52)
    text_pw2.SetTextSize(0.03)
    text_pw2.SetNDC()
    text_pw2.Draw()


    tlatex_probejetalgo = root.TLatex(coord.graph_to_pad_x(0.95), coord.graph_to_pad_y(0.95), 'AK8 PUPPI')
    tlatex_probejetalgo.SetTextAlign(33) # right top
    tlatex_probejetalgo.SetTextFont(62)
    tlatex_probejetalgo.SetTextSize(0.05)
    tlatex_probejetalgo.SetTextColor(root.kGray+1)
    tlatex_probejetalgo.SetNDC()
    tlatex_probejetalgo.Draw()


    y_shift = 0.031

    tlatex_pt1 = '#tau_{3}/#tau_{2} < ' + str(wp.get_cut_value())# + '(' + wp.alt_name + ')'
    # tlatex_pt1 = 'HOTVR t tag incl. #tau_{3}/#tau_{2} < 0.56'
    tlatex_pt1 = root.TLatex(coord.graph_to_pad_x(0.95), coord.graph_to_pad_y(0.867 - y_shift), tlatex_pt1)
    tlatex_pt1.SetTextAlign(31) # left top
    tlatex_pt1.SetTextFont(42)
    tlatex_pt1.SetTextSize(0.035)
    tlatex_pt1.SetNDC()
    tlatex_pt1.Draw()

    tlatex_pt2 = '(' + wp.alt_name + ')'
    tlatex_pt2 = root.TLatex(coord.graph_to_pad_x(0.95), coord.graph_to_pad_y(0.801 - y_shift), tlatex_pt2)
    tlatex_pt2.SetTextAlign(31) # left top
    tlatex_pt2.SetTextFont(42)
    tlatex_pt2.SetTextSize(0.035)
    tlatex_pt2.SetNDC()
    tlatex_pt2.Draw()

    tlatex_pt3 = '|#eta| < 2.5'
    tlatex_pt3 = root.TLatex(coord.graph_to_pad_x(0.95), coord.graph_to_pad_y(0.735 - y_shift), tlatex_pt3)
    tlatex_pt3.SetTextAlign(31) # left top
    tlatex_pt3.SetTextFont(42)
    tlatex_pt3.SetTextSize(0.035)
    tlatex_pt3.SetNDC()
    tlatex_pt3.Draw()


    text_top_left = 'T&P e/#mu+jets'
    # text_top_right = '#sqrt{#it{s}} = 13 TeV'
    text_top_right = str(_YEARS['run2']['lumi_fb_display']) + ' fb^{#minus1} (13 TeV)'

    if text_top_left:
        tlatex_top_left = root.TLatex(margin_l, 1. - margin_t + 0.01, text_top_left)
        tlatex_top_left.SetTextAlign(11)  # left bottom
        tlatex_top_left.SetTextFont(42)
        tlatex_top_left.SetTextSize(text_size)
        tlatex_top_left.SetNDC()
        tlatex_top_left.Draw()

    if text_top_right:
        tlatex_top_right = root.TLatex(1. - margin_r, 1. - margin_t + 0.01, text_top_right)
        tlatex_top_right.SetTextAlign(31)  # right bottom
        tlatex_top_right.SetTextFont(42)
        tlatex_top_right.SetTextSize(text_size)
        tlatex_top_right.SetNDC()
        tlatex_top_right.Draw()


    canvas.Update()
    tickScaleX = (canvas.GetUxmax() - canvas.GetUxmin()) / (canvas.GetX2() - canvas.GetX1()) * (canvas.GetWh() * canvas.GetAbsHNDC())
    tickScaleY = (canvas.GetUymax() - canvas.GetUymin()) / (canvas.GetY2() - canvas.GetY1()) * (canvas.GetWw() * canvas.GetAbsWNDC())
    histo_pt_300to400_nom.GetXaxis().SetTickLength(canvas.GetWh() * tick_length / tickScaleX)
    histo_pt_300to400_nom.GetYaxis().SetTickLength(canvas.GetWw() * tick_length / tickScaleY)

    root.gPad.Update()
    root.gPad.RedrawAxis()

    output_path = 'plots_pass_region'
    os.system('mkdir -p '+output_path)

    output_file_name = 'plotPassRegion_ak8_' + wp.name
    if normalized:
        output_file_name += '_norm'
    output_file_name += '.pdf'
    output_file_path = os.path.join(output_path, output_file_name)

    canvas.SaveAs(output_file_path)



def create_plot_hotvr(wp=None, normalized=False):

    tagger_names = [
        'hotvr_t__tau',
    ]
    pt_bins = _PT_INTERVALS_TANDP_HOTVR

    histos = OrderedDict()
    for tagger_name in tagger_names:
        the_tagger = taggers[tagger_name]
        for pt_bin in pt_bins.values():
            histos.setdefault(tagger_name, OrderedDict())[pt_bin.name] = add_histos(
                tagger=the_tagger,
                wp=wp,
                pt_bin=pt_bin,
            )

    # print(
    #     histos['ak8_t__tau']['pt_300toInf'].GetBinContent(20)
    # )


    canvas = root.TCanvas('canvas', '', 600, 600)
    canvas.cd()

    root.gStyle.SetOptStat(0)

    max_value = 300

    text_size = 0.035

    margin_l = 0.15
    margin_r = 0.05
    margin_b = 0.12
    margin_t = 0.08
    tick_length = 0.015 # fraction of canvas width/height
    canvas.SetMargin(margin_l, margin_r, margin_b, margin_t) #lrbt
    coord = CoordinateConverter(margin_l, margin_r, margin_b, margin_t)
    canvas.SetTickx(1)
    canvas.SetTicky(1)


    linestyles = True
    linestyles = False

    histo_pt_200to250_nom = deepcopy(histos['hotvr_t__tau']['pt_200to250'])
    histo_pt_200to250_nom.SetLineWidth(2)
    if linestyles: histo_pt_200to250_nom.SetLineStyle(1)
    histo_pt_200to250_nom.SetLineColor(root.kBlack)
    histo_pt_200to250_nom.Draw('hist')

    histo_pt_250to300_nom = deepcopy(histos['hotvr_t__tau']['pt_250to300'])
    histo_pt_250to300_nom.SetLineWidth(2)
    if linestyles: histo_pt_250to300_nom.SetLineStyle(2)
    histo_pt_250to300_nom.SetLineColor(root.kGreen)
    histo_pt_250to300_nom.Draw('hist same')

    histo_pt_300to350_nom = deepcopy(histos['hotvr_t__tau']['pt_300to350'])
    histo_pt_300to350_nom.SetLineWidth(2)
    if linestyles: histo_pt_300to350_nom.SetLineStyle(3)
    histo_pt_300to350_nom.SetLineColor(root.kOrange)
    histo_pt_300to350_nom.Draw('hist same')

    histo_pt_350to400_nom = deepcopy(histos['hotvr_t__tau']['pt_350to400'])
    histo_pt_350to400_nom.SetLineWidth(2)
    if linestyles: histo_pt_350to400_nom.SetLineStyle(4)
    histo_pt_350to400_nom.SetLineColor(root.kRed)
    histo_pt_350to400_nom.Draw('hist same')

    histo_pt_400to480_nom = deepcopy(histos['hotvr_t__tau']['pt_400to480'])
    histo_pt_400to480_nom.SetLineWidth(2)
    if linestyles: histo_pt_400to480_nom.SetLineStyle(5)
    histo_pt_400to480_nom.SetLineColor(root.kMagenta)
    histo_pt_400to480_nom.Draw('hist same')

    histo_pt_480to600_nom = deepcopy(histos['hotvr_t__tau']['pt_480to600'])
    histo_pt_480to600_nom.SetLineWidth(2)
    if linestyles: histo_pt_480to600_nom.SetLineStyle(6)
    histo_pt_480to600_nom.SetLineColor(root.kBlue)
    histo_pt_480to600_nom.Draw('hist same')

    histo_pt_600toInf_nom = deepcopy(histos['hotvr_t__tau']['pt_600toInf'])
    histo_pt_600toInf_nom.SetLineWidth(2)
    if linestyles: histo_pt_600toInf_nom.SetLineStyle(7)
    histo_pt_600toInf_nom.SetLineColor(root.kCyan)
    histo_pt_600toInf_nom.Draw('hist same')

    all_histos = [
        histo_pt_200to250_nom,
        histo_pt_250to300_nom,
        histo_pt_300to350_nom,
        histo_pt_350to400_nom,
        histo_pt_400to480_nom,
        histo_pt_480to600_nom,
        histo_pt_600toInf_nom,
    ]

    if normalized:
        for histo in all_histos:
            normalize_histo(histo)

    maximum = get_maximum(all_histos)
    histo_pt_200to250_nom.SetMaximum(maximum * 1.15)
    histo_pt_200to250_nom.SetMinimum(0)

    if normalized:
        histo_pt_200to250_nom.GetYaxis().SetTitle('Fraction of events / 5 GeV')
    else:
        histo_pt_200to250_nom.GetYaxis().SetTitle('Events / 5 GeV')
    histo_pt_200to250_nom.GetYaxis().SetTitleSize(0.045)
    histo_pt_200to250_nom.GetYaxis().SetTitleOffset(1.55)
    # histo_pt_200to250_nom.GetYaxis().CenterTitle()

    histo_pt_200to250_nom.GetXaxis().SetRangeUser(0, max_value)
    # histo_pt_200to250_nom.GetXaxis().SetLimits(0, max_value)

    histo_pt_200to250_nom.GetXaxis().SetTitle('Probe jet #it{m}_{jet} [GeV]')
    histo_pt_200to250_nom.GetXaxis().SetTitleSize(0.045)
    histo_pt_200to250_nom.GetXaxis().SetTitleOffset(1.2)
    # histo_pt_200to250_nom.GetXaxis().CenterTitle()

    histo_pt_200to250_nom.GetXaxis().ChangeLabel(7,-1,-1,-1,-1,-1,'#geq 300');


    leg_x1 = 0.44 - 0.25
    leg_y1 = 0.25
    leg_x2 = 0.77 - 0.25
    leg_y2 = 0.69

    legend = root.TLegend(leg_x1, leg_y1, leg_x2, leg_y2)

    legend.AddEntry(histo_pt_200to250_nom, "200 < #it{p}_{T}/GeV < 250", "l")
    legend.AddEntry(histo_pt_250to300_nom, "250 < #it{p}_{T}/GeV < 300", "l")
    legend.AddEntry(histo_pt_300to350_nom, "300 < #it{p}_{T}/GeV < 350", "l")
    legend.AddEntry(histo_pt_350to400_nom, "350 < #it{p}_{T}/GeV < 400", "l")
    legend.AddEntry(histo_pt_400to480_nom, "400 < #it{p}_{T}/GeV < 480", "l")
    legend.AddEntry(histo_pt_480to600_nom, "480 < #it{p}_{T}/GeV < 600", "l")
    legend.AddEntry(histo_pt_600toInf_nom, "#it{p}_{T} > 600 GeV", "l")

    legend.SetTextSize(0.03)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    legend.Draw()


    # header_leg2 = '#bf{Fully merged top quark, t#bar{t} MC}'
    header_leg2 = '#bf{t#bar{t} MC, fully merged t}'
    leg2_x1 = 0.44 - 0.25
    leg2_y1 = 0.72 - 0.08 + 0.1
    leg2_x2 = 0.77 - 0.25
    leg2_y2 = 0.87 - 0.08
    legend2 = root.TLegend(leg2_x1, leg2_y1, leg2_x2, leg2_y2)

    histo_pt_200to250_nom_leg2 = deepcopy(histo_pt_200to250_nom)
    histo_pt_200to250_nom_leg2.SetLineColor(root.kWhite)

    legend2.SetHeader(header_leg2)
    #legend2.AddEntry(histo_pt_200to250_nom_leg2, "")
    #legend2.AddEntry(histo_pt_200to250_nom_leg2, "")

    legend2.SetTextSize(0.03)
    legend2.SetBorderSize(0)
    legend2.SetFillStyle(0)

    legend2.Draw()


    text_pw1 = root.TLatex(coord.graph_to_pad_x(0.05), coord.graph_to_pad_y(0.95), 'Private work')
    text_pw1.SetTextAlign(13)
    text_pw1.SetTextFont(52)
    text_pw1.SetTextSize(0.03)
    text_pw1.SetNDC()
    text_pw1.Draw()

    text_pw2 = root.TLatex(coord.graph_to_pad_x(0.05), coord.graph_to_pad_y(0.90), '(CMS data/simulation)')
    text_pw2.SetTextAlign(13)
    text_pw2.SetTextFont(52)
    text_pw2.SetTextSize(0.03)
    text_pw2.SetNDC()
    text_pw2.Draw()


    tlatex_probejetalgo = root.TLatex(coord.graph_to_pad_x(0.95), coord.graph_to_pad_y(0.95), 'HOTVR PUPPI')
    tlatex_probejetalgo.SetTextAlign(33) # right top
    tlatex_probejetalgo.SetTextFont(62)
    tlatex_probejetalgo.SetTextSize(0.05)
    tlatex_probejetalgo.SetTextColor(root.kGray+1)
    tlatex_probejetalgo.SetNDC()
    tlatex_probejetalgo.Draw()


    y_shift = 0.031

    # tlatex_pt1 = '#tau_{3}/#tau_{2} < ' + str(wp.get_cut_value())# + '(' + wp.alt_name + ')'
    tlatex_pt1 = 'HOTVR t tag'
    tlatex_pt1 = root.TLatex(coord.graph_to_pad_x(0.95), coord.graph_to_pad_y(0.867 - y_shift), tlatex_pt1)
    tlatex_pt1.SetTextAlign(31) # left top
    tlatex_pt1.SetTextFont(42)
    tlatex_pt1.SetTextSize(0.035)
    tlatex_pt1.SetNDC()
    tlatex_pt1.Draw()

    tlatex_pt2 = '#tau_{3}/#tau_{2} < 0.56'
    tlatex_pt2 = root.TLatex(coord.graph_to_pad_x(0.95), coord.graph_to_pad_y(0.801 - y_shift), tlatex_pt2)
    tlatex_pt2.SetTextAlign(31) # left top
    tlatex_pt2.SetTextFont(42)
    tlatex_pt2.SetTextSize(0.035)
    tlatex_pt2.SetNDC()
    tlatex_pt2.Draw()

    tlatex_pt3 = '|#eta| < 2.5'
    tlatex_pt3 = root.TLatex(coord.graph_to_pad_x(0.95), coord.graph_to_pad_y(0.735 - y_shift), tlatex_pt3)
    tlatex_pt3.SetTextAlign(31) # left top
    tlatex_pt3.SetTextFont(42)
    tlatex_pt3.SetTextSize(0.035)
    tlatex_pt3.SetNDC()
    tlatex_pt3.Draw()


    text_top_left = 'T&P e/#mu+jets'
    # text_top_right = '#sqrt{#it{s}} = 13 TeV'
    text_top_right = str(_YEARS['run2']['lumi_fb_display']) + ' fb^{#minus1} (13 TeV)'

    if text_top_left:
        tlatex_top_left = root.TLatex(margin_l, 1. - margin_t + 0.01, text_top_left)
        tlatex_top_left.SetTextAlign(11)  # left bottom
        tlatex_top_left.SetTextFont(42)
        tlatex_top_left.SetTextSize(text_size)
        tlatex_top_left.SetNDC()
        tlatex_top_left.Draw()

    if text_top_right:
        tlatex_top_right = root.TLatex(1. - margin_r, 1. - margin_t + 0.01, text_top_right)
        tlatex_top_right.SetTextAlign(31)  # right bottom
        tlatex_top_right.SetTextFont(42)
        tlatex_top_right.SetTextSize(text_size)
        tlatex_top_right.SetNDC()
        tlatex_top_right.Draw()


    canvas.Update()
    tickScaleX = (canvas.GetUxmax() - canvas.GetUxmin()) / (canvas.GetX2() - canvas.GetX1()) * (canvas.GetWh() * canvas.GetAbsHNDC())
    tickScaleY = (canvas.GetUymax() - canvas.GetUymin()) / (canvas.GetY2() - canvas.GetY1()) * (canvas.GetWw() * canvas.GetAbsWNDC())
    histo_pt_200to250_nom.GetXaxis().SetTickLength(canvas.GetWh() * tick_length / tickScaleX)
    histo_pt_200to250_nom.GetYaxis().SetTickLength(canvas.GetWw() * tick_length / tickScaleY)

    root.gPad.Update()
    root.gPad.RedrawAxis()

    output_path = 'plots_pass_region'
    os.system('mkdir -p '+output_path)

    output_file_name = 'plotPassRegion_hotvr_' + wp.name
    if normalized:
        output_file_name += '_norm'
    output_file_name += '.pdf'
    output_file_path = os.path.join(output_path, output_file_name)

    canvas.SaveAs(output_file_path)


if __name__ == "__main__":
    #create_plot_ak8(
    #    wp=_NULL_WP,
    #)

    for i_wp in range(5):
        create_plot_ak8(
            wp=taggers['ak8_t__tau'].get_wp(wp_index=i_wp),
            # normalized=True,
        )

    create_plot_hotvr(
        wp=taggers['hotvr_t__tau'].get_wp(wp_index=0),
        # normalized=True,
    )
