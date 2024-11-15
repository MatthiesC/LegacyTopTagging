#! /usr/bin/env python
from __future__ import print_function

import os
import sys
import ROOT as root
from collections import OrderedDict
import numpy as np
import csv

sys.path.append(os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/Analysis'))
from constants import _YEARS, _TAGGERS, _NULL_WP

sys.path.append(os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/NicePlots/python'))
from plotter import CoordinateConverter

years = [
'UL16preVFP',
'UL16postVFP',
'UL17',
'UL18',
]
# regions = [
# 'Pass',
# 'Fail',
# ]
# channels = [
# 'muo',
# 'ele',
# ]
taggers = [
# 'ak8_t__tau',
# 'ak8_t_btagDJet__tau',
# 'ak8_t_btagDCSV__tau',
# 'hotvr_t__tau',
# 'ak8_w__partnet',
# 'ak8_t__MDdeepak8',
'ak8_t__leonidasTvsQCD',
'ak8_w__leonidasWvsQCD',
'ak8_w__leonidasWvsQCDMD',
]
taggers = {k: _TAGGERS[k] for k in taggers}

graph_types = [
'tot',
# 'stat',
# 'syst',
# 'uncorr',
# 'corr',
]

golden_ratio = .61803398875

# y_scaling_msc = {
#     'FullyMerged': 1.0,
#     'SemiMerged': ,
#     'NotMerged': ,
# }

class ValueAndAsymmErrors():

    def __init__(self, value, up=None, down=None):

        self.value = value
        self.up = up or 0.
        self.down = down or 0.

dummyScaleFactor = ValueAndAsymmErrors(value=1.)

class ScaleFactor():

    def __init__(self, tot=None, stat=None, syst=None, uncorr=None, corr=None):

        self.tot = tot or dummyScaleFactor
        self.stat = stat or dummyScaleFactor
        self.syst = syst or dummyScaleFactor
        self.uncorr = uncorr or dummyScaleFactor
        self.corr = corr or dummyScaleFactor

class ScaleFactorCollector():

    def __init__(self, taggerName, workingPoint=None):

        self.tagger = taggers[taggerName]

        self.years = years

        # self.workdirBasePath = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel/combine/', tagger.name)
        # self.outputBasePath = os.path.join(self.workdirBasePath, 'scale_factors')

        self.workdirBasePath = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/Analysis/Combine/scaleFactors_Leonidas')
        self.outputBasePath = self.workdirBasePath

        # if self.tagger.name.startswith('ak8_t'):
        #     self.mscSplitting = 'mscTop3'
        # elif self.tagger.name.startswith('ak8_w'):
        #     self.mscSplitting = 'mscW3'
        # elif self.tagger.name.startswith('hotvr_t'):
        #     self.mscSplitting = 'mscTop3'

        # if self.mscSplitting == 'mscTop3':
        #     self.mscSplits = [
        #         'FullyMerged',
        #         'SemiMerged',
        #         'NotMerged',
        #     ]
        # elif self.mscSplitting == 'mscW3':
        #     self.mscSplits = [
        #         'FullyMerged',
        #         'WMerged',
        #         'NotTopOrWMerged',
        #     ]
        # else:
        #     sys.exit()

        # get pt bins for tagger
        pt_bins = self.tagger.var_intervals
        self.sorted_pt_bins = []
        for pt_bin in pt_bins.values():
            if not pt_bin.fit == True: continue
            self.sorted_pt_bins.append(pt_bin)
        self.sorted_pt_bins = sorted(self.sorted_pt_bins, key=lambda x: x.var_min, reverse=False)

        self.sf_task_name = '-'.join(['scaleFactors', self.tagger.name])
        self.outputPath = os.path.join(self.outputBasePath, self.sf_task_name)
        print('Writing to', self.outputPath)
        os.system('mkdir -p '+self.outputPath)

        # init scale factor dict
        self.scale_factors = OrderedDict()
        for wp in self.tagger.get_wp():
            for pt_bin in self.sorted_pt_bins:
                for year in self.years:
                    self.scale_factors.setdefault(wp.name, OrderedDict()).setdefault(pt_bin.name, OrderedDict())[year] = ScaleFactor()

    def write_root_files(self):

        # outfile_name_combined = self.sf_task_name+'.root'
        # outFilePath_combined = os.path.join(self.outputPath, outfile_name_combined)
        # print(outFilePath_combined)
        # outFile_combined = root.TFile.Open(outFilePath_combined, 'RECREATE')

        infile_name = 'infileNameNotFound'
        if self.tagger.name == 'ak8_t__leonidasTvsQCD':
            infile_name = 'scaleFactors_Leonidas_TvsQCD.csv'
        elif self.tagger.name == 'ak8_w__leonidasWvsQCD':
            infile_name = 'scaleFactors_Leonidas_WvsQCD.csv'
        elif self.tagger.name == 'ak8_w__leonidasWvsQCDMD':
            infile_name = 'scaleFactors_Leonidas_WvsQCDMD.csv'

        inFilePath = os.path.join(self.outputBasePath, infile_name)

        # Reading the CSV file
        with open(inFilePath, mode='r') as file:
            csv_reader = csv.reader(file)
            data = list(csv_reader)

        header = data[0]
        data = data[1:]

        for year in self.years:

            value_Year = _YEARS.get(year)['leonidasNaming']
            # print(value_Year)

            # outfile_name = self.sf_task_name+'-'+year+'.root'
            # outFilePath = os.path.join(self.outputPath, outfile_name)
            # print(outFilePath)

            # outFile = root.TFile.Open(outFilePath, 'RECREATE')
            # outFolderName = '-'.join([self.tagger.name, year])
            # outFolder = outFile.mkdir(outFolderName)
            # outFolder_combined = outFile_combined.mkdir(outFolderName)

            for wp in self.tagger.get_wp():

                value_BkgEff = '{:.1f}'.format(wp.bkg_eff * 100)
                # print(value_BkgEff)

                for graph_type in graph_types:

                    # graph_name = '_'.join([wp.name, year, graph_type])

                    # input_graphs = []

                    # Create an empty TGraphAsymmErrors object to store the concatenated data
                    # concatenated_graph = root.TGraphAsymmErrors()

                    for pt_bin in self.sorted_pt_bins:

                        # combine_task_name_suffix = '-'.join([('PtTotal' if pt_bin.total_range else pt_bin.name), year]) # HACK
                        # combine_task_name = '-'.join(['combineTask', self.tagger.name, self.wp.name]) + ('-'+combine_task_name_suffix if len(combine_task_name_suffix) else '')
                        # combine_task_path = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel/combine', self.tagger.name, 'workdirs', combine_task_name)
                        # infile_name = 'scale_factors_obs.root'
                        # inFilePath = os.path.join(combine_task_path, infile_name)
                        
                        value_PtLow = str(int(pt_bin.var_min))
                        # print(value_PtLow)

                        for row in data:
                            this_row = True
                            if row[3] != value_BkgEff:
                                this_row = False
                            if row[4] != value_PtLow:
                                this_row = False
                            if row[1] != value_Year:
                                this_row = False
                            if not this_row:
                                continue
                            # code here
                            y = float(row[-3])
                            y_err_low = float(row[-2])
                            y_err_high = float(row[-1])
                            setattr(self.scale_factors[wp.name][pt_bin.name][year], graph_type, ValueAndAsymmErrors(value=y, up=y_err_high, down=y_err_low))
                            break
                            


                        # inFile = root.TFile.Open(inFilePath, 'READ')
                        # inFile.SetDirectory(0)

                        # input_graph = inFile.Get(graph_name)
                        # input_graphs.append(input_graph)

                        # Get the number of points in the input graph
                        # num_points = input_graph.GetN()
                        # if num_points != 1:
                        #     sys.exit('input graph should only have one entry')

                        # Loop through each point in the input graph and add it to the concatenated graph
                        # for i in range(num_points):
                        #     x = root.Double()
                        #     y = root.Double()
                        #     x_err_low = root.Double()
                        #     x_err_high = root.Double()
                        #     y_err_low = root.Double()
                        #     y_err_high = root.Double()

                        #     # Get the values and errors of the point
                        #     input_graph.GetPoint(i, x, y)
                        #     x_err_low = input_graph.GetErrorXlow(i)
                        #     x_err_high = input_graph.GetErrorXhigh(i)
                        #     y_err_low = input_graph.GetErrorYlow(i)
                        #     y_err_high = input_graph.GetErrorYhigh(i)

                        #     # setattr(x, 'y', v) is equivalent to x.y = v
                        #     setattr(self.scale_factors[wp.name][pt_bin.name][year], graph_type, ValueAndAsymmErrors(value=y, up=y_err_high, down=y_err_low))
                        #     # self.scale_factors[msc][pt_bin.name][year].tot.value = y
                        #     # self.scale_factors[msc][pt_bin.name][year].tot.up = y_err_high
                        #     # self.scale_factors[msc][pt_bin.name][year].tot.down = y_err_low

                        #     # Add the point with errors to the concatenated graph
                        #     # concatenated_graph.SetPoint(concatenated_graph.GetN(), x, y)
                        #     # concatenated_graph.SetPointError(concatenated_graph.GetN() - 1, x_err_low, x_err_high, y_err_low, y_err_high)

                        # inFile.Close()

                    # # Create an empty TGraphAsymmErrors object to store the concatenated data
                    # concatenated_graph = root.TGraphAsymmErrors()
                    #
                    # # Loop through each input TGraphAsymmErrors object and retrieve its data points and errors
                    # for input_graph in input_graphs:
                    #     # Get the number of points in the input graph
                    #     num_points = input_graph.GetN()
                    #     if num_points != 1:
                    #         sys.exit('input graph should only have one entry')
                    #
                    #     # Loop through each point in the input graph and add it to the concatenated graph
                    #     for i in range(num_points):
                    #         x = root.Double()
                    #         y = root.Double()
                    #         x_err_low = root.Double()
                    #         x_err_high = root.Double()
                    #         y_err_low = root.Double()
                    #         y_err_high = root.Double()
                    #
                    #         # setattr(x, 'y', v) is equivalent to x.y = v
                    #         setattr(self.scale_factors[msc][pt_bin.name][year], graph_type, ValueAndAsymmErrors(value=y, up=y_err_high, down=y_err_low))
                    #         # self.scale_factors[msc][pt_bin.name][year].tot.value = y
                    #         # self.scale_factors[msc][pt_bin.name][year].tot.up = y_err_high
                    #         # self.scale_factors[msc][pt_bin.name][year].tot.down = y_err_low
                    #
                    #         # Get the values and errors of the point
                    #         input_graph.GetPoint(i, x, y)
                    #         x_err_low = input_graph.GetErrorXlow(i)
                    #         x_err_high = input_graph.GetErrorXhigh(i)
                    #         y_err_low = input_graph.GetErrorYlow(i)
                    #         y_err_high = input_graph.GetErrorYhigh(i)
                    #
                    #         # Add the point with errors to the concatenated graph
                    #         concatenated_graph.SetPoint(concatenated_graph.GetN(), x, y)
                    #         concatenated_graph.SetPointError(concatenated_graph.GetN() - 1, x_err_low, x_err_high, y_err_low, y_err_high)

                    # outFolder.cd()
                    # concatenated_graph.Write(graph_name)
                    # outFolder_combined.cd()
                    # concatenated_graph.Write(graph_name)

            # print('Wrote', outFilePath)

        # print('Wrote', outFilePath_combined)

    def plot_scale_factors(self, dynamic_width=False):

        self.text_size = 0.035

        margin_b = 0.15
        margin_t = 0.05
        margin_l = 0.15
        margin_r = 0.05

        canvas_name = 'canvas-'+self.sf_task_name
        canvas_title = self.sf_task_name
        canvas_height = 600
        if dynamic_width:
            canvas_width_per_pt_bin = 160 # 200
            canvas_width = int(canvas_height * (margin_b + margin_t) + canvas_width_per_pt_bin * len(self.sorted_pt_bins))
        else:
            canvas_width = 750
        self.canvas = root.TCanvas(canvas_name, canvas_title, canvas_width, canvas_height)
        # hw_correction = self.canvas.GetWh() / self.canvas.GetWw()
        hw_correction = float(canvas_height) / float(canvas_width) # if not using float(), it will be an integer division which can lead to margin_correction = 0
        margin_l = margin_l * hw_correction
        margin_r = margin_r * hw_correction
        self.canvas.SetMargin(margin_l, margin_r, margin_b, margin_t); # lrbt

        self.coord = CoordinateConverter(margin_l, margin_r, margin_b, margin_t)

        header_y_fraction = 0.2
        y_sf_low = 0.001
        y_sf_high = 1.999
        y_sf_range = y_sf_high - y_sf_low
        y_gap = y_sf_range * 0.5 * golden_ratio
        y_low = y_sf_low
        y_high = (y_sf_range + y_gap) * len(self.tagger.get_wp()) / (1. - header_y_fraction)
        y_range = y_high - y_low

        y_sf_offset = OrderedDict()
        for index, wp in enumerate(reversed(self.tagger.get_wp())):
            y_sf_offset[wp.name] = index * (y_sf_range + y_gap)

        x_pt_low = 0.
        x_pt_high = 1.
        x_pt_range = x_pt_high - x_pt_low
        x_pt_margin = x_pt_range * 0.1
        x_pt_nonmargin = x_pt_range - 2. * x_pt_margin
        x_sf_distance = x_pt_nonmargin / len(self.years)
        x_sf_position = OrderedDict()
        for index, year in enumerate(self.years):
            x_sf_position[year] = x_pt_margin + (0.5 + index) * x_sf_distance
        x_sf_width = x_pt_range * 0.1
        x_low = x_pt_low
        x_high = x_pt_range * len(self.sorted_pt_bins)
        x_range = x_high - x_low

        x_pt_offset = OrderedDict()
        for index, pt_bin in enumerate(self.sorted_pt_bins):
            x_pt_offset[pt_bin.name] = index * x_pt_range

        #HACK for W tagger:
        arr_x__w_tagger = np.array([], dtype=np.float32)
        arr_y__w_tagger = np.array([], dtype=np.float32)
        arr_zeroes__w_tagger = np.array([], dtype=np.float32)

        # one graph for each year / each graph_type
        graphs = OrderedDict()
        points_per_graph = len(self.tagger.get_wp()) * len(self.sorted_pt_bins)
        for year in self.years:
            for graph_type in graph_types:
                arr_x = np.array([], dtype=np.float32)
                arr_y = np.array([], dtype=np.float32)
                arr_exl = np.array([], dtype=np.float32)
                arr_exh = np.array([], dtype=np.float32)
                arr_eyl = np.array([], dtype=np.float32)
                arr_eyh = np.array([], dtype=np.float32)
                for wp in self.tagger.get_wp():

                    # SCALE HERE
                    scale_ey = 0.7
                    # if msc == 'FullyMerged':
                    #     scale_ey = 0.4
                    # elif msc == 'SemiMerged':
                    #     scale_ey = 0.7
                    # elif self.mscSplitting == 'mscW3':
                    #     scale_ey = 0.4

                    for pt_bin in self.sorted_pt_bins:
                        arr_x = np.append(arr_x, x_sf_position[year] + x_pt_offset[pt_bin.name])
                        value_y = 0.5 * y_sf_range + ( getattr(self.scale_factors[wp.name][pt_bin.name][year], 'tot').value - 0.5 * y_sf_range ) / scale_ey + y_sf_offset[wp.name] # always use central value of tot # SCALE HERE
                        arr_y = np.append(arr_y, value_y)
                        if graph_type == 'stat':
                            exl = 0.6 * x_sf_width
                        elif graph_type == 'uncorr':
                            exl = 0.0 * x_sf_width
                        else:
                            exl = 0.4 * x_sf_width
                        arr_exl = np.append(arr_exl, exl)
                        arr_exh = np.append(arr_exh, exl)
                        value_eyl = getattr(self.scale_factors[wp.name][pt_bin.name][year], graph_type).down / scale_ey # SCALE HERE
                        arr_eyl = np.append(arr_eyl, value_eyl)
                        value_eyh = getattr(self.scale_factors[wp.name][pt_bin.name][year], graph_type).up / scale_ey # SCALE HERE
                        arr_eyh = np.append(arr_eyh, value_eyh)
                graphs.setdefault(year, OrderedDict())[graph_type] = root.TGraphAsymmErrors(points_per_graph, arr_x, arr_y, arr_exl, arr_exh, arr_eyl, arr_eyh)

        self.canvas.cd()
        # pad_name = 'pad-'+self.sf_task_name
        # pad_title = self.sf_task_name
        # pad = root.TPad(pad_name, pad_title, 0., 0., 1., 1.)
        # pad.SetBottomMargin(margin_b)
        # pad.SetTopMargin(margin_t)
        # pad.SetLeftMargin(margin_l)
        # pad.SetRightMargin(margin_r)
        # pad.SetTickx(1);
        # pad.SetTicky(1);
        # pad.Draw()
        # root.gStyle.SetLineWidth(1)
        # root.gStyle.SetOptStat(root.kFALSE)
        # pad.cd()


        # leg_x1 = 0.46
        # leg_y1 = 0.76
        # leg_x2 = 0.948
        # leg_y2 = 0.81
        # legend = root.TLegend(leg_x1, leg_y1, leg_x2, leg_y2)
        # legend.SetNColumns(len(self.years))
        leg_x1 = 0.015
        leg_y1 = 0.81
        leg_x2 = 0.18
        leg_y2 = 0.98
        legend = root.TLegend(leg_x1, leg_y1, leg_x2, leg_y2, '', 'brNDC')
        # legend.SetNColumns(len(self.years))
        legend.SetTextSize(0.03)
        # legend.SetBorderSize(0)
        # legend.SetFillStyle(0)

        # same color setup as in LegacyTopTagging/Analysis/ControlPlots/plots_HOTVR_ptResponse_v2.cxx
        transparency1 = 1.0
        transparency2 = 0.5

        mg = root.TMultiGraph()
        mg.SetTitle('')
        for year in self.years:
            graph_tot = graphs[year]['tot']
            graph_tot.SetFillColorAlpha(_YEARS.get(year)['tcolor'], transparency1)
            graph_tot.SetMarkerStyle(7)
            graph_tot.SetLineWidth(0)
            legend.AddEntry(graph_tot, _YEARS.get(year)['year'])
            mg.Add(graph_tot)

            # graph_uncorr = graphs[year]['uncorr']
            # graph_uncorr.SetFillColor(root.kWhite)
            # graph_uncorr.SetMarkerStyle(7)
            # mg.Add(graph_uncorr)

            # graph_stat = graphs[year]['stat']
            # graph_stat.SetFillColorAlpha(_YEARS.get(year)['tcolor'], transparency1)
            # graph_stat.SetMarkerStyle(7)
            # graph_stat.SetLineWidth(0)
            # legend.AddEntry(graph_stat, _YEARS.get(year)['year'])
            # mg.Add(graph_stat)
        mg.Draw('a2p')
        mg.GetHistogram().GetXaxis().SetLimits(x_low, x_high)
        mg.GetHistogram().GetXaxis().SetTickLength(0.)
        mg.GetHistogram().GetXaxis().SetLabelOffset(-2.)
        mg.GetHistogram().SetMinimum(y_low)
        mg.GetHistogram().SetMaximum(y_high)
        mg.GetHistogram().GetYaxis().SetTickLength(0.)
        mg.GetHistogram().GetYaxis().SetLabelOffset(-2.)

        lines = {}
        for index, wp in enumerate(self.tagger.get_wp()):
            shifted_one = 1. + y_sf_offset[wp.name]
            lines[wp.name] = root.TLine(x_low, shifted_one, x_high, shifted_one)
            lines[wp.name].SetLineColor(root.kGray + 1)
            lines[wp.name].SetLineStyle(3)
            lines[wp.name].Draw()


        end_error_size = 5 # this is a number of pixels. This could be dynamic depending on the size of the canvas and number of pt bins
        root.gStyle.SetEndErrorSize(end_error_size)

        # mg2 = root.TMultiGraph()
        # mg2.SetTitle('')
        # for year in self.years:
        #     graph_uncorr = graphs[year]['uncorr']
        #     # graph_uncorr.SetFillColor(root.kWhite)
        #     graph_uncorr.SetLineColor(root.kBlack)
        #     graph_uncorr.SetMarkerStyle(7)
        #     mg2.Add(graph_uncorr)
        # mg2.Draw('||p')

        legend.Draw()

        #HACK W tagger:
        if self.tagger.name == 'ak8_w__partnet':
            graph__w_tagger = root.TGraphAsymmErrors(len(arr_x__w_tagger), arr_x__w_tagger, arr_y__w_tagger, arr_zeroes__w_tagger, arr_zeroes__w_tagger, arr_zeroes__w_tagger, arr_zeroes__w_tagger)
            graph__w_tagger.SetMarkerStyle(47)
            graph__w_tagger.SetMarkerSize(1.2)
            graph__w_tagger.Draw('p')

        

        axes = {}
        ndiv = 505
        tick_length = 0.05
        w_low_base = y_sf_low
        w_high_base = y_sf_high
        for index_wp, wp in enumerate(self.tagger.get_wp()):

            # SCALE HERE
            scale_ey = 0.7
            # if msc == 'FullyMerged':
            #     scale_ey = 0.4
            # elif msc == 'SemiMerged':
            #     scale_ey = 0.7
            # elif self.mscSplitting == 'mscW3':
            #     scale_ey = 0.4

            w_low = 1 + ( w_low_base - 1 ) * scale_ey # SCALE HERE
            w_high = 1 + ( w_high_base - 1 ) * scale_ey # SCALE HERE
            x_axis_low = x_low
            y_axis_low = y_sf_low + y_sf_offset[wp.name]
            x_axis_high = x_low
            y_axis_high = y_sf_high + y_sf_offset[wp.name]
            axes.setdefault(wp.name, {})['left'] = root.TGaxis(x_axis_low, y_axis_low, x_axis_high, y_axis_high, w_low, w_high, ndiv, '-S')
            axes[wp.name]['left'].SetTickLength(tick_length)
            axes[wp.name]['left'].SetLabelSize(self.text_size)
            axes[wp.name]['left'].SetLabelFont(42)
            axes[wp.name]['left'].Draw()
            for index_pt_bin, pt_bin in enumerate(self.sorted_pt_bins):
                if index_pt_bin == 0:
                    continue
                x_axis_low = x_pt_low + x_pt_offset[pt_bin.name]
                x_axis_high = x_pt_low + x_pt_offset[pt_bin.name]
                axes.setdefault(wp.name, {})[pt_bin.name] = root.TGaxis(x_axis_low, y_axis_low, x_axis_high, y_axis_high, w_low, w_high, ndiv, 'U+-S')
                axes[wp.name][pt_bin.name].SetTickLength(tick_length)
                axes[wp.name][pt_bin.name].Draw()
            x_axis_low = x_high
            x_axis_high = x_high
            axes.setdefault(wp.name, {})['right'] = root.TGaxis(x_axis_low, y_axis_low, x_axis_high, y_axis_high, w_low, w_high, ndiv, 'U+S')
            axes[wp.name]['right'].SetTickLength(tick_length)
            axes[wp.name]['right'].Draw()

        deltaX = 0.03
        hw_corrected_deltaX = deltaX * hw_correction

        text_cms_offset_x = 0.055

        self.tlatex_cms = root.TLatex(self.coord.graph_to_pad_x(hw_corrected_deltaX) + text_cms_offset_x, self.coord.graph_to_pad_y(1. - deltaX), 'CMS')
        self.tlatex_cms.SetTextAlign(13) # left top
        self.tlatex_cms.SetTextFont(62)
        self.tlatex_cms.SetTextSize(0.07)
        self.tlatex_cms.SetNDC()
        self.tlatex_cms.Draw()

        # self.text_prelim = 'Work in Progress'
        self.text_prelim = 'Preliminary'
        # #HACK
        # if self.tagger.name == 'ak8_w__partnet' or self.tagger.name == 'ak8_t__MDdeepak8':
        #     self.text_prelim = 'Private Work'
        #HACK end
        self.tlatex_prelim = root.TLatex(self.coord.graph_to_pad_x(hw_corrected_deltaX + 0.14) + text_cms_offset_x, self.coord.graph_to_pad_y(0.948), self.text_prelim)
        self.tlatex_prelim.SetTextAlign(13) # left top
        self.tlatex_prelim.SetTextFont(52)
        self.tlatex_prelim.SetTextSize(0.045)
        self.tlatex_prelim.SetNDC()
        self.tlatex_prelim.Draw()

        probejetalgo = 'AK8'
        # if tagger.name.startswith('ak8_'):
        #     probejetalgo = 'AK8'
        # elif tagger.name.startswith('hotvr_'):
        #     probejetalgo = 'HOTVR'

        # if self.wp.alt_name is not None:
        #     pt_text_string = self.tagger.label.replace('{WP_VALUE}', '{} ({})'.format(self.wp.get_cut_value(year), self.wp.alt_name))
        # else:
        #     pt_text_string = self.tagger.label.replace('{WP_VALUE}', '{}'.format(self.wp.get_cut_value(year)))
        # pt_text_string2 = '|#eta| < 2.4'
        # pt_text_string += ', ' + pt_text_string2

        # #HACK for ak8_w__partnet and ak8_t__MDdeepak8:
        # if self.tagger.name == 'ak8_w__partnet':
        #     pt_text_string = 'ParticleNet WvsQCD (#varepsilon_{B} = 3%), |#eta| < 2.5'
        # if self.tagger.name == 'ak8_t__MDdeepak8':
        #     pt_text_string = 'MD-DeepAK8 TvsQCD (#varepsilon_{B} = 1%), |#eta| < 2.5'

        if self.tagger.name == 'ak8_t__leonidasTvsQCD':
            # pt_text_string = 'ParticleNet TvsQCD, fully merged t decays, |#eta| < 2.4'
            pt_text_string = 'ParticleNet TvsQCD'
            pt_text_string2 = 'Fully merged t decays, |#eta| < 2.4'
        elif self.tagger.name == 'ak8_w__leonidasWvsQCD':
            # pt_text_string = 'ParticleNet WvsQCD, fully merged W decays, |#eta| < 2.4'
            pt_text_string = 'ParticleNet WvsQCD'
            pt_text_string2 = 'Fully merged W decays, |#eta| < 2.4'
        elif self.tagger.name == 'ak8_w__leonidasWvsQCDMD':
            # pt_text_string = 'ParticleNet-MD WvsQCD-MD, fully merged W decays, |#eta| < 2.4'
            pt_text_string = 'ParticleNet-MD WvsQCD-MD'
            pt_text_string2 = 'Fully merged W decays, |#eta| < 2.4'

        #HACK end

        tlatex_pt = root.TLatex(self.coord.graph_to_pad_x(1. - hw_corrected_deltaX), self.coord.graph_to_pad_y(0.848), pt_text_string)
        tlatex_pt.SetTextAlign(31) # left top
        tlatex_pt.SetTextFont(42)
        tlatex_pt.SetTextSize(0.03)
        tlatex_pt.SetNDC()
        tlatex_pt.Draw()

        tlatex_pt2 = root.TLatex(self.coord.graph_to_pad_x(1. - hw_corrected_deltaX), self.coord.graph_to_pad_y(0.79), pt_text_string2)
        tlatex_pt2.SetTextAlign(31) # left top
        tlatex_pt2.SetTextFont(42)
        tlatex_pt2.SetTextSize(0.03)
        tlatex_pt2.SetNDC()
        tlatex_pt2.Draw()

        tlatex_probejetalgo = root.TLatex(self.coord.graph_to_pad_x(1. - hw_corrected_deltaX), self.coord.graph_to_pad_y(1. - deltaX), probejetalgo+' PUPPI')
        tlatex_probejetalgo.SetTextAlign(33) # right top
        tlatex_probejetalgo.SetTextFont(62)
        tlatex_probejetalgo.SetTextSize(0.07)
        tlatex_probejetalgo.SetTextColor(root.kGray+1)
        tlatex_probejetalgo.SetNDC()
        tlatex_probejetalgo.Draw()

        pt_labels = {}
        for index, pt_bin in enumerate(self.sorted_pt_bins):
            if pt_bin.max_set:
                pt_text = '#it{p}_{T} #in ({pt_min}, {pt_max}] GeV'.replace('{pt_min}', '{}'.format(pt_bin.var_min)).replace('{pt_max}', '{}'.format(pt_bin.var_max))
            else:
                pt_text = '#it{p}_{T} > {pt_min} GeV'.replace('{pt_min}', '{}'.format(pt_bin.var_min))
            # pt_label_position_x = margin_l + (((len(self.years) + 1) * index) / x_high) * (1 - margin_l - margin_r) + 0.015
            pt_label_position_x = margin_l + 0.015 + (1. - margin_l - margin_r) / len(self.sorted_pt_bins) * index
            # print(pt_label_position_x)
            # print(pt_text)
            pt_label_position_y = 0.9 * margin_b
            pt_labels[pt_bin.name] = root.TLatex(pt_label_position_x, pt_label_position_y, pt_text)
            pt_labels[pt_bin.name].SetTextAlign(13) # left top aligned
            if len(self.sorted_pt_bins) > 4:
                pt_labels[pt_bin.name].SetTextAngle(-20)
            else:
                pt_labels[pt_bin.name].SetTextAngle(-15)
            pt_labels[pt_bin.name].SetTextFont(42)
            pt_labels[pt_bin.name].SetTextSize(self.text_size)
            pt_labels[pt_bin.name].SetNDC()
            pt_labels[pt_bin.name].Draw()

        y_label = root.TLatex(
            # margin_l * 0.3, # margin_l * 0.45
            margin_l * 0.45, # margin_l * 0.45
            (1. - (y_gap / y_range + header_y_fraction)) / 2. * (1. - margin_t - margin_b) + margin_b,
            'Data-to-simulation scale factor'
        )
        y_label.SetTextAlign(21)
        y_label.SetTextAngle(90)
        y_label.SetTextFont(42)
        y_label.SetTextSize(self.text_size)
        y_label.SetNDC()
        y_label.Draw()

        # y_label2 = root.TLatex(
        #     margin_l * 0.5,
        #     (1. - (y_gap / y_range + header_y_fraction)) / 2. * (1. - margin_t - margin_b) + margin_b,
        #     '[ inner (outer) colored area = stat. (tot.) unc. / black bar = era-uncorrelated unc. ]' #FIXME
        # )
        # y_label2.SetTextAlign(21)
        # y_label2.SetTextAngle(90)
        # y_label2.SetTextFont(42)
        # y_label2.SetTextSize(0.02)
        # y_label2.SetNDC()
        # y_label2.Draw()




        msc_paves = {}
        msc_labels = {}
        for wp in self.tagger.get_wp():
            y_center = 1. + y_sf_offset[wp.name]
            y_center += 0.03 # small correction for shadow of pave
            pave_distance_to_scale_factor_range = 0.1
            msc_paves[wp.name] = root.TPave(
                # self.coord.graph_to_pad_x(0.04 * hw_correction),
                # self.coord.graph_to_pad_y((y_center + 1. + y_gap * (1. - pave_distance_to_scale_factor_range)) / y_range),
                # self.coord.graph_to_pad_x(0.34 * hw_correction),
                # self.coord.graph_to_pad_y((y_center + 1. + y_gap * pave_distance_to_scale_factor_range) / y_range),
                (0.04 * hw_correction) * x_range,
                ((y_center + 1. + y_gap * (1. - pave_distance_to_scale_factor_range)) / y_range) * y_range,
                (0.44 * hw_correction) * x_range,
                ((y_center + 1. + y_gap * pave_distance_to_scale_factor_range) / y_range) * y_range,
                3,
                ''
            )
            msc_paves[wp.name].SetLineWidth(0)
            msc_paves[wp.name].Draw()

            msc_text = 'Mistag rate = {:.1f}%'.format(wp.bkg_eff * 100)

            y_center += 0.05

            msc_labels[wp.name] = root.TText(
                0.24 * hw_correction * x_range,
                ((y_center + 1. + y_gap * pave_distance_to_scale_factor_range + 0.06) / y_range) * y_range,
                msc_text
            )
            msc_labels[wp.name].SetTextAlign(21)
            msc_labels[wp.name].SetTextFont(52)
            msc_labels[wp.name].SetTextSize(0.03)
            msc_labels[wp.name].Draw()


        text_top_left_offset_x = 0.075
        text_top_left = 'Run II Legacy'
        self.tlatex_top_left = root.TLatex(margin_l + text_top_left_offset_x, 1. - margin_t + 0.01, text_top_left)
        self.tlatex_top_left.SetTextAlign(11) # left bottom
        self.tlatex_top_left.SetTextFont(42)
        self.tlatex_top_left.SetTextSize(self.text_size)
        self.tlatex_top_left.SetNDC()
        self.tlatex_top_left.Draw()

        text_top_right = ', '.join([_YEARS.get(year).get('lumi_fb_display')+' fb^{#minus1}' for year in self.years])
        text_top_right += ' (13 TeV)'
        self.tlatex_top_right = root.TLatex(1. - margin_r, 1. - margin_t + 0.01, text_top_right)
        self.tlatex_top_right.SetTextAlign(31) # right bottom
        self.tlatex_top_right.SetTextFont(42)
        self.tlatex_top_right.SetTextSize(self.text_size)
        self.tlatex_top_right.SetNDC()
        self.tlatex_top_right.Draw()

        for fileType in ['.pdf', '.png']:
            plotFileName = self.sf_task_name+'-RescaledYAxis'+fileType
            plotFilePath = os.path.join(self.outputPath, plotFileName)
            self.canvas.SaveAs(plotFilePath)

if __name__ == '__main__':

    for tagger in taggers.values():

        sf = ScaleFactorCollector(taggerName=tagger.name)
        sf.write_root_files()
        sf.plot_scale_factors()
        # print(sf.scale_factors)

        # for wp_index in range(len(tagger.get_wp())):

        #     sf = ScaleFactorCollector(taggerName=tagger.name, workingPoint=wp_index)
        #     sf.write_root_files()
        #     sf.plot_scale_factors()
        #     # print(sf.scale_factors)
