from __future__ import print_function

import os
import sys
import ROOT as root
from collections import OrderedDict
import numpy as np

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
'ak8_t__MDdeepak8',
]
taggers = {k: _TAGGERS[k] for k in taggers}

graph_types = [
'tot',
'stat',
'syst',
'uncorr',
'corr',
]

golden_ratio = .61803398875

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

        self.workdirBasePath = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel/combine/', tagger.name)
        self.outputBasePath = os.path.join(self.workdirBasePath, 'scale_factors')

        if self.tagger.name.startswith('ak8_t'):
            self.mscSplitting = 'mscTop3'
        elif self.tagger.name.startswith('ak8_w'):
            self.mscSplitting = 'mscW3'
        elif self.tagger.name.startswith('hotvr_t'):
            self.mscSplitting = 'mscTop3'

        if self.mscSplitting == 'mscTop3':
            self.mscSplits = [
                'FullyMerged',
                'SemiMerged',
                'NotMerged',
            ]
        elif self.mscSplitting == 'mscW3':
            self.mscSplits = [
                'FullyMerged',
                'WMerged',
                'NotTopOrWMerged',
            ]
        else:
            sys.exit()

        # get pt bins for tagger
        pt_bins = self.tagger.var_intervals
        self.sorted_pt_bins = []
        for pt_bin in pt_bins.values():
            if not pt_bin.fit == True: continue
            self.sorted_pt_bins.append(pt_bin)
        self.sorted_pt_bins = sorted(self.sorted_pt_bins, key=lambda x: x.var_min, reverse=False)

        # get working point
        if isinstance(workingPoint, str):
            self.wp = self.tagger.get_wp(wp_name=workingPoint)
        elif isinstance(workingPoint, int):
            self.wp = self.tagger.get_wp(wp_index=workingPoint)
        else:
            sys.exit()

        self.sf_task_name = '-'.join(['scaleFactors', self.tagger.name, self.wp.name])
        self.outputPath = os.path.join(self.outputBasePath, self.sf_task_name)
        print('Writing to', self.outputPath)
        os.system('mkdir -p '+self.outputPath)

        # init scale factor dict
        self.scale_factors = OrderedDict()
        for msc in self.mscSplits:
            for pt_bin in self.sorted_pt_bins:
                for year in self.years:
                    self.scale_factors.setdefault(msc, OrderedDict()).setdefault(pt_bin.name, OrderedDict())[year] = ScaleFactor()

    def write_root_files(self):

        outfile_name_combined = self.sf_task_name+'.root'
        outFilePath_combined = os.path.join(self.outputPath, outfile_name_combined)
        # print(outFilePath_combined)
        outFile_combined = root.TFile.Open(outFilePath_combined, 'RECREATE')

        for year in self.years:

            outfile_name = self.sf_task_name+'-'+year+'.root'
            outFilePath = os.path.join(self.outputPath, outfile_name)
            # print(outFilePath)

            outFile = root.TFile.Open(outFilePath, 'RECREATE')
            outFolderName = '-'.join([self.tagger.name, self.wp.name, year])
            outFolder = outFile.mkdir(outFolderName)
            outFolder_combined = outFile_combined.mkdir(outFolderName)

            for msc in self.mscSplits:

                for graph_type in graph_types:

                    graph_name = '_'.join([msc, year, graph_type])

                    input_graphs = []

                    # Create an empty TGraphAsymmErrors object to store the concatenated data
                    concatenated_graph = root.TGraphAsymmErrors()

                    for pt_bin in self.sorted_pt_bins:

                        combine_task_name_suffix = '-'.join([('PtTotal' if pt_bin.total_range else pt_bin.name), year]) # HACK
                        combine_task_name = '-'.join(['combineTask', self.tagger.name, self.wp.name]) + ('-'+combine_task_name_suffix if len(combine_task_name_suffix) else '')
                        combine_task_path = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel/combine', self.tagger.name, 'workdirs', combine_task_name)
                        infile_name = 'scale_factors_obs.root'
                        inFilePath = os.path.join(combine_task_path, infile_name)

                        inFile = root.TFile.Open(inFilePath, 'READ')
                        # inFile.SetDirectory(0)

                        input_graph = inFile.Get(graph_name)
                        input_graphs.append(input_graph)

                        # Get the number of points in the input graph
                        num_points = input_graph.GetN()
                        if num_points != 1:
                            sys.exit('input graph should only have one entry')

                        # Loop through each point in the input graph and add it to the concatenated graph
                        for i in range(num_points):
                            x = root.Double()
                            y = root.Double()
                            x_err_low = root.Double()
                            x_err_high = root.Double()
                            y_err_low = root.Double()
                            y_err_high = root.Double()

                            # Get the values and errors of the point
                            input_graph.GetPoint(i, x, y)
                            x_err_low = input_graph.GetErrorXlow(i)
                            x_err_high = input_graph.GetErrorXhigh(i)
                            y_err_low = input_graph.GetErrorYlow(i)
                            y_err_high = input_graph.GetErrorYhigh(i)

                            # setattr(x, 'y', v) is equivalent to x.y = v
                            setattr(self.scale_factors[msc][pt_bin.name][year], graph_type, ValueAndAsymmErrors(value=y, up=y_err_high, down=y_err_low))
                            # self.scale_factors[msc][pt_bin.name][year].tot.value = y
                            # self.scale_factors[msc][pt_bin.name][year].tot.up = y_err_high
                            # self.scale_factors[msc][pt_bin.name][year].tot.down = y_err_low

                            # Add the point with errors to the concatenated graph
                            concatenated_graph.SetPoint(concatenated_graph.GetN(), x, y)
                            concatenated_graph.SetPointError(concatenated_graph.GetN() - 1, x_err_low, x_err_high, y_err_low, y_err_high)

                        inFile.Close()

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

                    outFolder.cd()
                    concatenated_graph.Write(graph_name)
                    outFolder_combined.cd()
                    concatenated_graph.Write(graph_name)

            print('Wrote', outFilePath)

        print('Wrote', outFilePath_combined)

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
        y_high = (y_sf_range + y_gap) * len(self.mscSplits) / (1. - header_y_fraction)
        y_range = y_high - y_low

        y_sf_offset = OrderedDict()
        for index, msc in enumerate(reversed(self.mscSplits)):
            y_sf_offset[msc] = index * (y_sf_range + y_gap)

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

        # one graph for each year / each graph_type
        graphs = OrderedDict()
        points_per_graph = len(self.mscSplits) * len(self.sorted_pt_bins)
        for year in self.years:
            for graph_type in graph_types:
                arr_x = np.array([], dtype=np.float32)
                arr_y = np.array([], dtype=np.float32)
                arr_exl = np.array([], dtype=np.float32)
                arr_exh = np.array([], dtype=np.float32)
                arr_eyl = np.array([], dtype=np.float32)
                arr_eyh = np.array([], dtype=np.float32)
                for msc in self.mscSplits:
                    for pt_bin in self.sorted_pt_bins:
                        arr_x = np.append(arr_x, x_sf_position[year] + x_pt_offset[pt_bin.name])
                        arr_y = np.append(arr_y, getattr(self.scale_factors[msc][pt_bin.name][year], 'tot').value + y_sf_offset[msc]) # always use central value of tot
                        arr_exl = np.append(arr_exl, 0.5 * x_sf_width)
                        arr_exh = np.append(arr_exh, 0.5 * x_sf_width)
                        arr_eyl = np.append(arr_eyl, getattr(self.scale_factors[msc][pt_bin.name][year], graph_type).down)
                        arr_eyh = np.append(arr_eyh, getattr(self.scale_factors[msc][pt_bin.name][year], graph_type).up)
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

        mg = root.TMultiGraph()
        mg.SetTitle('')
        for year in self.years:
            graph_tot = graphs[year]['tot']
            graph_tot.SetFillColor(root.kRed)
            graph_tot.SetMarkerStyle(7)
            mg.Add(graph_tot)
            graph_uncorr = graphs[year]['uncorr']
            graph_uncorr.SetFillColor(root.kBlue)
            graph_uncorr.SetMarkerStyle(7)
            mg.Add(graph_uncorr)
            graph_stat = graphs[year]['stat']
            graph_stat.SetFillColor(root.kGreen)
            graph_stat.SetMarkerStyle(7)
            mg.Add(graph_stat)
        mg.Draw('a2p')
        mg.GetHistogram().GetXaxis().SetLimits(x_low, x_high)
        mg.GetHistogram().GetXaxis().SetTickLength(0.)
        mg.GetHistogram().GetXaxis().SetLabelOffset(-2.)
        mg.GetHistogram().SetMinimum(y_low)
        mg.GetHistogram().SetMaximum(y_high)
        mg.GetHistogram().GetYaxis().SetTickLength(0.)
        mg.GetHistogram().GetYaxis().SetLabelOffset(-2.)

        lines = {}
        for index, msc in enumerate(self.mscSplits):
            shifted_one = 1. + y_sf_offset[msc]
            lines[msc] = root.TLine(x_low, shifted_one, x_high, shifted_one)
            lines[msc].SetLineColor(root.kGray + 2)
            lines[msc].SetLineStyle(2)
            lines[msc].Draw()

        axes = {}
        ndiv = 505
        tick_length = 0.05
        w_low = y_sf_low
        w_high = y_sf_high
        for index_msc, msc in enumerate(self.mscSplits):
            x_axis_low = x_low
            y_axis_low = y_sf_low + y_sf_offset[msc]
            x_axis_high = x_low
            y_axis_high = y_sf_high + y_sf_offset[msc]
            axes.setdefault(msc, {})['left'] = root.TGaxis(x_axis_low, y_axis_low, x_axis_high, y_axis_high, w_low, w_high, ndiv, '-S')
            axes[msc]['left'].SetTickLength(tick_length)
            axes[msc]['left'].SetLabelSize(self.text_size)
            axes[msc]['left'].SetLabelFont(42)
            axes[msc]['left'].Draw()
            for index_pt_bin, pt_bin in enumerate(self.sorted_pt_bins):
                if index_pt_bin == 0:
                    continue
                x_axis_low = x_pt_low + x_pt_offset[pt_bin.name]
                x_axis_high = x_pt_low + x_pt_offset[pt_bin.name]
                axes.setdefault(msc, {})[pt_bin.name] = root.TGaxis(x_axis_low, y_axis_low, x_axis_high, y_axis_high, w_low, w_high, ndiv, 'U+-S')
                axes[msc][pt_bin.name].SetTickLength(tick_length)
                axes[msc][pt_bin.name].Draw()
            x_axis_low = x_high
            x_axis_high = x_high
            axes.setdefault(msc, {})['right'] = root.TGaxis(x_axis_low, y_axis_low, x_axis_high, y_axis_high, w_low, w_high, ndiv, 'U+S')
            axes[msc]['right'].SetTickLength(tick_length)
            axes[msc]['right'].Draw()

        deltaX = 0.03
        hw_corrected_deltaX = deltaX * hw_correction

        self.tlatex_cms = root.TLatex(self.coord.graph_to_pad_x(hw_corrected_deltaX), self.coord.graph_to_pad_y(1. - deltaX), 'CMS')
        self.tlatex_cms.SetTextAlign(13) # left top
        self.tlatex_cms.SetTextFont(62)
        self.tlatex_cms.SetTextSize(0.07)
        self.tlatex_cms.SetNDC()
        self.tlatex_cms.Draw()

        self.text_prelim = 'Work in Progress'
        self.tlatex_prelim = root.TLatex(self.coord.graph_to_pad_x(hw_corrected_deltaX + 0.14), self.coord.graph_to_pad_y(0.948), self.text_prelim)
        self.tlatex_prelim.SetTextAlign(13) # left top
        self.tlatex_prelim.SetTextFont(52)
        self.tlatex_prelim.SetTextSize(0.045)
        self.tlatex_prelim.SetNDC()
        self.tlatex_prelim.Draw()

        if tagger.name.startswith('ak8_'):
            probejetalgo = 'AK8'
        elif tagger.name.startswith('hotvr_'):
            probejetalgo = 'HOTVR'

        if self.wp.alt_name is not None:
            pt_text_string = self.tagger.label.replace('{WP_VALUE}', '{} ({})'.format(self.wp.get_cut_value(year), self.wp.alt_name))
        else:
            pt_text_string = self.tagger.label.replace('{WP_VALUE}', '{}'.format(self.wp.get_cut_value(year)))
        pt_text_string2 = '|#eta| < 2.5'
        pt_text_string += ', ' + pt_text_string2

        tlatex_pt = root.TLatex(self.coord.graph_to_pad_x(1. - hw_corrected_deltaX), self.coord.graph_to_pad_y(0.848), pt_text_string)
        tlatex_pt.SetTextAlign(31) # left top
        tlatex_pt.SetTextFont(42)
        tlatex_pt.SetTextSize(0.03)
        tlatex_pt.SetNDC()
        tlatex_pt.Draw()

        # tlatex_pt2 = root.TLatex(self.coord.graph_to_pad_x(1. - hw_corrected_deltaX), self.coord.graph_to_pad_y(0.79), pt_text_string2)
        # tlatex_pt2.SetTextAlign(31) # left top
        # tlatex_pt2.SetTextFont(42)
        # tlatex_pt2.SetTextSize(0.03)
        # tlatex_pt2.SetNDC()
        # tlatex_pt2.Draw()

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
            margin_l * 0.45,
            (1. - (y_gap / y_range + header_y_fraction)) / 2. * (1. - margin_t - margin_b) + margin_b,
            'Data-to-simulation scale factor'
        )
        y_label.SetTextAlign(21)
        y_label.SetTextAngle(90)
        y_label.SetTextFont(42)
        y_label.SetTextSize(self.text_size)
        y_label.SetNDC()
        y_label.Draw()

        msc_paves = {}
        msc_labels = {}
        for msc in self.mscSplits:
            y_center = 1. + y_sf_offset[msc]
            y_center += 0.03 # small correction for shadow of pave
            pave_distance_to_scale_factor_range = 0.1
            msc_paves[msc] = root.TPave(
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
            msc_paves[msc].SetLineWidth(0)
            msc_paves[msc].Draw()

            if msc == 'FullyMerged':
                msc_text = 'Fully merged top quark'
            elif msc == 'SemiMerged':
                msc_text = 'Semi-merged top quark'
            elif msc == 'NotMerged':
                msc_text = 'Not merged top quark'
            elif msc == 'WMerged':
                msc_text = 'Fully merged W boson'
            elif msc == 'NotTopOrWMerged':
                msc_text = 'Not merged'
            else:
                msc_text = 'not defined'

            y_center += 0.05

            msc_labels[msc] = root.TText(
                0.24 * hw_correction * x_range,
                ((y_center + 1. + y_gap * pave_distance_to_scale_factor_range + 0.06) / y_range) * y_range,
                msc_text
            )
            msc_labels[msc].SetTextAlign(21)
            msc_labels[msc].SetTextFont(52)
            msc_labels[msc].SetTextSize(0.03)
            msc_labels[msc].Draw()

            text_top_left = 'Run II Ultra Legacy'
            self.tlatex_top_left = root.TLatex(margin_l, 1. - margin_t + 0.01, text_top_left)
            self.tlatex_top_left.SetTextAlign(11) # left bottom
            self.tlatex_top_left.SetTextFont(42)
            self.tlatex_top_left.SetTextSize(self.text_size)
            self.tlatex_top_left.SetNDC()
            self.tlatex_top_left.Draw()

            text_top_right = ' + '.join([_YEARS.get(year).get('lumi_fb_display')+' fb^{#minus1}' for year in self.years])
            text_top_right += ' (13 TeV)'
            self.tlatex_top_right = root.TLatex(1. - margin_r, 1. - margin_t + 0.01, text_top_right)
            self.tlatex_top_right.SetTextAlign(31) # right bottom
            self.tlatex_top_right.SetTextFont(42)
            self.tlatex_top_right.SetTextSize(self.text_size)
            self.tlatex_top_right.SetNDC()
            self.tlatex_top_right.Draw()

        for fileType in ['.pdf', '.png']:
            plotFileName = self.sf_task_name+fileType
            plotFilePath = os.path.join(self.outputPath, plotFileName)
            self.canvas.SaveAs(plotFilePath)

if __name__ == '__main__':

    for tagger in taggers.values():

        for wp_index in range(len(tagger.get_wp())):

            sf = ScaleFactorCollector(taggerName=tagger.name, workingPoint=wp_index)
            sf.write_root_files()
            sf.plot_scale_factors()
            # print(sf.scale_factors)
