import ROOT as root
import os
import sys
import numpy as np
from copy import deepcopy

sys.path.append(os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/NicePlots/python'))
from plotter import CoordinateConverter

sys.path.append(os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/Analysis'))
from constants import _YEARS

def plots_comparison_EOY_and_UL(
    wp_index=0,
    probejetalgo='AK8',
    subjetbtag=False,
    mscSplitting='FullyMerged',
    year='18'
    ):

    # wp_index = 0
    # # wp_index = 1
    # # wp_index = 2
    # # wp_index = 3
    # # wp_index = 4
    # # probejetalgo = 'AK8'
    # probejetalgo = 'HOTVR'
    # subjetbtag = True
    # # subjetbtag = False
    if probejetalgo == 'HOTVR':
        subjetbtag = False
        wp_index = 0

    wp_name_map = {
        0: 'wp0p38_vt',
        1: 'wp0p47_t',
        2: 'wp0p52_m',
        3: 'wp0p61_l',
        4: 'wp0p69_vl',
    }
    wp_name_map_long = {
        0: 'very tight',
        1: 'tight',
        2: 'medium',
        3: 'loose',
        4: 'very loose',
    }
    wp_mis_map = {
        0: 'mis0p001',
        1: 'mis0p005',
        2: 'mis0p010',
        3: 'mis0p025',
        4: 'mis0p050',
    }

    # mscSplitting = 'FullyMerged'
    # # mscSplitting = 'SemiMerged'
    # # mscSplitting = 'NotMerged'

    # # year = '17'
    # year = '18'

    year_color = year_color = _YEARS['UL'+year]['tcolor']

    range_around_one = 0.4

    if mscSplitting == 'FullyMerged':
        mscSplitting_EOY = 'mergedTop'
        range_around_one = 0.4
    elif mscSplitting == 'SemiMerged':
        mscSplitting_EOY = 'semimerged'
        range_around_one = 0.4
    elif mscSplitting == 'NotMerged':
        mscSplitting_EOY = 'notmerged'
        # range_around_one = 1.0
        range_around_one = 0.4

    def modify_last_x_value(graph, new_x_value):
        if not graph or graph.GetN() == 0:
            return

        last_point_index = graph.GetN() - 1
        x = root.Double(0)
        y = root.Double(0)
        exl = graph.GetErrorXlow(last_point_index)
        exh = graph.GetErrorXhigh(last_point_index)
        eyl = graph.GetErrorYlow(last_point_index)
        eyh = graph.GetErrorYhigh(last_point_index)

        graph.GetPoint(last_point_index, x, y)  # Get the current coordinates of the last point

        graph.RemovePoint(last_point_index)  # Remove the last point
        graph.SetPoint(last_point_index, new_x_value, y)  # Set the new x value for the last point
        graph.SetPointError(last_point_index, exl, exh, eyl, eyh)  # Set back the errors


    def set_xy_errors_to_zero(graph, zero_x=True, zero_y=True):
        if not graph:
            return

        n_points = graph.GetN()
        for i in range(n_points):
            if zero_x:
                graph.SetPointEXlow(i, 0.0)    # Set negative x-error to zero
                graph.SetPointEXhigh(i, 0.0)   # Set positive x-error to zero
            if zero_y:
                graph.SetPointEYlow(i, 0.0)    # Set negative x-error to zero
                graph.SetPointEYhigh(i, 0.0)   # Set positive x-error to zero



    x_min = 150.
    x_max = 800.
    y_min = 1. - range_around_one
    y_max = 1. + range_around_one
    y_max_modified = (y_max - y_min) * 1.5 + y_min



    file_path_EOY_base = 'TopTaggingScaleFactors/preUL/scaleFactors'
    file_path_UL_base = 'TopTaggingScaleFactors/RunIISummer20UL/scaleFactors'

    file_name_EOY = '20'+year+'TopTaggingScaleFactors.root'
    file_name_UL = 'TopTaggingScaleFactors_RunIISummer20UL'+year+'_PUPPIv15.root'

    file_path_EOY = os.path.join(file_path_EOY_base, file_name_EOY)
    file_path_UL = os.path.join(file_path_UL_base, file_name_UL)

    file_EOY = root.TFile.Open(file_path_EOY, 'READ')
    file_UL = root.TFile.Open(file_path_UL, 'READ')

    if probejetalgo == 'HOTVR':
        hist_folder_EOY = 'HOTVR'
        graph_folder_UL = 'HOTVR_PUPPI_default'
    else:
        hist_folder_EOY = 'PUPPI_wp'+str(wp_index+1)+('_btag' if subjetbtag else '')
        graph_folder_UL = 'AK8_PUPPI_'
        if subjetbtag:
            graph_folder_UL += 'subjetbtag_'
        graph_folder_UL += wp_name_map[wp_index]
        if not subjetbtag:
            graph_folder_UL += '_'+wp_mis_map[wp_index]
    print(graph_folder_UL)

    hist_name_EOY_nominal = 'sf_'+mscSplitting_EOY+'_nominal'
    hist_name_EOY_up = 'sf_'+mscSplitting_EOY+'_up'
    hist_name_EOY_down = 'sf_'+mscSplitting_EOY+'_down'
    graph_name_UL = mscSplitting+'_tot'

    hist_path_EOY_nominal = os.path.join(hist_folder_EOY, hist_name_EOY_nominal)
    hist_path_EOY_up = os.path.join(hist_folder_EOY, hist_name_EOY_up)
    hist_path_EOY_down = os.path.join(hist_folder_EOY, hist_name_EOY_down)
    graph_path_UL = os.path.join(graph_folder_UL, graph_name_UL)

    hist_EOY_nominal = file_EOY.Get(hist_path_EOY_nominal)
    hist_EOY_up = file_EOY.Get(hist_path_EOY_up)
    hist_EOY_down = file_EOY.Get(hist_path_EOY_down)
    graph_UL = file_UL.Get(graph_path_UL)

    #set all UL x errors to zero:
    set_xy_errors_to_zero(graph_UL, zero_y=False)
    modify_last_x_value(graph_UL, 700.)

    # max_pt_EOY = hist_EOY_nominal.GetBinUpEdge(hist_EOY_nominal.GetNbinsX())
    # max_pt_UL = 1000. # should get this value from TGraphAsymmErrors directly, but who cares
    # max_pt = min(max_pt_EOY, max_pt_UL)

    values_EOY_x = []
    values_EOY_y = []
    values_EOY_exl = []
    values_EOY_exh = []
    values_EOY_eyl = []
    values_EOY_eyh = []

    for i_bin in range(1, hist_EOY_nominal.GetNbinsX()): # skip the last bin of the EOY scale factors: nonsense pt bin!

        is_last_bin = hist_EOY_nominal.GetBinCenter(i_bin) > 600. # kinda hacky

        if is_last_bin:
            values_EOY_x.append(700.)
            ex = 100.
            values_EOY_exl.append(ex)
            values_EOY_exh.append(ex)
        else:
            values_EOY_x.append(hist_EOY_nominal.GetBinCenter(i_bin))
            values_EOY_exl.append(hist_EOY_nominal.GetBinWidth(i_bin)*0.5)
            values_EOY_exh.append(hist_EOY_nominal.GetBinWidth(i_bin)*0.5)

        values_EOY_y.append(hist_EOY_nominal.GetBinContent(i_bin))
        values_EOY_eyl.append(hist_EOY_nominal.GetBinContent(i_bin) - hist_EOY_down.GetBinContent(i_bin))
        values_EOY_eyh.append(hist_EOY_up.GetBinContent(i_bin) - hist_EOY_nominal.GetBinContent(i_bin))

    values_EOY_x = np.array(values_EOY_x, dtype=np.float32)
    values_EOY_y = np.array(values_EOY_y, dtype=np.float32)
    values_EOY_exl = np.array(values_EOY_exl, dtype=np.float32)
    values_EOY_exh = np.array(values_EOY_exh, dtype=np.float32)
    values_EOY_eyl = np.array(values_EOY_eyl, dtype=np.float32)
    values_EOY_eyh = np.array(values_EOY_eyh, dtype=np.float32)

    graph_EOY = root.TGraphAsymmErrors(len(values_EOY_x), values_EOY_x, values_EOY_y, values_EOY_exl, values_EOY_exh, values_EOY_eyl, values_EOY_eyh)



    canvas = root.TCanvas('canvas', '', 600, 600)
    canvas.cd()


    graph_EOY_markers = deepcopy(graph_EOY)
    set_xy_errors_to_zero(graph_EOY_markers, zero_x=False)

    graph_EOY.SetTitle('')
    margin_l = 0.15
    margin_r = 0.05
    margin_b = 0.12
    margin_t = 0.08
    tick_length = 0.015 # fraction of canvas width/height
    canvas.SetMargin(margin_l, margin_r, margin_b, margin_t) #lrbt
    coord = CoordinateConverter(margin_l, margin_r, margin_b, margin_t)
    canvas.SetTickx(1);
    canvas.SetTicky(1);


    legend_sf_x1 = 0.04
    legend_sf_y1 = 0.57
    legend_sf_x2 = 0.5
    legend_sf_y2 = 0.8
    legend_sf = root.TLegend(coord.graph_to_pad_x(legend_sf_x1), coord.graph_to_pad_y(legend_sf_y1), coord.graph_to_pad_x(legend_sf_x2), coord.graph_to_pad_y(legend_sf_y2))
    legend_sf.SetTextSize(0.025)
    legend_sf.SetBorderSize(0)
    legend_sf.SetFillStyle(0)
    header_text = ''
    if mscSplitting == 'FullyMerged':
        header_text = 'Fully merged'
    elif mscSplitting == 'SemiMerged':
        header_text = 'Semi-merged'
    elif mscSplitting == 'NotMerged':
        header_text = 'Not merged'
    header_text += ' top quark'
    legend_sf.SetHeader('#bf{'+header_text+'}')

    graph_EOY.SetFillColorAlpha(root.kGray, 1.)
    graph_EOY.Draw('a2p')
    graph_EOY.GetHistogram().GetXaxis().SetLimits(x_min, x_max)
    # graph_EOY.GetHistogram().GetXaxis().SetNdivisions(1006)
    graph_EOY.GetHistogram().GetXaxis().SetTitle('Probe jet #it{p}_{T} [GeV]')
    graph_EOY.GetHistogram().GetYaxis().SetTitle('Data-to-simulation scale factor')
    graph_EOY.GetHistogram().GetXaxis().SetTitleOffset(1.3)
    graph_EOY.GetHistogram().GetXaxis().CenterTitle()
    graph_EOY.GetHistogram().GetYaxis().SetTitleOffset(1.5)
    graph_EOY.GetHistogram().GetYaxis().CenterTitle()
    graph_EOY.GetHistogram().SetMinimum(y_min)
    graph_EOY.GetHistogram().SetMaximum(y_max_modified)
    graph_EOY.GetHistogram().GetXaxis().ChangeLabel(7,-1,-1,-1,-1,-1,'#geq 800');

    line_one = root.TLine(x_min, 1., x_max, 1.)
    line_one.SetLineColor(root.kGray + 1)
    line_one.SetLineStyle(3)
    line_one.Draw()

    graph_EOY_markers.SetMarkerStyle(8)
    graph_EOY_markers.SetMarkerColor(root.kBlack)
    graph_EOY_markers.Draw('same pz')

    # graph_UL.SetMarkerStyle(25) # hollow square
    graph_UL.SetMarkerStyle(21) # solid square
    graph_UL.SetMarkerColor(year_color)
    graph_UL.SetLineColor(year_color)
    graph_UL.SetLineWidth(2)
    # graph_UL.Draw('same pz')



    initial_value_A = 1.
    initial_value_lambda = 0.

    # Limit the x-axis range of the fit
    fit_x_min = 200. if probejetalgo == 'HOTVR' else 300.  # Minimum x value for the fit range
    fit_x_max = 800. # Maximum x value for the fit range

    exp_fit_func = root.TF1("exp_fit_func", "[0] * TMath::Exp([1] * x)")
    exp_fit_func.SetParameter(0, initial_value_A)      # Initial value for A
    exp_fit_func.SetParameter(1, initial_value_lambda) # Initial value for lambda
    exp_fit_func.SetRange(fit_x_min, fit_x_max)
    # Perform the fit
    graph_UL_for_fit = deepcopy(graph_UL) # make a deepcopy that is not drawn. If we perform the following fit directly on the graph that is drawn, then the plot has a red line that i cannot remove
    fit_result = graph_UL_for_fit.Fit(exp_fit_func, "S")  # "S" indicates to store the fit result


    fit_amplitude = exp_fit_func.GetParameter(0)
    fit_amplitude_error = exp_fit_func.GetParError(0)

    fit_lambda = exp_fit_func.GetParameter(1)
    fit_lambda_error = exp_fit_func.GetParError(1)

    # # Create a TGraphErrors to store the error band
    # error_band_graph = root.TGraphErrors(graph_UL.GetN())
    #
    # # Generate the error band by varying fit parameters
    # for i in range(graph_UL.GetN()):
    #     x = graph_UL.GetX()[i]
    #     y_fit = exp_fit_func.Eval(x)
    #     y_err_fit = fit_result.Chi2() / (graph_UL.GetN() - 2) * (1 + (x - fit_lambda) ** 2) * fit_lambda_error ** 2  # Error from fit parameters
    #     error_band_graph.SetPoint(i, x, y_fit)
    #     error_band_graph.SetPointError(i, 0, y_err_fit)


    # Get the covariance matrix from the fit result
    cov_matrix = fit_result.GetCovarianceMatrix()

    # # Create a TGraphAsymmErrors to store the error band
    # error_band_graph = root.TGraphAsymmErrors(graph_UL.GetN())
    #
    # for i in range(graph_UL.GetN()):
    #     x = graph_UL.GetX()[i]
    #     y_fit = exp_fit_func.Eval(x)
    #
    #     # Calculate the error band using the covariance matrix
    #     err_fit = root.TMath.Sqrt(cov_matrix(1,1) * (x**2) + cov_matrix(0,0) + 2 * x * cov_matrix(0,1))
    #
    #     # Fill the TGraphAsymmErrors
    #     error_band_graph.SetPoint(i, x, y_fit)
    #     error_band_graph.SetPointError(i, 0, 0, err_fit, err_fit)


    n_points_fit = 100
    error_band_graph = root.TGraphAsymmErrors(n_points_fit)
    for i_point in range(n_points_fit+1):
        x_value = fit_x_min + (fit_x_max - fit_x_min) / n_points_fit * i_point
        y_fit = exp_fit_func.Eval(x_value)
        # Calculate the error band using the covariance matrix
        err_fit = root.TMath.Sqrt(cov_matrix(1,1) * (x_value**2) + cov_matrix(0,0) + 2 * x_value * cov_matrix(0,1))
        # Fill the TGraphAsymmErrors
        error_band_graph.SetPoint(i_point, x_value, y_fit)
        error_band_graph.SetPointError(i_point, 0, 0, err_fit, err_fit)


    error_band_graph.SetFillColorAlpha(year_color, 0.4)
    # error_band_graph.SetFillStyle(3356)
    error_band_graph.SetLineColor(year_color)
    error_band_graph.SetLineWidth(2)
    error_band_graph.Draw("3 same")
    exp_fit_func.Draw("same")
    exp_fit_func.SetLineColor(year_color)
    exp_fit_func.SetLineWidth(2)
    graph_UL.Draw('same pz')

    # Access the chi-squared and number of degrees of freedom
    chi2 = fit_result.Chi2()
    ndf = fit_result.Ndf()
    print('chi2', chi2)
    print('ndf', ndf)



    graph_EOY_markers_for_legend = deepcopy(graph_EOY_markers)
    graph_EOY_markers_for_legend.SetFillColorAlpha(root.kGray, 1.)

    legend_sf.AddEntry(graph_UL, 'UL 20'+year+' (this study)', 'ep')
    legend_sf.AddEntry(error_band_graph, 'Exponential fit to UL 20'+year+' (#chi^{2}/n.d.f. = '+"{:.2f}".format(chi2)+'/'+str(ndf)+')', 'lf')
    legend_sf.AddEntry(graph_EOY_markers_for_legend, 'EOY 20'+year+' (CMS DP 2020/025)', 'lpf')






    text_size = 0.035

    tlatex_probejetalgo = root.TLatex(coord.graph_to_pad_x(0.95), coord.graph_to_pad_y(0.95), probejetalgo+' PUPPI')
    tlatex_probejetalgo.SetTextAlign(33) # right top
    tlatex_probejetalgo.SetTextFont(62)
    tlatex_probejetalgo.SetTextSize(0.05)
    tlatex_probejetalgo.SetTextColor(root.kGray+1)
    tlatex_probejetalgo.SetNDC()
    tlatex_probejetalgo.Draw()

    pt_text_string = ''
    if probejetalgo == 'AK8':
        pt_text_string += '#tau_{3}/#tau_{2} ('+wp_name_map_long[wp_index]+' working point)'
        if subjetbtag:
            pt_text_string += ' + subjet b tag'
    elif probejetalgo == 'HOTVR':
        pt_text_string += '#tau_{3}/#tau_{2} < 0.56 + HOTVR cuts'

    if pt_text_string:
        tlatex_pt = root.TLatex(coord.graph_to_pad_x(0.95), coord.graph_to_pad_y(0.84), pt_text_string)
        tlatex_pt.SetTextAlign(31) # left top
        tlatex_pt.SetTextFont(42)
        tlatex_pt.SetTextSize(0.025)
        tlatex_pt.SetNDC()
        tlatex_pt.Draw()

    tlatex_cms = root.TLatex(coord.graph_to_pad_x(0.05), coord.graph_to_pad_y(0.95), 'CMS')
    tlatex_cms.SetTextAlign(13) # left top
    tlatex_cms.SetTextFont(62)
    tlatex_cms.SetTextSize(0.05)
    tlatex_cms.SetNDC()
    tlatex_cms.Draw()

    text_prelim = 'Private Work'

    tlatex_prelim = root.TLatex(coord.graph_to_pad_x(0.05), coord.graph_to_pad_y(0.87), text_prelim)
    tlatex_prelim.SetTextAlign(13) # left top
    tlatex_prelim.SetTextFont(52)
    tlatex_prelim.SetTextSize(text_size)
    tlatex_prelim.SetNDC()
    tlatex_prelim.Draw()

    # text_top_left = '20'+year+' data-taking era'
    text_top_left = ''
    text_top_right = '#sqrt{#it{s}} = 13 TeV'

    if text_top_left:
        tlatex_top_left = root.TLatex(margin_l, 1. - margin_t + 0.01, text_top_left)
        tlatex_top_left.SetTextAlign(11) # left bottom
        tlatex_top_left.SetTextFont(42)
        tlatex_top_left.SetTextSize(text_size)
        tlatex_top_left.SetNDC()
        tlatex_top_left.Draw()

    if text_top_right:
        tlatex_top_right = root.TLatex(1. - margin_r, 1. - margin_t + 0.01, text_top_right)
        tlatex_top_right.SetTextAlign(31) # right bottom
        tlatex_top_right.SetTextFont(42)
        tlatex_top_right.SetTextSize(text_size)
        tlatex_top_right.SetNDC()
        tlatex_top_right.Draw()


    legend_sf.Draw()


    canvas.Update()
    tickScaleX = (canvas.GetUxmax() - canvas.GetUxmin()) / (canvas.GetX2() - canvas.GetX1()) * (canvas.GetWh() * canvas.GetAbsHNDC())
    tickScaleY = (canvas.GetUymax() - canvas.GetUymin()) / (canvas.GetY2() - canvas.GetY1()) * (canvas.GetWw() * canvas.GetAbsWNDC())
    graph_EOY.GetXaxis().SetTickLength(canvas.GetWh() * tick_length / tickScaleX)
    graph_EOY.GetYaxis().SetTickLength(canvas.GetWw() * tick_length / tickScaleY)

    root.gPad.Update()
    root.gPad.RedrawAxis()


    output_path = 'plots_comparison_EOY_and_UL'
    os.system('mkdir -p '+output_path)

    # output_file_name = 'test.pdf'
    output_file_name = 'ScaleFactorComparisonEOYandUL_20{}_{}_wp{}_{}'.format(year, probejetalgo, str(wp_index), mscSplitting)
    if subjetbtag:
        output_file_name += '_subjetbtag'
    output_file_name += '.pdf'
    output_file_path = os.path.join(output_path, output_file_name)

    canvas.SaveAs(output_file_path)


if __name__=='__main__':

    for year in ['17', '18']:
        for mscSplitting in ['FullyMerged', 'SemiMerged', 'NotMerged']:

            plots_comparison_EOY_and_UL(
                wp_index=0,
                probejetalgo='HOTVR',
                subjetbtag=False,
                mscSplitting=mscSplitting,
                year=year
            )

            for wp_index in range(5):
                for subjetbtag in [True, False]:
                    plots_comparison_EOY_and_UL(
                        wp_index=wp_index,
                        probejetalgo='AK8',
                        subjetbtag=subjetbtag,
                        mscSplitting=mscSplitting,
                        year=year
                    )
