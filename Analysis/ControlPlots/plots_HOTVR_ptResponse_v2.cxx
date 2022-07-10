#include "../constants.h"

using namespace macros;


double binc(TH1F* hist, const int i_bin) {
  return hist->GetBinContent(i_bin);
}


void do_plot(const Year & year) {
  const string inputDir = (string)getenv("CMSSW_BASE")+"/src/UHH2/LegacyTopTagging/output/WorkingPointStudy/"+kYears.at(year).name+"/nominal/";
  const string inputFileName = "HOTVR_JetPtResponse.root";
  const string inputFilePath = inputDir + inputFileName;
  TFile *inputFile = TFile::Open(inputFilePath.c_str(), "READ");

  vector<string> inputHistNames = {
    "response_raw_mean",
    "response_raw_mean_mc_unc",
    "response_raw_std",
    "response_raw_std_mc_unc",
    "response_corr_mean",
    "response_corr_mean_mc_unc",
    "response_corr_mean_jes_up",
    "response_corr_mean_jes_down",
    "response_corr_mean_jer_up",
    "response_corr_mean_jer_down",
    "response_corr_std",
    "response_corr_std_mc_unc",
    "response_corr_std_jes_up",
    "response_corr_std_jes_down",
    "response_corr_std_jer_up",
    "response_corr_std_jer_down",
  };

  map<string, TH1F*> inputHists;
  for (const string & inputHistName : inputHistNames) {
    inputHists[inputHistName] = (TH1F*)inputFile->Get(inputHistName.c_str());
  }

  TH1F *helperHist = inputHists["response_raw_mean"];
  int n_bins = helperHist->GetNbinsX();

  vector<double> x_center;
  vector<double> x_low;
  vector<double> x_high;

  vector<double> y_raw_mean;
  vector<double> y_raw_mean_mc_unc;

  vector<double> y_raw_std_plus;
  vector<double> y_raw_std_minus;
  vector<double> y_raw_std_mc_unc;

  vector<double> y_raw_std_norm;
  vector<double> y_raw_std_norm_mc_unc;

  vector<double> y_corr_mean;
  vector<double> y_corr_mean_mc_unc;
  vector<double> y_corr_mean_mc_jer_plus;
  vector<double> y_corr_mean_mc_jer_minus;
  vector<double> y_corr_mean_mc_jer_jes_plus;
  vector<double> y_corr_mean_mc_jer_jes_minus;

  vector<double> y_corr_std_plus;
  vector<double> y_corr_std_minus;
  vector<double> y_corr_std_mc_unc;
  vector<double> y_corr_std_mc_jer_plus;
  vector<double> y_corr_std_mc_jer_minus;
  vector<double> y_corr_std_mc_jer_jes_plus;
  vector<double> y_corr_std_mc_jer_jes_minus;

  vector<double> y_corr_std_norm;
  vector<double> y_corr_std_norm_mc_unc;
  vector<double> y_corr_std_norm_mc_jer_plus;
  vector<double> y_corr_std_norm_mc_jer_minus;
  vector<double> y_corr_std_norm_mc_jer_jes_plus;
  vector<double> y_corr_std_norm_mc_jer_jes_minus;

  for (int i_bin = 1; i_bin <= n_bins; i_bin++) {
    const int i_entry = i_bin - 1;

    x_center.push_back(helperHist->GetBinCenter(i_bin));
    x_low.push_back(helperHist->GetBinWidth(i_bin) / 2.);
    x_high.push_back(helperHist->GetBinWidth(i_bin) / 2.);

    double mean = binc(inputHists["response_raw_mean"], i_bin);
    y_raw_mean.push_back(mean);
    double mc_unc = binc(inputHists["response_raw_mean_mc_unc"], i_bin);
    y_raw_mean_mc_unc.push_back(mc_unc);

    double std = binc(inputHists["response_raw_std"], i_bin);
    y_raw_std_plus.push_back(mean + std);
    y_raw_std_minus.push_back(mean - std);
    y_raw_std_norm.push_back(std / mean);
    mc_unc = binc(inputHists["response_raw_std_mc_unc"], i_bin);
    y_raw_std_mc_unc.push_back(mc_unc);
    // y_raw_std_norm_mc_unc.push_back(mc_unc / mean);
    y_raw_std_norm_mc_unc.push_back(sqrt( pow(y_raw_std_mc_unc.at(i_entry), 2) + pow(std/mean, 2) * pow(y_raw_mean_mc_unc.at(i_entry), 2) ) / mean);

    mean = binc(inputHists["response_corr_mean"], i_bin);
    y_corr_mean.push_back(mean);
    mc_unc = binc(inputHists["response_corr_mean_mc_unc"], i_bin);
    double jes_up = binc(inputHists["response_corr_mean_jes_up"], i_bin) - mean;
    double jes_down = binc(inputHists["response_corr_mean_jes_down"], i_bin) - mean;
    double jes_plus = max(0., max(jes_up, jes_down));
    double jes_minus = max(0., -min(jes_up, jes_down));
    double jer_up = binc(inputHists["response_corr_mean_jer_up"], i_bin) - mean;
    double jer_down = binc(inputHists["response_corr_mean_jer_down"], i_bin) - mean;
    double jer_plus = max(0., max(jer_up, jer_down));
    double jer_minus = max(0., -min(jer_up, jer_down));
    double mc_jer_plus = sqrt(jer_plus*jer_plus);
    double mc_jer_minus = sqrt(jer_minus*jer_minus);
    double mc_jer_jes_plus = sqrt(mc_unc*mc_unc + jer_plus*jer_plus + jes_plus*jes_plus);
    double mc_jer_jes_minus = sqrt(mc_unc*mc_unc + jer_minus*jer_minus + jes_minus*jes_minus);
    y_corr_mean_mc_unc.push_back(mc_unc);
    y_corr_mean_mc_jer_plus.push_back(mc_jer_plus);
    y_corr_mean_mc_jer_minus.push_back(mc_jer_minus);
    y_corr_mean_mc_jer_jes_plus.push_back(mc_jer_jes_plus);
    y_corr_mean_mc_jer_jes_minus.push_back(mc_jer_jes_minus);

    std = binc(inputHists["response_corr_std"], i_bin);
    y_corr_std_plus.push_back(mean + std);
    y_corr_std_minus.push_back(mean - std);
    y_corr_std_norm.push_back(std / mean);
    mc_unc = binc(inputHists["response_corr_std_mc_unc"], i_bin);
    jes_up = binc(inputHists["response_corr_std_jes_up"], i_bin) - std;
    jes_down = binc(inputHists["response_corr_std_jes_down"], i_bin) - std;
    jes_plus = max(0., max(jes_up, jes_down));
    jes_minus = max(0., -min(jes_up, jes_down));
    jer_up = binc(inputHists["response_corr_std_jer_up"], i_bin) - std;
    jer_down = binc(inputHists["response_corr_std_jer_down"], i_bin) - std;
    jer_plus = max(0., max(jer_up, jer_down));
    jer_minus = max(0., -min(jer_up, jer_down));
    mc_jer_plus = sqrt(jer_plus*jer_plus);
    mc_jer_minus = sqrt(jer_minus*jer_minus);
    mc_jer_jes_plus = sqrt(mc_unc*mc_unc + jer_plus*jer_plus + jes_plus*jes_plus);
    mc_jer_jes_minus = sqrt(mc_unc*mc_unc + jer_minus*jer_minus + jes_minus*jes_minus);
    y_corr_std_mc_unc.push_back(mc_unc);
    // y_corr_std_norm_mc_unc.push_back(mc_unc / mean);
    y_corr_std_norm_mc_unc.push_back(sqrt( pow(y_corr_std_mc_unc.at(i_entry), 2) + pow(std/mean, 2) * pow(y_corr_mean_mc_unc.at(i_entry), 2) ) / mean);
    y_corr_std_mc_jer_plus.push_back(mc_jer_plus);
    // y_corr_std_norm_mc_jer_plus.push_back(mc_jer_plus / mean);
    y_corr_std_norm_mc_jer_plus.push_back(sqrt( pow(y_corr_std_mc_jer_plus.at(i_entry), 2) + pow(std/mean, 2) * pow(y_corr_mean_mc_jer_plus.at(i_entry), 2) ) / mean);
    y_corr_std_mc_jer_minus.push_back(mc_jer_minus);
    // y_corr_std_norm_mc_jer_minus.push_back(mc_jer_minus / mean);
    y_corr_std_norm_mc_jer_minus.push_back(sqrt( pow(y_corr_std_mc_jer_minus.at(i_entry), 2) + pow(std/mean, 2) * pow(y_corr_mean_mc_jer_minus.at(i_entry), 2) ) / mean);
    y_corr_std_mc_jer_jes_plus.push_back(mc_jer_jes_plus);
    // y_corr_std_norm_mc_jer_jes_plus.push_back(mc_jer_jes_plus / mean);
    y_corr_std_norm_mc_jer_jes_plus.push_back(sqrt( pow(y_corr_std_mc_jer_jes_plus.at(i_entry), 2) + pow(std/mean, 2) * pow(y_corr_mean_mc_jer_jes_plus.at(i_entry), 2) ) / mean);
    y_corr_std_mc_jer_jes_minus.push_back(mc_jer_jes_minus);
    // y_corr_std_norm_mc_jer_jes_minus.push_back(mc_jer_jes_minus / mean);
    y_corr_std_norm_mc_jer_jes_minus.push_back(sqrt( pow(y_corr_std_mc_jer_jes_minus.at(i_entry), 2) + pow(std/mean, 2) * pow(y_corr_mean_mc_jer_jes_minus.at(i_entry), 2) ) / mean);
  }

  TGraphAsymmErrors *g_raw_mean_unc_mc = new TGraphAsymmErrors(n_bins, &x_center[0], &y_raw_mean[0], &x_low[0], &x_high[0], &y_raw_mean_mc_unc[0], &y_raw_mean_mc_unc[0]);
  TGraphAsymmErrors *g_raw_mean_line = new TGraphAsymmErrors(*g_raw_mean_unc_mc);

  TGraphAsymmErrors *g_raw_std_plus_unc_mc = new TGraphAsymmErrors(n_bins, &x_center[0], &y_raw_std_plus[0], &x_low[0], &x_high[0], &y_raw_std_mc_unc[0], &y_raw_std_mc_unc[0]);
  TGraphAsymmErrors *g_raw_std_plus_line = new TGraphAsymmErrors(*g_raw_std_plus_unc_mc);

  TGraphAsymmErrors *g_raw_std_minus_unc_mc = new TGraphAsymmErrors(n_bins, &x_center[0], &y_raw_std_minus[0], &x_low[0], &x_high[0], &y_raw_std_mc_unc[0], &y_raw_std_mc_unc[0]);
  TGraphAsymmErrors *g_raw_std_minus_line = new TGraphAsymmErrors(*g_raw_std_minus_unc_mc);

  TGraphAsymmErrors *g_raw_std_norm_unc_mc = new TGraphAsymmErrors(n_bins, &x_center[0], &y_raw_std_norm[0], &x_low[0], &x_high[0], &y_raw_std_norm_mc_unc[0], &y_raw_std_norm_mc_unc[0]);
  TGraphAsymmErrors *g_raw_std_norm_line = new TGraphAsymmErrors(*g_raw_std_norm_unc_mc);

  TGraphAsymmErrors *g_corr_mean_unc_mc = new TGraphAsymmErrors(n_bins, &x_center[0], &y_corr_mean[0], &x_low[0], &x_high[0], &y_corr_mean_mc_unc[0], &y_corr_mean_mc_unc[0]);
  TGraphAsymmErrors *g_corr_mean_line = new TGraphAsymmErrors(*g_corr_mean_unc_mc);
  TGraphAsymmErrors *g_corr_mean_unc_mc_jer = new TGraphAsymmErrors(n_bins, &x_center[0], &y_corr_mean[0], &x_low[0], &x_high[0], &y_corr_mean_mc_jer_minus[0], &y_corr_mean_mc_jer_plus[0]);
  TGraphAsymmErrors *g_corr_mean_unc_mc_jer_jes = new TGraphAsymmErrors(n_bins, &x_center[0], &y_corr_mean[0], &x_low[0], &x_high[0], &y_corr_mean_mc_jer_jes_minus[0], &y_corr_mean_mc_jer_jes_plus[0]);

  TGraphAsymmErrors *g_corr_std_plus_unc_mc = new TGraphAsymmErrors(n_bins, &x_center[0], &y_corr_std_plus[0], &x_low[0], &x_high[0], &y_corr_std_mc_unc[0], &y_corr_std_mc_unc[0]);
  TGraphAsymmErrors *g_corr_std_plus_line = new TGraphAsymmErrors(*g_corr_std_plus_unc_mc);
  TGraphAsymmErrors *g_corr_std_plus_unc_mc_jer = new TGraphAsymmErrors(n_bins, &x_center[0], &y_corr_std_plus[0], &x_low[0], &x_high[0], &y_corr_std_mc_jer_minus[0], &y_corr_std_mc_jer_plus[0]);
  TGraphAsymmErrors *g_corr_std_plus_unc_mc_jer_jes = new TGraphAsymmErrors(n_bins, &x_center[0], &y_corr_std_plus[0], &x_low[0], &x_high[0], &y_corr_std_mc_jer_jes_minus[0], &y_corr_std_mc_jer_jes_plus[0]);

  TGraphAsymmErrors *g_corr_std_minus_unc_mc = new TGraphAsymmErrors(n_bins, &x_center[0], &y_corr_std_minus[0], &x_low[0], &x_high[0], &y_corr_std_mc_unc[0], &y_corr_std_mc_unc[0]);
  TGraphAsymmErrors *g_corr_std_minus_line = new TGraphAsymmErrors(*g_corr_std_minus_unc_mc);
  TGraphAsymmErrors *g_corr_std_minus_unc_mc_jer = new TGraphAsymmErrors(n_bins, &x_center[0], &y_corr_std_minus[0], &x_low[0], &x_high[0], &y_corr_std_mc_jer_minus[0], &y_corr_std_mc_jer_plus[0]);
  TGraphAsymmErrors *g_corr_std_minus_unc_mc_jer_jes = new TGraphAsymmErrors(n_bins, &x_center[0], &y_corr_std_minus[0], &x_low[0], &x_high[0], &y_corr_std_mc_jer_jes_minus[0], &y_corr_std_mc_jer_jes_plus[0]);

  TGraphAsymmErrors *g_corr_std_norm_unc_mc = new TGraphAsymmErrors(n_bins, &x_center[0], &y_corr_std_norm[0], &x_low[0], &x_high[0], &y_corr_std_norm_mc_unc[0], &y_corr_std_norm_mc_unc[0]);
  TGraphAsymmErrors *g_corr_std_norm_line = new TGraphAsymmErrors(*g_corr_std_norm_unc_mc);
  TGraphAsymmErrors *g_corr_std_norm_unc_mc_jer = new TGraphAsymmErrors(n_bins, &x_center[0], &y_corr_std_norm[0], &x_low[0], &x_high[0], &y_corr_std_norm_mc_jer_minus[0], &y_corr_std_norm_mc_jer_plus[0]);
  TGraphAsymmErrors *g_corr_std_norm_unc_mc_jer_jes = new TGraphAsymmErrors(n_bins, &x_center[0], &y_corr_std_norm[0], &x_low[0], &x_high[0], &y_corr_std_norm_mc_jer_jes_minus[0], &y_corr_std_norm_mc_jer_jes_plus[0]);

  Int_t color_raw = kYears.at(year).tcolor;
  Int_t color_corr = kYears.at(year).tcolor;
  const double transparency1 = 1.0;
  const double transparency2 = 0.5;

  // MAIN

  // g_raw_mean_unc_mc->SetFillColorAlpha(color_raw, transparency1);
  g_raw_mean_unc_mc->SetFillColor(kGray+2);
  g_raw_mean_unc_mc->SetFillStyle(3356);
  g_raw_mean_unc_mc->SetLineWidth(0);

  g_raw_mean_line->SetLineColor(kBlack);
  g_raw_mean_line->SetLineStyle(kDotted);
  g_raw_mean_line->SetFillColorAlpha(0, 0);
  g_raw_mean_line->SetLineWidth(1);


  // g_corr_mean_unc_mc_jer_jes->SetFillColorAlpha(color_corr, 1);
  g_corr_mean_unc_mc_jer_jes->SetFillColorAlpha(color_corr, transparency2);
  g_corr_mean_unc_mc_jer_jes->SetLineWidth(0);

  // g_corr_mean_unc_mc_jer->SetFillColorAlpha(kWhite, 1);
  g_corr_mean_unc_mc_jer->SetFillColorAlpha(color_corr, transparency1);
  g_corr_mean_unc_mc_jer->SetLineWidth(0);

  // g_corr_mean_unc_mc->SetFillColorAlpha(color_corr, transparency1);
  g_corr_mean_unc_mc->SetFillColor(kGray+2);
  g_corr_mean_unc_mc->SetFillStyle(3356);
  g_corr_mean_unc_mc->SetLineWidth(0);

  g_corr_mean_line->SetLineColor(kBlack);
  g_corr_mean_line->SetLineStyle(kSolid);
  g_corr_mean_line->SetFillColorAlpha(0, 0);
  g_corr_mean_line->SetLineWidth(1);


  // RATIO

  // g_raw_std_norm_unc_mc->SetFillColorAlpha(color_raw, transparency1);
  g_raw_std_norm_unc_mc->SetFillColor(kGray+2);
  g_raw_std_norm_unc_mc->SetFillStyle(3356);
  g_raw_std_norm_unc_mc->SetLineWidth(0);

  g_raw_std_norm_line->SetLineColor(kBlack);
  g_raw_std_norm_line->SetLineStyle(kDotted);
  g_raw_std_norm_line->SetFillColorAlpha(0, 0);
  g_raw_std_norm_line->SetLineWidth(1);


  // g_corr_std_norm_unc_mc_jer_jes->SetFillColorAlpha(color_corr, 1);
  g_corr_std_norm_unc_mc_jer_jes->SetFillColorAlpha(color_corr, transparency2);
  g_corr_std_norm_unc_mc_jer_jes->SetLineWidth(0);

  // g_corr_std_norm_unc_mc_jer->SetFillColorAlpha(kWhite, 1);
  g_corr_std_norm_unc_mc_jer->SetFillColorAlpha(color_corr, transparency1);
  g_corr_std_norm_unc_mc_jer->SetLineWidth(0);

  // g_corr_std_norm_unc_mc->SetFillColorAlpha(color_corr, transparency1);
  g_corr_std_norm_unc_mc->SetFillColor(kGray+2);
  g_corr_std_norm_unc_mc->SetFillStyle(3356);
  g_corr_std_norm_unc_mc->SetLineWidth(0);

  g_corr_std_norm_line->SetLineColor(kBlack);
  g_corr_std_norm_line->SetLineStyle(kSolid);
  g_corr_std_norm_line->SetFillColorAlpha(0, 0);
  g_corr_std_norm_line->SetLineWidth(1);




  const float canvas_height = 600;
  const float canvas_width = 600;

  const float canvas_margin_l = 0.15;
  const float canvas_margin_r = 0.05;
  const float canvas_margin_b = 0.12 * canvas_width / canvas_height;
  const float canvas_margin_t = 0.08 * canvas_width / canvas_height;

  auto coord = new CoordinateConverter();
  coord->init(canvas_margin_l, canvas_margin_r, canvas_margin_b, canvas_margin_t);

  const float tick_length = 0.015; // fraction of canvas width/height
  const float border_y = 1. - kGoldenRatio; // access here, where the canvas is split between main plot and ratio plot
  const float border_margin = 0.015 * canvas_width / canvas_height;

  TCanvas *canvas = new TCanvas("canvas", "canvas title", canvas_width, canvas_height);

  canvas->cd();
  TPad *p_main = new TPad("pad_main", "pad title", 0, border_y, 1, 1);
  const float p_main_relative_height = 1. - border_y;
  p_main->SetTopMargin(canvas_margin_t / p_main_relative_height);
  p_main->SetBottomMargin(border_margin / p_main_relative_height);
  p_main->SetLeftMargin(canvas_margin_l);
  p_main->SetRightMargin(canvas_margin_r);
  // p_main->SetMargin();
  // p_main->SetFrameFillColor(kGreen);
  p_main->SetTickx(1);
  p_main->SetTicky(1);
  // p_main->SetLogx();
  // p_main->SetGrid();
  p_main->Draw();

  canvas->cd();
  TPad *p_ratio = new TPad("pad_ratio", "pad title", 0, 0, 1, border_y);
  const float p_ratio_relative_height = border_y;
  p_ratio->SetTopMargin(border_margin / p_ratio_relative_height);
  p_ratio->SetBottomMargin(canvas_margin_b / p_ratio_relative_height);
  p_ratio->SetLeftMargin(canvas_margin_l);
  p_ratio->SetRightMargin(canvas_margin_r);
  // p_ratio->SetMargin();
  // p_ratio->SetFrameFillColor(kOrange);
  p_ratio->SetTickx(1);
  p_ratio->SetTicky(1);
  // p_ratio->SetLogx();
  // p_ratio->SetGrid();
  p_ratio->Draw();

  p_main->cd();

  const float text_size = 0.035;
  const float small_text_size = 0.018;


  TMultiGraph *mg = new TMultiGraph("mg", "mg title");

  mg->Add(g_raw_mean_unc_mc);
  mg->Add(g_raw_mean_line);
  // mg->Add(g_raw_std_plus_unc_mc);
  // mg->Add(g_raw_std_plus_line);
  // mg->Add(g_raw_std_minus_unc_mc);
  // mg->Add(g_raw_std_minus_line);
  mg->Add(g_corr_mean_unc_mc_jer_jes);
  mg->Add(g_corr_mean_unc_mc_jer);
  mg->Add(g_corr_mean_unc_mc);
  mg->Add(g_corr_mean_line);
  // mg->Add(g_corr_std_plus_unc_mc);
  // mg->Add(g_corr_std_plus_line);
  // mg->Add(g_corr_std_minus_unc_mc);
  // mg->Add(g_corr_std_minus_line);

  mg->Draw("AL3");
  mg->SetTitle("");
  mg->GetHistogram()->GetXaxis()->SetTitle("#it{p}_{T}^{gen} [GeV]");
  // mg->GetHistogram()->GetYaxis()->SetTitle("#LT #it{p}_{T}^{rec} / #it{p}_{T}^{gen} #GT");
  mg->GetHistogram()->GetYaxis()->SetTitle("#LT #it{R} #GT, #it{R} = #it{p}_{T}^{rec} / #it{p}_{T}^{gen}");
  // mg->GetHistogram()->GetXaxis()->CenterTitle();
  mg->GetHistogram()->GetXaxis()->SetTitleOffset(100);
  mg->GetHistogram()->GetXaxis()->SetLabelOffset(100);
  // mg->GetHistogram()->GetYaxis()->CenterTitle();

  mg->GetHistogram()->GetXaxis()->SetLimits(150., 1250.);

  mg->GetHistogram()->GetYaxis()->SetTitleSize(text_size / p_main_relative_height); // restore original font size
  mg->GetHistogram()->GetYaxis()->SetTitleOffset(1.3);
  mg->GetHistogram()->GetYaxis()->SetLabelSize(text_size / p_main_relative_height); // restore original font size
  mg->GetHistogram()->GetYaxis()->SetLabelOffset(0.0077 / p_main_relative_height);

  double maximum_saved = mg->GetHistogram()->GetMaximum();
  double minimum_saved = mg->GetHistogram()->GetMinimum();
  mg->GetHistogram()->SetMaximum(1.08);
  mg->GetHistogram()->SetMinimum(0.88);

  // https://root-forum.cern.ch/t/inconsistent-tick-length/18563/8
  p_main->Update();
  double tickScaleX = (p_main->GetUxmax() - p_main->GetUxmin()) / (p_main->GetX2() - p_main->GetX1()) * (p_main->GetWh() * p_main->GetAbsHNDC());
  double tickScaleY = (p_main->GetUymax() - p_main->GetUymin()) / (p_main->GetY2() - p_main->GetY1()) * (p_main->GetWw() * p_main->GetAbsWNDC());
  mg->GetHistogram()->GetYaxis()->SetTickLength(p_main->GetWw() * tick_length / tickScaleY);
  mg->GetHistogram()->GetXaxis()->SetTickLength(p_main->GetWh() * tick_length / tickScaleX);

  gPad->Update();
  gPad->RedrawAxis();

  p_ratio->cd();

  TMultiGraph *mg2 = new TMultiGraph("mg2", "mg2 title");


  mg2->Add(g_corr_std_norm_unc_mc_jer_jes);
  mg2->Add(g_corr_std_norm_unc_mc_jer);
  mg2->Add(g_corr_std_norm_unc_mc);
  mg2->Add(g_corr_std_norm_line);

  mg2->Add(g_raw_std_norm_unc_mc);
  mg2->Add(g_raw_std_norm_line);

  mg2->Draw("AL3");
  mg2->SetTitle("");
  mg2->GetHistogram()->GetXaxis()->SetTitle("#it{p}_{T}^{gen} [GeV]");
  mg2->GetHistogram()->GetYaxis()->SetTitle("#sigma_{#it{R}} / #LT #it{R} #GT");
  mg2->GetHistogram()->GetXaxis()->SetTitleSize(text_size / p_ratio_relative_height); // restore original font size
  mg2->GetHistogram()->GetXaxis()->SetLabelSize(text_size / p_ratio_relative_height); // restore original font size
  mg2->GetHistogram()->GetXaxis()->CenterTitle();
  mg2->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
  mg2->GetHistogram()->GetXaxis()->SetLabelOffset(0.007 / p_ratio_relative_height);
  mg2->GetHistogram()->GetYaxis()->CenterTitle();


  mg2->GetHistogram()->GetXaxis()->SetLimits(150., 1250.);

  mg2->GetHistogram()->GetYaxis()->SetTitleSize(text_size / p_ratio_relative_height); // restore original font size
  mg2->GetHistogram()->GetYaxis()->SetTitleOffset(0.81);
  mg2->GetHistogram()->GetYaxis()->SetLabelSize(text_size / p_ratio_relative_height); // restore original font size
  mg2->GetHistogram()->GetYaxis()->SetLabelOffset(0.005 / p_ratio_relative_height);

  double mg2_maximum_saved = mg2->GetHistogram()->GetMaximum();
  double mg2_minimum_saved = mg2->GetHistogram()->GetMinimum();
  mg2->GetHistogram()->SetMaximum(0.17);
  mg2->GetHistogram()->SetMinimum(0.07);
  mg2->GetHistogram()->GetYaxis()->SetNdivisions(506);

  p_ratio->Update();
  tickScaleX = (p_ratio->GetUxmax() - p_ratio->GetUxmin()) / (p_ratio->GetX2() - p_ratio->GetX1()) * (p_ratio->GetWh() * p_ratio->GetAbsHNDC());
  tickScaleY = (p_ratio->GetUymax() - p_ratio->GetUymin()) / (p_ratio->GetY2() - p_ratio->GetY1()) * (p_ratio->GetWw() * p_ratio->GetAbsWNDC());
  mg2->GetHistogram()->GetYaxis()->SetTickLength(p_ratio->GetWw() * tick_length / tickScaleY);
  mg2->GetHistogram()->GetXaxis()->SetTickLength(p_ratio->GetWh() * tick_length / tickScaleX);

  gPad->Update();
  gPad->RedrawAxis();

  // TLatex *down_arrow2 = new TLatex(250, 0.155, "#downarrow JEC applied");
  // down_arrow2->SetTextAlign(13); // left top
  // down_arrow2->SetTextFont(42);
  // down_arrow2->SetTextSize(text_size / p_ratio_relative_height * 0.8);
  // // down_arrow2->SetNDC();
  // down_arrow2->Draw();



  p_main->cd();

  const double leg_x1 = 0.5;
  const double leg_y1 = 0.19;
  const double leg_x2 = 0.87;
  const double leg_y2 = 0.44;
  TLegend *legend = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2);
  // legend->SetNColumns(3);
  legend->SetTextSize(text_size / p_main_relative_height * 0.8);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(g_corr_mean_unc_mc_jer_jes, "JES #oplus JER #oplus MC stat. uncert.", "f");
  legend->AddEntry(g_corr_mean_unc_mc_jer, "JER uncert.", "f");
  legend->AddEntry(g_corr_mean_unc_mc, "MC stat. uncert.", "f");
  legend->Draw();

  TGraphAsymmErrors *g_null = new TGraphAsymmErrors(*g_corr_mean_line);
  g_null->SetLineWidth(0);


  TLegend *legend2 = new TLegend(leg_x1-0.3, leg_y1-0.045, leg_x2-0.3, leg_y2-0.045);
  // legend2->SetNColumns(3);
  legend2->SetTextSize(text_size / p_main_relative_height * 0.8);
  legend2->SetBorderSize(0);
  legend2->SetFillStyle(0);
  legend2->AddEntry(g_corr_mean_line, "JEC applied", "l");
  legend2->AddEntry(g_raw_mean_line, "Uncorrected", "l");
  legend2->AddEntry(g_null, " ", "l");
  legend2->Draw();

  // TLatex *up_arrow = new TLatex(250, 0.977, "#uparrow JEC applied");
  // up_arrow->SetTextAlign(13); // left top
  // up_arrow->SetTextFont(42);
  // up_arrow->SetTextSize(text_size / p_main_relative_height * 0.8);
  // // up_arrow->SetNDC();
  // up_arrow->Draw();
  // TLatex *down_arrow = new TLatex(250, 0.921, "#downarrow Uncorrected");
  // down_arrow->SetTextAlign(11); // left bottom
  // down_arrow->SetTextFont(42);
  // down_arrow->SetTextSize(text_size / p_main_relative_height * 0.8);
  // // down_arrow->SetNDC();
  // down_arrow->Draw();

  TLatex *cms = new TLatex(coord->ConvertGraphXToPadX((0.05*p_main->GetWh()/p_main->GetWw())), coord->ConvertGraphYToPadY(0.87), "CMS");
  cms->SetTextAlign(13); // left top
  cms->SetTextFont(62);
  cms->SetTextSize(0.05 / p_main_relative_height);
  cms->SetNDC();
  cms->Draw();

  TLatex *prelim = new TLatex(coord->ConvertGraphXToPadX((0.05*p_main->GetWh()/p_main->GetWw())), coord->ConvertGraphYToPadY(0.75), "Simulation, Private Work");
  prelim->SetTextAlign(13); // left top
  prelim->SetTextFont(52);
  prelim->SetTextSize(text_size / p_main_relative_height);
  prelim->SetNDC();
  prelim->Draw();

  // const string text_top_right_string = macros::kYears.at(year).lumi_fb_display+" fb^{#minus1} (13 TeV)";
  const string text_top_right_string = "t#bar{t} #rightarrow all-jets @ #sqrt{#it{s}} = 13 TeV";
  TLatex *text_top_right = new TLatex(1. - canvas_margin_r, 1. - canvas_margin_t / p_main_relative_height + 0.015, text_top_right_string.c_str());
  text_top_right->SetTextAlign(31); // right bottom aligned
  text_top_right->SetTextFont(42);
  text_top_right->SetTextSize(text_size / p_main_relative_height);
  text_top_right->SetNDC();
  text_top_right->Draw();

  TLatex *text_top_left = new TLatex(canvas_margin_l, 1. - canvas_margin_t / p_main_relative_height + 0.015, (string("Ultra Legacy ")+kYears.at(year).nice_name).c_str());
  text_top_left->SetTextAlign(11); // left bottom aligned
  text_top_left->SetTextFont(42);
  text_top_left->SetTextSize(text_size / p_main_relative_height);
  text_top_left->SetNDC();
  text_top_left->Draw();

  TLatex *algo_label = new TLatex(coord->ConvertGraphXToPadX(1-(0.05*p_main->GetWh()/p_main->GetWw())), coord->ConvertGraphYToPadY(0.87), (kProbeJetAlgos.at(ProbeJetAlgo::isHOTVR).name+" PUPPI").c_str());
  algo_label->SetTextAlign(33); // right top
  algo_label->SetTextFont(62);
  algo_label->SetTextSize(0.05 / p_main_relative_height); // 0.05
  algo_label->SetTextColor(kGray+1);
  algo_label->SetNDC();
  algo_label->Draw();

  TLatex *eta_text = new TLatex(coord->ConvertGraphXToPadX(1-(0.05*p_main->GetWh()/p_main->GetWw())), coord->ConvertGraphYToPadY(0.75), "#left|#eta^{rec}#right| < 2.5");
  eta_text->SetTextAlign(33); // right top
  eta_text->SetTextFont(42);
  eta_text->SetTextSize(text_size / p_main_relative_height);
  eta_text->SetNDC();
  eta_text->Draw();



  // Save to disk
  const string plotName = (string)"plot_HOTVR_ptResponse_v2_"+kYears.at(year).name + ".pdf";
  const string plotBasePath = "plots";
  gSystem->Exec(((string)"mkdir -p "+plotBasePath).c_str());
  string plotPath = plotBasePath+"/"+plotName;
  canvas->SaveAs(plotPath.c_str());
  // delete c;
}



void plots_HOTVR_ptResponse_v2() {

  const vector<Year> years = { Year::isUL16preVFP, Year::isUL16postVFP, Year::isUL17, Year::isUL18 };
  // const vector<Year> years = { Year::isUL18 };

  for(const Year & year : years) {
    do_plot(year);
  }
}
