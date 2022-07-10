#include "../../constants.h"

using namespace macros;

// class Plotter {
// public:
//   Plotter(
//     const string & tagger_base,
//     const Year & year,
//     const Channel & channel,
//     const bool log_y = false
//   ):
//     fTaggerBase(tagger_base),
//     fYear(year),
//     fChannel(channel),
//     fLogY(log_y)
//   {}
//
// private:
//   string inputDir;
//   string infileName
// };




void plot_postfit(const string & tagger_base, const Year & year, const Channel & channel, const bool log_y = false) {

  const string inputDir = (string)getenv("CMSSW_BASE")+"/src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel/"+kYears.at(year).name+"/"+kChannels.at(channel).name+"/qcd_normalization_study/";
  const string infileName = "prepostfitshapes_"+tagger_base+"_rebinned.root";
  const string infilePath = inputDir+infileName;
  TFile *infile = TFile::Open(infilePath.c_str(), "READ");

  TH1F *hist_Data   = (TH1F*)infile->Get("postfit/data_obs");
  TH1F *hist_Fit    = (TH1F*)infile->Get("postfit/TotalProcs");
  TH1F *hist_nonQCD = (TH1F*)infile->Get("postfit/nonQCD");
  TH1F *hist_QCD_DD = (TH1F*)infile->Get("postfit/QCD_DD");

  hist_Data->SetMarkerStyle(8);
  hist_Data->SetMarkerColor(kBlack);

}





void plot_prefit_QCD_MC(const string & tagger_base, const Year & year, const Channel & channel, const bool log_y = false) {

  const string inputDir = (string)getenv("CMSSW_BASE")+"/src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel/"+kYears.at(year).name+"/"+kChannels.at(channel).name+"/qcd_normalization_study/";
  const string infileName = "qcd_normalization_histograms__"+tagger_base+".root";
  const string infilePath = inputDir+infileName;
  TFile *infile = TFile::Open(infilePath.c_str(), "READ");

  TH1F *hist_Main = (TH1F*)infile->Get("Main/MC.QCD");
  double integral_Main = 0.;
  TH1F *hist_QCD = (TH1F*)infile->Get("QCD/MC.QCD");
  double integral_QCD = 0.;
  const int n_bins = hist_Main->GetNbinsX();
  if(n_bins != hist_QCD->GetNbinsX()) throw runtime_error("Number of bins of QCD and Main hist are not the same!");
  // get integrals
  for (int i_bin = 0; i_bin <= n_bins + 1; i_bin++) {
    integral_Main += hist_Main->GetBinContent(i_bin);
    integral_QCD += hist_QCD->GetBinContent(i_bin);
  }

  // now normalize
  double maximum = 0.;
  for (int i_bin = 0; i_bin <= n_bins + 1; i_bin++) {
    hist_Main->SetBinContent(i_bin, hist_Main->GetBinContent(i_bin) / integral_Main);
    hist_Main->SetBinError(i_bin, hist_Main->GetBinError(i_bin) / integral_Main);
    hist_QCD->SetBinContent(i_bin, hist_QCD->GetBinContent(i_bin) / integral_QCD);
    hist_QCD->SetBinError(i_bin, hist_QCD->GetBinError(i_bin) / integral_QCD);
    maximum = max(maximum, max(hist_Main->GetBinContent(i_bin), hist_QCD->GetBinContent(i_bin)));
  }

  TGraphAsymmErrors *graph_Main = new TGraphAsymmErrors(hist_Main);
  TGraphAsymmErrors *graph_QCD = new TGraphAsymmErrors(hist_QCD);

  // shift the graph markers so they do not overlap
  for (int i_bin = 0; i_bin < graph_Main->GetN(); i_bin++) {
    double x, y, x_err_low, x_err_high;
    graph_Main->GetPoint(i_bin, x, y);
    x_err_low = graph_Main->GetErrorXlow(i_bin);
    x_err_high = graph_Main->GetErrorXhigh(i_bin);
    double x_offset = x_err_low / 3.;
    graph_Main->SetPoint(i_bin, x - x_offset, y);
    graph_Main->SetPointEXlow(i_bin, x_err_low - x_offset);
    graph_Main->SetPointEXhigh(i_bin, x_err_high + x_offset);
    graph_QCD->GetPoint(i_bin, x, y);
    graph_QCD->SetPoint(i_bin, x + x_offset, y);
    graph_QCD->SetPointEXlow(i_bin, x_err_low + x_offset);
    graph_QCD->SetPointEXhigh(i_bin, x_err_high - x_offset);
  }


  float marker_size = 0.8;
  graph_Main->SetLineColor(kColors.at("pyplot_blue"));
  graph_Main->SetMarkerColor(kColors.at("pyplot_blue"));
  graph_Main->SetMarkerStyle(24);
  graph_Main->SetMarkerSize(marker_size);
  // graph_Main->SetLineStyle(1);
  graph_QCD->SetLineColor(kColors.at("pyplot_orange"));
  graph_QCD->SetMarkerColor(kColors.at("pyplot_orange"));
  graph_QCD->SetMarkerStyle(25);
  graph_QCD->SetMarkerSize(marker_size);
  // graph_QCD->SetLineStyle(2);

  TCanvas *canvas = new TCanvas("canvas", "canvas title", 600, 600);
  canvas->cd();
  float margin_l = 0.15;
  float margin_r = 0.05;
  float margin_b = 0.12;
  float margin_t = 0.08;
  float tick_length = 0.015; // fraction of canvas width/height
  canvas->SetMargin(margin_l, margin_r, margin_b, margin_t); //lrbt
  auto coord = new CoordinateConverter();
  coord->init(margin_l, margin_r, margin_b, margin_t);
  if(log_y) canvas->SetLogy();
  canvas->SetTickx(1);
  canvas->SetTicky(1);

  TMultiGraph *mg = new TMultiGraph("mg", "");

  mg->Add(graph_Main);
  mg->Add(graph_QCD);
  mg->Draw("apz");

  TH1F *draw_hist = mg->GetHistogram();

  draw_hist->SetMaximum(log_y ? 1.: maximum*1.2);
  draw_hist->SetMinimum(log_y ? 0.0001 : 0.);

  draw_hist->GetXaxis()->SetTitle("#it{M}_{T, W} [GeV]");
  draw_hist->GetXaxis()->SetTitleOffset(1.3);
  draw_hist->GetXaxis()->SetLimits(hist_Main->GetXaxis()->GetBinLowEdge(1), hist_Main->GetXaxis()->GetBinUpEdge(n_bins));

  draw_hist->GetYaxis()->SetTitle("Normalized number of events / 5 GeV");
  draw_hist->GetYaxis()->SetTitleOffset(1.7);

  canvas->Update();
  double tickScaleX = (canvas->GetUxmax() - canvas->GetUxmin()) / (canvas->GetX2() - canvas->GetX1()) * (canvas->GetWh() * canvas->GetAbsHNDC());
  double tickScaleY = (canvas->GetUymax() - canvas->GetUymin()) / (canvas->GetY2() - canvas->GetY1()) * (canvas->GetWw() * canvas->GetAbsWNDC());
  draw_hist->GetXaxis()->SetTickLength(canvas->GetWh() * tick_length / tickScaleX);
  draw_hist->GetYaxis()->SetTickLength(canvas->GetWw() * tick_length / tickScaleY);
  gPad->Update();
  gPad->RedrawAxis();

  TLatex *text_top_left = new TLatex(margin_l, 1-(margin_t-0.01), (string("Ultra Legacy ")+kYears.at(year).nice_name).c_str());
  text_top_left->SetTextAlign(11); // left bottom aligned
  text_top_left->SetTextFont(42);
  text_top_left->SetTextSize(0.035);
  text_top_left->SetNDC();
  text_top_left->Draw();

  string string_text_top_right = kYears.at(year).lumi_fb_display+" fb^{#minus1} (13 TeV)";
  TLatex *text_top_right = new TLatex(1-margin_r, 1-(margin_t-0.01), string_text_top_right.c_str());
  text_top_right->SetTextAlign(31); // right bottom aligned
  text_top_right->SetTextFont(42);
  text_top_right->SetTextSize(0.035);
  text_top_right->SetNDC();
  text_top_right->Draw();

  TLatex *cms = new TLatex(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.95), "CMS");
  cms->SetTextAlign(13); // left top
  cms->SetTextFont(62);
  cms->SetTextSize(0.05);
  cms->SetNDC();
  cms->Draw();

  // TText *prelim = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.87), "Simulation, Private Work");
  // prelim->SetTextAlign(13); // left top
  TText *prelim = new TText(coord->ConvertGraphXToPadX(0.2), coord->ConvertGraphYToPadY(0.905), "Simulation, Private Work");
  prelim->SetTextAlign(11); // left bottom
  prelim->SetTextFont(52);
  prelim->SetTextSize(0.035);
  prelim->SetNDC();
  prelim->Draw();

  string pt_text_string;
  string pt_text_string2;
  if(tagger_base == "ak8_t") {
    // pt_text_string = "AK8 probe #it{p}_{T} > 300 GeV";
    pt_text_string = "T&P selection w/ AK8 probe";
    pt_text_string2 = "#it{p}_{T, probe} > 300 GeV, |#eta_{probe}| < 2.5";
  }
  else if(tagger_base == "ak8_w") {
    // pt_text_string = "AK8 probe #it{p}_{T} > 200 GeV";
    pt_text_string = "T&P selection w/ AK8 probe";
    pt_text_string2 = "#it{p}_{T, probe} > 200 GeV, |#eta_{probe}| < 2.5";
  }
  else if(tagger_base == "hotvr_t") {
    // pt_text_string = "HOTVR probe #it{p}_{T} > 200 GeV";
    pt_text_string = "T&P selection w/ HOTVR probe";
    pt_text_string2 = "#it{p}_{T, probe} > 200 GeV, |#eta_{probe}| < 2.5";
  }

  TLatex *pt_text = new TLatex(coord->ConvertGraphXToPadX(0.95), coord->ConvertGraphYToPadY(0.735), pt_text_string.c_str()); //.945
  pt_text->SetTextAlign(33); // left bottom
  pt_text->SetTextFont(42);
  pt_text->SetTextSize(0.03);
  pt_text->SetNDC();
  pt_text->Draw();

  TLatex *pt_text2 = new TLatex(coord->ConvertGraphXToPadX(0.95), coord->ConvertGraphYToPadY(0.65), pt_text_string2.c_str()); //.885
  pt_text2->SetTextAlign(33); // left bottom
  pt_text2->SetTextFont(42);
  pt_text2->SetTextSize(0.03);
  pt_text2->SetNDC();
  pt_text2->Draw();

  // TGraphAsymmErrors *null_graph = new TGraphAsymmErrors(*graph_Main);
  // null_graph->SetLineColor(kWhite);

  float legend_x1 = coord->ConvertGraphXToPadX(0.05);
  float legend_y1 = coord->ConvertGraphYToPadY(0.05);
  float legend_x2 = coord->ConvertGraphXToPadX(0.45);
  float legend_y2 = coord->ConvertGraphYToPadY(0.3);
  TLegend *legend = new TLegend(legend_x1, legend_y1, legend_x2, legend_y2); // x1 y1 x2 y2
  legend->SetHeader("#bf{Pre-fit QCD multijet MC}");
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(graph_Main, "signal region, isolated muon"); //passing muon isolation criterion
  legend->AddEntry(graph_QCD, "sideband region, non-iso. muon");
  // legend->AddEntry(null_graph, "\\ldots custom muon isolation", "l");
  legend->Draw();


  double chi2;
  int ndf;
  int igood;
  hist_Main->Chi2TestX(hist_QCD, chi2, ndf, igood, "WW");
  cout << chi2 << endl;
  cout << ndf << endl;
  cout << igood << endl;
  // hist_QCD->Chi2TestX(hist_Main, chi2, ndf, igood, "WW");
  // cout << chi2 << endl;
  // cout << ndf << endl;
  // cout << igood << endl;

  string chi2_string = "#chi^{2} / n.d.f. = ";
  stringstream ss;
  ss << fixed << setprecision(1) << chi2;
  chi2_string += ss.str() + " / " + to_string(ndf);
  TLatex *chi2_text = new TLatex(coord->ConvertGraphXToPadX(0.95), coord->ConvertGraphYToPadY(.5), chi2_string.c_str());
  chi2_text->SetTextAlign(33); // left bottom
  chi2_text->SetTextFont(42);
  chi2_text->SetTextSize(0.03);
  chi2_text->SetNDC();
  chi2_text->Draw();

  // Save to disk
  string plotName = "plot_QCDestimate_CompareShapes__"+tagger_base+"__"+kYears.at(year).name+"_"+kChannels.at(channel).name+"__";
  // for(bool do_x : {do_raw, do_msd, do_msd_btag}) { do_x ? plotName += "1" : plotName += "0"; }
  plotName += (log_y ? string("log") : string("lin"))+".pdf";
  const string plotDir = inputDir+"plots/";
  gSystem->Exec(((string)"mkdir -p "+plotDir).c_str());
  const string plotPath = plotDir+plotName;
  canvas->SaveAs(plotPath.c_str());
  delete canvas;
}


void plot_CompareShapes() {
  const set<string> taggers = { "ak8_t", "hotvr_t", "ak8_w" };
  for (const string & tagger : taggers) {
    for (const Year & year : kAllYears) {
      plot_prefit_QCD_MC(tagger, year, Channel::isMuo, true);
    }
  }
  // plot_prefit_QCD_MC("hotvr_t", Year::isUL18, Channel::isMuo, true);
}
