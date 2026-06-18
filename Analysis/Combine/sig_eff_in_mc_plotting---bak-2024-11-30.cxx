#include <iomanip>
#include <sstream>
#include <array>

#include "../constants.h"

using namespace macros;

const string data_path = "sig_eff_in_mc_BasicHists";




TGraph * new_null_graph() {
  int n = 1001; // should be the same number + 1 as the variable "nx" in analyze.py
  double x[n], y[n];
  for(int i = 0; i <= n; i++) {
    x[i] = i/double(n);
    y[i] = i/double(n);
  }
  TGraph * graph = new TGraph(n, x, y);
  graph->SetLineColorAlpha(kBlack, 0); // 0 = fully transparent
  return graph;
}



// {
//    double x[100], y[100];
//    int n = 20;
//    for (int i=0;i<n;i++) {
//      x[i] = i*0.1;
//      y[i] = 10*sin(x[i]+0.2);
//    }
//    auto g = new TGraph(n,x,y);
//    g->SetTitle("Graph title;X title;Y title");
//    g->Draw("AC*");
// }

//   double y[6] = {3, 8, 1, 10, 5, 7};


//________________________________________________________________________________
// Leonidas' values from his CSV Excel sheets:
// pt bin edges W: 200/300/400/800
// pt bin edges top: 300/400/480/600/1200
static const double leo_ptCenter_W[3] = {250, 350, 600};
static const double leo_ptCenter_top[4] = {350, 440, 540, 700};
// WvsQCD
typedef map<Year, std::array<double, 3>> values_leo_type_3;
typedef map<Year, std::array<double, 4>> values_leo_type_4;
static const values_leo_type_3 values_leo_WvsQCD_BkgEff0p050 = {
  { Year::isUL16preVFP, { 0.8425, 0.7367, 0.7009 } },
  { Year::isUL16postVFP, { 0.8426, 0.7348, 0.7009 } },
  { Year::isUL17, { 0.8377, 0.7297, 0.6967 } },
  { Year::isUL18, { 0.8391, 0.7350, 0.6997 } },
};
static const values_leo_type_3 values_leo_WvsQCD_BkgEff0p010 = {
  { Year::isUL16preVFP, { 0.6880, 0.5918, 0.5630 } },
  { Year::isUL16postVFP, { 0.6851, 0.5875, 0.5591 } },
  { Year::isUL17, { 0.6784, 0.5829, 0.5554 } },
  { Year::isUL18, { 0.6854, 0.5920, 0.5646 } },
};
static const values_leo_type_3 values_leo_WvsQCD_BkgEff0p005 = {
  { Year::isUL16preVFP, { 0.5493, 0.4649, 0.4436 } },
  { Year::isUL16postVFP, { 0.5410, 0.4609, 0.4410 } },
  { Year::isUL17, { 0.5369, 0.4586, 0.4418 } },
  { Year::isUL18, { 0.5143, 0.4394, 0.4259 } },
};
// WvsQCDMD
static const values_leo_type_3 values_leo_WvsQCDMD_BkgEff0p025 = {
  { Year::isUL16preVFP, { 0.7231, 0.6797, 0.6459 } },
  { Year::isUL16postVFP, { 0.7124, 0.6696, 0.6419 } },
  { Year::isUL17, { 0.7249, 0.6723, 0.6395 } },
  { Year::isUL18, { 0.7269, 0.6812, 0.6485 } },
};
static const values_leo_type_3 values_leo_WvsQCDMD_BkgEff0p010 = {
  { Year::isUL16preVFP, { 0.5139, 0.5127, 0.4929 } },
  { Year::isUL16postVFP, { 0.5079, 0.5107, 0.4900 } },
  { Year::isUL17, { 0.5218, 0.5122, 0.4889 } },
  { Year::isUL18, { 0.5191, 0.5172, 0.4980 } },
};
static const values_leo_type_3 values_leo_WvsQCDMD_BkgEff0p005 = {
  { Year::isUL16preVFP, { 0.3701, 0.3895, 0.3855 } },
  { Year::isUL16postVFP, { 0.3677, 0.3903, 0.3844 } },
  { Year::isUL17, { 0.3744, 0.3865, 0.3797 } },
  { Year::isUL18, { 0.3647, 0.3845, 0.3819 } },
};
// TvsQCD
static const values_leo_type_4 values_leo_TvsQCD_BkgEff0p010 = {
  { Year::isUL16preVFP, { 0.8942, 0.9219, 0.9050, 0.8546 } },
  { Year::isUL16postVFP, { 0.8875, 0.9217, 0.9045, 0.8572 } },
  { Year::isUL17, { 0.8747, 0.9209, 0.9025, 0.8672 } },
  { Year::isUL18, { 0.8802, 0.9182, 0.9022, 0.8617 } },
};
static const values_leo_type_4 values_leo_TvsQCD_BkgEff0p005 = {
  { Year::isUL16preVFP, { 0.8575, 0.8952, 0.8812, 0.8326 } },
  { Year::isUL16postVFP, { 0.8502, 0.8946, 0.8770, 0.8358 } },
  { Year::isUL17, { 0.8329, 0.8936, 0.8779, 0.8449 } },
  { Year::isUL18, { 0.8401, 0.8907, 0.8770, 0.8402 } },
};
static const values_leo_type_4 values_leo_TvsQCD_BkgEff0p001 = {
  { Year::isUL16preVFP, { 0.7279, 0.7912, 0.7869, 0.7436 } },
  { Year::isUL16postVFP, { 0.7189, 0.7913, 0.7824, 0.7462 } },
  { Year::isUL17, { 0.7018, 0.7899, 0.7833, 0.7581 } },
  { Year::isUL18, { 0.7106, 0.7882, 0.7837, 0.7559 } },
};


void do_plot_nsubjettiness(const Year & year, const bool bool_prelim) {

  int debug = 0;

  const string inputGraphName = "graph_sig_eff_withWindow_both";

  const string infilePath_nominal = data_path + "/sig_eff_in_mc-ak8_t__tau.root";
  const string infilePath_btag = data_path + "/sig_eff_in_mc-ak8_t_btagDJet__tau.root";
  const string infilePath_hotvr = data_path + "/sig_eff_in_mc-hotvr_t__tau.root";

  TFile* infile_nominal = TFile::Open(infilePath_nominal.c_str(), "READ");
  TFile* infile_btag = TFile::Open(infilePath_btag.c_str(), "READ");
  TFile* infile_hotvr = TFile::Open(infilePath_hotvr.c_str(), "READ");

  const int line_width = 2;
  const int marker_size = 1;
  const int line_style_nominal = 1;
  const int line_style_btag = 2;

  TMultiGraph *mg = new TMultiGraph("mg", "mg title");

  cout << debug++ << endl;

  const string graphName_nominal_vt = kYears.at(year).name + "-BkgEff0p001/"+inputGraphName;
  TGraph graph_nominal_vt = *(TGraph*)infile_nominal->Get(graphName_nominal_vt.c_str());
  graph_nominal_vt.SetLineWidth(line_width);
  graph_nominal_vt.SetLineColor(kGreen);
  graph_nominal_vt.SetLineStyle(line_style_nominal);
  graph_nominal_vt.SetMarkerColor(kGreen);
  graph_nominal_vt.SetMarkerSize(marker_size);
  graph_nominal_vt.SetMarkerStyle(20);
  mg->Add(&graph_nominal_vt);

  cout << debug++ << endl;
  
  const string graphName_nominal_ti = kYears.at(year).name + "-BkgEff0p005/"+inputGraphName;
  TGraph graph_nominal_ti = *(TGraph*)infile_nominal->Get(graphName_nominal_ti.c_str());
  graph_nominal_ti.SetLineWidth(line_width);
  graph_nominal_ti.SetLineColor(kCyan);
  graph_nominal_ti.SetLineStyle(line_style_nominal);
  graph_nominal_ti.SetMarkerColor(kCyan);
  graph_nominal_ti.SetMarkerSize(marker_size);
  graph_nominal_ti.SetMarkerStyle(21);
  mg->Add(&graph_nominal_ti);

  cout << debug++ << endl;

  const string graphName_nominal_me = kYears.at(year).name + "-BkgEff0p010/"+inputGraphName;
  TGraph graph_nominal_me = *(TGraph*)infile_nominal->Get(graphName_nominal_me.c_str());
  graph_nominal_me.SetLineWidth(line_width);
  graph_nominal_me.SetLineColor(kBlue);
  graph_nominal_me.SetLineStyle(line_style_nominal);
  graph_nominal_me.SetMarkerColor(kBlue);
  graph_nominal_me.SetMarkerSize(marker_size);
  graph_nominal_me.SetMarkerStyle(22);
  mg->Add(&graph_nominal_me);

  cout << debug++ << endl;

  const string graphName_nominal_lo = kYears.at(year).name + "-BkgEff0p025/"+inputGraphName;
  TGraph graph_nominal_lo = *(TGraph*)infile_nominal->Get(graphName_nominal_lo.c_str());
  graph_nominal_lo.SetLineWidth(line_width);
  graph_nominal_lo.SetLineColor(kMagenta);
  graph_nominal_lo.SetLineStyle(line_style_nominal);
  graph_nominal_lo.SetMarkerColor(kMagenta);
  graph_nominal_lo.SetMarkerSize(marker_size);
  graph_nominal_lo.SetMarkerStyle(23);
  mg->Add(&graph_nominal_lo);

  cout << debug++ << endl;

  const string graphName_nominal_vl = kYears.at(year).name + "-BkgEff0p050/"+inputGraphName;
  TGraph graph_nominal_vl = *(TGraph*)infile_nominal->Get(graphName_nominal_vl.c_str());
  graph_nominal_vl.SetLineWidth(line_width);
  graph_nominal_vl.SetLineColor(kRed);
  graph_nominal_vl.SetLineStyle(line_style_nominal);
  graph_nominal_vl.SetMarkerColor(kRed);
  graph_nominal_vl.SetMarkerSize(marker_size);
  graph_nominal_vl.SetMarkerStyle(47);
  mg->Add(&graph_nominal_vl);

  cout << debug++ << endl;




  const string graphName_btag_vt = kYears.at(year).name + "-BkgEff0p001/"+inputGraphName;
  TGraph graph_btag_vt = *(TGraph*)infile_btag->Get(graphName_btag_vt.c_str());
  graph_btag_vt.SetLineWidth(line_width);
  graph_btag_vt.SetLineColor(kGreen);
  graph_btag_vt.SetLineStyle(line_style_btag);
  graph_btag_vt.SetMarkerColor(kGreen);
  graph_btag_vt.SetMarkerSize(marker_size);
  graph_btag_vt.SetMarkerStyle(24);
  mg->Add(&graph_btag_vt);

  cout << debug++ << endl;
  
  const string graphName_btag_ti = kYears.at(year).name + "-BkgEff0p005/"+inputGraphName;
  TGraph graph_btag_ti = *(TGraph*)infile_btag->Get(graphName_btag_ti.c_str());
  graph_btag_ti.SetLineWidth(line_width);
  graph_btag_ti.SetLineColor(kCyan);
  graph_btag_ti.SetLineStyle(line_style_btag);
  graph_btag_ti.SetMarkerColor(kCyan);
  graph_btag_ti.SetMarkerSize(marker_size);
  graph_btag_ti.SetMarkerStyle(25);
  mg->Add(&graph_btag_ti);

  cout << debug++ << endl;

  const string graphName_btag_me = kYears.at(year).name + "-BkgEff0p010/"+inputGraphName;
  TGraph graph_btag_me = *(TGraph*)infile_btag->Get(graphName_btag_me.c_str());
  graph_btag_me.SetLineWidth(line_width);
  graph_btag_me.SetLineColor(kBlue);
  graph_btag_me.SetLineStyle(line_style_btag);
  graph_btag_me.SetMarkerColor(kBlue);
  graph_btag_me.SetMarkerSize(marker_size);
  graph_btag_me.SetMarkerStyle(26);
  mg->Add(&graph_btag_me);

  cout << debug++ << endl;

  const string graphName_btag_lo = kYears.at(year).name + "-BkgEff0p025/"+inputGraphName;
  TGraph graph_btag_lo = *(TGraph*)infile_btag->Get(graphName_btag_lo.c_str());
  graph_btag_lo.SetLineWidth(line_width);
  graph_btag_lo.SetLineColor(kMagenta);
  graph_btag_lo.SetLineStyle(line_style_btag);
  graph_btag_lo.SetMarkerColor(kMagenta);
  graph_btag_lo.SetMarkerSize(marker_size);
  graph_btag_lo.SetMarkerStyle(32);
  mg->Add(&graph_btag_lo);

  cout << debug++ << endl;

  const string graphName_btag_vl = kYears.at(year).name + "-BkgEff0p050/"+inputGraphName;
  TGraph graph_btag_vl = *(TGraph*)infile_btag->Get(graphName_btag_vl.c_str());
  graph_btag_vl.SetLineWidth(line_width);
  graph_btag_vl.SetLineColor(kRed);
  graph_btag_vl.SetLineStyle(line_style_btag);
  graph_btag_vl.SetMarkerColor(kRed);
  graph_btag_vl.SetMarkerSize(marker_size);
  graph_btag_vl.SetMarkerStyle(46);
  mg->Add(&graph_btag_vl);

  cout << debug++ << endl;




  // HOTVR

  const string graphName_hotvr = kYears.at(year).name + "-Standard/"+inputGraphName;
  TGraph graph_hotvr = *(TGraph*)infile_hotvr->Get(graphName_hotvr.c_str());
  graph_hotvr.SetLineWidth(line_width);
  graph_hotvr.SetLineColor(kBlack);
  graph_hotvr.SetLineStyle(line_style_nominal);
  graph_hotvr.SetMarkerColor(kBlack);
  graph_hotvr.SetMarkerSize(marker_size);
  graph_hotvr.SetMarkerStyle(34);
  mg->Add(&graph_hotvr);

  cout << debug++ << endl;


  // Canvas

  TCanvas *c = new TCanvas("canvas", "canvas title", 600, 600);
  c->cd();
  float margin_l = 0.15;
  float margin_r = 0.05;
  float margin_b = 0.12;
  float margin_t = 0.08;
  float tick_length = 0.015; // fraction of canvas width/height
  c->SetMargin(margin_l, margin_r, margin_b, margin_t); //lrbt
  auto coord = new CoordinateConverter();
  coord->init(margin_l, margin_r, margin_b, margin_t);
  // if(log_y) c->SetLogy();
  c->SetTickx(1);
  c->SetTicky(1);


  bool plot_at_bottom = true;


  mg->Draw("alp");
  mg->SetTitle("");
  // mg->GetHistogram()->GetXaxis()->SetRangeUser(300., 1000.);
  mg->GetHistogram()->GetXaxis()->SetLimits(200., 1200.);
  double maximum_saved = mg->GetHistogram()->GetMaximum();
  double minimum_saved = mg->GetHistogram()->GetMinimum();
  double new_maximum = maximum_saved * 1.618;
  // mg->GetHistogram()->SetMaximum(1.);
  // mg->GetHistogram()->SetMaximum(new_maximum);
  mg->GetHistogram()->SetMaximum(1.4);
  mg->GetHistogram()->SetMinimum(0.);
  // if(is_eff_bkg || is_eff_sig) mg->GetHistogram()->GetXaxis()->SetTitle(tagger.scan_var_tlatex_axis.c_str());
  // else if(is_roc) mg->GetHistogram()->GetXaxis()->SetTitle("True positive rate, #varepsilon_{S}");
  // if(is_eff_bkg || is_roc) mg->GetHistogram()->GetYaxis()->SetTitle("False positive rate, #varepsilon_{B}");
  // if(is_eff_sig) mg->GetHistogram()->GetYaxis()->SetTitle("True positive rate, #varepsilon_{S}");
  mg->GetHistogram()->GetXaxis()->SetTitle("Probe jet #it{p}_{T} [GeV]");
  mg->GetHistogram()->GetYaxis()->SetTitle("Signal efficiency");
  mg->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
  mg->GetHistogram()->GetXaxis()->CenterTitle();
  mg->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);
  mg->GetHistogram()->GetYaxis()->CenterTitle();
  // mg->GetHistogram()->GetXaxis()->SetNdivisions(405);
  // if(!log_y) mg->GetHistogram()->GetYaxis()->SetNdivisions(405);
  gStyle->SetLineWidth(1);

  // legend->Draw();





  // Lines
  // TLine(Double_t x1,Double_t y1,Double_t x2,Double_t y2)

  int lineX_color = kGray+2;
  int lineX_style = 3;
  int lineX_width = 1;

  // double line_max_y = maximum_saved * 1.05;
  double line_max_y = 1.;

  TLine* line1 = new TLine(250, 0.4, 250, 0.5);
  line1->SetLineColor(lineX_color);
  line1->SetLineStyle(lineX_style);
  line1->SetLineWidth(lineX_width);
  line1->Draw();
  TLine* line2 = new TLine(300, 0, 300, line_max_y);
  line2->SetLineColor(lineX_color);
  line2->SetLineStyle(lineX_style);
  line2->SetLineWidth(lineX_width);
  line2->Draw();
  TLine* line3 = new TLine(350, 0.4, 350, 0.5);
  line3->SetLineColor(lineX_color);
  line3->SetLineStyle(lineX_style);
  line3->SetLineWidth(lineX_width);
  line3->Draw();
  TLine* line4 = new TLine(400, 0, 400, line_max_y);
  line4->SetLineColor(lineX_color);
  line4->SetLineStyle(lineX_style);
  line4->SetLineWidth(lineX_width);
  line4->Draw();
  TLine* line5 = new TLine(480, 0, 480, line_max_y);
  line5->SetLineColor(lineX_color);
  line5->SetLineStyle(lineX_style);
  line5->SetLineWidth(lineX_width);
  line5->Draw();
  TLine* line6 = new TLine(600, 0, 600, line_max_y);
  line6->SetLineColor(lineX_color);
  line6->SetLineStyle(lineX_style);
  line6->SetLineWidth(lineX_width);
  line6->Draw();



  // Legend

  float leg_x1 = 0.62;
  float leg_y1 = 0.15;
  float leg_x2 = 0.83;
  float leg_y2 = 0.55;

  // TLegend *legend = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2, '', 'brNDC');
  TLegend *legend = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2);

  legend->AddEntry(&graph_nominal_vl, "AK8, #tau_{3}/#tau_{2} < 0.69", "p");
  legend->AddEntry(&graph_nominal_lo, "AK8, #tau_{3}/#tau_{2} < 0.61", "p");
  legend->AddEntry(&graph_hotvr, "HOTVR, #tau_{3}/#tau_{2} < 0.56", "p");
  legend->AddEntry(&graph_nominal_me, "AK8, #tau_{3}/#tau_{2} < 0.52", "p");
  legend->AddEntry(&graph_nominal_ti, "AK8, #tau_{3}/#tau_{2} < 0.47", "p");
  legend->AddEntry(&graph_nominal_vt, "AK8, #tau_{3}/#tau_{2} < 0.38", "p");

  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0); // Set background to transparent

  legend->Draw();


  // Legend nominal / subjet b tag

  float leg2_x1 = 0.44;
  float leg2_y1 = 0.72;
  float leg2_x2 = 0.77;
  float leg2_y2 = 0.87;

  // TLegend *legend2 = new TLegend(leg2_x1, leg2_y1, leg2_x2, leg2_y2, '', 'brNDC');
  TLegend *legend2 = new TLegend(leg2_x1, leg2_y1, leg2_x2, leg2_y2);

  TGraph *graph_nominal_dummy = new TGraph();
  graph_nominal_dummy->SetLineWidth(line_width);
  graph_nominal_dummy->SetLineColor(kBlack);
  graph_nominal_dummy->SetLineStyle(line_style_nominal);
  graph_nominal_dummy->SetMarkerColor(kBlack);
  graph_nominal_dummy->SetMarkerSize(marker_size);
  graph_nominal_dummy->SetMarkerStyle(20);

  TGraph *graph_btag_dummy = new TGraph();
  graph_btag_dummy->SetLineWidth(line_width);
  graph_btag_dummy->SetLineColor(kBlack);
  graph_btag_dummy->SetLineStyle(line_style_btag);
  graph_btag_dummy->SetMarkerColor(kBlack);
  graph_btag_dummy->SetMarkerSize(marker_size);
  graph_btag_dummy->SetMarkerStyle(24);

  legend2->SetHeader("#bf{t tagging using N-subjettiness}");
  legend2->AddEntry(graph_nominal_dummy, "w/o subjet b tagging");
  legend2->AddEntry(graph_btag_dummy, "w/ loose DeepJet subjet b tag");
  
  legend2->SetTextSize(0.03);
  legend2->SetBorderSize(0);
  legend2->SetFillStyle(0); // Set background to transparent

  legend2->Draw();



  // TLatex *text_top_left = new TLatex(margin_l, 1-(margin_t-0.01), (string("Ultra Legacy ")+kYears.at(year).nice_name).c_str());
  // text_top_left->SetTextAlign(11); // left bottom aligned
  // text_top_left->SetTextFont(42);
  // text_top_left->SetTextSize(0.035);
  // text_top_left->SetNDC();
  // text_top_left->Draw();

  // TLatex *algo_label = new TLatex(coord->ConvertGraphXToPadX(1-(0.05*c->GetWh()/c->GetWw())), coord->ConvertGraphYToPadY(0.95), "AK8 PUPPI");
  // algo_label->SetTextAlign(33); // right top
  // algo_label->SetTextFont(62);
  // algo_label->SetTextSize(0.05); // 0.05
  // algo_label->SetTextColor(kGray+1);
  // algo_label->SetNDC();
  // algo_label->Draw();

  // const double eta_pt_text_x = coord->ConvertGraphXToPadX(1-(0.05*c->GetWh()/c->GetWw()));
  // double eta_text_y = coord->ConvertGraphYToPadY(0.801);
  // // double pt_text_y = coord->ConvertGraphYToPadY(0.867);
  // if(!plot_at_bottom) {
  //   eta_text_y = coord->ConvertGraphYToPadY((0.05*c->GetWh()/c->GetWw()) + 0.066);
  //   // pt_text_y = coord->ConvertGraphYToPadY((0.05*c->GetWh()/c->GetWw()));
  // }

  // TLatex *eta_text = new TLatex(eta_pt_text_x, eta_text_y, "|#eta| < 2.5");
  // eta_text->SetTextAlign(plot_at_bottom ? 33 : 31); // right top
  // eta_text->SetTextFont(42);
  // eta_text->SetTextSize(0.035);
  // eta_text->SetNDC();
  // eta_text->Draw();

  // TLatex *pt_text = new TLatex(eta_pt_text_x, pt_text_y, pt_bin.text.c_str());
  // pt_text->SetTextAlign(plot_at_bottom ? 33 : 31); // right top
  // pt_text->SetTextFont(42);
  // pt_text->SetTextSize(0.035);
  // pt_text->SetNDC();
  // pt_text->Draw();

  // string string_text_top_right = "t#bar{t} #rightarrow all-jets @ #sqrt{#it{s}} = 13 TeV";
  // string string_text_top_right = "#sqrt{#it{s}} = 13 TeV";
  string string_text_top_right = kYears.at(year).lumi_fb_display+" fb^{#minus1} (13 TeV)";
  TLatex *text_top_right = new TLatex(1-margin_r, 1-(margin_t-0.01), string_text_top_right.c_str());
  text_top_right->SetTextAlign(31); // right bottom aligned
  text_top_right->SetTextFont(42);
  text_top_right->SetTextSize(0.035);
  text_top_right->SetNDC();
  text_top_right->Draw();

  TText *prelimPW1 = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.95), "Private work");
  prelimPW1->SetTextAlign(13); // left top
  prelimPW1->SetTextFont(52);
  prelimPW1->SetTextSize(0.03);
  prelimPW1->SetNDC();
  // prelimPW1->SetTextColor(kWhite);
  if(!bool_prelim) prelimPW1->Draw();

  TText *prelimPW2 = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.90), "(CMS data/simulation)");
  prelimPW2->SetTextAlign(13); // left top
  prelimPW2->SetTextFont(52);
  prelimPW2->SetTextSize(0.03);
  prelimPW2->SetNDC();
  // prelimPW2->SetTextColor(kWhite);
  if(!bool_prelim) prelimPW2->Draw();

  TLatex *cms = new TLatex(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.95), "CMS");
  cms->SetTextAlign(13); // left top
  cms->SetTextFont(62);
  cms->SetTextSize(0.05);
  cms->SetNDC();
  if(bool_prelim) cms->Draw();

  if(!plot_at_bottom) {
    TString prelim_text;
    if(bool_prelim) prelim_text = "Simulation, Preliminary";
    else prelim_text = "Simulation, Private Work";
    TText *prelim = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.87), prelim_text);
    prelim->SetTextAlign(13); // left top
    prelim->SetTextFont(52);
    prelim->SetTextSize(0.035);
    prelim->SetNDC();
    if(bool_prelim) prelim->Draw();
  }
  else {
    TText *prelim = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.87), "Simulation");
    prelim->SetTextAlign(13); // left top
    prelim->SetTextFont(52);
    prelim->SetTextSize(0.035);
    prelim->SetNDC();
    // prelim->SetTextColor(kWhite);
    if(bool_prelim) prelim->Draw();

    TString prelim_text;
    if(bool_prelim) prelim_text = "Preliminary";
    else prelim_text = "Private Work";
    TText *prelim2 = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.804), prelim_text);
    prelim2->SetTextAlign(13); // left top
    prelim2->SetTextFont(52);
    prelim2->SetTextSize(0.035);
    prelim2->SetNDC();
    // prelim2->SetTextColor(kWhite);
    if(bool_prelim) prelim2->Draw();
  }


  // https://root-forum.cern.ch/t/inconsistent-tick-length/18563/8
  c->Update();
  double tickScaleX = (c->GetUxmax() - c->GetUxmin()) / (c->GetX2() - c->GetX1()) * (c->GetWh() * c->GetAbsHNDC());
  double tickScaleY = (c->GetUymax() - c->GetUymin()) / (c->GetY2() - c->GetY1()) * (c->GetWw() * c->GetAbsWNDC());
  mg->GetHistogram()->GetYaxis()->SetTickLength(c->GetWw() * tick_length / tickScaleY);
  mg->GetHistogram()->GetXaxis()->SetTickLength(c->GetWh() * tick_length / tickScaleX);

  gPad->Update();
  gPad->RedrawAxis();

  // Save to disk
  string plotName;
  if(bool_prelim) plotName = "plot_sig_eff_in_mc__ak8__"+kYears.at(year).name+".pdf";
  else plotName = "plot_sig_eff_in_mc__ak8__"+kYears.at(year).name+"-PrivateWork.pdf";
  // string plotName = "plot__"+tagger.name_base+"__"+pt_bin.name+"__"+graph_base_name+"__";
  // for(bool do_x : {do_raw, do_msd, do_msd_btag}) { do_x ? plotName += "1" : plotName += "0"; }
  // plotName += (log_y ? string("log") : string("lin"))+".pdf";
  string plotPath = data_path+"/"+plotName;
  c->SaveAs(plotPath.c_str());
  delete c;
}

void do_plot_partnet_w(const Year & year, const bool bool_prelim) {

  int debug = 0;

  const int line_width = 2;
  const int marker_size = 1;
  const int line_style_nominal = 1;
  const int line_style_btag = 2;

  TMultiGraph *mg = new TMultiGraph("mg", "mg title");

  cout << debug++ << endl;

  TGraph graph_nominal_BkgEff0p005 = TGraph(3, leo_ptCenter_W, values_leo_WvsQCD_BkgEff0p005.at(year).data());
  graph_nominal_BkgEff0p005.SetLineWidth(line_width);
  graph_nominal_BkgEff0p005.SetLineColor(kCyan);
  graph_nominal_BkgEff0p005.SetLineStyle(line_style_nominal);
  graph_nominal_BkgEff0p005.SetMarkerColor(kCyan);
  graph_nominal_BkgEff0p005.SetMarkerSize(marker_size);
  graph_nominal_BkgEff0p005.SetMarkerStyle(21);
  mg->Add(&graph_nominal_BkgEff0p005);

  cout << debug++ << endl;
  
  TGraph graph_nominal_BkgEff0p010 = TGraph(3, leo_ptCenter_W, values_leo_WvsQCD_BkgEff0p010.at(year).data());
  graph_nominal_BkgEff0p010.SetLineWidth(line_width);
  graph_nominal_BkgEff0p010.SetLineColor(kBlue);
  graph_nominal_BkgEff0p010.SetLineStyle(line_style_nominal);
  graph_nominal_BkgEff0p010.SetMarkerColor(kBlue);
  graph_nominal_BkgEff0p010.SetMarkerSize(marker_size);
  graph_nominal_BkgEff0p010.SetMarkerStyle(22);
  mg->Add(&graph_nominal_BkgEff0p010);

  cout << debug++ << endl;

  TGraph graph_nominal_BkgEff0p050 = TGraph(3, leo_ptCenter_W, values_leo_WvsQCD_BkgEff0p050.at(year).data());
  graph_nominal_BkgEff0p050.SetLineWidth(line_width);
  graph_nominal_BkgEff0p050.SetLineColor(kRed);
  graph_nominal_BkgEff0p050.SetLineStyle(line_style_nominal);
  graph_nominal_BkgEff0p050.SetMarkerColor(kRed);
  graph_nominal_BkgEff0p050.SetMarkerSize(marker_size);
  graph_nominal_BkgEff0p050.SetMarkerStyle(47);
  mg->Add(&graph_nominal_BkgEff0p050);

  cout << debug++ << endl;




  TGraph graph_btag_BkgEff0p005 = TGraph(3, leo_ptCenter_W, values_leo_WvsQCDMD_BkgEff0p005.at(year).data());
  graph_btag_BkgEff0p005.SetLineWidth(line_width);
  graph_btag_BkgEff0p005.SetLineColor(kCyan);
  graph_btag_BkgEff0p005.SetLineStyle(line_style_btag);
  graph_btag_BkgEff0p005.SetMarkerColor(kCyan);
  graph_btag_BkgEff0p005.SetMarkerSize(marker_size);
  graph_btag_BkgEff0p005.SetMarkerStyle(25);
  mg->Add(&graph_btag_BkgEff0p005);

  cout << debug++ << endl;
  
  TGraph graph_btag_BkgEff0p010 = TGraph(3, leo_ptCenter_W, values_leo_WvsQCDMD_BkgEff0p010.at(year).data());
  graph_btag_BkgEff0p010.SetLineWidth(line_width);
  graph_btag_BkgEff0p010.SetLineColor(kBlue);
  graph_btag_BkgEff0p010.SetLineStyle(line_style_btag);
  graph_btag_BkgEff0p010.SetMarkerColor(kBlue);
  graph_btag_BkgEff0p010.SetMarkerSize(marker_size);
  graph_btag_BkgEff0p010.SetMarkerStyle(26);
  mg->Add(&graph_btag_BkgEff0p010);

  cout << debug++ << endl;

  TGraph graph_btag_BkgEff0p025 = TGraph(3, leo_ptCenter_W, values_leo_WvsQCDMD_BkgEff0p025.at(year).data());
  graph_btag_BkgEff0p025.SetLineWidth(line_width);
  graph_btag_BkgEff0p025.SetLineColor(kMagenta);
  graph_btag_BkgEff0p025.SetLineStyle(line_style_btag);
  graph_btag_BkgEff0p025.SetMarkerColor(kMagenta);
  graph_btag_BkgEff0p025.SetMarkerSize(marker_size);
  graph_btag_BkgEff0p025.SetMarkerStyle(32);
  mg->Add(&graph_btag_BkgEff0p025);

  cout << debug++ << endl;


  // Canvas

  TCanvas *c = new TCanvas("canvas", "canvas title", 600, 600);
  c->cd();
  float margin_l = 0.15;
  float margin_r = 0.05;
  float margin_b = 0.12;
  float margin_t = 0.08;
  float tick_length = 0.015; // fraction of canvas width/height
  c->SetMargin(margin_l, margin_r, margin_b, margin_t); //lrbt
  auto coord = new CoordinateConverter();
  coord->init(margin_l, margin_r, margin_b, margin_t);
  // if(log_y) c->SetLogy();
  c->SetTickx(1);
  c->SetTicky(1);


  bool plot_at_bottom = true;


  mg->Draw("alp");
  mg->SetTitle("");
  // mg->GetHistogram()->GetXaxis()->SetRangeUser(300., 1000.);
  mg->GetHistogram()->GetXaxis()->SetLimits(200., 1200.);
  double maximum_saved = mg->GetHistogram()->GetMaximum();
  double minimum_saved = mg->GetHistogram()->GetMinimum();
  double new_maximum = maximum_saved * 1.618;
  // mg->GetHistogram()->SetMaximum(1.);
  // mg->GetHistogram()->SetMaximum(new_maximum);
  mg->GetHistogram()->SetMaximum(1.4);
  mg->GetHistogram()->SetMinimum(0.);
  // if(is_eff_bkg || is_eff_sig) mg->GetHistogram()->GetXaxis()->SetTitle(tagger.scan_var_tlatex_axis.c_str());
  // else if(is_roc) mg->GetHistogram()->GetXaxis()->SetTitle("True positive rate, #varepsilon_{S}");
  // if(is_eff_bkg || is_roc) mg->GetHistogram()->GetYaxis()->SetTitle("False positive rate, #varepsilon_{B}");
  // if(is_eff_sig) mg->GetHistogram()->GetYaxis()->SetTitle("True positive rate, #varepsilon_{S}");
  mg->GetHistogram()->GetXaxis()->SetTitle("Probe jet #it{p}_{T} [GeV]");
  mg->GetHistogram()->GetYaxis()->SetTitle("Signal efficiency");
  mg->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
  mg->GetHistogram()->GetXaxis()->CenterTitle();
  mg->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);
  mg->GetHistogram()->GetYaxis()->CenterTitle();
  // mg->GetHistogram()->GetXaxis()->SetNdivisions(405);
  // if(!log_y) mg->GetHistogram()->GetYaxis()->SetNdivisions(405);
  gStyle->SetLineWidth(1);

  // legend->Draw();





  // Lines
  // TLine(Double_t x1,Double_t y1,Double_t x2,Double_t y2)

  int lineX_color = kGray+2;
  int lineX_style = 3;
  int lineX_width = 1;

  // double line_max_y = maximum_saved * 1.05;
  double line_max_y = 1.;

  TLine* line2 = new TLine(300, 0, 300, line_max_y);
  line2->SetLineColor(lineX_color);
  line2->SetLineStyle(lineX_style);
  line2->SetLineWidth(lineX_width);
  line2->Draw();
  TLine* line5 = new TLine(400, 0, 400, line_max_y);
  line5->SetLineColor(lineX_color);
  line5->SetLineStyle(lineX_style);
  line5->SetLineWidth(lineX_width);
  line5->Draw();
  TLine* line6 = new TLine(800, 0, 800, line_max_y);
  line6->SetLineColor(lineX_color);
  line6->SetLineStyle(lineX_style);
  line6->SetLineWidth(lineX_width);
  line6->Draw();



  // Legend

  float leg_x1 = 0.62;
  float leg_y1 = 0.27;
  float leg_x2 = 0.83;
  float leg_y2 = 0.57;

  // TLegend *legend = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2, '', 'brNDC');
  TLegend *legend = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2);

  legend->AddEntry(&graph_nominal_BkgEff0p050, "Mistag rate = 5.0%", "p");
  legend->AddEntry(&graph_btag_BkgEff0p025, "Mistag rate = 2.5%", "p");
  legend->AddEntry(&graph_nominal_BkgEff0p010, "Mistag rate = 1.0%", "p");
  legend->AddEntry(&graph_nominal_BkgEff0p005, "Mistag rate = 0.5%", "p");

  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0); // Set background to transparent

  legend->Draw();


  // Legend nominal / subjet b tag

  float leg2_x1 = 0.44;
  float leg2_y1 = 0.72;
  float leg2_x2 = 0.77;
  float leg2_y2 = 0.87;

  // TLegend *legend2 = new TLegend(leg2_x1, leg2_y1, leg2_x2, leg2_y2, '', 'brNDC');
  TLegend *legend2 = new TLegend(leg2_x1, leg2_y1, leg2_x2, leg2_y2);

  TGraph *graph_nominal_dummy = new TGraph();
  graph_nominal_dummy->SetLineWidth(line_width);
  graph_nominal_dummy->SetLineColor(kBlack);
  graph_nominal_dummy->SetLineStyle(line_style_nominal);
  graph_nominal_dummy->SetMarkerColor(kBlack);
  graph_nominal_dummy->SetMarkerSize(marker_size);
  graph_nominal_dummy->SetMarkerStyle(20);

  TGraph *graph_btag_dummy = new TGraph();
  graph_btag_dummy->SetLineWidth(line_width);
  graph_btag_dummy->SetLineColor(kBlack);
  graph_btag_dummy->SetLineStyle(line_style_btag);
  graph_btag_dummy->SetMarkerColor(kBlack);
  graph_btag_dummy->SetMarkerSize(marker_size);
  graph_btag_dummy->SetMarkerStyle(24);

  legend2->SetHeader("#bf{W tagging using ParticleNet}");
  legend2->AddEntry(graph_nominal_dummy, "nominal: WvsQCD");
  legend2->AddEntry(graph_btag_dummy, "mass-decorrelated: WvsQCD-MD");
  
  legend2->SetTextSize(0.03);
  legend2->SetBorderSize(0);
  legend2->SetFillStyle(0); // Set background to transparent

  legend2->Draw();



  // TLatex *text_top_left = new TLatex(margin_l, 1-(margin_t-0.01), (string("Ultra Legacy ")+kYears.at(year).nice_name).c_str());
  // text_top_left->SetTextAlign(11); // left bottom aligned
  // text_top_left->SetTextFont(42);
  // text_top_left->SetTextSize(0.035);
  // text_top_left->SetNDC();
  // text_top_left->Draw();

  // TLatex *algo_label = new TLatex(coord->ConvertGraphXToPadX(1-(0.05*c->GetWh()/c->GetWw())), coord->ConvertGraphYToPadY(0.95), "AK8 PUPPI");
  // algo_label->SetTextAlign(33); // right top
  // algo_label->SetTextFont(62);
  // algo_label->SetTextSize(0.05); // 0.05
  // algo_label->SetTextColor(kGray+1);
  // algo_label->SetNDC();
  // algo_label->Draw();

  // const double eta_pt_text_x = coord->ConvertGraphXToPadX(1-(0.05*c->GetWh()/c->GetWw()));
  // double eta_text_y = coord->ConvertGraphYToPadY(0.801);
  // // double pt_text_y = coord->ConvertGraphYToPadY(0.867);
  // if(!plot_at_bottom) {
  //   eta_text_y = coord->ConvertGraphYToPadY((0.05*c->GetWh()/c->GetWw()) + 0.066);
  //   // pt_text_y = coord->ConvertGraphYToPadY((0.05*c->GetWh()/c->GetWw()));
  // }

  // TLatex *eta_text = new TLatex(eta_pt_text_x, eta_text_y, "|#eta| < 2.5");
  // eta_text->SetTextAlign(plot_at_bottom ? 33 : 31); // right top
  // eta_text->SetTextFont(42);
  // eta_text->SetTextSize(0.035);
  // eta_text->SetNDC();
  // eta_text->Draw();

  // TLatex *pt_text = new TLatex(eta_pt_text_x, pt_text_y, pt_bin.text.c_str());
  // pt_text->SetTextAlign(plot_at_bottom ? 33 : 31); // right top
  // pt_text->SetTextFont(42);
  // pt_text->SetTextSize(0.035);
  // pt_text->SetNDC();
  // pt_text->Draw();

  // string string_text_top_right = "t#bar{t} #rightarrow all-jets @ #sqrt{#it{s}} = 13 TeV";
  // string string_text_top_right = "#sqrt{#it{s}} = 13 TeV";
  string string_text_top_right = kYears.at(year).lumi_fb_display+" fb^{#minus1} (13 TeV)";
  TLatex *text_top_right = new TLatex(1-margin_r, 1-(margin_t-0.01), string_text_top_right.c_str());
  text_top_right->SetTextAlign(31); // right bottom aligned
  text_top_right->SetTextFont(42);
  text_top_right->SetTextSize(0.035);
  text_top_right->SetNDC();
  text_top_right->Draw();

  TText *prelimPW1 = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.95), "Private work");
  prelimPW1->SetTextAlign(13); // left top
  prelimPW1->SetTextFont(52);
  prelimPW1->SetTextSize(0.03);
  prelimPW1->SetNDC();
  // prelimPW1->SetTextColor(kWhite);
  if(!bool_prelim) prelimPW1->Draw();

  TText *prelimPW2 = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.90), "(CMS data/simulation)");
  prelimPW2->SetTextAlign(13); // left top
  prelimPW2->SetTextFont(52);
  prelimPW2->SetTextSize(0.03);
  prelimPW2->SetNDC();
  // prelimPW2->SetTextColor(kWhite);
  if(!bool_prelim) prelimPW2->Draw();

  TLatex *cms = new TLatex(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.95), "CMS");
  cms->SetTextAlign(13); // left top
  cms->SetTextFont(62);
  cms->SetTextSize(0.05);
  cms->SetNDC();
  if(bool_prelim) cms->Draw();

  if(!plot_at_bottom) {
    TString prelim_text;
    if(bool_prelim) prelim_text = "Simulation, Preliminary";
    else prelim_text = "Simulation, Private Work";
    TText *prelim = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.87), prelim_text);
    prelim->SetTextAlign(13); // left top
    prelim->SetTextFont(52);
    prelim->SetTextSize(0.035);
    prelim->SetNDC();
    if(bool_prelim) prelim->Draw();
  }
  else {
    TText *prelim = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.87), "Simulation");
    prelim->SetTextAlign(13); // left top
    prelim->SetTextFont(52);
    prelim->SetTextSize(0.035);
    prelim->SetNDC();
    // prelim->SetTextColor(kWhite);
    if(bool_prelim) prelim->Draw();

    TString prelim_text;
    if(bool_prelim) prelim_text = "Preliminary";
    else prelim_text = "Private Work";
    TText *prelim2 = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.804), prelim_text);
    prelim2->SetTextAlign(13); // left top
    prelim2->SetTextFont(52);
    prelim2->SetTextSize(0.035);
    prelim2->SetNDC();
    // prelim2->SetTextColor(kWhite);
    if(bool_prelim) prelim2->Draw();
  }


  // https://root-forum.cern.ch/t/inconsistent-tick-length/18563/8
  c->Update();
  double tickScaleX = (c->GetUxmax() - c->GetUxmin()) / (c->GetX2() - c->GetX1()) * (c->GetWh() * c->GetAbsHNDC());
  double tickScaleY = (c->GetUymax() - c->GetUymin()) / (c->GetY2() - c->GetY1()) * (c->GetWw() * c->GetAbsWNDC());
  mg->GetHistogram()->GetYaxis()->SetTickLength(c->GetWw() * tick_length / tickScaleY);
  mg->GetHistogram()->GetXaxis()->SetTickLength(c->GetWh() * tick_length / tickScaleX);

  gPad->Update();
  gPad->RedrawAxis();

  // Save to disk
  string plotName;
  if(bool_prelim) plotName = "plot_sig_eff_in_mc__partnet_w__"+kYears.at(year).name+".pdf";
  else plotName = "plot_sig_eff_in_mc__partnet_w__"+kYears.at(year).name+"-PrivateWork.pdf";
  string plotPath = data_path+"/"+plotName;
  c->SaveAs(plotPath.c_str());
  delete c;
}

void do_plot_partnet_top(const Year & year, const bool bool_prelim) {

  int debug = 0;

  const int line_width = 2;
  const int marker_size = 1;
  const int line_style_nominal = 1;
  const int line_style_btag = 2;

  TMultiGraph *mg = new TMultiGraph("mg", "mg title");

  cout << debug++ << endl;

  TGraph graph_nominal_BkgEff0p001 = TGraph(4, leo_ptCenter_top, values_leo_TvsQCD_BkgEff0p001.at(year).data());
  graph_nominal_BkgEff0p001.SetLineWidth(line_width);
  graph_nominal_BkgEff0p001.SetLineColor(kGreen);
  graph_nominal_BkgEff0p001.SetLineStyle(line_style_nominal);
  graph_nominal_BkgEff0p001.SetMarkerColor(kGreen);
  graph_nominal_BkgEff0p001.SetMarkerSize(marker_size);
  graph_nominal_BkgEff0p001.SetMarkerStyle(20);
  mg->Add(&graph_nominal_BkgEff0p001);

  cout << debug++ << endl;
  
  TGraph graph_nominal_BkgEff0p005 = TGraph(4, leo_ptCenter_top, values_leo_TvsQCD_BkgEff0p005.at(year).data());
  graph_nominal_BkgEff0p005.SetLineWidth(line_width);
  graph_nominal_BkgEff0p005.SetLineColor(kCyan);
  graph_nominal_BkgEff0p005.SetLineStyle(line_style_nominal);
  graph_nominal_BkgEff0p005.SetMarkerColor(kCyan);
  graph_nominal_BkgEff0p005.SetMarkerSize(marker_size);
  graph_nominal_BkgEff0p005.SetMarkerStyle(21);
  mg->Add(&graph_nominal_BkgEff0p005);

  cout << debug++ << endl;

  TGraph graph_nominal_BkgEff0p010 = TGraph(4, leo_ptCenter_top, values_leo_TvsQCD_BkgEff0p010.at(year).data());
  graph_nominal_BkgEff0p010.SetLineWidth(line_width);
  graph_nominal_BkgEff0p010.SetLineColor(kBlue);
  graph_nominal_BkgEff0p010.SetLineStyle(line_style_nominal);
  graph_nominal_BkgEff0p010.SetMarkerColor(kBlue);
  graph_nominal_BkgEff0p010.SetMarkerSize(marker_size);
  graph_nominal_BkgEff0p010.SetMarkerStyle(22);
  mg->Add(&graph_nominal_BkgEff0p010);

  cout << debug++ << endl;





  // Canvas

  TCanvas *c = new TCanvas("canvas", "canvas title", 600, 600);
  c->cd();
  float margin_l = 0.15;
  float margin_r = 0.05;
  float margin_b = 0.12;
  float margin_t = 0.08;
  float tick_length = 0.015; // fraction of canvas width/height
  c->SetMargin(margin_l, margin_r, margin_b, margin_t); //lrbt
  auto coord = new CoordinateConverter();
  coord->init(margin_l, margin_r, margin_b, margin_t);
  // if(log_y) c->SetLogy();
  c->SetTickx(1);
  c->SetTicky(1);


  bool plot_at_bottom = true;


  mg->Draw("alp");
  mg->SetTitle("");
  // mg->GetHistogram()->GetXaxis()->SetRangeUser(300., 1000.);
  mg->GetHistogram()->GetXaxis()->SetLimits(200., 1200.);
  double maximum_saved = mg->GetHistogram()->GetMaximum();
  double minimum_saved = mg->GetHistogram()->GetMinimum();
  double new_maximum = maximum_saved * 1.618;
  // mg->GetHistogram()->SetMaximum(1.);
  // mg->GetHistogram()->SetMaximum(new_maximum);
  mg->GetHistogram()->SetMaximum(1.4);
  mg->GetHistogram()->SetMinimum(0.);
  // if(is_eff_bkg || is_eff_sig) mg->GetHistogram()->GetXaxis()->SetTitle(tagger.scan_var_tlatex_axis.c_str());
  // else if(is_roc) mg->GetHistogram()->GetXaxis()->SetTitle("True positive rate, #varepsilon_{S}");
  // if(is_eff_bkg || is_roc) mg->GetHistogram()->GetYaxis()->SetTitle("False positive rate, #varepsilon_{B}");
  // if(is_eff_sig) mg->GetHistogram()->GetYaxis()->SetTitle("True positive rate, #varepsilon_{S}");
  mg->GetHistogram()->GetXaxis()->SetTitle("Probe jet #it{p}_{T} [GeV]");
  mg->GetHistogram()->GetYaxis()->SetTitle("Signal efficiency");
  mg->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
  mg->GetHistogram()->GetXaxis()->CenterTitle();
  mg->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);
  mg->GetHistogram()->GetYaxis()->CenterTitle();
  // mg->GetHistogram()->GetXaxis()->SetNdivisions(405);
  // if(!log_y) mg->GetHistogram()->GetYaxis()->SetNdivisions(405);
  gStyle->SetLineWidth(1);

  // legend->Draw();





  // Lines
  // TLine(Double_t x1,Double_t y1,Double_t x2,Double_t y2)

  int lineX_color = kGray+2;
  int lineX_style = 3;
  int lineX_width = 1;

  // double line_max_y = maximum_saved * 1.05;
  double line_max_y = 1.;

  TLine* line2 = new TLine(300, 0, 300, line_max_y);
  line2->SetLineColor(lineX_color);
  line2->SetLineStyle(lineX_style);
  line2->SetLineWidth(lineX_width);
  line2->Draw();
  TLine* line4 = new TLine(400, 0, 400, line_max_y);
  line4->SetLineColor(lineX_color);
  line4->SetLineStyle(lineX_style);
  line4->SetLineWidth(lineX_width);
  line4->Draw();
  TLine* line5 = new TLine(480, 0, 480, line_max_y);
  line5->SetLineColor(lineX_color);
  line5->SetLineStyle(lineX_style);
  line5->SetLineWidth(lineX_width);
  line5->Draw();
  TLine* line6 = new TLine(600, 0, 600, line_max_y);
  line6->SetLineColor(lineX_color);
  line6->SetLineStyle(lineX_style);
  line6->SetLineWidth(lineX_width);
  line6->Draw();



  // Legend

  float leg_x1 = 0.62;
  float leg_y1 = 0.48;
  float leg_x2 = 0.83;
  float leg_y2 = 0.68;

  // TLegend *legend = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2, '', 'brNDC');
  TLegend *legend = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2);

  legend->AddEntry(&graph_nominal_BkgEff0p010, "Mistag rate = 1.0%", "p");
  legend->AddEntry(&graph_nominal_BkgEff0p005, "Mistag rate = 0.5%", "p");
  legend->AddEntry(&graph_nominal_BkgEff0p001, "Mistag rate = 0.1%", "p");

  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0); // Set background to transparent

  legend->Draw();


  // Legend nominal / subjet b tag

  float leg2_x1 = 0.44;
  float leg2_y1 = 0.72;
  float leg2_x2 = 0.77;
  float leg2_y2 = 0.87;

  // TLegend *legend2 = new TLegend(leg2_x1, leg2_y1, leg2_x2, leg2_y2, '', 'brNDC');
  TLegend *legend2 = new TLegend(leg2_x1, leg2_y1, leg2_x2, leg2_y2);

  TGraph *graph_nominal_dummy = new TGraph();
  graph_nominal_dummy->SetLineWidth(line_width);
  graph_nominal_dummy->SetLineColor(kBlack);
  graph_nominal_dummy->SetLineStyle(line_style_nominal);
  graph_nominal_dummy->SetMarkerColor(kBlack);
  graph_nominal_dummy->SetMarkerSize(marker_size);
  graph_nominal_dummy->SetMarkerStyle(20);

  TGraph *graph_btag_dummy = new TGraph();
  graph_btag_dummy->SetLineWidth(line_width);
  graph_btag_dummy->SetLineColor(kBlack);
  graph_btag_dummy->SetLineStyle(line_style_btag);
  graph_btag_dummy->SetMarkerColor(kBlack);
  graph_btag_dummy->SetMarkerSize(marker_size);
  graph_btag_dummy->SetMarkerStyle(24);

  legend2->SetHeader("#bf{t tagging using ParticleNet TvsQCD}");
  legend2->AddEntry(graph_nominal_dummy, "", "");
  legend2->AddEntry(graph_btag_dummy, "", "");
  
  legend2->SetTextSize(0.03);
  legend2->SetBorderSize(0);
  legend2->SetFillStyle(0); // Set background to transparent

  legend2->Draw();



  // TLatex *text_top_left = new TLatex(margin_l, 1-(margin_t-0.01), (string("Ultra Legacy ")+kYears.at(year).nice_name).c_str());
  // text_top_left->SetTextAlign(11); // left bottom aligned
  // text_top_left->SetTextFont(42);
  // text_top_left->SetTextSize(0.035);
  // text_top_left->SetNDC();
  // text_top_left->Draw();

  // TLatex *algo_label = new TLatex(coord->ConvertGraphXToPadX(1-(0.05*c->GetWh()/c->GetWw())), coord->ConvertGraphYToPadY(0.95), "AK8 PUPPI");
  // algo_label->SetTextAlign(33); // right top
  // algo_label->SetTextFont(62);
  // algo_label->SetTextSize(0.05); // 0.05
  // algo_label->SetTextColor(kGray+1);
  // algo_label->SetNDC();
  // algo_label->Draw();

  // const double eta_pt_text_x = coord->ConvertGraphXToPadX(1-(0.05*c->GetWh()/c->GetWw()));
  // double eta_text_y = coord->ConvertGraphYToPadY(0.801);
  // // double pt_text_y = coord->ConvertGraphYToPadY(0.867);
  // if(!plot_at_bottom) {
  //   eta_text_y = coord->ConvertGraphYToPadY((0.05*c->GetWh()/c->GetWw()) + 0.066);
  //   // pt_text_y = coord->ConvertGraphYToPadY((0.05*c->GetWh()/c->GetWw()));
  // }

  // TLatex *eta_text = new TLatex(eta_pt_text_x, eta_text_y, "|#eta| < 2.5");
  // eta_text->SetTextAlign(plot_at_bottom ? 33 : 31); // right top
  // eta_text->SetTextFont(42);
  // eta_text->SetTextSize(0.035);
  // eta_text->SetNDC();
  // eta_text->Draw();

  // TLatex *pt_text = new TLatex(eta_pt_text_x, pt_text_y, pt_bin.text.c_str());
  // pt_text->SetTextAlign(plot_at_bottom ? 33 : 31); // right top
  // pt_text->SetTextFont(42);
  // pt_text->SetTextSize(0.035);
  // pt_text->SetNDC();
  // pt_text->Draw();

  // string string_text_top_right = "t#bar{t} #rightarrow all-jets @ #sqrt{#it{s}} = 13 TeV";
  // string string_text_top_right = "#sqrt{#it{s}} = 13 TeV";
  string string_text_top_right = kYears.at(year).lumi_fb_display+" fb^{#minus1} (13 TeV)";
  TLatex *text_top_right = new TLatex(1-margin_r, 1-(margin_t-0.01), string_text_top_right.c_str());
  text_top_right->SetTextAlign(31); // right bottom aligned
  text_top_right->SetTextFont(42);
  text_top_right->SetTextSize(0.035);
  text_top_right->SetNDC();
  text_top_right->Draw();

  TText *prelimPW1 = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.95), "Private work");
  prelimPW1->SetTextAlign(13); // left top
  prelimPW1->SetTextFont(52);
  prelimPW1->SetTextSize(0.03);
  prelimPW1->SetNDC();
  // prelimPW1->SetTextColor(kWhite);
  if(!bool_prelim) prelimPW1->Draw();

  TText *prelimPW2 = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.90), "(CMS data/simulation)");
  prelimPW2->SetTextAlign(13); // left top
  prelimPW2->SetTextFont(52);
  prelimPW2->SetTextSize(0.03);
  prelimPW2->SetNDC();
  // prelimPW2->SetTextColor(kWhite);
  if(!bool_prelim) prelimPW2->Draw();

  TLatex *cms = new TLatex(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.95), "CMS");
  cms->SetTextAlign(13); // left top
  cms->SetTextFont(62);
  cms->SetTextSize(0.05);
  cms->SetNDC();
  if(bool_prelim) cms->Draw();

  if(!plot_at_bottom) {
    TString prelim_text;
    if(bool_prelim) prelim_text = "Simulation, Preliminary";
    else prelim_text = "Simulation, Private Work";
    TText *prelim = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.87), prelim_text);
    prelim->SetTextAlign(13); // left top
    prelim->SetTextFont(52);
    prelim->SetTextSize(0.035);
    prelim->SetNDC();
    if(bool_prelim) prelim->Draw();
  }
  else {
    TText *prelim = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.87), "Simulation");
    prelim->SetTextAlign(13); // left top
    prelim->SetTextFont(52);
    prelim->SetTextSize(0.035);
    prelim->SetNDC();
    // prelim->SetTextColor(kWhite);
    if(bool_prelim) prelim->Draw();

    TString prelim_text;
    if(bool_prelim) prelim_text = "Preliminary";
    else prelim_text = "Private Work";
    TText *prelim2 = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.804), prelim_text);
    prelim2->SetTextAlign(13); // left top
    prelim2->SetTextFont(52);
    prelim2->SetTextSize(0.035);
    prelim2->SetNDC();
    // prelim2->SetTextColor(kWhite);
    if(bool_prelim) prelim2->Draw();
  }


  // https://root-forum.cern.ch/t/inconsistent-tick-length/18563/8
  c->Update();
  double tickScaleX = (c->GetUxmax() - c->GetUxmin()) / (c->GetX2() - c->GetX1()) * (c->GetWh() * c->GetAbsHNDC());
  double tickScaleY = (c->GetUymax() - c->GetUymin()) / (c->GetY2() - c->GetY1()) * (c->GetWw() * c->GetAbsWNDC());
  mg->GetHistogram()->GetYaxis()->SetTickLength(c->GetWw() * tick_length / tickScaleY);
  mg->GetHistogram()->GetXaxis()->SetTickLength(c->GetWh() * tick_length / tickScaleX);

  gPad->Update();
  gPad->RedrawAxis();

  // Save to disk
  string plotName;
  if(bool_prelim) plotName = "plot_sig_eff_in_mc__partnet_top__"+kYears.at(year).name+".pdf";
  else plotName = "plot_sig_eff_in_mc__partnet_top__"+kYears.at(year).name+"-PrivateWork.pdf";
  string plotPath = data_path+"/"+plotName;
  c->SaveAs(plotPath.c_str());
  delete c;
}



void sig_eff_in_mc_plotting() {

    for(const auto year : kAllYears) {
      do_plot_nsubjettiness(year, true);
      do_plot_nsubjettiness(year, false);
      do_plot_partnet_w(year, true);
      do_plot_partnet_w(year, false);
      do_plot_partnet_top(year, true);
      do_plot_partnet_top(year, false);
  }

}