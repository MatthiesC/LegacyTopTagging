#include "../constants.h"

using namespace macros;

// typedef struct {
//   string name;
//   string long_name;
//   // string lumi_fb;
//   Int_t color;
// } Year;

typedef struct {
  string name;
  string watermark;
} Process;

class CoordinateConverter {
public:
  CoordinateConverter() {};
  void init(const float l, const float r, const float b, const float t);
  float ConvertGraphXToPadX(const float graph_x);
  float ConvertGraphYToPadY(const float graph_y);
private:
  float pad_l, pad_r, pad_b, pad_t; // pad margins
  float graph_width, graph_height;
};

void CoordinateConverter::init(const float l, const float r, const float b, const float t) {
  pad_l = l;
  pad_r = r;
  pad_b = b;
  pad_t = t;
  graph_width = 1.-l-r;
  graph_height = 1.-b-t;
}

float CoordinateConverter::ConvertGraphXToPadX(const float graph_x) {
  return pad_l+graph_x*graph_width;
}

float CoordinateConverter::ConvertGraphYToPadY(const float graph_y) {
  return pad_b+graph_y*graph_height;
}

enum class EVariation {
  nominal,
  jec_up,
  jec_down,
  jer_up,
  jer_down,
};

const map<EVariation, string> kVariationToString = {
  {EVariation::nominal, "nominal"},
  {EVariation::jec_up, "syst_jec_up"},
  {EVariation::jec_down, "syst_jec_down"},
  {EVariation::jer_up, "syst_jer_up"},
  {EVariation::jer_down, "syst_jer_down"},
};


void do_plot(const Year & year, const ProbeJetAlgo & probejetalgo, const string & dr_max, const bool log_y) {

  // initial binning of response histograms: 2000 bins between 0 and 2000 (--> 1 GeV bins)
  vector<float> v_binning;
  // for(uint i = 0; i < 100/100; i++) v_binning.push_back(i*100+0);
  // for(uint i = 0; i < 700/10; i++) v_binning.push_back(i*10+100);
  // for(uint i = 0; i < 300/20; i++) v_binning.push_back(i*20+800);
  // for(uint i = 0; i < 400/50; i++) v_binning.push_back(i*50+1100);
  // for(uint i = 0; i <= 500/100; i++) v_binning.push_back(i*100+1500);
  for(uint i = 0; i < 100/100; i++) v_binning.push_back(i*100+0);
  for(uint i = 0; i < 90/10; i++) v_binning.push_back(i*10+100);
  for(uint i = 0; i < 20/20; i++) v_binning.push_back(i*20+190); // this gives a point with a bin center of 200
  for(uint i = 0; i < 190/10; i++) v_binning.push_back(i*10+210);
  for(uint i = 0; i < 300/15; i++) v_binning.push_back(i*15+400);
  for(uint i = 0; i < 80/20; i++) v_binning.push_back(i*20+700);
  for(uint i = 0; i < 120/40; i++) v_binning.push_back(i*40+780);
  for(uint i = 0; i < 150/75; i++) v_binning.push_back(i*75+900);
  for(uint i = 0; i < 200/100; i++) v_binning.push_back(i*100+1050); // this gives a point with a bin center of 1200
  for(uint i = 0; i <= 750/250; i++) v_binning.push_back(i*250+1250);
  // for(float f : v_binning) cout << f << endl;
  Double_t a_binning[v_binning.size()];
  copy(v_binning.begin(), v_binning.end(), a_binning);

  map<EVariation, TH1F*> hists_corr;
  map<EVariation, TH1F*> hists_raw;
  for(const auto & var : kVariationToString) {
    const string infileBasePath = (string)getenv("CMSSW_BASE")+"/src/UHH2/LegacyTopTagging/output/WorkingPointStudy/"+kYears.at(year).name+"/"+var.second+"/";
    const string infileName = (string)"uhh2.AnalysisModuleRunner.MC.TTbarToHadronic_"+kYears.at(year).name+".root";
    const string infilePath = infileBasePath+infileName;
    TFile *infile = TFile::Open(infilePath.c_str(), "READ");
    const string hist_collection = kProbeJetAlgos.at(probejetalgo).name+"Hists_1_after_corrections";
    TH1F* inhist_gen = (TH1F*)infile->Get((hist_collection+"/response_gen_"+dr_max).c_str());
    TH1F* inhist_corr = (TH1F*)infile->Get((hist_collection+"/response_corr_"+dr_max).c_str());
    TH1F* inhist_raw = (TH1F*)infile->Get((hist_collection+"/response_raw_"+dr_max).c_str());
    TH1F* inhist_rebin_gen = dynamic_cast<TH1F*>(inhist_gen->Rebin(v_binning.size()-1, inhist_gen->GetName(), a_binning));
    TH1F* inhist_rebin_corr = dynamic_cast<TH1F*>(inhist_corr->Rebin(v_binning.size()-1, inhist_corr->GetName(), a_binning));
    TH1F* inhist_rebin_raw = dynamic_cast<TH1F*>(inhist_raw->Rebin(v_binning.size()-1, inhist_raw->GetName(), a_binning));
    TH1F *hist_gen = new TH1F(*inhist_rebin_gen);
    TH1F *hist_corr = new TH1F(*inhist_rebin_corr);
    TH1F *hist_raw = new TH1F(*inhist_rebin_raw);
    for(unsigned int ibin = 0; ibin <= hist_gen->GetNbinsX()+1; ++ibin) {
      hist_corr->SetBinContent(ibin, hist_corr->GetBinContent(ibin) / hist_gen->GetBinContent(ibin));
      hist_corr->SetBinError(ibin, hist_corr->GetBinError(ibin) / hist_gen->GetBinContent(ibin));
      hist_raw->SetBinContent(ibin, hist_raw->GetBinContent(ibin) / hist_gen->GetBinContent(ibin));
      hist_raw->SetBinError(ibin, hist_raw->GetBinError(ibin) / hist_gen->GetBinContent(ibin));
    }
    hists_corr[var.first] = hist_corr;
    hists_raw[var.first] = hist_raw;
  }

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
  if(log_y) c->SetLogy();
  c->SetTickx(1);
  c->SetTicky(1);

  const double pt_min_in_plot = 200.;
  const double pt_max_in_plot = 1300.;
  vector<double> values_zeroes;
  vector<double> values_ptgen;
  vector<double> values_raw_central;
  vector<double> values_raw_mc_unc;
  vector<double> values_corr_central;
  vector<double> values_corr_mc_unc; // symmetric error, thus no differentiation between up/down
  vector<double> values_corr_jec_plus;
  vector<double> values_corr_jec_minus;
  vector<double> values_corr_jer_plus;
  vector<double> values_corr_jer_minus;
  vector<double> values_corr_tot_plus;
  vector<double> values_corr_tot_minus;
  const TH1F *dummy_hist = hists_corr.at(EVariation::nominal); // used to get the binnning information; for all hists, the binning is and needs to be the same in order for this to work
  for(unsigned int ibin = 1; ibin <= dummy_hist->GetNbinsX(); ++ibin) { // skip overflow and underflow bins
    const double central_ptgen = dummy_hist->GetXaxis()->GetBinCenter(ibin);
    if(central_ptgen < pt_min_in_plot) continue;
    if(central_ptgen >= pt_max_in_plot) continue;
    values_zeroes.push_back(0.);
    values_ptgen.push_back(central_ptgen);
    // raw:
    values_raw_central.push_back(hists_raw.at(EVariation::nominal)->GetBinContent(ibin));
    values_raw_mc_unc.push_back(hists_raw.at(EVariation::nominal)->GetBinError(ibin));
    // corrected:
    const double value_central = hists_corr.at(EVariation::nominal)->GetBinContent(ibin);
    values_corr_central.push_back(value_central);
    const double value_mc_unc = hists_corr.at(EVariation::nominal)->GetBinError(ibin);
    values_corr_mc_unc.push_back(value_mc_unc);
    const double value_jec_up = hists_corr.at(EVariation::jec_up)->GetBinContent(ibin);
    const double value_jec_down = hists_corr.at(EVariation::jec_down)->GetBinContent(ibin);
    const double value_jec_plus = max(value_central, max(value_jec_up, value_jec_down)) - value_central;
    const double value_jec_minus = value_central - min(value_central, min(value_jec_up, value_jec_down));
    values_corr_jec_plus.push_back(value_jec_plus);
    values_corr_jec_minus.push_back(value_jec_minus);
    const double value_jer_up = hists_corr.at(EVariation::jer_up)->GetBinContent(ibin);
    const double value_jer_down = hists_corr.at(EVariation::jer_down)->GetBinContent(ibin);
    const double value_jer_plus = max(value_central, max(value_jer_up, value_jer_down)) - value_central;
    const double value_jer_minus = value_central - min(value_central, min(value_jer_up, value_jer_down));
    values_corr_jer_plus.push_back(value_jer_plus);
    values_corr_jer_minus.push_back(value_jer_minus);
    const double value_tot_plus = sqrt(value_mc_unc*value_mc_unc + value_jec_plus*value_jec_plus + value_jer_plus*value_jer_plus);
    const double value_tot_minus = sqrt(value_mc_unc*value_mc_unc + value_jec_minus*value_jec_minus + value_jer_minus*value_jer_minus);
    values_corr_tot_plus.push_back(value_tot_plus);
    values_corr_tot_minus.push_back(value_tot_minus);
  }

  TGraphAsymmErrors *graph_raw_mc = new TGraphAsymmErrors(values_ptgen.size(), &values_ptgen[0], &values_raw_central[0], &values_zeroes[0], &values_zeroes[0], &values_raw_mc_unc[0], &values_raw_mc_unc[0]);
  TGraphAsymmErrors *graph_raw_line = new TGraphAsymmErrors(*graph_raw_mc);
  TGraphAsymmErrors *graph_corr_mc = new TGraphAsymmErrors(values_ptgen.size(), &values_ptgen[0], &values_corr_central[0], &values_zeroes[0], &values_zeroes[0], &values_corr_mc_unc[0], &values_corr_mc_unc[0]);
  TGraphAsymmErrors *graph_corr_tot = new TGraphAsymmErrors(values_ptgen.size(), &values_ptgen[0], &values_corr_central[0], &values_zeroes[0], &values_zeroes[0], &values_corr_tot_minus[0], &values_corr_tot_plus[0]);
  TGraphAsymmErrors *graph_corr_line = new TGraphAsymmErrors(*graph_corr_tot);


  Int_t color_raw = kYears.at(year).tcolor;
  Int_t color_corr = kYears.at(year).tcolor;
  const double transparency1 = 1.0;
  const double transparency2 = 0.5;

  graph_raw_line->SetFillColorAlpha(0, 0);
  graph_raw_line->SetLineStyle(kDashed);
  graph_raw_line->SetLineColor(kBlack);
  graph_raw_line->SetLineWidth(2);

  graph_corr_line->SetFillColorAlpha(0, 0);
  graph_corr_line->SetLineStyle(kSolid);
  graph_corr_line->SetLineColor(kBlack);
  graph_corr_line->SetLineWidth(2);

  graph_raw_mc->SetFillColorAlpha(color_raw, transparency1);
  // graph_raw_mc->SetFillStyle(3357);
  graph_raw_mc->SetLineWidth(0);
  graph_corr_mc->SetFillColorAlpha(color_corr, transparency1);
  // graph_corr_mc->SetFillStyle(3357);
  graph_corr_mc->SetLineWidth(0);
  graph_corr_tot->SetFillColorAlpha(color_corr, transparency2);
  graph_corr_tot->SetLineWidth(0);

  // do this after all settings for graph_raw_mc are fixed!
  TGraphAsymmErrors *graph_raw_mc_copy = new TGraphAsymmErrors(*graph_raw_mc);
  graph_raw_mc_copy->SetFillColorAlpha(color_raw, transparency2);

  TMultiGraph *mg = new TMultiGraph("mg", "mg title");

  mg->Add(graph_corr_tot);
  mg->Add(graph_corr_mc);
  mg->Add(graph_corr_line);
  mg->Add(graph_raw_mc);
  mg->Add(graph_raw_mc_copy); // plot this graph twice to make the transparency for MC unc on raw plot same as MC unc on corr plot which is overlapping with tot unc
  mg->Add(graph_raw_line);

  mg->Draw("AL3");
  mg->SetTitle("");
  mg->GetHistogram()->GetXaxis()->SetTitle("#it{p}_{T}^{gen} [GeV]");
  mg->GetHistogram()->GetYaxis()->SetTitle("#LT #it{p}_{T}^{rec} / #it{p}_{T}^{gen} #GT");
  mg->GetHistogram()->GetXaxis()->CenterTitle();
  mg->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
  mg->GetHistogram()->GetYaxis()->CenterTitle();
  mg->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);


  double maximum_saved = mg->GetHistogram()->GetMaximum();
  double minimum_saved = mg->GetHistogram()->GetMinimum();
  mg->GetHistogram()->SetMaximum(1.15);
  mg->GetHistogram()->SetMinimum(0.85);

  // // graph-dependent attributes, order: raw, msd, msd_btag
  // vector<int> linewidths = {2, 2, 2};
  // vector<int> linestyles = {3, 1, 2};
  // vector<int> linecolors = {kBlue, kRed, kMagenta};
  // vector<string> legends = {"No #it{m}_{SD} window nor b tagging", "105 < #it{m}_{SD} [GeV] < 210", "#it{m}_{SD} as above + subjet b tag"};
  // int iterator = -1;
  // for(auto graph : graphs) {
  //   // min0max1(graph);
  //   mg->Add(graph);
  //   graph->SetLineWidth(linewidths.at(++iterator));
  //   graph->SetLineStyle(linestyles.at(iterator));
  //   // graph->SetLineColor(linecolors.at(iterator));
  //   legend->AddEntry(graph, legends.at(iterator).c_str(), "l");
  // }
  //
  // mg->Draw("ac");
  // mg->SetTitle("");
  // mg->GetHistogram()->GetXaxis()->SetRangeUser(0., 1.);
  // double maximum_saved = mg->GetHistogram()->GetMaximum();
  // double minimum_saved = mg->GetHistogram()->GetMinimum();
  // mg->GetHistogram()->SetMaximum(1.);
  // if(!log_y) mg->GetHistogram()->SetMinimum(0.);
  // else mg->GetHistogram()->SetMinimum(0.0002);
  // if(is_eff_qcd || is_eff_ttbar) mg->GetHistogram()->GetXaxis()->SetTitle("#tau_{3}/#tau_{2} upper limit");
  // else if(is_roc) mg->GetHistogram()->GetXaxis()->SetTitle("#varepsilon_{S}");
  // if(is_eff_qcd || is_roc) mg->GetHistogram()->GetYaxis()->SetTitle("#varepsilon_{B}");
  // if(is_eff_ttbar) mg->GetHistogram()->GetYaxis()->SetTitle("#varepsilon_{S}");
  // mg->GetHistogram()->GetXaxis()->SetTitleOffset(1.2);
  // mg->GetHistogram()->GetYaxis()->SetTitleOffset(1.2);
  // // mg->GetHistogram()->GetXaxis()->SetNdivisions(405);
  // if(!log_y) mg->GetHistogram()->GetYaxis()->SetNdivisions(405);
  // gStyle->SetLineWidth(1);
  //
  // legend->Draw();
  //



  TLatex *legend1 = new TLatex(coord->ConvertGraphXToPadX(0.07), coord->ConvertGraphYToPadY(0.655), "JEC applied, with MC stat. (inner) #oplus JEC (outer) uncertainties:");
  legend1->SetTextAlign(13); // left top
  legend1->SetTextFont(52);
  legend1->SetTextSize(0.025);
  legend1->SetNDC();
  legend1->Draw();

  TLatex *legend2 = new TLatex(coord->ConvertGraphXToPadX(0.07), coord->ConvertGraphYToPadY(0.325), "Uncorrected, with MC stat. uncertainty:");
  legend2->SetTextAlign(13); // left top
  legend2->SetTextFont(52);
  legend2->SetTextSize(0.025);
  legend2->SetNDC();
  legend2->Draw();



  TLatex *text_top_left = new TLatex(margin_l, 1-(margin_t-0.01), (string("Ultra Legacy ")+kYears.at(year).nice_name).c_str());
  text_top_left->SetTextAlign(11); // left bottom aligned
  text_top_left->SetTextFont(42);
  text_top_left->SetTextSize(0.035);
  text_top_left->SetNDC();
  text_top_left->Draw();

  TLatex *algo_label = new TLatex(coord->ConvertGraphXToPadX(1-(0.05*c->GetWh()/c->GetWw())), coord->ConvertGraphYToPadY(0.95), (kProbeJetAlgos.at(probejetalgo).name+" PUPPI").c_str());
  algo_label->SetTextAlign(33); // right top
  algo_label->SetTextFont(62);
  algo_label->SetTextSize(0.05); // 0.05
  algo_label->SetTextColor(kGray+1);
  algo_label->SetNDC();
  algo_label->Draw();

  TLatex *eta_text = new TLatex(coord->ConvertGraphXToPadX(1-(0.05*c->GetWh()/c->GetWw())), coord->ConvertGraphYToPadY(0.867), "#left|#eta^{rec}#right| < 2.5");
  eta_text->SetTextAlign(33); // right top
  eta_text->SetTextFont(42);
  eta_text->SetTextSize(0.035);
  eta_text->SetNDC();
  eta_text->Draw();

  // string string_text_top_right = kYears.at(year).name + " (CMSSW 10.6.X)";
  // string string_text_top_right = "POWHEG #plus Pythia8, pp #rightarrow t#bar{t} #rightarrow all-jets @ #sqrt{#it{s}} = 13 TeV";
  string string_text_top_right = "t#bar{t} #rightarrow all-jets @ #sqrt{#it{s}} = 13 TeV";
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

  TLatex *prelim = new TLatex(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.87), "Simulation, Private Work");
  prelim->SetTextAlign(13); // left top
  prelim->SetTextFont(52);
  prelim->SetTextSize(0.035);
  prelim->SetNDC();
  prelim->Draw();

  // https://root-forum.cern.ch/t/inconsistent-tick-length/18563/8
  c->Update();
  double tickScaleX = (c->GetUxmax() - c->GetUxmin()) / (c->GetX2() - c->GetX1()) * (c->GetWh() * c->GetAbsHNDC());
  double tickScaleY = (c->GetUymax() - c->GetUymin()) / (c->GetY2() - c->GetY1()) * (c->GetWw() * c->GetAbsWNDC());
  mg->GetHistogram()->GetYaxis()->SetTickLength(c->GetWw() * tick_length / tickScaleY);
  mg->GetHistogram()->GetXaxis()->SetTickLength(c->GetWh() * tick_length / tickScaleX);
  // graph_raw->GetHistogram()->GetYaxis()->SetTickLength(c->GetWw() * tick_length / tickScaleY);
  // graph_raw->GetHistogram()->GetXaxis()->SetTickLength(c->GetWh() * tick_length / tickScaleX);

  gPad->Update();
  gPad->RedrawAxis();

  // Save to disk
  string plotName = (string)"plot_"+kProbeJetAlgos.at(probejetalgo).name+"_ptResponse_"+dr_max+"_"+kYears.at(year).name;
  plotName += (log_y ? string("_log") : string("_lin"))+".pdf";
  string plotBasePath = "plots";
  gSystem->Exec(((string)"mkdir -p "+plotBasePath).c_str());
  string plotPath = plotBasePath+"/"+plotName;
  c->SaveAs(plotPath.c_str());
  delete c;
}


void plots_HOTVR_ptResponse() {

  // vector<Year> years;
  // years.push_back(Year{"UL16preVFP", "Ultra Legacy 2016 early", kPink+8});
  // years.push_back(Year{"UL16postVFP", "Ultra Legacy 2016 late", kAzure-6});
  // years.push_back(Year{"UL17", "Ultra Legacy 2017", kOrange-3});
  // years.push_back(Year{"UL18", "Ultra Legacy 2018", kSpring-8});

  const vector<Year> years = { Year::isUL16preVFP, Year::isUL16postVFP, Year::isUL17, Year::isUL18 };

  for(const Year & year : years) {
    for(unsigned int i = 0; i <= 0; i++) {
      const string dr_string = (string)"dr"+to_string(i/10)+"p"+to_string(i-(i/10)*10); // e.g. i=8 will be converted to "dr0p8" and i=12 will be converted to "dr1p2"
      cout << dr_string << endl;
      do_plot(year, ProbeJetAlgo::isHOTVR, dr_string, false);
      // do_plot(year, ProbeJetAlgo::isAK8, dr_string, false);
    }
  }
}
