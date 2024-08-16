#include <iomanip>
#include <sstream>

#include "../constants.h"

using namespace macros;

typedef struct {
  string name;
  string text;
  int summary_plot_linecolor;
  int summary_plot_linestyle;
  int summary_plot_linewidth=2;
} PtBin;

typedef struct {
  float target_eff_bkg;
  float real_eff_bkg;
  float real_eff_sig;
  float scan_cut;
  int scan_cut_i; // index of point of TGraph at which the cut is
  // float real_eff_qcd_with_btag;
  // float real_eff_ttbar_with_btag;
} WorkingPoint;

typedef struct {
  ProbeJetAlgo probejetalgo;
  string name_base; //ak8_t__tau
  vector<string> name_variants; //ak8_t_incl__tau, ak8_t_btag__tau
  vector<float> target_effs_bkg;
  map<string, vector<WorkingPoint>> wps; //string: variant, vector: WPs ordered
  vector<int> linestyle_variants; //ak8_t_incl__tau, ak8_t_btag__tau
  vector<int> linecolor_variants; //ak8_t_incl__tau, ak8_t_btag__tau
  string legend_base;
  vector<string> legend_variants; //ak8_t_incl__tau, ak8_t_btag__tau
  string scan_var_tlatex;
  string scan_var_tlatex_axis;
} Tagger;

// class CoordinateConverter {
// public:
//   CoordinateConverter() {};
//   void init(const float l, const float r, const float b, const float t);
//   float ConvertGraphXToPadX(const float graph_x);
//   float ConvertGraphYToPadY(const float graph_y);
// private:
//   float pad_l, pad_r, pad_b, pad_t; // pad margins
//   float graph_width, graph_height;
// };
//
// void CoordinateConverter::init(const float l, const float r, const float b, const float t) {
//   pad_l = l;
//   pad_r = r;
//   pad_b = b;
//   pad_t = t;
//   graph_width = 1.-l-r;
//   graph_height = 1.-b-t;
// }
//
// float CoordinateConverter::ConvertGraphXToPadX(const float graph_x) {
//   return pad_l+graph_x*graph_width;
// }
//
// float CoordinateConverter::ConvertGraphYToPadY(const float graph_y) {
//   return pad_b+graph_y*graph_height;
// }

// void min0max1(TGraph * graph) {
//   double x, y;
//   for(int i = 0; i < graph->GetN(); i++) {
//     graph->GetPoint(i, x, y);
//     if(x < 0.) graph->SetPoint(i, 0., y);
//     if(x > 1.) graph->SetPoint(i, 1., y);
//     if(y < 0.) graph->SetPoint(i, x, 0.);
//     if(y > 1.) graph->SetPoint(i, x, 1.);
//   }
// }

WorkingPoint init_working_point(const float eff, const TGraph * graph_eff_bkg, const TGraph * graph_eff_sig) {
  WorkingPoint result = { .target_eff_bkg = eff };
  int scan_cut_i;
  for(int i = graph_eff_bkg->GetN()-1; i >= 0; i--) {
    double x, y;
    graph_eff_bkg->GetPoint(i, x, y);
    if(y < eff) {
      result.real_eff_bkg = y;
      result.scan_cut = x;
      scan_cut_i = i;
      if(i != graph_eff_bkg->GetN()-1) {
        double x_before, y_before;
        graph_eff_bkg->GetPoint(i+1, x_before, y_before);
        if(abs(y_before - eff) < abs(y - eff)) { result.real_eff_bkg = y_before; result.scan_cut = x_before; scan_cut_i = i+1; }
      }
      break;
    }
  }
  double x, y;
  graph_eff_sig->GetPoint(scan_cut_i, x, y);
  result.real_eff_sig = y;
  result.scan_cut_i = scan_cut_i;
  return result;
}

// WorkingPoint init_working_point(float eff, const TGraph * graph_eff_qcd, const TGraph * graph_eff_ttbar, const TGraph * graph_eff_qcd_btag, const TGraph * graph_eff_ttbar_btag) {
//   WorkingPoint result = { .target_eff_qcd = eff };
//   int index_tau32cut;
//   for(int i = graph_eff_qcd->GetN()-1; i >= 0; i--) {
//     double x, y;
//     graph_eff_qcd->GetPoint(i, x, y);
//     if(y < eff) {
//       result.real_eff_qcd = y;
//       result.tau32cut = x;
//       index_tau32cut = i;
//       if(i != graph_eff_qcd->GetN()-1) {
//         double x_before, y_before;
//         graph_eff_qcd->GetPoint(i+1, x_before, y_before);
//         if(abs(y_before - eff) < abs(y - eff)) { result.real_eff_qcd = y_before; result.tau32cut = x_before; index_tau32cut = i+1; }
//       }
//       break;
//     }
//   }
//   double x, y;
//   graph_eff_ttbar->GetPoint(index_tau32cut, x, y);
//   result.real_eff_ttbar = y;
//   graph_eff_qcd_btag->GetPoint(index_tau32cut, x, y);
//   result.real_eff_qcd_with_btag = y;
//   graph_eff_ttbar_btag->GetPoint(index_tau32cut, x, y);
//   result.real_eff_ttbar_with_btag = y;
//   result.tau32cut_i = index_tau32cut;
//   return result;
// }

void print_working_point(const WorkingPoint & wp) {
  cout << ">>> WP: <<<" << endl;
  cout << "Target eff BKG:    " << wp.target_eff_bkg << endl;
  cout << "Actual eff BKG:    " << wp.real_eff_bkg << endl;
  cout << "Actual eff SIG:    " << wp.real_eff_sig << endl;
  cout << "Scan variable cut: " << wp.scan_cut << endl;
}

void calculate_working_points(const Year & year, Tagger & tagger, const PtBin & pt_bin, const bool print=false) {
  const string infilePath = (string)getenv("CMSSW_BASE")+"/src/UHH2/LegacyTopTagging/output/WorkingPointStudy/"+kYears.at(year).name+"/nominal/"+tagger.name_base+"/"+pt_bin.name+"/graphs.root";
  TFile* infile = TFile::Open(infilePath.c_str(), "READ");
  const TGraph* graph_eff_bkg = (TGraph*)infile->Get("eff_bkg");
  const TGraph* graph_eff_sig = (TGraph*)infile->Get("eff_sig");

  for(const auto & eff : tagger.target_effs_bkg) {
    WorkingPoint wp = init_working_point(eff, graph_eff_bkg, graph_eff_sig);
    tagger.wps[tagger.name_base].push_back(wp);
    if(print) print_working_point(wp);
  }
  infile->Close();

  for(const auto & variant : tagger.name_variants) {
    const string infilePath_var = (string)getenv("CMSSW_BASE")+"/src/UHH2/LegacyTopTagging/output/WorkingPointStudy/"+kYears.at(year).name+"/nominal/"+variant+"/"+pt_bin.name+"/graphs.root";
    TFile* infile_var = TFile::Open(infilePath_var.c_str(), "READ");
    const TGraph* graph_eff_bkg_var = (TGraph*)infile_var->Get("eff_bkg");
    const TGraph* graph_eff_sig_var = (TGraph*)infile_var->Get("eff_sig");
    for(const auto & wp_base : tagger.wps.at(tagger.name_base)) {
      double x, y;
      WorkingPoint wp_var = { .target_eff_bkg = -1., .scan_cut = wp_base.scan_cut, .scan_cut_i = wp_base.scan_cut_i }; // target eff only valid for base tagger, e.g. ak8_t__tau
      graph_eff_bkg_var->GetPoint(wp_base.scan_cut_i, x, y);
      wp_var.real_eff_bkg = y;
      graph_eff_sig_var->GetPoint(wp_base.scan_cut_i, x, y);
      wp_var.real_eff_sig = y;
      tagger.wps[variant].push_back(wp_var);
      if(print) print_working_point(wp_var);
    }
    infile_var->Close();
  }
}


void get_reference_working_points(const Year & year, Tagger & tagger, const PtBin & pt_bin, const Tagger & reference_tagger, const bool print=false) {
  const string infilePath = (string)getenv("CMSSW_BASE")+"/src/UHH2/LegacyTopTagging/output/WorkingPointStudy/"+kYears.at(year).name+"/nominal/"+tagger.name_base+"/"+pt_bin.name+"/graphs.root";
  TFile* infile = TFile::Open(infilePath.c_str(), "READ");
  const TGraph* graph_eff_bkg = (TGraph*)infile->Get("eff_bkg");
  const TGraph* graph_eff_sig = (TGraph*)infile->Get("eff_sig");

  for(const auto & wp_ref : reference_tagger.wps.at(reference_tagger.name_base)) {
    double x, y;
    WorkingPoint wp = { .target_eff_bkg = -1., .scan_cut = wp_ref.scan_cut, .scan_cut_i = wp_ref.scan_cut_i };
    graph_eff_bkg->GetPoint(wp_ref.scan_cut_i, x, y);
    wp.real_eff_bkg = y;
    graph_eff_sig->GetPoint(wp_ref.scan_cut_i, x, y);
    wp.real_eff_sig = y;
    tagger.wps[tagger.name_base].push_back(wp);
    if(print) print_working_point(wp);
  }
  infile->Close();

  for(const auto & variant : tagger.name_variants) {
    const string infilePath_var = (string)getenv("CMSSW_BASE")+"/src/UHH2/LegacyTopTagging/output/WorkingPointStudy/"+kYears.at(year).name+"/nominal/"+variant+"/"+pt_bin.name+"/graphs.root";
    TFile* infile_var = TFile::Open(infilePath_var.c_str(), "READ");
    const TGraph* graph_eff_bkg_var = (TGraph*)infile_var->Get("eff_bkg");
    const TGraph* graph_eff_sig_var = (TGraph*)infile_var->Get("eff_sig");
    for(const auto & wp_ref : reference_tagger.wps.at(reference_tagger.name_base)) {
      double x, y;
      WorkingPoint wp_var = { .target_eff_bkg = -1., .scan_cut = wp_ref.scan_cut, .scan_cut_i = wp_ref.scan_cut_i }; // target eff only valid for base tagger, e.g. ak8_t__tau
      graph_eff_bkg_var->GetPoint(wp_ref.scan_cut_i, x, y);
      wp_var.real_eff_bkg = y;
      graph_eff_sig_var->GetPoint(wp_ref.scan_cut_i, x, y);
      wp_var.real_eff_sig = y;
      tagger.wps[variant].push_back(wp_var);
      if(print) print_working_point(wp_var);
    }
    infile_var->Close();
  }
}


// vector<WorkingPoint> get_reference_working_points(const string & year, const PtBin & pt_bin, const vector<WorkingPoint> & wps_reference, const bool print=false) {
//   string infileBasePath = (string)getenv("CMSSW_BASE")+"/src/UHH2/LegacyTopTagging/output/WorkingPointStudy/"+year+"/workdir_npy/"+pt_bin.name+"/";
//   string infileName = "root_"+pt_bin.name+".root";
//   string infilePath = infileBasePath+infileName;
//
//   TFile* infile = TFile::Open(infilePath.c_str(), "READ");
//   const TGraph* graph_eff_qcd_msd = (TGraph*)infile->Get("eff_qcd_msd");
//   const TGraph* graph_eff_qcd_msd_btag = (TGraph*)infile->Get("eff_qcd_msd_deepcsv");
//   const TGraph* graph_eff_ttbar_msd = (TGraph*)infile->Get("eff_ttbar_msd");
//   const TGraph* graph_eff_ttbar_msd_btag = (TGraph*)infile->Get("eff_ttbar_msd_deepcsv");
//
//   vector<WorkingPoint> working_points;
//   double x=0, y=0;
//   for(int i = 0; i < wps_reference.size(); i++) {
//     WorkingPoint wp_ref = wps_reference.at(i);
//     WorkingPoint wp_ref_this_bin;
//     wp_ref_this_bin.target_eff_qcd = wp_ref.target_eff_qcd;
//     wp_ref_this_bin.tau32cut = wp_ref.tau32cut;
//     wp_ref_this_bin.tau32cut_i = wp_ref.tau32cut_i;
//     graph_eff_qcd_msd->GetPoint(wp_ref.tau32cut_i, x, y); wp_ref_this_bin.real_eff_qcd = y;
//     graph_eff_qcd_msd_btag->GetPoint(wp_ref.tau32cut_i, x, y); wp_ref_this_bin.real_eff_qcd_with_btag = y;
//     graph_eff_ttbar_msd->GetPoint(wp_ref.tau32cut_i, x, y); wp_ref_this_bin.real_eff_ttbar = y;
//     graph_eff_ttbar_msd_btag->GetPoint(wp_ref.tau32cut_i, x, y); wp_ref_this_bin.real_eff_ttbar_with_btag = y;
//     working_points.push_back(wp_ref_this_bin);
//   }
//   if(print) for(auto wp : working_points) {
//     cout << ">>> WP derived from reference bin <<<" << endl;
//     cout << "Target QCD eff: " << wp.target_eff_qcd << endl;
//     cout << "Actual QCD eff: " << wp.real_eff_qcd << endl;
//     cout << "Actual tt eff:  " << wp.real_eff_ttbar << endl;
//     cout << "tau32 cut:      " << wp.tau32cut << endl;
//     cout << "Actual QCD eff with b-tag: " << wp.real_eff_qcd_with_btag << endl;
//     cout << "Actual tt eff with b-tag:  " << wp.real_eff_ttbar_with_btag << endl;
//   };
//   return working_points;
// }


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


void do_plot(
  const Year & year,
  const Tagger & tagger,
  const PtBin & pt_bin,
  const string & graph_base_name,
  const bool log_y
  // const vector<WorkingPoint> & wps_this_bin,
  // const vector<WorkingPoint> & wps_this_bin_ref,
  // const bool do_raw = true,
  // const bool do_msd = true,
  // const bool do_msd_btag = true
)
{
  bool is_eff_bkg = graph_base_name == "eff_bkg";
  bool is_eff_sig = graph_base_name == "eff_sig";
  bool is_roc = graph_base_name == "roc";

  vector<TGraph> graphs;
  TGraph graph_pipe;
  vector<string> legends;

  const string infilePath = (string)getenv("CMSSW_BASE")+"/src/UHH2/LegacyTopTagging/output/WorkingPointStudy/"+kYears.at(year).name+"/nominal/"+tagger.name_base+"/"+pt_bin.name+"/graphs.root";
  TFile* infile = TFile::Open(infilePath.c_str(), "READ");
  graph_pipe = *(TGraph*)infile->Get(graph_base_name.c_str());
  graph_pipe.SetLineWidth(2);
  graph_pipe.SetLineStyle(1);
  graph_pipe.SetLineColor(kRed);
  graphs.push_back(graph_pipe);
  legends.push_back(tagger.legend_base);
  infile->Close();

  int iterator = 0;
  for(const auto & variant : tagger.name_variants) {
    const string infilePath_var = (string)getenv("CMSSW_BASE")+"/src/UHH2/LegacyTopTagging/output/WorkingPointStudy/"+kYears.at(year).name+"/nominal/"+variant+"/"+pt_bin.name+"/graphs.root";
    TFile* infile_var = TFile::Open(infilePath_var.c_str(), "READ");
    graph_pipe = *(TGraph*)infile_var->Get(graph_base_name.c_str());
    graph_pipe.SetLineWidth(2);
    graph_pipe.SetLineStyle(tagger.linestyle_variants.at(iterator));
    graph_pipe.SetLineColor(tagger.linecolor_variants.at(iterator));
    graphs.push_back(graph_pipe);
    legends.push_back(tagger.legend_variants.at(iterator));
    infile_var->Close();
    iterator++;
  }


  // // if(!(do_raw || do_msd || do_msd_btag)) cout << "Nothing to plot." << endl;
  //
  // string infileBasePath = (string)getenv("CMSSW_BASE")+"/src/UHH2/LegacyTopTagging/output/WorkingPointStudy/"+year+"/workdir_npy/"+pt_bin.name+"/";
  // string infileName = "root_"+pt_bin.name+".root";
  // string infilePath = infileBasePath+infileName;
  //
  // TFile* infile = TFile::Open(infilePath.c_str(), "READ");
  // vector<TGraph*> graphs;
  // if(do_raw) graphs.push_back((TGraph*)infile->Get(graph_base_name.c_str()));
  // if(do_msd) graphs.push_back((TGraph*)infile->Get((graph_base_name+"_msd").c_str()));
  // if(do_msd_btag) graphs.push_back((TGraph*)infile->Get((graph_base_name+"_msd_deepcsv").c_str()));


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

  TMultiGraph *mg = new TMultiGraph("mg", "mg title");

  bool plot_at_bottom = false;
  if((is_eff_sig || is_roc) && log_y) plot_at_bottom = true;

  float legend_x1 = coord->ConvertGraphXToPadX(plot_at_bottom ? 0.45 : 0.05);
  float legend_y1 = coord->ConvertGraphYToPadY(plot_at_bottom ? 0.05 : 0.47);
  float legend_x2 = coord->ConvertGraphXToPadX(plot_at_bottom ? 0.95 : 0.55);
  float legend_y2 = coord->ConvertGraphYToPadY(plot_at_bottom ? 0.35 : 0.77);
  TLegend *legend = new TLegend(legend_x1, legend_y1, legend_x2, legend_y2); // x1 y1 x2 y2
  // const string header = ("#bf{Scan of "+tagger.scan_var_tlatex+" for "+pt_bin.text+"}");
  const string header = ("#bf{Scan of "+tagger.scan_var_tlatex+"}");
  legend->SetHeader(header.c_str());
  legend->SetTextSize(0.025);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  // graph-dependent attributes, order: raw, msd, msd_btag
  // vector<int> linewidths = {2, 2, 2};
  // vector<int> linestyles = {3, 1, 2};
  // vector<int> linecolors = {kBlue, kRed, kMagenta};
  // vector<string> legends = {"No #it{m}_{SD} window nor b tagging", "105 < #it{m}_{SD} [GeV] < 210", "#it{m}_{SD} as above + subjet b tag"};
  iterator = 0;
  for(TGraph & graph : graphs) {
    // min0max1(graph);
    mg->Add(&graph);
    // graph->SetLineWidth(linewidths.at(++iterator));
    // graph.SetLineWidth(2);
    // graph->SetLineStyle(linestyles.at(iterator));
    // graph->SetLineColor(linecolors.at(iterator));
    legend->AddEntry(&graph, legends.at(iterator).c_str(), "l");
    iterator++;
  }
  auto wp_graph = new_null_graph();
  wp_graph->SetMarkerColor(kBlack);
  wp_graph->SetMarkerSize(1.2);
  wp_graph->SetMarkerStyle(24);
  mg->Add(wp_graph);
  legend->AddEntry(wp_graph, "Working points");


  mg->Draw("ac");
  mg->SetTitle("");
  mg->GetHistogram()->GetXaxis()->SetRangeUser(0., 1.);
  double maximum_saved = mg->GetHistogram()->GetMaximum();
  double minimum_saved = mg->GetHistogram()->GetMinimum();
  mg->GetHistogram()->SetMaximum(1.);
  if(!log_y) mg->GetHistogram()->SetMinimum(0.);
  else mg->GetHistogram()->SetMinimum(0.0002);
  if(is_eff_bkg || is_eff_sig) mg->GetHistogram()->GetXaxis()->SetTitle(tagger.scan_var_tlatex_axis.c_str());
  else if(is_roc) mg->GetHistogram()->GetXaxis()->SetTitle("True positive rate, #varepsilon_{S}");
  if(is_eff_bkg || is_roc) mg->GetHistogram()->GetYaxis()->SetTitle("False positive rate, #varepsilon_{B}");
  if(is_eff_sig) mg->GetHistogram()->GetYaxis()->SetTitle("True positive rate, #varepsilon_{S}");
  mg->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
  mg->GetHistogram()->GetXaxis()->CenterTitle();
  mg->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);
  mg->GetHistogram()->GetYaxis()->CenterTitle();
  // mg->GetHistogram()->GetXaxis()->SetNdivisions(405);
  if(!log_y) mg->GetHistogram()->GetYaxis()->SetNdivisions(405);
  gStyle->SetLineWidth(1);

  legend->Draw();

  TLatex *text_top_left = new TLatex(margin_l, 1-(margin_t-0.01), (string("Ultra Legacy ")+kYears.at(year).nice_name).c_str());
  text_top_left->SetTextAlign(11); // left bottom aligned
  text_top_left->SetTextFont(42);
  text_top_left->SetTextSize(0.035);
  text_top_left->SetNDC();
  text_top_left->Draw();

  TLatex *algo_label = new TLatex(coord->ConvertGraphXToPadX(1-(0.05*c->GetWh()/c->GetWw())), coord->ConvertGraphYToPadY(0.95), (kProbeJetAlgos.at(tagger.probejetalgo).name+" PUPPI").c_str());
  algo_label->SetTextAlign(33); // right top
  algo_label->SetTextFont(62);
  algo_label->SetTextSize(0.05); // 0.05
  algo_label->SetTextColor(kGray+1);
  algo_label->SetNDC();
  algo_label->Draw();

  const double eta_pt_text_x = coord->ConvertGraphXToPadX(1-(0.05*c->GetWh()/c->GetWw()));
  double eta_text_y = coord->ConvertGraphYToPadY(0.801);
  double pt_text_y = coord->ConvertGraphYToPadY(0.867);
  if(!plot_at_bottom) {
    eta_text_y = coord->ConvertGraphYToPadY((0.05*c->GetWh()/c->GetWw()) + 0.066);
    pt_text_y = coord->ConvertGraphYToPadY((0.05*c->GetWh()/c->GetWw()));
  }

  TLatex *eta_text = new TLatex(eta_pt_text_x, eta_text_y, "|#eta| < 2.5");
  eta_text->SetTextAlign(plot_at_bottom ? 33 : 31); // right top
  eta_text->SetTextFont(42);
  eta_text->SetTextSize(0.035);
  eta_text->SetNDC();
  eta_text->Draw();

  TLatex *pt_text = new TLatex(eta_pt_text_x, pt_text_y, pt_bin.text.c_str());
  pt_text->SetTextAlign(plot_at_bottom ? 33 : 31); // right top
  pt_text->SetTextFont(42);
  pt_text->SetTextSize(0.035);
  pt_text->SetNDC();
  pt_text->Draw();

  // string string_text_top_right = "t#bar{t} #rightarrow all-jets @ #sqrt{#it{s}} = 13 TeV";
  string string_text_top_right = "#sqrt{#it{s}} = 13 TeV";
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

  if(!plot_at_bottom) {
    TText *prelim = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.87), "Simulation, Private Work");
    prelim->SetTextAlign(13); // left top
    prelim->SetTextFont(52);
    prelim->SetTextSize(0.035);
    prelim->SetNDC();
    prelim->Draw();
  }
  else {
    TText *prelim = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.87), "Simulation");
    prelim->SetTextAlign(13); // left top
    prelim->SetTextFont(52);
    prelim->SetTextSize(0.035);
    prelim->SetNDC();
    // prelim->SetTextColor(kWhite);
    prelim->Draw();

    TText *prelim2 = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.804), "Private Work");
    prelim2->SetTextAlign(13); // left top
    prelim2->SetTextFont(52);
    prelim2->SetTextSize(0.035);
    prelim2->SetNDC();
    // prelim2->SetTextColor(kWhite);
    prelim2->Draw();
  }

  // WP markers
  for(const auto & wp : tagger.wps.at(tagger.name_base)) {
    double x=0, y=0;
    if(is_roc) { x=wp.real_eff_sig; y=wp.real_eff_bkg; }
    else if(is_eff_bkg) { x=wp.scan_cut; y=wp.real_eff_bkg; }
    else if(is_eff_sig) { x=wp.scan_cut; y=wp.real_eff_sig; }
    TMarker *marker1 = new TMarker(x, y, 24);
    marker1->SetMarkerColor(kRed);
    marker1->SetMarkerSize(1.2);
    marker1->Draw();
  }
  iterator = 0;
  for(const auto & variant : tagger.name_variants) {
    if(variant.find("_incl_") != string::npos) {
      cout << "Skip printing WP markers for inclusive tagger: " << variant << endl;
      break;
    }
    for(const auto & wp : tagger.wps.at(variant)) {
      double x=0, y=0;
      if(is_roc) { x=wp.real_eff_sig; y=wp.real_eff_bkg; }
      else if(is_eff_bkg) { x=wp.scan_cut; y=wp.real_eff_bkg; }
      else if(is_eff_sig) { x=wp.scan_cut; y=wp.real_eff_sig; }
      TMarker *marker1 = new TMarker(x, y, 24);
      marker1->SetMarkerColor(tagger.linecolor_variants.at(iterator));
      marker1->SetMarkerSize(1.2);
      marker1->Draw();
    }
    iterator++;
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
  string plotName = "plot__"+tagger.name_base+"__"+pt_bin.name+"__"+graph_base_name+"__";
  // for(bool do_x : {do_raw, do_msd, do_msd_btag}) { do_x ? plotName += "1" : plotName += "0"; }
  plotName += (log_y ? string("log") : string("lin"))+".pdf";
  // string plotPath = infileBasePath+plotName;
  TString plotPath = (TString(infilePath)).ReplaceAll("graphs.root", "")+plotName;
  c->SaveAs(plotPath.Data());
  delete c;
}



void do_summary_plot(
  const Year & year,
  const Tagger & tagger,
  const string & variant,
  const vector<PtBin> & pt_bins,
  const string & graph_base_name,
  const bool log_y,
  // const vector<WorkingPoint> & wps_ref_pt480to600,
  // const vector<WorkingPoint> & wps_dpnote_pt480to600
  const Tagger & tagger_pt_480to600,
  const Tagger & tagger_pt_480to600_DPnote
  // const bool show_DPnote_WPs = true,
)
{
  const bool plot_WPs = year == Year::isUL18 && (variant == "ak8_t__tau" || variant == "ak8_t_btagDCSV__tau");

  const bool is_base_tagger = variant == tagger.name_base;

  const bool is_eff_bkg = graph_base_name == "eff_bkg";
  const bool is_eff_sig = graph_base_name == "eff_sig";
  const bool is_roc = graph_base_name == "roc";

  vector<TGraph> graphs;
  for(const auto & pt_bin : pt_bins) {
    const string infilePath = (string)getenv("CMSSW_BASE")+"/src/UHH2/LegacyTopTagging/output/WorkingPointStudy/"+kYears.at(year).name+"/nominal/"+variant+"/"+pt_bin.name+"/graphs.root";
    TFile* infile = TFile::Open(infilePath.c_str(), "READ");
    graphs.push_back(*(TGraph*)infile->Get(graph_base_name.c_str()));
  }

  // string infileBaseBasePath = (string)getenv("CMSSW_BASE")+"/src/UHH2/LegacyTopTagging/output/WorkingPointStudy/"+year+"/workdir_npy/";
  // vector<TGraph*> graphs;
  // for(auto pt_bin : pt_bins) {
  //   string infileBasePath = infileBaseBasePath+pt_bin.name+"/";
  //   string infileName = "root_"+pt_bin.name+".root";
  //   string infilePath = infileBasePath+infileName;
  //   TFile* infile = TFile::Open(infilePath.c_str(), "READ");
  //   graphs.push_back((TGraph*)infile->Get((graph_base_name+(is_raw ? "" : "_"+variant)).c_str()));
  // }

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

  TMultiGraph *mg = new TMultiGraph("mg", "mg title");

  bool plot_at_bottom = false;
  if((is_eff_sig || is_roc) && log_y) plot_at_bottom = true;

  float legend_x1 = coord->ConvertGraphXToPadX(plot_at_bottom ? 0.45 : 0.05);
  float legend_y1 = coord->ConvertGraphYToPadY(plot_at_bottom ? 0.05 : (plot_WPs ? 0.32 : 0.42)); //0.47
  float legend_x2 = coord->ConvertGraphXToPadX(plot_at_bottom ? 0.95 : 0.55);
  float legend_y2 = coord->ConvertGraphYToPadY(plot_at_bottom ? (plot_WPs ? 0.50 : 0.40) : 0.77); //0.35
  TLegend *legend = new TLegend(legend_x1, legend_y1, legend_x2, legend_y2); // x1 y1 x2 y2
  string legend_title = "FIXME";
  if(is_base_tagger) legend_title = tagger.legend_base;
  else {
    int legend_i = -1;
    for(const auto & variant_ : tagger.name_variants) {
      legend_i++;
      if(variant_ == variant) break;
    }
    legend_title = tagger.legend_variants.at(legend_i);
  }
  // if(is_raw) legend_title = "No #it{m}_{SD} window nor b tagging";
  // else if(is_msd) legend_title = "105 < #it{m}_{SD} [GeV] < 210";
  // else if(is_msd_btag) legend_title = "105 < #it{m}_{SD} [GeV] < 210 + subjet b tag";
  const string header = ("#bf{Scan of "+tagger.scan_var_tlatex+"}");
  // legend->SetHeader(("#bf{"+legend_title+"}").c_str());
  legend->SetHeader(header.c_str());
  legend->SetTextSize(0.025);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  for(int i = 0; i < graphs.size(); i++) {
    mg->Add(&graphs.at(i));
    graphs.at(i).SetLineWidth(pt_bins.at(i).summary_plot_linewidth);
    graphs.at(i).SetLineStyle(pt_bins.at(i).summary_plot_linestyle);
    graphs.at(i).SetLineColor(pt_bins.at(i).summary_plot_linecolor);
    legend->AddEntry(&graphs.at(i), pt_bins.at(i).text.c_str(), "l");
  }
  auto wp_graph = new_null_graph();
  wp_graph->SetMarkerColor(kBlack);
  wp_graph->SetMarkerSize(1.2);
  wp_graph->SetMarkerStyle(24);
  mg->Add(wp_graph);
  if(plot_WPs) legend->AddEntry(wp_graph, "WPs for 480 < #it{p}_{T}/GeV < 600");
  auto dp_graph = new_null_graph();
  dp_graph->SetMarkerColor(kRed);
  dp_graph->SetMarkerSize(1.2);
  dp_graph->SetMarkerStyle(3);
  mg->Add(dp_graph);
  if(plot_WPs) legend->AddEntry(dp_graph, "WPs from CMS DP 2020/025");
  // TGraph * null_graph1 = new_null_graph();
  // mg->Add(null_graph1);

  mg->Draw("ac");
  mg->SetTitle("");
  mg->GetHistogram()->GetXaxis()->SetRangeUser(0., 1.);
  double maximum_saved = mg->GetHistogram()->GetMaximum();
  double minimum_saved = mg->GetHistogram()->GetMinimum();
  mg->GetHistogram()->SetMaximum(1.);
  if(!log_y) mg->GetHistogram()->SetMinimum(0.);
  else mg->GetHistogram()->SetMinimum(0.0002);
  if(is_eff_bkg || is_eff_sig) mg->GetHistogram()->GetXaxis()->SetTitle(tagger.scan_var_tlatex_axis.c_str());
  else if(is_roc) mg->GetHistogram()->GetXaxis()->SetTitle("True positive rate, #varepsilon_{S}");
  if(is_eff_bkg || is_roc) mg->GetHistogram()->GetYaxis()->SetTitle("False positive rate, #varepsilon_{B}");
  if(is_eff_sig) mg->GetHistogram()->GetYaxis()->SetTitle("True positive rate, #varepsilon_{S}");
  mg->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
  mg->GetHistogram()->GetXaxis()->CenterTitle();
  mg->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);
  mg->GetHistogram()->GetYaxis()->CenterTitle();
  // mg->GetHistogram()->GetXaxis()->SetNdivisions(405);
  if(!log_y) mg->GetHistogram()->GetYaxis()->SetNdivisions(405);
  gStyle->SetLineWidth(1);

  legend->Draw();



  // WP markers for Pt480to600
  if(plot_WPs) {
    for(const auto & wp : tagger_pt_480to600.wps.at(variant)) {
      double x=0, y=0;
      if(is_roc) { x=wp.real_eff_sig; y=wp.real_eff_bkg; }
      else if(is_eff_bkg) { x=wp.scan_cut; y=wp.real_eff_bkg; }
      else if(is_eff_sig) { x=wp.scan_cut; y=wp.real_eff_sig; }
      TMarker *marker1 = new TMarker(x, y, 24);
      marker1->SetMarkerColor(kBlack);
      marker1->SetMarkerSize(1.2);
      marker1->Draw();
    }
    for(const auto & wp : tagger_pt_480to600_DPnote.wps.at(variant)) {
      double x=0, y=0;
      if(is_roc) { x=wp.real_eff_sig; y=wp.real_eff_bkg; }
      else if(is_eff_bkg) { x=wp.scan_cut; y=wp.real_eff_bkg; }
      else if(is_eff_sig) { x=wp.scan_cut; y=wp.real_eff_sig; }
      TMarker *marker1 = new TMarker(x, y, 3);
      marker1->SetMarkerColor(kRed);
      marker1->SetMarkerSize(1.2);
      marker1->Draw();
    }
  }




  TLatex *text_top_left = new TLatex(margin_l, 1-(margin_t-0.01), (string("Ultra Legacy ")+kYears.at(year).nice_name).c_str());
  text_top_left->SetTextAlign(11); // left bottom aligned
  text_top_left->SetTextFont(42);
  text_top_left->SetTextSize(0.035);
  text_top_left->SetNDC();
  text_top_left->Draw();

  TLatex *algo_label = new TLatex(coord->ConvertGraphXToPadX(1-(0.05*c->GetWh()/c->GetWw())), coord->ConvertGraphYToPadY(0.95), (kProbeJetAlgos.at(tagger.probejetalgo).name+" PUPPI").c_str());
  algo_label->SetTextAlign(33); // right top
  algo_label->SetTextFont(62);
  algo_label->SetTextSize(0.05); // 0.05
  algo_label->SetTextColor(kGray+1);
  algo_label->SetNDC();
  algo_label->Draw();

  const double eta_pt_text_x = coord->ConvertGraphXToPadX(1-(0.05*c->GetWh()/c->GetWw()));
  double eta_text_y = coord->ConvertGraphYToPadY(0.801);
  double pt_text_y = coord->ConvertGraphYToPadY(0.867);
  if(!plot_at_bottom) {
    eta_text_y = coord->ConvertGraphYToPadY((0.05*c->GetWh()/c->GetWw()) + 0.066);
    pt_text_y = coord->ConvertGraphYToPadY((0.05*c->GetWh()/c->GetWw()));
  }

  TLatex *eta_text = new TLatex(eta_pt_text_x, eta_text_y, "|#eta| < 2.5");
  eta_text->SetTextAlign(plot_at_bottom ? 33 : 31); // right top
  eta_text->SetTextFont(42);
  eta_text->SetTextSize(0.035);
  eta_text->SetNDC();
  eta_text->Draw();

  TLatex *pt_text = new TLatex(eta_pt_text_x, pt_text_y, legend_title.c_str());
  pt_text->SetTextAlign(plot_at_bottom ? 33 : 31); // right top
  pt_text->SetTextFont(42);
  pt_text->SetTextSize(0.035);
  pt_text->SetNDC();
  pt_text->Draw();

  // string string_text_top_right = "t#bar{t} #rightarrow all-jets @ #sqrt{#it{s}} = 13 TeV";
  string string_text_top_right = "#sqrt{#it{s}} = 13 TeV";
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

  if(!plot_at_bottom) {
    TText *prelim = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.87), "Simulation, Private Work");
    prelim->SetTextAlign(13); // left top
    prelim->SetTextFont(52);
    prelim->SetTextSize(0.035);
    prelim->SetNDC();
    prelim->Draw();
  }
  else {
    TText *prelim = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.87), "Simulation");
    prelim->SetTextAlign(13); // left top
    prelim->SetTextFont(52);
    prelim->SetTextSize(0.035);
    prelim->SetNDC();
    // prelim->SetTextColor(kWhite);
    prelim->Draw();

    TText *prelim2 = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.804), "Private Work");
    prelim2->SetTextAlign(13); // left top
    prelim2->SetTextFont(52);
    prelim2->SetTextSize(0.035);
    prelim2->SetNDC();
    // prelim2->SetTextColor(kWhite);
    prelim2->Draw();
  }


  //   for(int i = 0; i < wps_ref_pt480to600.size(); i++) {
  //     double x=0, y=0;
  //
  //     WorkingPoint wp = wps_ref_pt480to600.at(i);
  //
  //     if(is_msd) {
  //       if(is_roc) { x=wp.real_eff_ttbar; y=wp.real_eff_qcd; }
  //       else if(is_eff_bkg) { x=wp.tau32cut; y=wp.real_eff_qcd; }
  //       else if(is_eff_sig) { x=wp.tau32cut; y=wp.real_eff_ttbar; }
  //     }
  //     else if(is_msd_btag) {
  //       if(is_roc) { x=wp.real_eff_ttbar_with_btag; y=wp.real_eff_qcd_with_btag; }
  //       else if(is_eff_bkg) { x=wp.tau32cut; y=wp.real_eff_qcd_with_btag; }
  //       else if(is_eff_sig) { x=wp.tau32cut; y=wp.real_eff_ttbar_with_btag; }
  //     }
  //     TMarker *marker = new TMarker(x, y, 20);
  //     marker->SetMarkerColor(kBlack);
  //     marker->SetMarkerSize(1.2);
  //     marker->Draw();
  //   }
  //   for(int i = 0; i < wps_dpnote_pt480to600.size(); i++) {
  //     double x=0, y=0;
  //
  //     WorkingPoint wp = wps_dpnote_pt480to600.at(i);
  //
  //     if(is_msd) {
  //       if(is_roc) { x=wp.real_eff_ttbar; y=wp.real_eff_qcd; }
  //       else if(is_eff_bkg) { x=wp.tau32cut; y=wp.real_eff_qcd; }
  //       else if(is_eff_sig) { x=wp.tau32cut; y=wp.real_eff_ttbar; }
  //     }
  //     else if(is_msd_btag) {
  //       if(is_roc) { x=wp.real_eff_ttbar_with_btag; y=wp.real_eff_qcd_with_btag; }
  //       else if(is_eff_bkg) { x=wp.tau32cut; y=wp.real_eff_qcd_with_btag; }
  //       else if(is_eff_sig) { x=wp.tau32cut; y=wp.real_eff_ttbar_with_btag; }
  //     }
  //     TMarker *marker = new TMarker(x, y, 3);
  //     marker->SetMarkerColor(kRed);
  //     marker->SetMarkerSize(1.2);
  //     marker->Draw();
  //   }
  // }

  // https://root-forum.cern.ch/t/inconsistent-tick-length/18563/8
  c->Update();
  double tickScaleX = (c->GetUxmax() - c->GetUxmin()) / (c->GetX2() - c->GetX1()) * (c->GetWh() * c->GetAbsHNDC());
  double tickScaleY = (c->GetUymax() - c->GetUymin()) / (c->GetY2() - c->GetY1()) * (c->GetWw() * c->GetAbsWNDC());
  mg->GetHistogram()->GetYaxis()->SetTickLength(c->GetWw() * tick_length / tickScaleY);
  mg->GetHistogram()->GetXaxis()->SetTickLength(c->GetWh() * tick_length / tickScaleX);

  gPad->Update();
  gPad->RedrawAxis();

  // Save to disk
  string plotName = "plot__"+variant+"__comparePtBins__"+graph_base_name+"__";
  // for(bool do_x : {do_raw, do_msd, do_msd_btag}) { do_x ? plotName += "1" : plotName += "0"; }
  plotName += (log_y ? string("log") : string("lin"))+".pdf";
  // string plotPath = infileBasePath+plotName;
  TString plotPath = (string)getenv("CMSSW_BASE")+"/src/UHH2/LegacyTopTagging/output/WorkingPointStudy/"+kYears.at(year).name+"/nominal/"+variant+"/"+plotName;
  c->SaveAs(plotPath.Data());
  delete c;
}


// string my_format(float value, int prec=2, int factor=100) {
//   std::stringstream stream;
//   stream << std::fixed << std::setprecision(prec) << value*factor;
//   return stream.str();
// }
//
//
// void print_latex_table(const vector<WorkingPoint> & wps) {
//   cout << "WPs in LaTeX tabular format:" << endl;
//   for(auto wp : wps) {
//     cout << "$< " << my_format(wp.tau32cut, 3, 1) << "$ & & " << my_format(wp.real_eff_qcd) << " & " << my_format(wp.real_eff_ttbar) << " & & " << my_format(wp.real_eff_qcd_with_btag) << " & " << my_format(wp.real_eff_ttbar_with_btag) << "\\\\" << endl;
//   }
// }


// vector<WorkingPoint> init_DPnote_WPs() {
//   // values taken from CMS DP 2020-025
//   vector<WorkingPoint> wps;
//   wps.push_back(WorkingPoint{ .scan_cut=0.80, .real_eff_bkg=0.159, .real_eff_sig=0.62, .real_eff_qcd_with_btag=0.053, .real_eff_ttbar_with_btag=0.53 });
//   wps.push_back(WorkingPoint{ .scan_cut=0.65, .real_eff_bkg=0.051, .real_eff_sig=0.49, .real_eff_qcd_with_btag=0.018, .real_eff_ttbar_with_btag=0.43 });
//   wps.push_back(WorkingPoint{ .scan_cut=0.54, .real_eff_bkg=0.017, .real_eff_sig=0.37, .real_eff_qcd_with_btag=0.006, .real_eff_ttbar_with_btag=0.33 });
//   wps.push_back(WorkingPoint{ .scan_cut=0.46, .real_eff_bkg=0.005, .real_eff_sig=0.26, .real_eff_qcd_with_btag=0.003, .real_eff_ttbar_with_btag=0.23 });
//   wps.push_back(WorkingPoint{ .scan_cut=0.40, .real_eff_bkg=0.002, .real_eff_sig=0.17, .real_eff_qcd_with_btag=0.001, .real_eff_ttbar_with_btag=0.16 });
//   return wps;
// }


Tagger init_DPnote_Tagger(const Tagger & tagger) {
  // values taken from CMS DP 2020-025
  Tagger tagger_DPnote = tagger;
  vector<WorkingPoint> wps__ak8_t__tau;
  wps__ak8_t__tau.push_back(WorkingPoint{ .scan_cut=0.80, .real_eff_bkg=0.159, .real_eff_sig=0.62, });
  wps__ak8_t__tau.push_back(WorkingPoint{ .scan_cut=0.65, .real_eff_bkg=0.051, .real_eff_sig=0.49, });
  wps__ak8_t__tau.push_back(WorkingPoint{ .scan_cut=0.54, .real_eff_bkg=0.017, .real_eff_sig=0.37, });
  wps__ak8_t__tau.push_back(WorkingPoint{ .scan_cut=0.46, .real_eff_bkg=0.005, .real_eff_sig=0.26, });
  wps__ak8_t__tau.push_back(WorkingPoint{ .scan_cut=0.40, .real_eff_bkg=0.002, .real_eff_sig=0.17, });
  tagger_DPnote.wps["ak8_t__tau"] = wps__ak8_t__tau;
  vector<WorkingPoint> wps__ak8_t_btagDCSV__tau;
  wps__ak8_t_btagDCSV__tau.push_back(WorkingPoint{ .scan_cut=0.80, .real_eff_bkg=0.053, .real_eff_sig=0.53, });
  wps__ak8_t_btagDCSV__tau.push_back(WorkingPoint{ .scan_cut=0.65, .real_eff_bkg=0.018, .real_eff_sig=0.43, });
  wps__ak8_t_btagDCSV__tau.push_back(WorkingPoint{ .scan_cut=0.54, .real_eff_bkg=0.006, .real_eff_sig=0.33, });
  wps__ak8_t_btagDCSV__tau.push_back(WorkingPoint{ .scan_cut=0.46, .real_eff_bkg=0.003, .real_eff_sig=0.23, });
  wps__ak8_t_btagDCSV__tau.push_back(WorkingPoint{ .scan_cut=0.40, .real_eff_bkg=0.001, .real_eff_sig=0.16, });
  tagger_DPnote.wps["ak8_t_btagDCSV__tau"] = wps__ak8_t_btagDCSV__tau;
  return tagger_DPnote;
}



void do_plot_all_years(
  const Tagger & tagger,
  const PtBin & pt_bin,
  const string & graph_base_name,
  const bool log_y
  // const vector<WorkingPoint> & wps_this_bin,
  // const vector<WorkingPoint> & wps_this_bin_ref,
  // const bool do_raw = true,
  // const bool do_msd = true,
  // const bool do_msd_btag = true
)
{
  bool is_eff_bkg = graph_base_name == "eff_bkg";
  bool is_eff_sig = graph_base_name == "eff_sig";
  bool is_roc = graph_base_name == "roc";

  vector<TGraph> graphs;
  TGraph graph_pipe;
  vector<string> legends;

  const string infilePathBase = (string)getenv("CMSSW_BASE")+"/src/UHH2/LegacyTopTagging/output/WorkingPointStudy/";
  string infilePath;
  for(const auto & year : kYears) {
    infilePath = infilePathBase+year.second.name+"/nominal/"+tagger.name_base+"/"+pt_bin.name+"/graphs.root";
    TFile* infile = TFile::Open(infilePath.c_str(), "READ");
    graph_pipe = *(TGraph*)infile->Get(graph_base_name.c_str());
    graph_pipe.SetLineWidth(2);
    graph_pipe.SetLineStyle(year.second.linestyle);
    graph_pipe.SetLineColor(year.second.tcolor);
    graphs.push_back(graph_pipe);
    TString year_nice_name = year.second.nice_name;
    year_nice_name.ReplaceAll("Ultra Legacy ", "");
    legends.push_back(string(year_nice_name.Data()));
    infile->Close();
  }


  // WP markers
  vector<TMarker*> markers;
  Tagger tagger_pipe = tagger;
  vector<float> wps_for_legend; // HACK
  for(const auto & year : kYears) {
    tagger_pipe.wps.clear();
    calculate_working_points(year.first, tagger_pipe, pt_bin, false);
    for(const auto & wp : tagger_pipe.wps.at(tagger_pipe.name_base)) {
      double x=0, y=0;
      if(is_roc) { x=wp.real_eff_sig; y=wp.real_eff_bkg; }
      else if(is_eff_bkg) { x=wp.scan_cut; y=wp.real_eff_bkg; }
      else if(is_eff_sig) { x=wp.scan_cut; y=wp.real_eff_sig; }
      TMarker *marker1 = new TMarker(x, y, 24);
      marker1->SetMarkerColor(year.second.tcolor);
      marker1->SetMarkerSize(1.2);
      // marker1->Draw();
      markers.push_back(marker1);
      if(wps_for_legend.size() < kYears.size()) wps_for_legend.push_back(wp.scan_cut); // HACK
    }
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

  TMultiGraph *mg = new TMultiGraph("mg", "mg title");

  bool plot_at_bottom = false;
  if((is_eff_sig || is_roc) && log_y) plot_at_bottom = true;

  float legend_x1 = coord->ConvertGraphXToPadX(plot_at_bottom ? 0.45 : 0.05);
  float legend_y1 = coord->ConvertGraphYToPadY(plot_at_bottom ? 0.05 : 0.47);
  float legend_x2 = coord->ConvertGraphXToPadX(plot_at_bottom ? 0.95 : 0.55);
  float legend_y2 = coord->ConvertGraphYToPadY(plot_at_bottom ? 0.35 : 0.77);

  if(is_eff_bkg && log_y) {
    legend_y1 = coord->ConvertGraphYToPadY(0.2);
    legend_y2 = coord->ConvertGraphYToPadY(0.5);
  }

  TLegend *legend = new TLegend(legend_x1, legend_y1, legend_x2, legend_y2); // x1 y1 x2 y2
  // const string header = ("#bf{Scan of "+tagger.scan_var_tlatex+" for "+pt_bin.text+"}");
  const string header = ("#bf{Scan of "+tagger.scan_var_tlatex+"}");
  legend->SetHeader(header.c_str());
  legend->SetTextSize(0.025);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  // graph-dependent attributes, order: raw, msd, msd_btag
  // vector<int> linewidths = {2, 2, 2};
  // vector<int> linestyles = {3, 1, 2};
  // vector<int> linecolors = {kBlue, kRed, kMagenta};
  // vector<string> legends = {"No #it{m}_{SD} window nor b tagging", "105 < #it{m}_{SD} [GeV] < 210", "#it{m}_{SD} as above + subjet b tag"};
  int iterator = 0;
  for(TGraph & graph : graphs) {
    // min0max1(graph);
    mg->Add(&graph);
    // graph->SetLineWidth(linewidths.at(++iterator));
    // graph.SetLineWidth(2);
    // graph->SetLineStyle(linestyles.at(iterator));
    // graph->SetLineColor(linecolors.at(iterator));
    stringstream wp_value;
    wp_value << setprecision(3) << wps_for_legend.at(iterator);
    const string legend_text = legends.at(iterator) + ", WP(#varepsilon_{B} = 3\%) = "+wp_value.str();
    legend->AddEntry(&graph, legend_text.c_str(), "l");
    iterator++;
  }
  auto wp_graph = new_null_graph();
  wp_graph->SetMarkerColor(kBlack);
  wp_graph->SetMarkerSize(1.2);
  wp_graph->SetMarkerStyle(24);
  mg->Add(wp_graph);
  legend->AddEntry(wp_graph, "Working points (WPs)");


  mg->Draw("ac");
  mg->SetTitle("");
  mg->GetHistogram()->GetXaxis()->SetRangeUser(0., 1.);
  double maximum_saved = mg->GetHistogram()->GetMaximum();
  double minimum_saved = mg->GetHistogram()->GetMinimum();
  mg->GetHistogram()->SetMaximum(1.);
  if(!log_y) mg->GetHistogram()->SetMinimum(0.);
  else mg->GetHistogram()->SetMinimum(0.0002);
  if(is_eff_bkg || is_eff_sig) mg->GetHistogram()->GetXaxis()->SetTitle(tagger.scan_var_tlatex_axis.c_str());
  else if(is_roc) mg->GetHistogram()->GetXaxis()->SetTitle("True positive rate, #varepsilon_{S}");
  if(is_eff_bkg || is_roc) mg->GetHistogram()->GetYaxis()->SetTitle("False positive rate, #varepsilon_{B}");
  if(is_eff_sig) mg->GetHistogram()->GetYaxis()->SetTitle("True positive rate, #varepsilon_{S}");
  mg->GetHistogram()->GetXaxis()->SetTitleOffset(1.3);
  mg->GetHistogram()->GetXaxis()->CenterTitle();
  mg->GetHistogram()->GetYaxis()->SetTitleOffset(1.5);
  mg->GetHistogram()->GetYaxis()->CenterTitle();
  // mg->GetHistogram()->GetXaxis()->SetNdivisions(405);
  if(!log_y) mg->GetHistogram()->GetYaxis()->SetNdivisions(405);
  gStyle->SetLineWidth(1);

  legend->Draw();

  TLatex *text_top_left = new TLatex(margin_l, 1-(margin_t-0.01), "Ultra Legacy");
  text_top_left->SetTextAlign(11); // left bottom aligned
  text_top_left->SetTextFont(42);
  text_top_left->SetTextSize(0.035);
  text_top_left->SetNDC();
  text_top_left->Draw();

  TLatex *algo_label = new TLatex(coord->ConvertGraphXToPadX(1-(0.05*c->GetWh()/c->GetWw())), coord->ConvertGraphYToPadY(0.95), (kProbeJetAlgos.at(tagger.probejetalgo).name+" PUPPI").c_str());
  algo_label->SetTextAlign(33); // right top
  algo_label->SetTextFont(62);
  algo_label->SetTextSize(0.05); // 0.05
  algo_label->SetTextColor(kGray+1);
  algo_label->SetNDC();
  algo_label->Draw();

  const double eta_pt_text_x = coord->ConvertGraphXToPadX(1-(0.05*c->GetWh()/c->GetWw()));
  double pt_text_y = coord->ConvertGraphYToPadY(0.867);
  double eta_text_y = coord->ConvertGraphYToPadY(0.801);
  double tagrule_text_y = coord->ConvertGraphYToPadY(0.735);
  if(!plot_at_bottom) {
    pt_text_y = coord->ConvertGraphYToPadY((0.05*c->GetWh()/c->GetWw()));
    eta_text_y = coord->ConvertGraphYToPadY((0.05*c->GetWh()/c->GetWw()) + 0.066);
    tagrule_text_y = coord->ConvertGraphYToPadY((0.05*c->GetWh()/c->GetWw()) + 0.066*2);
  }

  TLatex *eta_text = new TLatex(eta_pt_text_x, tagrule_text_y, "|#eta| < 2.5");
  eta_text->SetTextAlign(plot_at_bottom ? 33 : 31); // right top
  eta_text->SetTextFont(42);
  eta_text->SetTextSize(0.035);
  eta_text->SetNDC();
  eta_text->Draw();

  TLatex *pt_text = new TLatex(eta_pt_text_x, eta_text_y, pt_bin.text.c_str());
  pt_text->SetTextAlign(plot_at_bottom ? 33 : 31); // right top
  pt_text->SetTextFont(42);
  pt_text->SetTextSize(0.035);
  pt_text->SetNDC();
  pt_text->Draw();

  TLatex *tagrule_text = new TLatex(eta_pt_text_x, pt_text_y, tagger.legend_base.c_str());
  tagrule_text->SetTextAlign(plot_at_bottom ? 33 : 31); // right top
  tagrule_text->SetTextFont(42);
  tagrule_text->SetTextSize(0.035);
  tagrule_text->SetNDC();
  tagrule_text->Draw();


  // string string_text_top_right = "t#bar{t} #rightarrow all-jets @ #sqrt{#it{s}} = 13 TeV";
  string string_text_top_right = "#sqrt{#it{s}} = 13 TeV";
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

  if(!plot_at_bottom) {
    TText *prelim = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.87), "Simulation, Private Work");
    prelim->SetTextAlign(13); // left top
    prelim->SetTextFont(52);
    prelim->SetTextSize(0.035);
    prelim->SetNDC();
    prelim->Draw();
  }
  else {
    TText *prelim = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.87), "Simulation");
    prelim->SetTextAlign(13); // left top
    prelim->SetTextFont(52);
    prelim->SetTextSize(0.035);
    prelim->SetNDC();
    // prelim->SetTextColor(kWhite);
    prelim->Draw();

    TText *prelim2 = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.804), "Private Work");
    prelim2->SetTextAlign(13); // left top
    prelim2->SetTextFont(52);
    prelim2->SetTextSize(0.035);
    prelim2->SetNDC();
    // prelim2->SetTextColor(kWhite);
    prelim2->Draw();
  }


  for(const auto & marker : markers) {
    marker->Draw();
  }

  // HACK:
  TLine *wp_line = new TLine(0, 0.03, 1, 0.03);
  wp_line->SetLineColor(kGray+1);
  wp_line->Draw();
  TLatex *wp_text = new TLatex(0.95, 0.035, "3\%");
  wp_text->SetTextAlign(31); // left top
  wp_text->SetTextFont(42);
  wp_text->SetTextSize(0.035);
  wp_text->SetTextColor(kGray+1);
  // wp_text->SetNDC();
  wp_text->Draw();

  // https://root-forum.cern.ch/t/inconsistent-tick-length/18563/8
  c->Update();
  double tickScaleX = (c->GetUxmax() - c->GetUxmin()) / (c->GetX2() - c->GetX1()) * (c->GetWh() * c->GetAbsHNDC());
  double tickScaleY = (c->GetUymax() - c->GetUymin()) / (c->GetY2() - c->GetY1()) * (c->GetWw() * c->GetAbsWNDC());
  mg->GetHistogram()->GetYaxis()->SetTickLength(c->GetWw() * tick_length / tickScaleY);
  mg->GetHistogram()->GetXaxis()->SetTickLength(c->GetWh() * tick_length / tickScaleX);

  gPad->Update();
  gPad->RedrawAxis();

  // Save to disk
  string plotName = "plot__"+tagger.name_base+"__"+pt_bin.name+"__"+graph_base_name+"__";
  // for(bool do_x : {do_raw, do_msd, do_msd_btag}) { do_x ? plotName += "1" : plotName += "0"; }
  plotName += (log_y ? string("log") : string("lin"))+"__allYears.pdf";
  // string plotPath = infileBasePath+plotName;
  TString plotPath = (TString(infilePath)).ReplaceAll("graphs.root", "")+plotName;
  c->SaveAs(plotPath.Data());
  delete c;
}



// void plots(const string & year) {
void plots(const string & arg_year, const string & arg_tagger) {

  // const vector<Year> years = { Year::isUL16preVFP, Year::isUL16postVFP, Year::isUL17, Year::isUL18 };
  // const Year year = Year::isUL16preVFP;
  Year year = Year::isUL16preVFP;

  if(arg_year == "UL16preVFP") year = Year::isUL16preVFP;
  else if(arg_year == "UL16postVFP") year = Year::isUL16postVFP;
  else if(arg_year == "UL17") year = Year::isUL17;
  else if(arg_year == "UL18") year = Year::isUL18;

  vector<PtBin> pt_bins_vec;
  pt_bins_vec.push_back(PtBin{"pt_300to400", "300 < #it{p}_{T}/GeV < 400", kOrange-3, 1});
  pt_bins_vec.push_back(PtBin{"pt_400toInf", "#it{p}_{T} > 400 GeV", kCyan, 1});
  pt_bins_vec.push_back(PtBin{"pt_400to480", "400 < #it{p}_{T}/GeV < 480", kBlue, 2});
  pt_bins_vec.push_back(PtBin{"pt_480to600", "480 < #it{p}_{T}/GeV < 600", kBlue, 3});
  pt_bins_vec.push_back(PtBin{"pt_600toInf", "#it{p}_{T} > 600 GeV", kBlue, 4});
  pt_bins_vec.push_back(PtBin{"pt_300toInf", "#it{p}_{T} > 300 GeV"});
  pt_bins_vec.push_back(PtBin{"pt_1000toInf", "#it{p}_{T} > 1000 GeV"});

  pt_bins_vec.push_back(PtBin{"pt_300to350", "300 < #it{p}_{T}/GeV < 350", kOrange-3, 5});
  pt_bins_vec.push_back(PtBin{"pt_350to400", "350 < #it{p}_{T}/GeV < 400", kOrange-3, 6});
  pt_bins_vec.push_back(PtBin{"pt_200toInf", "#it{p}_{T} > 200 GeV", kOrange-3, 7});
  pt_bins_vec.push_back(PtBin{"pt_200to250", "200 < #it{p}_{T}/GeV < 250", kOrange-3, 8});
  pt_bins_vec.push_back(PtBin{"pt_250to300", "250 < #it{p}_{T}/GeV < 300", kOrange-3, 9});

  map<string, PtBin> pt_bins;
  for(const auto & pt_bin : pt_bins_vec) {
    pt_bins[pt_bin.name] = pt_bin;
  }

  const vector<string> graph_base_names = {"eff_bkg", "eff_sig", "roc"};

  // access here the target background efficiencies (the code is adaptive - you can freely change these values here and all WPs will correctly be recalculated)
  // const vector<float> working_point_mistag_rates = {0.001, 0.003, 0.01, 0.03, 0.1}; // old eB values used for the working points, see https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetTopTagging#Previous_working_points
  // const vector<float> working_point_mistag_rates = {0.001, 0.005, 0.01, 0.025, 0.05}; // new WPs synchronized with DeepAK8, see https://twiki.cern.ch/twiki/bin/viewauth/CMS/DeepAK8Tagging2018WPsSFs
  // for(int i = 0; i < pt_bins.size(); i++) cout << pt_bins.at(i).name << endl;
  // const vector<WorkingPoint> wps_reference = calculate_working_points(year, pt_bins.at("pt_400toInf"), working_point_mistag_rates);

  const string ref_pt_ak8_t__tau = "pt_400toInf";
  const Tagger tagger_template_ak8_t__tau = {
    .probejetalgo = ProbeJetAlgo::isAK8,
    .name_base = "ak8_t__tau",
    .name_variants = {"ak8_t_btagDCSV__tau", "ak8_t_btagDJet__tau", "ak8_t_incl__tau"},
    .target_effs_bkg = {0.001, 0.005, 0.01, 0.025, 0.05},
    .linestyle_variants = {2, 2, 3},
    .linecolor_variants = {kCyan, kBlue, kBlack},
    .legend_base = "105 < #it{m}_{SD}/GeV < 210",
    .legend_variants = {"#it{m}_{SD} + loose DeepCSV b tag", "#it{m}_{SD} + loose DeepJet b tag", "Inclusive jet collection"},
    .scan_var_tlatex = "#tau_{3}/#tau_{2}",
    .scan_var_tlatex_axis = "#tau_{3}/#tau_{2} upper limit",
  };

  const string ref_pt_hotvr_t__tau = "pt_400toInf";
  const Tagger tagger_template_hotvr_t__tau = {
    .probejetalgo = ProbeJetAlgo::isHOTVR,
    .name_base = "hotvr_t__tau",
    .name_variants = {"hotvr_t_incl__tau"},
    .target_effs_bkg = {0.03},
    .linestyle_variants = {3},
    .linecolor_variants = {kBlack},
    .legend_base = "Default substructure cuts",
    .legend_variants = {"Inclusive jet collection"},
    .scan_var_tlatex = "#tau_{3}/#tau_{2}",
    .scan_var_tlatex_axis = "#tau_{3}/#tau_{2} upper limit",
  };

  const string ref_pt_ak8_w__partnet = "pt_200toInf";
  const Tagger tagger_template_ak8_w__partnet = {
    .probejetalgo = ProbeJetAlgo::isAK8,
    .name_base = "ak8_w__partnet",
    // .name_variants = {"ak8_w_incl__partnet"},
    .target_effs_bkg = {0.03},
    .linestyle_variants = {3},
    .linecolor_variants = {kBlack},
    .legend_base = "65 < #it{m}_{SD}/GeV < 105",
    .legend_variants = {"Inclusive jet collection"},
    .scan_var_tlatex = "ParticleNet \"WvsQCD\" score",
    .scan_var_tlatex_axis = "ParticleNet \"WvsQCD\" score lower limit",
  };

  string ref_pt;
  Tagger tagger_template;
  if(arg_tagger == "ak8_t__tau") {
    tagger_template = tagger_template_ak8_t__tau;
    ref_pt = ref_pt_ak8_t__tau;
  }
  else if(arg_tagger == "hotvr_t__tau") {
    tagger_template = tagger_template_hotvr_t__tau;
    ref_pt = ref_pt_hotvr_t__tau;
  }
  else if(arg_tagger == "ak8_w__partnet") {
    tagger_template = tagger_template_ak8_w__partnet;
    ref_pt = ref_pt_ak8_w__partnet;
  }

  Tagger reference_tagger = tagger_template;
  Tagger tagger_pt_480to600 = tagger_template;
  calculate_working_points(year, reference_tagger, pt_bins.at(ref_pt), true);
  get_reference_working_points(year, tagger_pt_480to600, pt_bins.at("pt_480to600"), reference_tagger, false);
  Tagger tagger_pt_480to600_DPnote = init_DPnote_Tagger(tagger_template);

  // Tagger tagger = {
  //   .name_base = "ak8_w__partnet",
  //   // .name_variants = {"ak8_t_incl__tau", "ak8_t_btag__tau"},
  //   // .name_variants = {"ak8_t_btag__tau"},
  //   .target_effs_bkg = {0.005, 0.01, 0.05, 0.03},
  // };
  //
  // Tagger reference_tagger = tagger;
  // calculate_working_points(year, tagger, pt_bins.at("pt_300toInf"), true);

  for(const auto & pt_bin : pt_bins) {
    Tagger temp_tagger = tagger_template; // need to use new temporary copy of tagger for each pt_bin since get_reference_working_points will append WP vectors in Tagger struct for each pt_bin
    cout << "Working on " << pt_bin.second.name << endl;
    get_reference_working_points(year, temp_tagger, pt_bin.second, reference_tagger, false);
    for(const auto & graph_base_name : graph_base_names) {
      do_plot(year, temp_tagger, pt_bin.second, graph_base_name, true);
      do_plot(year, temp_tagger, pt_bin.second, graph_base_name, false);
      do_plot_all_years(temp_tagger, pt_bin.second, graph_base_name, true);
      do_plot_all_years(temp_tagger, pt_bin.second, graph_base_name, false);
    }
  }

  vector<PtBin> pt_bins_summary_plots;
  pt_bins_summary_plots.push_back(pt_bins.at("pt_400toInf"));
  // pt_bins_summary_plots.push_back(pt_bins.at("pt_300to400"));
  pt_bins_summary_plots.push_back(pt_bins.at("pt_300to350"));
  pt_bins_summary_plots.push_back(pt_bins.at("pt_350to400"));
  pt_bins_summary_plots.push_back(pt_bins.at("pt_400to480"));
  pt_bins_summary_plots.push_back(pt_bins.at("pt_480to600"));
  pt_bins_summary_plots.push_back(pt_bins.at("pt_600toInf"));

  for(const auto & graph_base_name : graph_base_names) {
    do_summary_plot(year, reference_tagger, reference_tagger.name_base, pt_bins_summary_plots, graph_base_name, true, tagger_pt_480to600, tagger_pt_480to600_DPnote);
    do_summary_plot(year, reference_tagger, reference_tagger.name_base, pt_bins_summary_plots, graph_base_name, false, tagger_pt_480to600, tagger_pt_480to600_DPnote);
    for(const auto & variant : reference_tagger.name_variants) {
      do_summary_plot(year, reference_tagger, variant, pt_bins_summary_plots, graph_base_name, true, tagger_pt_480to600, tagger_pt_480to600_DPnote);
      do_summary_plot(year, reference_tagger, variant, pt_bins_summary_plots, graph_base_name, false, tagger_pt_480to600, tagger_pt_480to600_DPnote);
    }
  }
}
