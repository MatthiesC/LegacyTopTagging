#include <iomanip>
#include <sstream>

typedef struct {
  string name;
  string text;
  int summary_plot_linecolor;
  int summary_plot_linestyle;
  int summary_plot_linewidth=2;
} PtBin;

typedef struct {
  float target_eff_qcd;
  float real_eff_qcd;
  float real_eff_ttbar;
  float tau32cut;
  int tau32cut_i;
  float real_eff_qcd_with_btag;
  float real_eff_ttbar_with_btag;
} WorkingPoint;

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

WorkingPoint init_working_point(float eff, const TGraph * graph_eff_qcd, const TGraph * graph_eff_ttbar, const TGraph * graph_eff_qcd_btag, const TGraph * graph_eff_ttbar_btag) {
  WorkingPoint result = { .target_eff_qcd = eff };
  int index_tau32cut;
  for(int i = graph_eff_qcd->GetN()-1; i >= 0; i--) {
    double x, y;
    graph_eff_qcd->GetPoint(i, x, y);
    if(y < eff) {
      result.real_eff_qcd = y;
      result.tau32cut = x;
      index_tau32cut = i;
      if(i != graph_eff_qcd->GetN()-1) {
        double x_before, y_before;
        graph_eff_qcd->GetPoint(i+1, x_before, y_before);
        if(abs(y_before - eff) < abs(y - eff)) { result.real_eff_qcd = y_before; result.tau32cut = x_before; index_tau32cut = i+1; }
      }
      break;
    }
  }
  double x, y;
  graph_eff_ttbar->GetPoint(index_tau32cut, x, y);
  result.real_eff_ttbar = y;
  graph_eff_qcd_btag->GetPoint(index_tau32cut, x, y);
  result.real_eff_qcd_with_btag = y;
  graph_eff_ttbar_btag->GetPoint(index_tau32cut, x, y);
  result.real_eff_ttbar_with_btag = y;
  result.tau32cut_i = index_tau32cut;
  return result;
}


vector<WorkingPoint> calculate_working_points(const string & year, const PtBin & pt_bin, const vector<float> & target_effs_qcd, const bool print=false) {
  string infileBasePath = (string)getenv("CMSSW_BASE")+"/src/UHH2/LegacyTopTagging/output/WorkingPointStudy/"+year+"/workdir_npy/"+pt_bin.name+"/";
  string infileName = "root_"+pt_bin.name+".root";
  string infilePath = infileBasePath+infileName;

  TFile* infile = TFile::Open(infilePath.c_str(), "READ");
  const TGraph* graph_eff_qcd_msd = (TGraph*)infile->Get("eff_qcd_msd");
  const TGraph* graph_eff_qcd_msd_btag = (TGraph*)infile->Get("eff_qcd_msd_btag");
  const TGraph* graph_eff_ttbar_msd = (TGraph*)infile->Get("eff_ttbar_msd");
  const TGraph* graph_eff_ttbar_msd_btag = (TGraph*)infile->Get("eff_ttbar_msd_btag");

  vector<WorkingPoint> working_points;
  for(auto wp : target_effs_qcd) {
    working_points.push_back(init_working_point(wp, graph_eff_qcd_msd, graph_eff_ttbar_msd, graph_eff_qcd_msd_btag, graph_eff_ttbar_msd_btag));
  }
  if(print) for(auto wp : working_points) {
    cout << ">>> WP derived from this pt bin <<<" << endl;
    cout << "Target QCD eff: " << wp.target_eff_qcd << endl;
    cout << "Actual QCD eff: " << wp.real_eff_qcd << endl;
    cout << "Actual tt eff:  " << wp.real_eff_ttbar << endl;
    cout << "tau32 cut:      " << wp.tau32cut << endl;
    cout << "Actual QCD eff with b-tag: " << wp.real_eff_qcd_with_btag << endl;
    cout << "Actual tt eff with b-tag:  " << wp.real_eff_ttbar_with_btag << endl;
  };
  return working_points;
}


vector<WorkingPoint> get_reference_working_points(const string & year, const PtBin & pt_bin, const vector<WorkingPoint> & wps_reference, const bool print=false) {
  string infileBasePath = (string)getenv("CMSSW_BASE")+"/src/UHH2/LegacyTopTagging/output/WorkingPointStudy/"+year+"/workdir_npy/"+pt_bin.name+"/";
  string infileName = "root_"+pt_bin.name+".root";
  string infilePath = infileBasePath+infileName;

  TFile* infile = TFile::Open(infilePath.c_str(), "READ");
  const TGraph* graph_eff_qcd_msd = (TGraph*)infile->Get("eff_qcd_msd");
  const TGraph* graph_eff_qcd_msd_btag = (TGraph*)infile->Get("eff_qcd_msd_btag");
  const TGraph* graph_eff_ttbar_msd = (TGraph*)infile->Get("eff_ttbar_msd");
  const TGraph* graph_eff_ttbar_msd_btag = (TGraph*)infile->Get("eff_ttbar_msd_btag");

  vector<WorkingPoint> working_points;
  double x=0, y=0;
  for(int i = 0; i < wps_reference.size(); i++) {
    WorkingPoint wp_ref = wps_reference.at(i);
    WorkingPoint wp_ref_this_bin;
    wp_ref_this_bin.target_eff_qcd = wp_ref.target_eff_qcd;
    wp_ref_this_bin.tau32cut = wp_ref.tau32cut;
    wp_ref_this_bin.tau32cut_i = wp_ref.tau32cut_i;
    graph_eff_qcd_msd->GetPoint(wp_ref.tau32cut_i, x, y); wp_ref_this_bin.real_eff_qcd = y;
    graph_eff_qcd_msd_btag->GetPoint(wp_ref.tau32cut_i, x, y); wp_ref_this_bin.real_eff_qcd_with_btag = y;
    graph_eff_ttbar_msd->GetPoint(wp_ref.tau32cut_i, x, y); wp_ref_this_bin.real_eff_ttbar = y;
    graph_eff_ttbar_msd_btag->GetPoint(wp_ref.tau32cut_i, x, y); wp_ref_this_bin.real_eff_ttbar_with_btag = y;
    working_points.push_back(wp_ref_this_bin);
  }
  if(print) for(auto wp : working_points) {
    cout << ">>> WP derived from reference bin <<<" << endl;
    cout << "Target QCD eff: " << wp.target_eff_qcd << endl;
    cout << "Actual QCD eff: " << wp.real_eff_qcd << endl;
    cout << "Actual tt eff:  " << wp.real_eff_ttbar << endl;
    cout << "tau32 cut:      " << wp.tau32cut << endl;
    cout << "Actual QCD eff with b-tag: " << wp.real_eff_qcd_with_btag << endl;
    cout << "Actual tt eff with b-tag:  " << wp.real_eff_ttbar_with_btag << endl;
  };
  return working_points;
}


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


void do_plot(const string & year, const PtBin & pt_bin, const string & graph_base_name, const bool log_y, const vector<WorkingPoint> & wps_this_bin, const vector<WorkingPoint> & wps_this_bin_ref, const bool do_raw = true, const bool do_msd = true, const bool do_msd_btag = true) {
  if(!(do_raw || do_msd || do_msd_btag)) cout << "Nothing to plot." << endl;
  bool is_eff_qcd = graph_base_name == "eff_qcd";
  bool is_eff_ttbar = graph_base_name == "eff_ttbar";
  bool is_roc = graph_base_name == "roc";

  string infileBasePath = (string)getenv("CMSSW_BASE")+"/src/UHH2/LegacyTopTagging/output/WorkingPointStudy/"+year+"/workdir_npy/"+pt_bin.name+"/";
  string infileName = "root_"+pt_bin.name+".root";
  string infilePath = infileBasePath+infileName;

  TFile* infile = TFile::Open(infilePath.c_str(), "READ");
  vector<TGraph*> graphs;
  if(do_raw) graphs.push_back((TGraph*)infile->Get(graph_base_name.c_str()));
  if(do_msd) graphs.push_back((TGraph*)infile->Get((graph_base_name+"_msd").c_str()));
  if(do_msd_btag) graphs.push_back((TGraph*)infile->Get((graph_base_name+"_msd_btag").c_str()));

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
  if((is_eff_ttbar || is_roc) && log_y) plot_at_bottom = true;

  float legend_x1 = coord->ConvertGraphXToPadX(plot_at_bottom ? 0.45 : 0.05);
  float legend_y1 = coord->ConvertGraphYToPadY(plot_at_bottom ? 0.05 : 0.47);
  float legend_x2 = coord->ConvertGraphXToPadX(plot_at_bottom ? 0.95 : 0.55);
  float legend_y2 = coord->ConvertGraphYToPadY(plot_at_bottom ? 0.35 : 0.77);
  TLegend *legend = new TLegend(legend_x1, legend_y1, legend_x2, legend_y2); // x1 y1 x2 y2
  legend->SetHeader(("#bf{Scan of #tau_{3}/#tau_{2} for "+pt_bin.text+"}").c_str());
  legend->SetTextSize(0.025);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  // graph-dependent attributes, order: raw, msd, msd_btag
  vector<int> linewidths = {2, 2, 2};
  vector<int> linestyles = {3, 1, 2};
  vector<int> linecolors = {kBlue, kRed, kMagenta};
  vector<string> legends = {"No #it{m}_{SD} window nor b tagging", "105 < #it{m}_{SD} [GeV] < 210", "#it{m}_{SD} as above + subjet b tag"};
  int iterator = -1;
  for(auto graph : graphs) {
    // min0max1(graph);
    mg->Add(graph);
    graph->SetLineWidth(linewidths.at(++iterator));
    graph->SetLineStyle(linestyles.at(iterator));
    // graph->SetLineColor(linecolors.at(iterator));
    legend->AddEntry(graph, legends.at(iterator).c_str(), "l");
  }

  mg->Draw("ac");
  mg->SetTitle("");
  mg->GetHistogram()->GetXaxis()->SetRangeUser(0., 1.);
  double maximum_saved = mg->GetHistogram()->GetMaximum();
  double minimum_saved = mg->GetHistogram()->GetMinimum();
  mg->GetHistogram()->SetMaximum(1.);
  if(!log_y) mg->GetHistogram()->SetMinimum(0.);
  else mg->GetHistogram()->SetMinimum(0.0002);
  if(is_eff_qcd || is_eff_ttbar) mg->GetHistogram()->GetXaxis()->SetTitle("#tau_{3}/#tau_{2} upper limit");
  else if(is_roc) mg->GetHistogram()->GetXaxis()->SetTitle("#varepsilon_{S}");
  if(is_eff_qcd || is_roc) mg->GetHistogram()->GetYaxis()->SetTitle("#varepsilon_{B}");
  if(is_eff_ttbar) mg->GetHistogram()->GetYaxis()->SetTitle("#varepsilon_{S}");
  mg->GetHistogram()->GetXaxis()->SetTitleOffset(1.2);
  mg->GetHistogram()->GetYaxis()->SetTitleOffset(1.2);
  // mg->GetHistogram()->GetXaxis()->SetNdivisions(405);
  if(!log_y) mg->GetHistogram()->GetYaxis()->SetNdivisions(405);
  gStyle->SetLineWidth(1);

  legend->Draw();

  TText *text_top_left = new TText(margin_l, 1-(margin_t-0.01), "AK8 PUPPI(v15)");
  text_top_left->SetTextAlign(11); // left bottom aligned
  text_top_left->SetTextFont(42);
  text_top_left->SetTextSize(0.035);
  text_top_left->SetNDC();
  text_top_left->Draw();

  string string_text_top_right = year + " (CMSSW 10.6.X)";
  TText *text_top_right = new TText(1-margin_r, 1-(margin_t-0.01), string_text_top_right.c_str());
  text_top_right->SetTextAlign(31); // right bottom aligned
  text_top_right->SetTextFont(42);
  text_top_right->SetTextSize(0.035);
  text_top_right->SetNDC();
  text_top_right->Draw();

  TText *cms = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.95), "CMS");
  cms->SetTextAlign(13); // left top
  cms->SetTextFont(62);
  cms->SetTextSize(0.05);
  cms->SetNDC();
  cms->Draw();

  TText *prelim = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.87), "Simulation, Work in Progress");
  prelim->SetTextAlign(13); // left top
  prelim->SetTextFont(52);
  prelim->SetTextSize(0.035);
  prelim->SetNDC();
  prelim->Draw();

  // WP markers
  for(int i = 0; i < wps_this_bin.size(); i++) {
    double x=0, y=0;

    WorkingPoint wp = wps_this_bin.at(i);
    WorkingPoint wp_ref = wps_this_bin_ref.at(i);

    if(is_roc) { x=wp.real_eff_ttbar; y=wp.real_eff_qcd; }
    else if(is_eff_qcd) { x=wp.tau32cut; y=wp.real_eff_qcd; }
    else if(is_eff_ttbar) { x=wp.tau32cut; y=wp.real_eff_ttbar; }
    TMarker *marker1 = new TMarker(x, y, 24);
    marker1->SetMarkerColor(kRed);
    marker1->SetMarkerSize(1.2);
    marker1->Draw();

    if(is_roc) { x=wp_ref.real_eff_ttbar; y=wp_ref.real_eff_qcd; }
    else if(is_eff_qcd) { x=wp_ref.tau32cut; y=wp_ref.real_eff_qcd; }
    else if(is_eff_ttbar) { x=wp_ref.tau32cut; y=wp_ref.real_eff_ttbar; }
    TMarker *marker2 = new TMarker(x, y, 3);
    marker2->SetMarkerColor(kBlue);
    marker2->SetMarkerSize(1.2);
    marker2->Draw();

    if(is_roc) { x=wp.real_eff_ttbar_with_btag; y=wp.real_eff_qcd_with_btag; }
    else if(is_eff_qcd) { x=wp.tau32cut; y=wp.real_eff_qcd_with_btag; }
    else if(is_eff_ttbar) { x=wp.tau32cut; y=wp.real_eff_ttbar_with_btag; }
    TMarker *marker3 = new TMarker(x, y, 24);
    marker3->SetMarkerColor(kOrange-3);
    marker3->SetMarkerSize(1.2);
    marker3->Draw();

    if(is_roc) { x=wp_ref.real_eff_ttbar_with_btag; y=wp_ref.real_eff_qcd_with_btag; }
    else if(is_eff_qcd) { x=wp_ref.tau32cut; y=wp_ref.real_eff_qcd_with_btag; }
    else if(is_eff_ttbar) { x=wp_ref.tau32cut; y=wp_ref.real_eff_ttbar_with_btag; }
    TMarker *marker4 = new TMarker(x, y, 3);
    marker4->SetMarkerColor(kCyan);
    marker4->SetMarkerSize(1.2);
    marker4->Draw();
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
  string plotName = "plot_"+pt_bin.name+"__"+graph_base_name+"__";
  for(bool do_x : {do_raw, do_msd, do_msd_btag}) { do_x ? plotName += "1" : plotName += "0"; }
  plotName += (log_y ? string("_log") : string("_lin"))+".pdf";
  string plotPath = infileBasePath+plotName;
  c->SaveAs(plotPath.c_str());
  delete c;
}


void do_summary_plot(const string & year, const string & tag_type, const vector<PtBin> & pt_bins, const string & graph_base_name, const bool log_y, const vector<WorkingPoint> & wps_ref_pt480to600, const vector<WorkingPoint> & wps_dpnote_pt480to600) {
  bool is_eff_qcd = graph_base_name == "eff_qcd";
  bool is_eff_ttbar = graph_base_name == "eff_ttbar";
  bool is_roc = graph_base_name == "roc";

  bool is_raw = tag_type == "raw";
  bool is_msd = tag_type == "msd";
  bool is_msd_btag = tag_type == "msd_btag";

  string infileBaseBasePath = (string)getenv("CMSSW_BASE")+"/src/UHH2/LegacyTopTagging/output/WorkingPointStudy/"+year+"/workdir_npy/";
  vector<TGraph*> graphs;
  for(auto pt_bin : pt_bins) {
    string infileBasePath = infileBaseBasePath+pt_bin.name+"/";
    string infileName = "root_"+pt_bin.name+".root";
    string infilePath = infileBasePath+infileName;
    TFile* infile = TFile::Open(infilePath.c_str(), "READ");
    graphs.push_back((TGraph*)infile->Get((graph_base_name+(is_raw ? "" : "_"+tag_type)).c_str()));
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
  if((is_eff_ttbar || is_roc) && log_y) plot_at_bottom = true;

  float legend_x1 = coord->ConvertGraphXToPadX(plot_at_bottom ? 0.45 : 0.05);
  float legend_y1 = coord->ConvertGraphYToPadY(plot_at_bottom ? 0.05 : 0.32); //0.47
  float legend_x2 = coord->ConvertGraphXToPadX(plot_at_bottom ? 0.95 : 0.55);
  float legend_y2 = coord->ConvertGraphYToPadY(plot_at_bottom ? 0.50 : 0.77); //0.35
  TLegend *legend = new TLegend(legend_x1, legend_y1, legend_x2, legend_y2); // x1 y1 x2 y2
  string legend_title = "";
  if(is_raw) legend_title = "No #it{m}_{SD} window nor b tagging";
  else if(is_msd) legend_title = "105 < #it{m}_{SD} [GeV] < 210";
  else if(is_msd_btag) legend_title = "105 < #it{m}_{SD} [GeV] < 210 + subjet b tag";
  legend->SetHeader(("#bf{"+legend_title+"}").c_str());
  legend->SetTextSize(0.025);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  for(int i = 0; i < graphs.size(); i++) {
    mg->Add(graphs.at(i));
    graphs.at(i)->SetLineWidth(pt_bins.at(i).summary_plot_linewidth);
    graphs.at(i)->SetLineStyle(pt_bins.at(i).summary_plot_linestyle);
    graphs.at(i)->SetLineColor(pt_bins.at(i).summary_plot_linecolor);
    legend->AddEntry(graphs.at(i), pt_bins.at(i).text.c_str(), "l");
  }
  TGraph * null_graph = new_null_graph();
  mg->Add(null_graph);

  mg->Draw("ac");
  mg->SetTitle("");
  mg->GetHistogram()->GetXaxis()->SetRangeUser(0., 1.);
  double maximum_saved = mg->GetHistogram()->GetMaximum();
  double minimum_saved = mg->GetHistogram()->GetMinimum();
  mg->GetHistogram()->SetMaximum(1.);
  if(!log_y) mg->GetHistogram()->SetMinimum(0.);
  else mg->GetHistogram()->SetMinimum(0.0002);
  if(is_eff_qcd || is_eff_ttbar) mg->GetHistogram()->GetXaxis()->SetTitle("#tau_{3}/#tau_{2} upper limit");
  else if(is_roc) mg->GetHistogram()->GetXaxis()->SetTitle("#varepsilon_{S}");
  if(is_eff_qcd || is_roc) mg->GetHistogram()->GetYaxis()->SetTitle("#varepsilon_{B}");
  if(is_eff_ttbar) mg->GetHistogram()->GetYaxis()->SetTitle("#varepsilon_{S}");
  mg->GetHistogram()->GetXaxis()->SetTitleOffset(1.2);
  mg->GetHistogram()->GetYaxis()->SetTitleOffset(1.2);
  // mg->GetHistogram()->GetXaxis()->SetNdivisions(405);
  if(!log_y) mg->GetHistogram()->GetYaxis()->SetNdivisions(405);
  gStyle->SetLineWidth(1);

  legend->Draw();

  TText *text_top_left = new TText(margin_l, 1-(margin_t-0.01), "AK8 PUPPI(v15)");
  text_top_left->SetTextAlign(11); // left bottom aligned
  text_top_left->SetTextFont(42);
  text_top_left->SetTextSize(0.035);
  text_top_left->SetNDC();
  text_top_left->Draw();

  string string_text_top_right = year + " (CMSSW 10.6.X)";
  TText *text_top_right = new TText(1-margin_r, 1-(margin_t-0.01), string_text_top_right.c_str());
  text_top_right->SetTextAlign(31); // right bottom aligned
  text_top_right->SetTextFont(42);
  text_top_right->SetTextSize(0.035);
  text_top_right->SetNDC();
  text_top_right->Draw();

  TText *cms = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.95), "CMS");
  cms->SetTextAlign(13); // left top
  cms->SetTextFont(62);
  cms->SetTextSize(0.05);
  cms->SetNDC();
  cms->Draw();

  TText *prelim = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.87), "Simulation, Work in Progress");
  prelim->SetTextAlign(13); // left top
  prelim->SetTextFont(52);
  prelim->SetTextSize(0.035);
  prelim->SetNDC();
  prelim->Draw();

  // WP markers for Pt480to600
  if(is_msd || is_msd_btag) {
    for(int i = 0; i < wps_ref_pt480to600.size(); i++) {
      double x=0, y=0;

      WorkingPoint wp = wps_ref_pt480to600.at(i);

      if(is_msd) {
        if(is_roc) { x=wp.real_eff_ttbar; y=wp.real_eff_qcd; }
        else if(is_eff_qcd) { x=wp.tau32cut; y=wp.real_eff_qcd; }
        else if(is_eff_ttbar) { x=wp.tau32cut; y=wp.real_eff_ttbar; }
      }
      else if(is_msd_btag) {
        if(is_roc) { x=wp.real_eff_ttbar_with_btag; y=wp.real_eff_qcd_with_btag; }
        else if(is_eff_qcd) { x=wp.tau32cut; y=wp.real_eff_qcd_with_btag; }
        else if(is_eff_ttbar) { x=wp.tau32cut; y=wp.real_eff_ttbar_with_btag; }
      }
      TMarker *marker = new TMarker(x, y, 20);
      marker->SetMarkerColor(kBlack);
      marker->SetMarkerSize(1.2);
      marker->Draw();
    }
    for(int i = 0; i < wps_dpnote_pt480to600.size(); i++) {
      double x=0, y=0;

      WorkingPoint wp = wps_dpnote_pt480to600.at(i);

      if(is_msd) {
        if(is_roc) { x=wp.real_eff_ttbar; y=wp.real_eff_qcd; }
        else if(is_eff_qcd) { x=wp.tau32cut; y=wp.real_eff_qcd; }
        else if(is_eff_ttbar) { x=wp.tau32cut; y=wp.real_eff_ttbar; }
      }
      else if(is_msd_btag) {
        if(is_roc) { x=wp.real_eff_ttbar_with_btag; y=wp.real_eff_qcd_with_btag; }
        else if(is_eff_qcd) { x=wp.tau32cut; y=wp.real_eff_qcd_with_btag; }
        else if(is_eff_ttbar) { x=wp.tau32cut; y=wp.real_eff_ttbar_with_btag; }
      }
      TMarker *marker = new TMarker(x, y, 3);
      marker->SetMarkerColor(kRed);
      marker->SetMarkerSize(1.2);
      marker->Draw();
    }
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
  string plotName = "plot_summary__"+graph_base_name+"__"+tag_type+"__";
  plotName += (log_y ? string("log") : string("lin"))+".pdf";
  string plotPath = infileBaseBasePath+plotName;
  c->SaveAs(plotPath.c_str());
  delete c;
}


string my_format(float value, int prec=2, int factor=100) {
  std::stringstream stream;
  stream << std::fixed << std::setprecision(prec) << value*factor;
  return stream.str();
}


void print_latex_table(const vector<WorkingPoint> & wps) {
  cout << "WPs in LaTeX tabular format:" << endl;
  for(auto wp : wps) {
    cout << "$< " << my_format(wp.tau32cut, 3, 1) << "$ & & " << my_format(wp.real_eff_qcd) << " & " << my_format(wp.real_eff_ttbar) << " & & " << my_format(wp.real_eff_qcd_with_btag) << " & " << my_format(wp.real_eff_ttbar_with_btag) << "\\\\" << endl;
  }
}


vector<WorkingPoint> init_DPnote_WPs() {
  // values taken from CMS DP 2020-025
  vector<WorkingPoint> wps;
  wps.push_back(WorkingPoint{ .tau32cut=0.80, .real_eff_qcd=0.159, .real_eff_ttbar=0.62, .real_eff_qcd_with_btag=0.053, .real_eff_ttbar_with_btag=0.53 });
  wps.push_back(WorkingPoint{ .tau32cut=0.65, .real_eff_qcd=0.051, .real_eff_ttbar=0.49, .real_eff_qcd_with_btag=0.018, .real_eff_ttbar_with_btag=0.43 });
  wps.push_back(WorkingPoint{ .tau32cut=0.54, .real_eff_qcd=0.017, .real_eff_ttbar=0.37, .real_eff_qcd_with_btag=0.006, .real_eff_ttbar_with_btag=0.33 });
  wps.push_back(WorkingPoint{ .tau32cut=0.46, .real_eff_qcd=0.005, .real_eff_ttbar=0.26, .real_eff_qcd_with_btag=0.003, .real_eff_ttbar_with_btag=0.23 });
  wps.push_back(WorkingPoint{ .tau32cut=0.40, .real_eff_qcd=0.002, .real_eff_ttbar=0.17, .real_eff_qcd_with_btag=0.001, .real_eff_ttbar_with_btag=0.16 });
  return wps;
}


void plots(const string & year) {
  vector<PtBin> pt_bins;
  // DO NOT CHANGE THE ORDER !!! .at(i) is used later here to refer to specific pt bins!
  pt_bins.push_back(PtBin{"Pt300to400", "300 < #it{p}_{T}^{jet} [GeV] < 400", kOrange-3, 1});
  pt_bins.push_back(PtBin{"Pt400toInf", "#it{p}_{T}^{jet} [GeV] > 400", kCyan, 1});
  pt_bins.push_back(PtBin{"Pt400to480", "400 < #it{p}_{T}^{jet} [GeV] < 480", kBlue, 2});
  pt_bins.push_back(PtBin{"Pt480to600", "480 < #it{p}_{T}^{jet} [GeV] < 600", kBlue, 3});
  pt_bins.push_back(PtBin{"Pt600toInf", "#it{p}_{T}^{jet} [GeV] > 600", kBlue, 4});
  pt_bins.push_back(PtBin{"Pt300toInf", "#it{p}_{T}^{jet} [GeV] > 300"});
  pt_bins.push_back(PtBin{"Pt1000toInf", "#it{p}_{T}^{jet} [GeV] > 1000"});

  // vector<string> graph_base_names = {"roc"};
  const vector<string> graph_base_names = {"eff_qcd", "eff_ttbar", "roc"};

  // access here the target background efficiencies (the code is adaptive - you can freely change these values here and all WPs will correctly be recalculated)
  // const vector<float> working_point_mistag_rates = {0.001, 0.003, 0.01, 0.03, 0.1}; // old eB values used for the working points, see https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetTopTagging#Previous_working_points
  const vector<float> working_point_mistag_rates = {0.001, 0.005, 0.01, 0.025, 0.05}; // new WPs synchronized with DeepAK8, see https://twiki.cern.ch/twiki/bin/viewauth/CMS/DeepAK8Tagging2018WPsSFs
  // for(int i = 0; i < pt_bins.size(); i++) cout << pt_bins.at(i).name << endl;
  const vector<WorkingPoint> wps_reference = calculate_working_points(year, pt_bins.at(1), working_point_mistag_rates);

  for(const auto & pt_bin : pt_bins) {
    cout << "Working on " << pt_bin.name << endl;
    const vector<WorkingPoint> wps_this_bin = calculate_working_points(year, pt_bin, working_point_mistag_rates);
    const vector<WorkingPoint> wps_this_bin_ref = get_reference_working_points(year, pt_bin, wps_reference);
    print_latex_table(wps_this_bin);
    print_latex_table(wps_this_bin_ref);
    for(const auto & graph_base_name : graph_base_names) {
      do_plot(year, pt_bin, graph_base_name, true, wps_this_bin, wps_this_bin_ref);
      do_plot(year, pt_bin, graph_base_name, false, wps_this_bin, wps_this_bin_ref);
    }
  }

  vector<PtBin> pt_bins_summary_plots;
  pt_bins_summary_plots.push_back(pt_bins.at(1));
  pt_bins_summary_plots.push_back(pt_bins.at(0));
  pt_bins_summary_plots.push_back(pt_bins.at(2));
  pt_bins_summary_plots.push_back(pt_bins.at(3));
  pt_bins_summary_plots.push_back(pt_bins.at(4));
  const vector<string> tag_types = {"raw", "msd", "msd_btag"};
  const vector<WorkingPoint> wps_ref_pt480to600 = get_reference_working_points(year, pt_bins.at(3), wps_reference);
  const vector<WorkingPoint> wps_dpnote_pt480to600 = init_DPnote_WPs();
  cout << "WPs for Pt480to600: this study (tau32cuts from reference bin) vs. CMS DP 2020-025" << endl;
  print_latex_table(wps_ref_pt480to600);
  print_latex_table(wps_dpnote_pt480to600);
  for(const auto & tag_type : tag_types) {
    for(const auto & graph_base_name : graph_base_names) {
      do_summary_plot(year, tag_type, pt_bins_summary_plots, graph_base_name, true, wps_ref_pt480to600, wps_dpnote_pt480to600);
      do_summary_plot(year, tag_type, pt_bins_summary_plots, graph_base_name, false, wps_ref_pt480to600, wps_dpnote_pt480to600);
    }
  }
}
