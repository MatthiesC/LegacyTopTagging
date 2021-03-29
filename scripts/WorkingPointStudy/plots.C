typedef struct {
  string name;
  string text;
} PtBin;

typedef struct {
  float target_eff_qcd;
  float real_eff_qcd;
  float real_eff_ttbar;
  float tau32cut;
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
  return result;
}



void do_plot(const PtBin & pt_bin, const string & graph_base_name, const bool log_y, const bool do_raw = true, const bool do_msd = true, const bool do_msd_btag = true) {
  if(!(do_raw || do_msd || do_msd_btag)) cout << "Nothing to plot." << endl;
  bool is_eff_qcd = graph_base_name == "eff_qcd";
  bool is_eff_ttbar = graph_base_name == "eff_ttbar";
  bool is_roc = graph_base_name == "roc";

  string infileBasePath = "/nfs/dust/cms/user/matthies/uhh2-106X-v1-el7/CMSSW_10_6_13/src/UHH2/LegacyTopTagging/output/WorkingPointStudy/UL17/workdir_npy/"+pt_bin.name+"/";
  string infileName = "root_"+pt_bin.name+".root";
  string infilePath = infileBasePath+infileName;

  TFile* infile = TFile::Open(infilePath.c_str(), "READ");
  vector<TGraph*> graphs;
  if(do_raw) graphs.push_back((TGraph*)infile->Get(graph_base_name.c_str()));
  if(do_msd) graphs.push_back((TGraph*)infile->Get((graph_base_name+"_msd").c_str()));
  if(do_msd_btag) graphs.push_back((TGraph*)infile->Get((graph_base_name+"_msd_btag").c_str()));

  vector<float> working_point_mistag_rates = {0.001, 0.003, 0.01, 0.03, 0.1}; // Access definition of working point mistag rates here!
  vector<int> markerstyles = {47, 22, 34, 23, 20};
  vector<int> markercolors = {kOrange, kRed, kBlue, kGreen, kMagenta};
  const TGraph* graph_eff_qcd_msd = (TGraph*)infile->Get("eff_qcd_msd");
  const TGraph* graph_eff_qcd_msd_btag = (TGraph*)infile->Get("eff_qcd_msd_btag");
  const TGraph* graph_eff_ttbar_msd = (TGraph*)infile->Get("eff_ttbar_msd");
  const TGraph* graph_eff_ttbar_msd_btag = (TGraph*)infile->Get("eff_ttbar_msd_btag");
  vector<WorkingPoint> working_points;
  for(auto wp : working_point_mistag_rates) {
    working_points.push_back(init_working_point(wp, graph_eff_qcd_msd, graph_eff_ttbar_msd, graph_eff_qcd_msd_btag, graph_eff_ttbar_msd_btag));
  }
  for(auto wp : working_points) {
    cout << "WP:" << endl;
    cout << "Target QCD eff: " << wp.target_eff_qcd << endl;
    cout << "Actual QCD eff: " << wp.real_eff_qcd << endl;
    cout << "Actual tt eff:  " << wp.real_eff_ttbar << endl;
    cout << "tau32 cut:      " << wp.tau32cut << endl;
    cout << "Actual QCD eff with b-tag: " << wp.real_eff_qcd_with_btag << endl;
    cout << "Actual tt eff with b-tag:  " << wp.real_eff_ttbar_with_btag << endl;
  }

  TCanvas *c = new TCanvas("canvas", "canvas title", 600, 600);
  c->cd();
  float margin_l = 0.15;
  float margin_r = 0.05;
  float margin_b = 0.12;
  float margin_t = 0.08;
  c->SetMargin(margin_l, margin_r, margin_b, margin_t); //lrbt
  auto coord = new CoordinateConverter();
  coord->init(margin_l, margin_r, margin_b, margin_t);
  if(log_y) c->SetLogy();
  c->SetTickx(1);
  c->SetTicky(1);

  TMultiGraph *mg = new TMultiGraph("mg", "mg title");

  float legend_x1 = coord->ConvertGraphXToPadX(log_y ? 0.45 : 0.05);
  float legend_y1 = coord->ConvertGraphYToPadY(log_y ? 0.05 : 0.47);
  float legend_x2 = coord->ConvertGraphXToPadX(log_y ? 0.95 : 0.55);
  float legend_y2 = coord->ConvertGraphYToPadY(log_y ? 0.35 : 0.77);
  TLegend *legend = new TLegend(legend_x1, legend_y1, legend_x2, legend_y2); // x1 y1 x2 y2
  legend->SetHeader(("#bf{Scan of #tau_{32} for "+pt_bin.text+"}").c_str());
  legend->SetTextSize(0.025);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  // graph-dependent attributes, order: raw, msd, msd_btag
  vector<int> linewidths = {2, 2, 2};
  vector<int> linestyles = {3, 1, 2};
  vector<int> linecolors = {kBlue, kRed, kMagenta};
  vector<string> legends = {"No #it{m}_{SD} window nor b tagging", "105 < #it{m}_{SD} (GeV) < 210", "#it{m}_{SD} as above + subjet b tag"};
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
  // if(log_y) mg->GetHistogram()->SetMinimum(0.0003);
  if(is_eff_qcd || is_eff_ttbar) mg->GetHistogram()->GetXaxis()->SetTitle("#tau_{3}/#tau_{2} upper limit");
  else if(is_roc) mg->GetHistogram()->GetXaxis()->SetTitle("#varepsilon_{S}");
  if(is_eff_qcd || is_roc) mg->GetHistogram()->GetYaxis()->SetTitle("#varepsilon_{B}");
  if(is_eff_ttbar) mg->GetHistogram()->GetYaxis()->SetTitle("#varepsilon_{S}");
  mg->GetHistogram()->GetXaxis()->SetTitleOffset(1.2);
  mg->GetHistogram()->GetYaxis()->SetTitleOffset(1.2);
  // mg->GetHistogram()->GetXaxis()->SetNdivisions(405);
  if(!log_y) mg->GetHistogram()->GetYaxis()->SetNdivisions(405);

  legend->Draw();

  TText *text_top_left = new TText(margin_l, 1-(margin_t-0.01), "AK8 PUPPI(v14)");
  text_top_left->SetTextAlign(11); // left bottom aligned
  text_top_left->SetTextFont(42);
  text_top_left->SetTextSize(0.035);
  text_top_left->SetNDC();
  text_top_left->Draw();

  TText *text_top_right = new TText(1-margin_r, 1-(margin_t-0.01), "UL17 (106X)");
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

  int iterator_wp = -1;
  for(auto wp : working_points) {
    float x=0, y=0;
    if(is_roc) { x=wp.real_eff_ttbar; y=wp.real_eff_qcd; }
    else if(is_eff_qcd) { x=wp.tau32cut; y=wp.real_eff_qcd; }
    else if(is_eff_ttbar) { x=wp.tau32cut; y=wp.real_eff_ttbar; }
    TMarker *marker = new TMarker(x, y, 24);
    marker->SetMarkerColor(kRed);
    marker->SetMarkerSize(1.2);
    marker->Draw();
  }

  // TLatex *info = new TLatex(coord->ConvertGraphXToPadX(0.5), coord->ConvertGraphYToPadY(0.28), pt_bin.text.c_str());
  // info->SetTextAlign(11); // left bottom
  // info->SetTextFont(42);
  // info->SetTextSize(0.035);
  // info->SetNDC();
  // info->Draw();

  string plotName = "plot_"+pt_bin.name+"__"+graph_base_name+"__";
  for(bool do_x : {do_raw, do_msd, do_msd_btag}) { do_x ? plotName += "1" : plotName += "0"; }
  plotName += (log_y ? string("_log") : string("_lin"))+".pdf";
  string plotPath = infileBasePath+plotName;
  c->SaveAs(plotPath.c_str());
}

void plots() {
  vector<PtBin> pt_bins;
  pt_bins.push_back(PtBin{"Pt300to400", "300 < p_{T}^{jet} (GeV) < 400"});
  pt_bins.push_back(PtBin{"Pt400to480", "400 < p_{T}^{jet} (GeV) < 480"});
  pt_bins.push_back(PtBin{"Pt480to600", "480 < p_{T}^{jet} (GeV) < 600"});
  pt_bins.push_back(PtBin{"Pt600toInf", "p_{T}^{jet} (GeV) > 600"});
  pt_bins.push_back(PtBin{"Pt300toInf", "p_{T}^{jet} (GeV) > 300"});

  vector<string> graph_base_names = {"roc"};
  // vector<string> graph_base_names = {"eff_qcd", "eff_ttbar", "roc"};

  for(auto & pt_bin : pt_bins) {
    cout << "Working on " << pt_bin.name << endl;
    for(auto & graph_base_name : graph_base_names) {
      do_plot(pt_bin, graph_base_name, true);
      // do_plot(pt_bin, graph_base_name, false);
    }
  }
}
