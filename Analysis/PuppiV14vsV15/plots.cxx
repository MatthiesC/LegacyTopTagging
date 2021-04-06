typedef struct {
  string hist_name;
  string x_axis_label;
  double x_low;
  double x_high;
  string unit="";
} Variable;


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


string get_bin_width_string(const double bin_width, const Variable & var) {
  string str = to_string(bin_width);
  str.erase(str.find_last_not_of('0')+1, string::npos);
  if(str.back() == *".") str.pop_back(); // https://stackoverflow.com/questions/60956326/comparison-between-pointer-and-integer-in-string-class-c
  if(!var.unit.empty()) str += " "+var.unit;
  return str;
}


void do_plot(const Variable & var, const bool log_y=false) {
  string infileBasePath = (string)getenv("CMSSW_BASE")+"/src/UHH2/LegacyTopTagging/output/PuppiV14vsV15/UL17/";
  string infileName_v14 = "uhh2.AnalysisModuleRunner.MC.TTbarToHadronicPuppiV14_UL17.root";
  string infilePath_v14 = infileBasePath+infileName_v14;
  TFile *infile_v14 = TFile::Open(infilePath_v14.c_str(), "READ");
  string infileName_v15 = "uhh2.AnalysisModuleRunner.MC.TTbarToHadronicPuppiV15_UL17.root";
  string infilePath_v15 = infileBasePath+infileName_v15;
  TFile *infile_v15 = TFile::Open(infilePath_v15.c_str(), "READ");

  TH1F *hist_v14 = (TH1F*)infile_v14->Get(("AK8Hists_matched/"+var.hist_name).c_str());
  TH1F *hist_v15 = (TH1F*)infile_v15->Get(("AK8Hists_matched/"+var.hist_name).c_str());

  vector<TH1F*> hists = { hist_v14, hist_v15 };

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

  double max_maximum = 0;
  for(auto hist : hists) {
    hist->SetTitle("");
    hist->Rebin(5);
    hist->GetXaxis()->SetTitle((var.x_axis_label).c_str());
    hist->GetXaxis()->SetTitleOffset(1.2);
    double bin_width = hist->GetXaxis()->GetBinWidth(1);
    string bin_width_string = get_bin_width_string(bin_width, var);
    hist->GetYaxis()->SetTitle(((string)"Fraction of jets / "+bin_width_string).c_str());
    hist->GetYaxis()->SetTitleOffset(1.9);
    hist->SetMarkerStyle(0);
    hist->GetXaxis()->SetRangeUser(var.x_low, var.x_high);
    hist->Scale(1./hist->Integral());
    if(!log_y) hist->SetMinimum(0.);
    if(hist->GetMaximum() > max_maximum) max_maximum = hist->GetMaximum();
  }
  gStyle->SetLineWidth(2);
  gStyle->SetOptStat(0);
  hist_v14->SetFillColorAlpha(kRed, 0.6);
  hist_v15->SetFillColorAlpha(kBlue, 0.6);
  if(log_y) {
    hist_v14->SetMaximum(max_maximum*10);
    hist_v15->SetMaximum(max_maximum*10);
  }
  else {
    hist_v14->SetMaximum(max_maximum*1.3);
    hist_v15->SetMaximum(max_maximum*1.3);
  }
  hist_v14->Draw("e3");
  // hist_v14->Scale(1./hist_v14->Integral());
  hist_v15->Draw("e3 same");
  // normalize histograms!!!

  TText *text_top_left = new TText(margin_l, 1-(margin_t-0.01), "AK8 PUPPI");
  text_top_left->SetTextAlign(11); // left bottom aligned
  text_top_left->SetTextFont(42);
  text_top_left->SetTextSize(0.035);
  text_top_left->SetNDC();
  text_top_left->Draw();

  TText *text_top_right = new TText(1-margin_r, 1-(margin_t-0.01), "UL17 (CMSSW 10.6.X)");
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

  // Save to disk
  string plotName = "plot_"+var.hist_name;
  plotName += (log_y ? string("_log") : string("_lin"))+".pdf";
  string plotPath = infileBasePath+plotName;
  c->SaveAs(plotPath.c_str());
  delete c;
}


void plots() {
  vector<Variable> variables;
  variables.push_back(Variable{"ak8jets_pt", "#it{p}_{T}^{jet} [GeV]", 300, 1000, "GeV"});
  variables.push_back(Variable{"ak8jets_mSD", "#it{m}_{SD} [GeV]", 0, 250, "GeV"});
  variables.push_back(Variable{"ak8jets_tau32", "#tau_{3}/#tau_{2}", 0, 1});
  variables.push_back(Variable{"ak8jets_maxDeepCSV", "Max. subjet #it{O}_{DeepCSV}^{prob(b)+prob(bb)}", 0, 1});
  for(auto var : variables) {
    do_plot(var, true);
    do_plot(var, false);
  }
}
