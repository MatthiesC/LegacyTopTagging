typedef struct {
  string name;
  string long_name;
  string lumi_fb;
} Year;

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


void do_plot(const Year & year, const Process & process, const bool log_z) {
  string infileBasePath = (string)getenv("CMSSW_BASE")+"/src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel/"+year.name+"/nominal/";
  TFile *infile = TFile::Open((infileBasePath+"hadded/uhh2.AnalysisModuleRunner.MC."+process.name+"__AllMergeScenarios_"+year.name+".root").c_str(), "READ");
  TH2F *hist_raw = (TH2F*)infile->Get("4_BeforeTwoD_Common/twodselection");

  TCanvas *c = new TCanvas("canvas", "canvas title", 600, 600);
  c->cd();
  float margin_l = 0.15;
  float margin_r = 0.19; // 0.05
  float margin_b = 0.12;
  float margin_t = 0.08;
  float tick_length = 0.015; // fraction of canvas width/height
  c->SetMargin(margin_l, margin_r, margin_b, margin_t); //lrbt
  auto coord = new CoordinateConverter();
  coord->init(margin_l, margin_r, margin_b, margin_t);
  if(log_z) c->SetLogz();
  c->SetTickx(0);
  c->SetTicky(0);

  hist_raw->SetTitle("");
  hist_raw->Rebin2D(5, 5); // initial number of bins in x and y: 1000

  double total_sum_of_events(0.);
  double sum_of_events_in_twod_corner(0.);
  for(unsigned int ix = 0; ix <= hist_raw->GetNbinsX()+1; ix++) {
    for(unsigned int iy = 0; iy <= hist_raw->GetNbinsY()+1; iy++) {
      double binc = hist_raw->GetBinContent(ix, iy);
      total_sum_of_events += binc;
      if(ix <= 10 && iy <= 10) sum_of_events_in_twod_corner += binc;
    }
  }
  double rejected_fraction = sum_of_events_in_twod_corner / total_sum_of_events;

  for(unsigned int ix = 0; ix <= hist_raw->GetNbinsX(); ix++) {
    for(unsigned int iy = 0; iy <= hist_raw->GetNbinsY(); iy++) {
      double scale_factor = hist_raw->GetXaxis()->GetBinWidth(ix) * hist_raw->GetYaxis()->GetBinWidth(iy);
      // cout << to_string(scale_factor) << endl;
      hist_raw->SetBinContent(ix, iy, hist_raw->GetBinContent(ix, iy) / scale_factor);
      if(hist_raw->GetBinContent(ix, iy) < 1) hist_raw->SetBinContent(ix, iy, 1);
      hist_raw->SetBinError(ix, iy, hist_raw->GetBinError(ix, iy) / scale_factor); // only needed if we were to plot the uncertainty of each bin (how to do this in 2D plot???)
    }
  }

  TH2F *hist = new TH2F("hist", "", 80, 0, 2, 80, 0, 200);
  for(unsigned int ix = 1; ix <= 80; ix++) {
    for(unsigned int iy = 1; iy <= 80; iy++) {
      hist->SetBinContent(ix, iy, hist_raw->GetBinContent(ix, iy));
      hist->SetBinError(ix, iy, hist_raw->GetBinError(ix, iy));
    }
  }

  hist->GetXaxis()->SetTitle("#Delta#it{R}(#mu, nearest AK4 jet)");
  hist->GetXaxis()->SetTitleOffset(1.3);
  // hist->GetXaxis()->SetRangeUser(0, 2);
  hist->GetXaxis()->SetNdivisions(805, false);
  hist->GetXaxis()->SetLabelOffset(0.02);

  hist->GetYaxis()->SetTitle("#it{p}_{T}^{rel}(#mu, nearest AK4 jet) [GeV]");
  hist->GetYaxis()->SetTitleOffset(1.8);
  // hist->GetYaxis()->SetRangeUser(0, 200);
  hist->GetYaxis()->SetNdivisions(410, false);
  hist->GetYaxis()->SetLabelOffset(0.02);

  hist->GetZaxis()->SetTitle("d^{2}(Events) / (d(#Delta#it{R}) d#it{p}_{T}^{rel}) [GeV^{#minus1}]");
  hist->GetZaxis()->SetTitleOffset(1.9);
  hist->GetZaxis()->ChangeLabel(1,-1,-1,-1,-1,-1,"#leq 1");
  hist->GetZaxis()->SetLabelOffset(0.02);

  gStyle->SetLineWidth(1);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kColorPrintableOnGrey);
  // TColor::InvertPalette();
  // gStyle->SetPalette(kLightTemperature);

  hist->Draw("colz");

  double ratio_hist_width_to_height = (1. - margin_l - margin_r) / (1. - margin_t - margin_b);
  double five_percent_x = 0.05 / ratio_hist_width_to_height;

  // TText *text_top_left = new TText(margin_l, 1-(margin_t-0.01), ((string)"T&P selection ("+year.name+")").c_str());
  TText *text_top_left = new TText(margin_l, 1-(margin_t-0.01), year.long_name.c_str());
  text_top_left->SetTextAlign(11); // left bottom aligned
  text_top_left->SetTextFont(42);
  text_top_left->SetTextSize(0.035);
  text_top_left->SetNDC();
  text_top_left->Draw();

  TLatex *text_top_right = new TLatex(1-margin_r, 1-(margin_t-0.01), (year.lumi_fb+(string)" fb^{#minus1} (13 TeV)").c_str());
  text_top_right->SetTextAlign(31); // right bottom aligned
  text_top_right->SetTextFont(42);
  text_top_right->SetTextSize(0.035);
  text_top_right->SetNDC();
  text_top_right->Draw();

  TText *cms = new TText(coord->ConvertGraphXToPadX(five_percent_x), coord->ConvertGraphYToPadY(0.95), "CMS");
  cms->SetTextAlign(13); // left top
  cms->SetTextFont(62);
  cms->SetTextSize(0.05);
  cms->SetNDC();
  cms->SetTextColor(kWhite);
  cms->Draw();

  TText *prelim = new TText(coord->ConvertGraphXToPadX(five_percent_x), coord->ConvertGraphYToPadY(0.87), "Simulation");
  prelim->SetTextAlign(13); // left top
  prelim->SetTextFont(52);
  prelim->SetTextSize(0.035);
  prelim->SetNDC();
  prelim->SetTextColor(kWhite);
  prelim->Draw();

  TText *prelim2 = new TText(coord->ConvertGraphXToPadX(five_percent_x), coord->ConvertGraphYToPadY(0.807), "Work in Progress");
  prelim2->SetTextAlign(13); // left top
  prelim2->SetTextFont(52);
  prelim2->SetTextSize(0.035);
  prelim2->SetNDC();
  prelim2->SetTextColor(kWhite);
  prelim2->Draw();

  TLatex *watermark = new TLatex(coord->ConvertGraphXToPadX(1. - five_percent_x), coord->ConvertGraphYToPadY(0.113), (process.watermark + (string)" MC").c_str());
  watermark->SetTextAlign(31); // right bottom
  watermark->SetTextFont(42);
  watermark->SetTextSize(0.035);
  watermark->SetNDC();
  watermark->SetTextColor(kWhite);
  watermark->Draw();

  TLatex *watermark2 = new TLatex(coord->ConvertGraphXToPadX(1. - five_percent_x), coord->ConvertGraphYToPadY(0.05), "passing T&P selection");
  watermark2->SetTextAlign(31); // right bottom
  watermark2->SetTextFont(42);
  watermark2->SetTextSize(0.035);
  watermark2->SetNDC();
  watermark2->SetTextColor(kWhite);
  watermark2->Draw();

  // https://root-forum.cern.ch/t/inconsistent-tick-length/18563/8
  c->Update();
  double tickScaleX = (c->GetUxmax() - c->GetUxmin()) / (c->GetX2() - c->GetX1()) * (c->GetWh() * c->GetAbsHNDC());
  double tickScaleY = (c->GetUymax() - c->GetUymin()) / (c->GetY2() - c->GetY1()) * (c->GetWw() * c->GetAbsWNDC());
  hist->GetXaxis()->SetTickLength(-1. * c->GetWh() * tick_length / tickScaleX);
  hist->GetYaxis()->SetTickLength(-1. * c->GetWw() * tick_length / tickScaleY);
  hist->GetZaxis()->SetTickLength(hist->GetYaxis()->GetTickLength() / (1. - margin_t - margin_b)); // divisor must be height of palette in NDC units

  gPad->Update();
  gPad->RedrawAxis();

  TPaletteAxis *palette = (TPaletteAxis*)hist->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(1. - margin_r + 0.01);
  palette->SetX2NDC(1. - margin_r + 0.04);
  // cannot set border for palette, thus drawing some lines around it:
  // double half_pix_ndc_w = 0.5 / c->GetWw();
  // double half_pix_ndc_h = 0.5 / c->GetWh();
  // TLine *palette_line_b = new TLine();
  // palette_line_b->DrawLineNDC(palette->GetX1NDC()-half_pix_ndc_w, palette->GetY1NDC(), palette->GetX2NDC()+half_pix_ndc_w, palette->GetY1NDC());
  // TLine *palette_line_t = new TLine();
  // palette_line_t->DrawLineNDC(palette->GetX1NDC()-half_pix_ndc_w, palette->GetY2NDC(), palette->GetX2NDC()+half_pix_ndc_w, palette->GetY2NDC());
  // TLine *palette_line_l = new TLine();
  // palette_line_l->DrawLineNDC(palette->GetX1NDC(), palette->GetY1NDC()-half_pix_ndc_h, palette->GetX1NDC(), palette->GetY2NDC()+half_pix_ndc_h);

  // highlight the 2D window that is cut away:
  Int_t corner_color = kGreen;
  TLine *corner_t = new TLine();
  corner_t->SetLineColor(corner_color);
  corner_t->SetLineWidth(2);
  corner_t->DrawLine(0, 25, 0.4, 25);
  TArrow *arrow_t = new TArrow(0.2, 25, 0.2, 25+0.03*200/(1.-margin_t-margin_b), 0.02, "|>");
  arrow_t->SetLineColor(corner_color);
  arrow_t->SetFillColor(corner_color);
  arrow_t->SetLineWidth(2);
  arrow_t->Draw();
  TLine *corner_r = new TLine();
  corner_r->SetLineColor(corner_color);
  corner_r->SetLineWidth(2);
  corner_r->DrawLine(0.4, 0, 0.4, 25);
  TArrow *arrow_r = new TArrow(0.4, 12.5, 0.4+0.03*2/(1.-margin_l-margin_r), 12.5, 0.02, "|>");
  arrow_r->SetLineColor(corner_color);
  arrow_r->SetFillColor(corner_color);
  arrow_r->SetLineWidth(2);
  arrow_r->Draw();

  // Save to disk
  string plotName = "plot_TwoDCut_"+process.name+"_"+year.name;
  plotName += (log_z ? string("_log") : string("_lin"))+".pdf";
  string plotDir = infileBasePath+"plots/";
  gSystem->Exec(("mkdir -p "+plotDir).c_str());
  string plotPath = plotDir+plotName;
  c->SaveAs(plotPath.c_str());
  delete c;
}


void plots_TwoDCut() {

  cout << "Caveat: I assume that the TwoDSelection histograms are given in the form of 1000x1000 bins with x = DeltaR [0 to 5] and y = pTrel [0 to 500]" << endl;
  cout << "NB: I will rebin and reduce the ranges to [0 to 2] and [0 to 200]" << endl;

  vector<Year> years;
  years.push_back(Year{"UL17", "Ultra Legacy 2017", "41.5"});
  // years.push_back(Year{"UL18", "59.8"});
  vector<Process> processes;
  processes.push_back(Process{"TTbar", "t#bar{t}"});
  processes.push_back(Process{"QCD", "QCD multijet"});
  for(const Year & year : years) {
    for(const Process & process : processes) {
      do_plot(year, process, true);
      do_plot(year, process, false);
    }
  }
}
