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


THStack * create_stack() {

  THStack *stack = new THStack("stack", "");

  const vector<string> h_name = { "ttbar_merged", "ttbar_semimerged", "ttbar_unmerged", "st_merged", "st_semimerged", "st_ummerged", "wjets", "dyjets", "vv", "qcd" };
  // const vector<int> h_color = { kMagenta, kMagenta-3, kMagenta+3, kYellow, kYellow-3, kYellow+3, kCyan, kAzure };
  const vector<int> h_color = { kPink-3, kPink+4, kPink-7, kOrange, kOrange-3, kOrange+4, kSpring-3, kSpring+4, kSpring-7, kAzure+10 };

  for(int i = h_name.size()-1; i >= 0; i--) {
    TH1F *h = new TH1F(h_name.at(i).c_str(), "", 20, 0, 1);
    // int n_process = (int)(100000/(int)(1+TMath::Sqrt(i*10)));
    int n_process = 500;
    h->FillRandom("gaus", n_process);
    h->SetFillColorAlpha(h_color.at(i), 1.);
    h->SetLineColor(kWhite);
    h->SetLineWidth(i == 0 ? 0 : 0);
    stack->Add(h);
  }

  return stack;
}


TH1F * create_th1f(const int n_events) {

  TH1F *hist = new TH1F("hist", "", 20, 0, 1);
  hist->FillRandom("gaus", n_events);

  return hist;
}


TGraphAsymmErrors * get_stack_unc(const TH1F * last) {
  unsigned int n = last->GetNbinsX();
  double x[n], y[n], exl[n], eyl[n], exh[n], eyh[n];
  for(unsigned int i = 0; i < n; i++) {
    x[i] = last->GetBinCenter(i+1);
    y[i] = last->GetBinContent(i+1);
    exl[i] = last->GetBinWidth(i+1) / 2.;
    eyl[i] = last->GetBinError(i+1) * 2; // *2 entfernen!
    exh[i] = last->GetBinWidth(i+1) / 2.;
    eyh[i] = last->GetBinError(i+1) * 2; // *2 entfernen!
  }
  return new TGraphAsymmErrors(n, x, y, exl, exh, eyl, eyh);
}


TGraphAsymmErrors * create_ratio_totalunc(const TGraphAsymmErrors * total_unc) {
  TGraphAsymmErrors * ratio_totalunc = new TGraphAsymmErrors(*total_unc);
  for(unsigned int i = 0; i < ratio_totalunc->GetN(); i++) {
    double x, y;
    ratio_totalunc->GetPoint(i, x, y);
    ratio_totalunc->SetPointEYlow(i, ratio_totalunc->GetErrorYlow(i) / y);
    ratio_totalunc->SetPointEYhigh(i, ratio_totalunc->GetErrorYhigh(i) / y);
    ratio_totalunc->SetPoint(i, x, 1.);
  }
  return ratio_totalunc;
}


TH1F * create_ratio_mc_stat(const TH1F * last) {
  TH1F *ratio_mc_stat = new TH1F(*last);
  for(unsigned int i = 0; i <= ratio_mc_stat->GetNbinsX()+1; i++) {
    ratio_mc_stat->SetBinError(i, ratio_mc_stat->GetBinError(i) / ratio_mc_stat->GetBinContent(i));
    ratio_mc_stat->SetBinContent(i, 1.);
  }
  return ratio_mc_stat;
}


TH1F * create_ratio_data(const TH1F * hist1, const TH1F * hist2) {
  TH1F *ratio = new TH1F(*hist2);
  ratio->Reset();
  for(unsigned int i = 0; i <= ratio->GetNbinsX()+1; i++) {
    if(hist1->GetBinContent(i) <= 0 || hist2->GetBinContent(i) <= 0) continue;
    ratio->SetBinContent(i, hist2->GetBinContent(i) / hist1->GetBinContent(i));
    ratio->SetBinError(i, hist2->GetBinError(i) / hist1->GetBinContent(i));
  }
  return ratio;
}


void redrawBorder() { // https://root-forum.cern.ch/t/how-to-redraw-axis-and-plot-borders/28252
   gPad->Update();
   gPad->RedrawAxis();
   // TLine l;
   // l.DrawLine(gPad->GetUxmin(), gPad->GetUymax(), gPad->GetUxmax(), gPad->GetUymax());
   // l.DrawLine(gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax());
}


void plots() {

  float margin_l = 0.15;
  float margin_r = 0.05;
  float margin_b = 0.12;
  float margin_t = 0.08;

  float tick_length = 0.015; // fraction of canvas width/height

  TCanvas *c = new TCanvas("canvas", "canvas title", 600, 600);
  double border_y = 0.28; // access here, where the canvas is split between main plot and ratio plot
  c->cd();
  TPad *p_main = new TPad("pad_main", "pad title", 0, border_y, 1, 1);
  p_main->SetTopMargin(margin_t/(1.-border_y));
  p_main->SetBottomMargin(0.015/(1.-border_y));
  // p_main->SetMargin();
  // p_main->SetFrameFillColor(kGreen);
  p_main->Draw();
  c->cd();
  TPad *p_ratio = new TPad("pad_ratio", "pad title", 0, 0, 1, border_y);
  p_ratio->SetTopMargin(0.015/border_y);
  p_ratio->SetBottomMargin(margin_b/border_y);
  // p_ratio->SetMargin();
  // p_ratio->SetFrameFillColor(kOrange);
  p_ratio->Draw();

  gStyle->SetLineWidth(1);
  gStyle->SetOptStat(0);

  vector<TPad *> pads = { p_main, p_ratio };
  for(TPad *p : pads) {
    p->SetLeftMargin(margin_l);
    p->SetRightMargin(margin_r);
    p->SetTickx(1);
    p->SetTicky(1);
  }

  p_main->cd();
  THStack *stack = create_stack();
  stack->Draw();

  stack->GetHistogram()->GetYaxis()->SetTitle("Events / bin");
  stack->GetHistogram()->GetYaxis()->SetLabelSize(stack->GetHistogram()->GetYaxis()->GetLabelSize()/(1.-border_y));
  stack->GetHistogram()->GetYaxis()->SetTitleSize(stack->GetHistogram()->GetYaxis()->GetTitleSize()/(1.-border_y));
  // stack->GetHistogram()->GetYaxis()->SetTitleOffset(stack->GetHistogram()->GetYaxis()->GetTitleOffset()/(1.-border_y));
  stack->GetHistogram()->GetXaxis()->SetLabelOffset(5); // hack to let the x axis labels vanish under ratio pad

  // https://root-forum.cern.ch/t/inconsistent-tick-length/18563/8
  p_main->Update();
  double tickScaleX = (p_main->GetUxmax() - p_main->GetUxmin()) / (p_main->GetX2() - p_main->GetX1()) * (p_main->GetWh() * p_main->GetAbsHNDC());
  double tickScaleY = (p_main->GetUymax() - p_main->GetUymin()) / (p_main->GetY2() - p_main->GetY1()) * (p_main->GetWw() * p_main->GetAbsWNDC());
  stack->GetYaxis()->SetTickLength(c->GetWw() * tick_length / tickScaleY);
  stack->GetXaxis()->SetTickLength(c->GetWh() * tick_length / tickScaleX);


  TH1F *last = (TH1F*)stack->GetStack()->Last();

  TGraphAsymmErrors *stack_totalunc = get_stack_unc(last);
  stack_totalunc->Draw("2");
  stack_totalunc->SetLineWidth(0);
  stack_totalunc->SetFillColor(kGray+1);
  stack_totalunc->SetFillStyle(3254); //3354

  TH1F *data = create_th1f(5000);
  data->Draw("same e x0");
  data->SetLineWidth(1);
  data->SetLineColor(kBlack);
  data->SetMarkerStyle(8);

  const double maximum_stack_saved = last->GetBinContent(last->GetMaximumBin());
  const double maximum_data_saved = data->GetBinContent(data->GetMaximumBin());
  const double maximum = max(maximum_stack_saved, maximum_data_saved);
  last->SetMaximum(maximum*1.5);

  redrawBorder();


  p_ratio->cd();
  TH1F *null_hist = new TH1F(*data); null_hist->Reset();
  null_hist->Draw(); // option for data points, need to change this later

  null_hist->GetXaxis()->SetTitle("Probe jet #it{p}_{T} [GeV]");
  null_hist->GetXaxis()->SetLabelSize(null_hist->GetXaxis()->GetLabelSize()/border_y);
  null_hist->GetXaxis()->SetTitleSize(null_hist->GetXaxis()->GetTitleSize()/border_y);
  null_hist->GetXaxis()->SetTitleOffset(1.3);
  null_hist->GetYaxis()->SetTitle("Data / pred.");
  null_hist->GetYaxis()->SetLabelSize(null_hist->GetYaxis()->GetLabelSize()/border_y);
  null_hist->GetYaxis()->SetTitleSize(null_hist->GetYaxis()->GetTitleSize()/border_y);
  null_hist->GetYaxis()->SetTitleOffset(0.5);
  null_hist->GetYaxis()->CenterTitle();
  null_hist->SetMinimum(0.7);
  null_hist->SetMaximum(1.3);
  null_hist->GetYaxis()->SetNdivisions(403);

  // https://root-forum.cern.ch/t/inconsistent-tick-length/18563/8
  p_ratio->Update();
  tickScaleX = (p_ratio->GetUxmax() - p_ratio->GetUxmin()) / (p_ratio->GetX2() - p_ratio->GetX1()) * (p_ratio->GetWh() * p_ratio->GetAbsHNDC());
  tickScaleY = (p_ratio->GetUymax() - p_ratio->GetUymin()) / (p_ratio->GetY2() - p_ratio->GetY1()) * (p_ratio->GetWw() * p_ratio->GetAbsWNDC());
  null_hist->GetYaxis()->SetTickLength(c->GetWw() * tick_length / tickScaleY);
  null_hist->GetXaxis()->SetTickLength(c->GetWh() * tick_length / tickScaleX);

  TGraphAsymmErrors *ratio_totalunc = create_ratio_totalunc(stack_totalunc);
  ratio_totalunc->Draw("same 2");
  ratio_totalunc->SetFillColor(kGray+1);
  // ratio_totalunc->SetFillStyle(3254);
  ratio_totalunc->SetFillStyle(1001);

  TH1F *ratio_mc_stat = create_ratio_mc_stat(last);
  ratio_mc_stat->Draw("same e2");
  ratio_mc_stat->SetFillColor(kGray);
  ratio_mc_stat->SetFillStyle(1001);


  TH1F *ratio_data = create_ratio_data(last, data);
  ratio_data->Draw("same e0 x0");

  TLine *ratio_line = new TLine(null_hist->GetXaxis()->GetXmin(), 1., null_hist->GetXaxis()->GetXmax(), 1.);
  ratio_line->SetLineStyle(2);
  ratio_line->Draw();

  redrawBorder();


  c->cd();
  auto coord = new CoordinateConverter();
  coord->init(margin_l, margin_r, margin_b, margin_t);

  // TLatex *text_top_left = new TLatex(margin_l, 1-(margin_t-0.01), "#mu+jets, 1b0t1W region");
  TLatex *text_top_left = new TLatex(margin_l, 1-(margin_t-0.01), "AK8 PUPPI");
  text_top_left->SetTextAlign(11); // left bottom aligned
  text_top_left->SetTextFont(42);
  text_top_left->SetTextSize(0.035);
  text_top_left->SetNDC();
  text_top_left->Draw();

  string string_text_top_right = "41.5 fb^{#minus1} (13 TeV)";
  TLatex *text_top_right = new TLatex(1-margin_r, 1-(margin_t-0.01), string_text_top_right.c_str());
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

  TText *prelim = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.87), "Preliminary"); //Work in Progress
  prelim->SetTextAlign(13); // left top
  prelim->SetTextFont(52);
  prelim->SetTextSize(0.035);
  prelim->SetNDC();
  prelim->Draw();


  // float legend_x1 = coord->ConvertGraphXToPadX(0.43);
  // float legend_y1 = coord->ConvertGraphYToPadY(0.7);
  // float legend_x2 = coord->ConvertGraphXToPadX(0.65); //71
  // float legend_y2 = coord->ConvertGraphYToPadY(0.96);
  float legend_x1 = coord->ConvertGraphXToPadX(0.75);
  float legend_y1 = coord->ConvertGraphYToPadY(0.71);
  float legend_x2 = coord->ConvertGraphXToPadX(0.97); //71
  float legend_y2 = coord->ConvertGraphYToPadY(0.97);
  TLegend *legend = new TLegend(legend_x1, legend_y1, legend_x2, legend_y2); // x1 y1 x2 y2
  legend->SetTextSize(0.025);
  legend->SetBorderSize(0);
  // legend->SetFillColor(kGreen);
  legend->SetFillStyle(0);

  TH1F *leg_wjets = new TH1F();
  leg_wjets->SetFillColor(kSpring-3); leg_wjets->SetLineWidth(0); leg_wjets->SetLineColor(kWhite);
  TH1F *leg_dyjets = new TH1F();
  leg_dyjets->SetFillColor(kSpring+4); leg_dyjets->SetLineWidth(0); leg_dyjets->SetLineColor(kWhite);
  TH1F *leg_vv = new TH1F();
  leg_vv->SetFillColor(kSpring-7); leg_vv->SetLineWidth(0); leg_vv->SetLineColor(kWhite);
  TH1F *leg_qcd = new TH1F();
  leg_qcd->SetFillColor(kAzure+10); leg_qcd->SetLineWidth(0); leg_qcd->SetLineColor(kWhite);


  legend->AddEntry(leg_wjets, "W + jets", "f");
  legend->AddEntry(leg_dyjets, "Z/#gamma* + jets", "f");
  legend->AddEntry(leg_vv, "WW, WZ, ZZ", "f");
  legend->AddEntry(leg_qcd, "QCD multijet", "f");
  legend->AddEntry(stack_totalunc, "Uncertainty", "f");

  legend->Draw();


  // legend_x1 = coord->ConvertGraphXToPadX(0.67); //71
  // legend_y1 = coord->ConvertGraphYToPadY(0.7);
  // legend_x2 = coord->ConvertGraphXToPadX(0.89); //97
  // legend_y2 = coord->ConvertGraphYToPadY(0.96);
  legend_x1 = coord->ConvertGraphXToPadX(0.46); //71
  legend_y1 = coord->ConvertGraphYToPadY(0.71);
  legend_x2 = coord->ConvertGraphXToPadX(0.68); //97
  legend_y2 = coord->ConvertGraphYToPadY(0.97);
  TLegend *legend2 = new TLegend(legend_x1, legend_y1, legend_x2, legend_y2); // x1 y1 x2 y2
  legend2->SetTextSize(0.025);
  legend2->SetBorderSize(0);
  // legend2->SetFillColor(kYellow);
  legend2->SetFillStyle(0);
  // legend2->SetHeader("t#bar{t} / single t:");

  TH1F *leg_merged_ttbar = new TH1F();
  leg_merged_ttbar->SetFillColor(kPink-3); leg_merged_ttbar->SetLineWidth(0); leg_merged_ttbar->SetLineColor(kWhite);
  TH1F *leg_semimerged_ttbar = new TH1F();
  leg_semimerged_ttbar->SetFillColor(kPink+4); leg_semimerged_ttbar->SetLineWidth(0); leg_semimerged_ttbar->SetLineColor(kWhite);
  TH1F *leg_unmerged_ttbar = new TH1F();
  leg_unmerged_ttbar->SetFillColor(kPink-7); leg_unmerged_ttbar->SetLineWidth(0); leg_unmerged_ttbar->SetLineColor(kWhite);

  legend2->AddEntry(data, "Data", "ep");
  legend2->AddEntry((TObject*)0, "", "");
  legend2->AddEntry(leg_merged_ttbar, "", "f");
  legend2->AddEntry(leg_semimerged_ttbar, "", "f");
  legend2->AddEntry(leg_unmerged_ttbar, "", "f");

  legend2->Draw();


  // legend_x1 = coord->ConvertGraphXToPadX(0.72);
  // legend_y1 = coord->ConvertGraphYToPadY(0.7);
  // legend_x2 = coord->ConvertGraphXToPadX(0.94);
  // legend_y2 = coord->ConvertGraphYToPadY(0.96);
  legend_x1 = coord->ConvertGraphXToPadX(0.505);
  legend_y1 = coord->ConvertGraphYToPadY(0.71);
  legend_x2 = coord->ConvertGraphXToPadX(0.725);
  legend_y2 = coord->ConvertGraphYToPadY(0.97);
  TLegend *legend3 = new TLegend(legend_x1, legend_y1, legend_x2, legend_y2); // x1 y1 x2 y2
  legend3->SetTextSize(0.025);
  legend3->SetBorderSize(0);
  // legend3->SetFillColor(kYellow);
  legend3->SetFillStyle(0);
  // legend3->SetHeader("");

  TH1F *leg_merged_singlet = new TH1F();
  leg_merged_singlet->SetFillColor(kOrange); leg_merged_singlet->SetLineWidth(0); leg_merged_singlet->SetLineColor(kWhite);
  TH1F *leg_semimerged_singlet = new TH1F();
  leg_semimerged_singlet->SetFillColor(kOrange-3); leg_semimerged_singlet->SetLineWidth(0); leg_semimerged_singlet->SetLineColor(kWhite);
  TH1F *leg_unmerged_singlet = new TH1F();
  leg_unmerged_singlet->SetFillColor(kOrange+4); leg_unmerged_singlet->SetLineWidth(0); leg_unmerged_singlet->SetLineColor(kWhite);

  legend3->AddEntry((TObject*)0, "", "");
  legend3->AddEntry((TObject*)0, "", "");
  legend3->AddEntry(leg_merged_singlet, "Merged t", "f");
  legend3->AddEntry(leg_semimerged_singlet, "Semimerged t", "f");
  legend3->AddEntry(leg_unmerged_singlet, "Unmerged t", "f");

  legend3->Draw();


  legend_x1 = coord->ConvertGraphXToPadX(0.417);
  legend_y1 = coord->ConvertGraphYToPadY(0.71);
  legend_x2 = coord->ConvertGraphXToPadX(0.627);
  legend_y2 = coord->ConvertGraphYToPadY(0.97);
  TLegend *legend4 = new TLegend(legend_x1, legend_y1, legend_x2, legend_y2); // x1 y1 x2 y2
  legend4->SetTextSize(0.025);
  legend4->SetBorderSize(0);
  // legend3->SetFillColor(kYellow);
  legend4->SetFillStyle(0);
  // legend3->SetHeader("");

  legend4->AddEntry((TObject*)0, "", "");
  legend4->AddEntry((TObject*)0, "t#bar{t}, single t categories:", "");
  legend4->AddEntry((TObject*)0, "", "");
  legend4->AddEntry((TObject*)0, "", "");
  legend4->AddEntry((TObject*)0, "", "");

  legend4->Draw();





  // Save to disk
  string infileBasePath = (string)getenv("CMSSW_BASE")+"/src/UHH2/LegacyTopTagging/output/test/";
  string plotName = "plot_test.eps";
  // plotName += (log_y ? string("_log") : string("_lin"))+".pdf";
  string plotPath = infileBasePath+plotName;
  c->SaveAs(plotPath.c_str());
  delete c;
}
