#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#include <TBox.h>
#include <TLine.h>
#include <TText.h>
#include <TLatex.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <THStack.h>
#include <TMath.h>
#include <TMarker.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TEnv.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TObjString.h>
#include <THashList.h>
#include <iostream>

using namespace std;

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


typedef struct {
  string name;
  int color;
} Process;

const vector<Process> processes = {
{"QCD_Mu", kAzure+10},
{"TTbar_FullyMerged", kPink-3},
{"TTbar_SemiMerged", kPink+4},
// {"TTbar_WMerged", kPink+5},
// {"TTbar_QBMerged", kPink+4},
{"TTbar_BkgOrNotMerged", kPink-7},
// {"TTbar_NotMerged", kPink+5},
// {"TTbar_Background", kPink+4},
{"ST_FullyMerged", kOrange},
{"ST_SemiMerged", kOrange-3},
// {"ST_WMerged", kOrange-3},
// {"ST_QBMerged", kOrange+4},
{"ST_BkgOrNotMerged", kOrange+4},
// {"ST_NotMerged", kOrange-3},
// {"ST_Background", kOrange+4},
{"WJetsToLNu", kSpring-3},
{"DYJetsToLLAndDiboson", kSpring-7},
};

const vector<string> syst_names = {
  "btaggingbc",
  "btaggingudsg",
  "jec",
  "jer",
  "scale",
  "pileup",
  "fsr",
  "tageff3prong",
  "tageff2prong",
  "tageff1prong",
  "topptA",
  "topptB",
};

THStack * create_stack(TFile * rootFile, const string & folderName, const bool divide_by_bin_width) {

  THStack *stack = new THStack("stack", "");

  for(int i = processes.size()-1; i >= 0; i--) {
    TH1F *h = (TH1F*)rootFile->Get((folderName+"/"+processes.at(i).name).c_str());
    if(divide_by_bin_width) {
      for(unsigned int j = 1; j <= h->GetNbinsX(); j++) {
        if(h->GetBinContent(j) == 0) continue;
        h->SetBinContent(j, h->GetBinContent(j) / h->GetBinWidth(j));
        h->SetBinError(j, h->GetBinError(j) / h->GetBinWidth(j));
      }
    }
    h->SetFillColorAlpha(processes.at(i).color, 1.);
    h->SetLineColor(kWhite);
    h->SetLineWidth(i == 0 ? 0 : 0);
    h->SetMarkerSize(0);
    stack->Add(h);
  }

  return stack;
}


THStack * scale_stack(const THStack * stack, const double scale_factor) {

  // Caveat: Does not copy over any members like cosmetics settings of the input stack, only the histograms within the stack (with an applied scaling)

  THStack *new_stack = new THStack(((string)stack->GetName()+"_scaled").c_str(), "");
  for(unsigned int i = 0; i < stack->GetNhists(); i++) {
    TH1F *hist = (TH1F*)stack->GetHists()->At(i);
    TH1F *scaled_hist = new TH1F(*hist);
    scaled_hist->Scale(scale_factor);
    new_stack->Add(scaled_hist);
  }

  return new_stack;
}


TH1F * create_th1f(const int n_events) {

  TH1F *hist = new TH1F("hist", "", 20, 0, 1);
  hist->FillRandom("gaus", n_events);

  return hist;
}


double get_syst_unc_for_bin(const unsigned int ibin, TFile * rootFile, const string & folderName, const string pm) {
  double squarederr(0);
  // vector<Process> processes;
  // processes.push_back({"QCD_Mu", kAzure+10});
  // processes.push_back({"TTbar_FullyMerged", kPink-3});
  // processes.push_back({"TTbar_WMerged", kPink+5});
  // processes.push_back({"TTbar_QBMerged", kPink+4});
  // processes.push_back({"TTbar_BkgOrNotMerged", kPink-7});
  // processes.push_back({"ST_FullyMerged", kOrange});
  // processes.push_back({"ST_WMerged", kOrange-3});
  // processes.push_back({"ST_QBMerged", kOrange+4});
  // processes.push_back({"ST_BkgOrNotMerged", kOrange+3});
  // processes.push_back({"WJetsToLNu", kSpring-3});
  // processes.push_back({"DYJetsToLLAndDiboson", kSpring-7});

  for(Process p : processes) {
    TH1F *hNominal = (TH1F*)rootFile->Get((folderName+"/"+p.name).c_str());
    double nominal = hNominal->GetBinContent(ibin);
    for(auto s : syst_names) {
    // unc_hists.push_back((TH1F*)rootFile->Get((folderName+"/"+s+"Down")));
      TH1F *hSysDown = (TH1F*)rootFile->Get((folderName+"/"+p.name+"_"+s+"Down").c_str());
      TH1F *hSysUp = (TH1F*)rootFile->Get((folderName+"/"+p.name+"_"+s+"Up").c_str());
      const double errDown = hSysDown->GetBinContent(ibin) - nominal;
      const double errUp = hSysUp->GetBinContent(ibin) - nominal;
      const bool errDown_is_positive = errDown >= 0;
      const bool errUp_is_positive = errUp >= 0;
      if(pm == "plus") {
        if(errDown_is_positive && !errUp_is_positive) squarederr += errDown*errDown;
        else if(!errDown_is_positive && errUp_is_positive) squarederr += errUp*errUp;
        else if(errDown_is_positive && errUp_is_positive) {
          if(errDown < errUp) squarederr += errUp*errUp;
          else squarederr += errDown*errDown;
        }
        // else: nothing added to squarederr
      }
      else if(pm == "minus") {
        if(errDown_is_positive && !errUp_is_positive) squarederr += errUp*errUp;
        else if(!errDown_is_positive && errUp_is_positive) squarederr += errDown*errDown;
        else if(!errDown_is_positive && !errUp_is_positive) {
          if(errDown > errUp) squarederr += errUp*errUp;
          else squarederr += errDown*errDown;
        }
        // else: nothing added to squarederr
      }
    }
  }
  return TMath::Sqrt(squarederr);
}


TGraphAsymmErrors * get_stack_unc(const TH1F * last, TFile * rootFile, const string & folderName, const bool divide_by_bin_width, const double scale_factor) {
  unsigned int n = last->GetNbinsX();
  double x[n], y[n], exl[n], eyl[n], exh[n], eyh[n];
  for(unsigned int i = 0; i < n; i++) {
    unsigned int ibin = i+1;
    x[i] = last->GetBinCenter(ibin);
    y[i] = last->GetBinContent(ibin) * scale_factor;
    exl[i] = last->GetBinWidth(ibin) / 2.;
    eyl[i] = scale_factor * TMath::Sqrt(TMath::Power(last->GetBinError(ibin)*(divide_by_bin_width ? last->GetBinWidth(ibin) : 1.), 2) + TMath::Power(get_syst_unc_for_bin(ibin, rootFile, folderName, "minus"), 2)) / (divide_by_bin_width ? last->GetBinWidth(ibin) : 1.);
    exh[i] = last->GetBinWidth(ibin) / 2.;
    eyh[i] = scale_factor * TMath::Sqrt(TMath::Power(last->GetBinError(ibin)*(divide_by_bin_width ? last->GetBinWidth(ibin) : 1.), 2) + TMath::Power(get_syst_unc_for_bin(ibin, rootFile, folderName, "plus"), 2)) / (divide_by_bin_width ? last->GetBinWidth(ibin) : 1.);
  }
  return new TGraphAsymmErrors(n, x, y, exl, exh, eyl, eyh);
}

TGraphAsymmErrors * get_stack_unc_PrePostFitShapes(TFile * rootFile, const string & folderName, const bool divide_by_bin_width, const double scale_factor) {
  const TH1F * totalprocs = (TH1F*)rootFile->Get((folderName+"/TotalProcs").c_str());
  unsigned int n = totalprocs->GetNbinsX();
  double x[n], y[n], exl[n], eyl[n], exh[n], eyh[n];
  for(unsigned int i = 0; i < n; i++) {
    unsigned int ibin = i+1;
    x[i] = totalprocs->GetBinCenter(ibin);
    y[i] = totalprocs->GetBinContent(ibin) * scale_factor / (divide_by_bin_width ? totalprocs->GetBinWidth(ibin) : 1.);
    exl[i] = totalprocs->GetBinWidth(ibin) / 2.;
    eyl[i] = totalprocs->GetBinError(ibin) * scale_factor / (divide_by_bin_width ? totalprocs->GetBinWidth(ibin) : 1.);
    exh[i] = totalprocs->GetBinWidth(ibin) / 2.;
    eyh[i] = totalprocs->GetBinError(ibin) * scale_factor / (divide_by_bin_width ? totalprocs->GetBinWidth(ibin) : 1.);
  }
  return new TGraphAsymmErrors(n, x, y, exl, exh, eyl, eyh);
}


// TGraphAsymmErrors * create_ratio_totalunc(const TGraphAsymmErrors * total_unc) {
//   TGraphAsymmErrors * ratio_totalunc = new TGraphAsymmErrors(*total_unc);
//   ratio_totalunc->SetTitle("DummyTitleToAvoidMemoryLeak");
//   for(unsigned int i = 0; i < ratio_totalunc->GetN(); i++) {
//     double x, y;
//     ratio_totalunc->GetPoint(i, x, y);
//     ratio_totalunc->SetPointEYlow(i, ratio_totalunc->GetErrorYlow(i) / y);
//     ratio_totalunc->SetPointEYhigh(i, ratio_totalunc->GetErrorYhigh(i) / y);
//     ratio_totalunc->SetPoint(i, x, 1.);
//   }
//   return ratio_totalunc;
// }

TGraphAsymmErrors * create_ratio_totalunc(const TGraphAsymmErrors * total_unc) {
  TGraphAsymmErrors * ratio_totalunc = new TGraphAsymmErrors(total_unc->GetN());
  for(unsigned int i = 0; i < total_unc->GetN(); i++) {
    double x, y;
    total_unc->GetPoint(i, x, y);
    ratio_totalunc->SetPoint(i, x, 1.);
    if(y <= 0) ratio_totalunc->SetPointError(i, total_unc->GetErrorXlow(i), total_unc->GetErrorXhigh(i), 0, 0); // like SFramePlotter's IgnoreEmptyBins
    else ratio_totalunc->SetPointError(i, total_unc->GetErrorXlow(i), total_unc->GetErrorXhigh(i), total_unc->GetErrorYlow(i) / y, total_unc->GetErrorYhigh(i) / y);
  }
  return ratio_totalunc;
}


TH1F * create_ratio_mc_stat(const TH1F * last) {
  TH1F *ratio_mc_stat = new TH1F(*last);
  for(unsigned int i = 0; i <= ratio_mc_stat->GetNbinsX()+1; i++) {
    ratio_mc_stat->SetBinContent(i, 1.);
    if(last->GetBinContent(i) <= 0) ratio_mc_stat->SetBinError(i, 0); // like SFramePlotter's IgnoreEmptyBins
    else ratio_mc_stat->SetBinError(i, last->GetBinError(i) / last->GetBinContent(i));
  }
  return ratio_mc_stat;
}


TH1F * create_ratio_data(const TH1F * hist1, const TH1F * hist2) { // hist1 = last of stack, hist2 = data
  TH1F *ratio = new TH1F(*hist2);
  // ratio->Reset();
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


enum class PlotType {
  isPrefit,
  isCombPrefit,
  isCombPostfit,
};


const map<string, PlotType> kStringToPlotType = {
  {"Prefit", PlotType::isPrefit},
  {"prefitComb", PlotType::isCombPrefit},
  {"postfitComb", PlotType::isCombPostfit},
};


typedef struct {
  PlotType fPlotType;
  TString fInputFilePath;
  TString fFolderName;
  TString fPlotDir;
  TString fPlotName;
  TString fTextTopLeft = "???";
  TString fTextTopRight = "??? fb^{#minus1} (??? TeV)";
  TString fTextPrelim = "Work in Progress";
  bool fDivideByBinWidth = true;
} PlotterArguments;


typedef struct {
  TString fName;
  TString fXAxisTitle;
  TString fYAxisTitle;
  bool fDivideByBinWidth = true;
} Variable;


// const map<TString, Variable> kVariables = {
//   {"pt", Variable{"pt", "Probe jet #it{p}_{T} [GeV]", "d(Events) / d#it{p}_{T} [GeV^{#minus1}]"}},
//   {"drlepton", Variable{"drlepton", "#Delta#it{R}(probe jet, #mu)", "d(Events) / d(#Delta#it{R})"}},
//   {"eta", Variable{"eta", "Probe jet #eta", "d(Events) / d#eta"}},
//   {"phi", Variable{"phi", "Probe jet #phi [rad]", "d(Events) / d#phi [rad^{#minus1}]"}},
//   {"mass", Variable{"mass", "Probe jet #it{m}_{jet} [GeV]", "d(Events) / d#it{m}_{jet} [GeV^{#minus1}]"}},
//   {"mSD", Variable{"mSD", "Probe jet #it{m}_{SD} [GeV]", "d(Events) / d#it{m}_{SD} [GeV^{#minus1}]"}},
//   {"tau32", Variable{"tau32", "Probe jet #tau_{3}/#tau_{2}", "d(Events) / d(#tau_{3}/#tau_{2})"}},
//   {"maxDeepCSV", Variable{"maxDeepCSV", "Max. #it{O}_{DeepCSV}^{prob(b)+prob(bb)} of probe subjets", "d(Events) / d#it{O}_{DeepCSV}^{prob(b)+prob(bb)}"}},
//   {"mpair", Variable{"mpair", "Min. #it{m}_{ij} [GeV] of leading three probe subjets", "d(Events) / d#it{m}_{ij} [GeV^{#minus1}]"}},
//   {"fpt1", Variable{"fpt1", "#it{p}_{T} fraction of leading probe subjet", "d(Events) / d(#it{p}_{T} fraction)"}},
// };

const map<TString, Variable> kVariables = {
  {"pt", Variable{"pt", "Probe jet #it{p}_{T} [GeV]", "Events / GeV"}},
  {"drlepton", Variable{"drlepton", "#Delta#it{R}(probe jet, #mu)", "Events / unit"}},
  {"eta", Variable{"eta", "Probe jet #eta", "Events / unit"}},
  {"phi", Variable{"phi", "Probe jet #phi [rad]", "Events / rad"}},
  {"mass", Variable{"mass", "Probe jet #it{m}_{jet} [GeV]", "Events / GeV"}},
  {"mSD", Variable{"mSD", "Probe jet #it{m}_{SD} [GeV]", "Events / GeV"}},
  {"tau32", Variable{"tau32", "Probe jet #tau_{3}/#tau_{2}", "Events / unit"}},
  {"maxDeepCSV", Variable{"maxDeepCSV", "Max. #it{O}_{DeepCSV}^{prob(b)+prob(bb)} of probe subjets", "Events / unit"}},
  {"mpair", Variable{"mpair", "Min. #it{m}_{ij} [GeV] of leading three probe subjets", "Events / GeV"}},
  {"fpt1", Variable{"fpt1", "#it{p}_{T} fraction of leading probe subjet", "Events / unit"}},
  {"nsub", Variable{"nsub", "Number of probe subjets", "Events", false}},
};


// void plotter(TFile * rootFile, const TString & folderName, const string & workDir, const bool divide_by_bin_width=true) {
void plotter(const PlotterArguments & args) {

  TFile * rootFile = TFile::Open(args.fInputFilePath.Data(), "READ");
  const TString & folderName = args.fFolderName;
  TString variable = ((TObjString*)folderName.Tokenize("_")->Last())->GetString();
  // const bool divide_by_bin_width = args.fDivideByBinWidth;
  const bool divide_by_bin_width = kVariables.at(variable).fDivideByBinWidth;


  // bool log_y(false);
  // if(folderName.Contains("mSD") && folderName.Contains("Fail")) log_y = true;

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
  // if(log_y) p_main->SetLogy();
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
  THStack *stack = create_stack(rootFile, (string)folderName, divide_by_bin_width);
  // stack = scale_stack(stack, 1000);

  TH1F *data = (TH1F*)rootFile->Get(((string)folderName+"/data_obs").c_str());
  if(divide_by_bin_width) {
    for(unsigned int j = 1; j <= data->GetNbinsX(); j++) {
      if(data->GetBinContent(j) == 0) continue;
      data->SetBinContent(j, data->GetBinContent(j) / data->GetBinWidth(j));
      data->SetBinError(j, data->GetBinError(j) / data->GetBinWidth(j));
    }
  }
  data->SetLineWidth(1);
  data->SetLineColor(kBlack);
  data->SetMarkerStyle(8);

  const double maximum_stack_saved = ((TH1F*)stack->GetStack()->Last())->GetBinContent(((TH1F*)stack->GetStack()->Last())->GetMaximumBin());
  const double maximum_data_saved = data->GetBinContent(data->GetMaximumBin());
  const double maximum = max(maximum_stack_saved, maximum_data_saved) * 1.5;

  // Re-scale the stack such that the y axis labels do not have higher values than O(100)
  // Changes y axis title: Events -> 10^N Events (N = multiple of 3)
  // THIS PROBABLY ONLY WORKS AS INTENDED WHEN USING LINEAR SCALE
  // ONLY WORKS FOR POSITIVE EXPONENTS (BIG NUMBERS)
  // // stack->GetHistogram()->GetYaxis()->SetNoExponent(); // disable default 10^N at the top of y axis
  THStack *original_stack = new THStack(*stack);
  int exponent(0);
  const int exponent_increment(3);
  const double exp_result(TMath::Power(10, exponent_increment));
  const double scale_factor = 1./exp_result;
  double total_scale_factor(1.);
  double new_maximum = maximum;
  while(new_maximum >= exp_result*10) {
    exponent += exponent_increment;
    total_scale_factor *= scale_factor;
    new_maximum = new_maximum * scale_factor;
  }
  stack = scale_stack(stack, total_scale_factor);
  data->Scale(total_scale_factor);

  TGraphAsymmErrors *stack_totalunc = nullptr;
  if(args.fPlotType == PlotType::isPrefit) {
    stack_totalunc = get_stack_unc((TH1F*)original_stack->GetStack()->Last(), rootFile, (string)folderName, divide_by_bin_width, total_scale_factor);
  }
  else if (args.fPlotType == PlotType::isCombPrefit || args.fPlotType == PlotType::isCombPostfit) {
    stack_totalunc = get_stack_unc_PrePostFitShapes(rootFile, (string)folderName, divide_by_bin_width, total_scale_factor);
  }
  else {
    cout << "Invalid PlotType" << endl;
  }
  stack_totalunc->SetLineWidth(0);
  stack_totalunc->SetFillColor(kGray+1);
  stack_totalunc->SetFillStyle(3254); //3354

  stack->Draw("hist");
  stack_totalunc->Draw("2");
  data->Draw("same e x0");

  // Segmentation violation if stack is not drawn before calling the following cosmetics settings
  TH1F *last = (TH1F*)stack->GetStack()->Last();
  last->SetMaximum(new_maximum);
  if(divide_by_bin_width) stack->GetHistogram()->GetYaxis()->SetTitle(kVariables.at(variable).fYAxisTitle.Data());
  else stack->GetHistogram()->GetYaxis()->SetTitle("Events / bin");
  if(exponent > 0) {
    // stack->GetHistogram()->GetYaxis()->SetTitle(((TString)("10^{")+exponent+"} "+stack->GetHistogram()->GetYaxis()->GetTitle()).Data()); // Events -> 10^N Events
    TString old_title = stack->GetHistogram()->GetYaxis()->GetTitle();
    TString unit_name = ((TObjString*)old_title.Tokenize(" ")->Last())->GetString();
    TString new_title;
    if(divide_by_bin_width) new_title = old_title.ReplaceAll(unit_name, (TString)"10^{#minus"+exponent+"} "+unit_name).ReplaceAll("unit", "units");
    else new_title = (TString)"10^{"+exponent+"} #times "+old_title;
    stack->GetHistogram()->GetYaxis()->SetTitle(new_title.Data());
  }
  stack->GetHistogram()->GetYaxis()->SetLabelSize(stack->GetHistogram()->GetYaxis()->GetLabelSize()/(1.-border_y));
  stack->GetHistogram()->GetYaxis()->SetTitleSize(stack->GetHistogram()->GetYaxis()->GetTitleSize()/(1.-border_y));
  stack->GetHistogram()->GetXaxis()->SetTitle(kVariables.at(variable).fXAxisTitle.Data());
  stack->GetHistogram()->GetXaxis()->SetLabelOffset(5); // hack to let the x axis labels vanish under ratio pad
  if(kVariables.at(variable).fXAxisTitle.BeginsWith("Number")) {
    stack->GetHistogram()->GetXaxis()->SetNdivisions(stack->GetHistogram()->GetNbinsX());
  }


  // https://root-forum.cern.ch/t/inconsistent-tick-length/18563/8
  p_main->Update();
  double tickScaleX = (p_main->GetUxmax() - p_main->GetUxmin()) / (p_main->GetX2() - p_main->GetX1()) * (p_main->GetWh() * p_main->GetAbsHNDC());
  double tickScaleY = (p_main->GetUymax() - p_main->GetUymin()) / (p_main->GetY2() - p_main->GetY1()) * (p_main->GetWw() * p_main->GetAbsWNDC());
  stack->GetYaxis()->SetTickLength(c->GetWw() * tick_length / tickScaleY);
  stack->GetXaxis()->SetTickLength(c->GetWh() * tick_length / tickScaleX);
  redrawBorder();


  p_ratio->cd();
  TH1F *null_hist = new TH1F(*data); null_hist->Reset();
  null_hist->Draw(); // option for data points, need to change this later

  null_hist->GetXaxis()->SetTitle(null_hist->GetTitle());
  null_hist->SetTitle("");

  null_hist->GetXaxis()->SetLabelSize(null_hist->GetXaxis()->GetLabelSize()/border_y);
  null_hist->GetXaxis()->SetTitleSize(null_hist->GetXaxis()->GetTitleSize()/border_y);
  null_hist->GetXaxis()->SetTitle(stack->GetHistogram()->GetXaxis()->GetTitle());
  null_hist->GetXaxis()->SetTitleOffset(1.3);
  null_hist->GetXaxis()->SetNdivisions(stack->GetHistogram()->GetXaxis()->GetNdivisions());

  null_hist->GetYaxis()->SetTitle("Data / pred.");
  null_hist->GetYaxis()->SetLabelSize(null_hist->GetYaxis()->GetLabelSize()/border_y);
  null_hist->GetYaxis()->SetTitleSize(null_hist->GetYaxis()->GetTitleSize()/border_y);
  null_hist->GetYaxis()->SetTitleOffset(0.45);
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

  TH1F *ratio_mc_stat = nullptr;
  if(args.fPlotType == PlotType::isPrefit) {
    ratio_mc_stat = create_ratio_mc_stat(last);
    ratio_mc_stat->Draw("same e2");
    ratio_mc_stat->SetFillColor(kGray);
    ratio_mc_stat->SetFillStyle(1001);
    ratio_mc_stat->SetTitle("MC stat. unc.");
    ratio_mc_stat->SetLineWidth(0);
    ratio_mc_stat->SetMarkerStyle(0);
  }

  TGraphAsymmErrors *ratio_totalunc = create_ratio_totalunc(stack_totalunc);
  // ratio_totalunc->SetTitle("DummyTitleToAvoidMemoryLeak");
  ratio_totalunc->Draw("same 2");
  ratio_totalunc->SetFillColor(kGray+1);
  ratio_totalunc->SetFillStyle(3254);
  stack_totalunc->SetLineWidth(0);
  stack_totalunc->SetMarkerStyle(0);
  // ratio_totalunc->SetFillStyle(1001); // if using this fill style, the plot bugs out and many of the following texts are not drawn. Only God knows why.

  TH1F *ratio_data = create_ratio_data(last, data);
  ratio_data->Draw("same e x0");

  TLine *ratio_line = new TLine(null_hist->GetXaxis()->GetXmin(), 1., null_hist->GetXaxis()->GetXmax(), 1.);
  ratio_line->SetLineStyle(2);
  ratio_line->Draw();

  // Minimize the "outliers" of Data in ratio plot by adjusting the y axis range:
  int data_ratio_outliers(0);
  for(unsigned int i = 1; i <= ratio_data->GetNbinsX(); i++) {
    const double binc = ratio_data->GetBinContent(i);
    if(binc <= 0) continue;
    if(binc < null_hist->GetMinimum() || binc > null_hist->GetMaximum()) data_ratio_outliers++;
  }
  if(data_ratio_outliers > 3) {
    null_hist->SetMinimum(0.3);
    null_hist->SetMaximum(1.7);
    null_hist->GetYaxis()->SetNdivisions(503);
  }
  data_ratio_outliers = 0;
  for(unsigned int i = 1; i <= ratio_data->GetNbinsX(); i++) {
    const double binc = ratio_data->GetBinContent(i);
    if(binc <= 0) continue;
    if(binc < null_hist->GetMinimum() || binc > null_hist->GetMaximum()) data_ratio_outliers++;
  }
  if(data_ratio_outliers > 3) {
    null_hist->SetMinimum(0.0);
    null_hist->SetMaximum(2.0);
    null_hist->GetYaxis()->SetNdivisions(503);
  }

  // const double height_of_ratio_plot_in_pixels = p_ratio->GetWh()*p_ratio->GetAbsHNDC()*(1 - p_ratio->GetTopMargin() - p_ratio->GetBottomMargin());
  // const double ratio_plot_y_height_corresponding_to_half_marker_size = 8*(null_hist->GetMaximum() - null_hist->GetMinimum()) * (0.5*ratio_data->GetMarkerSize() / height_of_ratio_plot_in_pixels);
  // A marker size of 1 corresponds to 8 pixels (https://root.cern.ch/doc/master/classTAttMarker.html)
  for(unsigned int i = 1; i <= ratio_data->GetNbinsX(); i++) {
    const double binc = ratio_data->GetBinContent(i);
    if(binc <= 0) continue;
    if(binc > null_hist->GetMaximum()) {
      TMarker *marker = new TMarker();
      marker->SetMarkerStyle(22);
      marker->SetMarkerSize(ratio_data->GetMarkerSize());
      // marker->DrawMarker(ratio_data->GetBinCenter(i), null_hist->GetMaximum() - ratio_plot_y_height_corresponding_to_half_marker_size); //  - 4/p_ratio->GetWh()
      marker->DrawMarker(ratio_data->GetBinCenter(i), null_hist->GetMaximum());
    }
    else if(binc < null_hist->GetMinimum()) {
      TMarker *marker = new TMarker();
      marker->SetMarkerStyle(23);
      marker->SetMarkerSize(ratio_data->GetMarkerSize());
      // marker->DrawMarker(ratio_data->GetBinCenter(i), null_hist->GetMinimum() + ratio_plot_y_height_corresponding_to_half_marker_size);
      marker->DrawMarker(ratio_data->GetBinCenter(i), null_hist->GetMinimum());
    }
  }

  // TLegend *leg_mcstat = new TLegend();
  // leg_mcstat->SetTextSize(0.02);
  // leg_mcstat->SetBorderSize(0);
  // leg_mcstat->SetFillStyle(0);
  // leg_mcstat->AddEntry(ratio_mc_stat);

  redrawBorder();


  c->cd();
  auto coord = new CoordinateConverter();
  coord->init(margin_l, margin_r, margin_b, margin_t);

  // TLatex *text_top_left = new TLatex(margin_l, 1-(margin_t-0.01), "#mu+jets, 1b0t1W region");
  TLatex *text_top_left = new TLatex(margin_l, 1-(margin_t-0.01), args.fTextTopLeft.Data()); // todo: hotvr puppi
  text_top_left->SetTextAlign(11); // left bottom aligned
  text_top_left->SetTextFont(42);
  text_top_left->SetTextSize(0.035);
  text_top_left->SetNDC();
  text_top_left->Draw();

  // string string_text_top_right = "41.5 fb^{#minus1} (13 TeV)"; // todo: other years have other lumis!
  TLatex *text_top_right = new TLatex(1-margin_r, 1-(margin_t-0.01), args.fTextTopRight.Data());
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

  TText *prelim = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.87), args.fTextPrelim.Data()); //Work in Progress //Preliminary
  prelim->SetTextAlign(13); // left top
  prelim->SetTextFont(52);
  prelim->SetTextSize(0.035);
  prelim->SetNDC();
  prelim->Draw();

  if(args.fPlotType == PlotType::isPrefit) {
    TText *text_mcstat = new TText(1-(margin_r-0.02), border_y*0.75, ratio_mc_stat->GetTitle());
    text_mcstat->SetTextAlign(22); // center center
    text_mcstat->SetTextFont(42);
    text_mcstat->SetTextSize(0.018);
    text_mcstat->SetTextAngle(90);
    text_mcstat->SetNDC();
    text_mcstat->Draw();

    TBox *box_mcstat = new TBox(1-(margin_r-0.012), border_y*0.45, 1-(margin_r-0.030), border_y*0.53);
    box_mcstat->SetFillColor(ratio_mc_stat->GetFillColor());
    box_mcstat->SetLineWidth(ratio_mc_stat->GetLineWidth());
    // box_mcstat->SetNDC();
    box_mcstat->Draw();
  }


  // float legend_x1 = coord->ConvertGraphXToPadX(0.43);
  // float legend_y1 = coord->ConvertGraphYToPadY(0.7);
  // float legend_x2 = coord->ConvertGraphXToPadX(0.65); //71
  // float legend_y2 = coord->ConvertGraphYToPadY(0.96);
  float legend_x1 = coord->ConvertGraphXToPadX(0.74);
  float legend_y1 = coord->ConvertGraphYToPadY(0.71);
  float legend_x2 = coord->ConvertGraphXToPadX(0.96); //71
  float legend_y2 = coord->ConvertGraphYToPadY(0.97);
  TLegend *legend = new TLegend(legend_x1, legend_y1, legend_x2, legend_y2); // x1 y1 x2 y2
  legend->SetTextSize(0.025);
  legend->SetBorderSize(0);
  // legend->SetFillColor(kGreen);
  legend->SetFillStyle(0);

  TH1F *leg_wjets = new TH1F();
  leg_wjets->SetFillColor(kSpring-3); leg_wjets->SetLineWidth(0); leg_wjets->SetLineColor(kWhite);
  TH1F *leg_dyjetsvv = new TH1F();
  leg_dyjetsvv->SetFillColor(kSpring-7); leg_dyjetsvv->SetLineWidth(0); leg_dyjetsvv->SetLineColor(kWhite);
  // TH1F *leg_vv = new TH1F();
  // leg_vv->SetFillColor(kSpring-7); leg_vv->SetLineWidth(0); leg_vv->SetLineColor(kWhite);
  TH1F *leg_qcd = new TH1F();
  leg_qcd->SetFillColor(kAzure+10); leg_qcd->SetLineWidth(0); leg_qcd->SetLineColor(kWhite);

  legend->AddEntry(data, "Data", "ep");
  legend->AddEntry(leg_wjets, "W + jets", "f");
  legend->AddEntry(leg_dyjetsvv, "Z/#gamma* + jets, VV", "f");
  // legend->AddEntry(leg_vv, "WW, WZ, ZZ", "f");
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
  leg_merged_ttbar->SetFillColor(kPink); leg_merged_ttbar->SetLineWidth(0); leg_merged_ttbar->SetLineColor(kWhite);
  // TH1F *leg_wmerged_ttbar = new TH1F();
  // leg_wmerged_ttbar->SetFillColor(kPink+5); leg_wmerged_ttbar->SetLineWidth(0); leg_wmerged_ttbar->SetLineColor(kWhite);
  TH1F *leg_qbmerged_ttbar = new TH1F();
  leg_qbmerged_ttbar->SetFillColor(kPink+4); leg_qbmerged_ttbar->SetLineWidth(0); leg_qbmerged_ttbar->SetLineColor(kWhite);
  TH1F *leg_unmerged_ttbar = new TH1F();
  leg_unmerged_ttbar->SetFillColor(kPink-7); leg_unmerged_ttbar->SetLineWidth(0); leg_unmerged_ttbar->SetLineColor(kWhite);

  // legend2->AddEntry(data, "Data", "ep");
  // legend2->AddEntry((TObject*)0, "", "");
  // legend2->AddEntry(leg_merged_ttbar, "", "f");
  // legend2->AddEntry(leg_wmerged_ttbar, "", "f");
  // legend2->AddEntry(leg_qbmerged_ttbar, "", "f");
  // legend2->AddEntry(leg_unmerged_ttbar, "", "f");

  legend2->AddEntry((TObject*)0, "", "");
  legend2->AddEntry(leg_merged_ttbar, "", "f");
  // legend2->AddEntry(leg_wmerged_ttbar, "", "f");
  legend2->AddEntry(leg_qbmerged_ttbar, "", "f");
  legend2->AddEntry(leg_unmerged_ttbar, "", "f");
  legend2->AddEntry((TObject*)0, "", "");

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
  // TH1F *leg_wmerged_singlet = new TH1F();
  // leg_wmerged_singlet->SetFillColor(kOrange-3); leg_wmerged_singlet->SetLineWidth(0); leg_wmerged_singlet->SetLineColor(kWhite);
  TH1F *leg_qbmerged_singlet = new TH1F();
  leg_qbmerged_singlet->SetFillColor(kOrange-3); leg_qbmerged_singlet->SetLineWidth(0); leg_qbmerged_singlet->SetLineColor(kWhite);
  TH1F *leg_unmerged_singlet = new TH1F();
  leg_unmerged_singlet->SetFillColor(kOrange+4); leg_unmerged_singlet->SetLineWidth(0); leg_unmerged_singlet->SetLineColor(kWhite);

  // legend3->AddEntry((TObject*)0, "", "");
  // legend3->AddEntry(leg_merged_singlet, "Fully merged", "f");
  // legend3->AddEntry(leg_wmerged_singlet, "W merged", "f");
  // legend3->AddEntry(leg_qbmerged_singlet, "q + b merged", "f");
  // legend3->AddEntry(leg_unmerged_singlet, "Other / bkg.", "f");

  legend3->AddEntry((TObject*)0, "", "");
  legend3->AddEntry(leg_merged_singlet, "Fully merged", "f");
  // legend3->AddEntry(leg_wmerged_singlet, "W merged", "f");
  legend3->AddEntry(leg_qbmerged_singlet, "Semi-merged", "f");
  legend3->AddEntry(leg_unmerged_singlet, "Not merged", "f");
  legend3->AddEntry((TObject*)0, "", "");

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

  // legend4->AddEntry((TObject*)0, "", "");
  legend4->AddEntry((TObject*)0, "t#bar{t}, single t categories:", "");
  legend4->AddEntry((TObject*)0, "", "");
  legend4->AddEntry((TObject*)0, "", "");
  legend4->AddEntry((TObject*)0, "", "");
  legend4->AddEntry((TObject*)0, "", "");

  legend4->Draw();


  // Save to disk
  gSystem->Exec(((TString)"mkdir -p "+args.fPlotDir).Data());
  TString plotPath = args.fPlotDir+"/"+args.fPlotName;
  c->SaveAs(plotPath.Data());
  delete c;
}


PlotterArguments convert_arguments(int argc, char **argv) {
  unsigned int n(1);
  if(argc != 9) {
    throw invalid_argument("Number of arguments not correct! Please check!");
  }
  PlotterArguments args{
    .fPlotType=kStringToPlotType.at((string)argv[n++]),
    .fInputFilePath=argv[n++],
    .fFolderName=argv[n++],
    .fPlotDir=argv[n++],
    .fPlotName=argv[n++],
    .fTextTopLeft=argv[n++],
    .fTextTopRight=argv[n++],
    .fTextPrelim=argv[n++],
    .fDivideByBinWidth=true,
  };
  return args;
}


int main(int argc, char **argv) {
  plotter(convert_arguments(argc, argv));
  return 0;
}


// void plots() {
//   vector<string> pt_bins = {
//     "Pt300to400",
//     // "Pt400to480",
//     // "Pt480to600",
//     // "Pt600toInf",
//     // "Pt300toInf",
//   };
//   vector<string> jet_versions = {
//     "All",
//     "BTag",
//   };
//   vector<string> wps = {
//     // "BkgEff0p001",
//     // "BkgEff0p005",
//     // "BkgEff0p010",
//     // "BkgEff0p025",
//     "BkgEff0p050",
//     // "NoTau32Cut",
//   };
//   vector<string> regions = {
//     "Pass",
//     "Fail",
//   };
//
//   string outputBaseDir = (string)getenv("CMSSW_BASE")+"/src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel/UL17/combine/AK8/";
//   for(const auto & pt_bin : pt_bins) {
//     for(const auto & jet_version : jet_versions) {
//       for(const auto & wp : wps) {
//         string task_name = pt_bin+"_"+jet_version+"_"+wp;
//         string workDir = outputBaseDir+"/"+task_name;
//         string inputFilePath = workDir+"/"+task_name+"__ProbeJetHists_AK8.root"; // _mSD10
//         TFile *rootFile = TFile::Open(inputFilePath.c_str(), "READ");
//         for(const auto & region : regions) {
//           string folderName = task_name+"_"+region+"_mSD";
//           plotter(rootFile, folderName, workDir);
//         }
//       }
//     }
//   }
//   // string folderName = "Pt480to600_All_BkgEff0p050_Pass_mSD";
//
// }
