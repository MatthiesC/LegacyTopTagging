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
#include <TMultiGraph.h>
#include <TColor.h>
#include <TGaxis.h>
#include <TPaveText.h>
#include <iostream>
#include <sstream>
#include <iomanip>

using namespace std;

template < typename Type > std::string to_str (const Type & t)
{
  std::ostringstream os;
  os << t;
  return os.str ();
}

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


string my_format(float value, int prec=2, int factor=1) {
  std::stringstream stream;
  stream << std::fixed << std::setprecision(prec) << value*factor;
  return stream.str();
}


typedef struct {
  double central = 1.;
  double err_tot_up = 0.;
  double err_tot_down = 0.;
  double err_stat_up = 0.;
  double err_stat_down = 0.;
  double err_syst_up = 0.;
  double err_syst_down = 0.;
} SFInfo;

enum class Process {
  isTTbar,
  isST,
  isOther,
};

enum class ProbeJetAlgo {
  isHOTVR,
  isAK8,
  notValid,
};

typedef struct {
  std::string name = "";
  double mass_min = 0.;
  double mass_max = std::numeric_limits<double>::infinity();
} ProbeJetAlgoInfo;

const std::map<ProbeJetAlgo, ProbeJetAlgoInfo> kProbeJetAlgos = {
  {ProbeJetAlgo::isHOTVR, ProbeJetAlgoInfo{"HOTVR", 140., 220.}},
  {ProbeJetAlgo::isAK8, ProbeJetAlgoInfo{"AK8", 105., 210.}},
};

enum class MergeScenario {
  isFullyMerged,
  isWMerged,
  isQBMerged,
  isNotMerged,
  isBackground,
  isAll,
  isSemiMerged,
  isBkgOrNotMerged,
  isNotFullyOrWMerged,
  isYllufMerged,
};

typedef struct {
  TString name;
  TString print_name = "print name not set";
} MergeScenarioInfo;

const std::map<MergeScenario, MergeScenarioInfo> kMergeScenarios = {
  {MergeScenario::isFullyMerged, MergeScenarioInfo{"FullyMerged", "Fully merged"}},
  {MergeScenario::isWMerged, MergeScenarioInfo{"WMerged", "W merged"}},
  {MergeScenario::isQBMerged, MergeScenarioInfo{"QBMerged", "q+b merged"}},
  {MergeScenario::isNotMerged, MergeScenarioInfo{"NotMerged", "Not merged"}},
  {MergeScenario::isBackground, MergeScenarioInfo{"Background"}},
  {MergeScenario::isAll, MergeScenarioInfo{"AllMergeScenarios"}},
  {MergeScenario::isSemiMerged, MergeScenarioInfo{"SemiMerged", "Semi-merged"}},
  {MergeScenario::isBkgOrNotMerged, MergeScenarioInfo{"BkgOrNotMerged", "Not merged"}},
  {MergeScenario::isNotFullyOrWMerged, MergeScenarioInfo{"NotFullyOrWMerged", "Not merged"}},
  {MergeScenario::isYllufMerged, MergeScenarioInfo{"YllufMerged", "Not merged"}},
};

enum class PtBin {
  Pt200to250,
  Pt250to300,
  Pt300to400,
  Pt400to480,
  Pt480to600,
  Pt600toInf,
  Pt200toInf,
  Pt300toInf,
  Pt400toInf,
};

typedef struct {
  std::string name = "";
  double pt_min = 0.;
  double pt_max = std::numeric_limits<double>::infinity();
} PtBinInfo;

const std::map<PtBin, PtBinInfo> kPtBins = {
  {PtBin::Pt200to250, PtBinInfo{"Pt200to250", 200, 250}},
  {PtBin::Pt250to300, PtBinInfo{"Pt250to300", 250, 300}},
  {PtBin::Pt300to400, PtBinInfo{"Pt300to400", 300, 400}},
  {PtBin::Pt400to480, PtBinInfo{"Pt400to480", 400, 480}},
  {PtBin::Pt480to600, PtBinInfo{"Pt480to600", 480, 600}},
  {PtBin::Pt600toInf, PtBinInfo{"Pt600toInf", 600}},
  {PtBin::Pt200toInf, PtBinInfo{"Pt200toInf", 200}},
  {PtBin::Pt300toInf, PtBinInfo{"Pt300toInf", 300}},
  {PtBin::Pt400toInf, PtBinInfo{"Pt400toInf", 400}},
};

const std::vector<PtBin> kPtBinsHOTVR = {
  PtBin::Pt200to250,
  PtBin::Pt250to300,
  PtBin::Pt300to400,
  PtBin::Pt400to480,
  PtBin::Pt480to600,
  PtBin::Pt600toInf,
  // only used for control distributions etc.:
  // PtBin::Pt200toInf,
  // PtBin::Pt300toInf,
  // PtBin::Pt400toInf,
};

const std::vector<PtBin> kPtBinsAK8 = {
  PtBin::Pt300to400,
  // PtBin::Pt400to480,
  // PtBin::Pt480to600,
  // PtBin::Pt600toInf,
  // only used for control distributions etc.:
  // PtBin::Pt300toInf,
  PtBin::Pt400toInf,
};

enum class WorkingPoint {
  Standard,
  NoTau32Cut,
  BkgEff0p001,
  BkgEff0p003,
  BkgEff0p005,
  BkgEff0p010,
  BkgEff0p025,
  BkgEff0p030,
  BkgEff0p050,
  BkgEff0p100,
};

typedef struct {
  std::string name;
  int color = kGray;
  int color2 = kGray+2;
} WorkingPointInfo;

const std::map<WorkingPoint, WorkingPointInfo> kWorkingPoints = {
  {WorkingPoint::Standard, WorkingPointInfo{"Standard", kAzure+10, kAzure+9}},
  {WorkingPoint::NoTau32Cut, WorkingPointInfo{"NoTau32Cut"}},
  // {WorkingPoint::BkgEff0p001, WorkingPointInfo{"BkgEff0p001", kRed, kRed+2}},
  // {WorkingPoint::BkgEff0p003, WorkingPointInfo{"BkgEff0p003", kMagenta, kMagenta+2}},
  // {WorkingPoint::BkgEff0p005, WorkingPointInfo{"BkgEff0p005", kMagenta, kMagenta+2}},
  // {WorkingPoint::BkgEff0p010, WorkingPointInfo{"BkgEff0p010", kCyan, kCyan+2}},
  // {WorkingPoint::BkgEff0p025, WorkingPointInfo{"BkgEff0p025", kGreen, kGreen+2}},
  // {WorkingPoint::BkgEff0p030, WorkingPointInfo{"BkgEff0p030", kGreen, kGreen+2}},
  // {WorkingPoint::BkgEff0p050, WorkingPointInfo{"BkgEff0p050", kYellow, kYellow+2}},
  // {WorkingPoint::BkgEff0p100, WorkingPointInfo{"BkgEff0p100", kYellow, kYellow+2}},

  // {WorkingPoint::BkgEff0p001, WorkingPointInfo{"BkgEff0p001", kAzure, kAzure-1}},
  // {WorkingPoint::BkgEff0p003, WorkingPointInfo{"BkgEff0p003", kViolet, kViolet-1}},
  // {WorkingPoint::BkgEff0p005, WorkingPointInfo{"BkgEff0p005", kViolet, kViolet-1}},
  // {WorkingPoint::BkgEff0p010, WorkingPointInfo{"BkgEff0p010", kPink, kPink-1}},
  // {WorkingPoint::BkgEff0p025, WorkingPointInfo{"BkgEff0p025", kOrange, kOrange-1}},
  // {WorkingPoint::BkgEff0p030, WorkingPointInfo{"BkgEff0p030", kOrange, kOrange-1}},
  // {WorkingPoint::BkgEff0p050, WorkingPointInfo{"BkgEff0p050", kSpring, kSpring-1}},
  // {WorkingPoint::BkgEff0p100, WorkingPointInfo{"BkgEff0p100", kSpring, kSpring-1}},

  // {WorkingPoint::BkgEff0p001, WorkingPointInfo{"BkgEff0p001", kViolet, kViolet-1}},
  // {WorkingPoint::BkgEff0p003, WorkingPointInfo{"BkgEff0p003", kPink, kPink-1}},
  // {WorkingPoint::BkgEff0p005, WorkingPointInfo{"BkgEff0p005", kPink, kPink-1}},
  // {WorkingPoint::BkgEff0p010, WorkingPointInfo{"BkgEff0p010", kOrange, kOrange-1}},
  // {WorkingPoint::BkgEff0p025, WorkingPointInfo{"BkgEff0p025", kSpring, kSpring-1}},
  // {WorkingPoint::BkgEff0p030, WorkingPointInfo{"BkgEff0p030", kSpring, kSpring-1}},
  // {WorkingPoint::BkgEff0p050, WorkingPointInfo{"BkgEff0p050", kTeal, kTeal-1}},
  // {WorkingPoint::BkgEff0p100, WorkingPointInfo{"BkgEff0p100", kTeal, kTeal-1}},

  {WorkingPoint::BkgEff0p001, WorkingPointInfo{"BkgEff0p001", kPink, kPink-1}},
  {WorkingPoint::BkgEff0p003, WorkingPointInfo{"BkgEff0p003", kOrange, kOrange-1}},
  {WorkingPoint::BkgEff0p005, WorkingPointInfo{"BkgEff0p005", kOrange, kOrange-1}},
  {WorkingPoint::BkgEff0p010, WorkingPointInfo{"BkgEff0p010", kSpring, kSpring-1}},
  {WorkingPoint::BkgEff0p025, WorkingPointInfo{"BkgEff0p025", kTeal, kTeal-1}},
  {WorkingPoint::BkgEff0p030, WorkingPointInfo{"BkgEff0p030", kTeal, kTeal-1}},
  {WorkingPoint::BkgEff0p050, WorkingPointInfo{"BkgEff0p050", kViolet, kViolet-1}},
  {WorkingPoint::BkgEff0p100, WorkingPointInfo{"BkgEff0p100", kViolet, kViolet-1}},
};

typedef std::map<WorkingPoint, double> WorkingPointMap;

const WorkingPointMap kWorkingPointsHOTVR = {
  {WorkingPoint::Standard, 0.56},
  // only used for control distributions etc.:
  // {WorkingPoint::NoTau32Cut, std::numeric_limits<double>::infinity()},
};

// Working points as of https://indico.cern.ch/event/1033067/contributions/4338454/attachments/2238372/3794612/JMAR_update_2021-05-04.pdf
const WorkingPointMap kWorkingPointsAK8 = {
  {WorkingPoint::BkgEff0p001, 0.38},
  {WorkingPoint::BkgEff0p005, 0.47},
  {WorkingPoint::BkgEff0p010, 0.52},
  {WorkingPoint::BkgEff0p025, 0.61},
  {WorkingPoint::BkgEff0p050, 0.69},
  // only used for control distributions etc.:
  // {WorkingPoint::NoTau32Cut, std::numeric_limits<double>::infinity()},
};

const WorkingPointMap kWorkingPointsAK8_old = {
  {WorkingPoint::BkgEff0p001, 0.40},
  {WorkingPoint::BkgEff0p003, 0.46},
  {WorkingPoint::BkgEff0p010, 0.54},
  {WorkingPoint::BkgEff0p030, 0.65},
  {WorkingPoint::BkgEff0p100, 0.80},
  // only used for control distributions etc.:
  // {WorkingPoint::NoTau32Cut, std::numeric_limits<double>::infinity()},
};

enum class JetCategory {
  All,
  Mass,
  BTag,
  MassAndBTag,
  HOTVRCuts,
  HOTVRCutsAndMass,
};

std::map<JetCategory, std::string> kJetCategoryAsString = {
  {JetCategory::All, "All"},
  {JetCategory::Mass, "Mass"},
  {JetCategory::BTag, "BTag"},
  {JetCategory::MassAndBTag, "MassAndBTag"},
  {JetCategory::HOTVRCuts, "HOTVRCuts"},
  {JetCategory::HOTVRCutsAndMass, "HOTVRCutsAndMass"},
};

enum class ExpObs {
  isExp,
  isObs,
};

typedef struct {
  TString short_name;
  TString long_name;
} ExpObsInfo;

const std::map<ExpObs, ExpObsInfo> kExpObs = {
  {ExpObs::isExp, {
    "exp",
    "Expected",
  }},
  {ExpObs::isObs, {
    "obs",
    "Observed",
  }},
};

enum class Year {
  isUL17,
  isUL18,
};

typedef struct {
  TString short_name;
  TString long_name;
  double lumi_fb;
  double lumi_unc; // relative
  int color = kGray;
  int color2 = kGray+2;
} YearInfo;

const std::map<Year, YearInfo> kYears = {
  {Year::isUL17, YearInfo{
    "UL17",
    "Ultra Legacy 2017",
    41.480,
    0.023,
    kSpring,
    kSpring-1,
  }},
  {Year::isUL18, YearInfo{
    "UL18",
    "Ultra Legacy 2018",
    59.830,
    0.025,
    kTeal,
    kTeal-1,
  }},
};

typedef std::map<Year, std::map<MergeScenario, std::map<PtBin, std::map<WorkingPoint, SFInfo>>>> ScaleFactorMap;

class ScaleFactorPlotter {
public:
  ScaleFactorPlotter(const ProbeJetAlgo algo, const vector<Year> & years, const JetCategory jet_cat, const ExpObs expobs, const bool use_PtAll = true);
  ScaleFactorPlotter(const ProbeJetAlgo algo, const Year year,            const JetCategory jet_cat, const ExpObs expobs, const bool use_PtAll = true);
  void ReadScaleFactors();
  void PlotSingleYears();
  void PlotSingleYear(const Year & year);
private:
  void Setup();
  ProbeJetAlgo fAlgo;
  vector<WorkingPoint> fWPs;
  vector<PtBin> fPtBins;
  ExpObs fExpObs;
  vector<Year> fYears;
  JetCategory fJetCat;
  bool fUsePtAll;
  vector<MergeScenario> fMergeScenarios = {MergeScenario::isYllufMerged, MergeScenario::isFullyMerged};
  // vector<MergeScenario> fMergeScenarios = {MergeScenario::isBkgOrNotMerged, MergeScenario::isSemiMerged, MergeScenario::isFullyMerged};
  // vector<MergeScenario> fMergeScenarios = {MergeScenario::isNotFullyOrWMerged, MergeScenario::isWMerged, MergeScenario::isFullyMerged};
  TString fWorkDirBase;
  ScaleFactorMap fSFMap;
};

ScaleFactorPlotter::ScaleFactorPlotter(
  const ProbeJetAlgo algo,
  const vector<Year> & years,
  const JetCategory jet_cat,
  const ExpObs expobs,
  const bool use_PtAll
):
  fAlgo(algo),
  fYears(years),
  fJetCat(jet_cat),
  fExpObs(expobs),
  fUsePtAll(use_PtAll)
{
  Setup();
}

ScaleFactorPlotter::ScaleFactorPlotter(
  const ProbeJetAlgo algo,
  const Year year,
  const JetCategory jet_cat,
  const ExpObs expobs,
  const bool use_PtAll
):
  fAlgo(algo),
  fYears(vector<Year>({year})),
  fJetCat(jet_cat),
  fExpObs(expobs),
  fUsePtAll(use_PtAll)
{
  Setup();
}

void ScaleFactorPlotter::Setup()
{
  fWorkDirBase = (TString)getenv("CMSSW_BASE")+"/src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel/";
  if(fAlgo == ProbeJetAlgo::isAK8) {
    for(const auto & wp : kWorkingPointsAK8) {
      fWPs.push_back(wp.first);
    }
    fPtBins = kPtBinsAK8;
  }
  else if(fAlgo == ProbeJetAlgo::isHOTVR) {
    for(const auto & wp : kWorkingPointsHOTVR) {
      fWPs.push_back(wp.first);
    }
    fPtBins = kPtBinsHOTVR;
  }
}

void ScaleFactorPlotter::ReadScaleFactors() {
  for(const auto & year : fYears) {
    for(const auto & msc : fMergeScenarios) {
      for(const auto & pt_bin : fPtBins) {
        double pt_above_min = kPtBins.at(pt_bin).pt_min+0.001;
        double pt_below_max = kPtBins.at(pt_bin).pt_max-0.001;
        for(const auto & wp : fWPs) {
          TString file_path_pt_part;
          if(fUsePtAll) file_path_pt_part = "PtAll";
          else file_path_pt_part = kPtBins.at(pt_bin).name;
          TString file_path = fWorkDirBase+kYears.at(year).short_name+"/combine/"+kProbeJetAlgos.at(fAlgo).name+"/"
            +kProbeJetAlgos.at(fAlgo).name+"_"+file_path_pt_part+"_"+kJetCategoryAsString.at(fJetCat)+"_"+kWorkingPoints.at(wp).name
            +"/scale_factors_"+kExpObs.at(fExpObs).short_name+".root";
          TFile *file = TFile::Open(file_path.Data(), "READ");
          SFInfo sf;
          TString graph_name = kMergeScenarios.at(msc).name+"_tot";
          TGraphAsymmErrors *graph = (TGraphAsymmErrors*)file->Get(graph_name.Data());
          unsigned int ibin(0);
          double x;
          double y;
          graph->GetPoint(ibin, x, y);
          while(!(x <= pt_below_max && x >= pt_above_min)) {
            ++ibin;
            graph->GetPoint(ibin, x, y);
          }
          sf.central = y;
          sf.err_tot_up = graph->GetErrorYhigh(ibin);
          sf.err_tot_down = graph->GetErrorYlow(ibin);
          graph_name = kMergeScenarios.at(msc).name+"_stat";
          graph = (TGraphAsymmErrors*)file->Get(graph_name.Data());
          sf.err_stat_up = graph->GetErrorYhigh(ibin);
          sf.err_stat_down = graph->GetErrorYlow(ibin);
          graph_name = kMergeScenarios.at(msc).name+"_syst";
          graph = (TGraphAsymmErrors*)file->Get(graph_name.Data());
          sf.err_syst_up = graph->GetErrorYhigh(ibin);
          sf.err_syst_down = graph->GetErrorYlow(ibin);
          fSFMap[year][msc][pt_bin][wp] = sf;
        }
      }
    }
  }
}

void ScaleFactorPlotter::PlotSingleYears() {
  for(const Year year : fYears) {
    PlotSingleYear(year);
  }
}

void ScaleFactorPlotter::PlotSingleYear(const Year & year) {
  TCanvas *c = new TCanvas("canvas", "canvas title", 849, 600); // 900,600 // 849 approx= sqrt(2)*600 (DIN A4 format)
  float margin_l = 0.15*c->GetWh()/c->GetWw();
  float margin_r = 0.05*c->GetWh()/c->GetWw();
  float margin_b = 0.15;
  float margin_t = 0.05;
  double msc_y_gap = 0.75;
  double y_sf_low = 0;
  double y_sf_high = 2;
  double header_y_fraction = 0.2;
  CoordinateConverter *cc = new CoordinateConverter();
  cc->init(margin_l, margin_r, margin_b, margin_t);
  c->cd();
  TPad *p = new TPad("pad", "pad title", 0, 0, 1, 1);
  p->SetTopMargin(margin_t);
  p->SetBottomMargin(margin_b);
  p->SetLeftMargin(margin_l);
  p->SetRightMargin(margin_r);
  p->SetTickx(1);
  p->SetTicky(1);
  p->Draw();
  gStyle->SetLineWidth(1);
  gStyle->SetOptStat(0);
  p->cd();

  vector<double> y_centers;
  map<WorkingPoint, TGraphAsymmErrors*> sf_tot, sf_stat, sf_syst;
  unsigned int points_per_graph = fMergeScenarios.size()*fPtBins.size();
  double sf_y[points_per_graph];
  double sf_x[points_per_graph];
  double sf_x_err[points_per_graph]; for(unsigned int i = 0; i < points_per_graph; i++) sf_x_err[i] = 0.25;
  double sf_err_tot_up[points_per_graph];
  double sf_err_tot_down[points_per_graph];
  double sf_err_stat_up[points_per_graph];
  double sf_err_stat_down[points_per_graph];
  double sf_err_syst_up[points_per_graph];
  double sf_err_syst_down[points_per_graph];
  double x_offset_wp(1.);
  for(const auto & wp : fWPs) {
    unsigned int point_i(0);
    double y_offset(0.);
    for(const auto & msc : fMergeScenarios) {
      double x_offset = x_offset_wp;
      for(const auto & pt_bin : fPtBins) {
        sf_y[point_i] = y_offset+fSFMap.at(year).at(msc).at(pt_bin).at(wp).central;
        sf_x[point_i] = x_offset;
        sf_err_tot_up[point_i] = fSFMap.at(year).at(msc).at(pt_bin).at(wp).err_tot_up;
        sf_err_tot_down[point_i] = fSFMap.at(year).at(msc).at(pt_bin).at(wp).err_tot_down;
        sf_err_stat_up[point_i] = fSFMap.at(year).at(msc).at(pt_bin).at(wp).err_stat_up;
        sf_err_stat_down[point_i] = fSFMap.at(year).at(msc).at(pt_bin).at(wp).err_stat_down;
        sf_err_syst_up[point_i] = fSFMap.at(year).at(msc).at(pt_bin).at(wp).err_syst_up;
        sf_err_syst_down[point_i] = fSFMap.at(year).at(msc).at(pt_bin).at(wp).err_syst_down;
        x_offset += fWPs.size()+1.;
        ++point_i;
      }
      y_centers.push_back(y_offset+1);
      y_offset += (y_sf_high - y_sf_low + msc_y_gap);
    }
    sf_tot[wp] = new TGraphAsymmErrors(points_per_graph, sf_x, sf_y, sf_x_err, sf_x_err, sf_err_tot_down, sf_err_tot_up);
    sf_stat[wp] = new TGraphAsymmErrors(points_per_graph, sf_x, sf_y, sf_x_err, sf_x_err, sf_err_stat_down, sf_err_stat_up);
    sf_syst[wp] = new TGraphAsymmErrors(points_per_graph, sf_x, sf_y, sf_x_err, sf_x_err, sf_err_syst_down, sf_err_syst_up);
    x_offset_wp += 1.;
  }

  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("");
  for(const auto wp : fWPs) {
    auto graph_tot = sf_tot.at(wp);
    graph_tot->SetLineColor(kWorkingPoints.at(wp).color2);
    graph_tot->SetFillColor(kWorkingPoints.at(wp).color2);
    graph_tot->SetMarkerStyle(7);
    mg->Add(graph_tot);
    auto graph_stat = sf_stat.at(wp);
    graph_stat->SetLineColor(kWorkingPoints.at(wp).color);
    graph_stat->SetFillColor(kWorkingPoints.at(wp).color);
    graph_stat->SetMarkerStyle(7);
    mg->Add(graph_stat);
  }
  mg->Draw("a2p");
  double x_max = fPtBins.size()*(fWPs.size()+1.);
  mg->GetHistogram()->GetXaxis()->SetLimits(0., x_max);
  mg->GetHistogram()->GetXaxis()->SetTickLength(0);
  // mg->GetHistogram()->GetXaxis()->SetNdivisions(fPtBins.size(), false);
  mg->GetHistogram()->GetXaxis()->SetLabelOffset(-2.);
  mg->GetHistogram()->SetMinimum(y_sf_low);
  mg->GetHistogram()->SetMaximum((y_sf_high - y_sf_low + msc_y_gap)*fMergeScenarios.size()/(1-header_y_fraction));
  // mg->GetHistogram()->GetYaxis()->SetNdivisions(10021, false);
  mg->GetHistogram()->GetYaxis()->SetTickLength(0);
  mg->GetHistogram()->GetYaxis()->SetLabelOffset(-2.);

  for(unsigned int i = 0; i < y_centers.size(); i++) {
    TLine *line = new TLine(0, y_centers.at(i), x_max, y_centers.at(i));
    line->SetLineColor(kGray+2);
    line->SetLineStyle(2);
    line->Draw();
  }

  for(unsigned int i = 0; i < y_centers.size(); i++) {
    double xmin = 0;
    double ymin = y_centers.at(i)-(y_sf_high-y_sf_low)/2.;
    double xmax = 0;
    double ymax = y_centers.at(i)+(y_sf_high-y_sf_low)/2.;
    double wmin = y_sf_low;
    double wmax = y_sf_high;
    double tick_length = 0.05;
    int ndiv = 505;
    // int ndiv = 503;
    TGaxis *axis_left = new TGaxis(xmin, ymin, xmax, ymax, wmin, wmax, ndiv, "-S");
    axis_left->SetTickLength(tick_length);
    axis_left->SetLabelSize(0.035);
    axis_left->SetLabelFont(42);
    // axis_left->ChangeLabel(1);
    // axis_left->ChangeLabel(-1);
    axis_left->Draw();
    for(unsigned int j = 1; j < fPtBins.size(); j++) {
      xmin = j*(fWPs.size()+1.);
      xmax = j*(fWPs.size()+1.);
      TGaxis *axis = new TGaxis(xmin, ymin, xmax, ymax, wmin, wmax, ndiv, "U+-S");
      axis->SetTickLength(tick_length);
      axis->Draw();
    }
    xmin = x_max;
    xmax = x_max;
    TGaxis *axis_right = new TGaxis(xmin, ymin, xmax, ymax, wmin, wmax, ndiv, "U+S");
    axis_right->SetTickLength(tick_length);
    axis_right->Draw();
  }


  c->cd();

  // double leg_width = 0.08;
  // double leg_height = 0.04;
  // for(const auto & wp : fWPs) {
  //   TLegend *leg = new TLegend();
  //
  // }

  for(unsigned int i = 0; i < fPtBins.size(); i++) {
    PtBin pt_bin = fPtBins.at(i);
    TString pt_text = "";
    if(i == fPtBins.size()-1) pt_text = (TString)"#it{p}_{T}^{jet} > "+my_format(kPtBins.at(pt_bin).pt_min, 0)+" GeV";
    else pt_text = (TString)"#it{p}_{T}^{jet} #in ("+my_format(kPtBins.at(pt_bin).pt_min, 0)+", "+my_format(kPtBins.at(pt_bin).pt_max, 0)+"] GeV";
    TLatex *pt_label = new TLatex(margin_l+(((fWPs.size()+1)*i)/x_max)*(1-margin_l-margin_r)+0.015, 0.9*margin_b, pt_text.Data());
    pt_label->SetTextAlign(13); // left top aligned
    pt_label->SetTextAngle(-15);
    pt_label->SetTextFont(42);
    pt_label->SetTextSize(0.035);
    pt_label->SetNDC();
    pt_label->Draw();
  }

  TLatex *y_label = new TLatex(margin_l*0.45, (1-(msc_y_gap/mg->GetHistogram()->GetMaximum()+header_y_fraction))/2.*(1-margin_t-margin_b)+margin_b, "Data-to-simulation scale factor"); // (#varepsilon_{Data}/#varepsilon_{MC})
  y_label->SetTextAlign(21);
  y_label->SetTextAngle(90);
  y_label->SetTextFont(42);
  y_label->SetTextSize(0.035);
  y_label->SetNDC();
  y_label->Draw();

  // TLatex *text_top_left = new TLatex(margin_l, 1-(margin_t-0.01), kYears.at(year).long_name.Data()); // todo: hotvr puppi
  TLatex *text_top_left = new TLatex(margin_l, 1-(margin_t-0.01), ((TString)"RunIISummer19"+kYears.at(year).short_name+", PUPPI(v15)").Data()); // todo: hotvr puppi
  text_top_left->SetTextAlign(11); // left bottom aligned
  text_top_left->SetTextFont(42);
  text_top_left->SetTextSize(0.035);
  text_top_left->SetNDC();
  text_top_left->Draw();

  // string string_text_top_right = "41.5 fb^{#minus1} (13 TeV)"; // todo: other years have other lumis!
  TLatex *text_top_right = new TLatex(1-margin_r, 1-(margin_t-0.01), (my_format(kYears.at(year).lumi_fb, 1)+" fb^{#minus1} (13 TeV)").c_str());
  text_top_right->SetTextAlign(31); // right bottom aligned
  text_top_right->SetTextFont(42);
  text_top_right->SetTextSize(0.035);
  text_top_right->SetNDC();
  text_top_right->Draw();

  // TText *cms = new TText(cc->ConvertGraphXToPadX(0.03*c->GetWh()/c->GetWw()), cc->ConvertGraphYToPadY(0.97), "CMS");
  TText *cms = new TText(cc->ConvertGraphXToPadX(0.03*c->GetWh()/c->GetWw()), cc->ConvertGraphYToPadY(0.905), "CMS");
  cms->SetTextAlign(11); // left bottom
  cms->SetTextFont(62);
  cms->SetTextSize(0.07); // 0.05
  cms->SetNDC();
  cms->Draw();

  TString prelim_text = (fExpObs == ExpObs::isExp ? (TString)"Simulation, " : (TString)"");
  // prelim_text += "Work in Progress";
  prelim_text += "Internal";
  // TText *prelim = new TText(cc->ConvertGraphXToPadX(0.03*c->GetWh()/c->GetWw()), cc->ConvertGraphYToPadY(0.87), prelim_text.Data());
  TText *prelim = new TText(cc->ConvertGraphXToPadX(0.21*c->GetWh()/c->GetWw()), cc->ConvertGraphYToPadY(0.905), prelim_text.Data());
  prelim->SetTextAlign(11); // left bottom
  prelim->SetTextFont(52);
  // prelim->SetTextSize(0.035);
  prelim->SetTextSize(0.05);
  prelim->SetNDC();
  prelim->Draw();

  TText *algo_label = new TText(cc->ConvertGraphXToPadX(1-(0.03*c->GetWh()/c->GetWw())), cc->ConvertGraphYToPadY(0.97), (kProbeJetAlgos.at(fAlgo).name+" PUPPI").c_str());
  algo_label->SetTextAlign(33); // right top
  algo_label->SetTextFont(62);
  algo_label->SetTextSize(0.05); // 0.05
  algo_label->SetTextColor(kGray+1);
  algo_label->SetNDC();
  algo_label->Draw();

  TString btag_label_text = "";
  TString btag_label2_text = "";
  if(fAlgo == ProbeJetAlgo::isAK8 && (fJetCat == JetCategory::BTag || fJetCat == JetCategory::MassAndBTag)) {
    btag_label_text = "w/ subjet b tag";
    btag_label2_text = "DeepCSV loose";
  }
  else if(fAlgo == ProbeJetAlgo::isHOTVR) {
    btag_label_text = "in t tagging mode";
    btag_label2_text = "w/ #tau_{3}/#tau_{2} < 0.56";
  }
  TLatex *btag_label = new TLatex(cc->ConvertGraphXToPadX(1-(0.03*c->GetWh()/c->GetWw())), cc->ConvertGraphYToPadY(0.885), btag_label_text);
  btag_label->SetTextAlign(33); // right top
  btag_label->SetTextFont(62);
  btag_label->SetTextSize(0.035); // 0.05
  btag_label->SetTextColor(kGray+1);
  btag_label->SetNDC();
  btag_label->Draw();
  TLatex *btag_label2 = new TLatex(cc->ConvertGraphXToPadX(1-(0.03*c->GetWh()/c->GetWw())), cc->ConvertGraphYToPadY(0.822), btag_label2_text);
  btag_label2->SetTextAlign(33); // right top
  btag_label2->SetTextFont(62);
  btag_label2->SetTextSize(0.035); // 0.05
  btag_label2->SetTextColor(kGray+1);
  btag_label2->SetNDC();
  btag_label2->Draw();

  unsigned int msc_i = 0;//fMergeScenarios.size()-1;
  for(const auto & msc : fMergeScenarios) {
    const float pave_distance_to_scale_factor_range = 0.125;
    TPave *pave_msc = new TPave(
      cc->ConvertGraphXToPadX(0.04*c->GetWh()/c->GetWw()),
      cc->ConvertGraphYToPadY((y_centers.at(msc_i)+1.+msc_y_gap-pave_distance_to_scale_factor_range)/mg->GetHistogram()->GetMaximum()),
      // cc->ConvertGraphYToPadY((y_centers.at(msc_i)+1.725)/mg->GetHistogram()->GetMaximum()),
      cc->ConvertGraphXToPadX(0.34*c->GetWh()/c->GetWw()),
      cc->ConvertGraphYToPadY((y_centers.at(msc_i)+1.+pave_distance_to_scale_factor_range)/mg->GetHistogram()->GetMaximum()),
      // cc->ConvertGraphYToPadY((y_centers.at(msc_i)+1.125)/mg->GetHistogram()->GetMaximum()),
      3,
      ""
    );
    // pave_msc->AddText((kMergeScenarios.at(msc).print_name).Data());
    // pave_msc->SetTextAlign(21); // left top
    // pave_msc->SetTextFont(52);
    // pave_msc->SetTextSize(0.035);
    pave_msc->SetLineWidth(0);
    // pave_msc->SetNDC();
    pave_msc->Draw();

    TText *text_msc = new TText(cc->ConvertGraphXToPadX(0.19*c->GetWh()/c->GetWw()), cc->ConvertGraphYToPadY((y_centers.at(msc_i)+1.275)/mg->GetHistogram()->GetMaximum()), (kMergeScenarios.at(msc).print_name).Data());
    text_msc->SetTextAlign(21); // left top
    text_msc->SetTextFont(52);
    text_msc->SetTextSize(0.035);
    // text_msc->SetNDC();
    text_msc->Draw();
    ++msc_i;
  }

  // LEGEND
  if(fAlgo == ProbeJetAlgo::isAK8) {
    float legend_x1 = cc->ConvertGraphXToPadX(0.45*c->GetWh()/c->GetWw());
    float legend_y1 = cc->ConvertGraphYToPadY(0.72);
    float legend_x2 = cc->ConvertGraphXToPadX(1.05*c->GetWh()/c->GetWw());
    float legend_y2 = cc->ConvertGraphYToPadY(0.98);
    TLegend *legend = new TLegend(legend_x1, legend_y1, legend_x2, legend_y2); // x1 y1 x2 y2
    legend->SetTextSize(0.035);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetNColumns(2);

    legend->SetHeader("Working points:");
    for(const auto wp : fWPs) {
      TString leg_text = (TString)"#tau_{3}/#tau_{2} < "+to_str((float)kWorkingPointsAK8.at(wp));
      // legend->AddEntry(sf_tot.at(wp), leg_text);
      legend->AddEntry(sf_stat.at(wp), leg_text);
    }
    legend->Draw();
  }

  // float legend2_x1 = cc->ConvertGraphXToPadX(0.5*c->GetWh()/c->GetWw());
  // float legend2_y1 = cc->ConvertGraphYToPadY(0.72);
  // float legend2_x2 = cc->ConvertGraphXToPadX(1.1*c->GetWh()/c->GetWw());
  // float legend2_y2 = cc->ConvertGraphYToPadY(0.98);
  // TLegend *legend2 = new TLegend(legend_x1, legend_y1, legend_x2, legend_y2); // x1 y1 x2 y2
  // legend2->SetTextSize(0.035);
  // legend2->SetBorderSize(0);
  // legend2->SetFillStyle(0);
  // legend2->SetNColumns(2);
  //
  // // legend2->SetHeader("Working points:");
  // for(const auto wp : fWPs) {
  //   TString leg2_text = (TString)"#tau_{3}/#tau_{2} < "+to_str((float)kWorkingPointsAK8.at(wp));
  //   // legend2->AddEntry(sf_tot.at(wp), leg_text);
  //   legend2->AddEntry(sf_stat.at(wp), leg_text);
  // }
  // legend2->Draw();



  TString plot_path = fWorkDirBase+kYears.at(year).short_name+"/combine/"+kProbeJetAlgos.at(fAlgo).name+"/plots/"
    +"ScaleFactors__"+kProbeJetAlgos.at(fAlgo).name+"_"+kJetCategoryAsString.at(fJetCat)+"__"+kExpObs.at(fExpObs).short_name+"__"+(fUsePtAll ? "PtAll" : "individualFits")+"__"+kYears.at(year).short_name+".pdf";
  c->SaveAs(plot_path);

  delete c;
}

int main(int argc, char **argv) {
  ScaleFactorPlotter *plotter = new ScaleFactorPlotter(ProbeJetAlgo::isAK8, Year::isUL17, JetCategory::All, ExpObs::isObs);
  plotter->ReadScaleFactors();
  plotter->PlotSingleYears();

  ScaleFactorPlotter *plotter1 = new ScaleFactorPlotter(ProbeJetAlgo::isAK8, Year::isUL18, JetCategory::All, ExpObs::isObs);
  plotter1->ReadScaleFactors();
  plotter1->PlotSingleYears();

  ScaleFactorPlotter *plotter2 = new ScaleFactorPlotter(ProbeJetAlgo::isAK8, Year::isUL17, JetCategory::BTag, ExpObs::isObs);
  plotter2->ReadScaleFactors();
  plotter2->PlotSingleYears();

  ScaleFactorPlotter *plotter3 = new ScaleFactorPlotter(ProbeJetAlgo::isAK8, Year::isUL18, JetCategory::BTag, ExpObs::isObs);
  plotter3->ReadScaleFactors();
  plotter3->PlotSingleYears();

  ScaleFactorPlotter *plotter4 = new ScaleFactorPlotter(ProbeJetAlgo::isHOTVR, Year::isUL17, JetCategory::HOTVRCuts, ExpObs::isObs);
  plotter4->ReadScaleFactors();
  plotter4->PlotSingleYears();

  ScaleFactorPlotter *plotter5 = new ScaleFactorPlotter(ProbeJetAlgo::isHOTVR, Year::isUL18, JetCategory::HOTVRCuts, ExpObs::isObs);
  plotter5->ReadScaleFactors();
  plotter5->PlotSingleYears();
}

// int main(int argc, char **argv) {
//   ProbeJetAlgo algo = ProbeJetAlgo::isAK8;
//   // ProbeJetAlgo algo = ProbeJetAlgo::isHOTVR;
//   vector<Year> years = {Year::isUL17, Year::isUL18};
//   JetCategory jet_cat = JetCategory::All;
//   // JetCategory jet_cat = JetCategory::HOTVRCuts;
//   ScaleFactorPlotter *plotter = new ScaleFactorPlotter(algo, years, jet_cat);
//   plotter->ReadScaleFactors();
//   plotter->Plot(Year::isUL17);
//   plotter->Plot(Year::isUL18);
//
//   ProbeJetAlgo algo3 = ProbeJetAlgo::isAK8;
//   // ProbeJetAlgo algo3 = ProbeJetAlgo::isHOTVR;
//   vector<Year> years3 = {Year::isUL17, Year::isUL18};
//   JetCategory jet_cat3 = JetCategory::BTag;
//   // JetCategory jet_cat3 = JetCategory::HOTVRCuts;
//   ScaleFactorPlotter *plotter3 = new ScaleFactorPlotter(algo3, years3, jet_cat3);
//   plotter3->ReadScaleFactors();
//   plotter3->Plot(Year::isUL17);
//   plotter3->Plot(Year::isUL18);
//
//   // ProbeJetAlgo algo2 = ProbeJetAlgo::isAK8;
//   ProbeJetAlgo algo2 = ProbeJetAlgo::isHOTVR;
//   vector<Year> years2 = {Year::isUL17, Year::isUL18};
//   // JetCategory jet_cat2 = JetCategory::All;
//   JetCategory jet_cat2 = JetCategory::HOTVRCuts;
//   ScaleFactorPlotter *plotter2 = new ScaleFactorPlotter(algo2, years2, jet_cat2);
//   plotter2->ReadScaleFactors();
//   plotter2->Plot(Year::isUL17);
//   plotter2->Plot(Year::isUL18);
//   return 0;
// }
