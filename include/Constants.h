#pragma once

#include <set>

#include "UHH2/common/include/Utils.h"

namespace uhh2 { namespace ltt {

enum class Channel {
  isEle,
  isMuo,
  notValid,
};

typedef struct {
  int index;
  std::string name;
} ChannelInfo;

const std::map<Channel, ChannelInfo> kChannels = {
  { Channel::isEle,    ChannelInfo{1, "ele"} },
  { Channel::isMuo,    ChannelInfo{2, "muo"} },
  { Channel::notValid, ChannelInfo{0, "invalid"} },
};

typedef struct {
  int index;
  std::string name;
  std::string nice_name;
} YearInfo;

const std::map<Year, YearInfo> kYears = {
  {Year::isUL16preVFP, YearInfo{1, "UL16preVFP", "2016 early"} },
  {Year::isUL16postVFP, YearInfo{2, "UL16postVFP", "2016 late"} },
  {Year::isUL17, YearInfo{3, "UL17", "2017"} },
  {Year::isUL18, YearInfo{4, "UL18", "2018"} },
};

enum class Band {
  MAIN,
  QCD,
};

typedef struct {
  int index;
  std::string name;
} BandInfo;

const std::map<Band, BandInfo> kBands = {
  { Band::MAIN, { .index=0, .name="Main" } },
  { Band::QCD,  { .index=1, .name="QCD" } },
};

const std::set<Band> kRelevantBands = {
  Band::MAIN,
  Band::QCD,
};

const std::string kCollectionName_AK4CHS = "jetsAk4CHS";
// const std::string kHandleName_AK4pairs = "AK4pairs";
const std::string kHandleName_pairedPUPPIjets = "pairedPUPPIjets";
const std::string kHandleName_pairedPUPPIjets_hemi = "pairedPUPPIjets_hemi";
const std::string kHandleName_pairedCHSjets = "pairedCHSjets";
const std::string kHandleName_pairedCHSjets_hemi = "pairedCHSjets_hemi";
const double kDeltaRForPuppiCHSMatch = 0.2;
const double kAbsEtaBTagThreshold = 2.5;
const std::string kHandleName_forwardPUPPIjets = "forwardPUPPIjets";
const std::string kHandleName_uncleanedPUPPIjets = "uncleanedPUPPIjets";
const std::string kCollectionName_METCHS = "slimmedMETs";
const std::string kCollectionName_METPUPPI = "slimmedMETsPuppi";

const std::string kHandleName_weight_btag_njet_sf = "weight_btag_njet_sf";

const std::string kCollectionName_AK8_rec = "jetsAk8PuppiSubstructure_SoftDropPuppi";
const std::string kCollectionName_AK8_gen = "genjetsAk8SubstructureSoftDrop";

const std::string kHandleName_PrimaryLepton = "PrimaryLepton";
const double kDeltaRLeptonicHemisphere = M_PI*2./3.;

const std::string kHandleName_bJets = "bJets";
const std::string kHandleName_bJets_loose = "bJets_loose";
const std::string kHandleName_bJets_medium = "bJets_medium";
const std::string kHandleName_bJets_tight = "bJets_tight";

const std::string kHandleName_n_jets = "n_jets";
const std::string kHandleName_n_jets_central = "n_jets_central";
const std::string kHandleName_n_jets_forward = "n_jets_forward";

const std::string kHandleName_n_bJets = "n_bJets";
const std::string kHandleName_n_bJets_loose = "n_bJets_loose";
const std::string kHandleName_n_bJets_medium = "n_bJets_medium";
const std::string kHandleName_n_bJets_tight = "n_bJets_tight";

const std::string kHandleName_bJets_hemi = "bJets_hemi";
const std::string kHandleName_bJets_hemi_loose = "bJets_hemi_loose";
const std::string kHandleName_bJets_hemi_medium = "bJets_hemi_medium";
const std::string kHandleName_bJets_hemi_tight = "bJets_hemi_tight";

const std::string kHandleName_n_bJets_hemi = "n_bJets_hemi";
const std::string kHandleName_n_bJets_hemi_loose = "n_bJets_hemi_loose";
const std::string kHandleName_n_bJets_hemi_medium = "n_bJets_hemi_medium";
const std::string kHandleName_n_bJets_hemi_tight = "n_bJets_hemi_tight";

const std::string kHandleName_SingleTopGen_tWch = "SingleTopGen_tWch";


const std::string k_jec_tag_2016 = "Summer16_07Aug2017";
const std::string k_jec_ver_2016 = "11";
const std::string k_jer_tag_2016 = "Summer16_25nsV1";

const std::string k_jec_tag_2017 = "Fall17_17Nov2017";
const std::string k_jec_ver_2017 = "32";
const std::string k_jer_tag_2017 = "Fall17_V3";

const std::string k_jec_tag_2018 = "Autumn18";
const std::string k_jec_ver_2018 = "19";
const std::string k_jer_tag_2018 = "Autumn18_V7";

const std::string k_jec_tag_UL16preVFP = "Summer19UL16APV";
const std::string k_jec_ver_UL16preVFP = "7";
const std::string k_jer_tag_UL16preVFP = "Summer20UL16APV_JRV3";

const std::string k_jec_tag_UL16postVFP = "Summer19UL16";
const std::string k_jec_ver_UL16postVFP = "7";
const std::string k_jer_tag_UL16postVFP = "Summer20UL16_JRV3";

const std::string k_jec_tag_UL17 = "Summer19UL17";
const std::string k_jec_ver_UL17 = "5";
const std::string k_jer_tag_UL17 = "Summer19UL17_JRV2";

const std::string k_jec_tag_UL18 = "Summer19UL18";
const std::string k_jec_ver_UL18 = "5";
const std::string k_jer_tag_UL18 = "Summer19UL18_JRV2";

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
};

typedef struct {
  int index;
  std::string name;
} MergeScenarioInfo;

const std::map<MergeScenario, MergeScenarioInfo> kMergeScenarios = {
  {MergeScenario::isFullyMerged, { .index=1,  .name="FullyMerged"} },
  {MergeScenario::isWMerged,     { .index=2,  .name="WMerged"} },
  {MergeScenario::isQBMerged,    { .index=3,  .name="QBMerged"} },
  {MergeScenario::isNotMerged,   { .index=4,  .name="NotMerged"} },
  {MergeScenario::isBackground,  { .index=-1, .name="Background"} },
  {MergeScenario::isAll,         { .index=0,  .name="AllMergeScenarios"} },
};

enum class DecayChannel {
  isTTbarToHadronic,
  isTTbarToSemiLeptonic,
  isTTbarToDiLeptonic,
  isSTToHadronic,
  isSTToLeptonic,
  isNotTopQuarkMC,
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
  PtBin::Pt200toInf,
  PtBin::Pt300toInf,
  PtBin::Pt400toInf,
};

const std::vector<PtBin> kPtBinsAK8 = {
  PtBin::Pt300to400,
  PtBin::Pt400to480,
  PtBin::Pt480to600,
  PtBin::Pt600toInf,
  // only used for control distributions etc.:
  PtBin::Pt300toInf,
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

const std::map<WorkingPoint, std::string> kWorkingPointAsString = {
  {WorkingPoint::Standard, "Standard"},
  {WorkingPoint::NoTau32Cut, "NoTau32Cut"},
  {WorkingPoint::BkgEff0p001, "BkgEff0p001"},
  {WorkingPoint::BkgEff0p003, "BkgEff0p003"},
  {WorkingPoint::BkgEff0p005, "BkgEff0p005"},
  {WorkingPoint::BkgEff0p010, "BkgEff0p010"},
  {WorkingPoint::BkgEff0p025, "BkgEff0p025"},
  {WorkingPoint::BkgEff0p030, "BkgEff0p030"},
  {WorkingPoint::BkgEff0p050, "BkgEff0p050"},
  {WorkingPoint::BkgEff0p100, "BkgEff0p100"},
};

typedef std::map<WorkingPoint, double> WorkingPointMap;

const WorkingPointMap kWorkingPointsHOTVR = {
  {WorkingPoint::Standard, 0.56},
  // only used for control distributions etc.:
  {WorkingPoint::NoTau32Cut, std::numeric_limits<double>::infinity()},
};

// Working points as of https://indico.cern.ch/event/1033067/contributions/4338454/attachments/2238372/3794612/JMAR_update_2021-05-04.pdf
const WorkingPointMap kWorkingPointsAK8 = {
  {WorkingPoint::BkgEff0p001, 0.38},
  {WorkingPoint::BkgEff0p005, 0.47},
  {WorkingPoint::BkgEff0p010, 0.52},
  {WorkingPoint::BkgEff0p025, 0.61},
  {WorkingPoint::BkgEff0p050, 0.69},
  // only used for control distributions etc.:
  {WorkingPoint::NoTau32Cut, std::numeric_limits<double>::infinity()},
};

const WorkingPointMap kWorkingPointsAK8_old = {
  {WorkingPoint::BkgEff0p001, 0.40},
  {WorkingPoint::BkgEff0p003, 0.46},
  {WorkingPoint::BkgEff0p010, 0.54},
  {WorkingPoint::BkgEff0p030, 0.65},
  {WorkingPoint::BkgEff0p100, 0.80},
  // only used for control distributions etc.:
  {WorkingPoint::NoTau32Cut, std::numeric_limits<double>::infinity()},
};

const double kWorkingPointVariation = 0.02;

const double kTau21Cut = 0.45;
const double kTau21Variation = 0.02;

}}
