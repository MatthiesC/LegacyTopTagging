#pragma once

namespace uhh2 { namespace ltt {

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

enum class MergeScenario {
  isFullyMerged,
  isWMerged,
  isQBMerged,
  isNotMerged,
  isBackgroundMCorData,
};

const std::map<MergeScenario, std::string> kMergeScenarioAsString = {
  {MergeScenario::isFullyMerged, "FullyMerged"},
  {MergeScenario::isWMerged, "WMerged"},
  {MergeScenario::isQBMerged, "QBMerged"},
  {MergeScenario::isNotMerged, "NotMerged"},
  {MergeScenario::isBackgroundMCorData, "BackgroundMCorData"},
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

const std::map<PtBin, std::string> kPtBinAsString = {
  {PtBin::Pt200to250, "Pt200to250"},
  {PtBin::Pt250to300, "Pt250to300"},
  {PtBin::Pt300to400, "Pt300to400"},
  {PtBin::Pt400to480, "Pt400to480"},
  {PtBin::Pt480to600, "Pt480to600"},
  {PtBin::Pt600toInf, "Pt600toInf"},
  {PtBin::Pt200toInf, "Pt200toInf"},
  {PtBin::Pt300toInf, "Pt300toInf"},
  {PtBin::Pt400toInf, "Pt400toInf"},
};

typedef std::map<PtBin, std::pair<double, double>> PtBinMap;

const PtBinMap kPtBinsHOTVR = {
  {PtBin::Pt200to250, {200, 250}},
  {PtBin::Pt250to300, {250, 300}},
  {PtBin::Pt300to400, {300, 400}},
  {PtBin::Pt400to480, {400, 480}},
  {PtBin::Pt480to600, {480, 600}},
  {PtBin::Pt600toInf, {600, std::numeric_limits<double>::infinity()}},
  // only used for control distributions etc.:
  {PtBin::Pt200toInf, {200, std::numeric_limits<double>::infinity()}},
  {PtBin::Pt300toInf, {300, std::numeric_limits<double>::infinity()}},
  {PtBin::Pt400toInf, {400, std::numeric_limits<double>::infinity()}},
};

const PtBinMap kPtBinsAK8 = {
  {PtBin::Pt300to400, {300, 400}},
  {PtBin::Pt400to480, {400, 480}},
  {PtBin::Pt480to600, {480, 600}},
  {PtBin::Pt600toInf, {600, std::numeric_limits<double>::infinity()}},
  // only used for control distributions etc.:
  {PtBin::Pt300toInf, {300, std::numeric_limits<double>::infinity()}},
  {PtBin::Pt400toInf, {400, std::numeric_limits<double>::infinity()}},
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

}}
