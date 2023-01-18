#pragma once

#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/LegacyTopTagging/include/Constants.h"
#include "UHH2/LegacyTopTagging/include/Utils.h"

namespace uhh2 { namespace ltt {

//____________________________________________________________________________________________________
class MyTriggerSelection : public uhh2::Selection {
public:
  MyTriggerSelection(uhh2::Context & ctx, const bool low_pt, const bool simplerEleSetup);
  virtual bool passes(const uhh2::Event & event) override;
private:
  const Year fYear;
  const ltt::Channel fChannel;
  const bool fSimplerEleSetup;
  const bool fLowPt;
  enum class DataStream {
    isMC,
    isSingleMuon,
    isSingleElectron,
    isSinglePhoton,
    isEGamma,
  };
  DataStream fDataStream;
  // const double lumi_percentage_UL17_RunB = (41.5 - 36.7) / 41.5; // "note 1) this is only 36.7 out 41.5 fb-1 due to triggers not been included at start up" (https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTRunIISummary#2017)
  const float lumi_percentage_UL17_RunB = 4.803 / 41.480; // = 11.580% ; Christopher's calculation with brilcalc tool
  // brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt [--end 299329] (either with or without the last option)
  bool is_UL17_RunB(const Event & event) const { return fYear == Year::isUL17 && event.run <= 299329; };
  const float lumi_percentage_UL16preVFP_without_TkMu50 = 2.792 / 19.536; // = 14.290% ; Christopher's calculation with brilcalc tool
  // brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON_UL16preVFP.txt [--end 274888]
  bool is_UL16preVFP_without_TkMu50(const Event & event) const { return fYear == Year::isUL16preVFP && event.run < 274889; }; // < is correct, not <=, source: https://twiki.cern.ch/twiki/bin/view/CMS/MuonHLT2016#2016_Runs

  std::unique_ptr<uhh2::Selection> fTrigSel_IsoMu24;
  std::unique_ptr<uhh2::Selection> fTrigSel_IsoTkMu24;
  std::unique_ptr<uhh2::Selection> fTrigSel_IsoMu27;

  std::unique_ptr<uhh2::Selection> fTrigSel_Mu50;
  std::unique_ptr<uhh2::Selection> fTrigSel_TkMu50;
  std::unique_ptr<uhh2::Selection> fTrigSel_OldMu100;
  std::unique_ptr<uhh2::Selection> fTrigSel_TkMu100;

  std::unique_ptr<uhh2::Selection> fTrigSel_Ele27_WPTight_Gsf;
  std::unique_ptr<uhh2::Selection> fTrigSel_Ele35_WPTight_Gsf;
  std::unique_ptr<uhh2::Selection> fTrigSel_Ele32_WPTight_Gsf;

  std::unique_ptr<uhh2::Selection> fTrigSel_Ele115_CaloIdVT_GsfTrkIdT;
  std::unique_ptr<uhh2::Selection> fTrigSel_Photon175;
  std::unique_ptr<uhh2::Selection> fTrigSel_Photon200;
};

}}
