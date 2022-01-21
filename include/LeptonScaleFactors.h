#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h"
#include "UHH2/core/include/Electron.h"
#include "UHH2/core/include/Muon.h"

#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/YearRunSwitchers.h"


namespace uhh2 { namespace ltt {

// For reference on EGamma scale factors in UL:
// https://twiki.cern.ch/twiki/bin/view/CMS/EgammaUL2016To2018

//____________________________________________________________________________________________________
// Note that there are different ROOT files for ele_pt in interval 10-20 GeV; right now, this class always applies the scale factors for ele_pt > 20 GeV even if the actual ele_pt is below that value
class ElectronRecoScaleFactors: public uhh2::AnalysisModule {
public:
  ElectronRecoScaleFactors(uhh2::Context & ctx, const bool dummy = false);
  virtual bool process(uhh2::Event & event) override;
private:
  const bool fDummy;
  std::vector<uhh2::Event::Handle<float>> fHandles;
  std::unique_ptr<YearSwitcher> fYearSwitcher;
};

//____________________________________________________________________________________________________
class ElectronIdScaleFactors: public uhh2::AnalysisModule {
public:
  ElectronIdScaleFactors(uhh2::Context & ctx, const Electron::tag & tagID, const bool dummy = false);
  virtual bool process(uhh2::Event & event) override;
private:
  const bool fDummy;
  std::vector<uhh2::Event::Handle<float>> fHandles;
  std::unique_ptr<YearSwitcher> fYearSwitcher;
};


// Twiki: https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL201{6,7,8}
// SF files: https://gitlab.cern.ch/cms-muonPOG/muonefficiencies/-/tree/master/Run2/UL

//____________________________________________________________________________________________________
class MuonIdScaleFactors: public uhh2::AnalysisModule {
public:
  MuonIdScaleFactors(uhh2::Context & ctx, const Muon::Selector & selectorID, const bool dummy = false);
  virtual bool process(uhh2::Event & event) override;
private:
  const bool fDummy;
  std::vector<uhh2::Event::Handle<float>> fHandles;
  std::unique_ptr<YearSwitcher> fYearSwitcher;
};

//____________________________________________________________________________________________________
class MuonIsoScaleFactors: public uhh2::AnalysisModule {
public:
  MuonIsoScaleFactors(uhh2::Context & ctx, const Muon::Selector & selectorISO, const Muon::Selector & selectorID, const bool dummy = false);
  virtual bool process(uhh2::Event & event) override;
private:
  const bool fDummy;
  std::vector<uhh2::Event::Handle<float>> fHandles;
  std::unique_ptr<YearSwitcher> fYearSwitcher;
};

//____________________________________________________________________________________________________
class MuonTriggerScaleFactors: public uhh2::AnalysisModule {
public:
  MuonTriggerScaleFactors(uhh2::Context & ctx, const bool use_Mu50, const bool absolute_eta = false, const bool dummy = false);
  virtual bool process(uhh2::Event & event) override;
private:
  const bool fDummy;
  std::vector<uhh2::Event::Handle<float>> fHandles;
  std::unique_ptr<YearSwitcher> fYearSwitcher;
};

}}
