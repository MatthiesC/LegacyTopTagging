#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h"

#include "UHH2/common/include/YearRunSwitchers.h"
#include "UHH2/common/include/JetIds.h"

#include "UHH2/LegacyTopTagging/include/Constants.h"


namespace uhh2 { namespace ltt {

template<typename T, typename U>
double deltaEta(const T & p1, const U & p2) {
    return fabs(p1.eta() - p2.eta());
}

double tau32(const TopJet & topjet);

double tau32groomed(const TopJet & topjet);

double mSD(const TopJet & topjet);

double maxDeepCSVSubJetValue(const TopJet & topjet);

double HOTVR_mpair(const TopJet & topjet, const bool safe = true);

double HOTVR_fpt(const TopJet & topjet, const unsigned int subjet_i = 0);

class HOTVRTopTag {
public:
  explicit HOTVRTopTag(const double _mass_min = 140., const double _mass_max = 220., const double _fpt_max = 0.8, const double _mpair_min = 50.);
  bool operator()(const TopJet & jet, const uhh2::Event & event) const;
private:
  double mass_min;
  double mass_max;
  double fpt_max;
  double mpair_min;
};

const TopJet * nextTopJet(const Particle & p, const std::vector<TopJet> & topjets);

class METSelection: public uhh2::Selection {
public:
  METSelection(const double _met_min = 0., const double _met_max = std::numeric_limits<double>::infinity());
  virtual bool passes(const uhh2::Event & event) override;
private:
  double met_min;
  double met_max;
};

class PTWSelection: public uhh2::Selection {
public:
  PTWSelection(uhh2::Context & ctx, const double _ptw_min = 0., const double _ptw_max = std::numeric_limits<double>::infinity());
  virtual bool passes(const uhh2::Event & event) override;
private:
  double ptw_min;
  double ptw_max;
  uhh2::Event::Handle<FlavorParticle> h_primlep;
};

class TwoDSelection: public uhh2::Selection {
public:
  TwoDSelection(uhh2::Context & ctx, const double _ptrel_min = 0., const double _dr_min = 0.);
  virtual bool passes(const uhh2::Event & event) override;
private:
  double ptrel_min;
  double dr_min;
  uhh2::Event::Handle<FlavorParticle> h_primlep;
};

class BTagCloseToLeptonSelection: public uhh2::Selection {
public:
  BTagCloseToLeptonSelection(uhh2::Context & ctx, const double _dr_max, const JetId & _btagID);
  virtual bool passes(const uhh2::Event & event) override;
private:
  double dr_max;
  JetId btagID;
  uhh2::Event::Handle<FlavorParticle> h_primlep;
};

class MuonScaleFactors: public uhh2::AnalysisModule {
public:
  MuonScaleFactors(uhh2::Context & ctx);
  virtual bool process(uhh2::Event & event) override;
private:
  std::unique_ptr<YearSwitcher> sf_id;
  // std::unique_ptr<YearSwitcher> sf_iso;
};

class TriggerScaleFactors: public uhh2::AnalysisModule {
public:
  TriggerScaleFactors(uhh2::Context & ctx);
  virtual bool process(uhh2::Event & event) override;
private:
  std::unique_ptr<YearSwitcher> sf_trig;
};

class ProbeJetHandleSetter: public uhh2::AnalysisModule {
public:
  ProbeJetHandleSetter(uhh2::Context & ctx, const ProbeJetAlgo & _algo, const std::string & coll_rec = "");
  virtual bool process(uhh2::Event & event) override;
private:
  uhh2::Event::Handle<TopJet> h_probejet;
  uhh2::Event::Handle<std::vector<TopJet>> h_topjets;
};

class DecayChannelAndHadronicTopHandleSetter: public uhh2::AnalysisModule {
public:
  DecayChannelAndHadronicTopHandleSetter(uhh2::Context & ctx);
  virtual bool process(uhh2::Event & event) override;
private:
  uhh2::Event::Handle<DecayChannel> output_decay_channel;
  uhh2::Event::Handle<GenParticle> h_hadronictop;
};

class MergeScenarioHandleSetter: public uhh2::AnalysisModule {
public:
  MergeScenarioHandleSetter(uhh2::Context & ctx, const ProbeJetAlgo & _algo);
  virtual bool process(uhh2::Event & event) override;
private:
  ProbeJetAlgo algo;
  uhh2::Event::Handle<TopJet> h_probejet;
  uhh2::Event::Handle<GenParticle> h_hadronictop;
  uhh2::Event::Handle<bool> output_has_probejet;
  uhh2::Event::Handle<MergeScenario> output_merge_scenario;
};

// This class can be used to setup an AnalysisTree which contains all variables needed to fill probejet histograms at a later stage than SFrame
class MainOutputSetter: public uhh2::AnalysisModule {
public:
  MainOutputSetter(uhh2::Context & ctx);
  virtual bool process(uhh2::Event & event) override;
private:
  uhh2::Event::Handle<TopJet> h_probejet_hotvr;
  uhh2::Event::Handle<TopJet> h_probejet_ak8;
  std::vector<uhh2::Event::Handle<float>> h_mainoutput;
  uhh2::Event::Handle<int> h_probejet_hotvr_nsub_integer;
};

}}
