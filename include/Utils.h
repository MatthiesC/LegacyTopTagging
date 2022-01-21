#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h"

#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/common/include/Utils.h"

#include "UHH2/LegacyTopTagging/include/Constants.h"


namespace uhh2 { namespace ltt {

//____________________________________________________________________________________________________
template<typename T, typename U>
double deltaEta(const T & p1, const U & p2) {
    return fabs(p1.eta() - p2.eta());
}

//____________________________________________________________________________________________________
double tau32(const TopJet & topjet);

//____________________________________________________________________________________________________
double tau21(const TopJet & topjet);

//____________________________________________________________________________________________________
double tau32groomed(const TopJet & topjet);

//____________________________________________________________________________________________________
double tau21groomed(const TopJet & topjet);

//____________________________________________________________________________________________________
double mSD(const TopJet & topjet);

//____________________________________________________________________________________________________
double maxDeepCSVSubJetValue(const TopJet & topjet);

//____________________________________________________________________________________________________
double maxDeepJetSubJetValue(const TopJet & topjet);

//____________________________________________________________________________________________________
double HOTVR_mpair(const TopJet & topjet, const bool safe = true);

//____________________________________________________________________________________________________
double HOTVR_fpt(const TopJet & topjet, const unsigned int subjet_i = 0);

//____________________________________________________________________________________________________
class HOTVRTopTag {
public:
  explicit HOTVRTopTag(const double _mass_min = kProbeJetAlgos.at(ProbeJetAlgo::isHOTVR).mass_min, const double _mass_max = kProbeJetAlgos.at(ProbeJetAlgo::isHOTVR).mass_max, const double _fpt_max = 0.8, const double _mpair_min = 50.);
  bool operator()(const TopJet & jet, const uhh2::Event & event) const;
private:
  const double mass_min;
  const double mass_max;
  const double fpt_max;
  const double mpair_min;
};

//____________________________________________________________________________________________________
const TopJet * nextTopJet(const Particle & p, const std::vector<TopJet> & topjets);

//____________________________________________________________________________________________________
class METSelection: public uhh2::Selection {
public:
  METSelection(const double _met_min = 0., const double _met_max = std::numeric_limits<double>::infinity());
  virtual bool passes(const uhh2::Event & event) override;
private:
  const double met_min;
  const double met_max;
};

//____________________________________________________________________________________________________
class PTWSelection: public uhh2::Selection {
public:
  PTWSelection(uhh2::Context & ctx, const double _ptw_min = 0., const double _ptw_max = std::numeric_limits<double>::infinity());
  virtual bool passes(const uhh2::Event & event) override;
private:
  const double ptw_min;
  const double ptw_max;
  uhh2::Event::Handle<FlavorParticle> h_primlep;
};

//____________________________________________________________________________________________________
class TwoDSelection: public uhh2::Selection {
public:
  TwoDSelection(uhh2::Context & ctx, const double _ptrel_min = 0., const double _dr_min = 0.);
  virtual bool passes(const uhh2::Event & event) override;
private:
  const double ptrel_min;
  const double dr_min;
  uhh2::Event::Handle<FlavorParticle> h_primlep;
};

//____________________________________________________________________________________________________
class BTagCloseToLeptonSelection: public uhh2::Selection {
public:
  BTagCloseToLeptonSelection(uhh2::Context & ctx, const double _dr_max, const JetId & _btagID);
  virtual bool passes(const uhh2::Event & event) override;
private:
  const double dr_max;
  const JetId btagID;
  uhh2::Event::Handle<FlavorParticle> h_primlep;
};

//____________________________________________________________________________________________________
class PartonShowerVariation: public uhh2::AnalysisModule {
public:
  PartonShowerVariation(uhh2::Context & ctx);
  virtual bool process(uhh2::Event & event) override;
private:
  enum class PSVariation {
    None,
    ISRup_sqrt2,
    ISRup_2,
    ISRup_4,
    ISRdown_sqrt2,
    ISRdown_2,
    ISRdown_4,
    FSRup_sqrt2,
    FSRup_2,
    FSRup_4,
    FSRdown_sqrt2,
    FSRdown_2,
    FSRdown_4,
  };
  typedef struct {
    std::string name = "";
    unsigned int index = 0;
  } PSVariationInfo;
  const std::map<PSVariation, PSVariationInfo> kPSVariations = {
    {PSVariation::ISRup_sqrt2, PSVariationInfo{"ISRup_sqrt2", 2}},
    {PSVariation::ISRup_2, PSVariationInfo{"ISRup_2", 6}},
    {PSVariation::ISRup_4, PSVariationInfo{"ISRup_4", 10}},
    {PSVariation::ISRdown_sqrt2, PSVariationInfo{"ISRdown_sqrt2", 4}},
    {PSVariation::ISRdown_2, PSVariationInfo{"ISRdown_2", 8}},
    {PSVariation::ISRdown_4, PSVariationInfo{"ISRdown_4", 12}},
    {PSVariation::FSRup_sqrt2, PSVariationInfo{"FSRup_sqrt2", 3}},
    {PSVariation::FSRup_2, PSVariationInfo{"FSRup_2", 7}},
    {PSVariation::FSRup_4, PSVariationInfo{"FSRup_4", 11}},
    {PSVariation::FSRdown_sqrt2, PSVariationInfo{"FSRdown_sqrt2", 5}},
    {PSVariation::FSRdown_2, PSVariationInfo{"FSRdown_2", 9}},
    {PSVariation::FSRdown_4, PSVariationInfo{"FSRdown_4", 13}},
  };
  bool warning_thrown = false;
  PSVariation applied_variation = PSVariation::None;
  std::map<PSVariation, uhh2::Event::Handle<double>> weights_map;
};

class MuonScaleFactors: public uhh2::AnalysisModule {
public:
  MuonScaleFactors(uhh2::Context & ctx);
  virtual bool process(uhh2::Event & event) override;
private:
  std::unique_ptr<uhh2::AnalysisModule> sf_id;
  // std::unique_ptr<uhh2::AnalysisModule> sf_iso;
};

//____________________________________________________________________________________________________
class LttTriggerSelection: public uhh2::Selection {
public:
  LttTriggerSelection(const uhh2::Context & ctx);
  virtual bool passes(const uhh2::Event & event) override;
private:
  const Year fYear;
  std::unique_ptr<uhh2::Selection> slct_trigger_Mu50;
  std::unique_ptr<uhh2::Selection> slct_trigger_TkMu50;
  std::unique_ptr<uhh2::Selection> slct_trigger_OldMu100;
  // std::unique_ptr<uhh2::Selection> slct_trigger_Mu100;
  std::unique_ptr<uhh2::Selection> slct_trigger_TkMu100;
};

//____________________________________________________________________________________________________
class TriggerScaleFactors: public uhh2::AnalysisModule {
public:
  TriggerScaleFactors(uhh2::Context & ctx);
  virtual bool process(uhh2::Event & event) override;
private:
  std::unique_ptr<uhh2::AnalysisModule> sf_trig;
};

//____________________________________________________________________________________________________
class ProbeJetHandleSetter: public uhh2::AnalysisModule {
public:
  ProbeJetHandleSetter(uhh2::Context & ctx, const ProbeJetAlgo & _algo, const std::string & coll_rec = "");
  virtual bool process(uhh2::Event & event) override;
private:
  uhh2::Event::Handle<TopJet> h_probejet;
  uhh2::Event::Handle<std::vector<TopJet>> h_topjets;
};

//____________________________________________________________________________________________________
void get_Wb_daughters(GenParticle & w_from_top, GenParticle & b_from_top, const GenParticle & top, const std::vector<GenParticle> & genparticles);

//____________________________________________________________________________________________________
class DecayChannelAndHadronicTopHandleSetter: public uhh2::AnalysisModule {
public:
  DecayChannelAndHadronicTopHandleSetter(uhh2::Context & ctx);
  virtual bool process(uhh2::Event & event) override;
private:
  uhh2::Event::Handle<DecayChannel> output_decay_channel;
  uhh2::Event::Handle<GenParticle> h_hadronictop;
};

//____________________________________________________________________________________________________
class MergeScenarioHandleSetter: public uhh2::AnalysisModule {
public:
  MergeScenarioHandleSetter(uhh2::Context & ctx, const ProbeJetAlgo & _algo);
  virtual bool process(uhh2::Event & event) override;
private:
  const ProbeJetAlgo algo;
  uhh2::Event::Handle<TopJet> h_probejet;
  uhh2::Event::Handle<GenParticle> h_hadronictop;
  uhh2::Event::Handle<bool> output_has_probejet;
  uhh2::Event::Handle<MergeScenario> output_merge_scenario;
};

//____________________________________________________________________________________________________
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

//____________________________________________________________________________________________________
// https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/2000.html
class HEM2018Selection: public uhh2::Selection {
public:
  HEM2018Selection(uhh2::Context & ctx);
  virtual bool passes(const uhh2::Event & event) override;
  double GetAffectedLumiFraction() const { return fAffectedLumiFraction; };
private:
  const Year fYear;
  const int fRunNumber = 319077;
  const std::pair<double, double> fEtaRange = {-3.2, -1.3};
  const std::pair<double, double> fPhiRange = {-1.57, -0.87};
  const double fAffectedLumiFraction = 0.64844705699; // (Run 319077 (17.370008/pb) + Run C + Run D) / all 2018
};

//____________________________________________________________________________________________________
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1PrefiringWeightRecipe
// (superseeds https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1ECALPrefiringWeightRecipe)
class PrefiringWeights: public uhh2::AnalysisModule {
public:
  PrefiringWeights(uhh2::Context & ctx, const bool apply = true);
  virtual bool process(uhh2::Event & event) override;
private:
  const Year fYear;
  uhh2::Event::Handle<double> h_weight_nominal;
  uhh2::Event::Handle<double> h_weight_up;
  uhh2::Event::Handle<double> h_weight_down;
  void set_dummy_weights(uhh2::Event & event);
  const double fDummyWeight = 1.;
  const bool fApply;
  enum class PrefireVariation {
    nominal,
    up,
    down,
  };
  PrefireVariation applied_variation;
};

//____________________________________________________________________________________________________
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
class TopPtReweighting: public uhh2::AnalysisModule {
public:
  TopPtReweighting(uhh2::Context & ctx, const bool apply = true);
  virtual bool process(uhh2::Event & event) override;
private:
  uhh2::Event::Handle<double> h_weight_nominal;
  uhh2::Event::Handle<double> h_weight_a_up;
  uhh2::Event::Handle<double> h_weight_a_down;
  uhh2::Event::Handle<double> h_weight_b_up;
  uhh2::Event::Handle<double> h_weight_b_down;
  void set_dummy_weights(uhh2::Event & event);
  const double fDummyWeight = 1.;
  const bool fApply;
  enum class TopPtVariation {
    nominal,
    a_up,
    a_down,
    b_up,
    b_down,
  };
  TopPtVariation applied_variation;
  const bool fPtCutOff_b = false;
  const double fPtCutOff = 500.;
  const double fA = 0.0615;
  const double fA_up = fA*1.5;
  const double fA_down = fA*0.5;
  const double fB = -0.0005;
  const double fB_up = fB*1.5;
  const double fB_down = fB*0.5;
};

//____________________________________________________________________________________________________
// Copy of https://github.com/MatthiesC/HighPtSingleTop/blob/master/include/TheoryCorrections.h
/*
Taken from https://github.com/UHH2/VHResonances/blob/master/Analysis/python/PlotNLOCorrections.py
- EWK corrections are available for W+jets, Z+jets, gamma+jets samples in the "merged_kfactors_*.root" files under the name "kfactor_monojet_ewk"
- QCD NLO corrections are available for W+jets, Z+jets, gamma+jets in the same files under the name "kfactor_monojet_qcd". Those are calculated for 2016 samples.
- 2017 version of QCD NLO corrections are available for Z+jets (ll + nunu cases) in the "kfac_*_filter" files.
- QCD NNLO corrections are in the "lindert_qcd_nnlo_sf" file with the following convention:
    - eej -> Z(ll) +jets
    - vvj -> Z(nunu) +jets
    - evj -> W +jets
    - aj -> gamma +jets
- QCD NNLO corrections need to be applied on top of EWK corrections for NLO samples and on top of EWK + QCD NLO corrections for LO samples.
- According to Andreas "I do not apply the NNLO corrections. I have not seen any evidence that they actually improve data/MC agreement. I do not trust them."
- For W+Jets @LO for 2017 and 2018: wjet_dress_monojet or wjet_dress_inclusive in "2017_gen_v_pt_qcd_sf.root"
- In the "merged_kfactors_*.root" file, for Z and W + jets, the qcd_ewk histograms are also present: qcd_ewk = QCD * EWK
- taken from https://github.com/bu-cms/bucoffea/tree/master/bucoffea/data/sf/theory
- relative to those studies https://arxiv.org/abs/1705.04664
Example module can be found here: https://github.com/UHH2/VHResonances/blob/master/src/HiggsToWWModules.cxx
*/

class VJetsReweighting: public uhh2::AnalysisModule {
 public:
  explicit VJetsReweighting(uhh2::Context & ctx, const std::string& weight_name="weight_vjets");
  virtual bool process(uhh2::Event & event) override;

 private:
  std::unordered_map<std::string, std::unique_ptr<TH1F>> histos;

  void load_histo(TFile* file, const std::string& name, const std::string& histName);
  double get_v_pt(uhh2::Event & event);
  double evaluate(const std::string& name, double pt);

  const bool is_2016_nonUL;
  const bool is_WJets;
  const bool is_DYJets;
  const bool apply_EWK;
  const bool apply_QCD_EWK;
  const bool apply_QCD_NLO;
  const bool apply_QCD_NNLO;

  uhh2::Event::Handle<double> h_weight_EWK;
  uhh2::Event::Handle<double> h_weight_QCD_EWK;
  uhh2::Event::Handle<double> h_weight_QCD_NLO;
  uhh2::Event::Handle<double> h_weight_QCD_NNLO;
};

}}
