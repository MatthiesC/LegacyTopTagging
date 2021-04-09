#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Utils.h"

#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/JetCorrectionSets.h"
#include "UHH2/common/include/YearRunSwitchers.h"


namespace uhh2 { namespace ltt {

class TopJetCorrections: public uhh2::AnalysisModule {
public:
  TopJetCorrections(const std::string & coll_rec="", const std::string & coll_gen="");
  virtual bool process(uhh2::Event & event) override;
  void init(uhh2::Context & ctx);
  void disable_topjet_corrections() { fail_if_init_done(); correct_topjets = false; }
  void disable_subjet_corrections() { fail_if_init_done(); correct_subjets = false; }
  void enable_rebuilding_topjets_from_subjets() { fail_if_init_done(); do_rebuild_topjets_from_subjets = true; }

private:
  void set_subjet_handles(uhh2::Event & event);
  void reset_smeared_subjets(uhh2::Event & event);
  void rebuild_topjets_from_subjets(uhh2::Event & event);
  void fail_if_init_done() const { if(init_done) throw std::runtime_error("TopJetCorrections: Not allowed to call a configuration switch after TopJetCorrections::init() has already been called"); }

  std::unique_ptr<YearSwitcher> tjet_corrector_MC, tjet_corrector_data, subjet_corrector_MC, subjet_corrector_data;
  std::shared_ptr<RunSwitcher> jec_switcher_16_topjets, jec_switcher_16_subjets;
  std::shared_ptr<RunSwitcher> jec_switcher_17_topjets, jec_switcher_17_subjets;
  std::shared_ptr<RunSwitcher> jec_switcher_18_topjets, jec_switcher_18_subjets;
  std::shared_ptr<RunSwitcher> jec_switcher_UL17_topjets, jec_switcher_UL17_subjets;
  std::shared_ptr<RunSwitcher> jec_switcher_UL18_topjets, jec_switcher_UL18_subjets;
  std::unique_ptr<GenericJetResolutionSmearer> tjet_resolution_smearer;
  std::unique_ptr<GenericJetResolutionSmearer> subjet_resolution_smearer;

  bool is_mc;
  bool init_done = false;

  bool correct_topjets = true;
  bool correct_subjets = true;
  bool do_rebuild_topjets_from_subjets = false;

  Year year;

  std::string jec_tag_2016, jec_ver_2016, jer_tag_2016;
  std::string jec_tag_2017, jec_ver_2017, jer_tag_2017;
  std::string jec_tag_2018, jec_ver_2018, jer_tag_2018;
  std::string jec_tag_UL17, jec_ver_UL17, jer_tag_UL17;
  std::string jec_tag_UL18, jec_ver_UL18, jer_tag_UL18;

  std::string collection_rec = "topjets";
  bool use_additional_branch_for_rec = false;
  std::string collection_gen = "gentopjets";
  bool use_additional_branch_for_gen = false;

  uhh2::Event::Handle<std::vector<TopJet>> h_topjets;
  uhh2::Event::Handle<std::vector<GenTopJet>> h_gentopjets;

  uhh2::Event::Handle<std::vector<GenJet>> h_gensubjets;
  uhh2::Event::Handle<std::vector<Jet>> h_subjets;
  uhh2::Event::Handle<std::vector<std::pair<int, int>>> h_subjets_map;
};


class AK8Cleaning : public uhh2::AnalysisModule {
 public:
  AK8Cleaning(uhh2::Context & ctx, const double minPt, const double maxEta, const double minDR = 0.8, const std::string & coll_rec="", const std::string & h_name_primlep = "PrimaryLepton");
  virtual bool process(uhh2::Event&) override;
  void switch_jet_pt_sorter(const bool b=true) { jet_pt_sorter=b; };

 private:
  double _minPt;
  double _maxEta;
  double _minDR;
  bool jet_pt_sorter = true;

  uhh2::Event::Handle<FlavorParticle> h_primlep;
  uhh2::Event::Handle<std::vector<TopJet>> h_ak8jets_rec;
};

}}
