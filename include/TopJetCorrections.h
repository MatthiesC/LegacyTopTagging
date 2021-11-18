#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Utils.h"

#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/JetCorrectionSets.h"
#include "UHH2/common/include/YearRunSwitchers.h"


namespace uhh2 { namespace ltt {

class TopJetCorrections: public uhh2::AnalysisModule {
public:
  TopJetCorrections(const std::string & coll_rec = "", const std::string & coll_gen = "");
  virtual bool process(uhh2::Event & event) override;
  void init(uhh2::Context & ctx);
  void switch_topjet_corrections(const bool b = true) { fail_if_init_done(); correct_topjets = b; }
  void switch_subjet_corrections(const bool b = true) { fail_if_init_done(); correct_subjets = b; }
  void switch_rebuilding_topjets_from_subjets(const bool b = true) { fail_if_init_done(); do_rebuild_topjets_from_subjets = b; }

private:
  void set_subjet_handles(uhh2::Event & event);
  void reset_smeared_subjets(uhh2::Event & event);
  void rebuild_topjets_from_subjets(uhh2::Event & event);
  void fail_if_init_done() const { if(init_done) throw std::runtime_error("TopJetCorrections: Not allowed to call a configuration switch after TopJetCorrections::init() has already been called"); }

  std::unique_ptr<YearSwitcher> tjet_corrector_MC, tjet_corrector_data, subjet_corrector_MC, subjet_corrector_data;
  std::shared_ptr<RunSwitcher> jec_switcher_16_topjets, jec_switcher_16_subjets;
  std::shared_ptr<RunSwitcher> jec_switcher_17_topjets, jec_switcher_17_subjets;
  std::shared_ptr<RunSwitcher> jec_switcher_18_topjets, jec_switcher_18_subjets;
  std::shared_ptr<RunSwitcher> jec_switcher_UL16preVFP_topjets, jec_switcher_UL16preVFP_subjets;
  std::shared_ptr<RunSwitcher> jec_switcher_UL16postVFP_topjets, jec_switcher_UL16postVFP_subjets;
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
  std::string jec_tag_UL16preVFP, jec_ver_UL16preVFP, jer_tag_UL16preVFP;
  std::string jec_tag_UL16postVFP, jec_ver_UL16postVFP, jer_tag_UL16postVFP;
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


class TopJetCleaning : public uhh2::AnalysisModule {
 public:
  TopJetCleaning(uhh2::Context & ctx, const double pt_min, const double eta_max, const double dr_min, const std::string & coll_rec = "");
  virtual bool process(uhh2::Event & event) override;

 private:
  double pt_min;
  double eta_max;
  double dr_min;

  uhh2::Event::Handle<FlavorParticle> h_primlep;
  uhh2::Event::Handle<std::vector<TopJet>> h_topjets;
};

}}
