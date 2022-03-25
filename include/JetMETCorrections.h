#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Utils.h"

#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/JetCorrectionSets.h"
#include "UHH2/common/include/YearRunSwitchers.h"


namespace uhh2 { namespace ltt {

class JetMETCorrections: public uhh2::AnalysisModule {
public:
  JetMETCorrections(const std::string & coll_rec = "", const std::string & coll_gen = "");
  virtual bool process(uhh2::Event & event) override;
  void init(uhh2::Context & ctx);
  void switch_jlc(const bool b = true) { fail_if_init_done(); do_jlc = b; }
  void switch_jec(const bool b = true) { fail_if_init_done(); do_jec = b; }
  void switch_met_type1_correction(const bool b = true) { fail_if_init_done(); do_met_type1_correction = b; }
  void switch_met_xy_correction(const bool b = true) { fail_if_init_done(); do_met_xy_correction = b; }

private:
  void fail_if_init_done() const { if(init_done) throw std::runtime_error("JetMETCorrections: Not allowed to call a configuration switch after JetMETCorrections::init() has already been called"); }

  std::string jec_tag_2016, jec_ver_2016, jer_tag_2016;
  std::string jec_tag_2017, jec_ver_2017, jer_tag_2017;
  std::string jec_tag_2018, jec_ver_2018, jer_tag_2018;
  std::string jec_tag_UL16preVFP, jec_ver_UL16preVFP, jer_tag_UL16preVFP;
  std::string jec_tag_UL16postVFP, jec_ver_UL16postVFP, jer_tag_UL16postVFP;
  std::string jec_tag_UL17, jec_ver_UL17, jer_tag_UL17;
  std::string jec_tag_UL18, jec_ver_UL18, jer_tag_UL18;

  std::string collection_rec = "jets";
  bool use_additional_branch_for_rec = false;
  std::string collection_gen = "genjets";
  bool use_additional_branch_for_gen = false;

  bool init_done = false;

  bool do_jlc = true;
  bool do_jec = true;
  bool do_met_type1_correction = true;
  bool do_met_xy_correction = true;

  bool is_mc;
  Year year;

  uhh2::Event::Handle<std::vector<Jet>> h_jets;
  // uhh2::Event::Handle<std::vector<GenJet>> h_genjets;

  std::unique_ptr<AnalysisModule> clnr_jetpfid;

  std::unique_ptr<YearSwitcher> jlc_MC;
  std::unique_ptr<YearSwitcher> jet_corrector_MC;
  std::unique_ptr<GenericJetResolutionSmearer> jet_resolution_smearer;

  std::shared_ptr<RunSwitcher> jlc_switcher_16;
  std::shared_ptr<RunSwitcher> jec_switcher_16;
  std::shared_ptr<RunSwitcher> jlc_switcher_17;
  std::shared_ptr<RunSwitcher> jec_switcher_17;
  std::shared_ptr<RunSwitcher> jlc_switcher_18;
  std::shared_ptr<RunSwitcher> jec_switcher_18;
  std::shared_ptr<RunSwitcher> jlc_switcher_UL16preVFP;
  std::shared_ptr<RunSwitcher> jec_switcher_UL16preVFP;
  std::shared_ptr<RunSwitcher> jlc_switcher_UL16postVFP;
  std::shared_ptr<RunSwitcher> jec_switcher_UL16postVFP;
  std::shared_ptr<RunSwitcher> jlc_switcher_UL17;
  std::shared_ptr<RunSwitcher> jec_switcher_UL17;
  std::shared_ptr<RunSwitcher> jlc_switcher_UL18;
  std::shared_ptr<RunSwitcher> jec_switcher_UL18;
  std::unique_ptr<YearSwitcher> jlc_data;
  std::unique_ptr<YearSwitcher> jet_corrector_data;

  std::unique_ptr<AnalysisModule> met_xy_correction;
};

}}
