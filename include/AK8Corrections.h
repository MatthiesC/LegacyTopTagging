#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Utils.h"

#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/JetCorrectionSets.h"
#include "UHH2/common/include/YearRunSwitchers.h"


namespace ltt {

class AK8Corrections: public uhh2::AnalysisModule {
public:
  AK8Corrections(const std::string & coll_rec="", const std::string & coll_gen="");
  virtual bool process(uhh2::Event & event) override;
  void init(uhh2::Context & ctx);

private:
  std::unique_ptr<YearSwitcher> tjet_corrector_MC, tjet_corrector_data;
  std::shared_ptr<RunSwitcher> jec_switcher_16, jec_switcher_17, jec_switcher_18, jec_switcher_UL17, jec_switcher_UL18;
  std::unique_ptr<GenericJetResolutionSmearer> tjet_resolution_smearer;

  bool is_mc;
  bool init_done = false;

  Year year;

  std::string jec_tag_2016, jec_ver_2016;
  std::string jec_tag_2017, jec_ver_2017;
  std::string jec_tag_2018, jec_ver_2018;
  std::string jec_tag_UL17, jec_ver_UL17;
  std::string jec_tag_UL18, jec_ver_UL18;
  std::string jec_tjet_coll;

  std::string collection_rec, collection_gen;
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

}
