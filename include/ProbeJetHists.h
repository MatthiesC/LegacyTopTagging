#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/JetIds.h"

#include "UHH2/LegacyTopTagging/include/Constants.h"
#include "UHH2/LegacyTopTagging/include/Utils.h"


namespace uhh2 { namespace ltt {

class AK8ProbeJetHists: public uhh2::Hists {
public:
  AK8ProbeJetHists(uhh2::Context & ctx, const std::string & dirname, const MergeScenario & _msc);
  virtual void fill(const uhh2::Event & event) override;

protected:
  enum class PassCategory {
    Pass,
    Fail,
  };

  enum class JetCategory {
    All,
    Mass,
    BTag,
    MassAndBTag,
  };

  std::map<PtBin,
    std::map<JetCategory,
      std::map<WorkingPoint,
        std::map<PassCategory,
          std::vector<TH1F*>>>>> hists_map;

private:
  void fill_probe(const std::vector<TH1F*> & hists);

  uhh2::Event::Handle<FlavorParticle> h_primlep;
  uhh2::Event::Handle<TopJet> h_probejet;
  uhh2::Event::Handle<MergeScenario> h_merge_scenario;
  MergeScenario msc;
  const std::vector<PtBin> pt_bins = kPtBinsAK8;
  const WorkingPointMap wps = kWorkingPointsAK8;
  const JetId SubjetBTagID = BTag(BTag::DEEPCSV, BTag::WP_LOOSE);
  double w;
  FlavorParticle primlep;
  TopJet probejet;

  std::map<PassCategory, std::string> kPassCategoryAsString = {
    {PassCategory::Pass, "Pass"},
    {PassCategory::Fail, "Fail"},
  };

  std::map<JetCategory, std::string> kJetCategoryAsString = {
    {JetCategory::All, "All"},
    {JetCategory::Mass, "Mass"},
    {JetCategory::BTag, "BTag"},
    {JetCategory::MassAndBTag, "MassAndBTag"},
  };
};

class HOTVRProbeJetHists: public uhh2::Hists {
public:
  HOTVRProbeJetHists(uhh2::Context & ctx, const std::string & dirname, const MergeScenario & _msc);
  virtual void fill(const uhh2::Event & event) override;

protected:
  enum class PassCategory {
    Pass,
    Fail,
  };

  enum class JetCategory {
    All,
    Mass,
    HOTVRCuts,
    HOTVRCutsAndMass,
  };

  std::map<PtBin,
    std::map<JetCategory,
      std::map<WorkingPoint,
        std::map<PassCategory,
          std::vector<TH1F*>>>>> hists_map;

private:
  void fill_probe(const std::vector<TH1F*> & hists);

  uhh2::Event::Handle<FlavorParticle> h_primlep;
  uhh2::Event::Handle<TopJet> h_probejet;
  uhh2::Event::Handle<MergeScenario> h_merge_scenario;
  MergeScenario msc;
  const std::vector<PtBin> pt_bins = kPtBinsHOTVR;
  const WorkingPointMap wps = kWorkingPointsHOTVR;
  const TopJetId HOTVRTopTagID = ltt::HOTVRTopTag(0., std::numeric_limits<double>::infinity()); // standard HOTVR t-tag without mass cut and without tau32 cut
  double w;
  FlavorParticle primlep;
  TopJet probejet;

  std::map<PassCategory, std::string> kPassCategoryAsString = {
    {PassCategory::Pass, "Pass"},
    {PassCategory::Fail, "Fail"},
  };

  std::map<JetCategory, std::string> kJetCategoryAsString = {
    {JetCategory::All, "All"},
    {JetCategory::Mass, "Mass"},
    {JetCategory::HOTVRCuts, "HOTVRCuts"},
    {JetCategory::HOTVRCutsAndMass, "HOTVRCutsAndMass"},
  };
};

class ProbeJetHistsRunner: public uhh2::Hists {
public:
  ProbeJetHistsRunner(uhh2::Context & ctx, const std::string & dirname);
  virtual void fill(const uhh2::Event & event) override;

protected:
  std::vector<uhh2::Hists*> hists_vector;
};

}}
