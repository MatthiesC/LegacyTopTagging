#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h"

#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/Utils.h"

#include "UHH2/LegacyTopTagging/include/AK8Hists.h"
#include "UHH2/LegacyTopTagging/include/Utils.h"

using namespace std;
using namespace uhh2;
using namespace ltt;


namespace uhh2 { namespace ltt {

class RawAK8JetsModule: public AnalysisModule {
public:
    explicit RawAK8JetsModule(Context & ctx);
    virtual bool process(Event & event) override;

private:
  bool debug;
  string dataset_version;

  unique_ptr<AnalysisModule> sf_lumi;
  unique_ptr<AnalysisModule> sf_pileup;
  unique_ptr<AnalysisModule> ak8_cleaner;
  unique_ptr<Selection> slct_1ak8;

  unique_ptr<Hists> hists_ak8;
  unique_ptr<Hists> hists_ak8_matched;
};


RawAK8JetsModule::RawAK8JetsModule(Context & ctx) {

  debug = string2bool(ctx.get("debug"));
  dataset_version = ctx.get("dataset_version");

  ctx.undeclare_all_event_output();

  sf_lumi.reset(new MCLumiWeight(ctx));
  sf_pileup.reset(new MCPileupReweight(ctx));
  TopJetId ak8_id = PtEtaCut(300, 2.4);
  ak8_cleaner.reset(new TopJetCleaner(ctx, ak8_id));
  slct_1ak8.reset(new NTopJetSelection(1, -1));

  hists_ak8.reset(new ltt::AK8Hists(ctx, "AK8Hists"));
  hists_ak8_matched.reset(new ltt::AK8Hists(ctx, "AK8Hists_matched"));
}


bool RawAK8JetsModule::process(Event & event) {

  if(debug) {
    cout << endl;
    cout << "+-----------+" << endl;
    cout << "| NEW EVENT |" << endl;
    cout << "+-----------+" << endl;
  }

  sort_by_pt(*event.topjets);

  sf_lumi->process(event);
  sf_pileup->process(event);

  hists_ak8->fill(event);

  ak8_cleaner->process(event);
  if(!slct_1ak8->passes(event)) return false;

  GenParticle top, antitop;
  for(GenParticle gp : *event.genparticles) {
    if(gp.pdgId() == 6) top = gp;
    else if(gp.pdgId() == -6) antitop = gp;
  }
  const vector<TopJet> & topjets = *event.topjets;
  const TopJet *tnearestjet = nextTopJet(top, topjets);
  const TopJet *antitnearestjet = nextTopJet(antitop, topjets);
  const bool the_two_jets_are_the_same = (tnearestjet == antitnearestjet);
  if(the_two_jets_are_the_same) return false;
  vector<TopJet> matched_topjets;
  if(deltaR(tnearestjet->v4(), top.v4()) < 0.6) matched_topjets.push_back(*tnearestjet);
  if(deltaR(antitnearestjet->v4(), antitop.v4()) < 0.6) matched_topjets.push_back(*antitnearestjet);
  swap(matched_topjets, *event.topjets);

  hists_ak8_matched->fill(event);

  if(debug) cout << "End of RawAK8JetsModule" << endl;
  return false;
}


UHH2_REGISTER_ANALYSIS_MODULE(RawAK8JetsModule)

}}
