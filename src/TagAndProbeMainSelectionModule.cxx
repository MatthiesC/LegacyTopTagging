#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h"

#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/PrimaryLepton.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/common/include/JetHists.h"

#include "UHH2/LegacyTopTagging/include/Utils.h"
#include "UHH2/LegacyTopTagging/include/AndHists.h"

using namespace std;
using namespace uhh2;
using namespace ltt;


namespace uhh2 { namespace ltt {

class TagAndProbeMainSelectionModule: public AnalysisModule {
public:
    explicit TagAndProbeMainSelectionModule(Context & ctx);
    virtual bool process(Event & event) override;

private:
  bool debug;
  bool run_btag_sf;

  unique_ptr<AnalysisModule> primlep;
  unique_ptr<AnalysisModule> scale_variation;
  unique_ptr<AnalysisModule> sf_lumi;
  unique_ptr<AnalysisModule> sf_muon;
  unique_ptr<AnalysisModule> sf_btag;
  unique_ptr<AnalysisModule> sf_trigger;
  unique_ptr<Selection> slct_trigger;
  unique_ptr<Selection> slct_muon;
  unique_ptr<Selection> slct_btag;
  unique_ptr<Selection> slct_twod;

  unique_ptr<AndHists> hist_presel;
  unique_ptr<AndHists> hist_btag;
  unique_ptr<BTagMCEfficiencyHists> hist_btag_eff;
  unique_ptr<AndHists> hist_before_twod;
  unique_ptr<AndHists> hist_before_trigger;
  unique_ptr<AndHists> hist_after_trigger;
  unique_ptr<AndHists> hist_full;
};


TagAndProbeMainSelectionModule::TagAndProbeMainSelectionModule(Context & ctx) {

  debug = string2bool(ctx.get("debug"));
  const string xml_key_of_btag_eff_file = "BTagMCEffFile";
  run_btag_sf = ctx.has(xml_key_of_btag_eff_file);

  const double deltaR_leptonicHemisphere = M_PI*2./3.;

  const MuonId muonID_tag = AndId<Muon>(PtEtaCut(55., 2.4), MuonID(Muon::Selector::CutBasedIdTight));
  const BTag::algo btagALGO = BTag::DEEPJET;
  const BTag::wp btagWP = BTag::WP_MEDIUM;
  const JetId btagID = BTag(btagALGO, btagWP);

  scale_variation.reset(new MCScaleVariation(ctx));
  sf_lumi.reset(new MCLumiWeight(ctx));
  sf_muon.reset(new MuonScaleFactors(ctx));
  if(run_btag_sf) sf_btag.reset(new MCBTagScaleFactor(ctx, btagALGO, btagWP, "jets", ctx.get("SystDirection_BTag", "nominal"), "mujets", "incl", xml_key_of_btag_eff_file, "", "BTagScaleFactorFile"));
  sf_trigger.reset(new TriggerScaleFactors(ctx));

  slct_trigger.reset(new TriggerSelection("HLT_Mu50_v*"));
  slct_muon.reset(new NMuonSelection(1, 1, muonID_tag));
  slct_btag.reset(new BTagCloseToLeptonSelection(ctx, deltaR_leptonicHemisphere, btagID));
  slct_twod.reset(new TwoDSelection(ctx, 25., 0.4));

  primlep.reset(new PrimaryLepton(ctx));

  hist_presel.reset(new AndHists(ctx, "0_PreSel"));
  hist_btag_eff.reset(new BTagMCEfficiencyHists(ctx, "0_BTagMCEff", btagID));
  hist_btag.reset(new AndHists(ctx, "1_BTag"));
  hist_before_twod.reset(new AndHists(ctx, "2_BeforeTwoD"));
  hist_before_trigger.reset(new AndHists(ctx, "2_BeforeTrigger"));
  hist_after_trigger.reset(new AndHists(ctx, "2_AfterTrigger"));
  hist_full.reset(new AndHists(ctx, "3_TagAndProbeSelection"));
}


bool TagAndProbeMainSelectionModule::process(Event & event) {

  if(debug) {
    cout << endl;
    cout << "+-----------+" << endl;
    cout << "| NEW EVENT |" << endl;
    cout << "+-----------+" << endl;
  }

  if(debug) cout << "PreSelection scale factors" << endl; // else getting error for some data samples, e.g. "RunSwitcher cannot handle run number 275656 for year 2016"
  primlep->process(event);
  scale_variation->process(event);
  sf_lumi->process(event);
  sf_muon->process(event);
  // prefiring weights not yet available for UL (for updates on this, see https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1ECALPrefiringWeightRecipe)
  hist_presel->fill(event);

  if(debug) cout << "2018 HEM15/16 issue selection" << endl;
  //TODO

  if(debug) cout << "Booleans for further selections" << endl;
  const bool passes_trigger = slct_trigger->passes(event);
  const bool passes_muon = slct_muon->passes(event);
  const bool passes_btag = slct_btag->passes(event);
  const bool passes_twod = slct_twod->passes(event);

  if(passes_trigger) sf_trigger->process(event);
  if(passes_trigger && passes_muon && passes_twod) {
    hist_btag_eff->fill(event);
  }
  if(!passes_btag) return false;
  if(!(passes_trigger || passes_muon || passes_twod)) return false;
  if(run_btag_sf) sf_btag->process(event);
  hist_btag->fill(event);
  if(passes_muon && passes_trigger) {
    hist_before_twod->fill(event);
  }
  if(passes_twod) {
    hist_before_trigger->fill(event);
    if(passes_trigger) hist_after_trigger->fill(event);
  }
  if(!(passes_trigger && passes_muon && passes_twod)) return false;
  hist_full->fill(event);



  if(debug) cout << "End of TagAndProbeMainSelectionModule" << endl;
  return true;
}


UHH2_REGISTER_ANALYSIS_MODULE(TagAndProbeMainSelectionModule)

}}
