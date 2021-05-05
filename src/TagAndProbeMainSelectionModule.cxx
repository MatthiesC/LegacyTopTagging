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
#include "UHH2/LegacyTopTagging/include/TopJetCorrections.h"
#include "UHH2/LegacyTopTagging/include/Constants.h"
#include "UHH2/LegacyTopTagging/include/ProbeJetHists.h"

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
  Process proc;
  ProbeJetAlgo pjalgo;
  MergeScenario msc;

  unique_ptr<AnalysisModule> primlep;
  unique_ptr<AnalysisModule> scale_variation;
  unique_ptr<AnalysisModule> sf_lumi;
  unique_ptr<AnalysisModule> sf_pileup;
  unique_ptr<AnalysisModule> sf_muon;
  unique_ptr<AnalysisModule> sf_btag;
  unique_ptr<AnalysisModule> sf_trigger;
  unique_ptr<Selection> slct_trigger;
  unique_ptr<Selection> slct_muon;
  unique_ptr<Selection> slct_btag;
  unique_ptr<Selection> slct_twod;

  unique_ptr<TopJetCorrections> corrections_hotvr;
  unique_ptr<TopJetCorrections> corrections_ak8;
  unique_ptr<TopJetCleaning> cleaner_hotvr;
  unique_ptr<TopJetCleaning> cleaner_ak8;
  unique_ptr<Selection> slct_1hotvr;
  unique_ptr<Selection> slct_1ak8;
  unique_ptr<AnalysisModule> decay_channel_and_hadronic_top;
  unique_ptr<AnalysisModule> probejet_hotvr;
  unique_ptr<AnalysisModule> probejet_ak8;
  unique_ptr<AnalysisModule> merge_scenarios_hotvr;
  unique_ptr<AnalysisModule> merge_scenarios_ak8;
  unique_ptr<AnalysisModule> main_output;

  unique_ptr<AndHists> hist_presel;
  unique_ptr<AndHists> hist_btag;
  unique_ptr<BTagMCEfficiencyHists> hist_btag_eff;
  unique_ptr<AndHists> hist_before_twod;
  unique_ptr<AndHists> hist_before_trigger;
  unique_ptr<AndHists> hist_after_trigger;
  unique_ptr<AndHists> hist_full_before_corrections;
  unique_ptr<AndHists> hist_full_after_corrections;
  unique_ptr<AndHists> hist_full_after_cleaner;

  unique_ptr<Hists> probejethists;
};


TagAndProbeMainSelectionModule::TagAndProbeMainSelectionModule(Context & ctx) {

  debug = string2bool(ctx.get("debug"));
  const string dv = ctx.get("dataset_version");

  ctx.undeclare_all_event_output();

  const double deltaR_leptonicHemisphere = M_PI*2./3.;
  const double hotvr_pt_min = 200.;
  const double hotvr_eta_max = 2.5;
  const double ak8_pt_min = 300.;
  const double ak8_eta_max = 2.4;

  const MuonId muonID_tag = AndId<Muon>(PtEtaCut(55., 2.4), MuonID(Muon::Selector::CutBasedIdTight));
  const BTag::algo btagALGO = BTag::DEEPJET;
  const BTag::wp btagWP = BTag::WP_MEDIUM;
  const JetId btagID = BTag(btagALGO, btagWP);

  scale_variation.reset(new MCScaleVariation(ctx));
  sf_lumi.reset(new MCLumiWeight(ctx));
  sf_pileup.reset(new MCPileupReweight(ctx, ctx.get("SystDirection_Pileup", "nominal")));
  sf_muon.reset(new MuonScaleFactors(ctx));
  const string xml_key_of_btag_eff_file = "BTagMCEffFile";
  run_btag_sf = ctx.has(xml_key_of_btag_eff_file);
  if(run_btag_sf) sf_btag.reset(new MCBTagScaleFactor(ctx, btagALGO, btagWP, "jets", ctx.get("SystDirection_BTag", "nominal"), "mujets", "incl", xml_key_of_btag_eff_file, "", "BTagScaleFactorFile"));
  sf_trigger.reset(new TriggerScaleFactors(ctx));

  slct_trigger.reset(new TriggerSelection("HLT_Mu50_v*"));
  slct_muon.reset(new NMuonSelection(1, 1, muonID_tag));
  slct_btag.reset(new BTagCloseToLeptonSelection(ctx, deltaR_leptonicHemisphere, btagID));
  slct_twod.reset(new TwoDSelection(ctx, 25., 0.4));

  primlep.reset(new PrimaryLepton(ctx));

  corrections_hotvr.reset(new TopJetCorrections());
  corrections_hotvr->switch_topjet_corrections(false);
  corrections_hotvr->switch_subjet_corrections(true);
  corrections_hotvr->switch_rebuilding_topjets_from_subjets(true);
  corrections_hotvr->init(ctx);
  corrections_ak8.reset(new TopJetCorrections(ctx.get("AK8Collection_rec"), ctx.get("AK8Collection_gen")));
  corrections_ak8->init(ctx);

  cleaner_hotvr.reset(new TopJetCleaning(ctx, hotvr_pt_min, hotvr_eta_max, deltaR_leptonicHemisphere));
  cleaner_ak8.reset(new TopJetCleaning(ctx, ak8_pt_min, ak8_eta_max, deltaR_leptonicHemisphere, ctx.get("AK8Collection_rec")));

  slct_1hotvr.reset(new NTopJetSelection(1, -1));
  slct_1ak8.reset(new NTopJetSelection(1, -1, boost::none, ctx.get_handle<vector<TopJet>>(ctx.get("AK8Collection_rec"))));

  decay_channel_and_hadronic_top.reset(new DecayChannelAndHadronicTopHandleSetter(ctx));
  probejet_hotvr.reset(new ProbeJetHandleSetter(ctx, ProbeJetAlgo::isHOTVR));
  probejet_ak8.reset(new ProbeJetHandleSetter(ctx, ProbeJetAlgo::isAK8, ctx.get("AK8Collection_rec")));
  merge_scenarios_hotvr.reset(new MergeScenarioHandleSetter(ctx, ProbeJetAlgo::isHOTVR));
  merge_scenarios_ak8.reset(new MergeScenarioHandleSetter(ctx, ProbeJetAlgo::isAK8));
  main_output.reset(new MainOutputSetter(ctx));

  hist_presel.reset(new AndHists(ctx, "0_PreSel"));
  hist_btag_eff.reset(new BTagMCEfficiencyHists(ctx, "0_BTagMCEff", btagID));
  hist_btag.reset(new AndHists(ctx, "1_BTag"));
  hist_before_twod.reset(new AndHists(ctx, "2_BeforeTwoD"));
  hist_before_trigger.reset(new AndHists(ctx, "2_BeforeTrigger"));
  hist_after_trigger.reset(new AndHists(ctx, "2_AfterTrigger"));
  hist_full_before_corrections.reset(new AndHists(ctx, "3_TagAndProbeSelection_BeforeCorr", true));
  hist_full_after_corrections.reset(new AndHists(ctx, "3_TagAndProbeSelection_AfterCorr", true));
  hist_full_after_cleaner.reset(new AndHists(ctx, "3_TagAndProbeSelection_AfterCleaning", true));

  probejethists.reset(new ProbeJetHistsRunner(ctx, "ProbeJetHists"));
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
  sf_pileup->process(event);
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
  hist_full_before_corrections->fill(event);
  corrections_hotvr->process(event);
  corrections_ak8->process(event);
  hist_full_after_corrections->fill(event);
  cleaner_hotvr->process(event);
  cleaner_ak8->process(event);
  hist_full_after_cleaner->fill(event);

  if(debug) cout << "Identify probe jets" << endl;
  const bool has_hotvr_jet = slct_1hotvr->passes(event);
  const bool has_ak8_jet = slct_1ak8->passes(event);

  if(!(has_hotvr_jet || has_ak8_jet)) return false;
  decay_channel_and_hadronic_top->process(event);

  if(has_hotvr_jet) probejet_hotvr->process(event);
  if(has_ak8_jet) probejet_ak8->process(event);
  merge_scenarios_hotvr->process(event); // needs to be outside of the previous if statement!
  merge_scenarios_ak8->process(event); // needs to be outside of the previous if statement!

  main_output->process(event);

  probejethists->fill(event);

  if(debug) cout << "End of TagAndProbeMainSelectionModule" << endl;
  return false;
}


UHH2_REGISTER_ANALYSIS_MODULE(TagAndProbeMainSelectionModule)

}}
