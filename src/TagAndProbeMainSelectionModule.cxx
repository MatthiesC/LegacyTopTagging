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

  unique_ptr<AnalysisModule> sf_toppt;
  unique_ptr<AnalysisModule> sf_vjets;
  unique_ptr<AnalysisModule> primlep;
  unique_ptr<AnalysisModule> scale_variation;
  unique_ptr<AnalysisModule> ps_variation;
  unique_ptr<AnalysisModule> sf_lumi;
  unique_ptr<AnalysisModule> sf_prefire;
  unique_ptr<AnalysisModule> sf_pileup;
  unique_ptr<AnalysisModule> sf_muon;
  unique_ptr<AnalysisModule> sf_btag;
  unique_ptr<AnalysisModule> sf_trigger;
  unique_ptr<Selection> slct_ak4;
  unique_ptr<HEM2018Selection> slct_hem2018;
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
  unique_ptr<AndHists> hist_ak4;
  unique_ptr<AndHists> hist_hem2018;
  unique_ptr<AndHists> hist_btag;
  unique_ptr<BTagMCEfficiencyHists> hist_btag_eff;
  unique_ptr<AndHists> hist_before_twod;
  unique_ptr<AndHists> hist_before_trigger_incl_low_muon_pt;
  unique_ptr<AndHists> hist_before_trigger_excl_low_muon_pt;
  unique_ptr<AndHists> hist_after_trigger_incl_low_muon_pt;
  unique_ptr<AndHists> hist_after_trigger_excl_low_muon_pt;
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

  sf_toppt.reset(new ltt::TopPtReweighting(ctx));
  sf_vjets.reset(new ltt::VJetsReweighting(ctx));
  scale_variation.reset(new MCScaleVariation(ctx));
  ps_variation.reset(new PartonShowerVariation(ctx));
  sf_lumi.reset(new MCLumiWeight(ctx));
  sf_prefire.reset(new PrefiringWeights(ctx));
  sf_pileup.reset(new MCPileupReweight(ctx, ctx.get("SystDirection_Pileup", "nominal")));
  sf_muon.reset(new MuonScaleFactors(ctx));
  const string xml_key_of_btag_eff_file = "BTagMCEffFile";
  run_btag_sf = ctx.has(xml_key_of_btag_eff_file);
  if(run_btag_sf) sf_btag.reset(new MCBTagScaleFactor(ctx, btagALGO, btagWP, "jets", ctx.get("SystDirection_BTag", "nominal"), "mujets", "incl", xml_key_of_btag_eff_file, "", "BTagScaleFactorFile"));
  sf_trigger.reset(new TriggerScaleFactors(ctx));

  slct_ak4.reset(new NJetSelection(1, -1));
  slct_hem2018.reset(new HEM2018Selection(ctx));
  slct_trigger.reset(new TriggerSelection("HLT_Mu50_v*"));
  // slct_trigger.reset(new LttTriggerSelection(ctx));
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
  hist_ak4.reset(new AndHists(ctx, "1_AK4"));
  hist_hem2018.reset(new AndHists(ctx, "1_HEM2018"));
  hist_btag_eff.reset(new BTagMCEfficiencyHists(ctx, "2_BTagMCEff", btagID));
  hist_btag.reset(new AndHists(ctx, "3_BTag"));
  hist_before_twod.reset(new AndHists(ctx, "4_BeforeTwoD"));
  hist_before_trigger_incl_low_muon_pt.reset(new AndHists(ctx, "4_BeforeTriggerInclLowMuonPt"));
  hist_before_trigger_excl_low_muon_pt.reset(new AndHists(ctx, "4_BeforeTriggerExclLowMuonPt"));
  hist_after_trigger_incl_low_muon_pt.reset(new AndHists(ctx, "4_AfterTriggerInclLowMuonPt"));
  hist_after_trigger_excl_low_muon_pt.reset(new AndHists(ctx, "5_AfterTriggerExclLowMuonPt", true)); // include AK8 and HOTVR plots from here on
  hist_full_after_corrections.reset(new AndHists(ctx, "5_AfterTopJetCorrections", true));
  hist_full_after_cleaner.reset(new AndHists(ctx, "5_AfterTopJetCleaning", true));

  probejethists.reset(new ProbeJetHistsRunner(ctx, "ProbeJetHists"));
}


bool TagAndProbeMainSelectionModule::process(Event & event) {

  if(debug) {
    cout << endl;
    cout << "+-----------+" << endl;
    cout << "| NEW EVENT |" << endl;
    cout << "+-----------+" << endl;
  }

  if(debug) cout << "PreSelection scale factors" << endl;
  primlep->process(event);
  scale_variation->process(event);
  ps_variation->process(event);
  sf_lumi->process(event);
  sf_prefire->process(event);
  sf_pileup->process(event);
  sf_muon->process(event);
  sf_toppt->process(event);
  sf_vjets->process(event);
  hist_presel->fill(event);

   // It is required to have >= 1 AK4 jet to calculate values like dR(jet, lepton) for the TwoDSelection.
   // Later on, the BTagCloseToLeptonSelection will require >= 1 AK4 jet anyway.
  if(debug) cout << "AK4 selection" << endl;
  if(!slct_ak4->passes(event)) return false;
  hist_ak4->fill(event);

  if(debug) cout << "2018 HEM15/16 issue selection" << endl;
  if(slct_hem2018->passes(event)) {
    if(event.isRealData) return false;
    else event.weight *= (1. - slct_hem2018->GetAffectedLumiFraction());
  }
  hist_hem2018->fill(event);

  if(debug) cout << "Booleans for further selections" << endl;
  const bool passes_trigger = slct_trigger->passes(event);
  const bool passes_muon = slct_muon->passes(event);
  const bool passes_btag = slct_btag->passes(event);
  const bool passes_twod = slct_twod->passes(event);

  if(passes_trigger) {
    if(debug) cout << "Apply trigger scale factor" << endl;
    sf_trigger->process(event);
  }
  else {
    if(debug) cout << "Skipped trigger scale factor application" << endl;
  }
  if(passes_trigger && passes_muon && passes_twod) {
    if(debug) cout << "Fill b-tagging efficiency histograms" << endl;
    hist_btag_eff->fill(event);
  }
  else {
    if(debug) cout << "Skipped filling b-tagging efficiency histograms" << endl;
  }
  if(!passes_btag) return false;
  if(!(passes_trigger || passes_muon || passes_twod)) return false;

  if(debug) cout << "Apply b-tagging scale factors (if possible)" << endl;
  if(run_btag_sf) sf_btag->process(event);
  hist_btag->fill(event);

  if(debug) cout << "Fill some control histograms regarding trigger efficiency and 2D cut" << endl;
  if(passes_muon && passes_trigger) {
    hist_before_twod->fill(event);
  }
  if(passes_twod) {
    hist_before_trigger_incl_low_muon_pt->fill(event);
    if(passes_muon) hist_before_trigger_excl_low_muon_pt->fill(event);
    if(passes_trigger) hist_after_trigger_incl_low_muon_pt->fill(event);
  }
  if(!(passes_trigger && passes_muon && passes_twod)) return false;
  hist_after_trigger_excl_low_muon_pt->fill(event);

  if(debug) cout << "Apply corrections to top jet collections" << endl;
  corrections_hotvr->process(event);
  corrections_ak8->process(event);
  hist_full_after_corrections->fill(event);
  cleaner_hotvr->process(event);
  cleaner_ak8->process(event);
  hist_full_after_cleaner->fill(event);

  if(debug) cout << "Check whether event has probe jets or not" << endl;
  const bool has_hotvr_jet = slct_1hotvr->passes(event);
  const bool has_ak8_jet = slct_1ak8->passes(event);
  if(!(has_hotvr_jet || has_ak8_jet)) return false;

  if(debug) cout << "Find out decay channel and identify generated hadronic top quark (only valid for top-MC)" << endl;
  decay_channel_and_hadronic_top->process(event);

  if(debug) cout << "Set probe jet handles and find out merge scenario" << endl;
  if(has_hotvr_jet) probejet_hotvr->process(event);
  if(has_ak8_jet) probejet_ak8->process(event);

  // Following modules need to be outside of the previous if statements! MergeScenario for event w/o probe jet will be "isBackground"
  merge_scenarios_hotvr->process(event);
  merge_scenarios_ak8->process(event);

  if(debug) cout << "Add probe jet properties to output tree" << endl;
  main_output->process(event);

  if(debug) cout << "Fill probe jet histograms" << endl;
  probejethists->fill(event);

  if(debug) cout << "End of TagAndProbeMainSelectionModule" << endl;
  return false;
}


UHH2_REGISTER_ANALYSIS_MODULE(TagAndProbeMainSelectionModule)

}}
