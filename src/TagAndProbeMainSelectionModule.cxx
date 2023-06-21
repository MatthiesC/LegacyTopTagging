#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/JetHists.h"
#include "UHH2/common/include/PrimaryLepton.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/PSWeights.h"

#include "UHH2/LegacyTopTagging/include/AK8Hists.h"
#include "UHH2/LegacyTopTagging/include/AndHists.h"
#include "UHH2/LegacyTopTagging/include/Constants.h"
#include "UHH2/LegacyTopTagging/include/HOTVRHists.h"
#include "UHH2/LegacyTopTagging/include/LeptonScaleFactors.h"
#include "UHH2/LegacyTopTagging/include/JetMETCorrections.h"
#include "UHH2/LegacyTopTagging/include/SingleTopGen_tWch.h"
#include "UHH2/LegacyTopTagging/include/TopJetCorrections.h"
#include "UHH2/LegacyTopTagging/include/TriggerSelection.h"
#include "UHH2/LegacyTopTagging/include/Utils.h"


using namespace std;
using namespace uhh2;
using namespace uhh2::ltt;


namespace uhh2 { namespace ltt {

class TagAndProbeMainSelectionModule: public AnalysisModule {
public:
  explicit TagAndProbeMainSelectionModule(Context & ctx);
  virtual bool process(Event & event) override;

private:
  const bool debug;
  unsigned long long i_event = 0;
  const Channel fChannel;
  const Year fYear;
  const string fDatasetVersion;
  const string fDatasetVersion_without_year_suffix;

  /*
  Kinematic variables:
  */

  const double muon_eta_max = 2.4;
  const double muon_lowpt_pt_min = 55.0; // 30.0 // only use highpt muons for tag-and-probe SF measurement
  const double muon_highpt_pt_min = 55.0;

  const double electron_eta_max = 2.4;
  const double electron_lowpt_pt_min = 55.0; // 30.0
  const double electron_highpt_pt_min = 120.0;

  const double met_min = 50.0;

  const double ptw_min = 120.0;

  const double hotvr_pt_min = 200.;
  const double hotvr_eta_max = 2.5;
  const double hotvr_dr_lep_min = kDeltaRLeptonicHemisphere;

  const double ak8_pt_min = 200.;
  const double ak8_eta_max = 2.5;
  const double ak8_dr_lep_min = kDeltaRLeptonicHemisphere;

  const double ak4_pt_min = 30.;
  // const double ak4_eta_max = 5.0; // keep forward jets; will define different handles for central (= b-taggable) and forward jets
  const double ak4_eta_max = 2.5;
  const double ak4_dr_lep_min = 0.4;
  const double ak4_ptrel_max = 30.; // 25.

  /*
  All modules etc. in chronological order as used in event loop:
  */

  bool is_tW;
  Event::Handle<bool> fHandle_bool_reco_sel;
  unique_ptr<AnalysisModule> prod_SingleTopGen_tWch;

  bool is_syst = false;

  const Muon::Selector muonIDselector_lowpt = Muon::Selector::CutBasedIdTight;
  const MuonId muonID_lowpt = AndId<Muon>(MuonID(muonIDselector_lowpt), PtEtaCut(muon_lowpt_pt_min, muon_eta_max));
  const Muon::Selector muonIDselector_highpt = Muon::Selector::CutBasedIdGlobalHighPt;
  const MuonId muonID_highpt = AndId<Muon>(MuonID(muonIDselector_highpt), PtEtaCut(muon_highpt_pt_min, muon_eta_max));

  // can be increased from wp90 to wp80 later...
  const Electron::tag electronIDtag_lowpt = Electron::tag::mvaEleID_Fall17_noIso_V2_wp90;
  const ElectronId electronID_lowpt = AndId<Electron>(ElectronTagID(electronIDtag_lowpt), PtEtaCut(electron_lowpt_pt_min, electron_eta_max));
  const Electron::tag electronIDtag_highpt = Electron::tag::mvaEleID_Fall17_noIso_V2_wp90;
  const ElectronId electronID_highpt = AndId<Electron>(ElectronTagID(electronIDtag_highpt), PtEtaCut(electron_highpt_pt_min, electron_eta_max));

  unique_ptr<Selection> slct_muon_lowpt;
  unique_ptr<Selection> slct_muon_highpt;
  unique_ptr<Selection> slct_elec_lowpt;
  unique_ptr<Selection> slct_elec_highpt;

  unique_ptr<AnalysisModule> sf_muon_id_highpt;
  unique_ptr<AnalysisModule> sf_muon_id_lowpt;
  unique_ptr<AnalysisModule> sf_muon_id_dummy;
  unique_ptr<AnalysisModule> sf_elec_id_highpt;
  unique_ptr<AnalysisModule> sf_elec_id_lowpt;
  unique_ptr<AnalysisModule> sf_elec_id_dummy;
  unique_ptr<AnalysisModule> sf_elec_reco;
  unique_ptr<AnalysisModule> sf_elec_reco_dummy;

  unique_ptr<AnalysisModule> primlep;
  unique_ptr<AnalysisModule> scale_variation;
  unique_ptr<AnalysisModule> ps_variation;
  unique_ptr<AnalysisModule> sf_lumi;
  unique_ptr<AnalysisModule> sf_pileup;
  unique_ptr<AnalysisModule> sf_prefire;
  unique_ptr<AnalysisModule> weight_trickery;

  unique_ptr<ltt::JetMETCorrections> jetmet_corrections_puppi;
  unique_ptr<AnalysisModule> cleaner_ak4puppi;
  unique_ptr<ltt::JetMETCorrections> jetmet_corrections_chs;

  unique_ptr<ltt::TopJetCorrections> corrections_hotvr;
  unique_ptr<AnalysisModule> cleaner_hotvr;
  unique_ptr<ltt::TopJetCorrections> corrections_ak8;
  unique_ptr<AnalysisModule> cleaner_ak8;

  unique_ptr<AnalysisModule> object_pt_sorter;
  const BTag::algo btagAlgo = BTag::algo::DEEPJET;
  const BTag::wp btagWP = BTag::wp::WP_MEDIUM;
  unique_ptr<AnalysisModule> puppichs_matching;

  unique_ptr<Selection> slct_met;
  unique_ptr<Selection> slct_metfilter;

  // unique_ptr<OrSelection> slct_1largejet;
  unique_ptr<Selection> slct_1hotvr;
  unique_ptr<Selection> slct_1ak8;

  unique_ptr<Selection> slct_1ak4jet;

  unique_ptr<Selection> slct_ptw;

  unique_ptr<ltt::HEM2018Selection> slct_hem2018;

  unique_ptr<AnalysisModule> sf_toppt;

  unique_ptr<AnalysisModule> sf_vjets;

  unique_ptr<Selection> slct_trigger_highpt;
  unique_ptr<Selection> slct_trigger_lowpt;

  unique_ptr<AnalysisModule> sf_muon_trigger_highpt;
  unique_ptr<AnalysisModule> sf_muon_trigger_lowpt;
  unique_ptr<AnalysisModule> sf_muon_trigger_dummy;

  unique_ptr<Selection> slct_twod;

  unique_ptr<Selection> slct_btag;
  const map<Band, string> xml_key_of_btag_eff_file = {
    { Band::MAIN, "BTagMCEffFile" },
    { Band::QCD, "BTagMCEffFile_QCDSideband" },
  };
  map<Band, bool> run_btag_sf;
  map<Band, unique_ptr<Hists>> hist_btag_eff;
  map<Band, unique_ptr<AnalysisModule>> sf_btagging;

  unique_ptr<ltt::AndHists> hist_before2d;
  map<Band, unique_ptr<ltt::AndHists>> hist_presel;

  unique_ptr<AnalysisModule> decay_channel_and_hadronic_top;
  unique_ptr<AnalysisModule> probejet_hotvr;
  unique_ptr<AnalysisModule> probejet_ak8;
  unique_ptr<AnalysisModule> merge_scenarios_hotvr;
  unique_ptr<AnalysisModule> merge_scenarios_ak8;
  unique_ptr<AnalysisModule> main_output;

  Event::Handle<float> fHandle_weight;
  Event::Handle<int> fHandle_band;
  Event::Handle<int> fHandle_year;
  Event::Handle<string> fHandle_dataset;
};


TagAndProbeMainSelectionModule::TagAndProbeMainSelectionModule(Context & ctx):
  debug(string2bool(ctx.get("debug"))),
  fChannel(extract_channel(ctx)),
  fYear(extract_year(ctx)),
  fDatasetVersion(ctx.get("dataset_version")),
  fDatasetVersion_without_year_suffix(extract_dataset(ctx))
{
  ctx.undeclare_all_event_output(); // throw away all output trees (jet collections etc.) which are not needed in further steps of the analysis
  // empty_output_tree = string2bool(ctx.get("EmptyOutputTree")); // handy to not have output trees for systematics files, reduces root file size
  // is_QCDsideband = string2bool(ctx.get("QCD_sideband"));

  // unsigned int i_hist(0);

  cout << "This is the input dataset: \"" << fDatasetVersion_without_year_suffix.c_str() << "\"\n";

  fHandle_bool_reco_sel = ctx.get_handle<bool>("btw_bool_reco_sel"); // kHandleName_bool_reco_sel // I really should have merged HighPtSingleTop and LegacyTopTagging into one repo...
  is_tW = fDatasetVersion.find("ST_tW") == 0;
  prod_SingleTopGen_tWch.reset(new ltt::SingleTopGen_tWchProducer(ctx, kHandleName_SingleTopGen_tWch));

  if(string2bool(ctx.get("extra_syst"))) is_syst = true;
  if(ctx.get("jecsmear_direction") != "nominal") is_syst = true;
  if(ctx.get("jersmear_direction") != "nominal") is_syst = true;

  slct_muon_lowpt.reset(new NMuonSelection(1, 1, muonID_lowpt));
  slct_muon_highpt.reset(new NMuonSelection(1, 1, muonID_highpt));
  slct_elec_lowpt.reset(new NElectronSelection(1, 1, electronID_lowpt));
  slct_elec_highpt.reset(new NElectronSelection(1, 1, electronID_highpt));

  sf_muon_id_highpt.reset(new ltt::MuonIdScaleFactors(ctx, muonIDselector_highpt));
  sf_muon_id_lowpt.reset(new ltt::MuonIdScaleFactors(ctx, muonIDselector_lowpt));
  sf_muon_id_dummy.reset(new ltt::MuonIdScaleFactors(ctx, boost::none, boost::none, boost::none, boost::none, true));
  sf_elec_id_highpt.reset(new ltt::ElectronIdScaleFactors(ctx, electronIDtag_highpt));
  sf_elec_id_lowpt.reset(new ltt::ElectronIdScaleFactors(ctx, electronIDtag_lowpt));
  sf_elec_id_dummy.reset(new ltt::ElectronIdScaleFactors(ctx, boost::none, boost::none, boost::none, boost::none, true));
  sf_elec_reco.reset(new ltt::ElectronRecoScaleFactors(ctx));
  sf_elec_reco_dummy.reset(new ltt::ElectronRecoScaleFactors(ctx, boost::none, boost::none, boost::none, boost::none, true));

  primlep.reset(new PrimaryLepton(ctx));
  scale_variation.reset(new MCScaleVariation(ctx));
  // ps_variation.reset(new ltt::PartonShowerVariation(ctx));
  bool has_ps_weights = true;
  if(fDatasetVersion.find("QCD") == 0) has_ps_weights = false;
  ps_variation.reset(new PSWeights(ctx, true, !has_ps_weights));
  sf_lumi.reset(new MCLumiWeight(ctx));
  sf_pileup.reset(new MCPileupReweight(ctx, ctx.get("SystDirection_Pileup", "nominal")));
  sf_prefire.reset(new ltt::PrefiringWeights(ctx));
  weight_trickery.reset(new ltt::WeightTrickery(ctx, kHandleName_SingleTopGen_tWch, false));

  // const JetId jetID = AndId<Jet>(PtEtaCut(ak4_pt_min, ak4_eta_max), ltt::NoLeptonInJet("all", ak4_dr_lep_min)); // JetPFID is already part of JetMETCorrections
  const JetId jetID = PtEtaCut(ak4_pt_min, ak4_eta_max); // JetPFID is already part of JetMETCorrections

  const string met_name = ctx.get("METName");
  const bool puppi_met = (met_name == kCollectionName_METPUPPI);
  jetmet_corrections_puppi.reset(new ltt::JetMETCorrections(boost::none, boost::none, puppi_met ? boost::none : (boost::optional<std::string>)kCollectionName_METPUPPI));
  jetmet_corrections_puppi->init(ctx);
  cleaner_ak4puppi.reset(new JetCleaner(ctx, jetID));
  jetmet_corrections_chs.reset(new ltt::JetMETCorrections(kCollectionName_AK4CHS, boost::none, puppi_met ? (boost::optional<std::string>)kCollectionName_METCHS : boost::none));
  jetmet_corrections_chs->init(ctx);

  const TopJetId hotvrID = AndId<TopJet>(JetPFID(JetPFID::WP_TIGHT_PUPPI), PtEtaCut(hotvr_pt_min, hotvr_eta_max), ltt::NoLeptonInJet("all", hotvr_dr_lep_min));
  const TopJetId ak8ID = AndId<TopJet>(JetPFID(JetPFID::WP_TIGHT_PUPPI), PtEtaCut(ak8_pt_min, ak8_eta_max), ltt::NoLeptonInJet("all", ak8_dr_lep_min));

  corrections_hotvr.reset(new ltt::TopJetCorrections());
  corrections_hotvr->switch_topjet_corrections(false);
  corrections_hotvr->switch_subjet_corrections(true);
  corrections_hotvr->switch_rebuilding_topjets_from_subjets(true);
  corrections_hotvr->init(ctx);
  cleaner_hotvr.reset(new TopJetCleaner(ctx, hotvrID));
  corrections_ak8.reset(new ltt::TopJetCorrections(kCollectionName_AK8_rec, kCollectionName_AK8_gen));
  corrections_ak8->init(ctx);
  cleaner_ak8.reset(new TopJetCleaner(ctx, ak8ID, kCollectionName_AK8_rec));

  object_pt_sorter.reset(new ltt::ObjectPtSorter(ctx));
  puppichs_matching.reset(new ltt::MatchPuppiToCHSAndSetBTagHandles(ctx, btagAlgo, btagWP));

  slct_met.reset(new ltt::METSelection(ctx, met_min));
  slct_metfilter.reset(new ltt::METFilterSelection(ctx));

  // slct_1largejet.reset(new OrSelection());
  // slct_1largejet->add(make_shared<NTopJetSelection>(1, -1));
  // slct_1largejet->add(make_shared<NTopJetSelection>(1, -1, boost::none, ctx.get_handle<vector<TopJet>>(kCollectionName_AK8_rec)));
  slct_1hotvr.reset(new NTopJetSelection(1, -1));
  slct_1ak8.reset(new NTopJetSelection(1, -1, boost::none, ctx.get_handle<vector<TopJet>>(kCollectionName_AK8_rec)));

  slct_1ak4jet.reset(new NJetSelection(1, -1, boost::none, ctx.get_handle<vector<Jet>>(kHandleName_pairedPUPPIjets)));

  slct_ptw.reset(new ltt::PTWSelection(ctx, ptw_min));

  slct_hem2018.reset(new ltt::HEM2018Selection(ctx));

  sf_toppt.reset(new ltt::TopPtReweighting(ctx, string2bool(ctx.get("apply_TopPtReweighting"))));

  sf_vjets.reset(new ltt::VJetsReweighting(ctx));

  slct_trigger_highpt.reset(new ltt::MyTriggerSelection(ctx, false, true));
  slct_trigger_lowpt.reset(new ltt::MyTriggerSelection(ctx, true, true));

  sf_muon_trigger_highpt.reset(new ltt::MuonTriggerScaleFactors(ctx, true, false)); // avoid check since we do not request official muon isolation
  sf_muon_trigger_lowpt.reset(new ltt::MuonTriggerScaleFactors(ctx, false, false)); // --"--
  sf_muon_trigger_dummy.reset(new ltt::MuonTriggerScaleFactors(ctx, boost::none, boost::none, boost::none, boost::none, boost::none, true));

  slct_twod.reset(new ltt::TwoDSelection(ctx, ak4_ptrel_max, ak4_dr_lep_min, true));

  const JetId btagID = BTag(btagAlgo, btagWP);
  slct_btag.reset(new NJetSelection(1, -1, boost::none, ctx.get_handle<vector<Jet>>(kHandleName_bJets_hemi)));
  for(const auto & band : kRelevantBands) {
    run_btag_sf[band] = ctx.has(xml_key_of_btag_eff_file.at(band));
    hist_btag_eff[band].reset(new BTagMCEfficiencyHists(ctx, "BTagMCEff_"+kBands.at(band).name, btagID, kHandleName_pairedCHSjets_hemi));
    if(run_btag_sf.at(band)) sf_btagging[band].reset(new MCBTagScaleFactor(ctx, btagAlgo, btagWP, kHandleName_pairedCHSjets_hemi, "mujets", "incl", xml_key_of_btag_eff_file.at(band)));
  }

  hist_before2d.reset(new ltt::AndHists(ctx, "Before2D", true, true));
  for(const Band & band : kRelevantBands) {
    hist_presel[band].reset(new ltt::AndHists(ctx, "Presel_"+kBands.at(band).name, true, true));
  }

  decay_channel_and_hadronic_top.reset(new ltt::DecayChannelAndHadronicTopHandleSetter(ctx));
  probejet_hotvr.reset(new ltt::ProbeJetHandleSetter(ctx, ProbeJetAlgo::isHOTVR));
  probejet_ak8.reset(new ltt::ProbeJetHandleSetter(ctx, ProbeJetAlgo::isAK8, kCollectionName_AK8_rec));
  merge_scenarios_hotvr.reset(new ltt::MergeScenarioHandleSetter(ctx, ProbeJetAlgo::isHOTVR, kHandleName_SingleTopGen_tWch));
  merge_scenarios_ak8.reset(new ltt::MergeScenarioHandleSetter(ctx, ProbeJetAlgo::isAK8, kHandleName_SingleTopGen_tWch));
  main_output.reset(new ltt::MainOutputSetter(ctx));

  fHandle_weight = ctx.declare_event_output<float>("weight");
  fHandle_band = ctx.declare_event_output<int>("band");
  fHandle_year = ctx.declare_event_output<int>("year");
  fHandle_dataset = ctx.declare_event_output<string>("dataset");
}


//------------//
// EVENT LOOP //
//------------//

bool TagAndProbeMainSelectionModule::process(Event & event) {

  if(debug) {
    cout << endl;
    cout << "+-----------+" << endl;
    cout << "| NEW EVENT |" << endl;
    cout << "+-----------+" << endl;
    cout << "i_event = " << to_string(i_event++) << endl;
    cout << endl;
  }

  if(debug) cout << "Initial stuff after preselection" << endl;
  if(event.is_valid(fHandle_bool_reco_sel) && event.get(fHandle_bool_reco_sel) == false) return false; // reject (tW) events which only pass gen level selections during preselection from HighPtSingleTop
  if(fChannel == Channel::isMuo) {
    if(event.muons->size() != 1) return false;
  }
  else if(fChannel == Channel::isEle) {
    if(event.electrons->size() != 1) return false;
  }
  if(is_tW) {
    prod_SingleTopGen_tWch->process(event); // needed for WeightTrickery and setting correct merge scenario for tW
  }

  bool lowpt(false);
  if(fChannel == Channel::isMuo) {
    const Muon *muon = &event.muons->at(0);
    if(muon->pt() < muon_highpt_pt_min) {
      if(!slct_muon_lowpt->passes(event)) return false;
      sf_muon_id_lowpt->process(event);
      lowpt = true;
    }
    else {
      if(!slct_muon_highpt->passes(event)) return false;
      sf_muon_id_highpt->process(event);
      lowpt = false;
    }
    sf_elec_id_dummy->process(event);
    sf_elec_reco_dummy->process(event);
  }
  else if(fChannel == Channel::isEle) {
    const Electron *electron = &event.electrons->at(0);
    const float abseta_sc = fabs(electron->supercluster_eta());
    if(abseta_sc > 1.4442 && abseta_sc < 1.566) return false; // gap between ECAL barrel and endcap
    if(electron->pt() < electron_highpt_pt_min) {
      if(!slct_elec_lowpt->passes(event)) return false;
      sf_elec_id_lowpt->process(event);
      lowpt = true;
    }
    else {
      if(!slct_elec_highpt->passes(event)) return false;
      sf_elec_id_highpt->process(event);
      lowpt = false;
    }
    sf_elec_reco->process(event);
    sf_muon_id_dummy->process(event);
  }
  primlep->process(event);
  scale_variation->process(event);
  ps_variation->process(event);
  sf_lumi->process(event);
  sf_pileup->process(event);
  sf_prefire->process(event);
  weight_trickery->process(event);

  if(debug) cout << "JetMET corrections, pt sorting, and PUPPI-CHS matching" << endl;
  jetmet_corrections_puppi->process(event);
  cleaner_ak4puppi->process(event);
  jetmet_corrections_chs->process(event);

  corrections_hotvr->process(event); // needs to come already here because of subsequent object pt sorter
  cleaner_hotvr->process(event);
  corrections_ak8->process(event);
  cleaner_ak8->process(event);

  object_pt_sorter->process(event); // needs to come after jet corrections but before PUPPI-CHS matching
  puppichs_matching->process(event);

  if(debug) cout << "MET selection and MET filters" << endl;
  if(!slct_met->passes(event)) return false;
  if(!slct_metfilter->passes(event)) return false;

  // if(debug) cout << "At least one large-R jet (HOTVR or AK8)" << endl;
  // if(!slct_1largejet->passes(event)) return false;
  if(debug) cout << "Check whether event has probe jets or not" << endl;
  const bool has_hotvr_jet = slct_1hotvr->passes(event);
  const bool has_ak8_jet = slct_1ak8->passes(event);
  if(!(has_hotvr_jet || has_ak8_jet)) return false;

  if(debug) cout << "Select at least one AK4 jet" << endl;
  if(!slct_1ak4jet->passes(event)) return false; // Require at least one AK4 jet for computational reasons (dR(lepton, jet) etc.); this rejects only \mathcal{O}(0.01\%) of events in real data (tested in pre-UL 2017 muo, RunB)

  if(debug) cout << "Leptonic W boson pT selection" << endl;
  if(!slct_ptw->passes(event)) return false;

  if(debug) cout << "2018 HEM15/16 issue selection" << endl;
  if(slct_hem2018->passes(event)) {
    if(event.isRealData) return false;
    else event.weight *= (1. - slct_hem2018->GetAffectedLumiFraction());
  }

  if(debug) cout << "Apply top-pt reweighting for ttbar events" << endl;
  sf_toppt->process(event);

  if(debug) cout << "Apply (N)NLO QCD/EWK corrections to V+jets samples" << endl;
  sf_vjets->process(event);

  if(debug) cout << "Trigger selection" << endl;
  const bool passes_trigger = lowpt ? slct_trigger_lowpt->passes(event) : slct_trigger_highpt->passes(event);
  if(!passes_trigger) return false;
  if(fChannel == Channel::isMuo) {
    if(lowpt) sf_muon_trigger_lowpt->process(event);
    else sf_muon_trigger_highpt->process(event);
    // sf_elec_trigger_dummy->process(event);
  }
  else if(fChannel == Channel::isEle) {
    // sf_elec_trigger->process(event); // need to differentiate between 2017 Run B and Run C-F
    sf_muon_trigger_dummy->process(event);
  }

  if(debug) cout << "Booleans for further selections" << endl;
  const bool passes_twod = slct_twod->passes(event);
  Band band;
  if(passes_twod) band = Band::MAIN;
  else {
    if(is_syst) return false;
    band = Band::QCD;
  }

  hist_btag_eff[band]->fill(event);

  const bool passes_btag = slct_btag->passes(event);
  if(!passes_btag) return false;
  if(run_btag_sf[band]) sf_btagging[band]->process(event);

  hist_before2d->fill(event);
  hist_presel[band]->fill(event);

  if(debug) cout << "Find out decay channel and identify generated hadronic top quark (only valid for top-MC)" << endl;
  decay_channel_and_hadronic_top->process(event);

  if(debug) cout << "Set probe jet handles and find out merge scenario" << endl;
  if(has_hotvr_jet) probejet_hotvr->process(event);
  if(has_ak8_jet) probejet_ak8->process(event);

  // Following modules need to be outside of the previous if statements! MergeScenario for event w/o probe jet will be "isBackground"
  merge_scenarios_hotvr->process(event);
  merge_scenarios_ak8->process(event);

  if(debug) cout << "Write main output to AnalysisTree" << endl;
  main_output->process(event);

  if(debug) cout << "End of MainSelectionModule. Event passed" << endl;
  event.set(fHandle_weight, event.weight);
  event.set(fHandle_band, kBands.at(band).index);
  event.set(fHandle_year, kYears.at(fYear).index);
  // event.set(fHandle_dataset, fDatasetVersion_without_year_suffix);
  event.set(fHandle_dataset, fDatasetVersion); // with year suffix

  return true;
}


UHH2_REGISTER_ANALYSIS_MODULE(TagAndProbeMainSelectionModule)

}}
