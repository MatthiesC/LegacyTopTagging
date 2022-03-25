#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/JetIds.h"

#include "UHH2/LegacyTopTagging/include/JetMETCorrections.h"
#include "UHH2/LegacyTopTagging/include/METXYCorrection.h"
#include "UHH2/LegacyTopTagging/include/Constants.h"

using namespace std;
using namespace uhh2;
using namespace ltt;


namespace uhh2 { namespace ltt {

//____________________________________________________________________________________________________
// Type-I MET correction
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETRun2Corrections
void propagate_JEC_to_MET(const Event & event, const Event::Handle<vector<Jet>> & h_jets) {
  // we start from raw MET
  LorentzVector metv4 = event.met->uncorr_v4();
  for(const auto & jet : event.get(h_jets)) {
    const bool to_be_corrected = jet.v4().Pt() > 15. && (jet.neutralEmEnergyFraction() + jet.chargedEmEnergyFraction()) < 0.9;
    if(to_be_corrected) {
      // slimmed MET is corrected by L1FastJet
      const auto factor_raw = jet.JEC_factor_raw();
      const auto L1factor_raw = jet.JEC_L1factor_raw();

      const LorentzVector L1corr = (L1factor_raw * factor_raw) * jet.v4(); // L1 corrected jets
      const LorentzVector L123corr = jet.v4(); // L123 corrected jets (L23 in case of PUPPI)
      metv4 -= L123corr;

      // slimmed MET is corrected by L1FastJet, for PUPPI: L1factor_raw = 1 --> L1corr = raw-jet pT.
      metv4 += L1corr;
    }
  }
  event.met->set_pt(metv4.Pt());
  event.met->set_phi(metv4.Phi());
}

//____________________________________________________________________________________________________
JetMETCorrections::JetMETCorrections(const string & coll_rec, const string & coll_gen) {

  cout << "Hello World from JetMETCorrections!" << endl;

  jec_tag_2016 = k_jec_tag_2016;
  jec_ver_2016 = k_jec_ver_2016;
  jer_tag_2016 = k_jer_tag_2016;

  jec_tag_2017 = k_jec_tag_2017;
  jec_ver_2017 = k_jec_ver_2017;
  jer_tag_2017 = k_jer_tag_2017;

  jec_tag_2018 = k_jec_tag_2018;
  jec_ver_2018 = k_jec_ver_2018;
  jer_tag_2018 = k_jer_tag_2018;

  jec_tag_UL16preVFP = k_jec_tag_UL16preVFP;
  jec_ver_UL16preVFP = k_jec_ver_UL16preVFP;
  jer_tag_UL16preVFP = k_jer_tag_UL16preVFP;

  jec_tag_UL16postVFP = k_jec_tag_UL16postVFP;
  jec_ver_UL16postVFP = k_jec_ver_UL16postVFP;
  jer_tag_UL16postVFP = k_jer_tag_UL16postVFP;

  jec_tag_UL17 = k_jec_tag_UL17;
  jec_ver_UL17 = k_jec_ver_UL17;
  jer_tag_UL17 = k_jer_tag_UL17;

  jec_tag_UL18 = k_jec_tag_UL18;
  jec_ver_UL18 = k_jec_ver_UL18;
  jer_tag_UL18 = k_jer_tag_UL18;

  if(!coll_rec.empty()) {
    collection_rec = coll_rec;
    use_additional_branch_for_rec = true;
  }

  if(!coll_gen.empty()) {
    collection_gen = coll_gen;
    use_additional_branch_for_gen = true;
  }
}

void JetMETCorrections::init(Context & ctx) {

  if(init_done) {
    throw runtime_error("JetMETCorrections::init() called twice!");
  }
  init_done = true;

  if(do_jlc) {
    cout << "JetMETCorrections will clean jets in jet collection '" << collection_rec << "' from leptons (ensure to clean leptons properly beforehand!)" << endl;
  }
  if(do_jec) {
    cout << "JetMETCorrections will correct jets in jet collection '" << collection_rec << "'" << endl;
  }
  if(do_met_type1_correction) {
    cout << "JetMETCorrections will propagate JES corrections for jet collection '" << collection_rec << "' to MET" << endl;
  }
  if(do_met_xy_correction) {
    cout << "JetMETCorrections will correct MET XY" << endl;
  }

  is_mc = ctx.get("dataset_type") == "MC";
  year = extract_year(ctx);

  h_jets = ctx.get_handle<vector<Jet>>(collection_rec);
  // h_genjets = ctx.get_handle<vector<GenJet>>(collection_gen);

  string userJetColl = string2lowercase(use_additional_branch_for_rec ? collection_rec : ctx.get("JetCollection"));

  string algo = "";
  if(userJetColl.find("ak4") != string::npos) {
    algo = "AK4";
  }
  else if(userJetColl.find("ak8") != string::npos) {
    algo = "AK8";
  }
  else {
    cout << "JetMETCorrections::init(): Cannot determine jet cone + radius (neither AK4 nor AK8) - going to assume it is AK4 for identifying JEC files" << endl;
    algo = "AK4";
  }

  JetPFID::wp jetpfID_wp = JetPFID::WP_TIGHT_CHS;
  string pus = "PFchs"; // Pileup subtraction
  if(userJetColl.find("puppi") != string::npos) {
    pus = "PFPuppi";
    jetpfID_wp = JetPFID::WP_TIGHT_PUPPI;
  }
  else if(userJetColl.find("chs") == string::npos) {
    cout << "JetMETCorrections::init(): Cannot determine pile-up subtraction (neither CHS nor PUPPI) - going to assume it is CHS for identifying JEC files and jet PF ID" << endl;
  }
  string jec_jet_coll = algo + pus;

  clnr_jetpfid.reset(new JetCleaner(ctx, JetPFID(jetpfID_wp), collection_rec));

  if(is_mc) {
    jlc_MC.reset(new YearSwitcher(ctx));
    jlc_MC->setup2016(make_shared<JetLeptonCleaner_by_KEYmatching>(ctx, JERFiles::JECFilesMC(jec_tag_2016, jec_ver_2016, jec_jet_coll), collection_rec));
    jlc_MC->setup2017(make_shared<JetLeptonCleaner_by_KEYmatching>(ctx, JERFiles::JECFilesMC(jec_tag_2017, jec_ver_2017, jec_jet_coll), collection_rec));
    jlc_MC->setup2018(make_shared<JetLeptonCleaner_by_KEYmatching>(ctx, JERFiles::JECFilesMC(jec_tag_2018, jec_ver_2018, jec_jet_coll), collection_rec));
    jlc_MC->setupUL16preVFP(make_shared<JetLeptonCleaner_by_KEYmatching>(ctx, JERFiles::JECFilesMC(jec_tag_UL16preVFP, jec_ver_UL16preVFP, jec_jet_coll), collection_rec));
    jlc_MC->setupUL16postVFP(make_shared<JetLeptonCleaner_by_KEYmatching>(ctx, JERFiles::JECFilesMC(jec_tag_UL16postVFP, jec_ver_UL16postVFP, jec_jet_coll), collection_rec));
    jlc_MC->setupUL17(make_shared<JetLeptonCleaner_by_KEYmatching>(ctx, JERFiles::JECFilesMC(jec_tag_UL17, jec_ver_UL17, jec_jet_coll), collection_rec));
    jlc_MC->setupUL18(make_shared<JetLeptonCleaner_by_KEYmatching>(ctx, JERFiles::JECFilesMC(jec_tag_UL18, jec_ver_UL18, jec_jet_coll), collection_rec));

    jet_corrector_MC.reset(new YearSwitcher(ctx));
    jet_corrector_MC->setup2016(make_shared<GenericJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_2016, jec_ver_2016, jec_jet_coll), collection_rec));
    jet_corrector_MC->setup2017(make_shared<GenericJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_2017, jec_ver_2017, jec_jet_coll), collection_rec));
    jet_corrector_MC->setup2018(make_shared<GenericJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_2018, jec_ver_2018, jec_jet_coll), collection_rec));
    jet_corrector_MC->setupUL16preVFP(make_shared<GenericJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_UL16preVFP, jec_ver_UL16preVFP, jec_jet_coll), collection_rec));
    jet_corrector_MC->setupUL16postVFP(make_shared<GenericJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_UL16postVFP, jec_ver_UL16postVFP, jec_jet_coll), collection_rec));
    jet_corrector_MC->setupUL17(make_shared<GenericJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_UL17, jec_ver_UL17, jec_jet_coll), collection_rec));
    jet_corrector_MC->setupUL18(make_shared<GenericJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_UL18, jec_ver_UL18, jec_jet_coll), collection_rec));

    string jer_tag = "";
    if(year == Year::is2016v2 || year == Year::is2016v3) {
      jer_tag = jer_tag_2016;
    }
    else if(year == Year::is2017v1 || year == Year::is2017v2) {
      jer_tag = jer_tag_2017;
    }
    else if(year == Year::is2018) {
      jer_tag = jer_tag_2018;
    }
    else if(year == Year::isUL16preVFP) {
      jer_tag = jer_tag_UL16preVFP;
    }
    else if(year == Year::isUL16postVFP) {
      jer_tag = jer_tag_UL16postVFP;
    }
    else if(year == Year::isUL17) {
      jer_tag = jer_tag_UL17;
    }
    else if(year == Year::isUL18) {
      jer_tag = jer_tag_UL18;
    }
    else {
      throw runtime_error("Cannot find suitable jet resolution file & scale factors for this year for JetResolutionSmearer");
    }

    jet_resolution_smearer.reset(new GenericJetResolutionSmearer(ctx, collection_rec, collection_gen, JERFiles::JERPathStringMC(jer_tag, jec_jet_coll, "SF"), JERFiles::JERPathStringMC(jer_tag, jec_jet_coll, "PtResolution")));
  }
  else {
    jlc_switcher_16.reset(new RunSwitcher(ctx, "2016"));
    jec_switcher_16.reset(new RunSwitcher(ctx, "2016"));
    for(const auto & runItr : runPeriods2016) { // runPeriods defined in common/include/Utils.h
      jlc_switcher_16->setupRun(runItr, make_shared<JetLeptonCleaner_by_KEYmatching>(ctx, JERFiles::JECFilesDATA(jec_tag_2016, jec_ver_2016, jec_jet_coll, runItr), collection_rec));
      jec_switcher_16->setupRun(runItr, make_shared<GenericJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_2016, jec_ver_2016, jec_jet_coll, runItr), collection_rec));
    }

    jlc_switcher_17.reset(new RunSwitcher(ctx, "2017"));
    jec_switcher_17.reset(new RunSwitcher(ctx, "2017"));
    for(const auto & runItr : runPeriods2017) {
      jlc_switcher_17->setupRun(runItr, make_shared<JetLeptonCleaner_by_KEYmatching>(ctx, JERFiles::JECFilesDATA(jec_tag_2017, jec_ver_2017, jec_jet_coll, runItr), collection_rec));
      jec_switcher_17->setupRun(runItr, make_shared<GenericJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_2017, jec_ver_2017, jec_jet_coll, runItr), collection_rec));
    }

    jlc_switcher_18.reset(new RunSwitcher(ctx, "2018"));
    jec_switcher_18.reset(new RunSwitcher(ctx, "2018"));
    for(const auto & runItr : runPeriods2018) {
      jlc_switcher_18->setupRun(runItr, make_shared<JetLeptonCleaner_by_KEYmatching>(ctx, JERFiles::JECFilesDATA(jec_tag_2018, jec_ver_2018, jec_jet_coll, runItr), collection_rec));
      jec_switcher_18->setupRun(runItr, make_shared<GenericJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_2018, jec_ver_2018, jec_jet_coll, runItr), collection_rec));
    }

    jlc_switcher_UL16preVFP.reset(new RunSwitcher(ctx, "2016"));
    jec_switcher_UL16preVFP.reset(new RunSwitcher(ctx, "2016"));
    for(const auto & runItr : runPeriodsUL16preVFP) {
      jlc_switcher_UL16preVFP->setupRun(runItr, make_shared<JetLeptonCleaner_by_KEYmatching>(ctx, JERFiles::JECFilesDATA(jec_tag_UL16preVFP, jec_ver_UL16preVFP, jec_jet_coll, runItr), collection_rec));
      jec_switcher_UL16preVFP->setupRun(runItr, make_shared<GenericJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_UL16preVFP, jec_ver_UL16preVFP, jec_jet_coll, runItr), collection_rec));
    }

    jlc_switcher_UL16postVFP.reset(new RunSwitcher(ctx, "2016"));
    jec_switcher_UL16postVFP.reset(new RunSwitcher(ctx, "2016"));
    for(const auto & runItr : runPeriodsUL16postVFP) {
      jlc_switcher_UL16postVFP->setupRun(runItr, make_shared<JetLeptonCleaner_by_KEYmatching>(ctx, JERFiles::JECFilesDATA(jec_tag_UL16postVFP, jec_ver_UL16postVFP, jec_jet_coll, runItr), collection_rec));
      jec_switcher_UL16postVFP->setupRun(runItr, make_shared<GenericJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_UL16postVFP, jec_ver_UL16postVFP, jec_jet_coll, runItr), collection_rec));
    }

    jlc_switcher_UL17.reset(new RunSwitcher(ctx, "2017"));
    jec_switcher_UL17.reset(new RunSwitcher(ctx, "2017"));
    for(const auto & runItr : runPeriods2017) {
      jlc_switcher_UL17->setupRun(runItr, make_shared<JetLeptonCleaner_by_KEYmatching>(ctx, JERFiles::JECFilesDATA(jec_tag_UL17, jec_ver_UL17, jec_jet_coll, runItr), collection_rec));
      jec_switcher_UL17->setupRun(runItr, make_shared<GenericJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_UL17, jec_ver_UL17, jec_jet_coll, runItr), collection_rec));
    }

    jlc_switcher_UL18.reset(new RunSwitcher(ctx, "2018"));
    jec_switcher_UL18.reset(new RunSwitcher(ctx, "2018"));
    for(const auto & runItr : runPeriods2018) {
      jlc_switcher_UL18->setupRun(runItr, make_shared<JetLeptonCleaner_by_KEYmatching>(ctx, JERFiles::JECFilesDATA(jec_tag_UL18, jec_ver_UL18, jec_jet_coll, runItr), collection_rec));
      jec_switcher_UL18->setupRun(runItr, make_shared<GenericJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_UL18, jec_ver_UL18, jec_jet_coll, runItr), collection_rec));
    }

    jlc_data.reset(new YearSwitcher(ctx));
    jlc_data->setup2016(jlc_switcher_16);
    jlc_data->setup2017(jlc_switcher_17);
    jlc_data->setup2018(jlc_switcher_18);
    jlc_data->setupUL16preVFP(jlc_switcher_UL16preVFP);
    jlc_data->setupUL16postVFP(jlc_switcher_UL16postVFP);
    jlc_data->setupUL17(jlc_switcher_UL17);
    jlc_data->setupUL18(jlc_switcher_UL18);

    jet_corrector_data.reset(new YearSwitcher(ctx));
    jet_corrector_data->setup2016(jec_switcher_16);
    jet_corrector_data->setup2017(jec_switcher_17);
    jet_corrector_data->setup2018(jec_switcher_18);
    jet_corrector_data->setupUL16preVFP(jec_switcher_UL16preVFP);
    jet_corrector_data->setupUL16postVFP(jec_switcher_UL16postVFP);
    jet_corrector_data->setupUL17(jec_switcher_UL17);
    jet_corrector_data->setupUL18(jec_switcher_UL18);
  }

  met_xy_correction.reset(new ltt::METXYCorrector(ctx));
}

bool JetMETCorrections::process(Event & event) {

  if(!init_done) {
    throw runtime_error("JetMETCorrections::init() not called!");
  }

  clnr_jetpfid->process(event);
  if(is_mc) {
    if(do_jlc) jlc_MC->process(event);
    if(do_jec) {
      jet_corrector_MC->process(event);
      jet_resolution_smearer->process(event);
    }
  }
  else {
    if(do_jlc) jlc_data->process(event);
    if(do_jec) jet_corrector_data->process(event);
  }
  if(do_met_type1_correction) propagate_JEC_to_MET(event, h_jets); // needs to be done AFTER JLC and JLC needs to be done AFTER cleaning leptons
  if(do_met_xy_correction) met_xy_correction->process(event);

  return true;
}

}}
