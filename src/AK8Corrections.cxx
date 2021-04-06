#include "UHH2/LegacyTopTagging/include/AK8Corrections.h"

using namespace std;
using namespace uhh2;
using namespace ltt;


namespace uhh2 { namespace ltt {

/*
* AK8 CORRECTIONS SETUP
*/

AK8Corrections::AK8Corrections(const string & coll_rec, const string & coll_gen) {

  cout << "Hello World from AK8Corrections!" << endl;
  cout << "Caveat: The module AK8Corrections will not correct subjets!" << endl;

  jec_tag_2016 = "Summer16_07Aug2017";
  jec_ver_2016 = "11";

  jec_tag_2017 = "Fall17_17Nov2017";
  jec_ver_2017 = "32";

  jec_tag_2018 = "Autumn18";
  jec_ver_2018 = "19";

  jec_tag_UL17 = "Summer19UL17";
  jec_ver_UL17 = "5";

  jec_tag_UL18 = "Summer19UL18";
  jec_ver_UL18 = "5";

  collection_rec = coll_rec.empty() ? "topjets" : coll_rec;
  collection_gen = coll_gen.empty() ? "gentopjets" : coll_gen;
}

void AK8Corrections::init(Context & ctx) {

  if(init_done) {
    throw invalid_argument("AK8Corrections::init called twice!");
  }
  init_done = true;

  is_mc = ctx.get("dataset_type") == "MC";
  year = extract_year(ctx);

  string userTopJetColl = string2lowercase(collection_rec == "topjets" ? ctx.get("TopJetCollection") : collection_rec);

  string algo = "";
  if(userTopJetColl.find("ak4") != string::npos) {
    algo = "AK4";
  }
  else if(userTopJetColl.find("ak8") != string::npos) {
    algo = "AK8";
  }
  else if(userTopJetColl.find("ak8") == string::npos) {
    cout << "AK8Corrections.cxx: Cannot determine tjet cone + radius (neither AK4 nor AK8) - going to assume it is AK8 for JECs" << '\n';
    algo = "AK8";
  }

  string pus = "PFchs"; // Pileup subtraction
  if(userTopJetColl.find("puppi") != string::npos) {
    pus = "PFPuppi";
  }
  else if(userTopJetColl.find("chs") == string::npos) {
    cout << "Cannot determine pileup subtraction (neither CHS nor PUPPI) - going to assume it is CHS for JECs" << endl;
  }
  jec_tjet_coll = algo + pus;

  if(is_mc) {
    tjet_corrector_MC.reset(new YearSwitcher(ctx));
    tjet_corrector_MC->setup2016(make_shared<GenericTopJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_2016, jec_ver_2016, jec_tjet_coll), collection_rec));
    tjet_corrector_MC->setup2017(make_shared<GenericTopJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_2017, jec_ver_2017, jec_tjet_coll), collection_rec));
    tjet_corrector_MC->setup2018(make_shared<GenericTopJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_2018, jec_ver_2018, jec_tjet_coll), collection_rec));
    tjet_corrector_MC->setupUL17(make_shared<GenericTopJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_UL17, jec_ver_UL17, jec_tjet_coll), collection_rec));
    tjet_corrector_MC->setupUL18(make_shared<GenericTopJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_UL18, jec_ver_UL18, jec_tjet_coll), collection_rec));

    string jer_tag = "";
    if(year == Year::is2016v2 || year == Year::is2016v3) {
      jer_tag = "Summer16_25nsV1";
    }
    else if(year == Year::is2017v1 || year == Year::is2017v2) {
      jer_tag = "Fall17_V3";
    }
    else if(year == Year::is2018) {
      jer_tag = "Autumn18_V7";
    }
    else if(year == Year::isUL17) {
      jer_tag = "Summer19UL17_JRV2";
    }
    else if(year == Year::isUL18) {
      jer_tag = "Summer19UL18_JRV2";
    }
    else {
      throw runtime_error("Cannot find suitable jet resolution file & scale factors for this year for JetResolutionSmearer");
    }

    tjet_resolution_smearer.reset(new GenericJetResolutionSmearer(ctx, collection_rec, collection_gen, JERFiles::JERPathStringMC(jer_tag, jec_tjet_coll, "SF"), JERFiles::JERPathStringMC(jer_tag, jec_tjet_coll, "PtResolution")));
  }
  else {
    jec_switcher_16.reset(new RunSwitcher(ctx, "2016"));
    for(const auto & runItr : runPeriods2016) { // runPeriods defined in common/include/Utils.h
      jec_switcher_16->setupRun(runItr, make_shared<GenericTopJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_2016, jec_ver_2016, jec_tjet_coll, runItr), collection_rec));
    }

    jec_switcher_17.reset(new RunSwitcher(ctx, "2017"));
    for(const auto & runItr : runPeriods2017) {
      jec_switcher_17->setupRun(runItr, make_shared<GenericTopJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_2017, jec_ver_2017, jec_tjet_coll, runItr), collection_rec));
    }

    jec_switcher_18.reset(new RunSwitcher(ctx, "2018"));
    for(const auto & runItr : runPeriods2018) {
      jec_switcher_18->setupRun(runItr, make_shared<GenericTopJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_2018, jec_ver_2018, jec_tjet_coll, runItr), collection_rec));
    }

    jec_switcher_UL17.reset(new RunSwitcher(ctx, "2017"));
    for(const auto & runItr : runPeriods2017) {
      jec_switcher_UL17->setupRun(runItr, make_shared<GenericTopJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_UL17, jec_ver_UL17, jec_tjet_coll, runItr), collection_rec));
    }

    jec_switcher_UL18.reset(new RunSwitcher(ctx, "2018"));
    for(const auto & runItr : runPeriods2018) {
      jec_switcher_UL18->setupRun(runItr, make_shared<GenericTopJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_UL18, jec_ver_UL18, jec_tjet_coll, runItr), collection_rec));
    }

    tjet_corrector_data.reset(new YearSwitcher(ctx));
    tjet_corrector_data->setup2016(jec_switcher_16);
    tjet_corrector_data->setup2017(jec_switcher_17);
    tjet_corrector_data->setup2018(jec_switcher_18);
    tjet_corrector_data->setupUL17(jec_switcher_UL17);
    tjet_corrector_data->setupUL18(jec_switcher_UL18);
  }
}

bool AK8Corrections::process(uhh2::Event & event) {

  if(!init_done) {
    throw runtime_error("AK8Corrections::init not called (has to be called in AnalysisModule constructor)");
  }
  if(is_mc) {
    tjet_corrector_MC->process(event);
    tjet_resolution_smearer->process(event);
  }
  else {
    tjet_corrector_data->process(event);
  }
  return true;
}


/*
* AK8 CLEANING
*/

AK8Cleaning::AK8Cleaning(Context & ctx, const double minPt, const double maxEta, const double minDR, const string & coll_rec, const string & h_primlep_name) {

  _minPt = minPt;
  _maxEta = maxEta;
  _minDR = minDR;
  h_ak8jets_rec = ctx.get_handle<vector<TopJet>>(coll_rec.empty() ? "topjets" : coll_rec);
  h_primlep = ctx.get_handle<FlavorParticle>(h_primlep_name);
}

bool AK8Cleaning::process(uhh2::Event & event) {

  vector<TopJet> initial_topjets = event.get(h_ak8jets_rec);
  vector<TopJet> cleaned_topjets;
  const auto & primlep = event.get(h_primlep);

  for(const auto & tjet : initial_topjets) {
    if(uhh2::deltaR(tjet, primlep) < _minDR) continue;
    if(abs(tjet.v4().eta()) > _maxEta) continue;
    if(tjet.v4().pt() < _minPt) continue;
    cleaned_topjets.push_back(tjet);
  }
  if(jet_pt_sorter) sort_by_pt<TopJet>(cleaned_topjets);

  event.set(h_ak8jets_rec, cleaned_topjets);

  return true;
}

}}
