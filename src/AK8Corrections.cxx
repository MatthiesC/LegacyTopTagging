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

  if(!coll_rec.empty()) {
    collection_rec = coll_rec;
    use_additional_branch_for_rec = true;
  }

  if(!coll_gen.empty()) {
    collection_gen = coll_gen;
    use_additional_branch_for_gen = true;
  }
}


void AK8Corrections::init(Context & ctx) {

  if(init_done) {
    throw runtime_error("AK8Corrections::init() called twice!");
  }
  init_done = true;

  if(correct_topjets) {
    cout << "AK8Corrections will correct superjets in t jet collection '" << collection_rec << "'" << endl;
  }
  else {
    cout << "AK8Corrections will NOT correct superjets in t jet collection '" << collection_rec << "'" << endl;
  }
  if(correct_subjets) {
    cout << "AK8Corrections will correct subjets in t jet collection '" << collection_rec << "'" << endl;
  }
  else {
    cout << "AK8Corrections will NOT correct subjets in t jet collection '" << collection_rec << "'" << endl;
  }
  if(!(correct_topjets || correct_subjets)) {
    throw invalid_argument("AK8Corrections::init(): You are initializing AK8Corrections but you specified to not correct t jets nor their subjets. This seems to be unintended. Please check!");
  }

  is_mc = ctx.get("dataset_type") == "MC";
  year = extract_year(ctx);

  h_topjets = ctx.get_handle<vector<TopJet>>(collection_rec);
  h_gentopjets = ctx.get_handle<vector<GenTopJet>>(collection_gen);

  string gensubjets_handle_name = (string)"ak8corrections_gensubjets_handle_for_" + collection_gen;
  h_gensubjets = ctx.get_handle<vector<GenJet>>(gensubjets_handle_name);
  string subjets_handle_name = (string)"ak8corrections_subjets_handle_for_" + collection_rec;
  h_subjets = ctx.get_handle<vector<Jet>>(subjets_handle_name);
  string subjets_map_handle_name = (string)"ak8corrections_subjets_map_handle_for_" + collection_rec;
  h_subjets_map = ctx.get_handle<vector<pair<int, int>>>(subjets_map_handle_name);

  string userTopJetColl = string2lowercase(use_additional_branch_for_rec ? collection_rec : ctx.get("TopJetCollection"));

  string algo = "";
  if(userTopJetColl.find("ak4") != string::npos) {
    algo = "AK4";
    if(correct_subjets) {
      throw invalid_argument("AK8Corrections::init(): You specified to correct subjets but your t jet collection is an AK4 collection. AK4 jets do not have subjets. Please check!");
    }
  }
  else if(userTopJetColl.find("ak8") != string::npos) {
    algo = "AK8";
  }
  else {
    cout << "AK8Corrections::init(): Cannot determine t jet cone + radius (neither AK4 nor AK8) - going to assume it is AK8" << '\n';
    algo = "AK8";
  }

  string pus = "PFchs"; // Pileup subtraction
  if(userTopJetColl.find("puppi") != string::npos) {
    pus = "PFPuppi";
  }
  else if(userTopJetColl.find("chs") == string::npos) {
    cout << "AK8Corrections::init(): Cannot determine pile-up subtraction (neither CHS nor PUPPI) - going to assume it is CHS" << endl;
  }
  string jec_tjet_coll = algo + pus;
  string jec_subjet_coll = (string)"AK4" + pus; // going to assume that subjets will always be corrected like AK4 jets

  if(is_mc) {
    tjet_corrector_MC.reset(new YearSwitcher(ctx));
    tjet_corrector_MC->setup2016(make_shared<GenericTopJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_2016, jec_ver_2016, jec_tjet_coll), collection_rec));
    tjet_corrector_MC->setup2017(make_shared<GenericTopJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_2017, jec_ver_2017, jec_tjet_coll), collection_rec));
    tjet_corrector_MC->setup2018(make_shared<GenericTopJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_2018, jec_ver_2018, jec_tjet_coll), collection_rec));
    tjet_corrector_MC->setupUL17(make_shared<GenericTopJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_UL17, jec_ver_UL17, jec_tjet_coll), collection_rec));
    tjet_corrector_MC->setupUL18(make_shared<GenericTopJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_UL18, jec_ver_UL18, jec_tjet_coll), collection_rec));

    subjet_corrector_MC.reset(new YearSwitcher(ctx));
    subjet_corrector_MC->setup2016(make_shared<GenericSubJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_2016, jec_ver_2016, jec_subjet_coll), collection_rec));
    subjet_corrector_MC->setup2017(make_shared<GenericSubJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_2017, jec_ver_2017, jec_subjet_coll), collection_rec));
    subjet_corrector_MC->setup2018(make_shared<GenericSubJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_2018, jec_ver_2018, jec_subjet_coll), collection_rec));
    subjet_corrector_MC->setupUL17(make_shared<GenericSubJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_UL17, jec_ver_UL17, jec_subjet_coll), collection_rec));
    subjet_corrector_MC->setupUL18(make_shared<GenericSubJetCorrector>(ctx, JERFiles::JECFilesMC(jec_tag_UL18, jec_ver_UL18, jec_subjet_coll), collection_rec));

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
    subjet_resolution_smearer.reset(new GenericJetResolutionSmearer(ctx, subjets_handle_name, gensubjets_handle_name, JERFiles::JERPathStringMC(jer_tag, jec_subjet_coll, "SF"), JERFiles::JERPathStringMC(jer_tag, jec_subjet_coll, "PtResolution")));
  }
  else {
    jec_switcher_16_topjets.reset(new RunSwitcher(ctx, "2016"));
    jec_switcher_16_subjets.reset(new RunSwitcher(ctx, "2016"));
    for(const auto & runItr : runPeriods2016) { // runPeriods defined in common/include/Utils.h
      jec_switcher_16_topjets->setupRun(runItr, make_shared<GenericTopJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_2016, jec_ver_2016, jec_tjet_coll, runItr), collection_rec));
      jec_switcher_16_subjets->setupRun(runItr, make_shared<GenericSubJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_2016, jec_ver_2016, jec_subjet_coll, runItr), collection_rec));
    }

    jec_switcher_17_topjets.reset(new RunSwitcher(ctx, "2017"));
    jec_switcher_17_subjets.reset(new RunSwitcher(ctx, "2017"));
    for(const auto & runItr : runPeriods2017) {
      jec_switcher_17_topjets->setupRun(runItr, make_shared<GenericTopJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_2017, jec_ver_2017, jec_tjet_coll, runItr), collection_rec));
      jec_switcher_17_subjets->setupRun(runItr, make_shared<GenericSubJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_2017, jec_ver_2017, jec_subjet_coll, runItr), collection_rec));
    }

    jec_switcher_18_topjets.reset(new RunSwitcher(ctx, "2018"));
    jec_switcher_18_subjets.reset(new RunSwitcher(ctx, "2018"));
    for(const auto & runItr : runPeriods2018) {
      jec_switcher_18_topjets->setupRun(runItr, make_shared<GenericTopJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_2018, jec_ver_2018, jec_tjet_coll, runItr), collection_rec));
      jec_switcher_18_subjets->setupRun(runItr, make_shared<GenericSubJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_2018, jec_ver_2018, jec_subjet_coll, runItr), collection_rec));
    }

    jec_switcher_UL17_topjets.reset(new RunSwitcher(ctx, "2017"));
    jec_switcher_UL17_subjets.reset(new RunSwitcher(ctx, "2017"));
    for(const auto & runItr : runPeriods2017) {
      jec_switcher_UL17_topjets->setupRun(runItr, make_shared<GenericTopJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_UL17, jec_ver_UL17, jec_tjet_coll, runItr), collection_rec));
      jec_switcher_UL17_subjets->setupRun(runItr, make_shared<GenericSubJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_UL17, jec_ver_UL17, jec_subjet_coll, runItr), collection_rec));
    }

    jec_switcher_UL18_topjets.reset(new RunSwitcher(ctx, "2018"));
    jec_switcher_UL18_subjets.reset(new RunSwitcher(ctx, "2018"));
    for(const auto & runItr : runPeriods2018) {
      jec_switcher_UL18_topjets->setupRun(runItr, make_shared<GenericTopJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_UL18, jec_ver_UL18, jec_tjet_coll, runItr), collection_rec));
      jec_switcher_UL18_subjets->setupRun(runItr, make_shared<GenericSubJetCorrector>(ctx, JERFiles::JECFilesDATA(jec_tag_UL18, jec_ver_UL18, jec_subjet_coll, runItr), collection_rec));
    }

    tjet_corrector_data.reset(new YearSwitcher(ctx));
    tjet_corrector_data->setup2016(jec_switcher_16_topjets);
    tjet_corrector_data->setup2017(jec_switcher_17_topjets);
    tjet_corrector_data->setup2018(jec_switcher_18_topjets);
    tjet_corrector_data->setupUL17(jec_switcher_UL17_topjets);
    tjet_corrector_data->setupUL18(jec_switcher_UL18_topjets);

    subjet_corrector_data.reset(new YearSwitcher(ctx));
    subjet_corrector_data->setup2016(jec_switcher_16_subjets);
    subjet_corrector_data->setup2017(jec_switcher_17_subjets);
    subjet_corrector_data->setup2018(jec_switcher_18_subjets);
    subjet_corrector_data->setupUL17(jec_switcher_UL17_subjets);
    subjet_corrector_data->setupUL18(jec_switcher_UL18_subjets);
  }
}


void AK8Corrections::set_subjet_handles(Event & event) {

  if(event.isRealData) {
    throw runtime_error("AK8Corrections::set_subjet_handles() may only be called on MC events");
  }

  vector<GenJet> gensubjets;
  vector<GenTopJet> *gentopjets = use_additional_branch_for_gen ? &event.get(h_gentopjets) : event.gentopjets;
  for(const GenTopJet & genjet : *gentopjets) {
    for(const GenJet & gensubjet : genjet.subjets()) {
      gensubjets.push_back(gensubjet);
    }
  }
  event.set(h_gensubjets, gensubjets);

  vector<Jet> subjets;
  vector<TopJet> *topjets = use_additional_branch_for_rec ? &event.get(h_topjets) : event.topjets;
  int subjetItr = 0;
  vector<pair<int, int>> subjets_map;
  for(unsigned int topjetItr = 0; topjetItr < topjets->size(); ++topjetItr) {
    for(const Jet & subjet : topjets->at(topjetItr).subjets()) {
      subjets.push_back(subjet);
      subjets_map.push_back(pair<int, int>({subjetItr++, topjetItr}));
    }
  }
  event.set(h_subjets, subjets);
  event.set(h_subjets_map, subjets_map);
}


void AK8Corrections::reset_smeared_subjets(Event & event) {

  vector<Jet> subjets = event.get(h_subjets);
  vector<pair<int, int>> subjets_map = event.get(h_subjets_map);
  int topjetItr = 0;
  vector<TopJet> *topjets = use_additional_branch_for_rec ? &event.get(h_topjets) : event.topjets;
  for(auto & topjet : *topjets) {
    vector<Jet> new_subjets;
    for(auto & p : subjets_map) {
      if(p.second == topjetItr) {
        new_subjets.push_back(subjets.at(p.first));
      }
    }
    topjet.set_subjets(new_subjets);
    ++topjetItr;
  }
}


void AK8Corrections::rebuild_topjets_from_subjets(Event & event) {

  vector<TopJet> *topjets = use_additional_branch_for_rec ? &event.get(h_topjets) : event.topjets;
  for(auto & topjet : *topjets) {
    LorentzVector v4;
    for(const auto & subjet : topjet.subjets()) {
      v4 += subjet.v4();
    }
    double jec_factor_raw = topjet.v4().Pt() / v4.Pt();
    topjet.set_JEC_factor_raw(jec_factor_raw);
    topjet.set_v4(v4);
  }
}


bool AK8Corrections::process(Event & event) {

  if(!init_done) {
    throw runtime_error("AK8Corrections::init() not called!");
  }

  if(correct_topjets) {
    if(is_mc) {
      tjet_corrector_MC->process(event);
      tjet_resolution_smearer->process(event);
    }
    else {
      tjet_corrector_data->process(event);
    }
  }

  if(correct_subjets) {
    if(is_mc) {
      subjet_corrector_MC->process(event);
      set_subjet_handles(event);
      subjet_resolution_smearer->process(event);
      reset_smeared_subjets(event);
    }
    else {
      subjet_corrector_data->process(event);
    }
  }

  if(do_rebuild_topjets_from_subjets) {
    rebuild_topjets_from_subjets(event);
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


bool AK8Cleaning::process(Event & event) {

  vector<TopJet> initial_topjets = event.get(h_ak8jets_rec);
  vector<TopJet> cleaned_topjets;
  const auto & primlep = event.get(h_primlep);

  for(const auto & tjet : initial_topjets) {
    if(deltaR(tjet, primlep) < _minDR) continue;
    if(abs(tjet.v4().eta()) > _maxEta) continue;
    if(tjet.v4().pt() < _minPt) continue;
    cleaned_topjets.push_back(tjet);
  }
  if(jet_pt_sorter) sort_by_pt<TopJet>(cleaned_topjets);

  event.set(h_ak8jets_rec, cleaned_topjets);

  return true;
}

}}
