#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/MCWeight.h"

#include "UHH2/LegacyTopTagging/include/Utils.h"

using namespace std;
using namespace uhh2;
using namespace ltt;


namespace uhh2 { namespace ltt {

//____________________________________________________________________________________________________
double tau32(const TopJet & topjet) {
  return min((double)(topjet.tau3() / topjet.tau2()), 0.99999);
}

//____________________________________________________________________________________________________
double tau32groomed(const TopJet & topjet) {
  return min((double)(topjet.tau3_groomed() / topjet.tau2_groomed()), 0.99999);
}

//____________________________________________________________________________________________________
double mSD(const TopJet & topjet) {
  LorentzVector subjet_sum(0,0,0,0);
  for(auto subjet : topjet.subjets()) {
      subjet_sum += subjet.v4();
  }
  return inv_mass_safe(subjet_sum);
}

//____________________________________________________________________________________________________
double maxDeepCSVSubJetValue(const TopJet & topjet) {
  double result(-99999.);
  for(auto subjet : topjet.subjets()) {
    if(subjet.btag_DeepCSV() > result) result = subjet.btag_DeepCSV();
  }
  return min(result, 0.99999);
}

//____________________________________________________________________________________________________
double HOTVR_mpair(const TopJet & topjet, const bool safe) {
  vector<Jet> subjets = topjet.subjets();
  if(subjets.size() < 3) {
    if(safe) throw runtime_error("HOTVR jet has less than 3 subjets, cannot calculate minimum pairwise mass.");
    else return -1.;
  }
  sort_by_pt(subjets);
  double m01 = (subjets.at(0).v4() + subjets.at(1).v4()).M();
  double m02 = (subjets.at(0).v4() + subjets.at(2).v4()).M();
  double m12 = (subjets.at(1).v4() + subjets.at(2).v4()).M();
  return min(m01, min(m02, m12));
}

//____________________________________________________________________________________________________
double HOTVR_fpt(const TopJet & topjet, const unsigned int subjet_i) {
  vector<Jet> subjets = topjet.subjets();
  sort_by_pt(subjets);
  return subjets.at(subjet_i).v4().Pt() / topjet.v4().Pt();
}

//____________________________________________________________________________________________________
HOTVRTopTag::HOTVRTopTag(const double _mass_min, const double _mass_max, const double _fpt_max, const double _mpair_min):
  mass_min(_mass_min), mass_max(_mass_max), fpt_max(_fpt_max), mpair_min(_mpair_min) {}

bool HOTVRTopTag::operator()(const TopJet & jet, const Event & event) const {
  if(jet.subjets().size() < 3) return false;
  if(jet.v4().M() < mass_min || jet.v4().M() > mass_max) return false;
  if(HOTVR_fpt(jet) > fpt_max) return false;
  if(HOTVR_mpair(jet) < mpair_min) return false;
  return true;
}

//____________________________________________________________________________________________________
const TopJet * nextTopJet(const Particle & p, const vector<TopJet> & topjets) {
  return closestParticle(p, topjets);
}

//____________________________________________________________________________________________________
METSelection::METSelection(const double _met_min, const double _met_max): met_min(_met_min), met_max(_met_max) {}

bool METSelection::passes(const Event & event) {
  bool passed_lower_limit = event.met->v4().Pt() > met_min;
  bool passed_upper_limit = event.met->v4().Pt() < met_max;
  return passed_lower_limit && passed_upper_limit;
}

//____________________________________________________________________________________________________
PTWSelection::PTWSelection(Context & ctx, const double _ptw_min, const double _ptw_max): ptw_min(_ptw_min), ptw_max(_ptw_max), h_primlep(ctx.get_handle<FlavorParticle>("PrimaryLepton")) {}

bool PTWSelection::passes(const Event & event) {
  const FlavorParticle & primlep = event.get(h_primlep);
  LorentzVector v4 = event.met->v4() + primlep.v4();
  bool passed_lower_limit = v4.Pt() > ptw_min;
  bool passed_upper_limit = v4.Pt() < ptw_max;
  return passed_lower_limit && passed_upper_limit;
}

//____________________________________________________________________________________________________
TwoDSelection::TwoDSelection(Context & ctx, const double _ptrel_min, const double _dr_min): ptrel_min(_ptrel_min), dr_min(_dr_min), h_primlep(ctx.get_handle<FlavorParticle>("PrimaryLepton")) {}

bool TwoDSelection::passes(const Event & event) {
  const FlavorParticle & primlep = event.get(h_primlep);
  const Jet *nextjet = nextJet(primlep, *event.jets);
  bool passed_ptrel_cut = pTrel(primlep, nextjet) > ptrel_min;
  bool passed_dr_cut = deltaR(primlep.v4(), nextjet->v4()) > dr_min;
  return passed_ptrel_cut || passed_dr_cut;
}

//____________________________________________________________________________________________________
BTagCloseToLeptonSelection::BTagCloseToLeptonSelection(Context & ctx, const double _dr_max, const JetId & _btagID): dr_max(_dr_max), btagID(_btagID), h_primlep(ctx.get_handle<FlavorParticle>("PrimaryLepton")) {}

bool BTagCloseToLeptonSelection::passes(const Event & event) {
  const FlavorParticle & primlep = event.get(h_primlep);
  unsigned int btags_close_to_lepton(0);
  for(const Jet & jet : *event.jets) {
    if(deltaR(jet.v4(), primlep.v4()) > dr_max) continue;
    if(btagID(jet, event)) ++btags_close_to_lepton;
  }
  return btags_close_to_lepton > 0;
}

//____________________________________________________________________________________________________
PartonShowerVariation::PartonShowerVariation(Context & ctx) {
  const string config = ctx.get("SystDirection_PS", "nominal");
  for(const auto & v : kPSVariations) {
    weights_map[v.first] = ctx.declare_event_output<double>("weight_partonshower_"+v.second.name);
    if(v.second.name == config) applied_variation = v.first;
  }
  if(applied_variation != PSVariation::None) {
    cout << "PartonShowerVariation will apply this variation to event weights: " << kPSVariations.at(applied_variation).name << endl;
  }
  else {
    cout << "PartonShowerVariation will only write the weights to the AnalysisTree but not apply any weight." << endl;
  }
}

bool PartonShowerVariation::process(Event & event) {
  bool skip(false);
  if(event.isRealData) skip = true;
  else if(event.genInfo->weights().size() <= 1) {
    if(!warning_thrown) {
      cout << "PartonShowerVariation::process(): No parton shower weights stored for this sample. Doing nothing for all events." << endl;
      warning_thrown = true;
    }
    skip = true;
  }
  if(skip) {
    for(const auto & v : kPSVariations) {
      event.set(weights_map[v.first], 1.);
    }
    return true;
  }
  for(const auto & v : kPSVariations) {
    event.set(weights_map[v.first], event.genInfo->weights().at(v.second.index) / event.genInfo->weights().at(0));
  }
  if(applied_variation != PSVariation::None) event.weight *= event.get(weights_map.at(applied_variation));
  return true;
}

//____________________________________________________________________________________________________
// Twiki: https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL201{6,7,8}
// SF files: https://gitlab.cern.ch/cms-muonPOG/muonefficiencies/-/tree/master/Run2/UL
MuonScaleFactors::MuonScaleFactors(Context & ctx) {
  const string filepath = ctx.get("MuonIDScaleFactorFile");
  const string histname = ctx.get("MuonIDScaleFactorHist");
  const string direction = ctx.get("SystDirection_MuonID", "nominal");
  const string weight_name = "id_tight";
  sf_id.reset(new MCMuonScaleFactor(ctx, filepath, histname, 0.0, weight_name, false, direction));
  // No isolation scale factors applied since we do not use PF muon isolation but our custom 2D cut!
  // sf_iso.reset(new MCMuonScaleFactor(ctx, ...));
  // ...
}

bool MuonScaleFactors::process(Event & event) {
  sf_id->process(event);
  // sf_iso->process(event);
  return true;
}

//____________________________________________________________________________________________________
// Twiki: https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL201{6,7,8}
// SF files: https://gitlab.cern.ch/cms-muonPOG/muonefficiencies/-/tree/master/Run2/UL
TriggerScaleFactors::TriggerScaleFactors(Context & ctx) {
  const string filepath = ctx.get("TriggerScaleFactorFile");
  const string histname = ctx.get("TriggerScaleFactorHist");
  const string direction = ctx.get("SystDirection_MuonTrigger", "nominal");
  const string weight_name = "trigger";
  sf_trig.reset(new MCMuonScaleFactor(ctx, filepath, histname, 0.0, weight_name, false, direction));
}

bool TriggerScaleFactors::process(Event & event) {
  sf_trig->process(event);
  return true;
}

//____________________________________________________________________________________________________
ProbeJetHandleSetter::ProbeJetHandleSetter(Context & ctx, const ProbeJetAlgo & _algo, const string & coll_rec):
  h_probejet(ctx.get_handle<TopJet>("ProbeJet"+kProbeJetAlgos.at(_algo).name)),
  h_topjets(ctx.get_handle<vector<TopJet>>(coll_rec.empty() ? "topjets" : coll_rec)) {}

bool ProbeJetHandleSetter::process(Event & event) {
  vector<TopJet> topjets = event.get(h_topjets);
  sort_by_pt(topjets);
  event.set(h_probejet, topjets.at(0)); // take the leading jet
  return true;
}

//____________________________________________________________________________________________________
void get_Wb_daughters(GenParticle & w_from_top, GenParticle & b_from_top, const GenParticle & top, const vector<GenParticle> & genparticles) {
  w_from_top = *(top.daughter(&genparticles, 1));
  b_from_top = *(top.daughter(&genparticles, 2));
  if(abs(w_from_top.pdgId()) != 24) {
    swap(w_from_top, b_from_top);
  }
  // In rare cases a top-MC event contains more than two genparticles which list a specific top quark as their mother (e.g. due to an additional emission).
  // In order to correctly identify the W boson daughter of a specific top quark, this workaround is needed:
  if(abs(w_from_top.pdgId()) != 24) {
    for(const GenParticle & gp : genparticles) {
      const GenParticle *m1 = gp.mother(&genparticles, 1);
      const GenParticle *m2 = gp.mother(&genparticles, 2);
      const bool has_top_mother = (m1 && m1->index() == top.index()) || (m2 && m2->index() == top.index());
      if(has_top_mother && abs(gp.pdgId()) == 24) {
        w_from_top = gp;
        break;
      }
    }
  }
  if(abs(w_from_top.pdgId()) != 24) {
    throw runtime_error("get_Wb_daughters(): Not able to find correct W boson daughter of top quark.");
  }
  // Do a similar workaround for the b quark daughter:
  if(abs(b_from_top.pdgId()) != 5 && abs(b_from_top.pdgId()) != 3 && abs(b_from_top.pdgId()) != 1) {
    for(const GenParticle & gp : genparticles) {
      const GenParticle *m1 = gp.mother(&genparticles, 1);
      const GenParticle *m2 = gp.mother(&genparticles, 2);
      const bool has_top_mother = (m1 && m1->index() == top.index()) || (m2 && m2->index() == top.index());
      if(has_top_mother && (abs(gp.pdgId()) == 5 || abs(gp.pdgId()) == 3 || abs(gp.pdgId()) == 1)) {
        b_from_top = gp;
        break;
      }
    }
  }
  if(abs(b_from_top.pdgId()) != 5 && abs(b_from_top.pdgId()) != 3 && abs(b_from_top.pdgId()) != 1) {
    throw runtime_error("get_Wb_daughters(): Not able to find correct b quark daughter of top quark.");
  }
}

//____________________________________________________________________________________________________
DecayChannelAndHadronicTopHandleSetter::DecayChannelAndHadronicTopHandleSetter(Context & ctx):
  output_decay_channel(ctx.declare_event_output<DecayChannel>("output_decay_channel")),
  h_hadronictop(ctx.get_handle<GenParticle>("HadronicTopQuark")) // will be unset if process is neither ttbar->l+jets nor single t->hadronic
  {}

bool DecayChannelAndHadronicTopHandleSetter::process(Event & event) {
  DecayChannel dc = DecayChannel::isNotTopQuarkMC;
  if(event.isRealData) {
    event.set(output_decay_channel, dc);
    return true;
  }

  unsigned int n_tops(0);
  for(const GenParticle & gp : *event.genparticles) {
    if(abs(gp.pdgId()) == 6) ++n_tops;
  }
  Process proc = Process::isOther;
  if(n_tops == 2) proc = Process::isTTbar;
  else if(n_tops == 1) proc = Process::isST;

  // TODO: Currently it is not checked whether the two identified W bosons stem from top quarks or not. This will lead to problems when using e.g. ttW samples
  if(proc == Process::isTTbar) {
    GenParticle (*top)(nullptr), (*antitop)(nullptr);
    GenParticle (*w_plus)(nullptr), (*w_minus)(nullptr);
    for(GenParticle & gp : *event.genparticles) {
      if(gp.pdgId() == 6) top = &gp;
      else if(gp.pdgId() == -6) antitop = &gp;
      else if(gp.pdgId() == 24) w_plus = &gp;
      else if(gp.pdgId() == -24) w_minus = &gp;
    }
    const bool w_plus_is_hadronic = abs(w_plus->daughter(event.genparticles, 1)->pdgId()) <= 5;
    const bool w_minus_is_hadronic = abs(w_minus->daughter(event.genparticles, 1)->pdgId()) <= 5;
    if(w_plus_is_hadronic && w_minus_is_hadronic) {
      dc = DecayChannel::isTTbarToHadronic;
    }
    else if(w_plus_is_hadronic && !w_minus_is_hadronic) {
      dc = DecayChannel::isTTbarToSemiLeptonic;
      event.set(h_hadronictop, *top);
    }
    else if(!w_plus_is_hadronic && w_minus_is_hadronic) {
      dc = DecayChannel::isTTbarToSemiLeptonic;
      event.set(h_hadronictop, *antitop);
    }
    else if(!w_plus_is_hadronic && !w_minus_is_hadronic) {
      dc = DecayChannel::isTTbarToDiLeptonic;
    }
  }

  else if(proc == Process::isST) {
    GenParticle (*top)(nullptr);
    for(GenParticle & gp : *event.genparticles) {
      if(abs(gp.pdgId()) == 6) top = &gp;
    }
    GenParticle w_from_top;
    GenParticle b_from_top;
    get_Wb_daughters(w_from_top, b_from_top, *top, *event.genparticles);
    if(abs(w_from_top.daughter(event.genparticles, 1)->pdgId()) <= 5) {
      dc = DecayChannel::isSTToHadronic;
      event.set(h_hadronictop, *top);
    }
    else {
      dc = DecayChannel::isSTToLeptonic;
    }
  }

  event.set(output_decay_channel, dc);
  return true;
}

//____________________________________________________________________________________________________
MergeScenarioHandleSetter::MergeScenarioHandleSetter(Context & ctx, const ProbeJetAlgo & _algo): algo(_algo) {
  h_probejet = ctx.get_handle<TopJet>("ProbeJet"+kProbeJetAlgos.at(_algo).name);
  h_hadronictop = ctx.get_handle<GenParticle>("HadronicTopQuark"); // will be unset if process is neither ttbar->l+jets nor single t->hadronic

  output_has_probejet = ctx.declare_event_output<bool>("output_has_probejet_"+kProbeJetAlgos.at(_algo).name);
  output_merge_scenario = ctx.declare_event_output<MergeScenario>("output_merge_scenario_"+kProbeJetAlgos.at(_algo).name);
}

bool MergeScenarioHandleSetter::process(Event & event) {
  if(event.is_valid(h_probejet)) event.set(output_has_probejet, true);
  else event.set(output_has_probejet, false);

  MergeScenario msc = MergeScenario::isBackground;
  if(event.isRealData || !event.is_valid(h_hadronictop) || !event.is_valid(h_probejet)) {
    event.set(output_merge_scenario, msc);
    return true;
  }

  const TopJet probejet = event.get(h_probejet);
  const GenParticle hadronic_top = event.get(h_hadronictop);
  double dRmatch(-1.);
  if(algo == ProbeJetAlgo::isAK8) {
    dRmatch = 0.8;
  }
  else if(algo == ProbeJetAlgo::isHOTVR) {
    dRmatch = min(1.5, max(0.1, 600./(probejet.v4().pt()*probejet.JEC_factor_raw())));
  }
  GenParticle gen_w;
  GenParticle gen_b;
  get_Wb_daughters(gen_w, gen_b, hadronic_top, *event.genparticles);
  GenParticle gen_q1 = *(gen_w.daughter(event.genparticles, 1));
  GenParticle gen_q2 = *(gen_w.daughter(event.genparticles, 2));
  bool merged_b = deltaR(gen_b.v4(), probejet.v4()) < dRmatch;
  bool merged_q1 = deltaR(gen_q1.v4(), probejet.v4()) < dRmatch;
  bool merged_q2 = deltaR(gen_q2.v4(), probejet.v4()) < dRmatch;
  // now check the merge scenarios:
  if(merged_b && merged_q1 && merged_q2) msc = MergeScenario::isFullyMerged;
  else if(!merged_b && merged_q1 && merged_q2) msc = MergeScenario::isWMerged;
  else if(merged_b && merged_q1 && !merged_q2) msc = MergeScenario::isQBMerged;
  else if(merged_b && !merged_q1 && merged_q2) msc = MergeScenario::isQBMerged;
  else msc = MergeScenario::isNotMerged;

  event.set(output_merge_scenario, msc);
  return true;
}

//____________________________________________________________________________________________________
MainOutputSetter::MainOutputSetter(Context & ctx) {
  h_probejet_hotvr = ctx.get_handle<TopJet>("ProbeJet"+kProbeJetAlgos.at(ProbeJetAlgo::isHOTVR).name);
  h_probejet_ak8 = ctx.get_handle<TopJet>("ProbeJet"+kProbeJetAlgos.at(ProbeJetAlgo::isAK8).name);

  vector<string> output_names;
  const string preprefix = "output_probejet_"; string prefix = "";

  prefix = preprefix+kProbeJetAlgos.at(ProbeJetAlgo::isHOTVR).name+"_";
  output_names.push_back(prefix+"pt");
  output_names.push_back(prefix+"eta");
  output_names.push_back(prefix+"phi");
  output_names.push_back(prefix+"mass");
  output_names.push_back(prefix+"nsub");
  h_probejet_hotvr_nsub_integer = ctx.declare_event_output<int>(prefix+"nsub_integer");
  output_names.push_back(prefix+"mpair");
  output_names.push_back(prefix+"fpt1");
  output_names.push_back(prefix+"tau32");

  prefix = preprefix+kProbeJetAlgos.at(ProbeJetAlgo::isAK8).name+"_";
  output_names.push_back(prefix+"pt");
  output_names.push_back(prefix+"eta");
  output_names.push_back(prefix+"phi");
  output_names.push_back(prefix+"mass");
  output_names.push_back(prefix+"mSD");
  output_names.push_back(prefix+"tau32");
  output_names.push_back(prefix+"maxDeepCSV");

  for(unsigned int i = 0; i < output_names.size(); i++) {
    h_mainoutput.push_back(ctx.declare_event_output<float>(output_names.at(i)));
  }
}

bool MainOutputSetter::process(Event & event) {
  const bool has_hotvr_jet = event.is_valid(h_probejet_hotvr);
  const bool has_ak8_jet = event.is_valid(h_probejet_ak8);
  TopJet probejet_hotvr = TopJet();
  TopJet probejet_ak8 = TopJet();
  if(has_hotvr_jet) probejet_hotvr = event.get(h_probejet_hotvr);
  if(has_ak8_jet) probejet_ak8 = event.get(h_probejet_ak8);

  const double zero_padding(-999.);
  vector<double> values;
  values.resize(h_mainoutput.size(), zero_padding);
  unsigned int i(0);

  values.at(i++) = has_hotvr_jet ? probejet_hotvr.v4().pt() : zero_padding;
  values.at(i++) = has_hotvr_jet ? probejet_hotvr.v4().eta() : zero_padding;
  values.at(i++) = has_hotvr_jet ? probejet_hotvr.v4().phi() : zero_padding;
  values.at(i++) = has_hotvr_jet ? probejet_hotvr.v4().mass() : zero_padding;
  values.at(i++) = has_hotvr_jet ? probejet_hotvr.subjets().size() : zero_padding;
  event.set(h_probejet_hotvr_nsub_integer, has_hotvr_jet ? probejet_hotvr.subjets().size() : zero_padding);
  values.at(i++) = has_hotvr_jet ? HOTVR_mpair(probejet_hotvr, false) : zero_padding;
  values.at(i++) = has_hotvr_jet ? HOTVR_fpt(probejet_hotvr) : zero_padding;
  values.at(i++) = has_hotvr_jet ? tau32groomed(probejet_hotvr) : zero_padding;

  values.at(i++) = has_ak8_jet ? probejet_ak8.v4().pt() : zero_padding;
  values.at(i++) = has_ak8_jet ? probejet_ak8.v4().eta() : zero_padding;
  values.at(i++) = has_ak8_jet ? probejet_ak8.v4().phi() : zero_padding;
  values.at(i++) = has_ak8_jet ? probejet_ak8.v4().mass() : zero_padding;
  values.at(i++) = has_ak8_jet ? mSD(probejet_ak8) : zero_padding;
  values.at(i++) = has_ak8_jet ? tau32(probejet_ak8) : zero_padding;
  values.at(i++) = has_ak8_jet ? maxDeepCSVSubJetValue(probejet_ak8) : zero_padding;

  for(unsigned int i = 0; i < values.size(); i++) {
    event.set(h_mainoutput.at(i), values.at(i));
  }
  return true;
}

//____________________________________________________________________________________________________
HEM2018Selection::HEM2018Selection(Context & ctx): fYear(extract_year(ctx)) {}

// Caveat: Returns "true" if event is affected by the HEM issue.
bool HEM2018Selection::passes(const Event & event) {
  if(fYear == Year::isUL18 && ((event.isRealData && event.run >= fRunNumber) || !event.isRealData)) {
    for(const Jet & jet : *event.jets) {
      if(jet.v4().eta() > fEtaRange.first && jet.v4().eta() < fEtaRange.second && jet.v4().phi() > fPhiRange.first && jet.v4().phi() < fPhiRange.second) {
        return true;
      }
    }
  }
  return false;
}

//____________________________________________________________________________________________________
TopPtReweighting::TopPtReweighting(Context & ctx, const bool apply): fApply(apply) {
  const string config = ctx.get("SystDirection_TopPt", "nominal");
  h_weight_nominal = ctx.declare_event_output<double>("weight_toppt");
  h_weight_a_up    = ctx.declare_event_output<double>("weight_toppt_a_up");
  h_weight_a_down  = ctx.declare_event_output<double>("weight_toppt_a_down");
  h_weight_b_up    = ctx.declare_event_output<double>("weight_toppt_b_up");
  h_weight_b_down  = ctx.declare_event_output<double>("weight_toppt_b_down");
  if(config == "nominal") {
    applied_variation = TopPtVariation::nominal;
  }
  else if(config == "a_up") {
    applied_variation = TopPtVariation::a_up;
  }
  else if(config == "a_down") {
    applied_variation = TopPtVariation::a_down;
  }
  else if(config == "b_up") {
    applied_variation = TopPtVariation::b_up;
  }
  else if(config == "b_down") {
    applied_variation = TopPtVariation::b_down;
  }
  else {
    throw invalid_argument("TopPtReweighting: Invalid systematic variation given in XML config.");
  }
}

void TopPtReweighting::set_dummy_weights(Event & event) {
  event.set(h_weight_nominal, fDummyWeight);
  event.set(h_weight_a_up,    fDummyWeight);
  event.set(h_weight_a_down,  fDummyWeight);
  event.set(h_weight_b_up,    fDummyWeight);
  event.set(h_weight_b_down,  fDummyWeight);
}

bool TopPtReweighting::process(Event & event) {
  if(event.isRealData) {
    set_dummy_weights(event);
    return true;
  }

  unsigned int n_tops(0);
  double top_pt(-1.);
  double antitop_pt(-1.);
  for(const GenParticle & gp : *event.genparticles) {
    if(abs(gp.pdgId()) == 6) ++n_tops;
    if(gp.pdgId() == 6) {
      top_pt = gp.v4().pt();
    }
    else if(gp.pdgId() == -6) {
      antitop_pt = gp.v4().pt();
    }
  }
  if(fPtCutOff_b) {
    top_pt = min(top_pt, fPtCutOff);
    antitop_pt = min(antitop_pt, fPtCutOff);
  }

  Process proc = Process::isOther;
  if(n_tops == 2) proc = Process::isTTbar;
  else if(n_tops == 1) proc = Process::isST;
  if(proc != Process::isTTbar) {
    set_dummy_weights(event);
    return true;
  }
  if(top_pt < 0 || antitop_pt < 0) throw runtime_error("TopPtReweighting::process(): top or antitop pT is negative.");

  const double sf_top_nominal     = exp(fA+fB*top_pt);
  const double sf_top_a_up        = exp(fA_up+fB*top_pt);
  const double sf_top_a_down      = exp(fA_down+fB*top_pt);
  const double sf_top_b_up        = exp(fA+fB_up*top_pt);
  const double sf_top_b_down      = exp(fA+fB_down*top_pt);
  const double sf_antitop_nominal = exp(fA+fB*antitop_pt);
  const double sf_antitop_a_up    = exp(fA_up+fB*antitop_pt);
  const double sf_antitop_a_down  = exp(fA_down+fB*antitop_pt);
  const double sf_antitop_b_up    = exp(fA+fB_up*antitop_pt);
  const double sf_antitop_b_down  = exp(fA+fB_down*antitop_pt);

  const double sf_nominal = sqrt(sf_top_nominal*sf_antitop_nominal);
  const double sf_a_up    = sqrt(sf_top_a_up*sf_antitop_a_up);
  const double sf_a_down  = sqrt(sf_top_a_down*sf_antitop_a_down);
  const double sf_b_up    = sqrt(sf_top_b_up*sf_antitop_b_up);
  const double sf_b_down  = sqrt(sf_top_b_down*sf_antitop_b_down);

  event.set(h_weight_nominal, sf_nominal);
  event.set(h_weight_a_up,    sf_a_up);
  event.set(h_weight_a_down,  sf_a_down);
  event.set(h_weight_b_up,    sf_b_up);
  event.set(h_weight_b_down,  sf_b_down);

  if(fApply) {
    if(applied_variation == TopPtVariation::nominal) {
      event.weight *= sf_nominal;
    }
    else if(applied_variation == TopPtVariation::a_up) {
      event.weight *= sf_a_up;
    }
    else if(applied_variation == TopPtVariation::a_down) {
      event.weight *= sf_a_down;
    }
    else if(applied_variation == TopPtVariation::b_up) {
      event.weight *= sf_b_up;
    }
    else if(applied_variation == TopPtVariation::b_down) {
      event.weight *= sf_b_down;
    }
  }

  return true;
}

//____________________________________________________________________________________________________
// Copy of https://github.com/MatthiesC/HighPtSingleTop/blob/master/src/TheoryCorrections.cxx
VJetsReweighting::VJetsReweighting(Context & ctx, const string& weight_name):
  is_2016_nonUL(extract_year(ctx) == Year::is2016v3),
  is_WJets(ctx.get("dataset_version").find("WJets") == 0),
  is_DYJets(ctx.get("dataset_version").find("DYJets") == 0),
  apply_EWK(string2bool(ctx.get("VJetsReweighting_do_EWK"))),
  apply_QCD_EWK(string2bool(ctx.get("VJetsReweighting_do_QCD_EWK"))),
  apply_QCD_NLO(string2bool(ctx.get("VJetsReweighting_do_QCD_NLO"))),
  apply_QCD_NNLO(string2bool(ctx.get("VJetsReweighting_do_QCD_NNLO"))),
  h_weight_EWK(ctx.declare_event_output<double>(weight_name+"_EWK")),
  h_weight_QCD_EWK(ctx.declare_event_output<double>(weight_name+"_QCD_EWK")),
  h_weight_QCD_NLO(ctx.declare_event_output<double>(weight_name+"_QCD_NLO")),
  h_weight_QCD_NNLO(ctx.declare_event_output<double>(weight_name+"_QCD_NNLO"))
{
  if((apply_QCD_EWK && (apply_EWK || apply_QCD_NLO)) || (apply_QCD_NNLO && !(apply_QCD_EWK || (apply_EWK && apply_QCD_NLO)))) {
    throw invalid_argument("VJetsReweighting: You are not allowed to use the specified combination of correction scale factors.");
  }

  string filesDir = ctx.get("uhh2Dir")+"LegacyTopTagging/data/ScaleFactors/VJetsCorrections/";
  for(const string& proc : {"w", "z"}) {
    TFile* file = new TFile((filesDir+"merged_kfactors_"+proc+"jets.root").c_str());
    for(const string& corr : {"ewk", "qcd", "qcd_ewk"}) load_histo(file, (string)(proc+"_"+corr), (string)("kfactor_monojet_"+corr));
    file->Close();
  }
  for(const string& proc : {"dy", "znn"}) {
    TFile* file = new TFile((filesDir+"kfac_"+proc+"_filter.root").c_str());
    load_histo(file, (string)(proc+"_qcd_2017"), (string)("kfac_"+proc+"_filter"));
    file->Close();
  }
  TFile* file = new TFile((filesDir+"2017_gen_v_pt_qcd_sf.root").c_str());
  load_histo(file, "w_qcd_2017", "wjet_dress_inclusive");
  file->Close();
  file = new TFile((filesDir+"lindert_qcd_nnlo_sf.root").c_str());
  for(const string& proc : {"eej", "evj", "vvj"}) load_histo(file, (string)(proc+"_qcd_nnlo"), (string)(proc));
  file->Close();
}

void VJetsReweighting::load_histo(TFile* file, const string& name, const string& histName) {

  histos[name].reset((TH1F*)file->Get(histName.c_str()));
  histos[name]->SetDirectory(0);
}

double VJetsReweighting::get_v_pt(Event & event) {

  if(!(is_DYJets || is_WJets)) throw runtime_error("VJetsReweighting::get_v_pt(): Calling this function on non-WJets/DYJets sample makes no sense.");
  double pt(-1.);
  bool v_found(false);
  for(const GenParticle & gp : *event.genparticles) {
    if(is_WJets && gp.status() == 22 && abs(gp.pdgId()) == 24) {
      pt = gp.v4().Pt();
      v_found = true;
    }
    else if(is_DYJets && gp.status() == 22 && abs(gp.pdgId()) == 23) {
      pt = gp.v4().Pt();
      v_found = true;
    }
  }
  if(!v_found) {
    int n_status23_leptons(0);
    GenParticle d1, d2; // daughters of V boson
    for(const GenParticle & gp : *event.genparticles) {
      if(gp.status() == 23 && abs(gp.pdgId()) >= 11 && abs(gp.pdgId()) <= 16) {
        n_status23_leptons++;
        if(gp.pdgId() > 0) d1 = gp;
        else d2 = gp;
      }
    }
    if(n_status23_leptons != 2) throw runtime_error("VJetsReweighting::get_v_pt(): Did not find exactly two V daughter candidates.");
    pt = (d1.v4() + d2.v4()).Pt();
  }

  return pt;
}

double VJetsReweighting::evaluate(const string& name, double pt) {

  int firstBin = 1;
  int lastBin = histos[name]->GetNbinsX();
  double h_min = histos[name]->GetBinCenter(firstBin)-0.5*histos[name]->GetBinWidth(firstBin);
  double h_max = histos[name]->GetBinCenter(lastBin)+0.5*histos[name]->GetBinWidth(lastBin);
  double pt_for_eval = pt;
  pt_for_eval = (pt_for_eval > h_min) ? pt_for_eval : h_min+0.001;
  pt_for_eval = (pt_for_eval < h_max) ? pt_for_eval : h_max-0.001;

  return histos[name]->GetBinContent(histos[name]->FindBin(pt_for_eval));
}

bool VJetsReweighting::process(Event & event) {

  double weight_EWK(1.);
  double weight_QCD_EWK(1.);
  double weight_QCD_NLO(1.);
  double weight_QCD_NNLO(1.);

  if(is_WJets || is_DYJets) {
    double pt = get_v_pt(event);
    string process = "";

    // QCD EWK
    if(is_WJets) process = "w";
    if(is_DYJets) process = "z";
    weight_QCD_EWK = evaluate(process+"_qcd_ewk", pt);
    if(apply_QCD_EWK) event.weight *= weight_QCD_EWK;

    // EWK
    weight_EWK = evaluate(process+"_ewk", pt);
    if(apply_EWK) event.weight *= weight_EWK;

    // QCD NLO
    if(!is_2016_nonUL && is_DYJets) process = "dy";
    weight_QCD_NLO = evaluate(process+"_qcd"+(is_2016_nonUL ? "" : "_2017"), pt);
    if(apply_QCD_NLO) event.weight *= weight_QCD_NLO;

    // QCD NNLO
    if(is_DYJets) process = "eej";
    if(is_WJets) process = "evj";
    weight_QCD_NNLO = evaluate(process+"_qcd_nnlo", pt);
    if(apply_QCD_NNLO) event.weight *= weight_QCD_NNLO;
  }

  event.set(h_weight_EWK, weight_EWK);
  event.set(h_weight_QCD_EWK, weight_QCD_EWK);
  event.set(h_weight_QCD_NLO, weight_QCD_NLO);
  event.set(h_weight_QCD_NNLO, weight_QCD_NNLO);

  return true;
}

}}
