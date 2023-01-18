#include <iomanip>

#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/AdditionalSelections.h"

#include "UHH2/LegacyTopTagging/include/Utils.h"

#include <TVectorF.h>

using namespace std;
using namespace uhh2;
using namespace ltt;


namespace uhh2 { namespace ltt {

//____________________________________________________________________________________________________
Channel extract_channel(const Context & ctx) {
  Channel channel = Channel::notValid;
  const TString config = string2lowercase(ctx.get("analysis_channel"));
  if(config.Contains("ele")) channel = Channel::isEle;
  if(config.Contains("muo")) channel = Channel::isMuo;
  if(channel == Channel::notValid) {
    std::runtime_error("extract_channel(): Invalid channel string in xml file. Please check.");
  }
  return channel;
}

//____________________________________________________________________________________________________
float rapidity(const Particle & p) {
  return p.v4().Rapidity();
}

//____________________________________________________________________________________________________
Particle add_Particles(const Particle & p1, const Particle & p2) {
  Particle sum;
  sum.set_v4(p1.v4() + p2.v4());
  sum.set_charge(p1.charge() + p2.charge());
  return sum;
}

//____________________________________________________________________________________________________
void print_jet_info(const Jet & jet, const string & prefix) {
  cout << scientific << prefix+"pt = " << jet.v4().pt() << " ; JEC_factor_raw = " << jet.JEC_factor_raw() << endl;
}

//____________________________________________________________________________________________________
void print_topjet_info(const TopJet & jet, const string & prefix) {
  print_jet_info(jet, prefix);
  cout << prefix+" subjets:" << endl;
  for(const Jet & subjet : jet.subjets()) print_jet_info(subjet, prefix+"  ");
}

//____________________________________________________________________________________________________
double tau32(const TopJet & topjet) {
  return min((double)(topjet.tau3() / topjet.tau2()), 0.99999);
}

//____________________________________________________________________________________________________
double tau21(const TopJet & topjet) {
  return min((double)(topjet.tau2() / topjet.tau1()), 0.99999);
}

//____________________________________________________________________________________________________
double tau32groomed(const TopJet & topjet) {
  return min((double)(topjet.tau3_groomed() / topjet.tau2_groomed()), 0.99999);
}

//____________________________________________________________________________________________________
double tau21groomed(const TopJet & topjet) {
  return min((double)(topjet.tau2_groomed() / topjet.tau1_groomed()), 0.99999);
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
double mTW(const FlavorParticle & lepton_, const MET & met_) {
  LorentzVector lepton = lepton_.v4();
  TVector3 lepton_pT = uhh2::toVector(lepton); // [PxPyPz]
  lepton_pT.SetZ(0);
  LorentzVector met(met_.pt(), 0., met_.phi(), met_.pt()); // [PtEtaPhiE4D] // LorentzVector met = met_.v4(); // does not work with const MET since MET.h is missing a const qualifier for v4() member function, causing compiler error
  TVector3 met_pT = uhh2::toVector(met); // [PxPyPz]
  met_pT.SetZ(0);
  return sqrt(2.*(lepton_pT.Mag() * met_pT.Mag() - lepton_pT * met_pT));
}

//____________________________________________________________________________________________________
double pTW(const FlavorParticle & lepton_, const MET & met_) {
  return (lepton_.v4() + met_.v4()).Pt();
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
double maxDeepJetSubJetValue(const TopJet & topjet) {
  double result(-99999.);
  for(auto subjet : topjet.subjets()) {
    if(subjet.btag_DeepJet() > result) result = subjet.btag_DeepJet();
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
  const double m01 = (subjets.at(0).v4() + subjets.at(1).v4()).M();
  const double m02 = (subjets.at(0).v4() + subjets.at(2).v4()).M();
  const double m12 = (subjets.at(1).v4() + subjets.at(2).v4()).M();
  return min(m01, min(m02, m12));
}

//____________________________________________________________________________________________________
double HOTVR_fpt(const TopJet & topjet, const unsigned int subjet_i) {
  vector<Jet> subjets = topjet.subjets();
  sort_by_pt(subjets);
  return subjets.at(subjet_i).v4().Pt() / topjet.v4().Pt();
}

//____________________________________________________________________________________________________
double HOTVR_Reff(const TopJet & topjet) {
  return min(1.5, max(0.1, 600. / (topjet.v4().Pt() * topjet.JEC_factor_raw())));
}

//____________________________________________________________________________________________________
double particleNet_TvsWandQCD(const TopJet & topjet) {
  return min((double)(
    (topjet.btag_ParticleNetJetTags_probTbcq() + topjet.btag_ParticleNetJetTags_probTbqq()) /
    (topjet.btag_ParticleNetJetTags_probTbcq() + topjet.btag_ParticleNetJetTags_probTbqq() + topjet.btag_ParticleNetJetTags_probWcq() + topjet.btag_ParticleNetJetTags_probWqq() + topjet.btag_ParticleNetJetTags_probQCD())
  ), 0.99999);
}

//____________________________________________________________________________________________________
double particleNet_WvsTandQCD(const TopJet & topjet) {
  return min((double)(
    (topjet.btag_ParticleNetJetTags_probWcq() + topjet.btag_ParticleNetJetTags_probWqq()) /
    (topjet.btag_ParticleNetJetTags_probWcq() + topjet.btag_ParticleNetJetTags_probWqq() + topjet.btag_ParticleNetJetTags_probTbcq() + topjet.btag_ParticleNetJetTags_probTbqq() + topjet.btag_ParticleNetJetTags_probQCD())
  ), 0.99999);
}

//____________________________________________________________________________________________________
HOTVRTopTag::HOTVRTopTag(const double _mass_min, const double _mass_max, const double _fpt_max, const double _mpair_min):
  mass_min(_mass_min), mass_max(_mass_max), fpt_max(_fpt_max), mpair_min(_mpair_min) {}

bool HOTVRTopTag::operator()(const TopJet & jet, const Event & event) const {
  (void) event;
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
// Copy of NoLeptonInJet from common/src/JetIds.cxx but reduced to only asking for deltaR
NoLeptonInJet::NoLeptonInJet(const string & _lepton, const double _dr, const boost::optional<ElectronId> & _ele_id, const boost::optional<MuonId> & _muo_id):
  lepton(_lepton), dr(_dr), ele_id(_ele_id), muo_id(_muo_id) {}

bool NoLeptonInJet::operator()(const Jet & jet, const Event & event) const {

  const bool doMuons = event.muons && (lepton == "muon" || lepton == "all");
  const bool doElectrons = event.electrons && (lepton == "ele" || lepton == "all");
  if(doMuons) {
    for(const auto & muo : *event.muons) {
      if(muo_id && !(*muo_id)(muo, event)) continue;
      if(deltaR(jet, muo) < dr) return false;
    }
  }
  if(doElectrons) {
    for(const auto & ele : *event.electrons) {
      if(ele_id && !(*ele_id)(ele, event)) continue;
      if(deltaR(jet, ele) < dr) return false;
    }
  }
  return true;
}

//____________________________________________________________________________________________________
METSelection::METSelection(Context & ctx, const boost::optional<double> & _met_min, const boost::optional<double> & _met_max, const boost::optional<string> & _met_name): met_min(_met_min), met_max(_met_max), h_met(ctx.get_handle<MET>(_met_name ? *_met_name : "met")) {}

bool METSelection::passes(const Event & event) {
  const MET *met = &event.get(h_met);
  bool passed_lower_limit(true);
  bool passed_upper_limit(true);
  if(met_min) passed_lower_limit = met->pt() > *met_min;
  if(met_max) passed_upper_limit = met->pt() < *met_max;
  return passed_lower_limit && passed_upper_limit;
}

//____________________________________________________________________________________________________
PTWSelection::PTWSelection(Context & ctx, const double _ptw_min, const double _ptw_max): ptw_min(_ptw_min), ptw_max(_ptw_max), h_primlep(ctx.get_handle<FlavorParticle>(kHandleName_PrimaryLepton)) {}

bool PTWSelection::passes(const Event & event) {
  const FlavorParticle & primlep = event.get(h_primlep);
  const float ptw = pTW(primlep, *event.met);
  const bool passed_lower_limit = ptw > ptw_min;
  const bool passed_upper_limit = ptw < ptw_max;
  return passed_lower_limit && passed_upper_limit;
}

//____________________________________________________________________________________________________
TwoDSelection::TwoDSelection(Context & ctx, const double _ptrel_min, const double _dr_min, const bool _circular):
  ptrel_min(_ptrel_min),
  dr_min(_dr_min),
  circular(_circular),
  h_primlep(ctx.get_handle<FlavorParticle>(kHandleName_PrimaryLepton)),
  h_jets(ctx.get_handle<vector<Jet>>(kHandleName_pairedPUPPIjets))
{}

bool TwoDSelection::passes(const Event & event) {
  const FlavorParticle & primlep = event.get(h_primlep);
  // throw runtime_error("fix TwoDSelection with new PUPPI CHS matching setup"); // DONE
  const vector<Jet> & jets = event.get(h_jets);
  const Jet *nextjet = nextJet(primlep, jets);
  const float ptrel = pTrel(primlep, nextjet);
  const float dr = deltaR(primlep.v4(), nextjet->v4());
  const bool passed_ptrel_cut = ptrel > ptrel_min;
  const bool passed_dr_cut = dr > dr_min;
  const bool passes_circular = dr*dr / (dr_min*dr_min) + ptrel*ptrel / (ptrel_min*ptrel_min) > 1.f;
  if(circular) return passes_circular;
  else return passed_ptrel_cut || passed_dr_cut;
}

//____________________________________________________________________________________________________
BTagCloseToLeptonSelection::BTagCloseToLeptonSelection(Context & ctx, const double _dr_max):
  dr_max(_dr_max),
  h_primlep(ctx.get_handle<FlavorParticle>(kHandleName_PrimaryLepton)),
  h_bJets(ctx.get_handle<vector<Jet>>(kHandleName_bJets))
{}

bool BTagCloseToLeptonSelection::passes(const Event & event) {
  const FlavorParticle & primlep = event.get(h_primlep);
  const vector<Jet> & bjets = event.get(h_bJets);
  unsigned int btags_close_to_lepton(0);
  for(const Jet & jet : bjets) {
    if(deltaR(jet.v4(), primlep.v4()) < dr_max) ++btags_close_to_lepton;
  }
  return btags_close_to_lepton > 0;
}

//____________________________________________________________________________________________________
PartonShowerVariation::PartonShowerVariation(Context & ctx) {
  const string config = ctx.get("SystDirection_PS", "nominal");
  for(const auto & v : kPSVariations) {
    weights_map[v.first] = ctx.declare_event_output<float>("weight_partonshower_"+v.second.name);
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
LttTriggerSelection::LttTriggerSelection(const Context & ctx): fYear(extract_year(ctx)) {
  slct_trigger_Mu50.reset(new TriggerSelection("HLT_Mu50_v*"));
  slct_trigger_TkMu50.reset(new TriggerSelection("HLT_TkMu50_v*"));
  slct_trigger_OldMu100.reset(new TriggerSelection("HLT_OldMu100_v*"));
  // slct_trigger_Mu100.reset(new TriggerSelection("HLT_Mu100_v*"));
  slct_trigger_TkMu100.reset(new TriggerSelection("HLT_TkMu100_v*"));
}

bool LttTriggerSelection::passes(const Event & event) {
  if(fYear == Year::isUL16preVFP || fYear == Year::isUL16postVFP) {
    if(slct_trigger_Mu50->passes(event) || slct_trigger_TkMu50->passes(event)) {
      return true;
    }
  }
  else if(fYear == Year::isUL17 || fYear == Year::isUL18) {
    if(slct_trigger_Mu50->passes(event) || slct_trigger_OldMu100->passes(event) || slct_trigger_TkMu100->passes(event)) {
      return true;
    }
  }
  else {
    throw runtime_error("uhh2::ltt::LttTriggerSelection::passes(): No trigger paths given for this year");
  }
  return true;
}

//____________________________________________________________________________________________________
// Twiki: https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL201{6,7,8}
// SF files: https://gitlab.cern.ch/cms-muonPOG/muonefficiencies/-/tree/master/Run2/UL
// The centrally provided scale factors are derived for the trigger paths given in uhh2::ltt::TriggerSelection, using the HighPt muon ID. However, we use the tight muon ID in this analysis. Be aware of that!
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
  h_decay_channel(ctx.get_handle<DecayChannel>("decay_channel")),
  h_hadronictop(ctx.get_handle<GenParticle>("HadronicTopQuark")) // will be unset if process is neither ttbar->l+jets nor single t->hadronic
  {}

bool DecayChannelAndHadronicTopHandleSetter::process(Event & event) {
  DecayChannel dc = DecayChannel::isNotTopQuarkMC;
  if(event.isRealData) {
    event.set(h_decay_channel, dc);
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

  event.set(h_decay_channel, dc);
  return true;
}

//____________________________________________________________________________________________________
MergeScenarioHandleSetter::MergeScenarioHandleSetter(Context & ctx, const ProbeJetAlgo & _algo, const string & handle_name_GENtW):
  algo(_algo),
  fHandle_GENtW(ctx.get_handle<ltt::SingleTopGen_tWch>(handle_name_GENtW))
{
  h_probejet = ctx.get_handle<TopJet>("ProbeJet"+kProbeJetAlgos.at(_algo).name);
  h_hadronictop = ctx.get_handle<GenParticle>("HadronicTopQuark"); // will be unset if process is neither ttbar->l+jets nor single t->hadronic

  output_has_probejet = ctx.declare_event_output<bool>("output_has_probejet_"+kProbeJetAlgos.at(_algo).name);
  output_merge_scenario = ctx.declare_event_output<int>("output_merge_scenario_"+kProbeJetAlgos.at(_algo).name);
  h_merge_scenario = ctx.get_handle<MergeScenario>("merge_scenario_"+kProbeJetAlgos.at(_algo).name);
}

bool MergeScenarioHandleSetter::process(Event & event) {
  if(event.is_valid(h_probejet)) event.set(output_has_probejet, true);
  else event.set(output_has_probejet, false);

  MergeScenario msc = MergeScenario::isBackground;
  // if(event.isRealData || !event.is_valid(h_hadronictop) || !event.is_valid(h_probejet)) {
  //   event.set(output_merge_scenario, kMergeScenarios.at(msc).index);
  //   event.set(h_merge_scenario, msc);
  //   return true;
  // }
  if(event.isRealData || !event.is_valid(h_probejet)) {
    event.set(output_merge_scenario, kMergeScenarios.at(msc).index);
    event.set(h_merge_scenario, msc);
    return true;
  }

  const TopJet & probejet = event.get(h_probejet);
  double dRmatch(-1.);
  if(algo == ProbeJetAlgo::isAK8) {
    dRmatch = 0.8;
  }
  else if(algo == ProbeJetAlgo::isHOTVR) {
    dRmatch = HOTVR_Reff(probejet);
  }

  if(!event.is_valid(h_hadronictop)) {
    if(event.is_valid(fHandle_GENtW)) {
       const SingleTopGen_tWch & GENtW = event.get(fHandle_GENtW);
       if(GENtW.IsAssHadronicDecay()) {
         const GenParticle & d1 = GENtW.WAssDecay1();
         const GenParticle & d2 = GENtW.WAssDecay2();
         const GenParticle bAss = GENtW.HasAssociatedBottom() ? GENtW.bAss() : GenParticle();
         const bool merged_d1 = deltaR(d1.v4(), probejet.v4()) < dRmatch;
         const bool merged_d2 = deltaR(d2.v4(), probejet.v4()) < dRmatch;
         const bool merged_bAss = GENtW.HasAssociatedBottom() ? deltaR(bAss.v4(), probejet.v4()) < dRmatch : false;
         // now check the merge scenarios:
         if(merged_bAss && merged_d1 && merged_d2) msc = MergeScenario::isBackground; // we have three quarks merged but they are not all from a top, so we call it background
         else if(!merged_bAss && merged_d1 && merged_d2) msc = MergeScenario::isWMerged;
         else if(merged_bAss && merged_d1 && !merged_d2) msc = MergeScenario::isQBMerged;
         else if(merged_bAss && !merged_d1 && merged_d2) msc = MergeScenario::isQBMerged;
         else msc = MergeScenario::isNotMerged;
       }
    }
    event.set(output_merge_scenario, kMergeScenarios.at(msc).index);
    event.set(h_merge_scenario, msc);
    return true;
  }

  const GenParticle & hadronic_top = event.get(h_hadronictop);
  GenParticle gen_w;
  GenParticle gen_b;
  get_Wb_daughters(gen_w, gen_b, hadronic_top, *event.genparticles);
  GenParticle gen_q1 = *(gen_w.daughter(event.genparticles, 1));
  // GenParticle gen_q2 = *(gen_w.daughter(event.genparticles, 2));
  //___________________
  // There is a very weird bug in a CR2 TTbarToSemiLeptonic UL17 sample for some events where the index of the daughter2 of the W- boson is 65535 and the pythia status code is 0; catch this bug by manually setting v4 of daughter2:
  GenParticle gen_q2;
  if (gen_w.daughter(event.genparticles, 2) != nullptr) gen_q2 = *(gen_w.daughter(event.genparticles, 2));
  else {
    cout << "event.event: " << event.event << endl;
    cout << "Weird GenParticle daughter2. Reconstruct it from daughter1 and mother:" << endl;
    gen_q2.set_v4(gen_w.v4() - gen_q1.v4());
    cout << "pt = " << gen_q2.pt() << endl;
    cout << "eta = " << gen_q2.eta() << endl;
  }
  //___________________
  bool merged_b = deltaR(gen_b.v4(), probejet.v4()) < dRmatch;
  bool merged_q1 = deltaR(gen_q1.v4(), probejet.v4()) < dRmatch;
  bool merged_q2 = deltaR(gen_q2.v4(), probejet.v4()) < dRmatch;
  // now check the merge scenarios:
  if(merged_b && merged_q1 && merged_q2) msc = MergeScenario::isFullyMerged;
  else if(!merged_b && merged_q1 && merged_q2) msc = MergeScenario::isWMerged;
  else if(merged_b && merged_q1 && !merged_q2) msc = MergeScenario::isQBMerged;
  else if(merged_b && !merged_q1 && merged_q2) msc = MergeScenario::isQBMerged;
  else msc = MergeScenario::isNotMerged;

  event.set(output_merge_scenario, kMergeScenarios.at(msc).index);
  event.set(h_merge_scenario, msc);
  return true;
}

//____________________________________________________________________________________________________
MainOutputSetter::MainOutputSetter(Context & ctx):
  h_probejet_hotvr(ctx.get_handle<TopJet>("ProbeJet"+kProbeJetAlgos.at(ProbeJetAlgo::isHOTVR).name)),
  h_probejet_ak8(ctx.get_handle<TopJet>("ProbeJet"+kProbeJetAlgos.at(ProbeJetAlgo::isAK8).name)),
  h_primlep(ctx.get_handle<FlavorParticle>(kHandleName_PrimaryLepton)),
  h_jets(ctx.get_handle<vector<Jet>>(kHandleName_pairedPUPPIjets))
{

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
  output_names.push_back(prefix+"tau21");
  output_names.push_back(prefix+"maxDeepCSV");
  output_names.push_back(prefix+"maxDeepJet");
  output_names.push_back(prefix+"DeepAK8_TvsQCD");
  output_names.push_back(prefix+"DeepAK8_WvsQCD");
  output_names.push_back(prefix+"MDDeepAK8_TvsQCD");
  output_names.push_back(prefix+"MDDeepAK8_WvsQCD");
  output_names.push_back(prefix+"ParticleNet_TvsQCD");
  output_names.push_back(prefix+"ParticleNet_WvsQCD");

  // Other outputs:
  prefix = "output_";
  output_names.push_back(prefix+"lepton_pt");
  output_names.push_back(prefix+"mtw");
  output_names.push_back(prefix+"ptmiss");
  output_names.push_back(prefix+"ptw");
  output_names.push_back(prefix+"2d_ptrel");
  output_names.push_back(prefix+"2d_drjet");

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

  values.at(i++) = has_hotvr_jet ? probejet_hotvr.pt() : zero_padding;
  values.at(i++) = has_hotvr_jet ? probejet_hotvr.eta() : zero_padding;
  values.at(i++) = has_hotvr_jet ? probejet_hotvr.phi() : zero_padding;
  values.at(i++) = has_hotvr_jet ? probejet_hotvr.v4().mass() : zero_padding;
  values.at(i++) = has_hotvr_jet ? probejet_hotvr.subjets().size() : zero_padding;
  event.set(h_probejet_hotvr_nsub_integer, has_hotvr_jet ? probejet_hotvr.subjets().size() : zero_padding);
  values.at(i++) = has_hotvr_jet ? HOTVR_mpair(probejet_hotvr, false) : zero_padding;
  values.at(i++) = has_hotvr_jet ? HOTVR_fpt(probejet_hotvr) : zero_padding;
  values.at(i++) = has_hotvr_jet ? tau32groomed(probejet_hotvr) : zero_padding;

  values.at(i++) = has_ak8_jet ? probejet_ak8.pt() : zero_padding;
  values.at(i++) = has_ak8_jet ? probejet_ak8.eta() : zero_padding;
  values.at(i++) = has_ak8_jet ? probejet_ak8.phi() : zero_padding;
  values.at(i++) = has_ak8_jet ? probejet_ak8.v4().mass() : zero_padding;
  values.at(i++) = has_ak8_jet ? mSD(probejet_ak8) : zero_padding;
  values.at(i++) = has_ak8_jet ? tau32(probejet_ak8) : zero_padding;
  values.at(i++) = has_ak8_jet ? tau21(probejet_ak8) : zero_padding;
  values.at(i++) = has_ak8_jet ? maxDeepCSVSubJetValue(probejet_ak8) : zero_padding;
  values.at(i++) = has_ak8_jet ? maxDeepJetSubJetValue(probejet_ak8) : zero_padding;
  values.at(i++) = has_ak8_jet ? probejet_ak8.btag_DeepBoosted_TvsQCD() : zero_padding;
  values.at(i++) = has_ak8_jet ? probejet_ak8.btag_DeepBoosted_WvsQCD() : zero_padding;
  values.at(i++) = has_ak8_jet ? probejet_ak8.btag_MassDecorrelatedDeepBoosted_TvsQCD() : zero_padding;
  values.at(i++) = has_ak8_jet ? probejet_ak8.btag_MassDecorrelatedDeepBoosted_WvsQCD() : zero_padding;
  values.at(i++) = has_ak8_jet ? probejet_ak8.btag_ParticleNetDiscriminatorsJetTags_TvsQCD() : zero_padding;
  values.at(i++) = has_ak8_jet ? probejet_ak8.btag_ParticleNetDiscriminatorsJetTags_WvsQCD() : zero_padding;

  // Other outputs:
  const FlavorParticle & primlep = event.get(h_primlep);
  values.at(i++) = primlep.pt();
  values.at(i++) = mTW(primlep, *event.met);
  values.at(i++) = event.met->pt();
  values.at(i++) = pTW(primlep, *event.met);
  values.at(i++) = pTrel(primlep, nextJet(primlep, event.get(h_jets)));
  values.at(i++) = deltaR(primlep.v4(), nextJet(primlep, event.get(h_jets))->v4());

  for(unsigned int i = 0; i < values.size(); i++) {
    event.set(h_mainoutput.at(i), values.at(i));
  }
  return true;
}

//____________________________________________________________________________________________________
HEM2018Selection::HEM2018Selection(Context & ctx, const string & handle_name_jets):
  fYear(extract_year(ctx)),
  fHandle_jets(ctx.get_handle<vector<Jet>>(handle_name_jets))
{}

// Caveat: Returns "true" if event is affected by the HEM issue.
bool HEM2018Selection::passes(const Event & event) {
  if(fYear == Year::isUL18 && ((event.isRealData && event.run >= fRunNumber) || !event.isRealData)) {
    vector<Particle> relevant_objects;
    relevant_objects.insert(relevant_objects.end(), event.get(fHandle_jets).begin(), event.get(fHandle_jets).end());
    relevant_objects.insert(relevant_objects.end(), event.electrons->begin(), event.electrons->end());
    for(const Particle & obj : relevant_objects) {
      if(obj.v4().eta() > fEtaRange.first && obj.v4().eta() < fEtaRange.second && obj.v4().phi() > fPhiRange.first && obj.v4().phi() < fPhiRange.second) {
        return true;
      }
    }
  }
  return false;
}

//____________________________________________________________________________________________________
PrefiringWeights::PrefiringWeights(Context & ctx, const bool apply): fYear(extract_year(ctx)), fApply(apply) {
  const string config = ctx.get("SystDirection_Prefiring", "nominal");
  h_weight_nominal = ctx.declare_event_output<float>("weight_prefire");
  h_weight_up      = ctx.declare_event_output<float>("weight_prefire_up");
  h_weight_down    = ctx.declare_event_output<float>("weight_prefire_down");
  if(config == "nominal") {
    applied_variation = PrefireVariation::nominal;
  }
  else if(config == "up") {
    applied_variation = PrefireVariation::up;
  }
  else if(config == "down") {
    applied_variation = PrefireVariation::down;
  }
  else {
    throw invalid_argument("PrefiringWeights: Invalid systematic variation given in XML config.");
  }
}

void PrefiringWeights::set_dummy_weights(Event & event) {
  event.set(h_weight_nominal, fDummyWeight);
  event.set(h_weight_up,      fDummyWeight);
  event.set(h_weight_down,    fDummyWeight);
}

bool PrefiringWeights::process(Event & event) {
  if(event.isRealData || fYear == Year::is2018) { // UL18 still gets prefiring weights because of prefired muons, but not rereco 2018
    set_dummy_weights(event);
    return true;
  }

  const float w_nominal = event.prefiringWeight;
  const float w_up      = event.prefiringWeightUp;
  const float w_down    = event.prefiringWeightDown;

  event.set(h_weight_nominal, w_nominal);
  event.set(h_weight_up     , w_up);
  event.set(h_weight_down   , w_down);

  if(fApply) {
    if(applied_variation == PrefireVariation::nominal) {
      event.weight *= w_nominal;
    }
    else if(applied_variation == PrefireVariation::up) {
      event.weight *= w_up;
    }
    else if(applied_variation == PrefireVariation::down) {
      event.weight *= w_down;
    }
  }

  return true;
}

//____________________________________________________________________________________________________
TopPtReweighting::TopPtReweighting(Context & ctx, const bool apply): fApply(apply) {
  const string config = ctx.get("SystDirection_TopPt", "nominal");
  h_weight_nominal = ctx.declare_event_output<float>("weight_toppt");
  h_weight_a_up    = ctx.declare_event_output<float>("weight_toppt_a_up");
  h_weight_a_down  = ctx.declare_event_output<float>("weight_toppt_a_down");
  h_weight_b_up    = ctx.declare_event_output<float>("weight_toppt_b_up");
  h_weight_b_down  = ctx.declare_event_output<float>("weight_toppt_b_down");
  h_weight_applied = ctx.declare_event_output<float>("weight_toppt_applied");
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
  event.set(h_weight_applied, fDummyWeight);
}

bool TopPtReweighting::process(Event & event) {
  if(event.isRealData) {
    set_dummy_weights(event);
    return true;
  }

  unsigned int n_tops(0);
  float top_pt(-1.);
  float antitop_pt(-1.);
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

  const float sf_top_nominal     = exp(fA+fB*top_pt);
  const float sf_top_a_up        = exp(fA_up+fB*top_pt);
  const float sf_top_a_down      = exp(fA_down+fB*top_pt);
  const float sf_top_b_up        = exp(fA+fB_up*top_pt);
  const float sf_top_b_down      = exp(fA+fB_down*top_pt);
  const float sf_antitop_nominal = exp(fA+fB*antitop_pt);
  const float sf_antitop_a_up    = exp(fA_up+fB*antitop_pt);
  const float sf_antitop_a_down  = exp(fA_down+fB*antitop_pt);
  const float sf_antitop_b_up    = exp(fA+fB_up*antitop_pt);
  const float sf_antitop_b_down  = exp(fA+fB_down*antitop_pt);

  const float sf_nominal = sqrt(sf_top_nominal*sf_antitop_nominal);
  const float sf_a_up    = sqrt(sf_top_a_up*sf_antitop_a_up);
  const float sf_a_down  = sqrt(sf_top_a_down*sf_antitop_a_down);
  const float sf_b_up    = sqrt(sf_top_b_up*sf_antitop_b_up);
  const float sf_b_down  = sqrt(sf_top_b_down*sf_antitop_b_down);
  float sf_applied = 1.0f;

  if(fApply) {
    if(applied_variation == TopPtVariation::nominal) {
      sf_applied *= sf_nominal;
    }
    else if(applied_variation == TopPtVariation::a_up) {
      sf_applied *= sf_a_up;
    }
    else if(applied_variation == TopPtVariation::a_down) {
      sf_applied *= sf_a_down;
    }
    else if(applied_variation == TopPtVariation::b_up) {
      sf_applied *= sf_b_up;
    }
    else if(applied_variation == TopPtVariation::b_down) {
      sf_applied *= sf_b_down;
    }
  }
  event.weight *= sf_applied;

  event.set(h_weight_nominal, sf_nominal);
  event.set(h_weight_a_up,    sf_a_up);
  event.set(h_weight_a_down,  sf_a_down);
  event.set(h_weight_b_up,    sf_b_up);
  event.set(h_weight_b_down,  sf_b_down);
  event.set(h_weight_applied, sf_applied);

  return true;
}

//____________________________________________________________________________________________________
// Copy of https://github.com/MatthiesC/HighPtSingleTop/blob/master/src/TheoryCorrections.cxx
VJetsReweighting::VJetsReweighting(Context & ctx, const string & weight_name):
  is_2016_nonUL(extract_year(ctx) == Year::is2016v3 || extract_year(ctx) == Year::is2016v2),
  is_WJets(ctx.get("dataset_version").find("WJets") == 0),
  is_DYJets(ctx.get("dataset_version").find("DYJets") == 0),
  apply_EWK(string2bool(ctx.get("VJetsReweighting_do_EWK"))),
  apply_QCD_EWK(string2bool(ctx.get("VJetsReweighting_do_QCD_EWK"))),
  apply_QCD_NLO(string2bool(ctx.get("VJetsReweighting_do_QCD_NLO"))),
  apply_QCD_NNLO(string2bool(ctx.get("VJetsReweighting_do_QCD_NNLO"))),
  h_weight_applied(ctx.declare_event_output<float>(weight_name+"_applied")),
  h_weight_EWK(ctx.declare_event_output<float>(weight_name+"_EWK")),
  h_weight_QCD_EWK(ctx.declare_event_output<float>(weight_name+"_QCD_EWK")),
  h_weight_QCD_NLO(ctx.declare_event_output<float>(weight_name+"_QCD_NLO")),
  h_weight_QCD_NNLO(ctx.declare_event_output<float>(weight_name+"_QCD_NNLO"))
{
  if((apply_QCD_EWK && (apply_EWK || apply_QCD_NLO)) || (apply_QCD_NNLO && !(apply_QCD_EWK || (apply_EWK && apply_QCD_NLO)))) {
    throw invalid_argument("VJetsReweighting: You are not allowed to use the specified combination of correction scale factors.");
  }

  string filesDir = (string)getenv("CMSSW_BASE")+"/src/UHH2/LegacyTopTagging/data/ScaleFactors/VJetsCorrections/";
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

double VJetsReweighting::get_v_pt(const Event & event) {

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

double VJetsReweighting::evaluate(const string& name, const double pt) {

  const int firstBin = 1;
  const int lastBin = histos[name]->GetNbinsX();
  const double h_min = histos[name]->GetBinCenter(firstBin)-0.5*histos[name]->GetBinWidth(firstBin);
  const double h_max = histos[name]->GetBinCenter(lastBin)+0.5*histos[name]->GetBinWidth(lastBin);
  double pt_for_eval = pt;
  pt_for_eval = (pt_for_eval > h_min) ? pt_for_eval : h_min+0.001;
  pt_for_eval = (pt_for_eval < h_max) ? pt_for_eval : h_max-0.001;

  return histos[name]->GetBinContent(histos[name]->FindBin(pt_for_eval));
}

bool VJetsReweighting::process(Event & event) {

  float weight_applied(1.);
  float weight_EWK(1.);
  float weight_QCD_EWK(1.);
  float weight_QCD_NLO(1.);
  float weight_QCD_NNLO(1.);

  if(is_WJets || is_DYJets) {
    double pt = get_v_pt(event);
    string process = "";

    // QCD EWK
    if(is_WJets) process = "w";
    if(is_DYJets) process = "z";
    weight_QCD_EWK = evaluate(process+"_qcd_ewk", pt);
    if(apply_QCD_EWK) weight_applied *= weight_QCD_EWK;

    // EWK
    weight_EWK = evaluate(process+"_ewk", pt);
    if(apply_EWK) weight_applied *= weight_EWK;

    // QCD NLO
    if(!is_2016_nonUL && is_DYJets) process = "dy";
    weight_QCD_NLO = evaluate(process+"_qcd"+(is_2016_nonUL ? "" : "_2017"), pt);
    if(apply_QCD_NLO) weight_applied *= weight_QCD_NLO;

    // QCD NNLO
    if(is_DYJets) process = "eej";
    if(is_WJets) process = "evj";
    weight_QCD_NNLO = evaluate(process+"_qcd_nnlo", pt);
    if(apply_QCD_NNLO) weight_applied *= weight_QCD_NNLO;
  }

  event.weight *= weight_applied;

  event.set(h_weight_applied, weight_applied);
  event.set(h_weight_EWK, weight_EWK);
  event.set(h_weight_QCD_EWK, weight_QCD_EWK);
  event.set(h_weight_QCD_NLO, weight_QCD_NLO);
  event.set(h_weight_QCD_NNLO, weight_QCD_NNLO);

  return true;
}

//____________________________________________________________________________________________________
WeightTrickery::WeightTrickery(Context & ctx, const string & handle_name_GENtW, const bool doing_PDF_variations, const bool apply):
  fDoingPDFVariations(doing_PDF_variations),
  fApply(apply),
  fYear(extract_year(ctx)),
  fHandle_weight(ctx.declare_event_output<float>("weight_trickery")),
  fHandle_GENtW(ctx.get_handle<ltt::SingleTopGen_tWch>(handle_name_GENtW))
{
  slct_Mtt700to1000.reset(new MttbarGenSelection(700, 1000));
  slct_Mtt1000toInf.reset(new MttbarGenSelection(1000, -1));

  const string dataset_version = ctx.get("dataset_version");

  const bool is_extra_syst = string2bool(ctx.get("extra_syst"));

  is_TTbar = dataset_version.find("TTbar") == 0;
  is_TTbar_Mtt = is_TTbar && dataset_version.find("Mtt") != string::npos;
  is_TTbar_Mtt700to1000 = is_TTbar_Mtt && dataset_version.find("700to1000") != string::npos;
  is_TTbar_Mtt1000toInf = is_TTbar_Mtt && dataset_version.find("1000toInf") != string::npos;
  // is_TTbar_syst = is_TTbar && dataset_version.find("syst") != string::npos;
  is_TTbar_syst = is_TTbar && is_extra_syst;

  is_tW = dataset_version.find("ST_tW") == 0;
  is_tW_incl = is_tW && dataset_version.find("inclusiveDecays") != string::npos;
  is_tW_nfhd = is_tW && dataset_version.find("NoFullyHadronic") != string::npos;
  is_tW_nfhd_DS = is_tW_nfhd && dataset_version.find("_DS_") != string::npos;
  is_tW_nfhd_PDF = is_tW_nfhd && dataset_version.find("PDF") != string::npos;
  // is_tW_nfhd_syst = is_tW_nfhd && dataset_version.find("syst") != string::npos;
  is_tW_nfhd_syst = is_tW_nfhd && is_extra_syst;
}

bool WeightTrickery::process(Event & event) {

  float weight(1.0f);

  // The regular MiniAODv2 UL TT samples provide ca. 13-30 MC events per expected real TT event.
  // The TT_Mtt samples provide between 6-70 MC events per expected real TT event.
  // Guessing a contribution of 1% all-jets, 89% semilept, and 10% dilept TT after full event selection, I calculated appropriate weights in Excel
  if(is_TTbar) {
    if(!is_TTbar_syst) {
      if(is_TTbar_Mtt700to1000) {
        switch(fYear) {
          case Year::isUL16preVFP :
          weight = 0.446; break;
          case Year::isUL16postVFP :
          weight = 0.513; break;
          case Year::isUL17 :
          weight = 0.315; break;
          case Year::isUL18 :
          weight = 0.226; break;
          default : break;
        }
      }
      else if(is_TTbar_Mtt1000toInf) {
        switch(fYear) {
          case Year::isUL16preVFP :
          weight = 0.761; break;
          case Year::isUL16postVFP :
          weight = 0.749; break;
          case Year::isUL17 :
          weight = 0.534; break;
          case Year::isUL18 :
          weight = 0.470; break;
          default : break;
        }
      }
      else {
        if(slct_Mtt700to1000->passes(event)) {
          switch(fYear) {
            case Year::isUL16preVFP :
            weight = 0.554; break;
            case Year::isUL16postVFP :
            weight = 0.487; break;
            case Year::isUL17 :
            weight = 0.685; break;
            case Year::isUL18 :
            weight = 0.774; break;
            default : break;
          }
        }
        else if(slct_Mtt1000toInf->passes(event)) {
          switch(fYear) {
            case Year::isUL16preVFP :
            weight = 0.239; break;
            case Year::isUL16postVFP :
            weight = 0.251; break;
            case Year::isUL17 :
            weight = 0.466; break;
            case Year::isUL18 :
            weight = 0.530; break;
            default : break;
          }
        }
      }
    }
  }
  // Get the maximum MC statistics for NoFullyHadronic decays from the actual NoFullyHadronicDecays samples, the inclusiveDecays samples and the NoFullyHadronicDecays samples with PDFWeights stored in them.
  // The weights given here have been calculated manually and are based on the average number of MC events across years and top/antitop samples.
  // For events with fully hadronic decays, we rely on the inclusiveDecays samples and let weight = 1.
  // The value "MC events / real expected event" is way more stable across years than for TT samples. Thus, we do not need to distinguish between years in order to not make things too complicated.
  else if(is_tW) {
    const auto & GENtW = event.get(fHandle_GENtW);
    if(is_tW_incl && !(GENtW.IsTopHadronicDecay() && GENtW.IsAssHadronicDecay())) weight = 0.16;
    else if(is_tW_nfhd && !is_tW_nfhd_syst && !is_tW_nfhd_PDF && !is_tW_nfhd_DS) weight = 0.41;
    // Need to make sure not to adjust the event weights in PDFWeights samples if we actually rely solely on the PDFWeights samples to get the tW distrubtions with varied PDFs.
    else if(is_tW_nfhd_PDF && !fDoingPDFVariations) weight = 0.43;
  }

  if(fApply) event.weight *= weight;
  event.set(fHandle_weight, weight);

  return true;
}

//____________________________________________________________________________________________________
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
METFilterSelection::METFilterSelection(Context & ctx): fYear(extract_year(ctx)) {

  fAndSel.reset(new AndSelection(ctx, "metfilters"));

  if(is_UL(fYear)) {
    fAndSel->add<TriggerSelection>("goodVertices", "Flag_goodVertices");
    fAndSel->add<TriggerSelection>("globalSuperTightHalo2016Filter", "Flag_globalSuperTightHalo2016Filter");
    fAndSel->add<TriggerSelection>("HBHENoiseFilter", "Flag_HBHENoiseFilter");
    fAndSel->add<TriggerSelection>("HBHENoiseIsoFilter", "Flag_HBHENoiseIsoFilter");
    fAndSel->add<TriggerSelection>("EcalDeadCellTriggerPrimitiveFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter");
    fAndSel->add<TriggerSelection>("BadPFMuonFilter", "Flag_BadPFMuonFilter");
    fAndSel->add<TriggerSelection>("BadPFMuonDzFilter", "Flag_BadPFMuonDzFilter"); // right now not available in data ntuples
    // fAndSel->add<TriggerSelection>("BadChargedCandidateFilter", "Flag_BadChargedCandidateFilter"); // currently not recommended
    fAndSel->add<TriggerSelection>("eeBadScFilter", "Flag_eeBadScFilter");
    if(fYear == Year::isUL17 || fYear == Year::isUL18) {
      fAndSel->add<TriggerSelection>("ecalBadCalibFilter", "Flag_ecalBadCalibFilter");
      // fAndSel->add<TriggerSelection>("hfNoisyHitsFilter", "Flag_hfNoisyHitsFilter"); // currently not recommended
    }
  }
  else {
    throw runtime_error("METFilterSelection: Non-UL years not implemented");
  }
}

bool METFilterSelection::passes(const Event & event) {

  if(!fAndSel->passes(event)) return false;

  return true;
}

//____________________________________________________________________________________________________
// https://twiki.cern.ch/twiki/bin/view/CMS/PileupJetIDUL
JetPUID::JetPUID(const wp & working_point): fWP(working_point) {}

bool JetPUID::operator()(const Jet & jet, const Event & event) const {
  const string year = event.year;
  const double eta = fabs(jet.v4().eta());
  const double pt = jet.v4().pt();
  if(eta > 5.0 || pt < 10 || pt > 50) return true;

  double x(-2.); // ranges between -1 and +1 (BDT discriminator)
  switch(fWP) {
    case WP_LOOSE :
      if(year == "UL16preVFP" || year == "UL16postVFP") {
        if(eta < 2.5) {
          if(pt < 20) x = -0.95;
          else if(pt < 30) x = -0.90;
          else if(pt < 40) x = -0.71;
          else x = -0.42;
        }
        else if(eta < 2.75) {
          if(pt < 20) x = -0.70;
          else if(pt < 30) x = -0.57;
          else if(pt < 40) x = -0.36;
          else x = -0.09;
        }
        else if(eta < 3.0) {
          if(pt < 20) x = -0.52;
          else if(pt < 30) x = -0.43;
          else if(pt < 40) x = -0.29;
          else x = -0.14;
        }
        else {
          if(pt < 20) x = -0.49;
          else if(pt < 30) x = -0.42;
          else if(pt < 40) x = -0.23;
          else x = -0.02;
        }
      }
      else if(year == "UL17" || year == "UL18") {
        if(eta < 2.5) {
          if(pt < 20) x = -0.95;
          else if(pt < 30) x = -0.88;
          else if(pt < 40) x = -0.63;
          else x = -0.19;
        }
        else if(eta < 2.75) {
          if(pt < 20) x = -0.72;
          else if(pt < 30) x = -0.55;
          else if(pt < 40) x = -0.18;
          else x = 0.22;
        }
        else if(eta < 3.0) {
          if(pt < 20) x = -0.68;
          else if(pt < 30) x = -0.60;
          else if(pt < 40) x = -0.43;
          else x = -0.13;
        }
        else {
          if(pt < 20) x = -0.47;
          else if(pt < 30) x = -0.43;
          else if(pt < 40) x = -0.24;
          else x = -0.03;
        }
      }
      else throw runtime_error((string)"JetPUID::operator()(): Year '"+year+"' not implemented");
      break;
    case WP_MEDIUM :
      if(year == "UL16preVFP" || year == "UL16postVFP") {
        if(eta < 2.5) {
          if(pt < 20) x = 0.20;
          else if(pt < 30) x = 0.62;
          else if(pt < 40) x = 0.86;
          else x = 0.93;
        }
        else if(eta < 2.75) {
          if(pt < 20) x = -0.56;
          else if(pt < 30) x = -0.39;
          else if(pt < 40) x = -0.10;
          else x = 0.19;
        }
        else if(eta < 3.0) {
          if(pt < 20) x = -0.43;
          else if(pt < 30) x = -0.32;
          else if(pt < 40) x = -0.15;
          else x = 0.04;
        }
        else {
          if(pt < 20) x = -0.38;
          else if(pt < 30) x = -0.29;
          else if(pt < 40) x = -0.08;
          else x = 0.12;
        }
      }
      else if(year == "UL17" || year == "UL18") {
        if(eta < 2.5) {
          if(pt < 20) x = 0.26;
          else if(pt < 30) x = 0.68;
          else if(pt < 40) x = 0.90;
          else x = 0.96;
        }
        else if(eta < 2.75) {
          if(pt < 20) x = -0.33;
          else if(pt < 30) x = -0.04;
          else if(pt < 40) x = 0.36;
          else x = 0.61;
        }
        else if(eta < 3.0) {
          if(pt < 20) x = -0.54;
          else if(pt < 30) x = -0.43;
          else if(pt < 40) x = -0.16;
          else x = 0.14;
        }
        else {
          if(pt < 20) x = -0.37;
          else if(pt < 30) x = -0.30;
          else if(pt < 40) x = -0.09;
          else x = 0.12;
        }
      }
      else throw runtime_error((string)"JetPUID::operator()(): Year '"+year+"' not implemented");
      break;
    case WP_TIGHT :
      if(year == "UL16preVFP" || year == "UL16postVFP") {
        if(eta < 2.5) {
          if(pt < 20) x = 0.71;
          else if(pt < 30) x = 0.87;
          else if(pt < 40) x = 0.94;
          else x = 0.97;
        }
        else if(eta < 2.75) {
          if(pt < 20) x = -0.32;
          else if(pt < 30) x = -0.08;
          else if(pt < 40) x = 0.24;
          else x = 0.48;
        }
        else if(eta < 3.0) {
          if(pt < 20) x = -0.30;
          else if(pt < 30) x = -0.16;
          else if(pt < 40) x = 0.05;
          else x = 0.26;
        }
        else {
          if(pt < 20) x = -0.22;
          else if(pt < 30) x = -0.12;
          else if(pt < 40) x = 0.10;
          else x = 0.29;
        }
      }
      else if(year == "UL17" || year == "UL18") {
        if(eta < 2.5) {
          if(pt < 20) x = 0.77;
          else if(pt < 30) x = 0.90;
          else if(pt < 40) x = 0.96;
          else x = 0.98;
        }
        else if(eta < 2.75) {
          if(pt < 20) x = 0.38;
          else if(pt < 30) x = 0.60;
          else if(pt < 40) x = 0.82;
          else x = 0.92;
        }
        else if(eta < 3.0) {
          if(pt < 20) x = -0.31;
          else if(pt < 30) x = -0.12;
          else if(pt < 40) x = 0.20;
          else x = 0.47;
        }
        else {
          if(pt < 20) x = -0.21;
          else if(pt < 30) x = -0.13;
          else if(pt < 40) x = 0.09;
          else x = 0.29;
        }
      }
      else throw runtime_error((string)"JetPUID::operator()(): Year '"+year+"' not implemented");
      break;
    default :
      throw runtime_error((string)"JetPUID::operator()(): Unknown working point");
  }

  return jet.pileupID() > x;
}

//____________________________________________________________________________________________________
const Jet * getCHSmatch(const Jet & puppijet, const Event & event, const Event::Handle<vector<Jet>> & h_chsjets, const bool safe) {
  const Jet *chsjet = match(puppijet, event.get(h_chsjets), kDeltaRForPuppiCHSMatch);
  if(safe && chsjet == nullptr) throw runtime_error("getCHSmatch(): CHS jet not found! Did you properly clean your PUPPI collection by calling 'MatchPuppiToCHSAndSetBTagHandles::process()'?");
  else return chsjet;
}

//____________________________________________________________________________________________________
MatchPuppiToCHSAndSetBTagHandles::MatchPuppiToCHSAndSetBTagHandles(Context & ctx, const BTag::algo & btag_algo, const BTag::wp & btag_wp):
  fHandle_PUPPIjets(ctx.get_handle<vector<Jet>>("jets")), // after this module, will contain all forward PUPPI jets and all CHS-matched central PUPPI jets
  fHandle_CHSjets(ctx.get_handle<vector<Jet>>(kCollectionName_AK4CHS)), // all CHS jets, unmodified
  fHandle_pairedPUPPIjets(ctx.get_handle<vector<Jet>>(kHandleName_pairedPUPPIjets)), // all central PUPPI jets matched to CHS jets
  fHandle_pairedCHSjets(ctx.get_handle<vector<Jet>>(kHandleName_pairedCHSjets)), // can be used as handle for b-tagging discriminator reweighting class; double-counted CHS jets not strictly ruled out but we'll ignore this
  fHandle_forwardPUPPIjets(ctx.get_handle<vector<Jet>>(kHandleName_forwardPUPPIjets)), // all forward PUPPI jets
  fHandle_uncleanedPUPPIjets(ctx.get_handle<vector<Jet>>(kHandleName_uncleanedPUPPIjets)), // all PUPPI jets no matter if matched to CHS jet or not, and no matter if central or forward
  fBTagWP(btag_wp),
  fBTagID_loose(BTag(btag_algo, BTag::WP_LOOSE)),
  fBTagID_medium(BTag(btag_algo, BTag::WP_MEDIUM)),
  fBTagID_tight(BTag(btag_algo, BTag::WP_TIGHT)),
  fHandle_bJets(ctx.get_handle<vector<Jet>>(kHandleName_bJets)),
  fHandle_bJets_loose(ctx.get_handle<vector<Jet>>(kHandleName_bJets_loose)),
  fHandle_bJets_medium(ctx.get_handle<vector<Jet>>(kHandleName_bJets_medium)),
  fHandle_bJets_tight(ctx.get_handle<vector<Jet>>(kHandleName_bJets_tight)),
  // Write number of jets to AnalysisTree:
  fHandle_n_jets(ctx.declare_event_output<int>(kHandleName_n_jets)),
  fHandle_n_jets_central(ctx.declare_event_output<int>(kHandleName_n_jets_central)),
  fHandle_n_jets_forward(ctx.declare_event_output<int>(kHandleName_n_jets_forward)),
  fHandle_n_bJets(ctx.declare_event_output<int>(kHandleName_n_bJets)),
  fHandle_n_bJets_loose(ctx.declare_event_output<int>(kHandleName_n_bJets_loose)),
  fHandle_n_bJets_medium(ctx.declare_event_output<int>(kHandleName_n_bJets_medium)),
  fHandle_n_bJets_tight(ctx.declare_event_output<int>(kHandleName_n_bJets_tight)),
  fHandle_PrimaryLepton(ctx.get_handle<FlavorParticle>(kHandleName_PrimaryLepton)),
  fHandle_pairedPUPPIjets_hemi(ctx.get_handle<vector<Jet>>(kHandleName_pairedPUPPIjets_hemi)),
  fHandle_pairedCHSjets_hemi(ctx.get_handle<vector<Jet>>(kHandleName_pairedCHSjets_hemi)),
  fHandle_bJets_hemi(ctx.get_handle<vector<Jet>>(kHandleName_bJets_hemi)),
  fHandle_bJets_hemi_loose(ctx.get_handle<vector<Jet>>(kHandleName_bJets_hemi_loose)),
  fHandle_bJets_hemi_medium(ctx.get_handle<vector<Jet>>(kHandleName_bJets_hemi_medium)),
  fHandle_bJets_hemi_tight(ctx.get_handle<vector<Jet>>(kHandleName_bJets_hemi_tight)),
  fHandle_n_bJets_hemi(ctx.declare_event_output<int>(kHandleName_n_bJets_hemi)),
  fHandle_n_bJets_hemi_loose(ctx.declare_event_output<int>(kHandleName_n_bJets_hemi_loose)),
  fHandle_n_bJets_hemi_medium(ctx.declare_event_output<int>(kHandleName_n_bJets_hemi_medium)),
  fHandle_n_bJets_hemi_tight(ctx.declare_event_output<int>(kHandleName_n_bJets_hemi_tight))
{}

bool MatchPuppiToCHSAndSetBTagHandles::process(Event & event) {
  FlavorParticle primlep;
  const bool valid_primlep = event.is_valid(fHandle_PrimaryLepton);
  if(valid_primlep) primlep = event.get(fHandle_PrimaryLepton);
  //__________________________________________________
  const vector<Jet> uncleaned_puppijets = event.get(fHandle_PUPPIjets);
  event.set(fHandle_uncleanedPUPPIjets, uncleaned_puppijets);
  //__________________________________________________
  // Clean all central PUPPI jets which don't have a match to a CHS jet; keep all forward PUPPI jets no matter if matched or not
  vector<Jet> cleaned_puppijets;
  vector<Jet> paired_puppijets;
  vector<Jet> paired_puppijets_hemi;
  vector<Jet> paired_chsjets;
  vector<Jet> paired_chsjets_hemi;
  vector<Jet> forward_puppijets;
  for(const Jet & puppijet : uncleaned_puppijets) {
    if(fabs(puppijet.v4().eta()) >= kAbsEtaBTagThreshold) {
      forward_puppijets.push_back(puppijet);
      cleaned_puppijets.push_back(puppijet);
    }
    else {
      const Jet *chsjetPtr = getCHSmatch(puppijet, event, fHandle_CHSjets, false);
      if(chsjetPtr != nullptr) {
        Jet chsjet = *chsjetPtr;
        if(fabs(chsjet.v4().eta()) >= kAbsEtaBTagThreshold) chsjet.set_eta(chsjet.v4().eta() > 0 ? kAbsEtaBTagThreshold-0.0001 : -kAbsEtaBTagThreshold+0.0001); // do this so that b-tagging SF can still be used if CHS eta lies outside SF eta range (but maybe SF = 1 is better?)
        paired_puppijets.push_back(puppijet);
        paired_chsjets.push_back(chsjet);
        cleaned_puppijets.push_back(puppijet);
        if(valid_primlep && deltaR(primlep.v4(), puppijet.v4()) < kDeltaRLeptonicHemisphere) {
          paired_puppijets_hemi.push_back(puppijet);
          paired_chsjets_hemi.push_back(chsjet);
        }
      }
    }
  }
  event.set(fHandle_PUPPIjets, cleaned_puppijets);
  event.set(fHandle_pairedPUPPIjets, paired_puppijets);
  event.set(fHandle_pairedPUPPIjets_hemi, paired_puppijets_hemi);
  event.set(fHandle_pairedCHSjets, paired_chsjets);
  event.set(fHandle_pairedCHSjets_hemi, paired_chsjets_hemi);
  event.set(fHandle_forwardPUPPIjets, forward_puppijets);

  event.set(fHandle_n_jets, cleaned_puppijets.size());
  event.set(fHandle_n_jets_central, paired_puppijets.size());
  event.set(fHandle_n_jets_forward, forward_puppijets.size());

  //__________________________________________________
  // Retrieve the b-tagging information from the matched CHS jets in order to find b-tagged PUPPI jets in central region
  // At this point, getCHSmatch will *always* find a CHS jet since the central PUPPI jet collection has already been cleaned accordingly (no PUPPI without CHS match, see above)
  vector<Jet> bjets_loose;
  vector<Jet> bjets_medium;
  vector<Jet> bjets_tight;
  vector<Jet> bjets_hemi_loose;
  vector<Jet> bjets_hemi_medium;
  vector<Jet> bjets_hemi_tight;
  for(const Jet & puppijet : paired_puppijets) {
    const Jet *chsjetPtr = getCHSmatch(puppijet, event, fHandle_CHSjets);
    if(fBTagID_loose(*chsjetPtr, event)) {
      bjets_loose.push_back(puppijet);
      if(valid_primlep && deltaR(primlep.v4(), puppijet.v4()) < kDeltaRLeptonicHemisphere) bjets_hemi_loose.push_back(puppijet);
    }
    if(fBTagID_medium(*chsjetPtr, event)) {
      bjets_medium.push_back(puppijet);
      if(valid_primlep && deltaR(primlep.v4(), puppijet.v4()) < kDeltaRLeptonicHemisphere) bjets_hemi_medium.push_back(puppijet);
    }
    if(fBTagID_tight(*chsjetPtr, event)) {
      bjets_tight.push_back(puppijet);
      if(valid_primlep && deltaR(primlep.v4(), puppijet.v4()) < kDeltaRLeptonicHemisphere) bjets_hemi_tight.push_back(puppijet);
    }
  }
  switch(fBTagWP) {
    case BTag::WP_LOOSE :
    event.set(fHandle_bJets, bjets_loose);
    event.set(fHandle_n_bJets, bjets_loose.size());
    event.set(fHandle_bJets_hemi, bjets_hemi_loose);
    event.set(fHandle_n_bJets_hemi, bjets_hemi_loose.size());
    break;
    case BTag::WP_MEDIUM :
    event.set(fHandle_bJets, bjets_medium);
    event.set(fHandle_n_bJets, bjets_medium.size());
    event.set(fHandle_bJets_hemi, bjets_hemi_medium);
    event.set(fHandle_n_bJets_hemi, bjets_hemi_medium.size());
    break;
    case BTag::WP_TIGHT :
    event.set(fHandle_bJets, bjets_tight);
    event.set(fHandle_n_bJets, bjets_tight.size());
    event.set(fHandle_bJets_hemi, bjets_hemi_tight);
    event.set(fHandle_n_bJets_hemi, bjets_hemi_tight.size());
    break;
    default :
    throw invalid_argument("MatchPuppiToCHSAndSetBTagHandles::process(): Invalid b-tagging WP!");
  }
  event.set(fHandle_bJets_loose, bjets_loose);
  event.set(fHandle_bJets_medium, bjets_medium);
  event.set(fHandle_bJets_tight, bjets_tight);

  event.set(fHandle_n_bJets_loose, bjets_loose.size());
  event.set(fHandle_n_bJets_medium, bjets_medium.size());
  event.set(fHandle_n_bJets_tight, bjets_tight.size());

  event.set(fHandle_bJets_hemi_loose, bjets_hemi_loose);
  event.set(fHandle_bJets_hemi_medium, bjets_hemi_medium);
  event.set(fHandle_bJets_hemi_tight, bjets_hemi_tight);

  event.set(fHandle_n_bJets_hemi_loose, bjets_hemi_loose.size());
  event.set(fHandle_n_bJets_hemi_medium, bjets_hemi_medium.size());
  event.set(fHandle_n_bJets_hemi_tight, bjets_hemi_tight.size());

  return true;
}

//____________________________________________________________________________________________________
ObjectPtSorter::ObjectPtSorter(Context & ctx, const bool do_fatjets):
  bDoFatJets(do_fatjets),
  fHandle_CHSjets(ctx.get_handle<vector<Jet>>(kCollectionName_AK4CHS)),
  fHandle_AK8Collection_rec(ctx.get_handle<vector<TopJet>>(kCollectionName_AK8_rec)),
  fHandle_AK8Collection_gen(ctx.get_handle<vector<GenTopJet>>(kCollectionName_AK8_gen))
{}

bool ObjectPtSorter::process(Event & event) {
  // Sorts AK4 PUPPI jets
  sort_by_pt<Jet>(*event.jets);
  // Sorts AK4 CHS jets
  vector<Jet> &ak4chsjets = event.get(fHandle_CHSjets);
  sort_by_pt<Jet>(ak4chsjets);
  // Sorts AK4 gen jets
  if(!event.isRealData) sort_by_pt<GenJet>(*event.genjets);

  if(bDoFatJets) {
    // Sorts HOTVR jets
    sort_by_pt<TopJet>(*event.topjets);
    for(auto & j : *event.topjets) {
      vector<Jet> subjets = j.subjets();
      sort_by_pt<Jet>(subjets);
      j.set_subjets(move(subjets));
    }
    // Sorts HOTVR gen jets
    if(!event.isRealData) {
      sort_by_pt<GenTopJet>(*event.gentopjets);
      // for(auto & j : *event.gentopjets) { // set_subjets not available for GenTopJet
      //   vector<GenJet> subjets = j.subjets();
      //   sort_by_pt<GenJet>(subjets);
      //   j.set_subjets(move(subjets));
      // }
    }
    // Sorts AK8 jets
    vector<TopJet> &ak8jets = event.get(fHandle_AK8Collection_rec);
    sort_by_pt<TopJet>(ak8jets);
    for(auto & j : ak8jets) {
      vector<Jet> subjets = j.subjets();
      sort_by_pt<Jet>(subjets);
      j.set_subjets(move(subjets));
    }
    // Sorts AK8 gen jets
    if(!event.isRealData) {
      vector<GenTopJet> &ak8genjets = event.get(fHandle_AK8Collection_gen);
      sort_by_pt<GenTopJet>(ak8genjets);
      // for(auto & j : ak8genjets) { // set_subjets not available for GenTopJet
      //   vector<GenJet> subjets = j.subjets();
      //   sort_by_pt<GenJet>(subjets);
      //   j.set_subjets(move(subjets));
      // }
    }
  }
  return true;
}

//____________________________________________________________________________________________________
BTagNJetScaleFactor::BTagNJetScaleFactor(Context & ctx):
 fYear(extract_year(ctx)),
 fChannel(extract_channel(ctx)),
 fHandle_weight_btag_njet_sf(ctx.declare_event_output<float>(kHandleName_weight_btag_njet_sf)),
 fHandle_pairedPUPPIjets(ctx.get_handle<vector<Jet>>(kHandleName_pairedPUPPIjets))
{
  for(const auto & proc : kMCProcess_toString) {
    if(ctx.get("dataset_version").find(proc.second) == 0) {
      fMCProcess = proc.first;
    }
  }
  if(fMCProcess == MCProcess::other) return;

  TFile *file = TFile::Open(ctx.get("BTagSFNJetReweightFile").c_str(), "READ");
  for(int njets = 1; njets <= 8; njets++) {
    string sf_name = "sf_";
    if(fYear == Year::isUL16preVFP) sf_name += "UL16preVFP";
    else if(fYear == Year::isUL16postVFP) sf_name += "UL16postVFP";
    else if(fYear == Year::isUL17) sf_name += "UL17";
    else if(fYear == Year::isUL18) sf_name += "UL18";
    if(fChannel == Channel::isEle) sf_name += "_ele";
    else if(fChannel == Channel::isMuo) sf_name += "_muo";
    sf_name += "_"+kMCProcess_toString.at(fMCProcess);
    sf_name += "_njets"+to_string(njets);
    TVectorF *sf_vec = (TVectorF*)file->Get(sf_name.c_str());
    fSFMap[fYear][fChannel][fMCProcess][njets] = (*sf_vec)[0];
  }
}

bool BTagNJetScaleFactor::process(Event & event) {
  float weight = 1.f;
  if(fMCProcess == MCProcess::other || event.isRealData) {
    event.set(fHandle_weight_btag_njet_sf, weight);
    return true;
  }
  int njets = max(1, min(8, (int)event.get(fHandle_pairedPUPPIjets).size()));
  weight = fSFMap.at(fYear).at(fChannel).at(fMCProcess).at(njets);
  event.weight *= weight;
  event.set(fHandle_weight_btag_njet_sf, weight);
  return true;
}

}}
