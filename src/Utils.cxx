#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/MCWeight.h"

#include "UHH2/LegacyTopTagging/include/Utils.h"

using namespace std;
using namespace uhh2;
using namespace ltt;


namespace uhh2 { namespace ltt {

double tau32(const TopJet & topjet) {
  return min((double)(topjet.tau3() / topjet.tau2()), 0.99999);
}

double tau32groomed(const TopJet & topjet) {
  return min((double)(topjet.tau3_groomed() / topjet.tau2_groomed()), 0.99999);
}

double mSD(const TopJet & topjet) {
  LorentzVector subjet_sum(0,0,0,0);
  for(auto subjet : topjet.subjets()) {
      subjet_sum += subjet.v4();
  }
  return inv_mass_safe(subjet_sum);
}

double maxDeepCSVSubJetValue(const TopJet & topjet) {
  double result(-99999.);
  for(auto subjet : topjet.subjets()) {
    if(subjet.btag_DeepCSV() > result) result = subjet.btag_DeepCSV();
  }
  return min(result, 0.99999);
}

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

double HOTVR_fpt(const TopJet & topjet, const unsigned int subjet_i) {
  vector<Jet> subjets = topjet.subjets();
  sort_by_pt(subjets);
  return subjets.at(subjet_i).v4().Pt() / topjet.v4().Pt();
}

const TopJet * nextTopJet(const Particle & p, const vector<TopJet> & topjets) {
  return closestParticle(p, topjets);
}

METSelection::METSelection(const double _met_min, const double _met_max): met_min(_met_min), met_max(_met_max) {}

bool METSelection::passes(const Event & event) {
  bool passed_lower_limit = event.met->v4().Pt() > met_min;
  bool passed_upper_limit = event.met->v4().Pt() < met_max;
  return passed_lower_limit && passed_upper_limit;
}

PTWSelection::PTWSelection(Context & ctx, const double _ptw_min, const double _ptw_max): ptw_min(_ptw_min), ptw_max(_ptw_max), h_primlep(ctx.get_handle<FlavorParticle>("PrimaryLepton")) {}

bool PTWSelection::passes(const Event & event) {
  const FlavorParticle & primlep = event.get(h_primlep);
  LorentzVector v4 = event.met->v4() + primlep.v4();
  bool passed_lower_limit = v4.Pt() > ptw_min;
  bool passed_upper_limit = v4.Pt() < ptw_max;
  return passed_lower_limit && passed_upper_limit;
}

TwoDSelection::TwoDSelection(Context & ctx, const double _ptrel_min, const double _dr_min): ptrel_min(_ptrel_min), dr_min(_dr_min), h_primlep(ctx.get_handle<FlavorParticle>("PrimaryLepton")) {}

bool TwoDSelection::passes(const Event & event) {
  const FlavorParticle & primlep = event.get(h_primlep);
  const Jet *nextjet = nextJet(primlep, *event.jets);
  bool passed_ptrel_cut = pTrel(primlep, nextjet) > ptrel_min;
  bool passed_dr_cut = deltaR(primlep.v4(), nextjet->v4()) > dr_min;
  return passed_ptrel_cut || passed_dr_cut;
}

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

// Twiki: https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL201{6,7,8}
// SF files: https://gitlab.cern.ch/cms-muonPOG/muonefficiencies/-/tree/master/Run2/UL
MuonScaleFactors::MuonScaleFactors(Context & ctx) {
  const string filepath = ctx.get("MuonIDScaleFactorFile");
  const string direction = ctx.get("SystDirection_MuonID", "nominal");
  const string weight_name = "id_tight";
  sf_id.reset(new YearSwitcher(ctx));
  sf_id->setupUL16preVFP (make_shared<MCMuonScaleFactor>(ctx, filepath, "NUM_TightID_DEN_TrackerMuons_abseta_pt", 0.0, weight_name, false, direction));
  sf_id->setupUL16postVFP(make_shared<MCMuonScaleFactor>(ctx, filepath, "NUM_TightID_DEN_TrackerMuons_abseta_pt", 0.0, weight_name, false, direction));
  sf_id->setupUL17       (make_shared<MCMuonScaleFactor>(ctx, filepath, "NUM_TightID_DEN_TrackerMuons_abseta_pt", 0.0, weight_name, false, direction));
  sf_id->setupUL18       (make_shared<MCMuonScaleFactor>(ctx, filepath, "NUM_TightID_DEN_TrackerMuons_abseta_pt", 0.0, weight_name, false, direction));
  // No isolation scale factors applied since we do not use PF muon isolation but our custom 2D cut!
  // sf_iso.reset(new YearSwitcher(ctx));
  // ...
}

bool MuonScaleFactors::process(Event & event) {
  sf_id->process(event);
  // sf_iso->process(event);
  return true;
}

// Twiki: https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL201{6,7,8}
// SF files: https://gitlab.cern.ch/cms-muonPOG/muonefficiencies/-/tree/master/Run2/UL
TriggerScaleFactors::TriggerScaleFactors(Context & ctx) {
  const string filepath = ctx.get("TriggerScaleFactorFile");
  const string direction = ctx.get("SystDirection_MuonTrigger", "nominal");
  const string weight_name = "trigger";
  sf_trig.reset(new YearSwitcher(ctx));
  sf_trig->setupUL16preVFP (make_shared<MCMuonScaleFactor>(ctx, filepath, "NUM_Mu50_or_TkMu50_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose_abseta_pt", 0.0, weight_name, false, direction));
  sf_trig->setupUL16postVFP(make_shared<MCMuonScaleFactor>(ctx, filepath, "NUM_Mu50_or_TkMu50_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose_abseta_pt", 0.0, weight_name, false, direction));
  sf_trig->setupUL17       (make_shared<MCMuonScaleFactor>(ctx, filepath, "NUM_Mu50_or_OldMu100_or_TkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose_abseta_pt", 0.0, weight_name, false, direction));
  sf_trig->setupUL18       (make_shared<MCMuonScaleFactor>(ctx, filepath, "NUM_Mu50_or_OldMu100_or_TkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose_abseta_pt", 0.0, weight_name, false, direction));
}

bool TriggerScaleFactors::process(Event & event) {
  sf_trig->process(event);
  return true;
}

ProbeJetHandleSetter::ProbeJetHandleSetter(Context & ctx, const string & jetalgo_name, const string & coll_rec):
  h_probejet(ctx.get_handle<TopJet>("ProbeJet"+jetalgo_name)),
  h_topjets(ctx.get_handle<vector<TopJet>>(coll_rec.empty() ? "topjets" : coll_rec)) {}

bool ProbeJetHandleSetter::process(Event & event) {
  vector<TopJet> topjets = event.get(h_topjets);
  sort_by_pt(topjets);
  event.set(h_probejet, topjets.at(0)); // take the leading jet
  return true;
}

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
    GenParticle w_from_top = *(top->daughter(event.genparticles, 1));
    GenParticle b_from_top = *(top->daughter(event.genparticles, 2));
    if(abs(w_from_top.pdgId()) != 24) {
      swap(w_from_top, b_from_top);
    }
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

MergeScenarioHandleSetter::MergeScenarioHandleSetter(Context & ctx, const ProbeJetAlgo & _algo): algo(_algo) {
  string algo_name = "";
  if(algo == ProbeJetAlgo::isHOTVR) algo_name = "HOTVR";
  else if(algo == ProbeJetAlgo::isAK8) algo_name = "AK8";
  else throw runtime_error("MergeScenarioHandleSetter: Requested probejet algorithm not yet implemented in this class!");
  h_probejet = ctx.get_handle<TopJet>("ProbeJet"+algo_name);
  h_hadronictop = ctx.get_handle<GenParticle>("HadronicTopQuark"); // will be unset if process is neither ttbar->l+jets nor single t->hadronic

  output_has_probejet = ctx.declare_event_output<bool>("output_has_probejet_"+algo_name);
  output_merge_scenario = ctx.declare_event_output<MergeScenario>("output_merge_scenario_"+algo_name);
}

bool MergeScenarioHandleSetter::process(Event & event) {
  if(event.is_valid(h_probejet)) event.set(output_has_probejet, true);
  else event.set(output_has_probejet, false);

  MergeScenario msc = MergeScenario::isBackgroundMCorData;
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
  GenParticle gen_b = *(hadronic_top.daughter(event.genparticles, 1));
  GenParticle gen_w = *(hadronic_top.daughter(event.genparticles, 2));
  if(abs(gen_w.pdgId()) != 24) swap(gen_b, gen_w);
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

}}
