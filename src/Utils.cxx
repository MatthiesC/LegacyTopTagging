#include "UHH2/common/include/Utils.h"

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

}}
