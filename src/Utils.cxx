#include "UHH2/common/include/Utils.h"

#include "UHH2/LegacyTopTagging/include/Utils.h"

using namespace std;
using namespace uhh2;


namespace ltt {

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

double HOTVR_mpair(const TopJet & topjet) {
  vector<Jet> subjets = topjet.subjets();
  if(subjets.size() < 3) throw runtime_error("HOTVR jet has less than 3 subjets, cannot calculate minimum pairwise mass.");
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

}
