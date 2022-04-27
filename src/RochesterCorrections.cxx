#include "UHH2/LegacyTopTagging/include/RochesterCorrections.h"

#include "TRandom3.h"

using namespace std;
using namespace uhh2;
using namespace ltt;


namespace uhh2 { namespace ltt {

RochesterCorrections::RochesterCorrections(Context & ctx):
  fYear(extract_year(ctx))
{
  cout << "Hello World from RochesterCorrections!" << endl;

  string filename = (string)getenv("CMSSW_BASE")+"/src/UHH2/LegacyTopTagging/data/rochester/";
  switch(fYear) {
    case Year::is2016v2 :
    filename += "RoccoR2016.txt";
    break;
    case Year::is2016v3 :
    filename += "RoccoR2016.txt";
    break;
    case Year::is2017v1 :
    filename += "RoccoR2017.txt";
    break;
    case Year::is2017v2 :
    filename += "RoccoR2017.txt";
    break;
    case Year::is2018 :
    filename += "RoccoR2018.txt";
    break;
    case Year::isUL16preVFP :
    filename += "RoccoR2016aUL.txt";
    break;
    case Year::isUL16postVFP :
    filename += "RoccoR2016bUL.txt";
    break;
    case Year::isUL17 :
    filename += "RoccoR2017UL.txt";
    break;
    case Year::isUL18 :
    filename += "RoccoR2018UL.txt";
    break;
  }
  rc.reset(new RoccoR());
  rc->init(filename);

  cout << "Caveat: Systematic variations not implemented yet" << endl;
  variation_set = 0;
  variation_member = 0;
}

bool RochesterCorrections::process(Event & event) {
  // procedure similar to JER smearing
  vector<GenParticle> gen_muons;
  if(!event.isRealData) {
    for(const GenParticle & gp : *event.genparticles) {
      if(abs(gp.pdgId()) == 13) gen_muons.push_back(gp);
    }
  }
  for(Muon & reco_muon : *event.muons) {
    double sf(1.);
    if(event.isRealData) {
      sf = rc->kScaleDT((int)reco_muon.charge(), reco_muon.v4().pt(), reco_muon.v4().eta(), reco_muon.v4().phi(), variation_set, variation_member);
    }
    else {
      const GenParticle *closest_gen_muon = closestParticle(reco_muon, gen_muons);
      // scaling method:
      if(closest_gen_muon != nullptr && deltaR(reco_muon, *closest_gen_muon) < 0.1 && closest_gen_muon->v4().pt() > 15.) { // deltaR and pt threshold are educated guesses done by myself
        sf = rc->kSpreadMC((int)reco_muon.charge(), reco_muon.v4().pt(), reco_muon.v4().eta(), reco_muon.v4().phi(), closest_gen_muon->v4().pt(), variation_set, variation_member);
      }
      // stochastic method:
      else {
        // Use random number generator with eta-dependent seed for reproducibility
        TRandom3 *random = new TRandom3((int)(fabs(reco_muon.v4().eta()*1000))); // TRandom3 is the same generator as used for gRandom
        sf = rc->kSmearMC((int)reco_muon.charge(), reco_muon.v4().pt(), reco_muon.v4().eta(), reco_muon.v4().phi(), reco_muon.innerTrack_trackerLayersWithMeasurement(), random->Rndm(), variation_set, variation_member);
      }
    }
    const LorentzVector muon_v4_before = reco_muon.v4();
    const LorentzVector muon_v4_after = muon_v4_before * sf;
    reco_muon.set_v4(muon_v4_after);
    // Change of use of unused "ptRatio" member: Use it as a kind of "JEC_factor_raw"
    reco_muon.set_ptRatio(-muon_v4_before.pt() / muon_v4_after.pt());
  }
  return true;
}

}}
