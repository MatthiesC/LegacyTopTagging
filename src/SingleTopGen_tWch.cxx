// based on TTbarGen.cxx
// edited by Christopher Matthies, 02/2018

#include "UHH2/LegacyTopTagging/include/SingleTopGen_tWch.h"
// #include "UHH2/core/include/Utils.h"

using namespace std;
using namespace uhh2;
using namespace ltt;

namespace uhh2 { namespace ltt {

SingleTopGen_tWch::SingleTopGen_tWch(const vector<GenParticle> & genparticles, bool throw_on_failure): m_type(e_notfound) {

  int n_top = 0;
  int n_WAss = 0;
  GenParticle top, WAss, initial1, initial2;

  for(unsigned int i = 0; i < genparticles.size(); ++i) {
    const GenParticle & genp = genparticles[i];
    if (genp.index() == 0){ // find the initial state particles
      initial1 = genp;
      m_initial1 = initial1;
    }
    else if (genp.index() == 1){
      initial2 = genp;
      m_initial2 = initial2;
    }
    else if (abs(genp.pdgId()) == 6){ // We don't distinguish between top and antitop. May be changed later but for now it's okay
      top = genp;
      auto w = genp.daughter(&genparticles, 1);
      auto b = genp.daughter(&genparticles, 2);
      if(!w || !b) {
        if(throw_on_failure) throw runtime_error("SingleTopGen_tWch: top has not ==2 daughters");
        return;
      }
      if (abs(w->pdgId()) != 24) {
        std::swap(w, b);
      }
      /* It rarely happens that the list of genparticles contains 4 or more (instead of 2) particles which reckon the same top
         as their mother although each particle including the tops can just have two daughters. E.g. if the top emits a photon
         before decaying into b and W, this photon may split up into two leptons which reckon the top as their mother, too.
         Therefore, it may happen that those leptons are considered as the top daughters whereas b and W are "ignored" and cannot
         be found. This workaround fixes that issue: */
      if(abs(w->pdgId()) != 24) {
        for(unsigned int j = 0; j < genparticles.size(); ++j) {
          const GenParticle & gp = genparticles[j];
          auto m1 = gp.mother(&genparticles, 1);
          auto m2 = gp.mother(&genparticles, 2);
          bool has_top_mother = ((m1 && m1->index() == genp.index()) || (m2 && m2->index() == genp.index()));
          if(has_top_mother && (abs(gp.pdgId()) == 24)) {
            w = &gp;
            break;
          }
        }
      }
      if (abs(w->pdgId()) != 24) {
        if(throw_on_failure) throw runtime_error("SingleTopGen_tWch: top has no W daughter");
        return;
      }

      // NOTE: here, we could skip over intermediate W bosons. However,
      // this Pythia8-related problem is now fixed when creating ntuples already,
      // so this should not be necessary.

      /* Do a similar workaround as above if the expected b daughter has not been found yet */
      if(abs(b->pdgId()) != 5 && abs(b->pdgId()) != 3 && abs(b->pdgId()) != 1) {
        for(unsigned int j = 0; j < genparticles.size(); ++j) {
          const GenParticle & gp = genparticles[j];
          auto m1 = gp.mother(&genparticles, 1);
          auto m2 = gp.mother(&genparticles, 2);
          bool has_top_mother = ((m1 && m1->index() == genp.index()) || (m2 && m2->index() == genp.index()));
          if(has_top_mother && (abs(gp.pdgId()) == 5 || abs(gp.pdgId()) == 3 || abs(gp.pdgId()) == 1)) {
            b = &gp;
            break;
          }
        }
      }
      if (abs(b->pdgId()) != 5 && abs(b->pdgId()) != 3   && abs(b->pdgId()) != 1) {
        if(throw_on_failure) throw runtime_error("SingleTopGen_tWch: top has no b daughter");
        return;
      }
      // now get WTop daughters:

      int n_WTopDau = 0;

      auto wd1 = w->daughter(&genparticles, 1);
      auto wd2 = w->daughter(&genparticles, 2);

      while (n_WTopDau != 2){

        if(wd1 && !wd2){
          w = wd1;
          wd1 = w->daughter(&genparticles, 1);
          wd2 = w->daughter(&genparticles, 2);
        }
        else if(wd1 && wd2){
          n_WTopDau = 2;
        }

        else{
          if(throw_on_failure) throw runtime_error("SingleTopGen_tWch: WTop has no daughters");
          return;
        }

      }
      if(!wd1 || !wd2){
        if(throw_on_failure) throw runtime_error("SingleTopGen_tWch: WTop has not ==2 daughters");
        return;
      }

      // now that we collected everything, fill the member variables.
      // Unlike in TTbarGen.cxx, where we use different member variables according to top quark charge, we do not distinguish them here because there is just one top
      m_Top = top;
      m_WTop = *w;
      m_bTop = *b;
      m_WTopDecay1 = *wd1;
      m_WTopDecay2 = *wd2;
      ++n_top;

    }

    // Now get the associated W boson of the tW channel and its daughters

    else if (abs(genp.pdgId()) == 24){ // We don't distinguish between W+ and W- associated (depends on the presence of top resp. anti-top)
      WAss = genp;
      auto m1 = genp.mother(&genparticles, 1);
      auto m2 = genp.mother(&genparticles, 2);
      if (!(m1 && m2 && m1->index() + m2->index() == 1)) { // check if the indices of both mothers are 0+1==1, meaning this W is the associated one and not the top daughter
	continue;
      }

      int n_WAssDau = 0;

      auto wassd1 = genp.daughter(&genparticles, 1);
      auto wassd2 = genp.daughter(&genparticles, 2);

      while (n_WAssDau != 2){

	if(wassd1 && !wassd2){
	  WAss = *wassd1;
	  wassd1 = WAss.daughter(&genparticles, 1);
	  wassd2 = WAss.daughter(&genparticles, 2);
	}
	else if(wassd1 && wassd2){
	  n_WAssDau = 2;
	}
	else{
	  if(throw_on_failure) throw runtime_error("SingleTopGen_tWch: WAss has no daughters");
	  return;
	}

      }

      if(!wassd1 || !wassd2){
	if(throw_on_failure) throw runtime_error("SingleTopGen_tWch: WAss has not ==2 daughters");
	return;
      }

      // fill member variables
      m_WAss = WAss;
      m_WAssDecay1 = *wassd1;
      m_WAssDecay2 = *wassd2;
      ++n_WAss;

    }
    // find possible final state gluon if NLO process
    else if (genp.pdgId() == 21 && genp.index() > 1){
      GenParticle final_gluon;
      final_gluon = genp;
      m_gluonAss = final_gluon;
    }
    // find possible final state bottom (not from top decay, but from non-LO processes)
    else if (abs(genp.pdgId()) == 5 && genp.mother(&genparticles, 1)->index() < 2 && genp.mother(&genparticles, 2)->index() < 2){
      GenParticle bAss;
      bAss = genp;
      m_bAss = bAss;
      m_has_bAss = true;
    }
  }



  if(n_top != 1){
    if(throw_on_failure)  throw runtime_error("SingleTopGen_tWch: did not find exactly one (anti)top in the event");
    return;
  }

  if(n_WAss != 1){
    if(throw_on_failure)  throw runtime_error("SingleTopGen_tWch: did not find exactly one associated (!) W in the event");
    return;
  }



  // calculate decay channel by counting the number of charged leptons
  // in the WTop and WAss daughters:
  int n_top_e = 0, n_top_m = 0, n_top_t = 0;
  for(const auto & wd : {m_WTopDecay1, m_WTopDecay2}){
    int id_top = abs(wd.pdgId());
    if(id_top == 11) ++n_top_e;
    else if(id_top == 13) ++n_top_m;
    else if(id_top == 15) ++n_top_t;
  }
  int n_ass_e = 0, n_ass_m = 0, n_ass_t = 0;
  for(const auto & wd : {m_WAssDecay1, m_WAssDecay2}) {
    int id_ass = abs(wd.pdgId());
    if(id_ass == 11) ++n_ass_e;
    else if(id_ass == 13) ++n_ass_m;
    else if(id_ass == 15) ++n_ass_t;
  }
  int n_ass_lep = n_ass_e + n_ass_m + n_ass_t;
  int n_top_lep = n_top_e + n_top_m + n_top_t;

  // fully leptonic or semi-leptonic (lepton from WAss):
  if(n_ass_lep == 1){
    if(n_ass_e == 1){
      if(n_top_e == 1) m_type = e_assele_topele;
      else if(n_top_m == 1) m_type = e_assele_topmuo;
      else if(n_top_t == 1) m_type = e_assele_toptau;
      else if(n_top_lep == 0) m_type = e_assele_tophad;
    }
    else if(n_ass_m == 1){
      if(n_top_e == 1) m_type = e_assmuo_topele;
      else if(n_top_m == 1) m_type = e_assmuo_topmuo;
      else if(n_top_t == 1) m_type = e_assmuo_toptau;
      else if(n_top_lep == 0) m_type = e_assmuo_tophad;
    }
    else if(n_ass_t == 1){
      if(n_top_e == 1) m_type = e_asstau_topele;
      else if(n_top_m == 1) m_type = e_asstau_topmuo;
      else if(n_top_t == 1) m_type = e_asstau_toptau;
      else if(n_top_lep == 0) m_type = e_asstau_tophad;
    }
  }
  // semi-leptonic (lepton from WTop) or fully hadronic:
  else if(n_ass_lep == 0){
    if(n_top_e == 1) m_type = e_asshad_topele;
    else if(n_top_m == 1) m_type = e_asshad_topmuo;
    else if(n_top_t == 1) m_type = e_asshad_toptau;
    else if(n_top_lep == 0) m_type = e_asshad_tophad;
  }

}


GenParticle SingleTopGen_tWch::Top() const{
  return m_Top;
}

GenParticle SingleTopGen_tWch::WTop() const{
  return m_WTop;
}

GenParticle SingleTopGen_tWch::bTop() const{
  return m_bTop;
}

GenParticle SingleTopGen_tWch::WTopDecay1() const{
  return m_WTopDecay1;
}

GenParticle SingleTopGen_tWch::WTopDecay2() const{
  return m_WTopDecay2;
}

GenParticle SingleTopGen_tWch::WAss() const{
  return m_WAss;
}

GenParticle SingleTopGen_tWch::WAssDecay1() const{
  return m_WAssDecay1;
}

GenParticle SingleTopGen_tWch::WAssDecay2() const{
  return m_WAssDecay2;
}

GenParticle SingleTopGen_tWch::LeptonicW() const{
  if(IsAssLeptonicDecay() && IsTopHadronicDecay()) return m_WAss;
  else if(IsAssHadronicDecay() && IsTopLeptonicDecay()) return m_WTop;
  else throw runtime_error("SingleTopGen_tWch::LeptonicW(): This is not a l+jets event!");
}

GenParticle SingleTopGen_tWch::SingleLepton() const{
  if(IsAssLeptonicDecay() && IsTopHadronicDecay()) {
    if(abs(m_WAssDecay1.pdgId()) == 11 || abs(m_WAssDecay1.pdgId()) == 13 || abs(m_WAssDecay1.pdgId()) == 15) return m_WAssDecay1;
    else return m_WAssDecay2;
  }
  else if(IsAssHadronicDecay() && IsTopLeptonicDecay()) {
    if(abs(m_WTopDecay1.pdgId()) == 11 || abs(m_WTopDecay1.pdgId()) == 13 || abs(m_WTopDecay1.pdgId()) == 15) return m_WTopDecay1;
    else return m_WTopDecay2;
  }
  else throw runtime_error("SingleTopGen_tWch::SingleLepton(): This is not a l+jets event!");
}

GenParticle SingleTopGen_tWch::SingleNeutrino() const{
  if(IsAssLeptonicDecay() && IsTopHadronicDecay()) {
    if(abs(m_WAssDecay1.pdgId()) == 11 || abs(m_WAssDecay1.pdgId()) == 13 || abs(m_WAssDecay1.pdgId()) == 15) return m_WAssDecay2;
    else return m_WAssDecay1;
  }
  else if(IsAssHadronicDecay() && IsTopLeptonicDecay()) {
    if(abs(m_WTopDecay1.pdgId()) == 11 || abs(m_WTopDecay1.pdgId()) == 13 || abs(m_WTopDecay1.pdgId()) == 15) return m_WTopDecay2;
    else return m_WTopDecay1;
  }
  else throw runtime_error("SingleTopGen_tWch::SingleNeutrino(): This is not a l+jets event!");
}

GenParticle SingleTopGen_tWch::LeptAss() const{
  if(!(IsAssLeptonicDecay())) throw runtime_error("SingleTopGen_tWch: funciton LeptAss() called but this is no event with Wass decaying leptonically");
  if(abs(m_WAssDecay1.pdgId()) == 11 || abs(m_WAssDecay1.pdgId()) == 13 || abs(m_WAssDecay1.pdgId()) == 15) return m_WAssDecay1;
  else return m_WAssDecay2;
}

GenParticle SingleTopGen_tWch::NuAss() const{
  if(!(IsAssLeptonicDecay())) throw runtime_error("SingleTopGen_tWch: funciton NuAss() called but this is no event with Wass decaying leptonically");
  if(abs(m_WAssDecay1.pdgId()) == 11 || abs(m_WAssDecay1.pdgId()) == 13 || abs(m_WAssDecay1.pdgId()) == 15) return m_WAssDecay2;
  else return m_WAssDecay1;
}

GenParticle SingleTopGen_tWch::Initial1() const{
  return m_initial1;
}

GenParticle SingleTopGen_tWch::Initial2() const{
  return m_initial2;
}

GenParticle SingleTopGen_tWch::InitialBottomQuark() const{
  if(!(IsBottomGluonProcess())) throw runtime_error("SingleTopGen_tWch: function InitialBottomQuark() called but this is no bottom-gluon event");
  if(m_initial1.pdgId() == 21) return m_initial2;
  else return m_initial1;
}

GenParticle SingleTopGen_tWch::gluonAss() const{
  return m_gluonAss;
}

GenParticle SingleTopGen_tWch::bAss() const{
  return m_bAss;
}



SingleTopGen_tWch::E_DecayChannel SingleTopGen_tWch::DecayChannel() const{
  return m_type;
}



bool SingleTopGen_tWch::IsGluonGluonProcess() const{
  return m_initial1.pdgId() == 21 && m_initial2.pdgId() == 21;
}

bool SingleTopGen_tWch::IsBottomGluonProcess() const{
  bool initial1_is_bottom = (abs(m_initial1.pdgId()) == 5 || abs(m_initial1.pdgId()) == 3 || abs(m_initial1.pdgId()) == 1);
  bool initial2_is_bottom = (abs(m_initial2.pdgId()) == 5 || abs(m_initial2.pdgId()) == 3 || abs(m_initial2.pdgId()) == 1);
  return (m_initial1.pdgId() == 21 && initial2_is_bottom) || (initial1_is_bottom && m_initial2.pdgId() == 21);
}

bool SingleTopGen_tWch::HasAssociatedBottom() const{
  return m_has_bAss;
}

bool SingleTopGen_tWch::HasAssociatedGluon() const{
  return m_gluonAss.pdgId();
}



bool SingleTopGen_tWch::IsTopHadronicDecay() const{
  return abs(m_WTopDecay1.pdgId()) <= 5;
}

bool SingleTopGen_tWch::IsTopLeptonicDecay(const bool consider_taus) const{
  if(consider_taus) return IsTopToElectronDecay() || IsTopToMuonDecay() || IsTopToTauonDecay();
  else return IsTopToElectronDecay() || IsTopToMuonDecay();
}

bool SingleTopGen_tWch::IsTopToElectronDecay() const{
  return m_type == e_assele_topele || m_type == e_assmuo_topele || m_type == e_asstau_topele || m_type == e_asshad_topele;
}

bool SingleTopGen_tWch::IsTopToMuonDecay() const{
  return m_type == e_assele_topmuo || m_type == e_assmuo_topmuo || m_type == e_asstau_topmuo || m_type == e_asshad_topmuo;
}

bool SingleTopGen_tWch::IsTopToTauonDecay() const{
  return m_type == e_assele_toptau || m_type == e_assmuo_toptau || m_type == e_asstau_toptau || m_type == e_asshad_toptau;
}



bool SingleTopGen_tWch::IsAssHadronicDecay() const{
  return abs(m_WAssDecay1.pdgId()) <= 5;
}

bool SingleTopGen_tWch::IsAssLeptonicDecay(const bool consider_taus) const{
  if(consider_taus) return IsAssToElectronDecay() || IsAssToMuonDecay() || IsAssToTauonDecay();
  else return IsAssToElectronDecay() || IsAssToMuonDecay();
}

bool SingleTopGen_tWch::IsAssToElectronDecay() const{
  return m_type == e_assele_topele || m_type == e_assele_topmuo || m_type == e_assele_toptau || m_type == e_assele_tophad;
}

bool SingleTopGen_tWch::IsAssToMuonDecay() const{
  return m_type == e_assmuo_topele || m_type == e_assmuo_topmuo || m_type == e_assmuo_toptau || m_type == e_assmuo_tophad;
}

bool SingleTopGen_tWch::IsAssToTauonDecay() const{
  return m_type == e_asstau_topele || m_type == e_asstau_topmuo || m_type == e_asstau_toptau || m_type == e_asstau_tophad;
}


bool SingleTopGen_tWch::IsSemiLeptonic(const bool consider_taus) const{
  if(
    ( IsTopHadronicDecay() && IsAssLeptonicDecay(consider_taus) ) ||
    ( IsTopLeptonicDecay(consider_taus) && IsAssHadronicDecay() )
  ) return true;
  else return false;
}


SingleTopGen_tWchProducer::SingleTopGen_tWchProducer(uhh2::Context & ctx, const std::string & name, bool throw_on_failure_): throw_on_failure(throw_on_failure_){
  h_singletopgen_twch = ctx.get_handle<SingleTopGen_tWch>(name);
}

bool SingleTopGen_tWchProducer::process(Event & event){
  event.set(h_singletopgen_twch, SingleTopGen_tWch(*event.genparticles, throw_on_failure));
  return true;
}

}}
