// based on TTbarGen.h
// edited by Christopher Matthies, 02/2018

#pragma once

#include "UHH2/core/include/GenParticle.h"
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include <vector>


namespace uhh2 { namespace ltt {

/** \brief ttbar generator truth information
 *
 * Interpret a vector of GenParticle as ttbar event, providing easier access to
 * the various components in the ttbar decay chain.
 *
 * The class can either be used directly by passing the genparticles vector;
 * another option is to use the SingleTopGen_tWchProducer (see below), which writes the
 * result to the event container, where it can be accessed later.
 */
class SingleTopGen_tWch {
public:

  /** construct from genparticles.
   *
   * The event should be an actual ttbar event, i.e.:
   *   - there are exactly two tops: one top and one anti-top
   *   - each top must have exactl two daughters
   *   - one of the top daughters must be a W
   *   - the Ws must have exactly two daughters
   *
   * Note that it is allowed that more particles than those belonging to the ttbar system
   * are in genparts; those are ignored.
   *
   * In case one of these conditions is not fulfilled, the behavior
   * depends on the throw_on_failure parameter: if it is true (the default), a runtime_error
   * is thrown with an explanation what failed. If it is false, no exception is thrown, but
   * not all contents is valid; most will return a 0 vector. The one thing guaranteed is that the
   * decaychannel will be e_notfound. If using throw_on_failure = false, it is thus a good idea
   * to check the decaychannel.
   */
  explicit SingleTopGen_tWch(const std::vector<GenParticle> & genparts, bool throw_on_failure = true);

  enum E_DecayChannel{
    e_assele_topele,
    e_assele_topmuo,
    e_assele_toptau,
    e_assele_tophad,
    e_assmuo_topele,
    e_assmuo_topmuo,
    e_assmuo_toptau,
    e_assmuo_tophad,
    e_asstau_topele,
    e_asstau_topmuo,
    e_asstau_toptau,
    e_asstau_tophad,
    e_asshad_topele,
    e_asshad_topmuo,
    e_asshad_toptau,
    e_asshad_tophad,
    e_notfound
  };

  GenParticle Top() const;
  GenParticle WTop() const;
  GenParticle bTop() const;
  GenParticle WTopDecay1() const;
  GenParticle WTopDecay2() const;
  GenParticle WAss() const;
  GenParticle WAssDecay1() const;
  GenParticle WAssDecay2() const;
  GenParticle LeptonicW() const;
  GenParticle SingleLepton() const;
  GenParticle SingleNeutrino() const;
  GenParticle LeptAss() const;
  GenParticle NuAss() const;
  GenParticle Initial1() const;
  GenParticle Initial2() const;
  GenParticle InitialBottomQuark() const;
  GenParticle gluonAss() const;
  GenParticle bAss() const;
  E_DecayChannel DecayChannel() const;

  bool IsTopHadronicDecay() const;
  bool IsTopLeptonicDecay(const bool consider_taus = true) const;
  bool IsTopToElectronDecay() const;
  bool IsTopToMuonDecay() const;
  bool IsTopToTauonDecay() const;

  bool IsAssHadronicDecay() const;
  bool IsAssLeptonicDecay(const bool consider_taus = true) const;
  bool IsAssToElectronDecay() const;
  bool IsAssToMuonDecay() const;
  bool IsAssToTauonDecay() const;

  bool IsSemiLeptonic(const bool consider_taus = true) const;

  bool IsGluonGluonProcess() const;
  bool IsBottomGluonProcess() const;
  bool HasAssociatedBottom() const;
  bool HasAssociatedGluon() const;

 private:

  GenParticle m_Top;
  GenParticle m_WTop;
  GenParticle m_bTop;
  GenParticle m_WTopDecay1;
  GenParticle m_WTopDecay2;
  GenParticle m_WAss;
  GenParticle m_WAssDecay1;
  GenParticle m_WAssDecay2;
  GenParticle m_initial1;
  GenParticle m_initial2;
  GenParticle m_gluonAss = GenParticle(); // make sure that this member variable is initialized
  GenParticle m_bAss = GenParticle(); // same
  bool m_has_bAss = false;

  E_DecayChannel m_type;
};


class SingleTopGen_tWchProducer: public uhh2::AnalysisModule {
public:
    explicit SingleTopGen_tWchProducer(uhh2::Context & ctx, const std::string & name = "singletopgen_twch", bool throw_on_failure = true);
    virtual bool process(uhh2::Event & event) override;

private:
    uhh2::Event::Handle<SingleTopGen_tWch> h_singletopgen_twch;
    bool throw_on_failure;
};

}}
