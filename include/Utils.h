#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h"


namespace uhh2 { namespace ltt {

template<typename T, typename U>
double deltaEta(const T & p1, const U & p2) {
    return fabs(p1.eta() - p2.eta());
}

double tau32(const TopJet & topjet);

double tau32groomed(const TopJet & topjet);

double mSD(const TopJet & topjet);

double maxDeepCSVSubJetValue(const TopJet & topjet);

double HOTVR_mpair(const TopJet & topjet);

double HOTVR_fpt(const TopJet & topjet, const unsigned int subjet_i = 0);

const TopJet * nextTopJet(const Particle & p, const std::vector<TopJet> & topjets);

class METSelection: public uhh2::Selection {
public:
  METSelection(const double _met_min = 0., const double _met_max = std::numeric_limits<double>::infinity());
  virtual bool passes(const uhh2::Event & event) override;
private:
  double met_min;
  double met_max;
};

class PTWSelection: public uhh2::Selection {
public:
  PTWSelection(uhh2::Context & ctx, const double _ptw_min = 0., const double _ptw_max = std::numeric_limits<double>::infinity());
  virtual bool passes(const uhh2::Event & event) override;
private:
  double ptw_min;
  double ptw_max;
  uhh2::Event::Handle<FlavorParticle> h_primlep;
};

class TwoDSelection: public uhh2::Selection {
public:
  TwoDSelection(uhh2::Context & ctx, const double _ptrel_min = 0., const double _dr_min = 0.);
  virtual bool passes(const uhh2::Event & event) override;
private:
  double ptrel_min;
  double dr_min;
  uhh2::Event::Handle<FlavorParticle> h_primlep;
};

}}
