#pragma once

#include "UHH2/core/include/AnalysisModule.h"
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

}}
