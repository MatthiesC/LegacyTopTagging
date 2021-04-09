#include "UHH2/common/include/Utils.h"

#include "UHH2/LegacyTopTagging/include/HOTVRHists.h"
#include "UHH2/LegacyTopTagging/include/Utils.h"

#include "TH1F.h"
#include "TH2F.h"

using namespace std;
using namespace uhh2;
using namespace ltt;


namespace uhh2 { namespace ltt {

HOTVRHists::HOTVRHists(Context & ctx, const string & dirname, const string & coll_rec): Hists(ctx, dirname) {

  h_hotvrjets = ctx.get_handle<vector<TopJet>>(coll_rec.empty() ? "topjets" : coll_rec);

  hist_number = book<TH1F>("number", "Number of HOTVR jets", 11, -0.5, 10.5);

  hist_hotvrjets_pt = book<TH1F>("hotvrjets_pt", "HOTVR jets: #it{p}_{T} [GeV]", 1000, 0, 1000);
  hist_hotvrjets_ptlog = book<TH1F>("hotvrjets_ptlog", "HOTVR jets: log_{10}(#it{p}_{T} [GeV])", 1000, 0, 4);
  hist_hotvrjets_eta = book<TH1F>("hotvrjets_eta", "HOTVR jets: #eta", 1000, -5.0, 5.0);
  hist_hotvrjets_phi = book<TH1F>("hotvrjets_phi", "HOTVR jets: #phi [rad]", 1000, -M_PI, M_PI);
  hist_hotvrjets_mass = book<TH1F>("hotvrjets_mass", "HOTVR jets: #it{m}_{jet} [GeV]", 1000, 0, 500);
  hist_hotvrjets_mpair = book<TH1F>("hotvrjets_mpair", "HOTVR jets: min. #it{m}_{ij} [GeV] of leading three subjets", 1000, 0, 500);
  hist_hotvrjets_tau32 = book<TH1F>("hotvrjets_tau32", "HOTVR jets: #tau_{3}/#tau_{2}", 1000, 0, 1);
  hist_hotvrjets_fpt1 = book<TH1F>("hotvrjets_fpt1", "HOTVR jets: #it{p}_{T} fraction of leading subjet", 1000, 0, 1);

  hist_subjets_number = book<TH1F>("subjets_number", "HOTVR jets: #it{N} of subjets", 11, -0.5, 10.5);
  hist_subjets_pt = book<TH1F>("subjets_pt", "HOTVR subjets: #it{p}_{T} [GeV]", 1000, 0, 500);
  hist_subjets_eta = book<TH1F>("subjets_eta", "HOTVR subjets: #eta", 1000, -5.0, 5.0);
  hist_subjets_phi = book<TH1F>("subjets_phi", "HOTVR subjets: #phi [rad]", 1000, -M_PI, M_PI);
  hist_subjets_mass = book<TH1F>("subjets_mass", "HOTVR subjets: #it{m}_{jet} [GeV]", 1000, 0, 200);
  hist_subjets_fpt = book<TH1F>("subjets_fpt", "HOTVR subjets: #it{p}_{T} fraction", 1000, 0, 1);

  hist_hotvrjet1_pt = book<TH1F>("hotvrjet1_pt", "Leading HOTVR jet: #it{p}_{T} [GeV]", 1000, 0, 1000);
  hist_hotvrjet1_ptlog = book<TH1F>("hotvrjet1_ptlog", "Leading HOTVR jet: log_{10}(#it{p}_{T} [GeV])", 1000, 0, 4);
  hist_hotvrjet1_eta = book<TH1F>("hotvrjet1_eta", "Leading HOTVR jet: #eta", 1000, -5.0, 5.0);
  hist_hotvrjet1_phi = book<TH1F>("hotvrjet1_phi", "Leading HOTVR jet: #phi [rad]", 1000, -M_PI, M_PI);
  hist_hotvrjet1_mass = book<TH1F>("hotvrjet1_mass", "Leading HOTVR jet: #it{m}_{jet} [GeV]", 1000, 0, 500);
  hist_hotvrjet1_mpair = book<TH1F>("hotvrjet1_mpair", "Leading HOTVR jet: min. #it{m}_{ij} [GeV] of leading three subjets", 1000, 0, 500);
  hist_hotvrjet1_tau32 = book<TH1F>("hotvrjet1_tau32", "Leading HOTVR jet: #tau_{3}/#tau_{2}", 1000, 0, 1);
  hist_hotvrjet1_fpt1 = book<TH1F>("hotvrjet1_fpt1", "Leading HOTVR jet: #it{p}_{T} fraction of leading subjet", 1000, 0, 1);

  hist_hotvrjet2_pt = book<TH1F>("hotvrjet2_pt", "Subleading HOTVR jet: #it{p}_{T} [GeV]", 1000, 0, 1000);
  hist_hotvrjet2_ptlog = book<TH1F>("hotvrjet2_ptlog", "Subleading HOTVR jet: log_{10}(#it{p}_{T} [GeV])", 1000, 0, 4);
  hist_hotvrjet2_eta = book<TH1F>("hotvrjet2_eta", "Subleading HOTVR jet: #eta", 1000, -5.0, 5.0);
  hist_hotvrjet2_phi = book<TH1F>("hotvrjet2_phi", "Subleading HOTVR jet: #phi [rad]", 1000, -M_PI, M_PI);
  hist_hotvrjet2_mass = book<TH1F>("hotvrjet2_mass", "Subleading HOTVR jet: #it{m}_{jet} [GeV]", 1000, 0, 500);
  hist_hotvrjet2_mpair = book<TH1F>("hotvrjet2_mpair", "Subleading HOTVR jet: min. #it{m}_{ij} [GeV] of leading three subjets", 1000, 0, 500);
  hist_hotvrjet2_tau32 = book<TH1F>("hotvrjet2_tau32", "Subleading HOTVR jet: #tau_{3}/#tau_{2}", 1000, 0, 1);
  hist_hotvrjet2_fpt1 = book<TH1F>("hotvrjet2_fpt1", "Subleading HOTVR jet: #it{p}_{T} fraction of leading subjet", 1000, 0, 1);
}


void HOTVRHists::fill(const Event & event) {

  const double w = event.weight;

  vector<TopJet> hotvrjets = event.get(h_hotvrjets);
  sort_by_pt(hotvrjets);

  hist_number->Fill(hotvrjets.size(), w);

  for(const TopJet & hotvrjet : hotvrjets) {
    hist_hotvrjets_pt->Fill(hotvrjet.v4().Pt(), w);
    hist_hotvrjets_ptlog->Fill(log10(hotvrjet.v4().Pt()), w);
    hist_hotvrjets_eta->Fill(hotvrjet.v4().Eta(), w);
    hist_hotvrjets_phi->Fill(hotvrjet.v4().Phi(), w);
    hist_hotvrjets_mass->Fill(hotvrjet.v4().M(), w);
    hist_hotvrjets_mpair->Fill(HOTVR_mpair(hotvrjet, false), w);
    hist_hotvrjets_tau32->Fill(tau32groomed(hotvrjet), w);
    hist_hotvrjets_fpt1->Fill(HOTVR_fpt(hotvrjet), w);

    hist_subjets_number->Fill(hotvrjet.subjets().size(), w);
    for(const auto & subjet : hotvrjet.subjets()) {
      hist_subjets_pt->Fill(subjet.v4().Pt(), w);
      hist_subjets_eta->Fill(subjet.v4().Eta(), w);
      hist_subjets_phi->Fill(subjet.v4().Phi(), w);
      hist_subjets_mass->Fill(subjet.v4().M(), w);
      hist_subjets_fpt->Fill(subjet.v4().Pt() / hotvrjet.v4().Pt(), w);
    }
  }

  if(hotvrjets.size() >= 1) {
    hist_hotvrjet1_pt->Fill(hotvrjets.at(0).v4().Pt(), w);
    hist_hotvrjet1_ptlog->Fill(log10(hotvrjets.at(0).v4().Pt()), w);
    hist_hotvrjet1_eta->Fill(hotvrjets.at(0).v4().Eta(), w);
    hist_hotvrjet1_phi->Fill(hotvrjets.at(0).v4().Phi(), w);
    hist_hotvrjet1_mass->Fill(hotvrjets.at(0).v4().M(), w);
    hist_hotvrjet1_mpair->Fill(HOTVR_mpair(hotvrjets.at(0), false), w);
    hist_hotvrjet1_tau32->Fill(tau32groomed(hotvrjets.at(0)), w);
    hist_hotvrjet1_fpt1->Fill(HOTVR_fpt(hotvrjets.at(0)), w);
  }

  if(hotvrjets.size() >= 2) {
    hist_hotvrjet2_pt->Fill(hotvrjets.at(1).v4().Pt(), w);
    hist_hotvrjet2_ptlog->Fill(log10(hotvrjets.at(1).v4().Pt()), w);
    hist_hotvrjet2_eta->Fill(hotvrjets.at(1).v4().Eta(), w);
    hist_hotvrjet2_phi->Fill(hotvrjets.at(1).v4().Phi(), w);
    hist_hotvrjet2_mass->Fill(hotvrjets.at(1).v4().M(), w);
    hist_hotvrjet2_mpair->Fill(HOTVR_mpair(hotvrjets.at(1), false), w);
    hist_hotvrjet2_tau32->Fill(tau32groomed(hotvrjets.at(1)), w);
    hist_hotvrjet2_fpt1->Fill(HOTVR_fpt(hotvrjets.at(1)), w);
  }
}

}}
