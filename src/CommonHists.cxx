#include "UHH2/common/include/Utils.h"

#include "UHH2/LegacyTopTagging/include/CommonHists.h"
#include "UHH2/LegacyTopTagging/include/Utils.h"

#include "TH1F.h"
#include "TH2F.h"

using namespace std;
using namespace uhh2;
using namespace ltt;


namespace uhh2 { namespace ltt {

CommonHists::CommonHists(Context & ctx, const string & dirname): Hists(ctx, dirname) {

  h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");

  hist_count = book<TH1F>("count", "Event count", 1, 0, 1);
  hist_weights_log10 = book<TH1F>("weights_log10", "log_{10}(event weight)", 1000, -4, 6);

  hist_ntrue = book<TH1F>("ntrue", "Number of true interactions", 101, -0.5, 100.5);
  hist_npv = book<TH1F>("npv", "Number of primary vertices", 101, -0.5, 100.5);

  hist_met_pt = book<TH1F>("met_pt", "#it{p}_{T}^{miss} [GeV]", 1000, 0, 1000);
  hist_met_phi = book<TH1F>("met_phi", "MET #phi [rad]", 1000, -M_PI, M_PI);

  hist_nmuons = book<TH1F>("nmuons", "Number of muons", 11, -0.5, 10.5);
  hist_nelectrons = book<TH1F>("nelectrons", "Number of electrons", 11, -0.5, 10.5);

  hist_primlep_pt = book<TH1F>("primlep_pt", "Lepton #it{p}_{T} [GeV]", 1000, 0, 1000);
  hist_ptw = book<TH1F>("ptw", "#it{p}_{T}^{W} [GeV]", 1000, 0, 1000);
  hist_primlep_eta = book<TH1F>("primlep_eta", "Lepton #eta", 1000, -5.0, 5.0);
  hist_primlep_phi = book<TH1F>("primlep_phi", "Lepton #phi [rad]", 1000, -M_PI, M_PI);
  hist_primlep_drjet = book<TH1F>("primlep_drjet", "#Delta#it{R}(lepton, closest AK4 jet)", 1000, 0, 5);
  hist_primlep_ptrel = book<TH1F>("primlep_ptrel", "#it{p}_{T}^{rel}(lepton, closest AK4 jet) [GeV]", 1000, 0, 500);
  hist_twodselection = book<TH2F>("twodselection", "#it{p}_{T}^{rel} [GeV] vs. #Delta#it{R}", 1000, 0, 5, 1000, 0, 500);

  hist_nak4 = book<TH1F>("nak4", "Number of AK4 jets", 11, -0.5, 10.5);
  hist_ak4jets_pt = book<TH1F>("ak4jets_pt", "AK4 jets: #it{p}_{T} [GeV]", 1000, 0, 1000);
  hist_ak4jets_eta = book<TH1F>("ak4jets_eta", "AK4 jets: #eta", 1000, -5.0, 5.0);
  hist_ak4jets_phi = book<TH1F>("ak4jets_phi", "AK4 jets: #phi [rad]", 1000, -M_PI, M_PI);
  hist_ak4jets_mass = book<TH1F>("ak4jets_mass", "AK4 jets: #it{m}_{jet} [GeV]", 1000, 0, 250);
  hist_ak4jets_drlepton = book<TH1F>("ak4jets_drlepton", "AK4 jets: #Delta#it{R}(AK4 jet, lepton)", 1000, 0, 5);
  hist_ak4jets_deepCSV = book<TH1F>("ak4jets_deepCSV", "AK4 jets: #it{O}_{DeepCSV}^{prob(b)+prob(bb)}", 1000, 0, 1);
  hist_ak4jets_deepJet = book<TH1F>("ak4jets_deepJet", "AK4 jets: #it{O}_{DeepJet}^{prob(b)+prob(bb)+prob(lepb)}", 1000, 0, 1);

  hist_ht_had = book<TH1F>("ht_had", "#it{H}_{T}^{had} [GeV]", 1000, 0, 2000);
  hist_ht_lep = book<TH1F>("ht_lep", "#it{H}_{T}^{lep} [GeV]", 1000, 0, 2000);
}


void CommonHists::fill(const Event & event) {

  const double w = event.weight;

  hist_count->Fill(0.5, w);
  if(w > 0) hist_weights_log10->Fill(log10(w), 1);
  else if(w < 0) hist_weights_log10->Fill(log10(-w), -1);

  if(!event.isRealData) hist_ntrue->Fill(event.genInfo->pileup_TrueNumInteractions(), w);
  hist_npv->Fill(event.pvs->size(), w);

  hist_met_pt->Fill(event.met->v4().Pt(), w);
  hist_met_phi->Fill(event.met->v4().Phi(), w);

  hist_nmuons->Fill(event.muons->size(), w);
  hist_nelectrons->Fill(event.electrons->size(), w);

  if(event.is_valid(h_primlep)) {
    const FlavorParticle & primlep = event.get(h_primlep);
    hist_primlep_pt->Fill(primlep.v4().Pt(), w);
    hist_ptw->Fill((event.met->v4() + primlep.v4()).Pt(), w);
    hist_primlep_eta->Fill(primlep.v4().Eta(), w);
    hist_primlep_phi->Fill(primlep.v4().Phi(), w);
    if(event.jets->size() > 0) {
      double drjet = deltaR(primlep.v4(), nextJet(primlep, *event.jets)->v4());
      double ptrel = pTrel(primlep, nextJet(primlep, *event.jets));
      hist_primlep_drjet->Fill(drjet, w);
      hist_primlep_ptrel->Fill(ptrel, w);
      hist_twodselection->Fill(drjet, ptrel, w);
    }
  }

  hist_nak4->Fill(event.jets->size(), w);
  for(const Jet & jet : *event.jets) {
    hist_ak4jets_pt->Fill(jet.v4().Pt(), w);
    hist_ak4jets_eta->Fill(jet.v4().Eta(), w);
    hist_ak4jets_phi->Fill(jet.v4().Phi(), w);
    hist_ak4jets_mass->Fill(jet.v4().M(), w);
    if(event.is_valid(h_primlep)) hist_ak4jets_drlepton->Fill(deltaR(jet.v4(), event.get(h_primlep).v4()), w);
    hist_ak4jets_deepCSV->Fill(jet.btag_DeepCSV(), w);
    hist_ak4jets_deepJet->Fill(jet.btag_DeepJet(), w);
  }

  double ht_had(0.);
  for(const Jet & jet : *event.jets) ht_had += jet.v4().Pt();
  hist_ht_had->Fill(ht_had, w);
  double ht_lep = event.met->v4().Pt();
  for(const Muon & muon : *event.muons) ht_lep += muon.v4().Pt();
  for(const Electron & electron : *event.electrons) ht_lep += electron.v4().Pt();
  hist_ht_lep->Fill(ht_lep, w);
}

}}
