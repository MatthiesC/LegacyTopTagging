#include "UHH2/common/include/Utils.h"

#include "UHH2/LegacyTopTagging/include/AK4Hists.h"
#include "UHH2/LegacyTopTagging/include/Utils.h"
#include "UHH2/LegacyTopTagging/include/Constants.h"

#include "TH1F.h"
#include "TH2F.h"

using namespace std;
using namespace uhh2;
using namespace ltt;


namespace uhh2 { namespace ltt {

AK4Hists::AK4Hists(Context & ctx, const string & dirname, const unsigned int default_nbins):
  Hists(ctx, dirname),
  fHandle_PUPPIjets(ctx.get_handle<vector<Jet>>("jets")),
  fHandle_CHSjets(ctx.get_handle<vector<Jet>>(kCollectionName_AK4CHS)),
  fHandle_pairedPUPPIjets(ctx.get_handle<vector<Jet>>(kHandleName_pairedPUPPIjets)),
  fHandle_forwardPUPPIjets(ctx.get_handle<vector<Jet>>(kHandleName_forwardPUPPIjets)),
  fHandle_uncleanedPUPPIjets(ctx.get_handle<vector<Jet>>(kHandleName_uncleanedPUPPIjets)),
  fHandle_bJets(ctx.get_handle<vector<Jet>>(kHandleName_bJets)),
  fHandle_bJets_loose(ctx.get_handle<vector<Jet>>(kHandleName_bJets_loose)),
  fHandle_bJets_medium(ctx.get_handle<vector<Jet>>(kHandleName_bJets_medium)),
  fHandle_bJets_tight(ctx.get_handle<vector<Jet>>(kHandleName_bJets_tight)),
  fHandle_PrimaryLepton(ctx.get_handle<FlavorParticle>(kHandleName_PrimaryLepton)),
  fHandle_weight_btagdisc_central(ctx.get_handle<float>("weight_btagdisc__central")), // same handle name as defined in MCWeight.cxx --> MCBTagDiscriminantReweighting
  fHandle_weight_btag_njet_sf(ctx.get_handle<float>(kHandleName_weight_btag_njet_sf))
{
  hist_number_puppijets = book<TH1F>("number_puppijets", "Number of AK4 PUPPI jets", 11, -0.5, 10.5);
  hist_number_puppijets_uncleaned = book<TH1F>("number_puppijets_uncleaned", "Number of AK4 PUPPI jets before CHS-matching", 11, -0.5, 10.5);
  hist_number_puppijets_central = book<TH1F>("number_puppijets_central", "Number of central AK4 PUPPI jets", 11, -0.5, 10.5);
  hist_number_puppijets_central_wo_btag_sf = book<TH1F>("number_puppijets_central_wo_btag_sf", "Number of central AK4 PUPPI jets (w/o b-tag SF)", 11, -0.5, 10.5);
  hist_number_puppijets_central_wo_njet_sf = book<TH1F>("number_puppijets_central_wo_njet_sf", "Number of central AK4 PUPPI jets (w/o njet SF)", 11, -0.5, 10.5);
  hist_number_puppijets_central_wo_btag_sf_wo_njet_sf = book<TH1F>("number_puppijets_central_wo_btag_sf_wo_njet_sf", "Number of central AK4 PUPPI jets (w/o b-tag SF, w/o njet SF)", 11, -0.5, 10.5);
  hist_number_puppijets_forward = book<TH1F>("number_puppijets_forward", "Number of forward AK4 PUPPI jets", 11, -0.5, 10.5);

  hist_puppichs_dr = book<TH1F>("puppichs_dr", "#Delta#it{R}(PUPPI jet, matched CHS jet)", default_nbins, 0, 0.5);
  hist_puppichs_ptresponse = book<TH1F>("puppichs_ptresponse", "#it{p}_{T}^{PUPPI} / #it{p}_{T}^{CHS}", default_nbins, 0.5, 1.5);

  hist_number_bjets = book<TH1F>("number_bjets", "Number of b jets", 11, -0.5, 10.5);
  hist_number_bjets_loose = book<TH1F>("number_bjets_loose", "Number of b jets (loose)", 11, -0.5, 10.5);
  hist_number_bjets_medium = book<TH1F>("number_bjets_medium", "Number of b jets (medium)", 11, -0.5, 10.5);
  hist_number_bjets_tight = book<TH1F>("number_bjets_tight", "Number of b jets (tight)", 11, -0.5, 10.5);

  hist_number_bjets_loose_medium = book<TH2F>("number_bjets_loose_medium", "Number of b jets (x = loose, y = medium)", 11, -0.5, 10.5, 11, -0.5, 10.5);
  hist_number_bjets_loose_tight = book<TH2F>("number_bjets_loose_tight", "Number of b jets (x = loose, y = tight)", 11, -0.5, 10.5, 11, -0.5, 10.5);
  hist_number_bjets_medium_tight = book<TH2F>("number_bjets_medium_tight", "Number of b jets (x = medium, y = tight)", 11, -0.5, 10.5, 11, -0.5, 10.5);

  hist_paired_puppijets_pt = book<TH1F>("paired_puppijets_pt", "Central AK4 jets: #it{p}_{T} [GeV]", default_nbins, 0, 1000);
  hist_paired_puppijets_drlepton = book<TH1F>("paired_puppijets_drlepton", "Central AK4 jets: #Delta#it{R}(AK4 jet, lepton)", default_nbins, 0, 5);
  hist_paired_puppijets_eta = book<TH1F>("paired_puppijets_eta", "Central AK4 jets: #eta", default_nbins, -5.0, 5.0);
  hist_paired_puppijets_phi = book<TH1F>("paired_puppijets_phi", "Central AK4 jets: #phi [rad]", default_nbins, -M_PI, M_PI);
  hist_paired_puppijets_mass = book<TH1F>("paired_puppijets_mass", "Central AK4 jets: #it{m}_{jet} [GeV]", default_nbins, 0, 200);
  hist_paired_puppijets_deepjet = book<TH1F>("paired_puppijets_deepjet", "Central AK4 jets: #it{O}(DeepJet)", default_nbins, 0, 1);


  hist_puppijets_pt = book<TH1F>("puppijets_pt", "AK4 jets: #it{p}_{T} [GeV]", default_nbins, 0, 1000);
  hist_forward_puppijets_pt = book<TH1F>("forward_puppijets_pt", "Forward AK4 jets: #it{p}_{T} [GeV]", default_nbins, 0, 1000);

  hist_puppijets_drlepton = book<TH1F>("puppijets_drlepton", "AK4 jets: #Delta#it{R}(AK4 jet, lepton)", default_nbins, 0, 5);
  hist_forward_puppijets_drlepton = book<TH1F>("forward_puppijets_drlepton", "Forward AK4 jets: #Delta#it{R}(AK4 jet, lepton)", default_nbins, 0, 5);

  hist_puppijets_eta = book<TH1F>("puppijets_eta", "AK4 jets: #eta", default_nbins, -5.0, 5.0);
  hist_forward_puppijets_eta = book<TH1F>("forward_puppijets_eta", "Forward AK4 jets: #eta", default_nbins, -5.0, 5.0);

  hist_puppijets_phi = book<TH1F>("puppijets_phi", "AK4 jets: #phi [rad]", default_nbins, -M_PI, M_PI);
  hist_forward_puppijets_phi = book<TH1F>("forward_puppijets_phi", "Forward AK4 jets: #phi [rad]", default_nbins, -M_PI, M_PI);

  hist_puppijets_mass = book<TH1F>("puppijets_mass", "AK4 jets: #it{m}_{jet} [GeV]", default_nbins, 0, 200);
  hist_forward_puppijets_mass = book<TH1F>("forward_puppijets_mass", "Forward AK4 jets: #it{m}_{jet} [GeV]", default_nbins, 0, 200);


  hist_paired_puppijet1pt_pt = book<TH1F>("paired_puppijet1pt_pt", "#it{p}_{T}-leading central AK4 jet: #it{p}_{T} [GeV]", default_nbins, 0, 1000);
  hist_paired_puppijet1dj_pt = book<TH1F>("paired_puppijet1dj_pt", "DeepJet-leading central AK4 jet: #it{p}_{T} [GeV]", default_nbins, 0, 1000);

  hist_paired_puppijet1pt_drlepton = book<TH1F>("paired_puppijet1pt_drlepton", "#it{p}_{T}-leading central AK4 jet: #Delta#it{R}(AK4 jet, lepton)", default_nbins, 0, 5);
  hist_paired_puppijet1dj_drlepton = book<TH1F>("paired_puppijet1dj_drlepton", "DeepJet-leading central AK4 jet: #Delta#it{R}(AK4 jet, lepton)", default_nbins, 0, 5);

  hist_paired_puppijet1pt_eta = book<TH1F>("paired_puppijet1pt_eta", "#it{p}_{T}-leading central AK4 jet: #eta", default_nbins, -5.0, 5.0);
  hist_paired_puppijet1dj_eta = book<TH1F>("paired_puppijet1dj_eta", "DeepJet-leading central AK4 jet: #eta", default_nbins, -5.0, 5.0);

  hist_paired_puppijet1pt_phi = book<TH1F>("paired_puppijet1pt_phi", "#it{p}_{T}-leading central AK4 jet: #phi [rad]", default_nbins, -M_PI, M_PI);
  hist_paired_puppijet1dj_phi = book<TH1F>("paired_puppijet1dj_phi", "DeepJet-leading central AK4 jet: #phi [rad]", default_nbins, -M_PI, M_PI);

  hist_paired_puppijet1pt_mass = book<TH1F>("paired_puppijet1pt_mass", "#it{p}_{T}-leading central AK4 jet: #it{m}_{jet} [GeV]", default_nbins, 0, 200);
  hist_paired_puppijet1dj_mass = book<TH1F>("paired_puppijet1dj_mass", "DeepJet-leading central AK4 jet: #it{m}_{jet} [GeV]", default_nbins, 0, 200);

  hist_paired_puppijet1pt_deepjet = book<TH1F>("paired_puppijet1pt_deepjet", "#it{p}_{T}-leading central AK4 jet: #it{O}(DeepJet)", default_nbins, 0, 1);
  hist_paired_puppijet1dj_deepjet = book<TH1F>("paired_puppijet1dj_deepjet", "DeepJet-leading central AK4 jet: #it{O}(DeepJet)", default_nbins, 0, 1);
}


void AK4Hists::fill(const Event & event) {

  const double w = event.weight;
  FlavorParticle primlep = FlavorParticle();
  bool valid_primlep(false);
  if(event.is_valid(fHandle_PrimaryLepton)) {
    primlep = event.get(fHandle_PrimaryLepton);
    valid_primlep = true;
  }
  bool matching_done(false);
  vector<Jet> puppijets = event.get(fHandle_PUPPIjets);
  sort_by_pt<Jet>(puppijets);
  vector<Jet> paired_puppijets;
  vector<Jet> paired_puppijets_dj;
  vector<Jet> forward_puppijets;
  if(event.is_valid(fHandle_pairedPUPPIjets)) {
    paired_puppijets = event.get(fHandle_pairedPUPPIjets);
    sort_by_pt<Jet>(paired_puppijets);
    forward_puppijets = event.get(fHandle_forwardPUPPIjets);
    sort_by_pt<Jet>(forward_puppijets);
    paired_puppijets_dj = paired_puppijets;
    sort_by_deepjet_from_matches(paired_puppijets_dj, event, fHandle_CHSjets);
    matching_done = true;
  }

  hist_number_puppijets->Fill(puppijets.size(), w);
  if(matching_done) {
    hist_number_puppijets_uncleaned->Fill(event.get(fHandle_uncleanedPUPPIjets).size(), w);
    hist_number_puppijets_central->Fill(paired_puppijets.size(), w);

    const float divisor1 = event.is_valid(fHandle_weight_btagdisc_central) && event.get(fHandle_weight_btagdisc_central) != 0 ? event.get(fHandle_weight_btagdisc_central) : 1.f;
    hist_number_puppijets_central_wo_btag_sf->Fill(paired_puppijets.size(), w / divisor1);
    const float divisor2 = event.is_valid(fHandle_weight_btag_njet_sf) && event.get(fHandle_weight_btag_njet_sf) != 0 ? event.get(fHandle_weight_btag_njet_sf) : 1.f;
    hist_number_puppijets_central_wo_njet_sf->Fill(paired_puppijets.size(), w / divisor2);
    hist_number_puppijets_central_wo_btag_sf_wo_njet_sf->Fill(paired_puppijets.size(), w / (divisor1 * divisor2));

    hist_number_puppijets_forward->Fill(forward_puppijets.size(), w);
    for(const Jet & puppijet : paired_puppijets) {
      const Jet *chsjet = getCHSmatch(puppijet, event, fHandle_CHSjets);
      hist_puppichs_dr->Fill(deltaR(puppijet.v4(), chsjet->v4()), w);
      hist_puppichs_ptresponse->Fill(puppijet.v4().pt() / chsjet->v4().pt(), w);
    }

    hist_number_bjets->Fill(event.get(fHandle_bJets).size(), w);
    hist_number_bjets_loose->Fill(event.get(fHandle_bJets_loose).size(), w);
    hist_number_bjets_medium->Fill(event.get(fHandle_bJets_medium).size(), w);
    hist_number_bjets_tight->Fill(event.get(fHandle_bJets_tight).size(), w);

    hist_number_bjets_loose_medium->Fill(event.get(fHandle_bJets_loose).size(), event.get(fHandle_bJets_medium).size(), w);
    hist_number_bjets_loose_tight->Fill(event.get(fHandle_bJets_loose).size(), event.get(fHandle_bJets_tight).size(), w);
    hist_number_bjets_medium_tight->Fill(event.get(fHandle_bJets_medium).size(), event.get(fHandle_bJets_tight).size(), w);
  }

  if(matching_done) {
    for(const Jet & puppijet : paired_puppijets) {
      hist_paired_puppijets_pt->Fill(puppijet.v4().Pt(), w);
      if(valid_primlep) hist_paired_puppijets_drlepton->Fill(deltaR(puppijet.v4(), primlep.v4()), w);
      hist_paired_puppijets_eta->Fill(puppijet.v4().Eta(), w);
      hist_paired_puppijets_phi->Fill(puppijet.v4().Phi(), w);
      hist_paired_puppijets_mass->Fill(puppijet.v4().M(), w);
      hist_paired_puppijets_deepjet->Fill(getCHSmatch(puppijet, event, fHandle_CHSjets)->btag_DeepJet(), w);
    }

    for(const Jet & puppijet : forward_puppijets) {
      hist_forward_puppijets_pt->Fill(puppijet.v4().Pt(), w);
      if(valid_primlep) hist_forward_puppijets_drlepton->Fill(deltaR(puppijet.v4(), primlep.v4()), w);
      hist_forward_puppijets_eta->Fill(puppijet.v4().Eta(), w);
      hist_forward_puppijets_phi->Fill(puppijet.v4().Phi(), w);
      hist_forward_puppijets_mass->Fill(puppijet.v4().M(), w);
    }
  }

  for(const Jet & puppijet : puppijets) {
    hist_puppijets_pt->Fill(puppijet.v4().Pt(), w);
    if(valid_primlep) hist_puppijets_drlepton->Fill(deltaR(puppijet.v4(), primlep.v4()), w);
    hist_puppijets_eta->Fill(puppijet.v4().Eta(), w);
    hist_puppijets_phi->Fill(puppijet.v4().Phi(), w);
    hist_puppijets_mass->Fill(puppijet.v4().M(), w);
  }

  if(matching_done) {
    if(paired_puppijets.size() >= 1) {
      hist_paired_puppijet1pt_pt->Fill(paired_puppijets.at(0).v4().Pt(), w);
      if(valid_primlep) hist_paired_puppijet1pt_drlepton->Fill(deltaR(paired_puppijets.at(0).v4(), primlep.v4()), w);
      hist_paired_puppijet1pt_eta->Fill(paired_puppijets.at(0).v4().Eta(), w);
      hist_paired_puppijet1pt_phi->Fill(paired_puppijets.at(0).v4().Phi(), w);
      hist_paired_puppijet1pt_mass->Fill(paired_puppijets.at(0).v4().M(), w);
      hist_paired_puppijet1pt_deepjet->Fill(getCHSmatch(paired_puppijets.at(0), event, fHandle_CHSjets)->btag_DeepJet(), w);
    }

    if(paired_puppijets_dj.size() >= 1) {
      hist_paired_puppijet1dj_pt->Fill(paired_puppijets_dj.at(0).v4().Pt(), w);
      if(valid_primlep) hist_paired_puppijet1dj_drlepton->Fill(deltaR(paired_puppijets_dj.at(0).v4(), primlep.v4()), w);
      hist_paired_puppijet1dj_eta->Fill(paired_puppijets_dj.at(0).v4().Eta(), w);
      hist_paired_puppijet1dj_phi->Fill(paired_puppijets_dj.at(0).v4().Phi(), w);
      hist_paired_puppijet1dj_mass->Fill(paired_puppijets_dj.at(0).v4().M(), w);
      hist_paired_puppijet1dj_deepjet->Fill(getCHSmatch(paired_puppijets_dj.at(0), event, fHandle_CHSjets)->btag_DeepJet(), w);
    }
  }
}

}}
