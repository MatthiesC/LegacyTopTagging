#include "UHH2/common/include/Utils.h"

#include "UHH2/LegacyTopTagging/include/HOTVRHists.h"
#include "UHH2/LegacyTopTagging/include/Utils.h"

#include "TH1F.h"
#include "TH2F.h"

using namespace std;
using namespace uhh2;
using namespace ltt;


namespace uhh2 { namespace ltt {

HOTVRHists::HOTVRHists(Context & ctx, const string & dirname, const string & coll_rec, const string & coll_gen, const string & handle_name_tag, const unsigned int default_nbins): Hists(ctx, dirname) {

  h_hotvrjets = ctx.get_handle<vector<TopJet>>(coll_rec.empty() ? "topjets" : coll_rec);
  h_hotvrgenjets = ctx.get_handle<vector<GenTopJet>>(coll_gen.empty() ? "gentopjets" : coll_gen);
  h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
  h_taggedjet = ctx.get_handle<TopJet>(handle_name_tag);

  hist_number = book<TH1F>("number", "Number of HOTVR jets", 11, -0.5, 10.5);

  hist_hotvrjets_pt = book<TH1F>("hotvrjets_pt", "HOTVR jets: #it{p}_{T} [GeV]", default_nbins, 0, 1000);
  hist_hotvrjets_drlepton = book<TH1F>("hotvrjets_drlepton", "HOTVR jets: #Delta#it{R}(HOTVR jet, lepton)", default_nbins, 0, 5);
  hist_hotvrjets_eta = book<TH1F>("hotvrjets_eta", "HOTVR jets: #eta", default_nbins, -5.0, 5.0);
  hist_hotvrjets_phi = book<TH1F>("hotvrjets_phi", "HOTVR jets: #phi [rad]", default_nbins, -M_PI, M_PI);
  hist_hotvrjets_mass = book<TH1F>("hotvrjets_mass", "HOTVR jets: #it{m}_{jet} [GeV]", default_nbins, 0, 500);
  hist_hotvrjets_mpair = book<TH1F>("hotvrjets_mpair", "HOTVR jets: min. #it{m}_{ij} [GeV] of leading three subjets", default_nbins, 0, 250);
  hist_hotvrjets_tau32 = book<TH1F>("hotvrjets_tau32", "HOTVR jets: #tau_{3}/#tau_{2}", default_nbins, 0, 1);
  hist_hotvrjets_fpt1 = book<TH1F>("hotvrjets_fpt1", "HOTVR jets: #it{p}_{T} fraction of leading subjet", default_nbins, 0, 1);

  hist_subjets_number = book<TH1F>("subjets_number", "HOTVR jets: #it{N} of subjets", 11, -0.5, 10.5);
  hist_subjets_pt = book<TH1F>("subjets_pt", "HOTVR subjets: #it{p}_{T} [GeV]", default_nbins, 0, 500);
  hist_subjets_eta = book<TH1F>("subjets_eta", "HOTVR subjets: #eta", default_nbins, -5.0, 5.0);
  hist_subjets_phi = book<TH1F>("subjets_phi", "HOTVR subjets: #phi [rad]", default_nbins, -M_PI, M_PI);
  hist_subjets_mass = book<TH1F>("subjets_mass", "HOTVR subjets: #it{m}_{jet} [GeV]", default_nbins, 0, 250);
  hist_subjets_fpt = book<TH1F>("subjets_fpt", "HOTVR subjets: #it{p}_{T} fraction", default_nbins, 0, 1);

  const string label_prefix = handle_name_tag != "" ? "Tagged HOTVR jet" : "Leading HOTVR jet";
  hist_hotvrjet1_pt = book<TH1F>("hotvrjet1_pt", (label_prefix+": #it{p}_{T} [GeV]").c_str(), default_nbins, 0, 1000);
  hist_hotvrjet1_drlepton = book<TH1F>("hotvrjet1_drlepton", (label_prefix+": #Delta#it{R}(HOTVR jet, lepton)").c_str(), default_nbins, 0, 5);
  hist_hotvrjet1_eta = book<TH1F>("hotvrjet1_eta", (label_prefix+": #eta").c_str(), default_nbins, -5.0, 5.0);
  hist_hotvrjet1_phi = book<TH1F>("hotvrjet1_phi", (label_prefix+": #phi [rad]").c_str(), default_nbins, -M_PI, M_PI);
  hist_hotvrjet1_mass = book<TH1F>("hotvrjet1_mass", (label_prefix+": #it{m}_{jet} [GeV]").c_str(), default_nbins, 0, 500);
  hist_hotvrjet1_mpair = book<TH1F>("hotvrjet1_mpair", (label_prefix+": min. #it{m}_{ij} [GeV] of leading three subjets").c_str(), default_nbins, 0, 250);
  hist_hotvrjet1_tau32 = book<TH1F>("hotvrjet1_tau32", (label_prefix+": #tau_{3}/#tau_{2}").c_str(), default_nbins, 0, 1);
  hist_hotvrjet1_fpt1 = book<TH1F>("hotvrjet1_fpt1", (label_prefix+": #it{p}_{T} fraction of leading subjet").c_str(), default_nbins, 0, 1);

  // hist_hotvrjet2_pt = book<TH1F>("hotvrjet2_pt", "Subleading HOTVR jet: #it{p}_{T} [GeV]", default_nbins, 0, 1000);
  // hist_hotvrjet2_drlepton = book<TH1F>("hotvrjet2_drlepton", "Subleading HOTVR jet: #Delta#it{R}(HOTVR jet, lepton)", default_nbins, 0, 5);
  // hist_hotvrjet2_eta = book<TH1F>("hotvrjet2_eta", "Subleading HOTVR jet: #eta", default_nbins, -5.0, 5.0);
  // hist_hotvrjet2_phi = book<TH1F>("hotvrjet2_phi", "Subleading HOTVR jet: #phi [rad]", default_nbins, -M_PI, M_PI);
  // hist_hotvrjet2_mass = book<TH1F>("hotvrjet2_mass", "Subleading HOTVR jet: #it{m}_{jet} [GeV]", default_nbins, 0, 500);
  // hist_hotvrjet2_mpair = book<TH1F>("hotvrjet2_mpair", "Subleading HOTVR jet: min. #it{m}_{ij} [GeV] of leading three subjets", default_nbins, 0, 250);
  // hist_hotvrjet2_tau32 = book<TH1F>("hotvrjet2_tau32", "Subleading HOTVR jet: #tau_{3}/#tau_{2}", default_nbins, 0, 1);
  // hist_hotvrjet2_fpt1 = book<TH1F>("hotvrjet2_fpt1", "Subleading HOTVR jet: #it{p}_{T} fraction of leading subjet", default_nbins, 0, 1);

  for(unsigned int i = 0; i <= fDRbins; i++) {
    const string dr_string = (string)"dr"+to_string(i/10)+"p"+to_string(i-(i/10)*10); // e.g. i=8 will be converted to "dr0p8" and i=12 will be converted to "dr1p2"
    hist_response_gen.push_back(book<TH1F>(((string)"response_gen_"+dr_string).c_str(), "#it{p}_{T}^{gen} [GeV]", 2000, 0, 2000));
    hist_response_corr.push_back(book<TH1F>(((string)"response_corr_"+dr_string).c_str(), "#it{p}_{T}^{gen} [GeV]", 2000, 0, 2000));
    hist_response_raw.push_back(book<TH1F>(((string)"response_raw_"+dr_string).c_str(), "#it{p}_{T}^{gen} [GeV]", 2000, 0, 2000));
    hist_response_eta2p5_gen.push_back(book<TH1F>(((string)"response_eta2p5_gen_"+dr_string).c_str(), "#it{p}_{T}^{gen} [GeV]", 2000, 0, 2000));
    hist_response_eta2p5_corr.push_back(book<TH1F>(((string)"response_eta2p5_corr_"+dr_string).c_str(), "#it{p}_{T}^{gen} [GeV]", 2000, 0, 2000));
    hist_response_eta2p5_raw.push_back(book<TH1F>(((string)"response_eta2p5_raw_"+dr_string).c_str(), "#it{p}_{T}^{gen} [GeV]", 2000, 0, 2000));
  }
}


void HOTVRHists::fill(const Event & event) {

  const double w = event.weight;
  FlavorParticle primlep = FlavorParticle();
  bool valid_primlep(false);
  if(event.is_valid(h_primlep)) {
    primlep = event.get(h_primlep);
    valid_primlep = true;
  }

  vector<TopJet> hotvrjets = event.get(h_hotvrjets);
  sort_by_pt(hotvrjets);

  hist_number->Fill(hotvrjets.size(), w);

  for(const TopJet & hotvrjet : hotvrjets) {
    hist_hotvrjets_pt->Fill(hotvrjet.v4().Pt(), w);
    if(valid_primlep) hist_hotvrjets_drlepton->Fill(deltaR(hotvrjet.v4(), primlep.v4()), w);
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
    TopJet the_hotvr_jet = hotvrjets.at(0);
    if(event.is_valid(h_taggedjet)) the_hotvr_jet = event.get(h_taggedjet);

    hist_hotvrjet1_pt->Fill(the_hotvr_jet.v4().Pt(), w);
    if(valid_primlep) hist_hotvrjet1_drlepton->Fill(deltaR(the_hotvr_jet.v4(), primlep.v4()), w);
    hist_hotvrjet1_eta->Fill(the_hotvr_jet.v4().Eta(), w);
    hist_hotvrjet1_phi->Fill(the_hotvr_jet.v4().Phi(), w);
    hist_hotvrjet1_mass->Fill(the_hotvr_jet.v4().M(), w);
    hist_hotvrjet1_mpair->Fill(HOTVR_mpair(the_hotvr_jet, false), w);
    hist_hotvrjet1_tau32->Fill(tau32groomed(the_hotvr_jet), w);
    hist_hotvrjet1_fpt1->Fill(HOTVR_fpt(the_hotvr_jet), w);
  }

  // if(hotvrjets.size() >= 2) {
  //   hist_hotvrjet2_pt->Fill(hotvrjets.at(1).v4().Pt(), w);
  //   if(valid_primlep) hist_hotvrjet2_drlepton->Fill(deltaR(hotvrjets.at(1).v4(), primlep.v4()), w);
  //   hist_hotvrjet2_eta->Fill(hotvrjets.at(1).v4().Eta(), w);
  //   hist_hotvrjet2_phi->Fill(hotvrjets.at(1).v4().Phi(), w);
  //   hist_hotvrjet2_mass->Fill(hotvrjets.at(1).v4().M(), w);
  //   hist_hotvrjet2_mpair->Fill(HOTVR_mpair(hotvrjets.at(1), false), w);
  //   hist_hotvrjet2_tau32->Fill(tau32groomed(hotvrjets.at(1)), w);
  //   hist_hotvrjet2_fpt1->Fill(HOTVR_fpt(hotvrjets.at(1)), w);
  // }

  if(!event.isRealData) {
    vector<GenTopJet> hotvrgenjets = event.get(h_hotvrgenjets);
    sort_by_pt(hotvrgenjets);

    for(const TopJet & recjet : hotvrjets) {
      const GenTopJet *genjet = closestParticle(recjet, hotvrgenjets);
      if(genjet == nullptr) continue;

      const double pt_gen = genjet->v4().pt();
      const double pt_rec_corr = recjet.v4().pt();
      const double pt_rec_raw = recjet.v4().pt() * recjet.JEC_factor_raw();
      const double eta_rec = recjet.v4().eta();

      for(unsigned int i = 0; i <= fDRbins; i++) { // iterate over different deltaR values; if i = 0, then use HOTVR_Reff
        const double dR = deltaR(genjet->v4(), recjet.v4());
        const double threshold = (i == 0) ? HOTVR_Reff(recjet) : 0.1*i;
        if(dR > threshold) continue;
        hist_response_gen.at(i)->Fill(pt_gen, w);
        hist_response_corr.at(i)->Fill(pt_gen, w * pt_rec_corr / pt_gen);
        hist_response_raw.at(i)->Fill(pt_gen, w * pt_rec_raw / pt_gen);
        if(eta_rec > 2.5) continue;
        hist_response_eta2p5_gen.at(i)->Fill(pt_gen, w);
        hist_response_eta2p5_corr.at(i)->Fill(pt_gen, w * pt_rec_corr / pt_gen);
        hist_response_eta2p5_raw.at(i)->Fill(pt_gen, w * pt_rec_raw / pt_gen);
      }
    }
  }
}

}}
