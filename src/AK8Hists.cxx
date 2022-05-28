#include "UHH2/common/include/Utils.h"

#include "UHH2/LegacyTopTagging/include/AK8Hists.h"
#include "UHH2/LegacyTopTagging/include/Utils.h"

#include "TH1F.h"
#include "TH2F.h"

using namespace std;
using namespace uhh2;
using namespace ltt;


namespace uhh2 { namespace ltt {

AK8Hists::AK8Hists(Context & ctx, const string & dirname, const string & coll_rec, const string & coll_gen, const string & handle_name_tag, const bool doResponseHists_, const unsigned int default_nbins):
  Hists(ctx, dirname),
  doResponseHists(doResponseHists_)
{

  h_ak8jets = ctx.get_handle<vector<TopJet>>(coll_rec.empty() ? "topjets" : coll_rec);
  h_ak8genjets = ctx.get_handle<vector<GenTopJet>>(coll_gen.empty() ? "gentopjets" : coll_gen);
  h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
  // if(handle_name_tag != "") h_taggedjet = ctx.get_handle<TopJet>(handle_name_tag);
  h_taggedjet = ctx.get_handle<TopJet>(handle_name_tag);

  hist_number = book<TH1F>("number", "Number of AK8 jets", 11, -0.5, 10.5);

  hist_ak8jets_pt = book<TH1F>("ak8jets_pt", "AK8 jets: #it{p}_{T} [GeV]", default_nbins, 0, 1000);
  hist_ak8jets_drlepton = book<TH1F>("ak8jets_drlepton", "AK8 jets: #Delta#it{R}(AK8 jet, lepton)", default_nbins, 0, 5);
  hist_ak8jets_eta = book<TH1F>("ak8jets_eta", "AK8 jets: #eta", default_nbins, -5.0, 5.0);
  hist_ak8jets_phi = book<TH1F>("ak8jets_phi", "AK8 jets: #phi [rad]", default_nbins, -M_PI, M_PI);
  hist_ak8jets_mass = book<TH1F>("ak8jets_mass", "AK8 jets: #it{m}_{jet} [GeV]", default_nbins, 0, 500);
  hist_ak8jets_mSD = book<TH1F>("ak8jets_mSD", "AK8 jets: #it{m}_{SD} [GeV]", default_nbins, 0, 500);
  hist_ak8jets_tau32 = book<TH1F>("ak8jets_tau32", "AK8 jets: #tau_{3}/#tau_{2}", default_nbins, 0, 1);
  hist_ak8jets_tau21 = book<TH1F>("ak8jets_tau21", "AK8 jets: #tau_{2}/#tau_{1}", default_nbins, 0, 1);
  hist_ak8jets_maxDeepCSV = book<TH1F>("ak8jets_maxDeepCSV", "AK8 jets: max. #it{O}_{DeepCSV}^{prob(b)+prob(bb)} of subjets", default_nbins, 0, 1);
  hist_ak8jets_maxDeepJet = book<TH1F>("ak8jets_maxDeepJet", "AK8 jets: max. #it{O}_{DeepJet}^{prob(b)+prob(bb)+prob(lepb)} of subjets", default_nbins, 0, 1);
  hist_ak8jets_DeepAK8_TvsQCD = book<TH1F>("ak8jets_DeepAK8_TvsQCD", "AK8 jets: #it{O}_{DeepAK8}^{TvsQCD}", default_nbins, 0, 1);
  hist_ak8jets_DeepAK8_WvsQCD = book<TH1F>("ak8jets_DeepAK8_WvsQCD", "AK8 jets: #it{O}_{DeepAK8}^{WvsQCD}", default_nbins, 0, 1);
  hist_ak8jets_PNet_TvsQCD = book<TH1F>("ak8jets_PNet_TvsQCD", "AK8 jets: #it{O}_{ParticleNet}^{TvsQCD}", default_nbins, 0, 1);
  hist_ak8jets_PNet_WvsQCD = book<TH1F>("ak8jets_PNet_WvsQCD", "AK8 jets: #it{O}_{ParticleNet}^{WvsQCD}", default_nbins, 0, 1);

  hist_subjets_number = book<TH1F>("subjets_number", "AK8 jets: #it{N} of subjets", 11, -0.5, 10.5);
  hist_subjets_pt = book<TH1F>("subjets_pt", "AK8 subjets: #it{p}_{T} [GeV]", default_nbins, 0, 500);
  hist_subjets_eta = book<TH1F>("subjets_eta", "AK8 subjets: #eta", default_nbins, -5.0, 5.0);
  hist_subjets_phi = book<TH1F>("subjets_phi", "AK8 subjets: #phi [rad]", default_nbins, -M_PI, M_PI);
  hist_subjets_mass = book<TH1F>("subjets_mass", "AK8 subjets: #it{m}_{jet} [GeV]", default_nbins, 0, 250);
  hist_subjets_deepCSV = book<TH1F>("subjets_deepCSV", "AK8 subjets: #it{O}_{DeepCSV}^{prob(b)+prob(bb)}", default_nbins, 0, 1);
  hist_subjets_deepJet = book<TH1F>("subjets_deepJet", "AK8 subjets: #it{O}_{DeepJet}^{prob(b)+prob(bb)+prob(lepb)}", default_nbins, 0, 1);

  const string label_prefix = handle_name_tag != "dummy" ? "Tagged AK8 jet" : "Leading AK8 jet";
  hist_ak8jet1_pt = book<TH1F>("ak8jet1_pt", (label_prefix+": #it{p}_{T} [GeV]").c_str(), default_nbins, 0, 1000);
  hist_ak8jet1_drlepton = book<TH1F>("ak8jet1_drlepton", (label_prefix+": #Delta#it{R}(AK8 jet, lepton)").c_str(), default_nbins, 0, 5);
  hist_ak8jet1_eta = book<TH1F>("ak8jet1_eta", (label_prefix+": #eta").c_str(), default_nbins, -5.0, 5.0);
  hist_ak8jet1_phi = book<TH1F>("ak8jet1_phi", (label_prefix+": #phi [rad]").c_str(), default_nbins, -M_PI, M_PI);
  hist_ak8jet1_mass = book<TH1F>("ak8jet1_mass", (label_prefix+": #it{m}_{jet} [GeV]").c_str(), default_nbins, 0, 500);
  hist_ak8jet1_mSD = book<TH1F>("ak8jet1_mSD", (label_prefix+": #it{m}_{SD} [GeV]").c_str(), default_nbins, 0, 500);
  hist_ak8jet1_tau32 = book<TH1F>("ak8jet1_tau32", (label_prefix+": #tau_{3}/#tau_{2}").c_str(), default_nbins, 0, 1);
  hist_ak8jet1_tau21 = book<TH1F>("ak8jet1_tau21", (label_prefix+": #tau_{2}/#tau_{1}").c_str(), default_nbins, 0, 1);
  hist_ak8jet1_maxDeepCSV = book<TH1F>("ak8jet1_maxDeepCSV", (label_prefix+": max. #it{O}_{DeepCSV}^{prob(b)+prob(bb)} of subjets").c_str(), default_nbins, 0, 1);
  hist_ak8jet1_maxDeepJet = book<TH1F>("ak8jet1_maxDeepJet", (label_prefix+": max. #it{O}_{DeepJet}^{prob(b)+prob(bb)+prob(lepb)} of subjets").c_str(), default_nbins, 0, 1);
  hist_ak8jet1_DeepAK8_TvsQCD = book<TH1F>("ak8jet1_DeepAK8_TvsQCD", (label_prefix+": #it{O}_{DeepAK8}^{TvsQCD}").c_str(), default_nbins, 0, 1);
  hist_ak8jet1_DeepAK8_WvsQCD = book<TH1F>("ak8jet1_DeepAK8_WvsQCD", (label_prefix+": #it{O}_{DeepAK8}^{WvsQCD}").c_str(), default_nbins, 0, 1);
  hist_ak8jet1_PNet_TvsQCD = book<TH1F>("ak8jet1_PNet_TvsQCD", (label_prefix+": #it{O}_{ParticleNet}^{TvsQCD}").c_str(), default_nbins, 0, 1);
  hist_ak8jet1_PNet_WvsQCD = book<TH1F>("ak8jet1_PNet_WvsQCD", (label_prefix+": #it{O}_{ParticleNet}^{WvsQCD}").c_str(), default_nbins, 0, 1);

  // hist_ak8jet2_pt = book<TH1F>("ak8jet2_pt", "Subleading AK8 jet: #it{p}_{T} [GeV]", default_nbins, 0, 1000);
  // hist_ak8jet2_drlepton = book<TH1F>("ak8jet2_drlepton", "Subleading AK8 jet: #Delta#it{R}(AK8 jet, lepton)", default_nbins, 0, 5);
  // hist_ak8jet2_eta = book<TH1F>("ak8jet2_eta", "Subleading AK8 jet: #eta", default_nbins, -5.0, 5.0);
  // hist_ak8jet2_phi = book<TH1F>("ak8jet2_phi", "Subleading AK8 jet: #phi [rad]", default_nbins, -M_PI, M_PI);
  // hist_ak8jet2_mass = book<TH1F>("ak8jet2_mass", "Subleading AK8 jet: #it{m}_{jet} [GeV]", default_nbins, 0, 500);
  // hist_ak8jet2_mSD = book<TH1F>("ak8jet2_mSD", "Subleading AK8 jet: #it{m}_{SD} [GeV]", default_nbins, 0, 500);
  // hist_ak8jet2_tau32 = book<TH1F>("ak8jet2_tau32", "Subleading AK8 jet: #tau_{3}/#tau_{2}", default_nbins, 0, 1);
  // hist_ak8jet2_tau21 = book<TH1F>("ak8jet2_tau21", "Subleading AK8 jet: #tau_{2}/#tau_{1}", default_nbins, 0, 1);
  // hist_ak8jet2_maxDeepCSV = book<TH1F>("ak8jet2_maxDeepCSV", "Subleading AK8 jet: max. #it{O}_{DeepCSV}^{prob(b)+prob(bb)} of subjets", default_nbins, 0, 1);
  // hist_ak8jet2_maxDeepJet = book<TH1F>("ak8jet2_maxDeepJet", "Subleading AK8 jet: max. #it{O}_{DeepJet}^{prob(b)+prob(bb)+prob(lepb)} of subjets", default_nbins, 0, 1);
  // hist_ak8jet2_DeepAK8_TvsQCD = book<TH1F>("ak8jet2_DeepAK8_TvsQCD", "Subleading AK8 jet: #it{O}_{DeepAK8}^{TvsQCD}", default_nbins, 0, 1);
  // hist_ak8jet2_DeepAK8_WvsQCD = book<TH1F>("ak8jet2_DeepAK8_WvsQCD", "Subleading AK8 jet: #it{O}_{DeepAK8}^{WvsQCD}", default_nbins, 0, 1);
  // hist_ak8jet2_PNet_TvsQCD = book<TH1F>("ak8jet2_PNet_TvsQCD", "Subleading AK8 jet: #it{O}_{ParticleNet}^{TvsQCD}", default_nbins, 0, 1);
  // hist_ak8jet2_PNet_WvsQCD = book<TH1F>("ak8jet2_PNet_WvsQCD", "Subleading AK8 jet: #it{O}_{ParticleNet}^{WvsQCD}", default_nbins, 0, 1);

  if(doResponseHists) {
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
}


void AK8Hists::fill(const Event & event) {

  const double w = event.weight;
  FlavorParticle primlep = FlavorParticle();
  bool valid_primlep(false);
  if(event.is_valid(h_primlep)) {
    primlep = event.get(h_primlep);
    valid_primlep = true;
  }

  vector<TopJet> ak8jets = event.get(h_ak8jets);
  sort_by_pt(ak8jets);

  hist_number->Fill(ak8jets.size(), w);

  for(const TopJet & ak8jet : ak8jets) {
    hist_ak8jets_pt->Fill(ak8jet.v4().Pt(), w);
    if(valid_primlep) hist_ak8jets_drlepton->Fill(deltaR(ak8jet.v4(), primlep.v4()), w);
    hist_ak8jets_eta->Fill(ak8jet.v4().Eta(), w);
    hist_ak8jets_phi->Fill(ak8jet.v4().Phi(), w);
    hist_ak8jets_mass->Fill(ak8jet.v4().M(), w);
    hist_ak8jets_mSD->Fill(mSD(ak8jet), w);
    hist_ak8jets_tau32->Fill(tau32(ak8jet), w);
    hist_ak8jets_tau21->Fill(tau21(ak8jet), w);
    hist_ak8jets_maxDeepCSV->Fill(maxDeepCSVSubJetValue(ak8jet), w);
    hist_ak8jets_maxDeepJet->Fill(maxDeepJetSubJetValue(ak8jet), w);
    hist_ak8jets_DeepAK8_TvsQCD->Fill(ak8jet.btag_DeepBoosted_TvsQCD(), w);
    hist_ak8jets_DeepAK8_WvsQCD->Fill(ak8jet.btag_DeepBoosted_WvsQCD(), w);
    hist_ak8jets_PNet_TvsQCD->Fill(ak8jet.btag_ParticleNetDiscriminatorsJetTags_TvsQCD(), w);
    hist_ak8jets_PNet_WvsQCD->Fill(ak8jet.btag_ParticleNetDiscriminatorsJetTags_WvsQCD(), w);

    hist_subjets_number->Fill(ak8jet.subjets().size(), w);
    for(const auto & subjet : ak8jet.subjets()) {
      hist_subjets_pt->Fill(subjet.v4().Pt(), w);
      hist_subjets_eta->Fill(subjet.v4().Eta(), w);
      hist_subjets_phi->Fill(subjet.v4().Phi(), w);
      hist_subjets_mass->Fill(subjet.v4().M(), w);
      hist_subjets_deepCSV->Fill(subjet.btag_DeepCSV(), w);
      hist_subjets_deepJet->Fill(subjet.btag_DeepJet(), w);
    }
  }

  if(ak8jets.size() >= 1) {
    TopJet the_ak8_jet = ak8jets.at(0);
    if(event.is_valid(h_taggedjet)) the_ak8_jet = event.get(h_taggedjet);

    hist_ak8jet1_pt->Fill(the_ak8_jet.v4().Pt(), w);
    if(valid_primlep) hist_ak8jet1_drlepton->Fill(deltaR(the_ak8_jet.v4(), primlep.v4()), w);
    hist_ak8jet1_eta->Fill(the_ak8_jet.v4().Eta(), w);
    hist_ak8jet1_phi->Fill(the_ak8_jet.v4().Phi(), w);
    hist_ak8jet1_mass->Fill(the_ak8_jet.v4().M(), w);
    hist_ak8jet1_mSD->Fill(mSD(the_ak8_jet), w);
    hist_ak8jet1_tau32->Fill(tau32(the_ak8_jet), w);
    hist_ak8jet1_tau21->Fill(tau21(the_ak8_jet), w);
    hist_ak8jet1_maxDeepCSV->Fill(maxDeepCSVSubJetValue(the_ak8_jet), w);
    hist_ak8jet1_maxDeepJet->Fill(maxDeepJetSubJetValue(the_ak8_jet), w);
    hist_ak8jet1_DeepAK8_TvsQCD->Fill(the_ak8_jet.btag_DeepBoosted_TvsQCD(), w);
    hist_ak8jet1_DeepAK8_WvsQCD->Fill(the_ak8_jet.btag_DeepBoosted_WvsQCD(), w);
    hist_ak8jet1_PNet_TvsQCD->Fill(the_ak8_jet.btag_ParticleNetDiscriminatorsJetTags_TvsQCD(), w);
    hist_ak8jet1_PNet_WvsQCD->Fill(the_ak8_jet.btag_ParticleNetDiscriminatorsJetTags_WvsQCD(), w);
  }

  // if(ak8jets.size() >= 2) {
  //   hist_ak8jet2_pt->Fill(ak8jets.at(1).v4().Pt(), w);
  //   if(valid_primlep) hist_ak8jet2_drlepton->Fill(deltaR(ak8jets.at(1).v4(), primlep.v4()), w);
  //   hist_ak8jet2_eta->Fill(ak8jets.at(1).v4().Eta(), w);
  //   hist_ak8jet2_phi->Fill(ak8jets.at(1).v4().Phi(), w);
  //   hist_ak8jet2_mass->Fill(ak8jets.at(1).v4().M(), w);
  //   hist_ak8jet2_mSD->Fill(mSD(ak8jets.at(1)), w);
  //   hist_ak8jet2_tau32->Fill(tau32(ak8jets.at(1)), w);
  //   hist_ak8jet2_tau21->Fill(tau21(ak8jets.at(1)), w);
  //   hist_ak8jet2_maxDeepCSV->Fill(maxDeepCSVSubJetValue(ak8jets.at(1)), w);
  //   hist_ak8jet2_maxDeepJet->Fill(maxDeepJetSubJetValue(ak8jets.at(1)), w);
  //   hist_ak8jet2_DeepAK8_TvsQCD->Fill(ak8jets.at(1).btag_DeepBoosted_TvsQCD(), w);
  //   hist_ak8jet2_DeepAK8_WvsQCD->Fill(ak8jets.at(1).btag_DeepBoosted_WvsQCD(), w);
  //   hist_ak8jet2_PNet_TvsQCD->Fill(ak8jets.at(1).btag_ParticleNetDiscriminatorsJetTags_TvsQCD(), w);
  //   hist_ak8jet2_PNet_WvsQCD->Fill(ak8jets.at(1).btag_ParticleNetDiscriminatorsJetTags_WvsQCD(), w);
  // }

  if(doResponseHists) {
    if(!event.isRealData) {
      vector<GenTopJet> ak8genjets = event.get(h_ak8genjets);
      sort_by_pt(ak8genjets);

      for(const TopJet & recjet : ak8jets) {
        const GenTopJet *genjet = closestParticle(recjet, ak8genjets);
        if(genjet == nullptr) continue;

        const double pt_gen = genjet->v4().pt();
        const double pt_rec_corr = recjet.v4().pt();
        const double pt_rec_raw = recjet.v4().pt() * recjet.JEC_factor_raw();
        const double eta_rec = recjet.v4().eta();

        for(unsigned int i = 0; i <= fDRbins; i++) { // iterate over different deltaR values; if i = 0, then use 0.8
          const double dR = deltaR(genjet->v4(), recjet.v4());
          const double threshold = (i == 0) ? 0.8 : 0.1*i;
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
}

}}
