#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"


namespace uhh2 { namespace ltt {

class HOTVRHists: public uhh2::Hists {
public:
  HOTVRHists(uhh2::Context & ctx, const std::string & dirname, const std::string & coll_rec = "", const std::string & coll_gen = "", const std::string & handle_name_tag = "HOTVRHists_dummy_handle", const unsigned int default_nbins = 100);

  virtual void fill(const uhh2::Event & event) override;

protected:
  TH1F *hist_number;

  TH1F *hist_hotvrjets_pt;
  TH1F *hist_hotvrjets_drlepton;
  TH1F *hist_hotvrjets_eta;
  TH1F *hist_hotvrjets_phi;
  TH1F *hist_hotvrjets_mass;
  TH1F *hist_hotvrjets_mpair;
  TH1F *hist_hotvrjets_tau32;
  TH1F *hist_hotvrjets_fpt1;

  TH1F *hist_subjets_number;
  TH1F *hist_subjets_pt;
  TH1F *hist_subjets_eta;
  TH1F *hist_subjets_phi;
  TH1F *hist_subjets_mass;
  TH1F *hist_subjets_fpt;

  TH1F *hist_hotvrjet1_pt;
  TH1F *hist_hotvrjet1_drlepton;
  TH1F *hist_hotvrjet1_eta;
  TH1F *hist_hotvrjet1_phi;
  TH1F *hist_hotvrjet1_mass;
  TH1F *hist_hotvrjet1_mpair;
  TH1F *hist_hotvrjet1_tau32;
  TH1F *hist_hotvrjet1_fpt1;

  // TH1F *hist_hotvrjet2_pt;
  // TH1F *hist_hotvrjet2_drlepton;
  // TH1F *hist_hotvrjet2_eta;
  // TH1F *hist_hotvrjet2_phi;
  // TH1F *hist_hotvrjet2_mass;
  // TH1F *hist_hotvrjet2_mpair;
  // TH1F *hist_hotvrjet2_tau32;
  // TH1F *hist_hotvrjet2_fpt1;

  const unsigned int fDRbins = 15;

  std::vector<TH1F*> hist_response_gen;
  std::vector<TH1F*> hist_response_corr;
  std::vector<TH1F*> hist_response_raw;
  std::vector<TH1F*> hist_response_eta2p5_gen;
  std::vector<TH1F*> hist_response_eta2p5_corr;
  std::vector<TH1F*> hist_response_eta2p5_raw;

private:
  uhh2::Event::Handle<std::vector<TopJet>> h_hotvrjets;
  uhh2::Event::Handle<std::vector<GenTopJet>> h_hotvrgenjets;
  uhh2::Event::Handle<FlavorParticle> h_primlep;
  uhh2::Event::Handle<TopJet> h_taggedjet;
};

}}
