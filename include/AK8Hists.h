#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"


namespace uhh2 { namespace ltt {

class AK8Hists: public uhh2::Hists {
public:
  AK8Hists(uhh2::Context & ctx, const std::string & dirname, const std::string & coll_rec = "", const std::string & coll_gen = "", const unsigned int default_nbins = 1000);

  virtual void fill(const uhh2::Event & event) override;

protected:
  TH1F *hist_number;

  TH1F *hist_ak8jets_pt;
  TH1F *hist_ak8jets_drlepton;
  TH1F *hist_ak8jets_eta;
  TH1F *hist_ak8jets_phi;
  TH1F *hist_ak8jets_mass;
  TH1F *hist_ak8jets_mSD;
  TH1F *hist_ak8jets_tau32;
  TH1F *hist_ak8jets_tau21;
  TH1F *hist_ak8jets_maxDeepCSV;
  TH1F *hist_ak8jets_maxDeepJet;
  TH1F *hist_ak8jets_DeepAK8_TvsQCD;
  TH1F *hist_ak8jets_DeepAK8_WvsQCD;
  TH1F *hist_ak8jets_PNet_TvsQCD;
  TH1F *hist_ak8jets_PNet_WvsQCD;

  TH1F *hist_subjets_number;
  TH1F *hist_subjets_pt;
  TH1F *hist_subjets_eta;
  TH1F *hist_subjets_phi;
  TH1F *hist_subjets_mass;
  TH1F *hist_subjets_deepCSV;
  TH1F *hist_subjets_deepJet;

  TH1F *hist_ak8jet1_pt;
  TH1F *hist_ak8jet1_drlepton;
  TH1F *hist_ak8jet1_eta;
  TH1F *hist_ak8jet1_phi;
  TH1F *hist_ak8jet1_mass;
  TH1F *hist_ak8jet1_mSD;
  TH1F *hist_ak8jet1_tau32;
  TH1F *hist_ak8jet1_tau21;
  TH1F *hist_ak8jet1_maxDeepCSV;
  TH1F *hist_ak8jet1_maxDeepJet;
  TH1F *hist_ak8jet1_DeepAK8_TvsQCD;
  TH1F *hist_ak8jet1_DeepAK8_WvsQCD;
  TH1F *hist_ak8jet1_PNet_TvsQCD;
  TH1F *hist_ak8jet1_PNet_WvsQCD;

  TH1F *hist_ak8jet2_pt;
  TH1F *hist_ak8jet2_drlepton;
  TH1F *hist_ak8jet2_eta;
  TH1F *hist_ak8jet2_phi;
  TH1F *hist_ak8jet2_mass;
  TH1F *hist_ak8jet2_mSD;
  TH1F *hist_ak8jet2_tau32;
  TH1F *hist_ak8jet2_tau21;
  TH1F *hist_ak8jet2_maxDeepCSV;
  TH1F *hist_ak8jet2_maxDeepJet;
  TH1F *hist_ak8jet2_DeepAK8_TvsQCD;
  TH1F *hist_ak8jet2_DeepAK8_WvsQCD;
  TH1F *hist_ak8jet2_PNet_TvsQCD;
  TH1F *hist_ak8jet2_PNet_WvsQCD;

  const unsigned int fDRbins = 15;

  std::vector<TH1F*> hist_response_gen;
  std::vector<TH1F*> hist_response_corr;
  std::vector<TH1F*> hist_response_raw;
  std::vector<TH1F*> hist_response_eta2p5_gen;
  std::vector<TH1F*> hist_response_eta2p5_corr;
  std::vector<TH1F*> hist_response_eta2p5_raw;

private:
  uhh2::Event::Handle<std::vector<TopJet>> h_ak8jets;
  uhh2::Event::Handle<std::vector<GenTopJet>> h_ak8genjets;
  uhh2::Event::Handle<FlavorParticle> h_primlep;
};

}}
