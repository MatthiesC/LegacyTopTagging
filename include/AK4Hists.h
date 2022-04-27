#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"


namespace uhh2 { namespace ltt {

class AK4Hists: public uhh2::Hists {
public:
  AK4Hists(uhh2::Context & ctx, const std::string & dirname, const unsigned int default_nbins = 100);

  virtual void fill(const uhh2::Event & event) override;

protected:
  TH1F *hist_number_puppijets;
  TH1F *hist_number_puppijets_uncleaned;
  TH1F *hist_number_puppijets_central;
  TH1F *hist_number_puppijets_forward;

  TH1F *hist_puppichs_dr;
  TH1F *hist_puppichs_ptresponse;

  TH1F *hist_number_bjets;
  TH1F *hist_number_bjets_loose;
  TH1F *hist_number_bjets_medium;
  TH1F *hist_number_bjets_tight;

  TH2F *hist_number_bjets_loose_medium;
  TH2F *hist_number_bjets_loose_tight;
  TH2F *hist_number_bjets_medium_tight;

  TH1F *hist_paired_puppijets_pt;
  TH1F *hist_paired_puppijets_drlepton;
  TH1F *hist_paired_puppijets_eta;
  TH1F *hist_paired_puppijets_phi;
  TH1F *hist_paired_puppijets_mass;
  TH1F *hist_paired_puppijets_deepjet;


  TH1F *hist_puppijets_pt;
  TH1F *hist_forward_puppijets_pt;

  TH1F *hist_puppijets_drlepton;
  TH1F *hist_forward_puppijets_drlepton;

  TH1F *hist_puppijets_eta;
  TH1F *hist_forward_puppijets_eta;

  TH1F *hist_puppijets_phi;
  TH1F *hist_forward_puppijets_phi;

  TH1F *hist_puppijets_mass;
  TH1F *hist_forward_puppijets_mass;


  TH1F *hist_paired_puppijet1pt_pt;
  TH1F *hist_paired_puppijet1dj_pt;

  TH1F *hist_paired_puppijet1pt_drlepton;
  TH1F *hist_paired_puppijet1dj_drlepton;

  TH1F *hist_paired_puppijet1pt_eta;
  TH1F *hist_paired_puppijet1dj_eta;

  TH1F *hist_paired_puppijet1pt_phi;
  TH1F *hist_paired_puppijet1dj_phi;

  TH1F *hist_paired_puppijet1pt_mass;
  TH1F *hist_paired_puppijet1dj_mass;

  TH1F *hist_paired_puppijet1pt_deepjet;
  TH1F *hist_paired_puppijet1dj_deepjet;

private:
  const uhh2::Event::Handle<std::vector<Jet>> fHandle_PUPPIjets;
  const uhh2::Event::Handle<std::vector<Jet>> fHandle_CHSjets;
  const uhh2::Event::Handle<std::vector<Jet>> fHandle_pairedPUPPIjets;
  const uhh2::Event::Handle<std::vector<Jet>> fHandle_forwardPUPPIjets;
  const uhh2::Event::Handle<std::vector<Jet>> fHandle_uncleanedPUPPIjets;
  const uhh2::Event::Handle<std::vector<Jet>> fHandle_bJets;
  const uhh2::Event::Handle<std::vector<Jet>> fHandle_bJets_loose;
  const uhh2::Event::Handle<std::vector<Jet>> fHandle_bJets_medium;
  const uhh2::Event::Handle<std::vector<Jet>> fHandle_bJets_tight;
  const uhh2::Event::Handle<FlavorParticle> fHandle_PrimaryLepton;
};

}}
