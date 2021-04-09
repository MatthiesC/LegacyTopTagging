#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"


namespace uhh2 { namespace ltt {

class HOTVRHists: public uhh2::Hists {
public:
  HOTVRHists(uhh2::Context & ctx, const std::string & dirname, const std::string & coll_rec = "");

  virtual void fill(const uhh2::Event & event) override;

protected:
  TH1F *hist_number;

  TH1F *hist_hotvrjets_pt;
  TH1F *hist_hotvrjets_ptlog;
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
  TH1F *hist_hotvrjet1_ptlog;
  TH1F *hist_hotvrjet1_eta;
  TH1F *hist_hotvrjet1_phi;
  TH1F *hist_hotvrjet1_mass;
  TH1F *hist_hotvrjet1_mpair;
  TH1F *hist_hotvrjet1_tau32;
  TH1F *hist_hotvrjet1_fpt1;

  TH1F *hist_hotvrjet2_pt;
  TH1F *hist_hotvrjet2_ptlog;
  TH1F *hist_hotvrjet2_eta;
  TH1F *hist_hotvrjet2_phi;
  TH1F *hist_hotvrjet2_mass;
  TH1F *hist_hotvrjet2_mpair;
  TH1F *hist_hotvrjet2_tau32;
  TH1F *hist_hotvrjet2_fpt1;

private:
  uhh2::Event::Handle<std::vector<TopJet>> h_hotvrjets;
};

}}
