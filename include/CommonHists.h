#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"


namespace uhh2 { namespace ltt {

class CommonHists: public uhh2::Hists {
public:
  CommonHists(uhh2::Context & ctx, const std::string & dirname);

  virtual void fill(const uhh2::Event & event) override;

protected:
  TH1F *hist_count;
  TH1F *hist_weights_log10;

  TH1F *hist_ntrue;
  TH1F *hist_npv;

  TH1F *hist_met_pt;
  TH1F *hist_met_phi;

  TH1F *hist_nmuons;
  TH1F *hist_nelectrons;

  TH1F *hist_primlep_pt;
  TH1F *hist_ptw;
  TH1F *hist_primlep_eta;
  TH1F *hist_primlep_phi;
  TH1F *hist_primlep_drjet;
  TH1F *hist_primlep_ptrel;
  TH2F *hist_twodselection;

  TH1F *hist_nak4;
  TH1F *hist_ak4jets_pt;
  TH1F *hist_ak4jets_eta;
  TH1F *hist_ak4jets_phi;
  TH1F *hist_ak4jets_mass;
  TH1F *hist_ak4jets_drlepton;
  TH1F *hist_ak4jets_deepCSV;
  TH1F *hist_ak4jets_deepJet;

  TH1F *hist_ht_had;
  TH1F *hist_ht_lep;

private:
  uhh2::Event::Handle<FlavorParticle> h_primlep;
};

}}
