#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"


namespace uhh2 { namespace ltt {

class AK8Hists: public uhh2::Hists {
public:
    AK8Hists(uhh2::Context & ctx, const std::string & dirname, const std::string & coll_rec = "");

    virtual void fill(const uhh2::Event & event) override;

protected:
    TH1F *hist_number;

    TH1F *hist_ak8jets_pt;
    TH1F *hist_ak8jets_ptlog;
    TH1F *hist_ak8jets_eta;
    TH1F *hist_ak8jets_phi;
    TH1F *hist_ak8jets_mass;
    TH1F *hist_ak8jets_mSD;
    TH1F *hist_ak8jets_tau32;
    TH1F *hist_ak8jets_maxDeepCSV;

    TH1F *hist_subjets_number;
    TH1F *hist_subjets_pt;
    TH1F *hist_subjets_eta;
    TH1F *hist_subjets_phi;
    TH1F *hist_subjets_mass;
    TH1F *hist_subjets_deepCSV;

    TH1F *hist_ak8jet1_pt;
    TH1F *hist_ak8jet1_ptlog;
    TH1F *hist_ak8jet1_eta;
    TH1F *hist_ak8jet1_phi;
    TH1F *hist_ak8jet1_mass;
    TH1F *hist_ak8jet1_mSD;
    TH1F *hist_ak8jet1_tau32;
    TH1F *hist_ak8jet1_maxDeepCSV;

    TH1F *hist_ak8jet2_pt;
    TH1F *hist_ak8jet2_ptlog;
    TH1F *hist_ak8jet2_eta;
    TH1F *hist_ak8jet2_phi;
    TH1F *hist_ak8jet2_mass;
    TH1F *hist_ak8jet2_mSD;
    TH1F *hist_ak8jet2_tau32;
    TH1F *hist_ak8jet2_maxDeepCSV;

private:
    uhh2::Event::Handle<std::vector<TopJet>> h_ak8jets;
};

}}
