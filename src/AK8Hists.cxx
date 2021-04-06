#include "UHH2/common/include/Utils.h"

#include "UHH2/LegacyTopTagging/include/AK8Hists.h"
#include "UHH2/LegacyTopTagging/include/Utils.h"

#include "TH1F.h"
#include "TH2F.h"

using namespace std;
using namespace uhh2;
using namespace ltt;


namespace uhh2 { namespace ltt {

AK8Hists::AK8Hists(Context & ctx, const string & dirname, const string & ak8_coll_rec): Hists(ctx, dirname) {

    h_ak8_coll_rec = ctx.get_handle<vector<TopJet>>(ak8_coll_rec);

    hist_number = book<TH1F>("number", "Number of AK8 jets", 11, -0.5, 10.5);

    hist_ak8jets_pt = book<TH1F>("ak8jets_pt", "AK8 jets: p_{T} [GeV]", 1000, 0, 1000);
    hist_ak8jets_ptlog = book<TH1F>("ak8jets_ptlog", "AK8 jets: log_{10}(p_{T}) [GeV]", 1000, 0, 4);
    hist_ak8jets_eta = book<TH1F>("ak8jets_eta", "AK8 jets: #eta", 1000, -5.0, 5.0);
    hist_ak8jets_phi = book<TH1F>("ak8jets_phi", "AK8 jets: #phi [rad]", 1000, -M_PI, M_PI);
    hist_ak8jets_mass = book<TH1F>("ak8jets_mass", "AK8 jets: m_{jet} [GeV]", 1000, 0, 500);
    hist_ak8jets_mSD = book<TH1F>("ak8jets_mSD", "AK8 jets: m_{SD} [GeV]", 1000, 0, 500);
    hist_ak8jets_tau32 = book<TH1F>("ak8jets_tau32", "AK8 jets: #tau_{3}/#tau_{2}", 1000, 0, 1);
    hist_ak8jets_maxDeepCSV = book<TH1F>("ak8jets_maxDeepCSV", "AK8 jets: max. #mathcal{O}_{DeepCSV}^{prob(b)+prob(bb)} of subjets", 1000, 0, 1);

    hist_ak8jet1_pt = book<TH1F>("ak8jet1_pt", "Leading AK8 jet: p_{T} [GeV]", 1000, 0, 1000);
    hist_ak8jet1_ptlog = book<TH1F>("ak8jet1_ptlog", "Leading AK8 jet: log_{10}(p_{T}) [GeV]", 1000, 0, 4);
    hist_ak8jet1_eta = book<TH1F>("ak8jet1_eta", "Leading AK8 jet: #eta", 1000, -5.0, 5.0);
    hist_ak8jet1_phi = book<TH1F>("ak8jet1_phi", "Leading AK8 jet: #phi [rad]", 1000, -M_PI, M_PI);
    hist_ak8jet1_mass = book<TH1F>("ak8jet1_mass", "Leading AK8 jet: m_{jet} [GeV]", 1000, 0, 500);
    hist_ak8jet1_mSD = book<TH1F>("ak8jet1_mSD", "Leading AK8 jet: m_{SD} [GeV]", 1000, 0, 500);
    hist_ak8jet1_tau32 = book<TH1F>("ak8jet1_tau32", "Leading AK8 jet: #tau_{3}/#tau_{2}", 1000, 0, 1);
    hist_ak8jet1_maxDeepCSV = book<TH1F>("ak8jet1_maxDeepCSV", "Leading AK8 jet: max. #mathcal{O}_{DeepCSV}^{prob(b)+prob(bb)} of subjets", 1000, 0, 1);

    hist_ak8jet2_pt = book<TH1F>("ak8jet2_pt", "Subleading AK8 jet: p_{T} [GeV]", 1000, 0, 1000);
    hist_ak8jet2_ptlog = book<TH1F>("ak8jet2_ptlog", "Subleading AK8 jet: log_{10}(p_{T}) [GeV]", 1000, 0, 4);
    hist_ak8jet2_eta = book<TH1F>("ak8jet2_eta", "Subleading AK8 jet: #eta", 1000, -5.0, 5.0);
    hist_ak8jet2_phi = book<TH1F>("ak8jet2_phi", "Subleading AK8 jet: #phi [rad]", 1000, -M_PI, M_PI);
    hist_ak8jet2_mass = book<TH1F>("ak8jet2_mass", "Subleading AK8 jet: m_{jet} [GeV]", 1000, 0, 500);
    hist_ak8jet2_mSD = book<TH1F>("ak8jet2_mSD", "Subleading AK8 jet: m_{SD} [GeV]", 1000, 0, 500);
    hist_ak8jet2_tau32 = book<TH1F>("ak8jet2_tau32", "Subleading AK8 jet: #tau_{3}/#tau_{2}", 1000, 0, 1);
    hist_ak8jet2_maxDeepCSV = book<TH1F>("ak8jet2_maxDeepCSV", "Subleading AK8 jet: max. #mathcal{O}_{DeepCSV}^{prob(b)+prob(bb)} of subjets", 1000, 0, 1);
}


void AK8Hists::fill(const Event & event) {

    const double w = event.weight;

    vector<TopJet> ak8jets = event.get(h_ak8_coll_rec);
    sort_by_pt(ak8jets);

    hist_number->Fill(ak8jets.size(), w);

    for(unsigned int i = 0; i < ak8jets.size(); i++) {
        hist_ak8jets_pt->Fill(ak8jets.at(i).v4().Pt(), w);
        hist_ak8jets_ptlog->Fill(log10(ak8jets.at(i).v4().Pt()), w);
        hist_ak8jets_eta->Fill(ak8jets.at(i).v4().Eta(), w);
        hist_ak8jets_phi->Fill(ak8jets.at(i).v4().Phi(), w);
        hist_ak8jets_mass->Fill(ak8jets.at(i).v4().M(), w);
        hist_ak8jets_mSD->Fill(mSD(ak8jets.at(i)), w);
        hist_ak8jets_tau32->Fill(tau32(ak8jets.at(i)), w);
        hist_ak8jets_maxDeepCSV->Fill(maxDeepCSVSubJetValue(ak8jets.at(i)), w);
    }

    if(ak8jets.size() >= 1) {
        hist_ak8jet1_pt->Fill(ak8jets.at(0).v4().Pt(), w);
        hist_ak8jet1_ptlog->Fill(log10(ak8jets.at(0).v4().Pt()), w);
        hist_ak8jet1_eta->Fill(ak8jets.at(0).v4().Eta(), w);
        hist_ak8jet1_phi->Fill(ak8jets.at(0).v4().Phi(), w);
        hist_ak8jet1_mass->Fill(ak8jets.at(0).v4().M(), w);
        hist_ak8jet1_mSD->Fill(mSD(ak8jets.at(0)), w);
        hist_ak8jet1_tau32->Fill(tau32(ak8jets.at(0)), w);
        hist_ak8jet1_maxDeepCSV->Fill(maxDeepCSVSubJetValue(ak8jets.at(0)), w);
    }

    if(ak8jets.size() >= 2) {
        hist_ak8jet2_pt->Fill(ak8jets.at(1).v4().Pt(), w);
        hist_ak8jet2_ptlog->Fill(log10(ak8jets.at(1).v4().Pt()), w);
        hist_ak8jet2_eta->Fill(ak8jets.at(1).v4().Eta(), w);
        hist_ak8jet2_phi->Fill(ak8jets.at(1).v4().Phi(), w);
        hist_ak8jet2_mass->Fill(ak8jets.at(1).v4().M(), w);
        hist_ak8jet2_mSD->Fill(mSD(ak8jets.at(1)), w);
        hist_ak8jet2_tau32->Fill(tau32(ak8jets.at(1)), w);
        hist_ak8jet2_maxDeepCSV->Fill(maxDeepCSVSubJetValue(ak8jets.at(1)), w);
    }
}

}}
