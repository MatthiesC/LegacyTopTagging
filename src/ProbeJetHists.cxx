#include "UHH2/common/include/Utils.h"

#include "UHH2/LegacyTopTagging/include/ProbeJetHists.h"
#include "UHH2/LegacyTopTagging/include/Utils.h"

#include "TH1F.h"

using namespace std;
using namespace uhh2;
using namespace ltt;


namespace uhh2 { namespace ltt {

AK8ProbeJetHists::AK8ProbeJetHists(Context & ctx, const string & dirname, const MergeScenario & _msc): Hists(ctx, dirname), msc(_msc) {

  h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
  h_probejet = ctx.get_handle<TopJet>("ProbeJetAK8");
  h_merge_scenario = ctx.get_handle<MergeScenario>("output_merge_scenario_AK8");

  for(const auto & pt_bin : pt_bin_map) {
    const string & pt_bin_string = kPtBinAsString.at(pt_bin.first);
    for(const auto & jet_cat : kJetCategoryAsString) {
      const string & jet_cat_string = jet_cat.second;
      for(const auto & wp : wp_map) {
        const string & wp_string = kWorkingPointAsString.at(wp.first);
        for(const auto & pass_cat : kPassCategoryAsString) {
          const string & pass_cat_string = pass_cat.second;
          vector<TH1F*> hists;
          const string & prefix = pt_bin_string+"_"+jet_cat_string+"_"+wp_string+"_"+pass_cat_string+"_";
          hists.push_back(book<TH1F>((prefix+"pt").c_str(), "Probe jet #it{p}_{T} [GeV]", 1000, 0, 1000));
          hists.push_back(book<TH1F>((prefix+"drlepton").c_str(), "#Delta#it{R}(probe jet, lepton)", 1000, 0, 5));
          hists.push_back(book<TH1F>((prefix+"eta").c_str(), "Probe jet #eta", 1000, -5.0, 5.0));
          hists.push_back(book<TH1F>((prefix+"phi").c_str(), "Probe jet #phi [rad]", 1000, -M_PI, M_PI));
          hists.push_back(book<TH1F>((prefix+"mass").c_str(), "Probe jet #it{m}_{jet} [GeV]", 1000, 0, 500));
          hists.push_back(book<TH1F>((prefix+"mSD").c_str(), "Probe jet #it{m}_{SD} [GeV]", 1000, 0, 500));
          hists.push_back(book<TH1F>((prefix+"tau32").c_str(), "Probe jet #tau_{3}/#tau_{2}", 1000, 0, 1));
          hists.push_back(book<TH1F>((prefix+"maxDeepCSV").c_str(), "Max. #it{O}_{DeepCSV}^{prob(b)+prob(bb)} of probe subjets", 1000, 0, 1));
          hists_map[pt_bin.first][jet_cat.first][wp.first][pass_cat.first] = hists;
        }
      }
    }
  }
}

void AK8ProbeJetHists::fill_probe(const vector<TH1F*> & hists) {

  unsigned int i(0);
  hists.at(i++)->Fill(probejet.v4().pt(), w);
  hists.at(i++)->Fill(deltaR(probejet.v4(), primlep.v4()), w);
  hists.at(i++)->Fill(probejet.v4().eta(), w);
  hists.at(i++)->Fill(probejet.v4().phi(), w);
  hists.at(i++)->Fill(probejet.v4().M(), w);
  hists.at(i++)->Fill(mSD(probejet), w);
  hists.at(i++)->Fill(tau32(probejet), w);
  hists.at(i++)->Fill(maxDeepCSVSubJetValue(probejet), w);
  if(i != hists.size()) throw runtime_error("AK8ProbeJetHists::fill_probe(): Number of declared and filled histograms do not match. Please check!");
}

void AK8ProbeJetHists::fill(const Event & event) {

  if(!event.is_valid(h_probejet)) return;
  if(event.get(h_merge_scenario) != msc) return;

  w = event.weight;
  primlep = event.get(h_primlep);
  probejet = event.get(h_probejet);

  const bool passes_mass_cut = (mSD(probejet) > 105. && mSD(probejet) < 210.);
  bool passes_btagging(false);
  for(const auto & subjet : probejet.subjets()) {
    if(SubjetBTagID(subjet, event)) {
      passes_btagging = true;
      break;
    }
  }

  for(const auto & pt_bin : pt_bin_map) {
    if(probejet.v4().pt() < pt_bin.second.first || probejet.v4().pt() > pt_bin.second.second) continue;
    for(const auto & wp : wp_map) {
      const bool passes_tau32_cut = tau32(probejet) < wp.second;

      if(passes_tau32_cut) fill_probe(hists_map[pt_bin.first][JetCategory::All][wp.first][PassCategory::Pass]);
      else fill_probe(hists_map[pt_bin.first][JetCategory::All][wp.first][PassCategory::Fail]);

      if(passes_mass_cut) {
        if(passes_tau32_cut) fill_probe(hists_map[pt_bin.first][JetCategory::Mass][wp.first][PassCategory::Pass]);
        else fill_probe(hists_map[pt_bin.first][JetCategory::Mass][wp.first][PassCategory::Fail]);
      }

      if(passes_btagging) {
        if(passes_tau32_cut) fill_probe(hists_map[pt_bin.first][JetCategory::BTag][wp.first][PassCategory::Pass]);
        else fill_probe(hists_map[pt_bin.first][JetCategory::BTag][wp.first][PassCategory::Fail]);
      }

      if(passes_mass_cut && passes_btagging) {
        if(passes_tau32_cut) fill_probe(hists_map[pt_bin.first][JetCategory::MassAndBTag][wp.first][PassCategory::Pass]);
        else fill_probe(hists_map[pt_bin.first][JetCategory::MassAndBTag][wp.first][PassCategory::Fail]);
      }
    }
  }
}

HOTVRProbeJetHists::HOTVRProbeJetHists(Context & ctx, const string & dirname, const MergeScenario & _msc): Hists(ctx, dirname), msc(_msc) {

  h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
  h_probejet = ctx.get_handle<TopJet>("ProbeJetHOTVR");
  h_merge_scenario = ctx.get_handle<MergeScenario>("output_merge_scenario_HOTVR");

  for(const auto & pt_bin : pt_bin_map) {
    const string & pt_bin_string = kPtBinAsString.at(pt_bin.first);
    for(const auto & jet_cat : kJetCategoryAsString) {
      const string & jet_cat_string = jet_cat.second;
      for(const auto & wp : wp_map) {
        const string & wp_string = kWorkingPointAsString.at(wp.first);
        for(const auto & pass_cat : kPassCategoryAsString) {
          const string & pass_cat_string = pass_cat.second;
          vector<TH1F*> hists;
          const string & prefix = pt_bin_string+"_"+jet_cat_string+"_"+wp_string+"_"+pass_cat_string+"_";
          hists.push_back(book<TH1F>((prefix+"pt").c_str(), "Probe jet #it{p}_{T} [GeV]", 1000, 0, 1000));
          hists.push_back(book<TH1F>((prefix+"drlepton").c_str(), "#Delta#it{R}(probe jet, lepton)", 1000, 0, 5));
          hists.push_back(book<TH1F>((prefix+"eta").c_str(), "Probe jet #eta", 1000, -5.0, 5.0));
          hists.push_back(book<TH1F>((prefix+"phi").c_str(), "Probe jet #phi [rad]", 1000, -M_PI, M_PI));
          hists.push_back(book<TH1F>((prefix+"mass").c_str(), "Probe jet #it{m}_{jet} [GeV]", 1000, 0, 500));
          hists.push_back(book<TH1F>((prefix+"mpair").c_str(), "Min. #it{m}_{ij} [GeV] of leading three probe subjets", 1000, 0, 250));
          hists.push_back(book<TH1F>((prefix+"tau32").c_str(), "Probe jet #tau_{3}/#tau_{2}", 1000, 0, 1));
          hists.push_back(book<TH1F>((prefix+"fpt1").c_str(), "#it{p}_{T} fraction of leading probe subjet", 1000, 0, 1));
          hists_map[pt_bin.first][jet_cat.first][wp.first][pass_cat.first] = hists;
        }
      }
    }
  }
}

void HOTVRProbeJetHists::fill_probe(const vector<TH1F*> & hists) {

  unsigned int i(0);
  hists.at(i++)->Fill(probejet.v4().pt(), w);
  hists.at(i++)->Fill(deltaR(probejet.v4(), primlep.v4()), w);
  hists.at(i++)->Fill(probejet.v4().eta(), w);
  hists.at(i++)->Fill(probejet.v4().phi(), w);
  hists.at(i++)->Fill(probejet.v4().M(), w);
  hists.at(i++)->Fill(HOTVR_mpair(probejet, false), w);
  hists.at(i++)->Fill(tau32groomed(probejet), w);
  hists.at(i++)->Fill(HOTVR_fpt(probejet), w);
  if(i != hists.size()) throw runtime_error("HOTVRProbeJetHists::fill_probe(): Number of declared and filled histograms do not match. Please check!");
}

void HOTVRProbeJetHists::fill(const Event & event) {

  if(!event.is_valid(h_probejet)) return;
  if(event.get(h_merge_scenario) != msc) return;

  w = event.weight;
  primlep = event.get(h_primlep);
  probejet = event.get(h_probejet);

  const bool passes_mass_cut = (probejet.v4().M() > 140. && probejet.v4().M() < 220.);
  const bool passes_hotvr_cuts = HOTVRTopTagID(probejet, event);

  for(const auto & pt_bin : pt_bin_map) {
    if(probejet.v4().pt() < pt_bin.second.first || probejet.v4().pt() > pt_bin.second.second) continue;
    for(const auto & wp : wp_map) {
      const bool passes_tau32_cut = tau32groomed(probejet) < wp.second;

      if(passes_tau32_cut) fill_probe(hists_map[pt_bin.first][JetCategory::All][wp.first][PassCategory::Pass]);
      else fill_probe(hists_map[pt_bin.first][JetCategory::All][wp.first][PassCategory::Fail]);

      if(passes_mass_cut) {
        if(passes_tau32_cut) fill_probe(hists_map[pt_bin.first][JetCategory::Mass][wp.first][PassCategory::Pass]);
        else fill_probe(hists_map[pt_bin.first][JetCategory::Mass][wp.first][PassCategory::Fail]);
      }

      if(passes_hotvr_cuts) {
        if(passes_tau32_cut) fill_probe(hists_map[pt_bin.first][JetCategory::HOTVRCuts][wp.first][PassCategory::Pass]);
        else fill_probe(hists_map[pt_bin.first][JetCategory::HOTVRCuts][wp.first][PassCategory::Fail]);
      }

      if(passes_mass_cut && passes_hotvr_cuts) {
        if(passes_tau32_cut) fill_probe(hists_map[pt_bin.first][JetCategory::HOTVRCutsAndMass][wp.first][PassCategory::Pass]);
        else fill_probe(hists_map[pt_bin.first][JetCategory::HOTVRCutsAndMass][wp.first][PassCategory::Fail]);
      }
    }
  }
}

ProbeJetHistsRunner::ProbeJetHistsRunner(Context & ctx, const string & dirname): Hists(ctx, dirname) {

  for(const auto & msc : kMergeScenarioAsString) {
    hists_vector.push_back(new ltt::AK8ProbeJetHists(ctx, dirname+"_AK8_"+msc.second, msc.first));
    hists_vector.push_back(new ltt::HOTVRProbeJetHists(ctx, dirname+"_HOTVR_"+msc.second, msc.first));
  }
}

void ProbeJetHistsRunner::fill(const Event & event) {

  for(Hists *hist : hists_vector) {
    hist->fill(event);
  }
}

}}
