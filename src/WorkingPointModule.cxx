#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h"

#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/Utils.h"

#include "UHH2/LegacyTopTagging/include/TopJetCorrections.h"
#include "UHH2/LegacyTopTagging/include/Utils.h"

using namespace std;
using namespace uhh2;


namespace uhh2 { namespace ltt {

class WorkingPointModule: public AnalysisModule {
public:
    explicit WorkingPointModule(Context & ctx);
    virtual bool process(Event & event) override;

private:
  bool debug;
  string dataset_version;
  bool is_qcd, is_ttbar;

  unique_ptr<AnalysisModule> sf_lumi;
  unique_ptr<AnalysisModule> sf_pileup;
  unique_ptr<TopJetCorrections> ak8_corrections;
  unique_ptr<AnalysisModule> ak8_cleaner;
  unique_ptr<Selection> slct_1ak8;

  Event::Handle<float> event_weight;

  // for ttbar

  Event::Handle<float> tnearestjet_pt;
  Event::Handle<float> tnearestjet_msd;
  Event::Handle<float> tnearestjet_subjets_deepcsv_max;
  Event::Handle<float> tnearestjet_tau32;
  Event::Handle<float> tnearestjet_dr;

  Event::Handle<float> antitnearestjet_pt;
  Event::Handle<float> antitnearestjet_msd;
  Event::Handle<float> antitnearestjet_subjets_deepcsv_max;
  Event::Handle<float> antitnearestjet_tau32;
  Event::Handle<float> antitnearestjet_dr;

  Event::Handle<bool> the_two_jets_are_the_same;

  Event::Handle<float> t_pt;
  Event::Handle<float> t_eta;
  Event::Handle<float> t_phi;

  Event::Handle<float> antit_pt;
  Event::Handle<float> antit_eta;
  Event::Handle<float> antit_phi;

  Event::Handle<float> tt_deta;
  Event::Handle<float> tt_dphi;
  Event::Handle<float> tt_dr;

  // for qcd

  Event::Handle<vector<float>> jets_pt;
  Event::Handle<vector<float>> jets_msd;
  Event::Handle<vector<float>> jets_subjets_deepcsv_max;
  Event::Handle<vector<float>> jets_tau32;
};


WorkingPointModule::WorkingPointModule(Context & ctx) {

  debug = string2bool(ctx.get("debug"));
  dataset_version = ctx.get("dataset_version");
  is_ttbar = (dataset_version.find("TTbar") != string::npos);
  is_qcd = (dataset_version.find("QCD") != string::npos);

  ctx.undeclare_all_event_output();

  sf_lumi.reset(new MCLumiWeight(ctx));
  sf_pileup.reset(new MCPileupReweight(ctx));
  ak8_corrections.reset(new TopJetCorrections());
  ak8_corrections->init(ctx);
  TopJetId ak8_id = PtEtaCut(300, 2.4);
  ak8_cleaner.reset(new TopJetCleaner(ctx, ak8_id));
  slct_1ak8.reset(new NTopJetSelection(1, -1));

  event_weight = ctx.declare_event_output<float>("event_weight");

  if(is_ttbar) {
    tnearestjet_pt = ctx.declare_event_output<float>("tnearestjet_pt");
    tnearestjet_msd = ctx.declare_event_output<float>("tnearestjet_msd");
    tnearestjet_subjets_deepcsv_max = ctx.declare_event_output<float>("tnearestjet_subjets_deepcsv_max");
    tnearestjet_tau32 = ctx.declare_event_output<float>("tnearestjet_tau32");
    tnearestjet_dr = ctx.declare_event_output<float>("tnearestjet_dr");

    antitnearestjet_pt = ctx.declare_event_output<float>("antitnearestjet_pt");
    antitnearestjet_msd = ctx.declare_event_output<float>("antitnearestjet_msd");
    antitnearestjet_subjets_deepcsv_max = ctx.declare_event_output<float>("antitnearestjet_subjets_deepcsv_max");
    antitnearestjet_tau32 = ctx.declare_event_output<float>("antitnearestjet_tau32");
    antitnearestjet_dr = ctx.declare_event_output<float>("antitnearestjet_dr");

    the_two_jets_are_the_same = ctx.declare_event_output<bool>("the_two_jets_are_the_same");

    t_pt = ctx.declare_event_output<float>("t_pt");
    t_eta = ctx.declare_event_output<float>("t_eta");
    t_phi = ctx.declare_event_output<float>("t_phi");

    antit_pt = ctx.declare_event_output<float>("antit_pt");
    antit_eta = ctx.declare_event_output<float>("antit_eta");
    antit_phi = ctx.declare_event_output<float>("antit_phi");

    tt_deta = ctx.declare_event_output<float>("tt_deta");
    tt_dphi = ctx.declare_event_output<float>("tt_dphi");
    tt_dr = ctx.declare_event_output<float>("tt_dr");
  }
  else if(is_qcd) {
    jets_pt = ctx.declare_event_output<vector<float>>("jets_pt");
    jets_msd = ctx.declare_event_output<vector<float>>("jets_msd");
    jets_subjets_deepcsv_max = ctx.declare_event_output<vector<float>>("jets_subjets_deepcsv_max");
    jets_tau32 = ctx.declare_event_output<vector<float>>("jets_tau32");
  }
}


bool WorkingPointModule::process(Event & event) {

  if(debug) {
    cout << endl;
    cout << "+-----------+" << endl;
    cout << "| NEW EVENT |" << endl;
    cout << "+-----------+" << endl;
  }

  sf_lumi->process(event);
  sf_pileup->process(event);
  ak8_corrections->process(event);
  ak8_cleaner->process(event);
  if(!slct_1ak8->passes(event)) return false;
  sort_by_pt(*event.topjets);
  const vector<TopJet> & topjets = *event.topjets;

  if(debug) cout << "Set event output" << endl;
  event.set(event_weight, event.weight);
  if(is_ttbar) {
    GenParticle top, antitop;
    for(GenParticle gp : *event.genparticles) {
      if(gp.pdgId() == 6) top = gp;
      else if(gp.pdgId() == -6) antitop = gp;
    }
    event.set(t_pt, top.v4().Pt());
    event.set(t_eta, top.v4().Eta());
    event.set(t_phi, top.v4().Phi());
    event.set(antit_pt, antitop.v4().Pt());
    event.set(antit_eta, antitop.v4().Eta());
    event.set(antit_phi, antitop.v4().Phi());
    event.set(tt_deta, deltaEta(top.v4(), antitop.v4()));
    event.set(tt_dphi, deltaPhi(top.v4(), antitop.v4()));
    event.set(tt_dr, deltaR(top.v4(), antitop.v4()));

    const TopJet *tnearestjet = nextTopJet(top, topjets);
    event.set(tnearestjet_pt, tnearestjet->v4().Pt());
    event.set(tnearestjet_msd, mSD(*tnearestjet));
    event.set(tnearestjet_subjets_deepcsv_max, maxDeepCSVSubJetValue(*tnearestjet));
    event.set(tnearestjet_tau32, tau32(*tnearestjet));
    event.set(tnearestjet_dr, deltaR(tnearestjet->v4(), top.v4()));

    const TopJet *antitnearestjet = nextTopJet(antitop, topjets);
    event.set(antitnearestjet_pt, antitnearestjet->v4().Pt());
    event.set(antitnearestjet_msd, mSD(*antitnearestjet));
    event.set(antitnearestjet_subjets_deepcsv_max, maxDeepCSVSubJetValue(*antitnearestjet));
    event.set(antitnearestjet_tau32, tau32(*antitnearestjet));
    event.set(antitnearestjet_dr, deltaR(antitnearestjet->v4(), antitop.v4()));

    event.set(the_two_jets_are_the_same, tnearestjet == antitnearestjet);
  }
  else if(is_qcd) {
    vector<float> _jets_pt;
    vector<float> _jets_msd;
    vector<float> _jets_subjets_deepcsv_max;
    vector<float> _jets_tau32;
    for(auto j : topjets) {
      _jets_pt.push_back(j.v4().Pt());
      _jets_msd.push_back(mSD(j));
      _jets_subjets_deepcsv_max.push_back(maxDeepCSVSubJetValue(j));
      _jets_tau32.push_back(tau32(j));
    }
    event.set(jets_pt, _jets_pt);
    event.set(jets_msd, _jets_msd);
    event.set(jets_subjets_deepcsv_max, _jets_subjets_deepcsv_max);
    event.set(jets_tau32, _jets_tau32);
  }

  if(debug) cout << "End of WorkingPointModule" << endl;
  return true;
}


UHH2_REGISTER_ANALYSIS_MODULE(WorkingPointModule)

}}
