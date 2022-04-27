#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h"

#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/TTbarGen.h"

#include "UHH2/LegacyTopTagging/include/TopJetCorrections.h"
#include "UHH2/LegacyTopTagging/include/Utils.h"
#include "UHH2/LegacyTopTagging/include/AK8Hists.h"
#include "UHH2/LegacyTopTagging/include/HOTVRHists.h"


using namespace std;
using namespace uhh2;
using namespace ltt;


namespace uhh2 { namespace ltt {

class WorkingPointModule: public AnalysisModule {
public:
    explicit WorkingPointModule(Context & ctx);
    virtual bool process(Event & event) override;

private:
  const bool debug;
  const bool empty_output_tree;
  string dataset_version;
  bool is_qcd, is_ttbar, is_wjets;

  unique_ptr<AnalysisModule> sf_lumi;
  unique_ptr<AnalysisModule> sf_pileup;

  unique_ptr<TopJetCorrections> ak8_corrections;
  unique_ptr<TopJetCorrections> hotvr_corrections;
  unique_ptr<AnalysisModule> ak8_cleaner;
  unique_ptr<AnalysisModule> hotvr_cleaner;
  unique_ptr<Selection> slct_1ak8;
  unique_ptr<Selection> slct_1hotvr;

  unique_ptr<Hists> hists_ak8_before_corrections;
  unique_ptr<Hists> hists_ak8_after_corrections;
  unique_ptr<Hists> hists_ak8_after_cleaning;

  unique_ptr<Hists> hists_hotvr_before_corrections;
  unique_ptr<Hists> hists_hotvr_after_corrections;
  unique_ptr<Hists> hists_hotvr_after_cleaning;

  Event::Handle<float> event_weight;
  Event::Handle<vector<TopJet>> hotvrjets_handle;
  Event::Handle<unsigned int> n_ak8jets;
  Event::Handle<unsigned int> n_hotvrjets;

  // for ttbar

  unique_ptr<AnalysisModule> ttbargen_producer;
  Event::Handle<TTbarGen> ttbargen_handle;

  Event::Handle<float> t_pt;
  Event::Handle<float> t_eta;
  Event::Handle<float> t_phi;
  Event::Handle<float> antit_pt;
  Event::Handle<float> antit_eta;
  Event::Handle<float> antit_phi;
  Event::Handle<float> tt_deta;
  Event::Handle<float> tt_dphi;
  Event::Handle<float> tt_dr;

  Event::Handle<float> wplus_pt;
  Event::Handle<float> wplus_eta;
  Event::Handle<float> wplus_phi;
  Event::Handle<float> wminus_pt;
  Event::Handle<float> wminus_eta;
  Event::Handle<float> wminus_phi;
  Event::Handle<float> ww_deta;
  Event::Handle<float> ww_dphi;
  Event::Handle<float> ww_dr;

  Event::Handle<float> b_pt;
  Event::Handle<float> b_eta;
  Event::Handle<float> b_phi;
  Event::Handle<float> antib_pt;
  Event::Handle<float> antib_eta;
  Event::Handle<float> antib_phi;

  Event::Handle<float> tnearestak8jet_pt;
  Event::Handle<float> tnearestak8jet_msd;
  Event::Handle<float> tnearestak8jet_subjets_deepcsv_max;
  Event::Handle<float> tnearestak8jet_subjets_deepjet_max;
  Event::Handle<float> tnearestak8jet_tau32;
  Event::Handle<float> tnearestak8jet_tau21;
  Event::Handle<float> tnearestak8jet_deepak8_TvsQCD;
  Event::Handle<float> tnearestak8jet_deepak8_WvsQCD;
  Event::Handle<float> tnearestak8jet_partnet_TvsQCD;
  Event::Handle<float> tnearestak8jet_partnet_WvsQCD;
  Event::Handle<float> tnearestak8jet_dr;

  Event::Handle<float> antitnearestak8jet_pt;
  Event::Handle<float> antitnearestak8jet_msd;
  Event::Handle<float> antitnearestak8jet_subjets_deepcsv_max;
  Event::Handle<float> antitnearestak8jet_subjets_deepjet_max;
  Event::Handle<float> antitnearestak8jet_tau32;
  Event::Handle<float> antitnearestak8jet_tau21;
  Event::Handle<float> antitnearestak8jet_deepak8_TvsQCD;
  Event::Handle<float> antitnearestak8jet_deepak8_WvsQCD;
  Event::Handle<float> antitnearestak8jet_partnet_TvsQCD;
  Event::Handle<float> antitnearestak8jet_partnet_WvsQCD;
  Event::Handle<float> antitnearestak8jet_dr;

  Event::Handle<bool> the_two_t_ak8jets_are_the_same;

  Event::Handle<float> wplusnearestak8jet_pt;
  Event::Handle<float> wplusnearestak8jet_msd;
  Event::Handle<float> wplusnearestak8jet_subjets_deepcsv_max;
  Event::Handle<float> wplusnearestak8jet_subjets_deepjet_max;
  Event::Handle<float> wplusnearestak8jet_tau32;
  Event::Handle<float> wplusnearestak8jet_tau21;
  Event::Handle<float> wplusnearestak8jet_deepak8_TvsQCD;
  Event::Handle<float> wplusnearestak8jet_deepak8_WvsQCD;
  Event::Handle<float> wplusnearestak8jet_partnet_TvsQCD;
  Event::Handle<float> wplusnearestak8jet_partnet_WvsQCD;
  Event::Handle<float> wplusnearestak8jet_dr;
  Event::Handle<float> wplusnearestak8jet_dr_b;

  Event::Handle<float> wminusnearestak8jet_pt;
  Event::Handle<float> wminusnearestak8jet_msd;
  Event::Handle<float> wminusnearestak8jet_subjets_deepcsv_max;
  Event::Handle<float> wminusnearestak8jet_subjets_deepjet_max;
  Event::Handle<float> wminusnearestak8jet_tau32;
  Event::Handle<float> wminusnearestak8jet_tau21;
  Event::Handle<float> wminusnearestak8jet_deepak8_TvsQCD;
  Event::Handle<float> wminusnearestak8jet_deepak8_WvsQCD;
  Event::Handle<float> wminusnearestak8jet_partnet_TvsQCD;
  Event::Handle<float> wminusnearestak8jet_partnet_WvsQCD;
  Event::Handle<float> wminusnearestak8jet_dr;
  Event::Handle<float> wminusnearestak8jet_dr_antib;

  Event::Handle<bool> the_two_w_ak8jets_are_the_same;

  Event::Handle<float> tnearesthotvrjet_pt;
  Event::Handle<float> tnearesthotvrjet_mass;
  Event::Handle<int> tnearesthotvrjet_nsubjets;
  Event::Handle<float> tnearesthotvrjet_mpair;
  Event::Handle<float> tnearesthotvrjet_fpt1;
  Event::Handle<float> tnearesthotvrjet_tau32;
  Event::Handle<float> tnearesthotvrjet_dr;
  Event::Handle<float> tnearesthotvrjet_dr_b;
  Event::Handle<float> tnearesthotvrjet_dr_genWplus_d1;
  Event::Handle<float> tnearesthotvrjet_dr_genWplus_d2;

  Event::Handle<float> antitnearesthotvrjet_pt;
  Event::Handle<float> antitnearesthotvrjet_mass;
  Event::Handle<int> antitnearesthotvrjet_nsubjets;
  Event::Handle<float> antitnearesthotvrjet_mpair;
  Event::Handle<float> antitnearesthotvrjet_fpt1;
  Event::Handle<float> antitnearesthotvrjet_tau32;
  Event::Handle<float> antitnearesthotvrjet_dr;
  Event::Handle<float> antitnearesthotvrjet_dr_antib;
  Event::Handle<float> antitnearesthotvrjet_dr_genWminus_d1;
  Event::Handle<float> antitnearesthotvrjet_dr_genWminus_d2;

  Event::Handle<bool> the_two_t_hotvrjets_are_the_same;

  // for wjets

  Event::Handle<float> w_pt;
  Event::Handle<float> w_eta;
  Event::Handle<float> w_phi;

  Event::Handle<float> wnearestak8jet_pt;
  Event::Handle<float> wnearestak8jet_msd;
  Event::Handle<float> wnearestak8jet_tau21;
  Event::Handle<float> wnearestak8jet_deepak8_WvsQCD;
  Event::Handle<float> wnearestak8jet_partnet_WvsQCD;
  Event::Handle<float> wnearestak8jet_dr;

  // for qcd

  Event::Handle<vector<float>> ak8jets_pt;
  Event::Handle<vector<float>> ak8jets_msd;
  Event::Handle<vector<float>> ak8jets_subjets_deepcsv_max;
  Event::Handle<vector<float>> ak8jets_subjets_deepjet_max;
  Event::Handle<vector<float>> ak8jets_tau32;
  Event::Handle<vector<float>> ak8jets_tau21;
  Event::Handle<vector<float>> ak8jets_deepak8_TvsQCD;
  Event::Handle<vector<float>> ak8jets_deepak8_WvsQCD;
  Event::Handle<vector<float>> ak8jets_partnet_TvsQCD;
  Event::Handle<vector<float>> ak8jets_partnet_WvsQCD;

  Event::Handle<vector<float>> hotvrjets_pt;
  Event::Handle<vector<float>> hotvrjets_mass;
  Event::Handle<vector<int>> hotvrjets_nsubjets;
  Event::Handle<vector<float>> hotvrjets_mpair;
  Event::Handle<vector<float>> hotvrjets_fpt1;
  Event::Handle<vector<float>> hotvrjets_tau32;
};


WorkingPointModule::WorkingPointModule(Context & ctx):
  debug(string2bool(ctx.get("debug"))),
  empty_output_tree(string2bool(ctx.get("empty_output_tree")))
{
  dataset_version = ctx.get("dataset_version");
  is_ttbar = (dataset_version.find("TTbar") != string::npos);
  is_wjets = (dataset_version.find("WJets") != string::npos);
  is_qcd = (dataset_version.find("QCD") != string::npos);

  ctx.undeclare_all_event_output();

  sf_lumi.reset(new MCLumiWeight(ctx));
  sf_pileup.reset(new MCPileupReweight(ctx));

  ak8_corrections.reset(new TopJetCorrections());
  ak8_corrections->init(ctx);
  hotvr_corrections.reset(new TopJetCorrections(ctx.get("hotvrCollection_rec"), ctx.get("hotvrCollection_gen")));
  hotvr_corrections->switch_topjet_corrections(false);
  hotvr_corrections->switch_subjet_corrections(true);
  hotvr_corrections->switch_rebuilding_topjets_from_subjets(true);
  hotvr_corrections->init(ctx);
  hotvrjets_handle = ctx.get_handle<vector<TopJet>>(ctx.get("hotvrCollection_rec"));
  const TopJetId ak8_id = AndId<TopJet>(PtEtaCut(200, 2.5), JetPFID(JetPFID::WP_TIGHT_PUPPI));
  const TopJetId hotvr_id = AndId<TopJet>(PtEtaCut(200, 2.5), JetPFID(JetPFID::WP_TIGHT_PUPPI));
  ak8_cleaner.reset(new TopJetCleaner(ctx, ak8_id));
  hotvr_cleaner.reset(new TopJetCleaner(ctx, hotvr_id, ctx.get("hotvrCollection_rec")));
  slct_1ak8.reset(new NTopJetSelection(1, -1));
  slct_1hotvr.reset(new NTopJetSelection(1, -1, boost::none, hotvrjets_handle));

  event_weight = ctx.declare_event_output<float>("event_weight");
  n_ak8jets = ctx.declare_event_output<unsigned int>("n_ak8jets");
  n_hotvrjets = ctx.declare_event_output<unsigned int>("n_hotvrjets");

  if(is_ttbar) {
    ttbargen_producer.reset(new TTbarGenProducer(ctx));
    ttbargen_handle = ctx.get_handle<TTbarGen>("ttbargen");

    t_pt = ctx.declare_event_output<float>("t_pt");
    t_eta = ctx.declare_event_output<float>("t_eta");
    t_phi = ctx.declare_event_output<float>("t_phi");
    antit_pt = ctx.declare_event_output<float>("antit_pt");
    antit_eta = ctx.declare_event_output<float>("antit_eta");
    antit_phi = ctx.declare_event_output<float>("antit_phi");
    tt_deta = ctx.declare_event_output<float>("tt_deta");
    tt_dphi = ctx.declare_event_output<float>("tt_dphi");
    tt_dr = ctx.declare_event_output<float>("tt_dr");

    wplus_pt = ctx.declare_event_output<float>("wplus_pt");
    wplus_eta = ctx.declare_event_output<float>("wplus_eta");
    wplus_phi = ctx.declare_event_output<float>("wplus_phi");
    wminus_pt = ctx.declare_event_output<float>("wminus_pt");
    wminus_eta = ctx.declare_event_output<float>("wminus_eta");
    wminus_phi = ctx.declare_event_output<float>("wminus_phi");
    ww_deta = ctx.declare_event_output<float>("ww_deta");
    ww_dphi = ctx.declare_event_output<float>("ww_dphi");
    ww_dr = ctx.declare_event_output<float>("ww_dr");

    b_pt = ctx.declare_event_output<float>("b_pt");
    b_eta = ctx.declare_event_output<float>("b_eta");
    b_phi = ctx.declare_event_output<float>("b_phi");
    antib_pt = ctx.declare_event_output<float>("antib_pt");
    antib_eta = ctx.declare_event_output<float>("antib_eta");
    antib_phi = ctx.declare_event_output<float>("antib_phi");

    tnearestak8jet_pt = ctx.declare_event_output<float>("tnearestak8jet_pt");
    tnearestak8jet_msd = ctx.declare_event_output<float>("tnearestak8jet_msd");
    tnearestak8jet_subjets_deepcsv_max = ctx.declare_event_output<float>("tnearestak8jet_subjets_deepcsv_max");
    tnearestak8jet_subjets_deepjet_max = ctx.declare_event_output<float>("tnearestak8jet_subjets_deepjet_max");
    tnearestak8jet_tau32 = ctx.declare_event_output<float>("tnearestak8jet_tau32");
    tnearestak8jet_tau21 = ctx.declare_event_output<float>("tnearestak8jet_tau21");
    tnearestak8jet_deepak8_TvsQCD = ctx.declare_event_output<float>("tnearestak8jet_deepak8_TvsQCD");
    tnearestak8jet_deepak8_WvsQCD = ctx.declare_event_output<float>("tnearestak8jet_deepak8_WvsQCD");
    tnearestak8jet_partnet_TvsQCD = ctx.declare_event_output<float>("tnearestak8jet_partnet_TvsQCD");
    tnearestak8jet_partnet_WvsQCD = ctx.declare_event_output<float>("tnearestak8jet_partnet_WvsQCD");
    tnearestak8jet_dr = ctx.declare_event_output<float>("tnearestak8jet_dr");

    antitnearestak8jet_pt = ctx.declare_event_output<float>("antitnearestak8jet_pt");
    antitnearestak8jet_msd = ctx.declare_event_output<float>("antitnearestak8jet_msd");
    antitnearestak8jet_subjets_deepcsv_max = ctx.declare_event_output<float>("antitnearestak8jet_subjets_deepcsv_max");
    antitnearestak8jet_subjets_deepjet_max = ctx.declare_event_output<float>("antitnearestak8jet_subjets_deepjet_max");
    antitnearestak8jet_tau32 = ctx.declare_event_output<float>("antitnearestak8jet_tau32");
    antitnearestak8jet_tau21 = ctx.declare_event_output<float>("antitnearestak8jet_tau21");
    antitnearestak8jet_deepak8_TvsQCD = ctx.declare_event_output<float>("antitnearestak8jet_deepak8_TvsQCD");
    antitnearestak8jet_deepak8_WvsQCD = ctx.declare_event_output<float>("antitnearestak8jet_deepak8_WvsQCD");
    antitnearestak8jet_partnet_TvsQCD = ctx.declare_event_output<float>("antitnearestak8jet_partnet_TvsQCD");
    antitnearestak8jet_partnet_WvsQCD = ctx.declare_event_output<float>("antitnearestak8jet_partnet_WvsQCD");
    antitnearestak8jet_dr = ctx.declare_event_output<float>("antitnearestak8jet_dr");

    the_two_t_ak8jets_are_the_same = ctx.declare_event_output<bool>("the_two_t_ak8jets_are_the_same");

    wplusnearestak8jet_pt = ctx.declare_event_output<float>("wplusnearestak8jet_pt");
    wplusnearestak8jet_msd = ctx.declare_event_output<float>("wplusnearestak8jet_msd");
    wplusnearestak8jet_subjets_deepcsv_max = ctx.declare_event_output<float>("wplusnearestak8jet_subjets_deepcsv_max");
    wplusnearestak8jet_subjets_deepjet_max = ctx.declare_event_output<float>("wplusnearestak8jet_subjets_deepjet_max");
    wplusnearestak8jet_tau32 = ctx.declare_event_output<float>("wplusnearestak8jet_tau32");
    wplusnearestak8jet_tau21 = ctx.declare_event_output<float>("wplusnearestak8jet_tau21");
    wplusnearestak8jet_deepak8_TvsQCD = ctx.declare_event_output<float>("wplusnearestak8jet_deepak8_TvsQCD");
    wplusnearestak8jet_deepak8_WvsQCD = ctx.declare_event_output<float>("wplusnearestak8jet_deepak8_WvsQCD");
    wplusnearestak8jet_partnet_TvsQCD = ctx.declare_event_output<float>("wplusnearestak8jet_partnet_TvsQCD");
    wplusnearestak8jet_partnet_WvsQCD = ctx.declare_event_output<float>("wplusnearestak8jet_partnet_WvsQCD");
    wplusnearestak8jet_dr = ctx.declare_event_output<float>("wplusnearestak8jet_dr");
    wplusnearestak8jet_dr_b = ctx.declare_event_output<float>("wplusnearestak8jet_dr_b");

    wminusnearestak8jet_pt = ctx.declare_event_output<float>("wminusnearestak8jet_pt");
    wminusnearestak8jet_msd = ctx.declare_event_output<float>("wminusnearestak8jet_msd");
    wminusnearestak8jet_subjets_deepcsv_max = ctx.declare_event_output<float>("wminusnearestak8jet_subjets_deepcsv_max");
    wminusnearestak8jet_subjets_deepjet_max = ctx.declare_event_output<float>("wminusnearestak8jet_subjets_deepjet_max");
    wminusnearestak8jet_tau32 = ctx.declare_event_output<float>("wminusnearestak8jet_tau32");
    wminusnearestak8jet_tau21 = ctx.declare_event_output<float>("wminusnearestak8jet_tau21");
    wminusnearestak8jet_deepak8_TvsQCD = ctx.declare_event_output<float>("wminusnearestak8jet_deepak8_TvsQCD");
    wminusnearestak8jet_deepak8_WvsQCD = ctx.declare_event_output<float>("wminusnearestak8jet_deepak8_WvsQCD");
    wminusnearestak8jet_partnet_TvsQCD = ctx.declare_event_output<float>("wminusnearestak8jet_partnet_TvsQCD");
    wminusnearestak8jet_partnet_WvsQCD = ctx.declare_event_output<float>("wminusnearestak8jet_partnet_WvsQCD");
    wminusnearestak8jet_dr = ctx.declare_event_output<float>("wminusnearestak8jet_dr");
    wminusnearestak8jet_dr_antib = ctx.declare_event_output<float>("wminusnearestak8jet_dr_antib");

    the_two_w_ak8jets_are_the_same = ctx.declare_event_output<bool>("the_two_w_ak8jets_are_the_same");

    tnearesthotvrjet_pt = ctx.declare_event_output<float>("tnearesthotvrjet_pt");
    tnearesthotvrjet_mass = ctx.declare_event_output<float>("tnearesthotvrjet_mass");
    tnearesthotvrjet_nsubjets = ctx.declare_event_output<int>("tnearesthotvrjet_nsubjets");
    tnearesthotvrjet_mpair = ctx.declare_event_output<float>("tnearesthotvrjet_mpair");
    tnearesthotvrjet_fpt1 = ctx.declare_event_output<float>("tnearesthotvrjet_fpt1");
    tnearesthotvrjet_tau32 = ctx.declare_event_output<float>("tnearesthotvrjet_tau32");
    tnearesthotvrjet_dr = ctx.declare_event_output<float>("tnearesthotvrjet_dr");
    tnearesthotvrjet_dr_b = ctx.declare_event_output<float>("tnearesthotvrjet_dr_b");
    tnearesthotvrjet_dr_genWplus_d1 = ctx.declare_event_output<float>("tnearesthotvrjet_dr_genWplus_d1");
    tnearesthotvrjet_dr_genWplus_d2 = ctx.declare_event_output<float>("tnearesthotvrjet_dr_genWplus_d2");

    antitnearesthotvrjet_pt = ctx.declare_event_output<float>("antitnearesthotvrjet_pt");
    antitnearesthotvrjet_mass = ctx.declare_event_output<float>("antitnearesthotvrjet_mass");
    antitnearesthotvrjet_nsubjets = ctx.declare_event_output<int>("antitnearesthotvrjet_nsubjets");
    antitnearesthotvrjet_mpair = ctx.declare_event_output<float>("antitnearesthotvrjet_mpair");
    antitnearesthotvrjet_fpt1 = ctx.declare_event_output<float>("antitnearesthotvrjet_fpt1");
    antitnearesthotvrjet_tau32 = ctx.declare_event_output<float>("antitnearesthotvrjet_tau32");
    antitnearesthotvrjet_dr = ctx.declare_event_output<float>("antitnearesthotvrjet_dr");
    antitnearesthotvrjet_dr_antib = ctx.declare_event_output<float>("antitnearesthotvrjet_dr_antib");
    antitnearesthotvrjet_dr_genWminus_d1 = ctx.declare_event_output<float>("antitnearesthotvrjet_dr_genWminus_d1");
    antitnearesthotvrjet_dr_genWminus_d2 = ctx.declare_event_output<float>("antitnearesthotvrjet_dr_genWminus_d2");

    the_two_t_hotvrjets_are_the_same = ctx.declare_event_output<bool>("the_two_t_hotvrjets_are_the_same");
  }
  else if(is_wjets) {
    w_pt = ctx.declare_event_output<float>("w_pt");
    w_eta = ctx.declare_event_output<float>("w_eta");
    w_phi = ctx.declare_event_output<float>("w_phi");

    wnearestak8jet_pt = ctx.declare_event_output<float>("wnearestak8jet_pt");
    wnearestak8jet_msd = ctx.declare_event_output<float>("wnearestak8jet_msd");
    wnearestak8jet_tau21 = ctx.declare_event_output<float>("wnearestak8jet_tau21");
    wnearestak8jet_deepak8_WvsQCD = ctx.declare_event_output<float>("wnearestak8jet_deepak8_WvsQCD");
    wnearestak8jet_partnet_WvsQCD = ctx.declare_event_output<float>("wnearestak8jet_partnet_WvsQCD");
    wnearestak8jet_dr = ctx.declare_event_output<float>("wnearestak8jet_dr");
  }
  else if(is_qcd) {
    ak8jets_pt = ctx.declare_event_output<vector<float>>("ak8jets_pt");
    ak8jets_msd = ctx.declare_event_output<vector<float>>("ak8jets_msd");
    ak8jets_subjets_deepcsv_max = ctx.declare_event_output<vector<float>>("ak8jets_subjets_deepcsv_max");
    ak8jets_subjets_deepjet_max = ctx.declare_event_output<vector<float>>("ak8jets_subjets_deepjet_max");
    ak8jets_tau32 = ctx.declare_event_output<vector<float>>("ak8jets_tau32");
    ak8jets_tau21 = ctx.declare_event_output<vector<float>>("ak8jets_tau21");
    ak8jets_deepak8_TvsQCD = ctx.declare_event_output<vector<float>>("ak8jets_deepak8_TvsQCD");
    ak8jets_deepak8_WvsQCD = ctx.declare_event_output<vector<float>>("ak8jets_deepak8_WvsQCD");
    ak8jets_partnet_TvsQCD = ctx.declare_event_output<vector<float>>("ak8jets_partnet_TvsQCD");
    ak8jets_partnet_WvsQCD = ctx.declare_event_output<vector<float>>("ak8jets_partnet_WvsQCD");

    hotvrjets_pt = ctx.declare_event_output<vector<float>>("hotvrjets_pt");
    hotvrjets_mass = ctx.declare_event_output<vector<float>>("hotvrjets_mass");
    hotvrjets_nsubjets = ctx.declare_event_output<vector<int>>("hotvrjets_nsubjets");
    hotvrjets_mpair = ctx.declare_event_output<vector<float>>("hotvrjets_mpair");
    hotvrjets_fpt1 = ctx.declare_event_output<vector<float>>("hotvrjets_fpt1");
    hotvrjets_tau32 = ctx.declare_event_output<vector<float>>("hotvrjets_tau32");
  }

  hists_ak8_before_corrections.reset(new AK8Hists(ctx, "AK8Hists_0_before_corrections", "", "", "", 1000));
  hists_ak8_after_corrections.reset(new AK8Hists(ctx, "AK8Hists_1_after_corrections", "", "", "", 1000));
  hists_ak8_after_cleaning.reset(new AK8Hists(ctx, "AK8Hists_2_after_cleaning", "", "", "", 1000));

  hists_hotvr_before_corrections.reset(new HOTVRHists(ctx, "HOTVRHists_0_before_corrections", ctx.get("hotvrCollection_rec"), ctx.get("hotvrCollection_gen"), "", 1000));
  hists_hotvr_after_corrections.reset(new HOTVRHists(ctx, "HOTVRHists_1_after_corrections", ctx.get("hotvrCollection_rec"), ctx.get("hotvrCollection_gen"), "", 1000));
  hists_hotvr_after_cleaning.reset(new HOTVRHists(ctx, "HOTVRHists_2_after_cleaning", ctx.get("hotvrCollection_rec"), ctx.get("hotvrCollection_gen"), "", 1000));
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

  hists_ak8_before_corrections->fill(event);
  hists_hotvr_before_corrections->fill(event);
  ak8_corrections->process(event);
  hotvr_corrections->process(event);
  hists_ak8_after_corrections->fill(event);
  hists_hotvr_after_corrections->fill(event);
  ak8_cleaner->process(event);
  hotvr_cleaner->process(event);
  hists_ak8_after_cleaning->fill(event);
  hists_hotvr_after_cleaning->fill(event);

  sort_by_pt(*event.topjets);
  const vector<TopJet> & ak8jets = *event.topjets;
  sort_by_pt(event.get(hotvrjets_handle));
  const vector<TopJet> & hotvrjets = event.get(hotvrjets_handle);

  if(debug) cout << "Set event output" << endl;
  event.set(event_weight, event.weight);
  event.set(n_ak8jets, ak8jets.size());
  event.set(n_hotvrjets, hotvrjets.size());

  if(is_ttbar) {
    ttbargen_producer->process(event);
    TTbarGen *ttbargen = &event.get(ttbargen_handle);
    const GenParticle top = ttbargen->Top();
    const GenParticle genWplus = ttbargen->WTop();
    const GenParticle b = ttbargen->bTop();
    const GenParticle genWplus_d1 = ttbargen->Wdecay1();
    const GenParticle genWplus_d2 = ttbargen->Wdecay2();
    const GenParticle antitop = ttbargen->Antitop();
    const GenParticle genWminus = ttbargen->WAntitop();
    const GenParticle antib = ttbargen->bAntitop();
    const GenParticle genWminus_d1 = ttbargen->WMinusdecay1();
    const GenParticle genWminus_d2 = ttbargen->WMinusdecay2();

    event.set(t_pt, top.v4().Pt());
    event.set(t_eta, top.v4().Eta());
    event.set(t_phi, top.v4().Phi());
    event.set(antit_pt, antitop.v4().Pt());
    event.set(antit_eta, antitop.v4().Eta());
    event.set(antit_phi, antitop.v4().Phi());
    event.set(tt_deta, deltaEta(top.v4(), antitop.v4()));
    event.set(tt_dphi, deltaPhi(top.v4(), antitop.v4()));
    event.set(tt_dr, deltaR(top.v4(), antitop.v4()));

    event.set(wplus_pt, genWplus.v4().Pt());
    event.set(wplus_eta, genWplus.v4().Eta());
    event.set(wplus_phi, genWplus.v4().Phi());
    event.set(wminus_pt, genWminus.v4().Pt());
    event.set(wminus_eta, genWminus.v4().Eta());
    event.set(wminus_phi, genWminus.v4().Phi());
    event.set(ww_deta, deltaEta(genWplus.v4(), genWminus.v4()));
    event.set(ww_dphi, deltaPhi(genWplus.v4(), genWminus.v4()));
    event.set(ww_dr, deltaR(genWplus.v4(), genWminus.v4()));

    event.set(b_pt, b.v4().Pt());
    event.set(b_eta, b.v4().Eta());
    event.set(b_phi, b.v4().Phi());
    event.set(antib_pt, antib.v4().Pt());
    event.set(antib_eta, antib.v4().Eta());
    event.set(antib_phi, antib.v4().Phi());

    float _tnearestak8jet_pt = -1.;
    float _tnearestak8jet_msd = -1.;
    float _tnearestak8jet_subjets_deepcsv_max = -1.;
    float _tnearestak8jet_subjets_deepjet_max = -1.;
    float _tnearestak8jet_tau32 = -1.;
    float _tnearestak8jet_tau21 = -1.;
    float _tnearestak8jet_deepak8_TvsQCD = -1.;
    float _tnearestak8jet_deepak8_WvsQCD = -1.;
    float _tnearestak8jet_partnet_TvsQCD = -1.;
    float _tnearestak8jet_partnet_WvsQCD = -1.;
    float _tnearestak8jet_dr = -1.;

    float _antitnearestak8jet_pt = -1.;
    float _antitnearestak8jet_msd = -1.;
    float _antitnearestak8jet_subjets_deepcsv_max = -1.;
    float _antitnearestak8jet_subjets_deepjet_max = -1.;
    float _antitnearestak8jet_tau32 = -1.;
    float _antitnearestak8jet_tau21 = -1.;
    float _antitnearestak8jet_deepak8_TvsQCD = -1.;
    float _antitnearestak8jet_deepak8_WvsQCD = -1.;
    float _antitnearestak8jet_partnet_TvsQCD = -1.;
    float _antitnearestak8jet_partnet_WvsQCD = -1.;
    float _antitnearestak8jet_dr = -1.;

    bool _the_two_t_ak8jets_are_the_same = false;

    float _wplusnearestak8jet_pt = -1.;
    float _wplusnearestak8jet_msd = -1.;
    float _wplusnearestak8jet_subjets_deepcsv_max = -1.;
    float _wplusnearestak8jet_subjets_deepjet_max = -1.;
    float _wplusnearestak8jet_tau32 = -1.;
    float _wplusnearestak8jet_tau21 = -1.;
    float _wplusnearestak8jet_deepak8_TvsQCD = -1.;
    float _wplusnearestak8jet_deepak8_WvsQCD = -1.;
    float _wplusnearestak8jet_partnet_TvsQCD = -1.;
    float _wplusnearestak8jet_partnet_WvsQCD = -1.;
    float _wplusnearestak8jet_dr = -1.;
    float _wplusnearestak8jet_dr_b = -1.;

    float _wminusnearestak8jet_pt = -1.;
    float _wminusnearestak8jet_msd = -1.;
    float _wminusnearestak8jet_subjets_deepcsv_max = -1.;
    float _wminusnearestak8jet_subjets_deepjet_max = -1.;
    float _wminusnearestak8jet_tau32 = -1.;
    float _wminusnearestak8jet_tau21 = -1.;
    float _wminusnearestak8jet_deepak8_TvsQCD = -1.;
    float _wminusnearestak8jet_deepak8_WvsQCD = -1.;
    float _wminusnearestak8jet_partnet_TvsQCD = -1.;
    float _wminusnearestak8jet_partnet_WvsQCD = -1.;
    float _wminusnearestak8jet_dr = -1.;
    float _wminusnearestak8jet_dr_antib = -1.;

    bool _the_two_w_ak8jets_are_the_same = false;

    if(ak8jets.size() > 0) {
      const TopJet *tnearestak8jet = nextTopJet(top, ak8jets);
      _tnearestak8jet_pt = tnearestak8jet->v4().Pt();
      _tnearestak8jet_msd = mSD(*tnearestak8jet);
      _tnearestak8jet_subjets_deepcsv_max = maxDeepCSVSubJetValue(*tnearestak8jet);
      _tnearestak8jet_subjets_deepjet_max = maxDeepJetSubJetValue(*tnearestak8jet);
      _tnearestak8jet_tau32 = tau32(*tnearestak8jet);
      _tnearestak8jet_tau21 = tau21(*tnearestak8jet);
      _tnearestak8jet_deepak8_TvsQCD = tnearestak8jet->btag_DeepBoosted_TvsQCD();
      _tnearestak8jet_deepak8_WvsQCD = tnearestak8jet->btag_DeepBoosted_WvsQCD();
      _tnearestak8jet_partnet_TvsQCD = tnearestak8jet->btag_ParticleNetDiscriminatorsJetTags_TvsQCD();
      _tnearestak8jet_partnet_WvsQCD = tnearestak8jet->btag_ParticleNetDiscriminatorsJetTags_WvsQCD();
      _tnearestak8jet_dr = deltaR(tnearestak8jet->v4(), top.v4());

      const TopJet *antitnearestak8jet = nextTopJet(antitop, ak8jets);
      _antitnearestak8jet_pt = antitnearestak8jet->v4().Pt();
      _antitnearestak8jet_msd = mSD(*antitnearestak8jet);
      _antitnearestak8jet_subjets_deepcsv_max = maxDeepCSVSubJetValue(*antitnearestak8jet);
      _antitnearestak8jet_subjets_deepjet_max = maxDeepJetSubJetValue(*antitnearestak8jet);
      _antitnearestak8jet_tau32 = tau32(*antitnearestak8jet);
      _antitnearestak8jet_tau21 = tau21(*antitnearestak8jet);
      _antitnearestak8jet_deepak8_TvsQCD = antitnearestak8jet->btag_DeepBoosted_TvsQCD();
      _antitnearestak8jet_deepak8_WvsQCD = antitnearestak8jet->btag_DeepBoosted_WvsQCD();
      _antitnearestak8jet_partnet_TvsQCD = antitnearestak8jet->btag_ParticleNetDiscriminatorsJetTags_TvsQCD();
      _antitnearestak8jet_partnet_WvsQCD = antitnearestak8jet->btag_ParticleNetDiscriminatorsJetTags_WvsQCD();
      _antitnearestak8jet_dr = deltaR(antitnearestak8jet->v4(), antitop.v4());

      _the_two_t_ak8jets_are_the_same = tnearestak8jet == antitnearestak8jet;

      const TopJet *wplusnearestak8jet = nextTopJet(genWplus, ak8jets);
      _wplusnearestak8jet_pt = wplusnearestak8jet->v4().Pt();
      _wplusnearestak8jet_msd = mSD(*wplusnearestak8jet);
      _wplusnearestak8jet_subjets_deepcsv_max = maxDeepCSVSubJetValue(*wplusnearestak8jet);
      _wplusnearestak8jet_subjets_deepjet_max = maxDeepJetSubJetValue(*wplusnearestak8jet);
      _wplusnearestak8jet_tau32 = tau32(*wplusnearestak8jet);
      _wplusnearestak8jet_tau21 = tau21(*wplusnearestak8jet);
      _wplusnearestak8jet_deepak8_TvsQCD = wplusnearestak8jet->btag_DeepBoosted_TvsQCD();
      _wplusnearestak8jet_deepak8_WvsQCD = wplusnearestak8jet->btag_DeepBoosted_WvsQCD();
      _wplusnearestak8jet_partnet_TvsQCD = wplusnearestak8jet->btag_ParticleNetDiscriminatorsJetTags_TvsQCD();
      _wplusnearestak8jet_partnet_WvsQCD = wplusnearestak8jet->btag_ParticleNetDiscriminatorsJetTags_WvsQCD();
      _wplusnearestak8jet_dr = deltaR(wplusnearestak8jet->v4(), genWplus.v4());
      _wplusnearestak8jet_dr_b = deltaR(wplusnearestak8jet->v4(), b.v4());

      const TopJet *wminusnearestak8jet = nextTopJet(genWminus, ak8jets);
      _wminusnearestak8jet_pt = wminusnearestak8jet->v4().Pt();
      _wminusnearestak8jet_msd = mSD(*wminusnearestak8jet);
      _wminusnearestak8jet_subjets_deepcsv_max = maxDeepCSVSubJetValue(*wminusnearestak8jet);
      _wminusnearestak8jet_subjets_deepjet_max = maxDeepJetSubJetValue(*wminusnearestak8jet);
      _wminusnearestak8jet_tau32 = tau32(*wminusnearestak8jet);
      _wminusnearestak8jet_tau21 = tau21(*wminusnearestak8jet);
      _wminusnearestak8jet_deepak8_TvsQCD = wminusnearestak8jet->btag_DeepBoosted_TvsQCD();
      _wminusnearestak8jet_deepak8_WvsQCD = wminusnearestak8jet->btag_DeepBoosted_WvsQCD();
      _wminusnearestak8jet_partnet_TvsQCD = wminusnearestak8jet->btag_ParticleNetDiscriminatorsJetTags_TvsQCD();
      _wminusnearestak8jet_partnet_WvsQCD = wminusnearestak8jet->btag_ParticleNetDiscriminatorsJetTags_WvsQCD();
      _wminusnearestak8jet_dr = deltaR(wminusnearestak8jet->v4(), genWminus.v4());
      _wminusnearestak8jet_dr_antib = deltaR(wminusnearestak8jet->v4(), antib.v4());

      _the_two_w_ak8jets_are_the_same = wplusnearestak8jet == wminusnearestak8jet;
    }

    event.set(tnearestak8jet_pt, _tnearestak8jet_pt);
    event.set(tnearestak8jet_msd, _tnearestak8jet_msd);
    event.set(tnearestak8jet_subjets_deepcsv_max, _tnearestak8jet_subjets_deepcsv_max);
    event.set(tnearestak8jet_subjets_deepjet_max, _tnearestak8jet_subjets_deepjet_max);
    event.set(tnearestak8jet_tau32, _tnearestak8jet_tau32);
    event.set(tnearestak8jet_tau21, _tnearestak8jet_tau21);
    event.set(tnearestak8jet_deepak8_TvsQCD, _tnearestak8jet_deepak8_TvsQCD);
    event.set(tnearestak8jet_deepak8_WvsQCD, _tnearestak8jet_deepak8_WvsQCD);
    event.set(tnearestak8jet_partnet_TvsQCD, _tnearestak8jet_partnet_TvsQCD);
    event.set(tnearestak8jet_partnet_WvsQCD, _tnearestak8jet_partnet_WvsQCD);
    event.set(tnearestak8jet_dr, _tnearestak8jet_dr);

    event.set(antitnearestak8jet_pt, _antitnearestak8jet_pt);
    event.set(antitnearestak8jet_msd, _antitnearestak8jet_msd);
    event.set(antitnearestak8jet_subjets_deepcsv_max, _antitnearestak8jet_subjets_deepcsv_max);
    event.set(antitnearestak8jet_subjets_deepjet_max, _antitnearestak8jet_subjets_deepjet_max);
    event.set(antitnearestak8jet_tau32, _antitnearestak8jet_tau32);
    event.set(antitnearestak8jet_tau21, _antitnearestak8jet_tau21);
    event.set(antitnearestak8jet_deepak8_TvsQCD, _antitnearestak8jet_deepak8_TvsQCD);
    event.set(antitnearestak8jet_deepak8_WvsQCD, _antitnearestak8jet_deepak8_WvsQCD);
    event.set(antitnearestak8jet_partnet_TvsQCD, _antitnearestak8jet_partnet_TvsQCD);
    event.set(antitnearestak8jet_partnet_WvsQCD, _antitnearestak8jet_partnet_WvsQCD);
    event.set(antitnearestak8jet_dr, _antitnearestak8jet_dr);

    event.set(the_two_t_ak8jets_are_the_same, _the_two_t_ak8jets_are_the_same);

    event.set(wplusnearestak8jet_pt, _wplusnearestak8jet_pt);
    event.set(wplusnearestak8jet_msd, _wplusnearestak8jet_msd);
    event.set(wplusnearestak8jet_subjets_deepcsv_max, _wplusnearestak8jet_subjets_deepcsv_max);
    event.set(wplusnearestak8jet_subjets_deepjet_max, _wplusnearestak8jet_subjets_deepjet_max);
    event.set(wplusnearestak8jet_tau32, _wplusnearestak8jet_tau32);
    event.set(wplusnearestak8jet_tau21, _wplusnearestak8jet_tau21);
    event.set(wplusnearestak8jet_deepak8_TvsQCD, _wplusnearestak8jet_deepak8_TvsQCD);
    event.set(wplusnearestak8jet_deepak8_WvsQCD, _wplusnearestak8jet_deepak8_WvsQCD);
    event.set(wplusnearestak8jet_partnet_TvsQCD, _wplusnearestak8jet_partnet_TvsQCD);
    event.set(wplusnearestak8jet_partnet_WvsQCD, _wplusnearestak8jet_partnet_WvsQCD);
    event.set(wplusnearestak8jet_dr, _wplusnearestak8jet_dr);
    event.set(wplusnearestak8jet_dr_b, _wplusnearestak8jet_dr_b);

    event.set(wminusnearestak8jet_pt, _wminusnearestak8jet_pt);
    event.set(wminusnearestak8jet_msd, _wminusnearestak8jet_msd);
    event.set(wminusnearestak8jet_subjets_deepcsv_max, _wminusnearestak8jet_subjets_deepcsv_max);
    event.set(wminusnearestak8jet_subjets_deepjet_max, _wminusnearestak8jet_subjets_deepjet_max);
    event.set(wminusnearestak8jet_tau32, _wminusnearestak8jet_tau32);
    event.set(wminusnearestak8jet_tau21, _wminusnearestak8jet_tau21);
    event.set(wminusnearestak8jet_deepak8_TvsQCD, _wminusnearestak8jet_deepak8_TvsQCD);
    event.set(wminusnearestak8jet_deepak8_WvsQCD, _wminusnearestak8jet_deepak8_WvsQCD);
    event.set(wminusnearestak8jet_partnet_TvsQCD, _wminusnearestak8jet_partnet_TvsQCD);
    event.set(wminusnearestak8jet_partnet_WvsQCD, _wminusnearestak8jet_partnet_WvsQCD);
    event.set(wminusnearestak8jet_dr, _wminusnearestak8jet_dr);
    event.set(wminusnearestak8jet_dr_antib, _wminusnearestak8jet_dr_antib);

    event.set(the_two_w_ak8jets_are_the_same, _the_two_w_ak8jets_are_the_same);

    float _tnearesthotvrjet_pt = -1.;
    float _tnearesthotvrjet_mass = -1.;
    int _tnearesthotvrjet_nsubjets = -1;
    float _tnearesthotvrjet_mpair = -1.;
    float _tnearesthotvrjet_fpt1 = -1.;
    float _tnearesthotvrjet_tau32 = -1.;
    float _tnearesthotvrjet_dr = -1.;
    float _tnearesthotvrjet_dr_b = -1.;
    float _tnearesthotvrjet_dr_genWplus_d1 = -1.;
    float _tnearesthotvrjet_dr_genWplus_d2 = -1.;

    float _antitnearesthotvrjet_pt = -1.;
    float _antitnearesthotvrjet_mass = -1.;
    int _antitnearesthotvrjet_nsubjets = -1;
    float _antitnearesthotvrjet_mpair = -1.;
    float _antitnearesthotvrjet_fpt1 = -1.;
    float _antitnearesthotvrjet_tau32 = -1.;
    float _antitnearesthotvrjet_dr = -1.;
    float _antitnearesthotvrjet_dr_antib = -1.;
    float _antitnearesthotvrjet_dr_genWminus_d1 = -1.;
    float _antitnearesthotvrjet_dr_genWminus_d2 = -1.;

    bool _the_two_t_hotvrjets_are_the_same = false;

    if(hotvrjets.size() > 0) {
      const TopJet *tnearesthotvrjet = nextTopJet(top, hotvrjets);
      _tnearesthotvrjet_pt = tnearesthotvrjet->v4().Pt();
      _tnearesthotvrjet_mass = tnearesthotvrjet->v4().M();
      _tnearesthotvrjet_nsubjets = tnearesthotvrjet->subjets().size();
      _tnearesthotvrjet_mpair = HOTVR_mpair(*tnearesthotvrjet, false);
      _tnearesthotvrjet_fpt1 = HOTVR_fpt(*tnearesthotvrjet);
      _tnearesthotvrjet_tau32 = tau32groomed(*tnearesthotvrjet);
      _tnearesthotvrjet_dr = deltaR(tnearesthotvrjet->v4(), top);
      _tnearesthotvrjet_dr_b = deltaR(tnearesthotvrjet->v4(), b);
      _tnearesthotvrjet_dr_genWplus_d1 = deltaR(tnearesthotvrjet->v4(), genWplus_d1);
      _tnearesthotvrjet_dr_genWplus_d2 = deltaR(tnearesthotvrjet->v4(), genWplus_d2);

      const TopJet *antitnearesthotvrjet = nextTopJet(top, hotvrjets);
      _antitnearesthotvrjet_pt = antitnearesthotvrjet->v4().Pt();
      _antitnearesthotvrjet_mass = antitnearesthotvrjet->v4().M();
      _antitnearesthotvrjet_nsubjets = antitnearesthotvrjet->subjets().size();
      _antitnearesthotvrjet_mpair = HOTVR_mpair(*antitnearesthotvrjet, false);
      _antitnearesthotvrjet_fpt1 = HOTVR_fpt(*antitnearesthotvrjet);
      _antitnearesthotvrjet_tau32 = tau32groomed(*antitnearesthotvrjet);
      _antitnearesthotvrjet_dr = deltaR(antitnearesthotvrjet->v4(), antitop);
      _antitnearesthotvrjet_dr_antib = deltaR(antitnearesthotvrjet->v4(), antib);
      _antitnearesthotvrjet_dr_genWminus_d1 = deltaR(antitnearesthotvrjet->v4(), genWminus_d1);
      _antitnearesthotvrjet_dr_genWminus_d2 = deltaR(antitnearesthotvrjet->v4(), genWminus_d2);

      _the_two_t_hotvrjets_are_the_same = tnearesthotvrjet == antitnearesthotvrjet;
    }

    event.set(tnearesthotvrjet_pt, _tnearesthotvrjet_pt);
    event.set(tnearesthotvrjet_mass, _tnearesthotvrjet_mass);
    event.set(tnearesthotvrjet_nsubjets, _tnearesthotvrjet_nsubjets);
    event.set(tnearesthotvrjet_mpair, _tnearesthotvrjet_mpair);
    event.set(tnearesthotvrjet_fpt1, _tnearesthotvrjet_fpt1);
    event.set(tnearesthotvrjet_tau32, _tnearesthotvrjet_tau32);
    event.set(tnearesthotvrjet_dr, _tnearesthotvrjet_dr);
    event.set(tnearesthotvrjet_dr_b, _tnearesthotvrjet_dr_b);
    event.set(tnearesthotvrjet_dr_genWplus_d1, _tnearesthotvrjet_dr_genWplus_d1);
    event.set(tnearesthotvrjet_dr_genWplus_d2, _tnearesthotvrjet_dr_genWplus_d2);

    event.set(antitnearesthotvrjet_pt, _antitnearesthotvrjet_pt);
    event.set(antitnearesthotvrjet_mass, _antitnearesthotvrjet_mass);
    event.set(antitnearesthotvrjet_nsubjets, _antitnearesthotvrjet_nsubjets);
    event.set(antitnearesthotvrjet_mpair, _antitnearesthotvrjet_mpair);
    event.set(antitnearesthotvrjet_fpt1, _antitnearesthotvrjet_fpt1);
    event.set(antitnearesthotvrjet_tau32, _antitnearesthotvrjet_tau32);
    event.set(antitnearesthotvrjet_dr, _antitnearesthotvrjet_dr);
    event.set(antitnearesthotvrjet_dr_antib, _antitnearesthotvrjet_dr_antib);
    event.set(antitnearesthotvrjet_dr_genWminus_d1, _antitnearesthotvrjet_dr_genWminus_d1);
    event.set(antitnearesthotvrjet_dr_genWminus_d2, _antitnearesthotvrjet_dr_genWminus_d2);

    event.set(the_two_t_hotvrjets_are_the_same, _the_two_t_hotvrjets_are_the_same);
  }
  else if(is_wjets) {
    GenParticle genW;
    for(GenParticle gp : *event.genparticles) {
      if(abs(gp.pdgId()) == 24) genW = gp;
      break;
    }
    event.set(w_pt, genW.v4().Pt());
    event.set(w_eta, genW.v4().Eta());
    event.set(w_phi, genW.v4().Phi());

    float _wnearestak8jet_pt = -1.;
    float _wnearestak8jet_msd = -1.;
    float _wnearestak8jet_tau21 = -1.;
    float _wnearestak8jet_deepak8_WvsQCD = -1.;
    float _wnearestak8jet_partnet_WvsQCD = -1.;
    float _wnearestak8jet_dr = -1.;

    if(ak8jets.size() > 0) {
      const TopJet *wnearestak8jet = nextTopJet(genW, ak8jets);
      _wnearestak8jet_pt = wnearestak8jet->v4().Pt();
      _wnearestak8jet_msd = mSD(*wnearestak8jet);
      _wnearestak8jet_tau21 = tau21(*wnearestak8jet);
      _wnearestak8jet_deepak8_WvsQCD = wnearestak8jet->btag_DeepBoosted_WvsQCD();
      _wnearestak8jet_partnet_WvsQCD = wnearestak8jet->btag_ParticleNetDiscriminatorsJetTags_WvsQCD();
      _wnearestak8jet_dr = deltaR(wnearestak8jet->v4(), genW.v4());
    }

    event.set(wnearestak8jet_pt, _wnearestak8jet_pt);
    event.set(wnearestak8jet_msd, _wnearestak8jet_msd);
    event.set(wnearestak8jet_tau21, _wnearestak8jet_tau21);
    event.set(wnearestak8jet_deepak8_WvsQCD, _wnearestak8jet_deepak8_WvsQCD);
    event.set(wnearestak8jet_partnet_WvsQCD, _wnearestak8jet_partnet_WvsQCD);
    event.set(wnearestak8jet_dr, _wnearestak8jet_dr);
  }
  else if(is_qcd) {
    vector<float> _ak8jets_pt;
    vector<float> _ak8jets_msd;
    vector<float> _ak8jets_subjets_deepcsv_max;
    vector<float> _ak8jets_subjets_deepjet_max;
    vector<float> _ak8jets_tau32;
    vector<float> _ak8jets_tau21;
    vector<float> _ak8jets_deepak8_TvsQCD;
    vector<float> _ak8jets_deepak8_WvsQCD;
    vector<float> _ak8jets_partnet_TvsQCD;
    vector<float> _ak8jets_partnet_WvsQCD;
    for(auto j : ak8jets) {
      _ak8jets_pt.push_back(j.v4().Pt());
      _ak8jets_msd.push_back(mSD(j));
      _ak8jets_subjets_deepcsv_max.push_back(maxDeepCSVSubJetValue(j));
      _ak8jets_subjets_deepjet_max.push_back(maxDeepJetSubJetValue(j));
      _ak8jets_tau32.push_back(tau32(j));
      _ak8jets_tau21.push_back(tau21(j));
      _ak8jets_deepak8_TvsQCD.push_back(j.btag_DeepBoosted_TvsQCD());
      _ak8jets_deepak8_WvsQCD.push_back(j.btag_DeepBoosted_WvsQCD());
      _ak8jets_partnet_TvsQCD.push_back(j.btag_ParticleNetDiscriminatorsJetTags_TvsQCD());
      _ak8jets_partnet_WvsQCD.push_back(j.btag_ParticleNetDiscriminatorsJetTags_WvsQCD());
    }
    event.set(ak8jets_pt, _ak8jets_pt);
    event.set(ak8jets_msd, _ak8jets_msd);
    event.set(ak8jets_subjets_deepcsv_max, _ak8jets_subjets_deepcsv_max);
    event.set(ak8jets_subjets_deepjet_max, _ak8jets_subjets_deepjet_max);
    event.set(ak8jets_tau32, _ak8jets_tau32);
    event.set(ak8jets_tau21, _ak8jets_tau21);
    event.set(ak8jets_deepak8_TvsQCD, _ak8jets_deepak8_TvsQCD);
    event.set(ak8jets_deepak8_WvsQCD, _ak8jets_deepak8_WvsQCD);
    event.set(ak8jets_partnet_TvsQCD, _ak8jets_partnet_TvsQCD);
    event.set(ak8jets_partnet_WvsQCD, _ak8jets_partnet_WvsQCD);

    vector<float> _hotvrjets_pt;
    vector<float> _hotvrjets_mass;
    vector<int> _hotvrjets_nsubjets;
    vector<float> _hotvrjets_mpair;
    vector<float> _hotvrjets_fpt1;
    vector<float> _hotvrjets_tau32;
    for(auto j : hotvrjets) {
      _hotvrjets_pt.push_back(j.v4().Pt());
      _hotvrjets_mass.push_back(j.v4().M());
      _hotvrjets_nsubjets.push_back(j.subjets().size());
      _hotvrjets_mpair.push_back(HOTVR_mpair(j, false));
      _hotvrjets_fpt1.push_back(HOTVR_fpt(j));
      _hotvrjets_tau32.push_back(tau32groomed(j));
    }
    event.set(hotvrjets_pt, _hotvrjets_pt);
    event.set(hotvrjets_mass, _hotvrjets_mass);
    event.set(hotvrjets_nsubjets, _hotvrjets_nsubjets);
    event.set(hotvrjets_mpair, _hotvrjets_mpair);
    event.set(hotvrjets_fpt1, _hotvrjets_fpt1);
    event.set(hotvrjets_tau32, _hotvrjets_tau32);
  }

  if(debug) cout << "End of WorkingPointModule" << endl;
  return !empty_output_tree;
}


UHH2_REGISTER_ANALYSIS_MODULE(WorkingPointModule)

}}
