#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h"

#include "UHH2/common/include/LumiSelection.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/PrimaryLepton.h"
#include "UHH2/common/include/CleaningModules.h"

#include "UHH2/LegacyTopTagging/include/Utils.h"
#include "UHH2/LegacyTopTagging/include/AndHists.h"

using namespace std;
using namespace uhh2;
using namespace ltt;


namespace uhh2 { namespace ltt {

class TagAndProbePreSelectionModule: public AnalysisModule {
public:
    explicit TagAndProbePreSelectionModule(Context & ctx);
    virtual bool process(Event & event) override;

private:
  bool debug;

  unique_ptr<Selection> slct_lumi;
  unique_ptr<AnalysisModule> sf_lumi;
  unique_ptr<AnalysisModule> sf_prefire;
  unique_ptr<AnalysisModule> sf_muon;
  unique_ptr<CommonModules> common_modules;
  unique_ptr<Selection> slct_elec;
  unique_ptr<Selection> slct_muon;
  unique_ptr<AnalysisModule> clnr_muon;
  unique_ptr<Selection> slct_met;
  unique_ptr<Selection> slct_ptw;
  unique_ptr<AnalysisModule> primlep;

  unique_ptr<AndHists> hist_nocuts;
  unique_ptr<AndHists> hist_common;
  unique_ptr<AndHists> hist_muon;
  unique_ptr<AndHists> hist_met;
  unique_ptr<AndHists> hist_ptw;
};


TagAndProbePreSelectionModule::TagAndProbePreSelectionModule(Context & ctx) {

  debug = string2bool(ctx.get("debug"));

  const JetPFID::wp jetPFID = JetPFID::WP_TIGHT_CHS;
  const JetId jetID = PtEtaCut(30., 2.4);

  // const ElectronId elecID = AndId<Electron>(PtEtaCut(55., 2.4), ElectronID_Fall17_medium_noIso); // used in RunII EOY studies
  // const ElectronId elecID_veto = AndId<Electron>(PtEtaCut(30., 2.4), ElectronID_Fall17_veto_noIso);
  // const MuonId muonID_veto = AndId<Muon>(PtEtaCut(30., 2.4), MuonID(Muon::Selector::CutBasedIdLoose));
  // for more selection efficiency, use medium lepton IDs for the lepton vetoes:
  const ElectronId elecID_veto = AndId<Electron>(PtEtaCut(30., 2.4), ElectronID_Fall17_medium_noIso); // Christopher's Selection
  // const ElectronId elecID_veto = AndId<Electron>(PtEtaCut(55., 2.4), ElectronID_Fall17_medium_noIso);
  const MuonId muonID_veto = AndId<Muon>(PtEtaCut(30., 2.4), MuonID(Muon::Selector::CutBasedIdMedium)); // Christopher's Selection
  // const MuonId muonID_veto = AndId<Muon>(PtEtaCut(55., 2.4), MuonID(Muon::Selector::CutBasedIdTight));

  const MuonId muonID_tag = AndId<Muon>(PtEtaCut(40., 2.4), MuonID(Muon::Selector::CutBasedIdTight)); // Christopher's Selection
  // pT(tag muon) = 55 GeV should be the final value (will be required in the TagAndProbeMainSelectionModule);
  // in order to plot trigger efficiency histograms starting at values lower than the trigger threshold, we use a looser cut here

  slct_lumi.reset(new LumiSelection(ctx));
  sf_lumi.reset(new MCLumiWeight(ctx));
  sf_prefire.reset(new PrefiringWeights(ctx));
  sf_muon.reset(new MuonScaleFactors(ctx));

  common_modules.reset(new CommonModules());
  common_modules->change_pf_id(jetPFID);
  common_modules->set_jet_id(jetID);
  common_modules->set_muon_id(muonID_veto);
  common_modules->set_electron_id(elecID_veto);
  common_modules->switch_metcorrection(true);
  common_modules->switch_jetlepcleaner(true);
  common_modules->switch_jetPtSorter(true);
  common_modules->disable_lumisel();
  common_modules->disable_mclumiweight();
  common_modules->init(ctx);

  slct_elec.reset(new NElectronSelection(0, 0));
  slct_muon.reset(new NMuonSelection(1, 1));
  clnr_muon.reset(new MuonCleaner(muonID_tag));
  slct_met.reset(new METSelection(50.));
  slct_ptw.reset(new PTWSelection(ctx, 150.));

  primlep.reset(new PrimaryLepton(ctx));

  hist_nocuts.reset(new AndHists(ctx, "0_NoCuts"));
  hist_common.reset(new AndHists(ctx, "1_Common"));
  hist_muon.reset(new AndHists(ctx, "2_Muon"));
  hist_met.reset(new AndHists(ctx, "3_MET"));
  hist_ptw.reset(new AndHists(ctx, "4_PtW"));
}


bool TagAndProbePreSelectionModule::process(Event & event) {

  if(debug) {
    cout << endl;
    cout << "+-----------+" << endl;
    cout << "| NEW EVENT |" << endl;
    cout << "+-----------+" << endl;
  }

  if(debug) cout << "Lumi selection (need to do this manually before CommonModules)" << endl; // else getting error for some data samples, e.g. "RunSwitcher cannot handle run number 275656 for year 2016"
  if(event.isRealData && !slct_lumi->passes(event)) return false;
  sf_lumi->process(event);
  sf_prefire->process(event);
  hist_nocuts->fill(event);

  if(debug) cout << "CommonModules: jet, electron, muon id cleaning; jet-lepton cleaning; MET+PV filter; AK4+MET corrections; PU weights" << endl;
  if(!common_modules->process(event)) return false;
  hist_common->fill(event);

  if(debug) cout << "Lepton selection" << endl;
  if(!slct_elec->passes(event)) return false; // veto vetoID electrons
  if(!slct_muon->passes(event)) return false; // veto additional looseID muons
  clnr_muon->process(event);
  if(!slct_muon->passes(event)) return false; // require exactly one tightID muon
  primlep->process(event);
  sf_muon->process(event);
  hist_muon->fill(event);

  if(debug) cout << "MET selection" << endl;
  if(!slct_met->passes(event)) return false;
  hist_met->fill(event);

  if(debug) cout << "PTW selection" << endl;
  if(!slct_ptw->passes(event)) return false;
  hist_ptw->fill(event);

  if(debug) cout << "End of TagAndProbePreSelectionModule" << endl;
  return true;
}


UHH2_REGISTER_ANALYSIS_MODULE(TagAndProbePreSelectionModule)

}}
