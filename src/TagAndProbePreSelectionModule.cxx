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
#include "UHH2/common/include/TriggerSelection.h"

#include "UHH2/LegacyTopTagging/include/Utils.h"

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
  unique_ptr<CommonModules> common_modules;
  unique_ptr<Selection> slct_elec;
  unique_ptr<Selection> slct_muon;
  unique_ptr<Selection> slct_met;
  unique_ptr<Selection> slct_ptw;
  unique_ptr<Selection> slct_ak4;
  unique_ptr<Selection> slct_trigger;
  unique_ptr<AnalysisModule> primlep;
};


TagAndProbePreSelectionModule::TagAndProbePreSelectionModule(Context & ctx) {

  debug = string2bool(ctx.get("debug"));

  JetPFID::wp jetPFID = JetPFID::WP_TIGHT_CHS;
  JetId jetID = PtEtaCut(30., 2.4);
  ElectronId elecID = AndId<Electron>(PtEtaCut(55., 2.4), ElectronID_Fall17_medium_noIso);
  MuonId muonID = AndId<Muon>(PtEtaCut(55., 2.4), MuonID(Muon::Selector::CutBasedIdTight));

  slct_lumi.reset(new LumiSelection(ctx));
  sf_lumi.reset(new MCLumiWeight(ctx));

  common_modules.reset( new CommonModules());
  common_modules->change_pf_id(jetPFID);
  common_modules->set_jet_id(jetID);
  common_modules->set_muon_id(muonID);
  common_modules->set_electron_id(elecID);
  common_modules->switch_metcorrection(true);
  common_modules->switch_jetlepcleaner(true);
  common_modules->switch_jetPtSorter(true);
  common_modules->disable_lumisel();
  common_modules->disable_mclumiweight();
  common_modules->init(ctx);

  slct_elec.reset(new NElectronSelection(0, 0));
  slct_muon.reset(new NMuonSelection(1, 1));
  slct_met.reset(new METSelection(50.));
  slct_ptw.reset(new PTWSelection(ctx, 150.));
  slct_ak4.reset(new NJetSelection(2));
  slct_trigger.reset(new TriggerSelection("HLT_Mu50_v*"));

  primlep.reset(new PrimaryLepton(ctx));
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

  if(debug) cout << "CommonModules: jet, electron, muon id cleaning; jet-lepton cleaning; MET+PV filter; AK4+MET corrections; PU weights" << endl;
  if(!common_modules->process(event)) return false;

  if(debug) cout << "Lepton selection" << endl;
  if(!slct_elec->passes(event)) return false;
  if(!slct_muon->passes(event)) return false;
  primlep->process(event);
  // lepton SF

  if(debug) cout << "MET selection" << endl;
  if(!slct_met->passes(event)) return false;

  if(debug) cout << "PTW selection" << endl;
  if(!slct_ptw->passes(event)) return false;

  if(debug) cout << "AK4 selection" << endl;
  if(!slct_ak4->passes(event)) return false;

  if(debug) cout << "Trigger selection" << endl;
  if(!slct_trigger->passes(event)) return false;
  // trigger SF

  if(debug) cout << "End of TagAndProbePreSelectionModule" << endl;
  return true;
}


UHH2_REGISTER_ANALYSIS_MODULE(TagAndProbePreSelectionModule)

}}
