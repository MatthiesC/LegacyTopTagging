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
#include "UHH2/LegacyTopTagging/include/METXYCorrection.h"
#include "UHH2/LegacyTopTagging/include/LeptonScaleFactors.h"

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
  // unique_ptr<AnalysisModule> sf_muon;
  unique_ptr<AnalysisModule> sf_muon_id;
  unique_ptr<CommonModules> common_modules;
  unique_ptr<Selection> slct_elec;
  unique_ptr<Selection> slct_muon;
  unique_ptr<AnalysisModule> clnr_muon;
  unique_ptr<Selection> slct_met;
  unique_ptr<Selection> slct_ptw;
  unique_ptr<AnalysisModule> primlep;
  unique_ptr<AnalysisModule> met_xy_correction;

  unique_ptr<AndHists> hist_nocuts;
  unique_ptr<AndHists> hist_common;
  unique_ptr<AndHists> hist_muon;
  unique_ptr<AndHists> hist_met_xy_correction;
  unique_ptr<AndHists> hist_met;
  unique_ptr<AndHists> hist_ptw;

  // unique_ptr<Selection> slct_trigger;
};


TagAndProbePreSelectionModule::TagAndProbePreSelectionModule(Context & ctx) {

  debug = string2bool(ctx.get("debug"));

  const JetPFID::wp jetPFID = JetPFID::WP_TIGHT_CHS;
  const JetId jetID = PtEtaCut(30., 2.5);

  // const ElectronId elecID = AndId<Electron>(PtEtaCut(55., 2.4), ElectronID_Fall17_medium_noIso); // used in RunII EOY studies
  // const ElectronId elecID_veto = AndId<Electron>(PtEtaCut(30., 2.4), ElectronID_Fall17_veto_noIso);
  // const MuonId muonID_veto = AndId<Muon>(PtEtaCut(30., 2.4), MuonID(Muon::Selector::CutBasedIdLoose));
  // for more selection efficiency, use medium lepton IDs for the lepton vetoes:
  // const ElectronId elecID_veto = AndId<Electron>(PtEtaCut(30., 2.4), ElectronID_Fall17_medium_noIso); // Christopher's Selection (in 106X_v1)
  const ElectronId elecID_veto = AndId<Electron>(PtEtaCut(30., 2.4), ElectronTagID(Electron::mvaEleID_Fall17_noIso_V2_wpLoose)); // Christopher's Selection (in 106X_v2)
  // const ElectronId elecID_veto = AndId<Electron>(PtEtaCut(55., 2.4), ElectronID_Fall17_medium_noIso);
  const MuonId muonID_veto = AndId<Muon>(PtEtaCut(30., 2.4), MuonID(Muon::Selector::CutBasedIdMedium)); // Christopher's Selection
  // const MuonId muonID_veto = AndId<Muon>(PtEtaCut(55., 2.4), MuonID(Muon::Selector::CutBasedIdTight));

  // const MuonId muonID_tag = AndId<Muon>(PtEtaCut(40., 2.4), MuonID(Muon::Selector::CutBasedIdTight)); // Christopher's Selection
  // pT(tag muon) = 55 GeV should be the final value (will be required in the TagAndProbeMainSelectionModule);
  // in order to plot trigger efficiency histograms starting at values lower than the trigger threshold, we use a looser cut here

  // const Muon::Selector muonSelector = Muon::CutBasedIdTight;
  const Muon::Selector muonSelector = Muon::CutBasedIdGlobalHighPt; // recommended to use this, see https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL2017#High_pT_above_120_GeV (see also AN-2018/008 and in the corresponding paper MUO-17-001)
  // Nachteil dieser ID ist, dass sie nicht auf der loose/medium/tight ID aufbaut, d.h. wenn ich einen Electron channel einbauen sollte, dann kann ich dort nicht durch ein veto auf looseID garantieren, dass events mit GlobalHighPt muon auch vetoed sind !!!
  // Ausserdem wurde die HighPtID benutzt fuer die Bestimmung der HLT_Mu50 trigger-SFs und ist laut oben erwaehntem TWiki fuer die gleichzeitige Benutztung mit diesen sehr gut geeignet
  // Die Benutztung der Mu50 trigger anstatt der IsoMu24/27 etc. trigger in Kombiniation mit der TightID ist sinnvoll, weil die IsoMu24/27 trigger isolierte Myonen verlangen, aber wir wenden eine custom 2D isolation an
  const MuonId muonID_tag = AndId<Muon>(PtEtaCut(55., 2.4), MuonID(muonSelector)); // in order to save disk space, I cut at 55 GeV already in presel!

  slct_lumi.reset(new LumiSelection(ctx));
  sf_lumi.reset(new MCLumiWeight(ctx));
  sf_prefire.reset(new PrefiringWeights(ctx));
  // sf_muon.reset(new MuonScaleFactors(ctx));
  sf_muon_id.reset(new ltt::MuonIdScaleFactors(ctx, muonSelector));

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
  met_xy_correction.reset(new METXYCorrector(ctx));

  hist_nocuts.reset(new AndHists(ctx, "0_NoCuts"));
  hist_common.reset(new AndHists(ctx, "1_Common"));
  hist_muon.reset(new AndHists(ctx, "2_Muon"));
  hist_met_xy_correction.reset(new AndHists(ctx, "3_METXYcorr"));
  hist_met.reset(new AndHists(ctx, "3_MET"));
  hist_ptw.reset(new AndHists(ctx, "4_PtW"));

  // slct_trigger.reset(new LttTriggerSelection(ctx));
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
  // sf_muon->process(event);
  sf_muon_id->process(event);
  hist_muon->fill(event);

  if(debug) cout << "MET XY correction" << endl;
  met_xy_correction->process(event);
  hist_met_xy_correction->fill(event);

  if(debug) cout << "MET selection" << endl;
  if(!slct_met->passes(event)) return false;
  hist_met->fill(event);

  if(debug) cout << "PTW selection" << endl;
  if(!slct_ptw->passes(event)) return false;
  hist_ptw->fill(event);

  // just for trigger debugging
  // if(!slct_trigger->passes(event)) return false;

  if(debug) cout << "End of TagAndProbePreSelectionModule" << endl;
  return true;
}


UHH2_REGISTER_ANALYSIS_MODULE(TagAndProbePreSelectionModule)

}}
