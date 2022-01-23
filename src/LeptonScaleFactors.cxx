#include "UHH2/LegacyTopTagging/include/LeptonScaleFactors.h"

#include "UHH2/common/include/MCWeight.h"

using namespace std;
using namespace uhh2;
using namespace ltt;


namespace uhh2 { namespace ltt {

//____________________________________________________________________________________________________
ElectronRecoScaleFactors::ElectronRecoScaleFactors(Context & ctx, const bool dummy):
  fDummy(dummy)
{
  const string weight_postfix = "reco";
  const string syst_direction = ctx.get("SystDirection_ElectronReco", "nominal");
  if(fDummy) {
    fHandles.push_back(ctx.declare_event_output<float>((string)"weight_sfelec_"+weight_postfix));
    fHandles.push_back(ctx.declare_event_output<float>((string)"weight_sfelec_"+weight_postfix+"_up"));
    fHandles.push_back(ctx.declare_event_output<float>((string)"weight_sfelec_"+weight_postfix+"_down"));
  }
  else {
    const string base_file_path = ctx.get("uhh2Dir")+"LegacyTopTagging/data/ScaleFactors/egamma_SFs/";
    fYearSwitcher.reset(new YearSwitcher(ctx));
    fYearSwitcher->setupUL16preVFP(
      make_shared<MCElecScaleFactor>(
        ctx,
        base_file_path+"UL16preVFP/egammaEffi_ptAbove20.txt_EGM2D_UL2016preVFP.root",
        0.0,
        weight_postfix,
        syst_direction,
        "electrons",
        "EGamma_SF2D",
        false
      )
    );
    fYearSwitcher->setupUL16postVFP(
      make_shared<MCElecScaleFactor>(
        ctx,
        base_file_path+"UL16postVFP/egammaEffi_ptAbove20.txt_EGM2D_UL2016postVFP.root",
        0.0,
        weight_postfix,
        syst_direction,
        "electrons",
        "EGamma_SF2D",
        false
      )
    );
    fYearSwitcher->setupUL17(
      make_shared<MCElecScaleFactor>(
        ctx,
        base_file_path+"UL17/egammaEffi_ptAbove20.txt_EGM2D_UL2017.root",
        0.0,
        weight_postfix,
        syst_direction,
        "electrons",
        "EGamma_SF2D",
        false
      )
    );
    fYearSwitcher->setupUL18(
      make_shared<MCElecScaleFactor>(
        ctx,
        base_file_path+"UL18/egammaEffi_ptAbove20.txt_EGM2D_UL2018.root",
        0.0,
        weight_postfix,
        syst_direction,
        "electrons",
        "EGamma_SF2D",
        false
      )
    );
  }
}

bool ElectronRecoScaleFactors::process(Event & event) {
  if(fDummy) {
    for(uint i = 0; i < fHandles.size(); i++) {
      event.set(fHandles.at(i), -1.);
    }
  }
  else {
    fYearSwitcher->process(event);
  }
  return true;
}

//____________________________________________________________________________________________________
ElectronIdScaleFactors::ElectronIdScaleFactors(Context & ctx, const Electron::tag & tagID, const bool dummy):
  fDummy(dummy)
{
  const string weight_postfix = "id";
  const string syst_direction = ctx.get("SystDirection_ElectronID", "nominal");
  if(fDummy) {
    fHandles.push_back(ctx.declare_event_output<float>((string)"weight_sfelec_"+weight_postfix));
    fHandles.push_back(ctx.declare_event_output<float>((string)"weight_sfelec_"+weight_postfix+"_up"));
    fHandles.push_back(ctx.declare_event_output<float>((string)"weight_sfelec_"+weight_postfix+"_down"));
  }
  else {
    const string base_file_path = ctx.get("uhh2Dir")+"LegacyTopTagging/data/ScaleFactors/egamma_SFs/";
    string file_name_UL16preVFP;
    string file_name_UL16postVFP;
    string file_name_UL17;
    string file_name_UL18;
    switch(tagID) {
      case Electron::tag::cutBasedElectronID_Fall17_94X_V2_veto :
        file_name_UL16preVFP = "egammaEffi.txt_Ele_Veto_preVFP_EGM2D.root";
        file_name_UL16postVFP = "egammaEffi.txt_Ele_Veto_postVFP_EGM2D.root";
        file_name_UL17 = "egammaEffi.txt_EGM2D_Veto_UL17.root";
        file_name_UL18 = "egammaEffi.txt_Ele_Veto_EGM2D.root";
        break;
      case Electron::tag::cutBasedElectronID_Fall17_94X_V2_loose :
        file_name_UL16preVFP = "egammaEffi.txt_Ele_Loose_preVFP_EGM2D.root";
        file_name_UL16postVFP = "egammaEffi.txt_Ele_Loose_postVFP_EGM2D.root";
        file_name_UL17 = "egammaEffi.txt_EGM2D_Loose_UL17.root";
        file_name_UL18 = "egammaEffi.txt_Ele_Loose_EGM2D.root";
        break;
      case Electron::tag::cutBasedElectronID_Fall17_94X_V2_medium :
        file_name_UL16preVFP = "egammaEffi.txt_Ele_Medium_preVFP_EGM2D.root";
        file_name_UL16postVFP = "egammaEffi.txt_Ele_Medium_postVFP_EGM2D.root";
        file_name_UL17 = "egammaEffi.txt_EGM2D_Medium_UL17.root";
        file_name_UL18 = "egammaEffi.txt_Ele_Medium_EGM2D.root";
        break;
      case Electron::tag::cutBasedElectronID_Fall17_94X_V2_tight :
        file_name_UL16preVFP = "egammaEffi.txt_Ele_Tight_preVFP_EGM2D.root";
        file_name_UL16postVFP = "egammaEffi.txt_Ele_Tight_postVFP_EGM2D.root";
        file_name_UL17 = "egammaEffi.txt_EGM2D_Tight_UL17.root";
        file_name_UL18 = "egammaEffi.txt_Ele_Tight_EGM2D.root";
        break;
      case Electron::tag::mvaEleID_Fall17_noIso_V2_wp90 :
        file_name_UL16preVFP = "egammaEffi.txt_Ele_wp90noiso_preVFP_EGM2D.root";
        file_name_UL16postVFP = "egammaEffi.txt_Ele_wp90noiso_postVFP_EGM2D.root";
        file_name_UL17 = "egammaEffi.txt_EGM2D_MVA90noIso_UL17.root";
        file_name_UL18 = "egammaEffi.txt_Ele_wp90noiso_EGM2D.root";
        break;
      case Electron::tag::mvaEleID_Fall17_noIso_V2_wp80 :
        file_name_UL16preVFP = "egammaEffi.txt_Ele_wp80noiso_preVFP_EGM2D.root";
        file_name_UL16postVFP = "egammaEffi.txt_Ele_wp80noiso_postVFP_EGM2D.root";
        file_name_UL17 = "egammaEffi.txt_EGM2D_MVA80noIso_UL17.root";
        file_name_UL18 = "egammaEffi.txt_Ele_wp80noiso_EGM2D.root";
        break;
      case Electron::tag::mvaEleID_Fall17_iso_V2_wp90 :
        file_name_UL16preVFP = "egammaEffi.txt_Ele_wp90iso_preVFP_EGM2D.root";
        file_name_UL16postVFP = "egammaEffi.txt_Ele_wp90iso_postVFP_EGM2D.root";
        file_name_UL17 = "egammaEffi.txt_EGM2D_MVA90iso_UL17.root";
        file_name_UL18 = "egammaEffi.txt_Ele_wp90iso_EGM2D.root";
        break;
      case Electron::tag::mvaEleID_Fall17_iso_V2_wp80 :
        file_name_UL16preVFP = "egammaEffi.txt_Ele_wp80iso_preVFP_EGM2D.root";
        file_name_UL16postVFP = "egammaEffi.txt_Ele_wp80iso_postVFP_EGM2D.root";
        file_name_UL17 = "egammaEffi.txt_EGM2D_MVA80iso_UL17.root";
        file_name_UL18 = "egammaEffi.txt_Ele_wp80iso_EGM2D.root";
        break;
      default :
        throw invalid_argument("ElectronIdScaleFactors: No scale factors implemented for given electron ID");
    }
    fYearSwitcher.reset(new YearSwitcher(ctx));
    fYearSwitcher->setupUL16preVFP(
      make_shared<MCElecScaleFactor>(
        ctx,
        base_file_path+"UL16preVFP/"+file_name_UL16preVFP,
        0.0,
        weight_postfix,
        syst_direction,
        "electrons",
        "EGamma_SF2D",
        false
      )
    );
    fYearSwitcher->setupUL16postVFP(
      make_shared<MCElecScaleFactor>(
        ctx,
        base_file_path+"UL16postVFP/"+file_name_UL16postVFP,
        0.0,
        weight_postfix,
        syst_direction,
        "electrons",
        "EGamma_SF2D",
        false
      )
    );
    fYearSwitcher->setupUL17(
      make_shared<MCElecScaleFactor>(
        ctx,
        base_file_path+"UL17/"+file_name_UL17,
        0.0,
        weight_postfix,
        syst_direction,
        "electrons",
        "EGamma_SF2D",
        false
      )
    );
    fYearSwitcher->setupUL18(
      make_shared<MCElecScaleFactor>(
        ctx,
        base_file_path+"UL18/"+file_name_UL18,
        0.0,
        weight_postfix,
        syst_direction,
        "electrons",
        "EGamma_SF2D",
        false
      )
    );
  }
}

bool ElectronIdScaleFactors::process(Event & event) {
  if(fDummy) {
    for(uint i = 0; i < fHandles.size(); i++) {
      event.set(fHandles.at(i), -1.);
    }
  }
  else {
    fYearSwitcher->process(event);
  }
  return true;
}

//____________________________________________________________________________________________________
MuonIdScaleFactors::MuonIdScaleFactors(Context & ctx, const Muon::Selector & selectorID, const bool dummy):
  fDummy(dummy)
{
  const string weight_postfix = "id";
  const string syst_direction = ctx.get("SystDirection_MuonID", "nominal");
  if(fDummy) {
    fHandles.push_back(ctx.declare_event_output<float>((string)"weight_sfmu_"+weight_postfix));
    fHandles.push_back(ctx.declare_event_output<float>((string)"weight_sfmu_"+weight_postfix+"_up"));
    fHandles.push_back(ctx.declare_event_output<float>((string)"weight_sfmu_"+weight_postfix+"_down"));
  }
  else {
    const string base_file_path = ctx.get("uhh2Dir")+"LegacyTopTagging/data/ScaleFactors/muon_SFs/";
    string hist_name = "NUM_";
    switch(selectorID) {
      case Muon::Selector::CutBasedIdGlobalHighPt :
        hist_name += "HighPtID";
        break;
      case Muon::Selector::CutBasedIdLoose :
        hist_name += "LooseID";
        break;
      case Muon::Selector::CutBasedIdMedium :
        hist_name += "MediumID";
        break;
      case Muon::Selector::CutBasedIdMediumPrompt :
        hist_name += "MediumPromptID";
        break;
      case Muon::Selector::SoftCutBasedId :
        hist_name += "SoftID";
        break;
      case Muon::Selector::CutBasedIdTight :
        hist_name += "TightID";
        break;
      case Muon::Selector::CutBasedIdTrkHighPt :
        hist_name += "TrkHighPtID";
        break;
      default :
        throw invalid_argument("MuonIdScaleFactors: No scale factors implemented for given muon ID");
    }
    hist_name += "_DEN_TrackerMuons_abseta_pt";
    fYearSwitcher.reset(new YearSwitcher(ctx));
    fYearSwitcher->setupUL16preVFP(
      make_shared<MCMuonScaleFactor>(
        ctx,
        base_file_path+"UL16preVFP/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ID.root",
        hist_name,
        0.0,
        weight_postfix,
        false,
        syst_direction
      )
    );
    fYearSwitcher->setupUL16postVFP(
      make_shared<MCMuonScaleFactor>(
        ctx,
        base_file_path+"UL16postVFP/Efficiencies_muon_generalTracks_Z_Run2016_UL_ID.root",
        hist_name,
        0.0,
        weight_postfix,
        false,
        syst_direction
      )
    );
    fYearSwitcher->setupUL17(
      make_shared<MCMuonScaleFactor>(
        ctx,
        base_file_path+"UL17/Efficiencies_muon_generalTracks_Z_Run2017_UL_ID.root",
        hist_name,
        0.0,
        weight_postfix,
        false,
        syst_direction
      )
    );
    fYearSwitcher->setupUL18(
      make_shared<MCMuonScaleFactor>(
        ctx,
        base_file_path+"UL18/Efficiencies_muon_generalTracks_Z_Run2018_UL_ID.root",
        hist_name,
        0.0,
        weight_postfix,
        false,
        syst_direction
      )
    );
  }
}

bool MuonIdScaleFactors::process(Event & event) {
  if(fDummy) {
    for(uint i = 0; i < fHandles.size(); i++) {
      event.set(fHandles.at(i), -1.);
    }
  }
  else {
    fYearSwitcher->process(event);
  }
  return true;
}

//____________________________________________________________________________________________________
MuonIsoScaleFactors::MuonIsoScaleFactors(Context & ctx, const Muon::Selector & selectorISO, const Muon::Selector & selectorID, const bool dummy):
  fDummy(dummy)
{
  const string weight_postfix = "iso";
  const string syst_direction = ctx.get("SystDirection_MuonIso", "nominal");
  if(fDummy) {
    fHandles.push_back(ctx.declare_event_output<float>((string)"weight_sfmu_"+weight_postfix));
    fHandles.push_back(ctx.declare_event_output<float>((string)"weight_sfmu_"+weight_postfix+"_up"));
    fHandles.push_back(ctx.declare_event_output<float>((string)"weight_sfmu_"+weight_postfix+"_down"));
  }
  else {
    const string base_file_path = ctx.get("uhh2Dir")+"LegacyTopTagging/data/ScaleFactors/muon_SFs/";
    string hist_name = "NUM_";
    if(selectorISO == Muon::Selector::PFIsoLoose && selectorID == Muon::Selector::CutBasedIdLoose) {
      hist_name += "LooseRelIso_DEN_LooseID";
    }
    else if(selectorISO == Muon::Selector::PFIsoLoose && selectorID == Muon::Selector::CutBasedIdMedium) {
      hist_name += "LooseRelIso_DEN_MediumID";
    }
    else if(selectorISO == Muon::Selector::PFIsoLoose && selectorID == Muon::Selector::CutBasedIdMediumPrompt) {
      hist_name += "LooseRelIso_DEN_MediumPromptID";
    }
    else if(selectorISO == Muon::Selector::PFIsoLoose && selectorID == Muon::Selector::CutBasedIdTight) {
      hist_name += "LooseRelIso_DEN_TightIDandIPCut";
    }
    else if(selectorISO == Muon::Selector::TkIsoLoose && selectorID == Muon::Selector::CutBasedIdGlobalHighPt) {
      hist_name += "LooseRelTkIso_DEN_HighPtIDandIPCut";
    }
    else if(selectorISO == Muon::Selector::TkIsoLoose && selectorID == Muon::Selector::CutBasedIdTrkHighPt) {
      hist_name += "LooseRelTkIso_DEN_TrkHighPtIDandIPCut";
    }
    else if(selectorISO == Muon::Selector::PFIsoTight && selectorID == Muon::Selector::CutBasedIdMedium) {
      hist_name += "TightRelIso_DEN_MediumID";
    }
    else if(selectorISO == Muon::Selector::PFIsoTight && selectorID == Muon::Selector::CutBasedIdMediumPrompt) {
      hist_name += "TightRelIso_DEN_MediumPromptID";
    }
    else if(selectorISO == Muon::Selector::PFIsoTight && selectorID == Muon::Selector::CutBasedIdTight) {
      hist_name += "TightRelIso_DEN_TightIDandIPCut";
    }
    else if(selectorISO == Muon::Selector::TkIsoTight && selectorID == Muon::Selector::CutBasedIdGlobalHighPt) {
      hist_name += "TightRelTkIso_DEN_HighPtIDandIPCut";
    }
    else if(selectorISO == Muon::Selector::TkIsoTight && selectorID == Muon::Selector::CutBasedIdTrkHighPt) {
      hist_name += "TightRelTkIso_DEN_TrkHighPtIDandIPCut";
    }
    else {
      throw invalid_argument("MuonIsoScaleFactors: No scale factors implemented for given combination of muon ID + ISO");
    }
    hist_name += "_abseta_pt";
    fYearSwitcher.reset(new YearSwitcher(ctx));
    fYearSwitcher->setupUL16preVFP(
      make_shared<MCMuonScaleFactor>(
        ctx,
        base_file_path+"UL16preVFP/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ISO.root",
        hist_name,
        0.0,
        weight_postfix,
        false,
        syst_direction
      )
    );
    fYearSwitcher->setupUL16postVFP(
      make_shared<MCMuonScaleFactor>(
        ctx,
        base_file_path+"UL16postVFP/Efficiencies_muon_generalTracks_Z_Run2016_UL_ISO.root",
        hist_name,
        0.0,
        weight_postfix,
        false,
        syst_direction
      )
    );
    fYearSwitcher->setupUL17(
      make_shared<MCMuonScaleFactor>(
        ctx,
        base_file_path+"UL17/Efficiencies_muon_generalTracks_Z_Run2017_UL_ISO.root",
        hist_name,
        0.0,
        weight_postfix,
        false,
        syst_direction
      )
    );
    fYearSwitcher->setupUL18(
      make_shared<MCMuonScaleFactor>(
        ctx,
        base_file_path+"UL18/Efficiencies_muon_generalTracks_Z_Run2018_UL_ISO.root",
        hist_name,
        0.0,
        weight_postfix,
        false,
        syst_direction
      )
    );
  }
}

bool MuonIsoScaleFactors::process(Event & event) {
  if(fDummy) {
    for(uint i = 0; i < fHandles.size(); i++) {
      event.set(fHandles.at(i), -1.);
    }
  }
  else {
    fYearSwitcher->process(event);
  }
  return true;
}

//____________________________________________________________________________________________________
MuonTriggerScaleFactors::MuonTriggerScaleFactors(Context & ctx, const bool use_Mu50, const bool absolute_eta, const bool dummy):
  fDummy(dummy)
{
  const string weight_postfix = "trigger";
  const string syst_direction = ctx.get("SystDirection_MuonTrigger", "nominal");
  if(fDummy) {
    fHandles.push_back(ctx.declare_event_output<float>((string)"weight_sfmu_"+weight_postfix));
    fHandles.push_back(ctx.declare_event_output<float>((string)"weight_sfmu_"+weight_postfix+"_up"));
    fHandles.push_back(ctx.declare_event_output<float>((string)"weight_sfmu_"+weight_postfix+"_down"));
  }
  else {
    const string base_file_path = ctx.get("uhh2Dir")+"LegacyTopTagging/data/ScaleFactors/muon_SFs/";
    string hist_name_UL16preVFP;
    string hist_name_UL16postVFP;
    string hist_name_UL17;
    string hist_name_UL18;
    const string hist_name_extension = absolute_eta ? "_abseta_pt" : "_eta_pt";
    switch(use_Mu50) {
      case true :
        hist_name_UL16preVFP = "NUM_Mu50_or_TkMu50_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose";
        hist_name_UL16postVFP = "NUM_Mu50_or_TkMu50_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose";
        hist_name_UL17 = "NUM_Mu50_or_OldMu100_or_TkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose";
        hist_name_UL18 = "NUM_Mu50_or_OldMu100_or_TkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose";
        break;
      case false :
        hist_name_UL16preVFP = "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight";
        hist_name_UL16postVFP = "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight";
        hist_name_UL17 = "NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoTight";
        hist_name_UL18 = "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight";
        break;
    }
    fYearSwitcher.reset(new YearSwitcher(ctx));
    fYearSwitcher->setupUL16preVFP(
      make_shared<MCMuonScaleFactor>(
        ctx,
        base_file_path+"UL16preVFP/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_SingleMuonTriggers.root",
        hist_name_UL16preVFP+hist_name_extension,
        0.0,
        weight_postfix,
        false,
        syst_direction,
        "muons",
        absolute_eta
      )
    );
    fYearSwitcher->setupUL16postVFP(
      make_shared<MCMuonScaleFactor>(
        ctx,
        base_file_path+"UL16postVFP/Efficiencies_muon_generalTracks_Z_Run2016_UL_SingleMuonTriggers.root",
        hist_name_UL16postVFP+hist_name_extension,
        0.0,
        weight_postfix,
        false,
        syst_direction,
        "muons",
        absolute_eta
      )
    );
    fYearSwitcher->setupUL17(
      make_shared<MCMuonScaleFactor>(
        ctx,
        base_file_path+"UL17/Efficiencies_muon_generalTracks_Z_Run2017_UL_SingleMuonTriggers.root",
        hist_name_UL17+hist_name_extension,
        0.0,
        weight_postfix,
        false,
        syst_direction,
        "muons",
        absolute_eta
      )
    );
    fYearSwitcher->setupUL18(
      make_shared<MCMuonScaleFactor>(
        ctx,
        base_file_path+"UL18/Efficiencies_muon_generalTracks_Z_Run2018_UL_SingleMuonTriggers.root",
        hist_name_UL18+hist_name_extension,
        0.0,
        weight_postfix,
        false,
        syst_direction,
        "muons",
        absolute_eta
      )
    );
  }
}

bool MuonTriggerScaleFactors::process(Event & event) {
  if(fDummy) {
    for(uint i = 0; i < fHandles.size(); i++) {
      event.set(fHandles.at(i), -1.);
    }
  }
  else {
    fYearSwitcher->process(event);
  }
  return true;
}

}}
