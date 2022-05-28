#include "UHH2/common/include/TriggerSelection.h"

#include "UHH2/LegacyTopTagging/include/TriggerSelection.h"

#include "TRandom3.h"

using namespace std;
using namespace uhh2;
using namespace uhh2::ltt;

namespace uhh2 { namespace ltt {
//____________________________________________________________________________________________________
// Information on EGamma triggers taken from https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTRunIISummary
MyTriggerSelection::MyTriggerSelection(Context & ctx, const bool low_pt):
  fYear(extract_year(ctx)),
  fChannel(extract_channel(ctx)),
  fLowPt(low_pt)
{
  const TString dataset_version = ctx.get("dataset_version");
  if(dataset_version.Contains("SingleMuon")) fDataStream = DataStream::isSingleMuon;
  else if(dataset_version.Contains("SingleElectron")) fDataStream = DataStream::isSingleElectron;
  else if(dataset_version.Contains("SinglePhoton")) fDataStream = DataStream::isSinglePhoton;
  else if(dataset_version.Contains("EGamma")) fDataStream = DataStream::isEGamma;
  else fDataStream = DataStream::isMC;

  fTrigSel_IsoMu24.reset(new TriggerSelection("HLT_IsoMu24_v*"));
  fTrigSel_IsoTkMu24.reset(new TriggerSelection("HLT_IsoTkMu24_v*"));
  fTrigSel_IsoMu27.reset(new TriggerSelection("HLT_IsoMu27_v*"));

  fTrigSel_Mu50.reset(new TriggerSelection("HLT_Mu50_v*"));
  fTrigSel_TkMu50.reset(new TriggerSelection("HLT_TkMu50_v*"));
  fTrigSel_OldMu100.reset(new TriggerSelection("HLT_OldMu100_v*"));
  fTrigSel_TkMu100.reset(new TriggerSelection("HLT_TkMu100_v*"));

  fTrigSel_Ele27_WPTight_Gsf.reset(new TriggerSelection("HLT_Ele27_WPTight_Gsf_v*"));
  fTrigSel_Ele35_WPTight_Gsf.reset(new TriggerSelection("HLT_Ele35_WPTight_Gsf_v*"));
  fTrigSel_Ele32_WPTight_Gsf.reset(new TriggerSelection("HLT_Ele32_WPTight_Gsf_v*"));

  fTrigSel_Ele115_CaloIdVT_GsfTrkIdT.reset(new TriggerSelection("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*"));
  fTrigSel_Photon175.reset(new TriggerSelection("HLT_Photon175_v*"));
  fTrigSel_Photon200.reset(new TriggerSelection("HLT_Photon200_v*"));
}

bool MyTriggerSelection::passes(const Event & event) {
  if(fDataStream == DataStream::isMC && event.isRealData) throw runtime_error("BTWTriggerSelection::passes(): Conflict with event.isRealData and dataset_version");
  if(fYear == Year::isUL16preVFP || fYear == Year::isUL16postVFP) {
    if(fChannel == Channel::isMuo) {
      if(fYear == Year::isUL16postVFP) {
        if(fLowPt) return fTrigSel_IsoMu24->passes(event) || fTrigSel_IsoTkMu24->passes(event);
        else return fTrigSel_Mu50->passes(event) || fTrigSel_TkMu50->passes(event);
      }
      else if(fYear == Year::isUL16preVFP) {
        if(fDataStream == DataStream::isMC) {
          if(fLowPt) return fTrigSel_IsoMu24->passes(event) || fTrigSel_IsoTkMu24->passes(event);
          else {
            // Use random number generator with eta-dependent seed for reproducibility
            TRandom3 *random = new TRandom3((int)(fabs(event.muons->at(0).v4().eta()*1000))); // TRandom3 is the same generator as used for gRandom
            if(random->Rndm() > lumi_percentage_UL16preVFP_without_TkMu50) { // emulation of UL16preVFP Run B with run >= 274889
              return fTrigSel_Mu50->passes(event) || fTrigSel_TkMu50->passes(event);
            }
            else { // emulation of UL16preVFP Run B with run < 274889
              return fTrigSel_Mu50->passes(event);
            }
          }
        }
        else if(fDataStream == DataStream::isSingleMuon) {
          if(fLowPt) return fTrigSel_IsoMu24->passes(event) || fTrigSel_IsoTkMu24->passes(event);
          else {
            if(is_UL16preVFP_without_TkMu50(event)) return fTrigSel_Mu50->passes(event);
            else return fTrigSel_Mu50->passes(event) || fTrigSel_TkMu50->passes(event);
          }
        }
        else return false;
      }
      else return false;
    }
    else if(fChannel == Channel::isEle) {
      if(fDataStream == DataStream::isMC) {
        if(fLowPt) return fTrigSel_Ele27_WPTight_Gsf->passes(event);
        else return fTrigSel_Ele115_CaloIdVT_GsfTrkIdT->passes(event) || fTrigSel_Photon175->passes(event);
      }
      else if(fDataStream == DataStream::isSingleElectron) {
        if(fLowPt) return fTrigSel_Ele27_WPTight_Gsf->passes(event);
        else return fTrigSel_Ele115_CaloIdVT_GsfTrkIdT->passes(event);
      }
      else if(fDataStream == DataStream::isSinglePhoton) {
        // Veto Ele115 since those events will be in the SingleElectron stream already
        if(fLowPt) return false;
        else return !fTrigSel_Ele115_CaloIdVT_GsfTrkIdT->passes(event) && fTrigSel_Photon175->passes(event);
      }
      else return false;
    }
    else return false;
  }
  else if(fYear == Year::isUL17) {
    if(fChannel == Channel::isMuo) {
      // if(fLowPt) return fTrigSel_IsoMu27->passes(event);
      // else return fTrigSel_Mu50->passes(event) || fTrigSel_OldMu100->passes(event) || fTrigSel_TkMu100->passes(event);
      if(fDataStream == DataStream::isMC) {
        if(fLowPt) return fTrigSel_IsoMu27->passes(event);
        else {
          // Use random number generator with eta-dependent seed for reproducibility
          TRandom3 *random = new TRandom3((int)(fabs(event.muons->at(0).v4().eta()*1000))); // TRandom3 is the same generator as used for gRandom
          if(random->Rndm() > lumi_percentage_UL17_RunB) { // Run C-F emulation
            return fTrigSel_Mu50->passes(event) || fTrigSel_OldMu100->passes(event) || fTrigSel_TkMu100->passes(event);
          }
          else { // Run B emulation
            return fTrigSel_Mu50->passes(event);
          }
        }
      }
      else if(fDataStream == DataStream::isSingleMuon) {
        if(fLowPt) return fTrigSel_IsoMu27->passes(event);
        else {
          if(is_UL17_RunB(event)) return fTrigSel_Mu50->passes(event);
          else return fTrigSel_Mu50->passes(event) || fTrigSel_OldMu100->passes(event) || fTrigSel_TkMu100->passes(event);
        }
      }
      else return false;
    }
    else if(fChannel == Channel::isEle) {
      if(fDataStream == DataStream::isMC) {
        // Use random number generator with eta-dependent seed for reproducibility
        TRandom3 *random = new TRandom3((int)(fabs(event.electrons->at(0).v4().eta()*1000))); // TRandom3 is the same generator as used for gRandom
        if(random->Rndm() > lumi_percentage_UL17_RunB) { // Run C-F emulation
          if(fLowPt) return fTrigSel_Ele35_WPTight_Gsf->passes(event);
          else return fTrigSel_Ele115_CaloIdVT_GsfTrkIdT->passes(event) || fTrigSel_Photon200->passes(event);
        }
        else { // Run B emulation
          return fTrigSel_Ele35_WPTight_Gsf->passes(event) || fTrigSel_Photon200->passes(event);
        }
      }
      else if(fDataStream == DataStream::isSingleElectron) {
        if(is_UL17_RunB(event)) {
          return fTrigSel_Ele35_WPTight_Gsf->passes(event);
        }
        else {
          if(fLowPt) return fTrigSel_Ele35_WPTight_Gsf->passes(event);
          else return fTrigSel_Ele115_CaloIdVT_GsfTrkIdT->passes(event);
        }
      }
      else if(fDataStream == DataStream::isSinglePhoton) {
        if(is_UL17_RunB(event)) {
          // Veto Ele35 since those events will be in the SingleElectron stream already
          return !fTrigSel_Ele35_WPTight_Gsf->passes(event) && fTrigSel_Photon200->passes(event);
        }
        else {
          // Veto Ele115 since those events will be in the SingleElectron stream already
          if(fLowPt) return false;
          else return !fTrigSel_Ele115_CaloIdVT_GsfTrkIdT->passes(event) && fTrigSel_Photon200->passes(event);
        }
      }
      else return false;
    }
    else return false;
  }
  else if(fYear == Year::isUL18)  {
    if(fChannel == Channel::isMuo) {
      if(fLowPt) return fTrigSel_IsoMu24->passes(event);
      else return fTrigSel_Mu50->passes(event) || fTrigSel_OldMu100->passes(event) || fTrigSel_TkMu100->passes(event);
    }
    else if(fChannel == Channel::isEle) {
      // No need for differentiation between SingleElectron and SinglePhoton streams since we have EGamma in 2018
      // According to https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTRunIISummary#2018 there is no need for a photon trigger
      if(fLowPt) return fTrigSel_Ele32_WPTight_Gsf->passes(event);
      else return fTrigSel_Ele115_CaloIdVT_GsfTrkIdT->passes(event) || fTrigSel_Photon200->passes(event);
    }
    else return false;
  }
  else return false;
}

}}
