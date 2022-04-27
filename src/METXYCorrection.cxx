#include "UHH2/LegacyTopTagging/include/METXYCorrection.h"
#include "UHH2/LegacyTopTagging/include/XYMETCorrection_withUL17andUL18andUL16_corrected.h"

using namespace std;
using namespace uhh2;
using namespace ltt;


namespace uhh2 { namespace ltt {

METXYCorrector::METXYCorrector(Context & ctx, const string & met_name, const bool is_puppi):
  fHandleMET(ctx.get_handle<MET>(met_name)),
  // fMETType(uhh2::string2lowercase(ctx.get("METName")).find("puppi") != string::npos ? METType::puppi : METType::pf),
  fMETType(is_puppi ? METType::puppi : METType::pf),
  fYear(extract_year(ctx))
{
  switch(fYear) {
    case Year::is2016v2 :
    fYear_TString = "2016";
    fIsUL = false;
    break;
    case Year::is2016v3 :
    fYear_TString = "2016";
    fIsUL = false;
    break;
    case Year::is2017v1 :
    fYear_TString = "2017";
    fIsUL = false;
    break;
    case Year::is2017v2 :
    fYear_TString = "2017";
    fIsUL = false;
    break;
    case Year::is2018 :
    fYear_TString = "2018";
    fIsUL = false;
    break;
    case Year::isUL16preVFP :
    fYear_TString = "2016APV";
    fIsUL = true;
    break;
    case Year::isUL16postVFP :
    fYear_TString = "2016nonAPV";
    fIsUL = true;
    break;
    case Year::isUL17 :
    fYear_TString = "2017";
    fIsUL = true;
    break;
    case Year::isUL18 :
    fYear_TString = "2018";
    fIsUL = true;
    break;
  }
}

bool METXYCorrector::process(Event & event) {
  MET *met = &event.get(fHandleMET);
  auto new_MET_METphi = METXYCorr_Met_MetPhi(met->pt(), met->phi(), event.run, fYear_TString, !event.isRealData, event.pvs->size(), fIsUL, fMETType == METType::puppi);
  met->set_pt(new_MET_METphi.first);
  met->set_phi(new_MET_METphi.second);
  return true;
}

}}
