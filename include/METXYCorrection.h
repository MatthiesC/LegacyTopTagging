#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h"

#include "UHH2/common/include/Utils.h"


namespace uhh2 { namespace ltt {

class METXYCorrector: public uhh2::AnalysisModule {
public:
  METXYCorrector(uhh2::Context & ctx, const std::string & met_name = "met", const bool is_puppi = false);
  virtual bool process(uhh2::Event & event) override;
private:
  enum class METType {
    pf,
    puppi,
    notfound,
  };
  uhh2::Event::Handle<MET> fHandleMET;
  const METType fMETType;
  const Year fYear;
  TString fYear_TString;
  bool fIsUL;
};

}}
