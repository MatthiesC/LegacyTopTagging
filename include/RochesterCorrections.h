#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Utils.h"

#include "UHH2/common/include/Utils.h"

#include "UHH2/LegacyTopTagging/include/RoccoR.h"


namespace uhh2 { namespace ltt {

/*
https://gitlab.cern.ch/akhukhun/roccor
https://twiki.cern.ch/twiki/bin/viewauth/CMS/RochcorMuon
https://indico.cern.ch/event/981770/contributions/4135530/attachments/2157765/3639739/roccor.pdf

Cite: EPJC V72, 10.2194 (2012) (arXiv:1208.3710)
*/
class RochesterCorrections: public uhh2::AnalysisModule {
public:
  RochesterCorrections(uhh2::Context & ctx);
  virtual bool process(uhh2::Event & event) override;
private:
  const Year fYear;
  std::unique_ptr<RoccoR> rc;
  int variation_set;
  int variation_member;
};

}}
