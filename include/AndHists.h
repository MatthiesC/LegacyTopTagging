#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"


namespace uhh2 { namespace ltt {

class AndHists: public uhh2::Hists {
public:
  AndHists(uhh2::Context & ctx, const std::string & dirname, const bool b_topjethists = false, const bool b_ak4hists = false);
  virtual void fill(const uhh2::Event & event) override;
  void add_hist(uhh2::Hists *hist);
  std::string dirname() const { return m_dirname; };

private:
  const std::string m_dirname;

protected:
  std::vector<uhh2::Hists*> hists_vector;
};

}}
