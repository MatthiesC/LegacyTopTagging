#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Utils.h"

#include "UHH2/common/include/MCWeight.h"

#include <TH1F.h>

using namespace std;
using namespace uhh2;
// using namespace ltt;


namespace uhh2 { namespace ltt {

class TT_Mtt_Hists: public uhh2::Hists {
public:
  TT_Mtt_Hists(uhh2::Context & ctx, const std::string & dirname);
  virtual void fill(const uhh2::Event & event) override;
protected:
  TH1F *hist_mtt;
};


TT_Mtt_Hists::TT_Mtt_Hists(Context & ctx, const string & dirname): Hists(ctx, dirname) {
  hist_mtt = book<TH1F>("mtt", "#it{m}_{t#bar{t}}", 2000, 0, 2000);
}


void TT_Mtt_Hists::fill(const Event & event) {
  const double w = event.weight;
  GenParticle *top(nullptr), *antitop(nullptr);
  for(GenParticle & gp : *event.genparticles) {
    if(top && antitop) {
      break;
    }
    else if(gp.pdgId() == 6) {
      top = &gp;
    }
    else if(gp.pdgId() == -6) {
      antitop = &gp;
    }
  }
  if(!(top && antitop)) cout << "was not able to find top and antitop" << endl;
  hist_mtt->Fill((top->v4() + antitop->v4()).M(), w);
}


class TT_Mtt_Module: public AnalysisModule {
public:
    explicit TT_Mtt_Module(Context & ctx);
    virtual bool process(Event & event) override;

private:
  const bool debug;
  unique_ptr<AnalysisModule> lumi_sf;
  unique_ptr<Hists> mtt_hist;
};


TT_Mtt_Module::TT_Mtt_Module(Context & ctx): debug(string2bool(ctx.get("debug"))) {
  ctx.undeclare_all_event_output();
  lumi_sf.reset(new MCLumiWeight(ctx));
  mtt_hist.reset(new TT_Mtt_Hists(ctx, "TT_Mtt_Hists"));
}


bool TT_Mtt_Module::process(Event & event) {

  if(debug) {
    cout << endl;
    cout << "+-----------+" << endl;
    cout << "| NEW EVENT |" << endl;
    cout << "+-----------+" << endl;
  }

  if(debug) cout << "Lumi Weight" << endl;
  lumi_sf->process(event);

  if(debug) cout << "Filling mtt histogram" << endl;
  mtt_hist->fill(event);

  if(debug) cout << "End of TT_Mtt_Module" << endl;
  return false;
}


UHH2_REGISTER_ANALYSIS_MODULE(TT_Mtt_Module)

}}
