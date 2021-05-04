#include "UHH2/common/include/LuminosityHists.h"

#include "UHH2/LegacyTopTagging/include/AndHists.h"
#include "UHH2/LegacyTopTagging/include/AK8Hists.h"
#include "UHH2/LegacyTopTagging/include/CommonHists.h"
#include "UHH2/LegacyTopTagging/include/HOTVRHists.h"

using namespace std;
using namespace uhh2;
using namespace ltt;


namespace uhh2 { namespace ltt {

AndHists::AndHists(Context & ctx, const string & dirname, const bool b_topjethists):
  Hists(ctx, dirname), m_dirname(dirname) {

  hists_vector.push_back(new LuminosityHists(ctx, dirname + "_Lumi"));
  hists_vector.push_back(new ltt::CommonHists(ctx, dirname + "_Common"));
  if(b_topjethists) {
    hists_vector.push_back(new ltt::AK8Hists(ctx, dirname + "_AK8", ctx.get("AK8Collection_rec")));
    hists_vector.push_back(new ltt::HOTVRHists(ctx, dirname + "_HOTVR"));
  }
}


void AndHists::fill(const Event & event) {

  for(Hists *hist : hists_vector) {
    hist->fill(event);
  }
}


void AndHists::add_hist(Hists *hist) {

  hists_vector.push_back(hist);
}

}}
