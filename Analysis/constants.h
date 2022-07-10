#include <set>
#include <string>

namespace macros {

static const float kGoldenRatio = 0.618033;

static const map<string, int> kColors = {
  // Colors close to matplotlib's default color cycle:
  {"pyplot_orange", kOrange + 7},
  {"pyplot_blue", kAzure + 9},
};

enum class Channel {
  isEle,
  isMuo,
  notValid,
  isBoth,
};

typedef struct {
  int index;
  std::string name;
  std::string tlatex;
} ChannelInfo;

const std::map<Channel, ChannelInfo> kChannels = {
  { Channel::isEle,    ChannelInfo{1, "ele", "e + jets"} },
  { Channel::isMuo,    ChannelInfo{2, "muo", "#mu + jets"} },
  { Channel::notValid, ChannelInfo{0, "invalid"} },
  { Channel::isBoth, ChannelInfo{3, "both", "e/#mu + jets"} },
};

enum class Year {
  isUL16preVFP,
  isUL16postVFP,
  isUL17,
  isUL18,

  isUL16,
  isRun2,
};

typedef struct {
  int index;
  std::string name;
  std::string nice_name;
  float lumi_fb;
  float lumi_pb;
  float lumi_unc;
  std::string lumi_fb_display;
  int tcolor;
  int linestyle;
} YearInfo;

const std::map<Year, YearInfo> kYears = {
  {Year::isUL16preVFP, YearInfo{1, "UL16preVFP", "2016 (early)", 19.536, 19536., 0.012, "19.5", 900+8, 1} },
  {Year::isUL16postVFP, YearInfo{2, "UL16postVFP", "2016 (late)", 16.797, 16797., 0.012, "16.8", 860-6, 2} },
  {Year::isUL17, YearInfo{3, "UL17", "2017", 41.480, 41480., 0.023, "41.5", 800-3, 3} },
  {Year::isUL18, YearInfo{4, "UL18", "2018", 59.832, 59832., 0.025, "59.8", 820-8, 4} },

  {Year::isRun2, YearInfo{4, "run2", "Run II", 137.645, 137645., 0.016, "138", 0, 1} },
};

const std::set<Year> kAllYears = { Year::isUL16preVFP, Year::isUL16postVFP, Year::isUL17, Year::isUL18 };

//kWhite  = 0,   kBlack  = 1,   kGray    = 920,  kRed    = 632,  kGreen  = 416,
//kBlue   = 600, kYellow = 400, kMagenta = 616,  kCyan   = 432,  kOrange = 800,
//kSpring = 820, kTeal   = 840, kAzure   =  860, kViolet = 880,  kPink   = 900

enum class ProbeJetAlgo {
  isHOTVR,
  isAK8,
  notValid,
};

typedef struct {
  std::string name = "";
  double mass_min = 0.;
  double mass_max = std::numeric_limits<double>::infinity();
} ProbeJetAlgoInfo;

const std::map<ProbeJetAlgo, ProbeJetAlgoInfo> kProbeJetAlgos = {
  {ProbeJetAlgo::isHOTVR, ProbeJetAlgoInfo{"HOTVR", 140., 220.}},
  {ProbeJetAlgo::isAK8, ProbeJetAlgoInfo{"AK8", 105., 210.}},
};



class CoordinateConverter {
public:
  CoordinateConverter() {};
  void init(const float l, const float r, const float b, const float t);
  float ConvertGraphXToPadX(const float graph_x);
  float ConvertGraphYToPadY(const float graph_y);
private:
  float pad_l, pad_r, pad_b, pad_t; // pad margins
  float graph_width, graph_height;
};

void CoordinateConverter::init(const float l, const float r, const float b, const float t) {
  pad_l = l;
  pad_r = r;
  pad_b = b;
  pad_t = t;
  graph_width = 1.-l-r;
  graph_height = 1.-b-t;
}

float CoordinateConverter::ConvertGraphXToPadX(const float graph_x) {
  return pad_l+graph_x*graph_width;
}

float CoordinateConverter::ConvertGraphYToPadY(const float graph_y) {
  return pad_b+graph_y*graph_height;
}




// Formula taken from last table in https://www.phenix.bnl.gov/WWW/publish/elke/EIC/Files-for-Wiki/lara.02-008.errors.pdf
// A is the fractional set; B is the full set
pair<double, double> fraction_uncertainty(const double eA, const double eB, const double sigmaA, const double sigmaB) {
  const double fraction = eA / eB;
  const double relative_uncertainty = sqrt(sigmaA * sigmaA / (eA * eA) + sigmaB * sigmaB / (eB * eB) - 2. * sigmaA * sigmaA / (eA * eB));
  return {fraction, fraction * relative_uncertainty};
}

}
