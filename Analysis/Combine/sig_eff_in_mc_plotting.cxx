#include <iomanip>
#include <sstream>
#include <array>

#include "../constants.h"

using namespace macros;

const string data_path = "sig_eff_in_mc_BasicHists";




TGraphAsymmErrors * create_sig_eff_graph_run2(
  TFile * infile,
  const string & graphNamePostfix,
  double x_offset,
  TFile * sf_file,
  const string & tagger_and_wp
) {
  map<Year, TGraph> graph_map;
  map<Year, TGraphAsymmErrors> sf_map_corr, sf_map_uncorr;
  for (auto year : kAllYears) {
    const string graphName = kYears.at(year).name + graphNamePostfix;
    graph_map[year] = *(TGraph*)infile->Get(graphName.c_str());
    const string sfName = tagger_and_wp + "-" + kYears.at(year).name + "/FullyMerged_" + kYears.at(year).name;
    sf_map_corr[year] = *(TGraphAsymmErrors*)sf_file->Get((sfName + "_corr").c_str());
    sf_map_uncorr[year] = *(TGraphAsymmErrors*)sf_file->Get((sfName + "_uncorr").c_str());
  }

  int n = graph_map[Year::isUL18].GetN();
  //double x[n], y[n], exl[n], exh[n], eyl[n], eyh[n];
  //double exl[n], exh[n], eyl[n], eyh[n];
  vector<double> x, y, exl, exh, eyl, eyh;
  //vector<double> x, y;
  for(int i = 0; i <= n; i++) {
    double px, py, sfx, sfy;  // px is the jet pt
    //x[i] = 0;
    //y[i] = 0;
    //x.push_back(0);
    //y.push_back(0);
    //map<Year, double> relErrors_corr_l, relErrors_uncorr_l;
    //map<Year, double> relErrors_corr_h, relErrors_uncorr_h;
    double temp_x(0.), temp_y(0.);
    double e_corr_l(0.), e_corr_h(0.), e_uncorr_l(0.), e_uncorr_h(0.);
    for (auto year : kAllYears) {
      double lumi_factor = kYears.at(year).lumi_fb / kYears.at(Year::isRun2).lumi_fb;
      graph_map[year].GetPoint(i, px, py);
      //cout << px << endl;
      temp_x = px + x_offset;
      //cout << x[i] << endl;

      double sf = sf_map_corr.at(year).GetPoint(i, sfx, sfy);
      temp_y += py * lumi_factor * sfy;

      //relErrors_corr_l[year] = sf_map_corr[year].GetErrorXlow() / sfy;
      //relErrors_corr_h[year] = sf_map_corr[year].GetErrorXhigh() / sfy;
      //relErrors_uncorr_l[year] = sf_map_uncorr[year].GetErrorXlow() / sfy;
      //relErrors_uncorr_h[year] = sf_map_uncorr[year].GetErrorXhigh() / sfy;

      e_corr_l += lumi_factor * sf_map_corr.at(year).GetErrorYlow(i) / sfy;
      e_corr_h += lumi_factor * sf_map_corr.at(year).GetErrorYhigh(i) / sfy;

      e_uncorr_l += lumi_factor * pow(sf_map_uncorr.at(year).GetErrorYlow(i) / sfy, 2);
      e_uncorr_h += lumi_factor * pow(sf_map_uncorr.at(year).GetErrorYhigh(i) / sfy, 2);
    }

    x.push_back(temp_x);
    y.push_back(temp_y);

    //eyl[i] = y[i] * sqrt( e_corr_l * e_corr_l + e_uncorr_l );
    //eyh[i] = y[i] * sqrt( e_corr_h * e_corr_h + e_uncorr_h );
    //exl[i] = 0;
    //exh[i] = 0;

    eyl.push_back(y.at(i) * sqrt( e_corr_l * e_corr_l + e_uncorr_l ));
    eyh.push_back(y.at(i) * sqrt( e_corr_h * e_corr_h + e_uncorr_h ));
    exl.push_back(0);
    exh.push_back(0);
  }

  //if ( graphNamePostfix.find("BkgEff") != string::npos ) x[0] = 350 + x_offset; // FIX for some unexplainable bug where x[0] is always some small number albeit px was found correctly

  double ax[n], ay[n];
  double aexl[n], aexh[n], aeyl[n], aeyh[n];
  for (int i = 0; i < n ; i++) {
    ax[i] = x.at(i);
    ay[i] = y.at(i);
    aexl[i] = exl.at(i);
    aexh[i] = exh.at(i);
    aeyl[i] = eyl.at(i);
    aeyh[i] = eyh.at(i);
  }

  TGraphAsymmErrors * result = new TGraphAsymmErrors(n, ax, ay, aexl, aexh, aeyl, aeyh);
  return result;

  // const string graphName = kYears.at(Year::isUL16preVFP).name + graphNamePostfix;
  // TGraph graph_UL16preVFP = *(TGraph*)infile->Get(graphName.c_str());
}




TGraph * new_null_graph() {
  int n = 1001; // should be the same number + 1 as the variable "nx" in analyze.py
  double x[n], y[n];
  for(int i = 0; i <= n; i++) {
    x[i] = i/double(n);
    y[i] = i/double(n);
  }
  TGraph * graph = new TGraph(n, x, y);
  graph->SetLineColorAlpha(kBlack, 0); // 0 = fully transparent
  return graph;
}



// {
//    double x[100], y[100];
//    int n = 20;
//    for (int i=0;i<n;i++) {
//      x[i] = i*0.1;
//      y[i] = 10*sin(x[i]+0.2);
//    }
//    auto g = new TGraph(n,x,y);
//    g->SetTitle("Graph title;X title;Y title");
//    g->Draw("AC*");
// }

//   double y[6] = {3, 8, 1, 10, 5, 7};


//________________________________________________________________________________
// Leonidas' values from his CSV Excel sheets:
// pt bin edges W: 200/300/400/800
// pt bin edges top: 300/400/480/600/1200
static const double leo_ptCenter_W[3] = {250, 350, 600};
static const double leo_ptCenter_top[4] = {350, 440, 540, 700};
// WvsQCD
typedef map<Year, std::array<double, 3>> values_leo_type_3;
typedef map<Year, std::array<double, 4>> values_leo_type_4;
static const values_leo_type_3 values_leo_WvsQCD_BkgEff0p050 = {
  { Year::isUL16preVFP, { 0.8425, 0.7367, 0.7009 } },
  { Year::isUL16postVFP, { 0.8426, 0.7348, 0.7009 } },
  { Year::isUL17, { 0.8377, 0.7297, 0.6967 } },
  { Year::isUL18, { 0.8391, 0.7350, 0.6997 } },
};
static const values_leo_type_3 values_leo_WvsQCD_BkgEff0p010 = {
  { Year::isUL16preVFP, { 0.6880, 0.5918, 0.5630 } },
  { Year::isUL16postVFP, { 0.6851, 0.5875, 0.5591 } },
  { Year::isUL17, { 0.6784, 0.5829, 0.5554 } },
  { Year::isUL18, { 0.6854, 0.5920, 0.5646 } },
};
static const values_leo_type_3 values_leo_WvsQCD_BkgEff0p005 = {
  { Year::isUL16preVFP, { 0.5493, 0.4649, 0.4436 } },
  { Year::isUL16postVFP, { 0.5410, 0.4609, 0.4410 } },
  { Year::isUL17, { 0.5369, 0.4586, 0.4418 } },
  { Year::isUL18, { 0.5143, 0.4394, 0.4259 } },
};
// WvsQCDMD
static const values_leo_type_3 values_leo_WvsQCDMD_BkgEff0p025 = {
  { Year::isUL16preVFP, { 0.7231, 0.6797, 0.6459 } },
  { Year::isUL16postVFP, { 0.7124, 0.6696, 0.6419 } },
  { Year::isUL17, { 0.7249, 0.6723, 0.6395 } },
  { Year::isUL18, { 0.7269, 0.6812, 0.6485 } },
};
static const values_leo_type_3 values_leo_WvsQCDMD_BkgEff0p010 = {
  { Year::isUL16preVFP, { 0.5139, 0.5127, 0.4929 } },
  { Year::isUL16postVFP, { 0.5079, 0.5107, 0.4900 } },
  { Year::isUL17, { 0.5218, 0.5122, 0.4889 } },
  { Year::isUL18, { 0.5191, 0.5172, 0.4980 } },
};
static const values_leo_type_3 values_leo_WvsQCDMD_BkgEff0p005 = {
  { Year::isUL16preVFP, { 0.3701, 0.3895, 0.3855 } },
  { Year::isUL16postVFP, { 0.3677, 0.3903, 0.3844 } },
  { Year::isUL17, { 0.3744, 0.3865, 0.3797 } },
  { Year::isUL18, { 0.3647, 0.3845, 0.3819 } },
};
// TvsQCD
static const values_leo_type_4 values_leo_TvsQCD_BkgEff0p010 = {
  { Year::isUL16preVFP, { 0.8942, 0.9219, 0.9050, 0.8546 } },
  { Year::isUL16postVFP, { 0.8875, 0.9217, 0.9045, 0.8572 } },
  { Year::isUL17, { 0.8747, 0.9209, 0.9025, 0.8672 } },
  { Year::isUL18, { 0.8802, 0.9182, 0.9022, 0.8617 } },
};
static const values_leo_type_4 values_leo_TvsQCD_BkgEff0p005 = {
  { Year::isUL16preVFP, { 0.8575, 0.8952, 0.8812, 0.8326 } },
  { Year::isUL16postVFP, { 0.8502, 0.8946, 0.8770, 0.8358 } },
  { Year::isUL17, { 0.8329, 0.8936, 0.8779, 0.8449 } },
  { Year::isUL18, { 0.8401, 0.8907, 0.8770, 0.8402 } },
};
static const values_leo_type_4 values_leo_TvsQCD_BkgEff0p001 = {
  { Year::isUL16preVFP, { 0.7279, 0.7912, 0.7869, 0.7436 } },
  { Year::isUL16postVFP, { 0.7189, 0.7913, 0.7824, 0.7462 } },
  { Year::isUL17, { 0.7018, 0.7899, 0.7833, 0.7581 } },
  { Year::isUL18, { 0.7106, 0.7882, 0.7837, 0.7559 } },
};

enum class E_leo_Tagger {
  isTvsQCD,
  isWvsQCD,
  isWvsQCDMD,
};
enum class E_leo_WP {
  isBkgEff0p001,
  isBkgEff0p005,
  isBkgEff0p010,
  isBkgEff0p025,
  isBkgEff0p050,
};

static const map<E_leo_WP, values_leo_type_4> values_leo_partnet_top = {
  { E_leo_WP::isBkgEff0p010, values_leo_TvsQCD_BkgEff0p010 },
  { E_leo_WP::isBkgEff0p005, values_leo_TvsQCD_BkgEff0p005 },
  { E_leo_WP::isBkgEff0p001, values_leo_TvsQCD_BkgEff0p001 },
};

static const map<E_leo_Tagger, map<E_leo_WP, values_leo_type_3>> values_leo_partnet_w = {
  {
    E_leo_Tagger::isWvsQCD, {
      { E_leo_WP::isBkgEff0p050, values_leo_WvsQCD_BkgEff0p050 },
      { E_leo_WP::isBkgEff0p010, values_leo_WvsQCD_BkgEff0p010 },
      { E_leo_WP::isBkgEff0p005, values_leo_WvsQCD_BkgEff0p005 },
    }
  },
  {
    E_leo_Tagger::isWvsQCDMD, {
      { E_leo_WP::isBkgEff0p025, values_leo_WvsQCDMD_BkgEff0p025 },
      { E_leo_WP::isBkgEff0p010, values_leo_WvsQCDMD_BkgEff0p010 },
      { E_leo_WP::isBkgEff0p005, values_leo_WvsQCDMD_BkgEff0p005 },
    }
  },
};

typedef struct {
  double pt_low;
  double pt_high;
  double sf;
  double sf_down;
  double sf_up;
} leoSF;

static const map<E_leo_Tagger, map<Year, map<E_leo_WP, vector<leoSF>>>> leoScaleFactors = {
  {
    E_leo_Tagger::isWvsQCD, {
      {
        Year::isUL18, {
          {
            E_leo_WP::isBkgEff0p050, {
              {200,300,1.015232832331004,0.07299516180427412,0.0329949782319105},
              {300,400,1.0262724517712463,0.0796128335692613,0.09061680739568656},
              {400,800,1.1618406081440518,0.1427165764604552,0.1623263594732257},
            }
          },
          {
            E_leo_WP::isBkgEff0p010, {
              {200,300,1.0207719507098578,0.02690925331582794,0.022482432399033303},
              {300,400,0.9662561264360947,0.05594758980276515,0.05535242117725969},
              {400,800,1.0342599879066012,0.1560764591233591,0.059018250935544414},
            }
          },
          {
            E_leo_WP::isBkgEff0p005, {
              {200,300,0.8916898767041578,0.025075497400440905,0.022158744403629638},
              {300,400,0.8739086463525401,0.054830386704272205,0.025863559827948746},
              {400,800,0.874348326268582,0.08838743871532617,0.06205314407461171},
            }
          },
        }
      },
      {
        Year::isUL17, {
          {
            E_leo_WP::isBkgEff0p050, {
              {200,300,1.0594467407510653,0.06450862631724785,0.05141082528499452},
              {300,400,1.0681023994893237,0.07090597397908116,0.10857980986795207},
              {400,800,1.1088731271809773,0.1712376434049938,0.16189605518799743},
            }
          },
          {
            E_leo_WP::isBkgEff0p010, {
              {200,300,0.9587246533669983,0.040422897033849936,0.04680109244747821},
              {300,400,0.9439519374013672,0.052956710009177455,0.08460017084174382},
              {400,800,0.7885701617038146,0.09584342583997385,0.1244672628782863},
            }
          },
          {
            E_leo_WP::isBkgEff0p005, {
              {200,300,0.9312602988581649,0.03757637686976645,0.025369301998965865},
              {300,400,0.8888791161455953,0.03204540181548643,0.02916984299483044},
              {400,800,0.7874248132193202,0.11770157229453493,0.11614097356312864},
            }
          },
        }
      },
      {
        Year::isUL16postVFP, {
          {
            E_leo_WP::isBkgEff0p050, {
              {200,300,1.3252393934339752,0.2637682305439748,0.3690533462901563},
              {300,400,1.03602530145928,0.13461939426429326,0.2335249972758235},
              {400,800,1.057893614667379,0.14366556189273427,0.13293171918031427},
            }
          },
          {
            E_leo_WP::isBkgEff0p010, {
              {200,300,1.1277146862875185,0.06347493639596635,0.06355513497258214},
              {300,400,1.0836717298969585,0.049694041648482346,0.0456811845094115},
              {400,800,1.155322907991434,0.13837723466132212,0.07985775972577314},
            }
          },
          {
            E_leo_WP::isBkgEff0p005, {
              {200,300,1.0070754979118601,0.03924686225998919,0.03744092024580048},
              {300,400,0.9544834920137664,0.05977861242522742,0.04202920805173388},
              {400,800,0.9520494290102449,0.09900686900643507,0.09659807429823458},
            }
          },
        }
      },
      {
        Year::isUL16preVFP, {
          {
            E_leo_WP::isBkgEff0p050, {
              {200,300,1.068682130166015,0.03152041939219363,0.03},
              {300,400,0.9633918974993291,0.07491282006829247,0.08394047863056131},
              {400,800,1.1669175779068541,0.069,0.07},
            }
          },
          {
            E_leo_WP::isBkgEff0p010, {
              {200,300,1.0112160124733622,0.07,0.07},
              {300,400,1.0030071418307132,0.0938784945034159,0.10313623868572142},
              {400,800,1.0529715854028163,0.10026392476050128,0.1004798066677628},
            }
          },
          {
            E_leo_WP::isBkgEff0p005, {
              {200,300,0.9885688949491379,0.06403079626501329,0.06146208378962281},
              {300,400,1.029636050453152,0.10405381164675886,0.044860843269695805},
              {400,800,1.0374200373620732,0.09092040667580814,0.07030884889378247},
            }
          },
        }
      },
    }
  },
  {
    E_leo_Tagger::isWvsQCDMD, {
      {
        Year::isUL18, {
          {
            E_leo_WP::isBkgEff0p025, {
              {200,300,0.9190298705570141,0.019133418313060258,0.022213209889155705},
              {300,400,0.91497012004329,0.021042142455296897,0.02020293492707881},
              {400,800,0.8647360917703466,0.03861458078723401,0.04112924159099074},
            }
          },
          {
            E_leo_WP::isBkgEff0p010, {
              {200,300,0.8634666100405313,0.01976610829273484,0.020531964646055834},
              {300,400,0.8614614438439079,0.02132661802511293,0.021720285208615353},
              {400,800,0.8152646656428022,0.03804387486150085,0.038157492538026916},
            }
          },
          {
            E_leo_WP::isBkgEff0p005, {
              {200,300,0.8142381072562699,0.02080028476344986,0.025445994517837023},
              {300,400,0.8144892282975787,0.02307827769971449,0.023090889513834978},
              {400,800,0.7686201753237074,0.03886600828228404,0.04079137354527396},
            }
          },
        }
      },
      {
        Year::isUL17, {
          {
            E_leo_WP::isBkgEff0p025, {
              {200,300,0.982730281923659,0.036673732392177194,0.03593481145009797},
              {300,400,0.9594886915963721,0.02437641178667438,0.02474377023123797},
              {400,800,0.9892794381973279,0.04497765154075517,0.04720321554228557},
            }
          },
          {
            E_leo_WP::isBkgEff0p010, {
              {200,300,0.918743136510155,0.04068105977223624,0.04253825676422884},
              {300,400,0.898412171557027,0.02500432426821242,0.024991547234687606},
              {400,800,0.9039402189100576,0.045117445347443996,0.04548813841444799},
            }
          },
          {
            E_leo_WP::isBkgEff0p005, {
              {200,300,0.8555180974523378,0.0283034531162365,0.032316100577830575},
              {300,400,0.8440061389098098,0.025807274531479574,0.02635991726279796},
              {400,800,0.8608276370198319,0.04413138554471052,0.04597389622181508},
            }
          },
        }
      },
      {
        Year::isUL16postVFP, {
          {
            E_leo_WP::isBkgEff0p025, {
              {200,300,0.9450578338842212,0.0636582927237117,0.05026579254779012},
              {300,400,0.9092346655321338,0.03703680267053133,0.03772561131817631},
              {400,800,0.8925475131701243,0.06871418729146117,0.07248030367171265},
            }
          },
          {
            E_leo_WP::isBkgEff0p010, {
              {200,300,0.887064825215807,0.04387845567911064,0.04742499030716374},
              {300,400,0.8662333616667055,0.03701362814765874,0.038059141188278967},
              {400,800,0.7866962271032931,0.06577430373422988,0.06748947013573098},
            }
          },
          {
            E_leo_WP::isBkgEff0p005, {
              {200,300,0.8578770720810618,0.041861699988094636,0.023774482530758756},
              {300,400,0.8313985780254843,0.0389008458789295,0.03997533404797682},
              {400,800,0.7322054452349709,0.06039902846704037,0.06300043645628861},
            }
          },
        }
      },
      {
        Year::isUL16preVFP, {
          {
            E_leo_WP::isBkgEff0p025, {
              {200,300,0.9195379553371,0.06554066913349732,0.036269103542342196},
              {300,400,0.9463957277711702,0.03469390108094872,0.03527766629510892},
              {400,800,0.9398161345034564,0.0598211031246495,0.06321678642886641},
            }
          },
          {
            E_leo_WP::isBkgEff0p010, {
              {200,300,0.9042614655651481,0.0316215048707319,0.03218926479843087},
              {300,400,0.8757183251530393,0.03414284963301184,0.03527721784960958},
              {400,800,0.9022723963132271,0.06549020119568982,0.06770280642640181},
            }
          },
          {
            E_leo_WP::isBkgEff0p005, {
              {200,300,0.8451073945199086,0.03740932701944466,0.03608965667673614},
              {300,400,0.8659212176684473,0.03819408101916044,0.038788281650652956},
              {400,800,0.8702821494336126,0.06541945721213815,0.06759452062068222},
            }
          },
        }
      },
    }
  },
  {
    E_leo_Tagger::isTvsQCD, {
      {
        Year::isUL18, {
          {
            E_leo_WP::isBkgEff0p010, {
              {300,400,1.1163204630587051,0.08589322774156138,0.1145103421879804},
              {400,480,0.9878718467695988,0.03394913833483326,0.03536992757274554},
              {480,600,0.988999891725833,0.03135142808228031,0.03254552277081374},
              {600,1200,0.9805739808831712,0.05931408212949585,0.051187894522590094},
            }
          },
          {
            E_leo_WP::isBkgEff0p005, {
              {300,400,1.1920909462006473,0.11764549348223086,0.11908708002963475},
              {400,480,0.9774648125322577,0.035550840455008226,0.0398907134046641},
              {480,600,0.9610178961036971,0.031596130357776864,0.03270118031346542},
              {600,1200,0.9719547878688892,0.04718604414578065,0.04838027950930873},
            }
          },
          {
            E_leo_WP::isBkgEff0p001, {
              {300,400,1.0334005773459902,0.07650971737263357,0.09380051325239541},
              {400,480,0.949257944706414,0.045225626829771715,0.046893611197207585},
              {480,600,0.9107008737018066,0.032422008884733344,0.03424149006948135},
              {600,1200,0.9527002945910971,0.05328536808253126,0.06053691794789251},
            }
          },
        }
      },
      {
        Year::isUL17, {
          {
            E_leo_WP::isBkgEff0p010, {
              {300,400,1.2665423171877026,0.12594170644336278,0.13075885577178425},
              {400,480,1.01651719047938,0.04200062224205936,0.03718799071755147},
              {480,600,1.0563219429599855,0.04406769273143074,0.12450427752619087},
              {600,1200,1.004343034237564,0.05644574156991433,0.06164302682382439},
            }
          },
          {
            E_leo_WP::isBkgEff0p005, {
              {300,400,1.1062034119256086,0.07529397680618755,0.08045105877312297},
              {400,480,1.0066113887295043,0.0369232520042021,0.03919635614723771},
              {480,600,1.049061790460786,0.04147167493389903,0.04338946817230571},
              {600,1200,0.9973377822040789,0.050104016301182774,0.05405382873318071},
            }
          },
          {
            E_leo_WP::isBkgEff0p001, {
              {300,400,1.1196092503548305,0.09802412387572934,0.11200766928244099},
              {400,480,0.9580004039531106,0.03816914560278717,0.04072105441844143},
              {480,600,1.0045270269743276,0.04642644457267664,0.04962975757396476},
              {600,1200,0.9323482894176377,0.05638889259445923,0.051473939960285375},
            }
          },
        }
      },
      {
        Year::isUL16postVFP, {
          {
            E_leo_WP::isBkgEff0p010, {
              {300,400,1.0286204672478074,0.07738336621550035,0.0837428626798436},
              {400,480,0.9861968638399904,0.046,0.05284760599583793},
              {480,600,1.1295417328901183,0.14011597977479784,0.17885647825191325},
              {600,1200,1.3035872134528272,0.2775637935767241,0.31002043724845657},
            }
          },
          {
            E_leo_WP::isBkgEff0p005, {
              {300,400,1.0712531531907075,0.09389896914842288,0.10059110760657686},
              {400,480,0.9672581863799624,0.04762944565463889,0.05182789160227019},
              {480,600,1.0340496715510143,0.0406073589552004,0.05083348138255023},
              {600,1200,1.3396199462217941,0.2711476744604726,0.23822635606740838},
            }
          },
          {
            E_leo_WP::isBkgEff0p001, {
              {300,400,0.9693713610029029,0.08151633318067164,0.08689368400132591},
              {400,480,0.896687995785267,0.04715255343674973,0.04840402903850832},
              {480,600,0.9968234331208976,0.05852845489084546,0.05677047796763779},
              {600,1200,1.0103636762075094,0.08543865285456187,0.09037365986555956},
            }
          },
        }
      },
      {
        Year::isUL16preVFP, {
          {
            E_leo_WP::isBkgEff0p010, {
              {300,400,1.2529632568286877,0.1401589562980461,0.14893548637024145},
              {400,480,1.140886941405853,0.11013160312991532,0.1453378248525864},
              {480,600,1.2321784406468166,0.16845345175872173,0.17834081774322008},
              {600,1200,1.3823424969039724,0.3203322566182154,0.35739488023541155},
            }
          },
          {
            E_leo_WP::isBkgEff0p005, {
              {300,400,1.2309559151349918,0.14489036980725922,0.13490834643442517},
              {400,480,1.0738515100011004,0.04694235506050637,0.1291071030184015},
              {480,600,1.057574178518796,0.05238859810504026,0.05645882027015936},
              {600,1200,1.0994421241009649,0.0981719413385509,0.28653083056383055},
            }
          },
          {
            E_leo_WP::isBkgEff0p001, {
              {300,400,1.0952484275692702,0.07172482413510917,0.07631575391928647},
              {400,480,1.0567220114984588,0.05890151615878081,0.06400976396644487},
              {480,600,1.0423524947175915,0.05886351490829689,0.06347414682865082},
              {600,1200,0.9978899695223552,0.0906911555419252,0.10535484120276511},
            }
          },
        }
      },
    }
  },
};


TGraphAsymmErrors create_sig_eff_graph_run2_partnet_w(
  E_leo_Tagger e_leo_Tagger,
  E_leo_WP e_leo_WP,
  double x_offset
) {
  const int n = 3;

  vector<double> x, y, exl, exh, eyl, eyh;

  for (int i = 0; i < n; i++) {
    double temp_x = leo_ptCenter_W[i];
    double temp_y = 0.;
    double temp_eyl = 0.;
    double temp_eyh = 0.;
    for (const auto year : kAllYears) {
      double lumi_factor = kYears.at(year).lumi_fb / kYears.at(Year::isRun2).lumi_fb;

      double sf = 1.; // this year's SF
      double sf_eyl = 0;
      double sf_eyh = 0;
      const auto scaleFactors = leoScaleFactors.at(e_leo_Tagger).at(year).at(e_leo_WP);
      for (const auto & scaleFactor : scaleFactors) {
        if (scaleFactor.pt_low > temp_x || scaleFactor.pt_high < temp_x) continue;
        sf = scaleFactor.sf;
        sf_eyl = scaleFactor.sf_down;
        sf_eyh = scaleFactor.sf_up;
        break;
      }

      temp_y += sf * lumi_factor * values_leo_partnet_w.at(e_leo_Tagger).at(e_leo_WP).at(year)[i];

      temp_eyl += lumi_factor * pow(sf_eyl / sf, 2);
      temp_eyh += lumi_factor * pow(sf_eyh / sf, 2);
    }
    x.push_back(temp_x + x_offset);
    y.push_back(temp_y);
    exl.push_back(0);
    exh.push_back(0);
    temp_eyl = temp_y * sqrt(temp_eyl);
    temp_eyh = temp_y * sqrt(temp_eyh);
    eyl.push_back(temp_eyl);
    eyh.push_back(temp_eyh);
  }

  double ax[n], ay[n], aexl[n], aexh[n], aeyl[n], aeyh[n];

  for (int i = 0; i < n; i++) {
    ax[i] = x.at(i);
    ay[i] = y.at(i);
    aexl[i] = exl.at(i);
    aexh[i] = exh.at(i);
    aeyl[i] = eyl.at(i);
    aeyh[i] = eyh.at(i);
  }

  TGraphAsymmErrors result = TGraphAsymmErrors(n, ax, ay, aexl, aexh, aeyl, aeyh);
  return result;
}
  /* map<Year, TGraph> graph_map;
  map<Year, TGraphAsymmErrors> sf_map_corr, sf_map_uncorr;
  for (auto year : kAllYears) {
    const string graphName = kYears.at(year).name + graphNamePostfix;
    graph_map[year] = *(TGraph*)infile->Get(graphName.c_str());
    const string sfName = tagger_and_wp + "-" + kYears.at(year).name + "/FullyMerged_" + kYears.at(year).name;
    sf_map_corr[year] = *(TGraphAsymmErrors*)sf_file->Get((sfName + "_corr").c_str());
    sf_map_uncorr[year] = *(TGraphAsymmErrors*)sf_file->Get((sfName + "_uncorr").c_str());
  }

  int n = graph_map[Year::isUL18].GetN();
  //double x[n], y[n], exl[n], exh[n], eyl[n], eyh[n];
  //double exl[n], exh[n], eyl[n], eyh[n];
  vector<double> x, y, exl, exh, eyl, eyh;
  //vector<double> x, y;
  for(int i = 0; i <= n; i++) {
    double px, py, sfx, sfy;  // px is the jet pt
    //x[i] = 0;
    //y[i] = 0;
    //x.push_back(0);
    //y.push_back(0);
    //map<Year, double> relErrors_corr_l, relErrors_uncorr_l;
    //map<Year, double> relErrors_corr_h, relErrors_uncorr_h;
    double temp_x(0.), temp_y(0.);
    double e_corr_l(0.), e_corr_h(0.), e_uncorr_l(0.), e_uncorr_h(0.);
    for (auto year : kAllYears) {
      double lumi_factor = kYears.at(year).lumi_fb / kYears.at(Year::isRun2).lumi_fb;
      graph_map[year].GetPoint(i, px, py);
      //cout << px << endl;
      temp_x = px + x_offset;
      //cout << x[i] << endl;

      double sf = sf_map_corr.at(year).GetPoint(i, sfx, sfy);
      temp_y += py * lumi_factor * sfy;

      //relErrors_corr_l[year] = sf_map_corr[year].GetErrorXlow() / sfy;
      //relErrors_corr_h[year] = sf_map_corr[year].GetErrorXhigh() / sfy;
      //relErrors_uncorr_l[year] = sf_map_uncorr[year].GetErrorXlow() / sfy;
      //relErrors_uncorr_h[year] = sf_map_uncorr[year].GetErrorXhigh() / sfy;

      e_corr_l += lumi_factor * sf_map_corr.at(year).GetErrorYlow(i) / sfy;
      e_corr_h += lumi_factor * sf_map_corr.at(year).GetErrorYhigh(i) / sfy;

      e_uncorr_l += lumi_factor * pow(sf_map_uncorr.at(year).GetErrorYlow(i) / sfy, 2);
      e_uncorr_h += lumi_factor * pow(sf_map_uncorr.at(year).GetErrorYhigh(i) / sfy, 2);
    }

    x.push_back(temp_x);
    y.push_back(temp_y);

    //eyl[i] = y[i] * sqrt( e_corr_l * e_corr_l + e_uncorr_l );
    //eyh[i] = y[i] * sqrt( e_corr_h * e_corr_h + e_uncorr_h );
    //exl[i] = 0;
    //exh[i] = 0;

    eyl.push_back(y.at(i) * sqrt( e_corr_l * e_corr_l + e_uncorr_l ));
    eyh.push_back(y.at(i) * sqrt( e_corr_h * e_corr_h + e_uncorr_h ));
    exl.push_back(0);
    exh.push_back(0);
  }

  //if ( graphNamePostfix.find("BkgEff") != string::npos ) x[0] = 350 + x_offset; // FIX for some unexplainable bug where x[0] is always some small number albeit px was found correctly

  double ax[n], ay[n];
  double aexl[n], aexh[n], aeyl[n], aeyh[n];
  for (int i = 0; i < n ; i++) {
    ax[i] = x.at(i);
    ay[i] = y.at(i);
    aexl[i] = exl.at(i);
    aexh[i] = exh.at(i);
    aeyl[i] = eyl.at(i);
    aeyh[i] = eyh.at(i);
  }

  TGraphAsymmErrors * result = new TGraphAsymmErrors(n, ax, ay, aexl, aexh, aeyl, aeyh);
  return result;

  // const string graphName = kYears.at(Year::isUL16preVFP).name + graphNamePostfix;
  // TGraph graph_UL16preVFP = *(TGraph*)infile->Get(graphName.c_str());
} */



//void do_plot_nsubjettiness(const Year & year, const bool bool_prelim) {
void do_plot_nsubjettiness(const bool bool_prelim) {

  int debug = 0;

  const string inputGraphName = "graph_sig_eff_withWindow_both";

  const string infilePath_nominal = data_path + "/sig_eff_in_mc-ak8_t__tau.root";
  const string infilePath_btag = data_path + "/sig_eff_in_mc-ak8_t_btagDJet__tau.root";
  const string infilePath_hotvr = data_path + "/sig_eff_in_mc-hotvr_t__tau.root";

  TFile* infile_nominal = TFile::Open(infilePath_nominal.c_str(), "READ");
  TFile* infile_btag = TFile::Open(infilePath_btag.c_str(), "READ");
  TFile* infile_hotvr = TFile::Open(infilePath_hotvr.c_str(), "READ");

  const int line_width = 2;
  const int marker_size = 1;
  const int line_style_nominal = 1;
  const int line_style_btag = 2;

  TMultiGraph *mg = new TMultiGraph("mg", "mg title");

  cout << debug++ << endl;

  double x_offset_nominal = 4;
  double x_offset_btag = -x_offset_nominal;

  double x_offset_general = 16;

  string tagger_and_wp = "";

  //const string graphName_nominal_vt = kYears.at(year).name + "-BkgEff0p001/"+inputGraphName;
  //TGraph graph_nominal_vt = *(TGraph*)infile_nominal->Get(graphName_nominal_vt.c_str());
  const string graphNamePostfix_nominal_vt = "-BkgEff0p001/"+inputGraphName;
  TFile* sf_file_nominal_vt = TFile::Open("scaleFactors_2023-08-10/scaleFactors-ak8_t__tau-BkgEff0p001/scaleFactors-ak8_t__tau-BkgEff0p001.root", "READ");
  tagger_and_wp = "ak8_t__tau-BkgEff0p001";
  auto graph_nominal_vt = *(create_sig_eff_graph_run2(infile_nominal, graphNamePostfix_nominal_vt, x_offset_nominal-2*x_offset_general, sf_file_nominal_vt, tagger_and_wp));
  graph_nominal_vt.SetLineWidth(line_width);
  graph_nominal_vt.SetLineColor(kGreen);
  graph_nominal_vt.SetLineStyle(line_style_nominal);
  graph_nominal_vt.SetMarkerColor(kGreen);
  graph_nominal_vt.SetMarkerSize(marker_size);
  graph_nominal_vt.SetMarkerStyle(20);
  mg->Add(&graph_nominal_vt);

  cout << debug++ << endl;
  
  //const string graphName_nominal_ti = kYears.at(year).name + "-BkgEff0p005/"+inputGraphName;
  //TGraph graph_nominal_ti = *(TGraph*)infile_nominal->Get(graphName_nominal_ti.c_str());
  const string graphNamePostfix_nominal_ti = "-BkgEff0p005/"+inputGraphName;
  TFile* sf_file_nominal_ti = TFile::Open("scaleFactors_2023-08-10/scaleFactors-ak8_t__tau-BkgEff0p005/scaleFactors-ak8_t__tau-BkgEff0p005.root", "READ");
  tagger_and_wp = "ak8_t__tau-BkgEff0p005";
  auto graph_nominal_ti = *(create_sig_eff_graph_run2(infile_nominal, graphNamePostfix_nominal_ti, x_offset_nominal-x_offset_general, sf_file_nominal_ti, tagger_and_wp));
  graph_nominal_ti.SetLineWidth(line_width);
  graph_nominal_ti.SetLineColor(kCyan);
  graph_nominal_ti.SetLineStyle(line_style_nominal);
  graph_nominal_ti.SetMarkerColor(kCyan);
  graph_nominal_ti.SetMarkerSize(marker_size);
  graph_nominal_ti.SetMarkerStyle(21);
  mg->Add(&graph_nominal_ti);

  cout << debug++ << endl;

  //const string graphName_nominal_me = kYears.at(year).name + "-BkgEff0p010/"+inputGraphName;
  //TGraph graph_nominal_me = *(TGraph*)infile_nominal->Get(graphName_nominal_me.c_str());
  const string graphNamePostfix_nominal_me = "-BkgEff0p010/"+inputGraphName;
  TFile* sf_file_nominal_me = TFile::Open("scaleFactors_2023-08-10/scaleFactors-ak8_t__tau-BkgEff0p010/scaleFactors-ak8_t__tau-BkgEff0p010.root", "READ");
  tagger_and_wp = "ak8_t__tau-BkgEff0p010";
  auto graph_nominal_me = *(create_sig_eff_graph_run2(infile_nominal, graphNamePostfix_nominal_me, x_offset_nominal, sf_file_nominal_me, tagger_and_wp));
  graph_nominal_me.SetLineWidth(line_width);
  graph_nominal_me.SetLineColor(kBlue);
  graph_nominal_me.SetLineStyle(line_style_nominal);
  graph_nominal_me.SetMarkerColor(kBlue);
  graph_nominal_me.SetMarkerSize(marker_size);
  graph_nominal_me.SetMarkerStyle(22);
  mg->Add(&graph_nominal_me);

  cout << debug++ << endl;

  //const string graphName_nominal_lo = kYears.at(year).name + "-BkgEff0p025/"+inputGraphName;
  //TGraph graph_nominal_lo = *(TGraph*)infile_nominal->Get(graphName_nominal_lo.c_str());
  const string graphNamePostfix_nominal_lo = "-BkgEff0p025/"+inputGraphName;
  TFile* sf_file_nominal_lo = TFile::Open("scaleFactors_2023-08-10/scaleFactors-ak8_t__tau-BkgEff0p025/scaleFactors-ak8_t__tau-BkgEff0p025.root", "READ");
  tagger_and_wp = "ak8_t__tau-BkgEff0p025";
  auto graph_nominal_lo = *(create_sig_eff_graph_run2(infile_nominal, graphNamePostfix_nominal_lo, x_offset_nominal+x_offset_general, sf_file_nominal_lo, tagger_and_wp));
  graph_nominal_lo.SetLineWidth(line_width);
  graph_nominal_lo.SetLineColor(kMagenta);
  graph_nominal_lo.SetLineStyle(line_style_nominal);
  graph_nominal_lo.SetMarkerColor(kMagenta);
  graph_nominal_lo.SetMarkerSize(marker_size);
  graph_nominal_lo.SetMarkerStyle(23);
  mg->Add(&graph_nominal_lo);

  cout << debug++ << endl;

  //const string graphName_nominal_vl = kYears.at(year).name + "-BkgEff0p050/"+inputGraphName;
  //TGraph graph_nominal_vl = *(TGraph*)infile_nominal->Get(graphName_nominal_vl.c_str());
  const string graphNamePostfix_nominal_vl = "-BkgEff0p050/"+inputGraphName;
  TFile* sf_file_nominal_vl = TFile::Open("scaleFactors_2023-08-10/scaleFactors-ak8_t__tau-BkgEff0p050/scaleFactors-ak8_t__tau-BkgEff0p050.root", "READ");
  tagger_and_wp = "ak8_t__tau-BkgEff0p050";
  auto graph_nominal_vl = *(create_sig_eff_graph_run2(infile_nominal, graphNamePostfix_nominal_vl, x_offset_nominal+2*x_offset_general, sf_file_nominal_vl, tagger_and_wp));
  graph_nominal_vl.SetLineWidth(line_width);
  graph_nominal_vl.SetLineColor(kRed);
  graph_nominal_vl.SetLineStyle(line_style_nominal);
  graph_nominal_vl.SetMarkerColor(kRed);
  graph_nominal_vl.SetMarkerSize(marker_size);
  graph_nominal_vl.SetMarkerStyle(47);
  mg->Add(&graph_nominal_vl);

  cout << debug++ << endl;




  //const string graphName_btag_vt = kYears.at(year).name + "-BkgEff0p001/"+inputGraphName;
  //TGraph graph_btag_vt = *(TGraph*)infile_btag->Get(graphName_btag_vt.c_str());
  const string graphNamePostfix_btag_vt = "-BkgEff0p001/"+inputGraphName;
  TFile* sf_file_btag_vt = TFile::Open("scaleFactors_2023-08-10/scaleFactors-ak8_t_btagDJet__tau-BkgEff0p001/scaleFactors-ak8_t_btagDJet__tau-BkgEff0p001.root", "READ");
  tagger_and_wp = "ak8_t_btagDJet__tau-BkgEff0p001";
  auto graph_btag_vt = *(create_sig_eff_graph_run2(infile_btag, graphNamePostfix_btag_vt, x_offset_btag-2*x_offset_general, sf_file_btag_vt, tagger_and_wp));
  graph_btag_vt.SetLineWidth(line_width);
  graph_btag_vt.SetLineColor(kGreen);
  graph_btag_vt.SetLineStyle(line_style_btag);
  graph_btag_vt.SetMarkerColor(kGreen);
  graph_btag_vt.SetMarkerSize(marker_size);
  graph_btag_vt.SetMarkerStyle(24);
  mg->Add(&graph_btag_vt);

  cout << debug++ << endl;
  
  //const string graphName_btag_ti = kYears.at(year).name + "-BkgEff0p005/"+inputGraphName;
  //TGraph graph_btag_ti = *(TGraph*)infile_btag->Get(graphName_btag_ti.c_str());
  const string graphNamePostfix_btag_ti = "-BkgEff0p005/"+inputGraphName;
  TFile* sf_file_btag_ti = TFile::Open("scaleFactors_2023-08-10/scaleFactors-ak8_t_btagDJet__tau-BkgEff0p005/scaleFactors-ak8_t_btagDJet__tau-BkgEff0p005.root", "READ");
  tagger_and_wp = "ak8_t_btagDJet__tau-BkgEff0p005";
  auto graph_btag_ti = *(create_sig_eff_graph_run2(infile_btag, graphNamePostfix_btag_ti, x_offset_btag-x_offset_general, sf_file_btag_ti, tagger_and_wp));
  graph_btag_ti.SetLineWidth(line_width);
  graph_btag_ti.SetLineColor(kCyan);
  graph_btag_ti.SetLineStyle(line_style_btag);
  graph_btag_ti.SetMarkerColor(kCyan);
  graph_btag_ti.SetMarkerSize(marker_size);
  graph_btag_ti.SetMarkerStyle(25);
  mg->Add(&graph_btag_ti);

  cout << debug++ << endl;

  //const string graphName_btag_me = kYears.at(year).name + "-BkgEff0p010/"+inputGraphName;
  //TGraph graph_btag_me = *(TGraph*)infile_btag->Get(graphName_btag_me.c_str());
  const string graphNamePostfix_btag_me = "-BkgEff0p010/"+inputGraphName;
  TFile* sf_file_btag_me = TFile::Open("scaleFactors_2023-08-10/scaleFactors-ak8_t_btagDJet__tau-BkgEff0p010/scaleFactors-ak8_t_btagDJet__tau-BkgEff0p010.root", "READ");
  tagger_and_wp = "ak8_t_btagDJet__tau-BkgEff0p010";
  auto graph_btag_me = *(create_sig_eff_graph_run2(infile_btag, graphNamePostfix_btag_me, x_offset_btag, sf_file_btag_me, tagger_and_wp));
  graph_btag_me.SetLineWidth(line_width);
  graph_btag_me.SetLineColor(kBlue);
  graph_btag_me.SetLineStyle(line_style_btag);
  graph_btag_me.SetMarkerColor(kBlue);
  graph_btag_me.SetMarkerSize(marker_size);
  graph_btag_me.SetMarkerStyle(26);
  mg->Add(&graph_btag_me);

  cout << debug++ << endl;

  //const string graphName_btag_lo = kYears.at(year).name + "-BkgEff0p025/"+inputGraphName;
  //TGraph graph_btag_lo = *(TGraph*)infile_btag->Get(graphName_btag_lo.c_str());
  const string graphNamePostfix_btag_lo = "-BkgEff0p025/"+inputGraphName;
  TFile* sf_file_btag_lo = TFile::Open("scaleFactors_2023-08-10/scaleFactors-ak8_t_btagDJet__tau-BkgEff0p025/scaleFactors-ak8_t_btagDJet__tau-BkgEff0p025.root", "READ");
  tagger_and_wp = "ak8_t_btagDJet__tau-BkgEff0p025";
  auto graph_btag_lo = *(create_sig_eff_graph_run2(infile_btag, graphNamePostfix_btag_lo, x_offset_btag+x_offset_general, sf_file_btag_lo, tagger_and_wp));
  graph_btag_lo.SetLineWidth(line_width);
  graph_btag_lo.SetLineColor(kMagenta);
  graph_btag_lo.SetLineStyle(line_style_btag);
  graph_btag_lo.SetMarkerColor(kMagenta);
  graph_btag_lo.SetMarkerSize(marker_size);
  graph_btag_lo.SetMarkerStyle(32);
  mg->Add(&graph_btag_lo);

  cout << debug++ << endl;

  //const string graphName_btag_vl = kYears.at(year).name + "-BkgEff0p050/"+inputGraphName;
  //TGraph graph_btag_vl = *(TGraph*)infile_btag->Get(graphName_btag_vl.c_str());
  const string graphNamePostfix_btag_vl = "-BkgEff0p050/"+inputGraphName;
  TFile* sf_file_btag_vl = TFile::Open("scaleFactors_2023-08-10/scaleFactors-ak8_t_btagDJet__tau-BkgEff0p050/scaleFactors-ak8_t_btagDJet__tau-BkgEff0p050.root", "READ");
  tagger_and_wp = "ak8_t_btagDJet__tau-BkgEff0p050";
  auto graph_btag_vl = *(create_sig_eff_graph_run2(infile_btag, graphNamePostfix_btag_vl, x_offset_btag+2*x_offset_general, sf_file_btag_vl, tagger_and_wp));
  graph_btag_vl.SetLineWidth(line_width);
  graph_btag_vl.SetLineColor(kRed);
  graph_btag_vl.SetLineStyle(line_style_btag);
  graph_btag_vl.SetMarkerColor(kRed);
  graph_btag_vl.SetMarkerSize(marker_size);
  graph_btag_vl.SetMarkerStyle(46);
  mg->Add(&graph_btag_vl);

  cout << debug++ << endl;




  // HOTVR

  //const string graphName_hotvr = kYears.at(year).name + "-Standard/"+inputGraphName;
  //TGraph graph_hotvr = *(TGraph*)infile_hotvr->Get(graphName_hotvr.c_str());
  const string graphNamePostfix_hotvr = "-Standard/"+inputGraphName;
  TFile* sf_file_hotvr = TFile::Open("scaleFactors_2023-08-10/scaleFactors-hotvr_t__tau-Standard/scaleFactors-hotvr_t__tau-Standard.root", "READ");
  tagger_and_wp = "hotvr_t__tau-Standard";
  auto graph_hotvr = *(create_sig_eff_graph_run2(infile_hotvr, graphNamePostfix_hotvr, 0., sf_file_hotvr, tagger_and_wp));
  graph_hotvr.SetLineWidth(line_width);
  graph_hotvr.SetLineColor(kBlack);
  graph_hotvr.SetLineStyle(line_style_nominal);
  graph_hotvr.SetMarkerColor(kBlack);
  graph_hotvr.SetMarkerSize(marker_size);
  graph_hotvr.SetMarkerStyle(34);
  mg->Add(&graph_hotvr);

  cout << debug++ << endl;


  // Canvas

  TCanvas *c = new TCanvas("canvas", "canvas title", 600, 600);
  c->cd();
  float margin_l = 0.15;
  float margin_r = 0.05;
  float margin_b = 0.12;
  float margin_t = 0.08;
  float tick_length = 0.015; // fraction of canvas width/height
  c->SetMargin(margin_l, margin_r, margin_b, margin_t); //lrbt
  auto coord = new CoordinateConverter();
  coord->init(margin_l, margin_r, margin_b, margin_t);
  // if(log_y) c->SetLogy();
  c->SetTickx(1);
  c->SetTicky(1);


  bool plot_at_bottom = true;


  mg->Draw("alp");
  mg->SetTitle("");
  // mg->GetHistogram()->GetXaxis()->SetRangeUser(300., 1000.);
  mg->GetHistogram()->GetXaxis()->SetLimits(200., 1200.);
  double maximum_saved = mg->GetHistogram()->GetMaximum();
  double minimum_saved = mg->GetHistogram()->GetMinimum();
  double new_maximum = maximum_saved * 1.618;
  // mg->GetHistogram()->SetMaximum(1.);
  // mg->GetHistogram()->SetMaximum(new_maximum);
  mg->GetHistogram()->SetMaximum(1.4);
  mg->GetHistogram()->SetMinimum(0.);
  // if(is_eff_bkg || is_eff_sig) mg->GetHistogram()->GetXaxis()->SetTitle(tagger.scan_var_tlatex_axis.c_str());
  // else if(is_roc) mg->GetHistogram()->GetXaxis()->SetTitle("True positive rate, #varepsilon_{S}");
  // if(is_eff_bkg || is_roc) mg->GetHistogram()->GetYaxis()->SetTitle("False positive rate, #varepsilon_{B}");
  // if(is_eff_sig) mg->GetHistogram()->GetYaxis()->SetTitle("True positive rate, #varepsilon_{S}");
  mg->GetHistogram()->GetXaxis()->SetTitle("Probe jet #it{p}_{T} [GeV]");
  mg->GetHistogram()->GetYaxis()->SetTitle("Signal efficiency");
  mg->GetHistogram()->GetXaxis()->SetTitleSize(0.045);
  mg->GetHistogram()->GetYaxis()->SetTitleSize(0.045);
  mg->GetHistogram()->GetXaxis()->SetTitleOffset(1.2);
  mg->GetHistogram()->GetXaxis()->CenterTitle();
  mg->GetHistogram()->GetYaxis()->SetTitleOffset(1.4);
  mg->GetHistogram()->GetYaxis()->CenterTitle();
  // mg->GetHistogram()->GetXaxis()->SetNdivisions(405);
  // if(!log_y) mg->GetHistogram()->GetYaxis()->SetNdivisions(405);
  gStyle->SetLineWidth(1);

  // legend->Draw();





  // Lines
  // TLine(Double_t x1,Double_t y1,Double_t x2,Double_t y2)

  int lineX_color = kGray+2;
  int lineX_style = 3;
  int lineX_width = 1;

  // double line_max_y = maximum_saved * 1.05;
  double line_max_y = 1.;

  TLine* line1 = new TLine(250, 0.4, 250, 0.5);
  line1->SetLineColor(lineX_color);
  line1->SetLineStyle(lineX_style);
  line1->SetLineWidth(lineX_width);
  line1->Draw();
  TLine* line2 = new TLine(300, 0, 300, line_max_y);
  line2->SetLineColor(lineX_color);
  line2->SetLineStyle(lineX_style);
  line2->SetLineWidth(lineX_width);
  line2->Draw();
  TLine* line3 = new TLine(350, 0.4, 350, 0.5);
  line3->SetLineColor(lineX_color);
  line3->SetLineStyle(lineX_style);
  line3->SetLineWidth(lineX_width);
  line3->Draw();
  TLine* line4 = new TLine(400, 0, 400, line_max_y);
  line4->SetLineColor(lineX_color);
  line4->SetLineStyle(lineX_style);
  line4->SetLineWidth(lineX_width);
  line4->Draw();
  TLine* line5 = new TLine(480, 0, 480, line_max_y);
  line5->SetLineColor(lineX_color);
  line5->SetLineStyle(lineX_style);
  line5->SetLineWidth(lineX_width);
  line5->Draw();
  TLine* line6 = new TLine(600, 0, 600, line_max_y);
  line6->SetLineColor(lineX_color);
  line6->SetLineStyle(lineX_style);
  line6->SetLineWidth(lineX_width);
  line6->Draw();



  // Legend

  float leg_x1 = 0.62;
  float leg_y1 = 0.15;
  float leg_x2 = 0.83;
  float leg_y2 = 0.55;

  // TLegend *legend = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2, '', 'brNDC');
  TLegend *legend = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2);

  legend->AddEntry(&graph_nominal_vl, "AK8, #tau_{3}/#tau_{2} < 0.69", "p");
  legend->AddEntry(&graph_nominal_lo, "AK8, #tau_{3}/#tau_{2} < 0.61", "p");
  legend->AddEntry(&graph_hotvr, "HOTVR, #tau_{3}/#tau_{2} < 0.56", "p");
  legend->AddEntry(&graph_nominal_me, "AK8, #tau_{3}/#tau_{2} < 0.52", "p");
  legend->AddEntry(&graph_nominal_ti, "AK8, #tau_{3}/#tau_{2} < 0.47", "p");
  legend->AddEntry(&graph_nominal_vt, "AK8, #tau_{3}/#tau_{2} < 0.38", "p");

  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0); // Set background to transparent

  legend->Draw();


  // Legend nominal / subjet b tag

  float leg2_x1 = 0.44;
  float leg2_y1 = 0.72;
  float leg2_x2 = 0.77;
  float leg2_y2 = 0.87;

  // TLegend *legend2 = new TLegend(leg2_x1, leg2_y1, leg2_x2, leg2_y2, '', 'brNDC');
  TLegend *legend2 = new TLegend(leg2_x1, leg2_y1, leg2_x2, leg2_y2);

  TGraph *graph_nominal_dummy = new TGraph();
  graph_nominal_dummy->SetLineWidth(line_width);
  graph_nominal_dummy->SetLineColor(kBlack);
  graph_nominal_dummy->SetLineStyle(line_style_nominal);
  graph_nominal_dummy->SetMarkerColor(kBlack);
  graph_nominal_dummy->SetMarkerSize(marker_size);
  graph_nominal_dummy->SetMarkerStyle(20);

  TGraph *graph_btag_dummy = new TGraph();
  graph_btag_dummy->SetLineWidth(line_width);
  graph_btag_dummy->SetLineColor(kBlack);
  graph_btag_dummy->SetLineStyle(line_style_btag);
  graph_btag_dummy->SetMarkerColor(kBlack);
  graph_btag_dummy->SetMarkerSize(marker_size);
  graph_btag_dummy->SetMarkerStyle(24);

  legend2->SetHeader("#bf{t tagging using N-subjettiness}");
  legend2->AddEntry(graph_nominal_dummy, "w/o subjet b tagging");
  legend2->AddEntry(graph_btag_dummy, "w/ loose DeepJet subjet b tag");
  
  legend2->SetTextSize(0.03);
  legend2->SetBorderSize(0);
  legend2->SetFillStyle(0); // Set background to transparent

  legend2->Draw();



  // TLatex *text_top_left = new TLatex(margin_l, 1-(margin_t-0.01), (string("Ultra Legacy ")+kYears.at(year).nice_name).c_str());
  // text_top_left->SetTextAlign(11); // left bottom aligned
  // text_top_left->SetTextFont(42);
  // text_top_left->SetTextSize(0.035);
  // text_top_left->SetNDC();
  // text_top_left->Draw();

  // TLatex *algo_label = new TLatex(coord->ConvertGraphXToPadX(1-(0.05*c->GetWh()/c->GetWw())), coord->ConvertGraphYToPadY(0.95), "AK8 PUPPI");
  // algo_label->SetTextAlign(33); // right top
  // algo_label->SetTextFont(62);
  // algo_label->SetTextSize(0.05); // 0.05
  // algo_label->SetTextColor(kGray+1);
  // algo_label->SetNDC();
  // algo_label->Draw();

  // const double eta_pt_text_x = coord->ConvertGraphXToPadX(1-(0.05*c->GetWh()/c->GetWw()));
  // double eta_text_y = coord->ConvertGraphYToPadY(0.801);
  // // double pt_text_y = coord->ConvertGraphYToPadY(0.867);
  // if(!plot_at_bottom) {
  //   eta_text_y = coord->ConvertGraphYToPadY((0.05*c->GetWh()/c->GetWw()) + 0.066);
  //   // pt_text_y = coord->ConvertGraphYToPadY((0.05*c->GetWh()/c->GetWw()));
  // }

  // TLatex *eta_text = new TLatex(eta_pt_text_x, eta_text_y, "|#eta| < 2.5");
  // eta_text->SetTextAlign(plot_at_bottom ? 33 : 31); // right top
  // eta_text->SetTextFont(42);
  // eta_text->SetTextSize(0.035);
  // eta_text->SetNDC();
  // eta_text->Draw();

  // TLatex *pt_text = new TLatex(eta_pt_text_x, pt_text_y, pt_bin.text.c_str());
  // pt_text->SetTextAlign(plot_at_bottom ? 33 : 31); // right top
  // pt_text->SetTextFont(42);
  // pt_text->SetTextSize(0.035);
  // pt_text->SetNDC();
  // pt_text->Draw();

  // string string_text_top_right = "t#bar{t} #rightarrow all-jets @ #sqrt{#it{s}} = 13 TeV";
  // string string_text_top_right = "#sqrt{#it{s}} = 13 TeV";
  string string_text_top_right = kYears.at(Year::isRun2).lumi_fb_display+" fb^{#minus1} (13 TeV)";
  TLatex *text_top_right = new TLatex(1-margin_r, 1-(margin_t-0.01), string_text_top_right.c_str());
  text_top_right->SetTextAlign(31); // right bottom aligned
  text_top_right->SetTextFont(42);
  text_top_right->SetTextSize(0.035);
  text_top_right->SetNDC();
  text_top_right->Draw();

  TText *prelimPW1 = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.95), "Private work");
  prelimPW1->SetTextAlign(13); // left top
  prelimPW1->SetTextFont(52);
  prelimPW1->SetTextSize(0.03);
  prelimPW1->SetNDC();
  // prelimPW1->SetTextColor(kWhite);
  if(!bool_prelim) prelimPW1->Draw();

  TText *prelimPW2 = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.90), "(CMS data/simulation)");
  prelimPW2->SetTextAlign(13); // left top
  prelimPW2->SetTextFont(52);
  prelimPW2->SetTextSize(0.03);
  prelimPW2->SetNDC();
  // prelimPW2->SetTextColor(kWhite);
  if(!bool_prelim) prelimPW2->Draw();

  TLatex *cms = new TLatex(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.95), "CMS");
  cms->SetTextAlign(13); // left top
  cms->SetTextFont(62);
  cms->SetTextSize(0.05);
  cms->SetNDC();
  if(bool_prelim) cms->Draw();

  if(!plot_at_bottom) {
    TString prelim_text;
    if(bool_prelim) prelim_text = "Simulation, Preliminary";
    else prelim_text = "Simulation, Private Work";
    TText *prelim = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.87), prelim_text);
    prelim->SetTextAlign(13); // left top
    prelim->SetTextFont(52);
    prelim->SetTextSize(0.035);
    prelim->SetNDC();
    if(bool_prelim) prelim->Draw();
  }
  else {
    TText *prelim = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.87), "Simulation");
    prelim->SetTextAlign(13); // left top
    prelim->SetTextFont(52);
    prelim->SetTextSize(0.035);
    prelim->SetNDC();
    // prelim->SetTextColor(kWhite);
    if(bool_prelim) prelim->Draw();

    TString prelim_text;
    if(bool_prelim) prelim_text = "Preliminary";
    else prelim_text = "Private Work";
    TText *prelim2 = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.804), prelim_text);
    prelim2->SetTextAlign(13); // left top
    prelim2->SetTextFont(52);
    prelim2->SetTextSize(0.035);
    prelim2->SetNDC();
    // prelim2->SetTextColor(kWhite);
    if(bool_prelim) prelim2->Draw();
  }


  // https://root-forum.cern.ch/t/inconsistent-tick-length/18563/8
  c->Update();
  double tickScaleX = (c->GetUxmax() - c->GetUxmin()) / (c->GetX2() - c->GetX1()) * (c->GetWh() * c->GetAbsHNDC());
  double tickScaleY = (c->GetUymax() - c->GetUymin()) / (c->GetY2() - c->GetY1()) * (c->GetWw() * c->GetAbsWNDC());
  mg->GetHistogram()->GetYaxis()->SetTickLength(c->GetWw() * tick_length / tickScaleY);
  mg->GetHistogram()->GetXaxis()->SetTickLength(c->GetWh() * tick_length / tickScaleX);

  gPad->Update();
  gPad->RedrawAxis();

  // Save to disk
  string plotName;
  if(bool_prelim) plotName = "plot_sig_eff_in_mc__ak8__"+kYears.at(Year::isRun2).name+".pdf";
  else plotName = "plot_sig_eff_in_mc__ak8__"+kYears.at(Year::isRun2).name+"-PrivateWork.pdf";
  // string plotName = "plot__"+tagger.name_base+"__"+pt_bin.name+"__"+graph_base_name+"__";
  // for(bool do_x : {do_raw, do_msd, do_msd_btag}) { do_x ? plotName += "1" : plotName += "0"; }
  // plotName += (log_y ? string("log") : string("lin"))+".pdf";
  string plotPath = data_path+"/"+plotName;
  c->SaveAs(plotPath.c_str());
  delete c;
}

//void do_plot_partnet_w(const Year & year, const bool bool_prelim) {
void do_plot_partnet_w(const bool bool_prelim) {

  int debug = 0;

  const int line_width = 2;
  const int marker_size = 1;
  const int line_style_nominal = 1;
  const int line_style_btag = 2;

  TMultiGraph *mg = new TMultiGraph("mg", "mg title");

  cout << debug++ << endl;

  double x_offset_nominal = 4;
  double x_offset_btag = -x_offset_nominal;
  double x_offset_general = 16;

  //auto graph_nominal_BkgEff0p005 = TGraph(3, leo_ptCenter_W, values_leo_WvsQCD_BkgEff0p005.at(year).data());
  auto graph_nominal_BkgEff0p005 = create_sig_eff_graph_run2_partnet_w(E_leo_Tagger::isWvsQCD, E_leo_WP::isBkgEff0p005, x_offset_nominal);
  graph_nominal_BkgEff0p005.SetLineWidth(line_width);
  graph_nominal_BkgEff0p005.SetLineColor(kCyan);
  graph_nominal_BkgEff0p005.SetLineStyle(line_style_nominal);
  graph_nominal_BkgEff0p005.SetMarkerColor(kCyan);
  graph_nominal_BkgEff0p005.SetMarkerSize(marker_size);
  graph_nominal_BkgEff0p005.SetMarkerStyle(21);
  mg->Add(&graph_nominal_BkgEff0p005);

  cout << debug++ << endl;
  
  //auto graph_nominal_BkgEff0p010 = TGraph(3, leo_ptCenter_W, values_leo_WvsQCD_BkgEff0p010.at(year).data());
  auto graph_nominal_BkgEff0p010 = create_sig_eff_graph_run2_partnet_w(E_leo_Tagger::isWvsQCD, E_leo_WP::isBkgEff0p010, x_offset_nominal);
  graph_nominal_BkgEff0p010.SetLineWidth(line_width);
  graph_nominal_BkgEff0p010.SetLineColor(kBlue);
  graph_nominal_BkgEff0p010.SetLineStyle(line_style_nominal);
  graph_nominal_BkgEff0p010.SetMarkerColor(kBlue);
  graph_nominal_BkgEff0p010.SetMarkerSize(marker_size);
  graph_nominal_BkgEff0p010.SetMarkerStyle(22);
  mg->Add(&graph_nominal_BkgEff0p010);

  cout << debug++ << endl;

  //auto graph_nominal_BkgEff0p050 = TGraph(3, leo_ptCenter_W, values_leo_WvsQCD_BkgEff0p050.at(year).data());
  auto graph_nominal_BkgEff0p050 = create_sig_eff_graph_run2_partnet_w(E_leo_Tagger::isWvsQCD, E_leo_WP::isBkgEff0p050, x_offset_nominal);
  graph_nominal_BkgEff0p050.SetLineWidth(line_width);
  graph_nominal_BkgEff0p050.SetLineColor(kRed);
  graph_nominal_BkgEff0p050.SetLineStyle(line_style_nominal);
  graph_nominal_BkgEff0p050.SetMarkerColor(kRed);
  graph_nominal_BkgEff0p050.SetMarkerSize(marker_size);
  graph_nominal_BkgEff0p050.SetMarkerStyle(47);
  mg->Add(&graph_nominal_BkgEff0p050);

  cout << debug++ << endl;




  //auto graph_btag_BkgEff0p005 = TGraph(3, leo_ptCenter_W, values_leo_WvsQCDMD_BkgEff0p005.at(year).data());
  auto graph_btag_BkgEff0p005 = create_sig_eff_graph_run2_partnet_w(E_leo_Tagger::isWvsQCDMD, E_leo_WP::isBkgEff0p005, x_offset_btag);
  graph_btag_BkgEff0p005.SetLineWidth(line_width);
  graph_btag_BkgEff0p005.SetLineColor(kCyan);
  graph_btag_BkgEff0p005.SetLineStyle(line_style_btag);
  graph_btag_BkgEff0p005.SetMarkerColor(kCyan);
  graph_btag_BkgEff0p005.SetMarkerSize(marker_size);
  graph_btag_BkgEff0p005.SetMarkerStyle(25);
  mg->Add(&graph_btag_BkgEff0p005);

  cout << debug++ << endl;
  
  //auto graph_btag_BkgEff0p010 = TGraph(3, leo_ptCenter_W, values_leo_WvsQCDMD_BkgEff0p010.at(year).data());
  auto graph_btag_BkgEff0p010 = create_sig_eff_graph_run2_partnet_w(E_leo_Tagger::isWvsQCDMD, E_leo_WP::isBkgEff0p010, x_offset_btag);
  graph_btag_BkgEff0p010.SetLineWidth(line_width);
  graph_btag_BkgEff0p010.SetLineColor(kBlue);
  graph_btag_BkgEff0p010.SetLineStyle(line_style_btag);
  graph_btag_BkgEff0p010.SetMarkerColor(kBlue);
  graph_btag_BkgEff0p010.SetMarkerSize(marker_size);
  graph_btag_BkgEff0p010.SetMarkerStyle(26);
  mg->Add(&graph_btag_BkgEff0p010);

  cout << debug++ << endl;

  //auto graph_btag_BkgEff0p025 = TGraph(3, leo_ptCenter_W, values_leo_WvsQCDMD_BkgEff0p025.at(year).data());
  auto graph_btag_BkgEff0p025 = create_sig_eff_graph_run2_partnet_w(E_leo_Tagger::isWvsQCDMD, E_leo_WP::isBkgEff0p025, x_offset_btag);
  graph_btag_BkgEff0p025.SetLineWidth(line_width);
  graph_btag_BkgEff0p025.SetLineColor(kMagenta);
  graph_btag_BkgEff0p025.SetLineStyle(line_style_btag);
  graph_btag_BkgEff0p025.SetMarkerColor(kMagenta);
  graph_btag_BkgEff0p025.SetMarkerSize(marker_size);
  graph_btag_BkgEff0p025.SetMarkerStyle(32);
  mg->Add(&graph_btag_BkgEff0p025);

  cout << debug++ << endl;


  // Canvas

  TCanvas *c = new TCanvas("canvas", "canvas title", 600, 600);
  c->cd();
  float margin_l = 0.15;
  float margin_r = 0.05;
  float margin_b = 0.12;
  float margin_t = 0.08;
  float tick_length = 0.015; // fraction of canvas width/height
  c->SetMargin(margin_l, margin_r, margin_b, margin_t); //lrbt
  auto coord = new CoordinateConverter();
  coord->init(margin_l, margin_r, margin_b, margin_t);
  // if(log_y) c->SetLogy();
  c->SetTickx(1);
  c->SetTicky(1);


  bool plot_at_bottom = true;


  mg->Draw("alp");
  mg->SetTitle("");
  // mg->GetHistogram()->GetXaxis()->SetRangeUser(300., 1000.);
  mg->GetHistogram()->GetXaxis()->SetLimits(200., 1200.);
  double maximum_saved = mg->GetHistogram()->GetMaximum();
  double minimum_saved = mg->GetHistogram()->GetMinimum();
  double new_maximum = maximum_saved * 1.618;
  // mg->GetHistogram()->SetMaximum(1.);
  // mg->GetHistogram()->SetMaximum(new_maximum);
  mg->GetHistogram()->SetMaximum(1.4);
  mg->GetHistogram()->SetMinimum(0.);
  // if(is_eff_bkg || is_eff_sig) mg->GetHistogram()->GetXaxis()->SetTitle(tagger.scan_var_tlatex_axis.c_str());
  // else if(is_roc) mg->GetHistogram()->GetXaxis()->SetTitle("True positive rate, #varepsilon_{S}");
  // if(is_eff_bkg || is_roc) mg->GetHistogram()->GetYaxis()->SetTitle("False positive rate, #varepsilon_{B}");
  // if(is_eff_sig) mg->GetHistogram()->GetYaxis()->SetTitle("True positive rate, #varepsilon_{S}");
  mg->GetHistogram()->GetXaxis()->SetTitle("Probe jet #it{p}_{T} [GeV]");
  mg->GetHistogram()->GetYaxis()->SetTitle("Signal efficiency");
  mg->GetHistogram()->GetXaxis()->SetTitleSize(0.045);
  mg->GetHistogram()->GetYaxis()->SetTitleSize(0.045);
  mg->GetHistogram()->GetXaxis()->SetTitleOffset(1.2);
  mg->GetHistogram()->GetXaxis()->CenterTitle();
  mg->GetHistogram()->GetYaxis()->SetTitleOffset(1.4);
  mg->GetHistogram()->GetYaxis()->CenterTitle();
  // mg->GetHistogram()->GetXaxis()->SetNdivisions(405);
  // if(!log_y) mg->GetHistogram()->GetYaxis()->SetNdivisions(405);
  gStyle->SetLineWidth(1);

  // legend->Draw();





  // Lines
  // TLine(Double_t x1,Double_t y1,Double_t x2,Double_t y2)

  int lineX_color = kGray+2;
  int lineX_style = 3;
  int lineX_width = 1;

  // double line_max_y = maximum_saved * 1.05;
  double line_max_y = 1.;

  TLine* line2 = new TLine(300, 0, 300, line_max_y);
  line2->SetLineColor(lineX_color);
  line2->SetLineStyle(lineX_style);
  line2->SetLineWidth(lineX_width);
  line2->Draw();
  TLine* line5 = new TLine(400, 0, 400, line_max_y);
  line5->SetLineColor(lineX_color);
  line5->SetLineStyle(lineX_style);
  line5->SetLineWidth(lineX_width);
  line5->Draw();
  TLine* line6 = new TLine(800, 0, 800, line_max_y);
  line6->SetLineColor(lineX_color);
  line6->SetLineStyle(lineX_style);
  line6->SetLineWidth(lineX_width);
  line6->Draw();



  // Legend

  float leg_x1 = 0.62;
  float leg_y1 = 0.27;
  float leg_x2 = 0.83;
  float leg_y2 = 0.57;

  // TLegend *legend = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2, '', 'brNDC');
  TLegend *legend = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2);

  legend->AddEntry(&graph_nominal_BkgEff0p050, "Mistag rate = 5.0%", "p");
  legend->AddEntry(&graph_btag_BkgEff0p025, "Mistag rate = 2.5%", "p");
  legend->AddEntry(&graph_nominal_BkgEff0p010, "Mistag rate = 1.0%", "p");
  legend->AddEntry(&graph_nominal_BkgEff0p005, "Mistag rate = 0.5%", "p");

  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0); // Set background to transparent

  legend->Draw();


  // Legend nominal / subjet b tag

  float leg2_x1 = 0.44;
  float leg2_y1 = 0.72;
  float leg2_x2 = 0.77;
  float leg2_y2 = 0.87;

  // TLegend *legend2 = new TLegend(leg2_x1, leg2_y1, leg2_x2, leg2_y2, '', 'brNDC');
  TLegend *legend2 = new TLegend(leg2_x1, leg2_y1, leg2_x2, leg2_y2);

  TGraph *graph_nominal_dummy = new TGraph();
  graph_nominal_dummy->SetLineWidth(line_width);
  graph_nominal_dummy->SetLineColor(kBlack);
  graph_nominal_dummy->SetLineStyle(line_style_nominal);
  graph_nominal_dummy->SetMarkerColor(kBlack);
  graph_nominal_dummy->SetMarkerSize(marker_size);
  graph_nominal_dummy->SetMarkerStyle(20);

  TGraph *graph_btag_dummy = new TGraph();
  graph_btag_dummy->SetLineWidth(line_width);
  graph_btag_dummy->SetLineColor(kBlack);
  graph_btag_dummy->SetLineStyle(line_style_btag);
  graph_btag_dummy->SetMarkerColor(kBlack);
  graph_btag_dummy->SetMarkerSize(marker_size);
  graph_btag_dummy->SetMarkerStyle(24);

  legend2->SetHeader("#bf{W tagging using ParticleNet}");
  legend2->AddEntry(graph_nominal_dummy, "nominal: WvsQCD");
  legend2->AddEntry(graph_btag_dummy, "mass-decorrelated: WvsQCD-MD");
  
  legend2->SetTextSize(0.03);
  legend2->SetBorderSize(0);
  legend2->SetFillStyle(0); // Set background to transparent

  legend2->Draw();



  // TLatex *text_top_left = new TLatex(margin_l, 1-(margin_t-0.01), (string("Ultra Legacy ")+kYears.at(year).nice_name).c_str());
  // text_top_left->SetTextAlign(11); // left bottom aligned
  // text_top_left->SetTextFont(42);
  // text_top_left->SetTextSize(0.035);
  // text_top_left->SetNDC();
  // text_top_left->Draw();

  // TLatex *algo_label = new TLatex(coord->ConvertGraphXToPadX(1-(0.05*c->GetWh()/c->GetWw())), coord->ConvertGraphYToPadY(0.95), "AK8 PUPPI");
  // algo_label->SetTextAlign(33); // right top
  // algo_label->SetTextFont(62);
  // algo_label->SetTextSize(0.05); // 0.05
  // algo_label->SetTextColor(kGray+1);
  // algo_label->SetNDC();
  // algo_label->Draw();

  // const double eta_pt_text_x = coord->ConvertGraphXToPadX(1-(0.05*c->GetWh()/c->GetWw()));
  // double eta_text_y = coord->ConvertGraphYToPadY(0.801);
  // // double pt_text_y = coord->ConvertGraphYToPadY(0.867);
  // if(!plot_at_bottom) {
  //   eta_text_y = coord->ConvertGraphYToPadY((0.05*c->GetWh()/c->GetWw()) + 0.066);
  //   // pt_text_y = coord->ConvertGraphYToPadY((0.05*c->GetWh()/c->GetWw()));
  // }

  // TLatex *eta_text = new TLatex(eta_pt_text_x, eta_text_y, "|#eta| < 2.5");
  // eta_text->SetTextAlign(plot_at_bottom ? 33 : 31); // right top
  // eta_text->SetTextFont(42);
  // eta_text->SetTextSize(0.035);
  // eta_text->SetNDC();
  // eta_text->Draw();

  // TLatex *pt_text = new TLatex(eta_pt_text_x, pt_text_y, pt_bin.text.c_str());
  // pt_text->SetTextAlign(plot_at_bottom ? 33 : 31); // right top
  // pt_text->SetTextFont(42);
  // pt_text->SetTextSize(0.035);
  // pt_text->SetNDC();
  // pt_text->Draw();

  // string string_text_top_right = "t#bar{t} #rightarrow all-jets @ #sqrt{#it{s}} = 13 TeV";
  // string string_text_top_right = "#sqrt{#it{s}} = 13 TeV";
  string string_text_top_right = kYears.at(Year::isRun2).lumi_fb_display+" fb^{#minus1} (13 TeV)";
  TLatex *text_top_right = new TLatex(1-margin_r, 1-(margin_t-0.01), string_text_top_right.c_str());
  text_top_right->SetTextAlign(31); // right bottom aligned
  text_top_right->SetTextFont(42);
  text_top_right->SetTextSize(0.035);
  text_top_right->SetNDC();
  text_top_right->Draw();

  TText *prelimPW1 = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.95), "Private work");
  prelimPW1->SetTextAlign(13); // left top
  prelimPW1->SetTextFont(52);
  prelimPW1->SetTextSize(0.03);
  prelimPW1->SetNDC();
  // prelimPW1->SetTextColor(kWhite);
  if(!bool_prelim) prelimPW1->Draw();

  TText *prelimPW2 = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.90), "(CMS data/simulation)");
  prelimPW2->SetTextAlign(13); // left top
  prelimPW2->SetTextFont(52);
  prelimPW2->SetTextSize(0.03);
  prelimPW2->SetNDC();
  // prelimPW2->SetTextColor(kWhite);
  if(!bool_prelim) prelimPW2->Draw();

  TLatex *cms = new TLatex(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.95), "CMS");
  cms->SetTextAlign(13); // left top
  cms->SetTextFont(62);
  cms->SetTextSize(0.05);
  cms->SetNDC();
  if(bool_prelim) cms->Draw();

  if(!plot_at_bottom) {
    TString prelim_text;
    if(bool_prelim) prelim_text = "Simulation, Preliminary";
    else prelim_text = "Simulation, Private Work";
    TText *prelim = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.87), prelim_text);
    prelim->SetTextAlign(13); // left top
    prelim->SetTextFont(52);
    prelim->SetTextSize(0.035);
    prelim->SetNDC();
    if(bool_prelim) prelim->Draw();
  }
  else {
    TText *prelim = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.87), "Simulation");
    prelim->SetTextAlign(13); // left top
    prelim->SetTextFont(52);
    prelim->SetTextSize(0.035);
    prelim->SetNDC();
    // prelim->SetTextColor(kWhite);
    if(bool_prelim) prelim->Draw();

    TString prelim_text;
    if(bool_prelim) prelim_text = "Preliminary";
    else prelim_text = "Private Work";
    TText *prelim2 = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.804), prelim_text);
    prelim2->SetTextAlign(13); // left top
    prelim2->SetTextFont(52);
    prelim2->SetTextSize(0.035);
    prelim2->SetNDC();
    // prelim2->SetTextColor(kWhite);
    if(bool_prelim) prelim2->Draw();
  }


  // https://root-forum.cern.ch/t/inconsistent-tick-length/18563/8
  c->Update();
  double tickScaleX = (c->GetUxmax() - c->GetUxmin()) / (c->GetX2() - c->GetX1()) * (c->GetWh() * c->GetAbsHNDC());
  double tickScaleY = (c->GetUymax() - c->GetUymin()) / (c->GetY2() - c->GetY1()) * (c->GetWw() * c->GetAbsWNDC());
  mg->GetHistogram()->GetYaxis()->SetTickLength(c->GetWw() * tick_length / tickScaleY);
  mg->GetHistogram()->GetXaxis()->SetTickLength(c->GetWh() * tick_length / tickScaleX);

  gPad->Update();
  gPad->RedrawAxis();

  // Save to disk
  string plotName;
  if(bool_prelim) plotName = "plot_sig_eff_in_mc__partnet_w__"+kYears.at(Year::isRun2).name+".pdf";
  else plotName = "plot_sig_eff_in_mc__partnet_w__"+kYears.at(Year::isRun2).name+"-PrivateWork.pdf";
  string plotPath = data_path+"/"+plotName;
  c->SaveAs(plotPath.c_str());
  delete c;
}

void do_plot_partnet_top(const Year & year, const bool bool_prelim) {
//void do_plot_partnet_top(const bool bool_prelim) {

  int debug = 0;

  const int line_width = 2;
  const int marker_size = 1;
  const int line_style_nominal = 1;
  const int line_style_btag = 2;

  TMultiGraph *mg = new TMultiGraph("mg", "mg title");

  cout << debug++ << endl;

  TGraph graph_nominal_BkgEff0p001 = TGraph(4, leo_ptCenter_top, values_leo_TvsQCD_BkgEff0p001.at(year).data());
  graph_nominal_BkgEff0p001.SetLineWidth(line_width);
  graph_nominal_BkgEff0p001.SetLineColor(kGreen);
  graph_nominal_BkgEff0p001.SetLineStyle(line_style_nominal);
  graph_nominal_BkgEff0p001.SetMarkerColor(kGreen);
  graph_nominal_BkgEff0p001.SetMarkerSize(marker_size);
  graph_nominal_BkgEff0p001.SetMarkerStyle(20);
  mg->Add(&graph_nominal_BkgEff0p001);

  cout << debug++ << endl;
  
  TGraph graph_nominal_BkgEff0p005 = TGraph(4, leo_ptCenter_top, values_leo_TvsQCD_BkgEff0p005.at(year).data());
  graph_nominal_BkgEff0p005.SetLineWidth(line_width);
  graph_nominal_BkgEff0p005.SetLineColor(kCyan);
  graph_nominal_BkgEff0p005.SetLineStyle(line_style_nominal);
  graph_nominal_BkgEff0p005.SetMarkerColor(kCyan);
  graph_nominal_BkgEff0p005.SetMarkerSize(marker_size);
  graph_nominal_BkgEff0p005.SetMarkerStyle(21);
  mg->Add(&graph_nominal_BkgEff0p005);

  cout << debug++ << endl;

  TGraph graph_nominal_BkgEff0p010 = TGraph(4, leo_ptCenter_top, values_leo_TvsQCD_BkgEff0p010.at(year).data());
  graph_nominal_BkgEff0p010.SetLineWidth(line_width);
  graph_nominal_BkgEff0p010.SetLineColor(kBlue);
  graph_nominal_BkgEff0p010.SetLineStyle(line_style_nominal);
  graph_nominal_BkgEff0p010.SetMarkerColor(kBlue);
  graph_nominal_BkgEff0p010.SetMarkerSize(marker_size);
  graph_nominal_BkgEff0p010.SetMarkerStyle(22);
  mg->Add(&graph_nominal_BkgEff0p010);

  cout << debug++ << endl;





  // Canvas

  TCanvas *c = new TCanvas("canvas", "canvas title", 600, 600);
  c->cd();
  float margin_l = 0.15;
  float margin_r = 0.05;
  float margin_b = 0.12;
  float margin_t = 0.08;
  float tick_length = 0.015; // fraction of canvas width/height
  c->SetMargin(margin_l, margin_r, margin_b, margin_t); //lrbt
  auto coord = new CoordinateConverter();
  coord->init(margin_l, margin_r, margin_b, margin_t);
  // if(log_y) c->SetLogy();
  c->SetTickx(1);
  c->SetTicky(1);


  bool plot_at_bottom = true;


  mg->Draw("alp");
  mg->SetTitle("");
  // mg->GetHistogram()->GetXaxis()->SetRangeUser(300., 1000.);
  mg->GetHistogram()->GetXaxis()->SetLimits(200., 1200.);
  double maximum_saved = mg->GetHistogram()->GetMaximum();
  double minimum_saved = mg->GetHistogram()->GetMinimum();
  double new_maximum = maximum_saved * 1.618;
  // mg->GetHistogram()->SetMaximum(1.);
  // mg->GetHistogram()->SetMaximum(new_maximum);
  mg->GetHistogram()->SetMaximum(1.4);
  mg->GetHistogram()->SetMinimum(0.);
  // if(is_eff_bkg || is_eff_sig) mg->GetHistogram()->GetXaxis()->SetTitle(tagger.scan_var_tlatex_axis.c_str());
  // else if(is_roc) mg->GetHistogram()->GetXaxis()->SetTitle("True positive rate, #varepsilon_{S}");
  // if(is_eff_bkg || is_roc) mg->GetHistogram()->GetYaxis()->SetTitle("False positive rate, #varepsilon_{B}");
  // if(is_eff_sig) mg->GetHistogram()->GetYaxis()->SetTitle("True positive rate, #varepsilon_{S}");
  mg->GetHistogram()->GetXaxis()->SetTitle("Probe jet #it{p}_{T} [GeV]");
  mg->GetHistogram()->GetYaxis()->SetTitle("Signal efficiency");
  mg->GetHistogram()->GetXaxis()->SetTitleSize(0.045);
  mg->GetHistogram()->GetYaxis()->SetTitleSize(0.045);
  mg->GetHistogram()->GetXaxis()->SetTitleOffset(1.2);
  mg->GetHistogram()->GetXaxis()->CenterTitle();
  mg->GetHistogram()->GetYaxis()->SetTitleOffset(1.4);
  mg->GetHistogram()->GetYaxis()->CenterTitle();
  // mg->GetHistogram()->GetXaxis()->SetNdivisions(405);
  // if(!log_y) mg->GetHistogram()->GetYaxis()->SetNdivisions(405);
  gStyle->SetLineWidth(1);

  // legend->Draw();





  // Lines
  // TLine(Double_t x1,Double_t y1,Double_t x2,Double_t y2)

  int lineX_color = kGray+2;
  int lineX_style = 3;
  int lineX_width = 1;

  // double line_max_y = maximum_saved * 1.05;
  double line_max_y = 1.;

  TLine* line2 = new TLine(300, 0, 300, line_max_y);
  line2->SetLineColor(lineX_color);
  line2->SetLineStyle(lineX_style);
  line2->SetLineWidth(lineX_width);
  line2->Draw();
  TLine* line4 = new TLine(400, 0, 400, line_max_y);
  line4->SetLineColor(lineX_color);
  line4->SetLineStyle(lineX_style);
  line4->SetLineWidth(lineX_width);
  line4->Draw();
  TLine* line5 = new TLine(480, 0, 480, line_max_y);
  line5->SetLineColor(lineX_color);
  line5->SetLineStyle(lineX_style);
  line5->SetLineWidth(lineX_width);
  line5->Draw();
  TLine* line6 = new TLine(600, 0, 600, line_max_y);
  line6->SetLineColor(lineX_color);
  line6->SetLineStyle(lineX_style);
  line6->SetLineWidth(lineX_width);
  line6->Draw();



  // Legend

  float leg_x1 = 0.62;
  float leg_y1 = 0.48;
  float leg_x2 = 0.83;
  float leg_y2 = 0.68;

  // TLegend *legend = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2, '', 'brNDC');
  TLegend *legend = new TLegend(leg_x1, leg_y1, leg_x2, leg_y2);

  legend->AddEntry(&graph_nominal_BkgEff0p010, "Mistag rate = 1.0%", "p");
  legend->AddEntry(&graph_nominal_BkgEff0p005, "Mistag rate = 0.5%", "p");
  legend->AddEntry(&graph_nominal_BkgEff0p001, "Mistag rate = 0.1%", "p");

  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0); // Set background to transparent

  legend->Draw();


  // Legend nominal / subjet b tag

  float leg2_x1 = 0.44;
  float leg2_y1 = 0.72;
  float leg2_x2 = 0.77;
  float leg2_y2 = 0.87;

  // TLegend *legend2 = new TLegend(leg2_x1, leg2_y1, leg2_x2, leg2_y2, '', 'brNDC');
  TLegend *legend2 = new TLegend(leg2_x1, leg2_y1, leg2_x2, leg2_y2);

  TGraph *graph_nominal_dummy = new TGraph();
  graph_nominal_dummy->SetLineWidth(line_width);
  graph_nominal_dummy->SetLineColor(kBlack);
  graph_nominal_dummy->SetLineStyle(line_style_nominal);
  graph_nominal_dummy->SetMarkerColor(kBlack);
  graph_nominal_dummy->SetMarkerSize(marker_size);
  graph_nominal_dummy->SetMarkerStyle(20);

  TGraph *graph_btag_dummy = new TGraph();
  graph_btag_dummy->SetLineWidth(line_width);
  graph_btag_dummy->SetLineColor(kBlack);
  graph_btag_dummy->SetLineStyle(line_style_btag);
  graph_btag_dummy->SetMarkerColor(kBlack);
  graph_btag_dummy->SetMarkerSize(marker_size);
  graph_btag_dummy->SetMarkerStyle(24);

  legend2->SetHeader("#bf{t tagging using ParticleNet TvsQCD}");
  legend2->AddEntry(graph_nominal_dummy, "", "");
  legend2->AddEntry(graph_btag_dummy, "", "");
  
  legend2->SetTextSize(0.03);
  legend2->SetBorderSize(0);
  legend2->SetFillStyle(0); // Set background to transparent

  legend2->Draw();



  // TLatex *text_top_left = new TLatex(margin_l, 1-(margin_t-0.01), (string("Ultra Legacy ")+kYears.at(year).nice_name).c_str());
  // text_top_left->SetTextAlign(11); // left bottom aligned
  // text_top_left->SetTextFont(42);
  // text_top_left->SetTextSize(0.035);
  // text_top_left->SetNDC();
  // text_top_left->Draw();

  // TLatex *algo_label = new TLatex(coord->ConvertGraphXToPadX(1-(0.05*c->GetWh()/c->GetWw())), coord->ConvertGraphYToPadY(0.95), "AK8 PUPPI");
  // algo_label->SetTextAlign(33); // right top
  // algo_label->SetTextFont(62);
  // algo_label->SetTextSize(0.05); // 0.05
  // algo_label->SetTextColor(kGray+1);
  // algo_label->SetNDC();
  // algo_label->Draw();

  // const double eta_pt_text_x = coord->ConvertGraphXToPadX(1-(0.05*c->GetWh()/c->GetWw()));
  // double eta_text_y = coord->ConvertGraphYToPadY(0.801);
  // // double pt_text_y = coord->ConvertGraphYToPadY(0.867);
  // if(!plot_at_bottom) {
  //   eta_text_y = coord->ConvertGraphYToPadY((0.05*c->GetWh()/c->GetWw()) + 0.066);
  //   // pt_text_y = coord->ConvertGraphYToPadY((0.05*c->GetWh()/c->GetWw()));
  // }

  // TLatex *eta_text = new TLatex(eta_pt_text_x, eta_text_y, "|#eta| < 2.5");
  // eta_text->SetTextAlign(plot_at_bottom ? 33 : 31); // right top
  // eta_text->SetTextFont(42);
  // eta_text->SetTextSize(0.035);
  // eta_text->SetNDC();
  // eta_text->Draw();

  // TLatex *pt_text = new TLatex(eta_pt_text_x, pt_text_y, pt_bin.text.c_str());
  // pt_text->SetTextAlign(plot_at_bottom ? 33 : 31); // right top
  // pt_text->SetTextFont(42);
  // pt_text->SetTextSize(0.035);
  // pt_text->SetNDC();
  // pt_text->Draw();

  // string string_text_top_right = "t#bar{t} #rightarrow all-jets @ #sqrt{#it{s}} = 13 TeV";
  // string string_text_top_right = "#sqrt{#it{s}} = 13 TeV";
  string string_text_top_right = kYears.at(Year::isRun2).lumi_fb_display+" fb^{#minus1} (13 TeV)";
  TLatex *text_top_right = new TLatex(1-margin_r, 1-(margin_t-0.01), string_text_top_right.c_str());
  text_top_right->SetTextAlign(31); // right bottom aligned
  text_top_right->SetTextFont(42);
  text_top_right->SetTextSize(0.035);
  text_top_right->SetNDC();
  text_top_right->Draw();

  TText *prelimPW1 = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.95), "Private work");
  prelimPW1->SetTextAlign(13); // left top
  prelimPW1->SetTextFont(52);
  prelimPW1->SetTextSize(0.03);
  prelimPW1->SetNDC();
  // prelimPW1->SetTextColor(kWhite);
  if(!bool_prelim) prelimPW1->Draw();

  TText *prelimPW2 = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.90), "(CMS data/simulation)");
  prelimPW2->SetTextAlign(13); // left top
  prelimPW2->SetTextFont(52);
  prelimPW2->SetTextSize(0.03);
  prelimPW2->SetNDC();
  // prelimPW2->SetTextColor(kWhite);
  if(!bool_prelim) prelimPW2->Draw();

  TLatex *cms = new TLatex(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.95), "CMS");
  cms->SetTextAlign(13); // left top
  cms->SetTextFont(62);
  cms->SetTextSize(0.05);
  cms->SetNDC();
  if(bool_prelim) cms->Draw();

  if(!plot_at_bottom) {
    TString prelim_text;
    if(bool_prelim) prelim_text = "Simulation, Preliminary";
    else prelim_text = "Simulation, Private Work";
    TText *prelim = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.87), prelim_text);
    prelim->SetTextAlign(13); // left top
    prelim->SetTextFont(52);
    prelim->SetTextSize(0.035);
    prelim->SetNDC();
    if(bool_prelim) prelim->Draw();
  }
  else {
    TText *prelim = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.87), "Simulation");
    prelim->SetTextAlign(13); // left top
    prelim->SetTextFont(52);
    prelim->SetTextSize(0.035);
    prelim->SetNDC();
    // prelim->SetTextColor(kWhite);
    if(bool_prelim) prelim->Draw();

    TString prelim_text;
    if(bool_prelim) prelim_text = "Preliminary";
    else prelim_text = "Private Work";
    TText *prelim2 = new TText(coord->ConvertGraphXToPadX(0.05), coord->ConvertGraphYToPadY(0.804), prelim_text);
    prelim2->SetTextAlign(13); // left top
    prelim2->SetTextFont(52);
    prelim2->SetTextSize(0.035);
    prelim2->SetNDC();
    // prelim2->SetTextColor(kWhite);
    if(bool_prelim) prelim2->Draw();
  }


  // https://root-forum.cern.ch/t/inconsistent-tick-length/18563/8
  c->Update();
  double tickScaleX = (c->GetUxmax() - c->GetUxmin()) / (c->GetX2() - c->GetX1()) * (c->GetWh() * c->GetAbsHNDC());
  double tickScaleY = (c->GetUymax() - c->GetUymin()) / (c->GetY2() - c->GetY1()) * (c->GetWw() * c->GetAbsWNDC());
  mg->GetHistogram()->GetYaxis()->SetTickLength(c->GetWw() * tick_length / tickScaleY);
  mg->GetHistogram()->GetXaxis()->SetTickLength(c->GetWh() * tick_length / tickScaleX);

  gPad->Update();
  gPad->RedrawAxis();

  // Save to disk
  string plotName;
  if(bool_prelim) plotName = "plot_sig_eff_in_mc__partnet_top__"+kYears.at(Year::isRun2).name+".pdf";
  else plotName = "plot_sig_eff_in_mc__partnet_top__"+kYears.at(Year::isRun2).name+"-PrivateWork.pdf";
  string plotPath = data_path+"/"+plotName;
  c->SaveAs(plotPath.c_str());
  delete c;
}



void sig_eff_in_mc_plotting() {

    //do_plot_nsubjettiness(true);
    //do_plot_nsubjettiness(false);
    do_plot_partnet_w(true);
    do_plot_partnet_w(false);
    //do_plot_partnet_top(true);
    //do_plot_partnet_top(false);

    for(const auto year : kAllYears) {
      //do_plot_nsubjettiness(year, true);
      //do_plot_nsubjettiness(year, false);
      //do_plot_partnet_w(year, true);
      //do_plot_partnet_w(year, false);
      //do_plot_partnet_top(year, true);
      //do_plot_partnet_top(year, false);
  }

}