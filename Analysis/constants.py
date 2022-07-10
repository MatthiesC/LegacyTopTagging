from __future__ import print_function

from collections import OrderedDict
import numpy as np


# Source: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopSystematics#Luminosity
# See also: https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM#SummaryTable (includes links to CMS publications which you can cite)
_YEARS = OrderedDict([
    ('UL16', {
        'short_name': 'UL16',
        'year': '2016',
        'long_name': 'Ultra Legacy 2016',
        'lumi_fb': 36.333, # TopSystematics source and my own brilcalc give the same results (36.333380073.916976929 with PHYSICS normtag)
        'lumi_pb': 36333.,
        'lumi_unc': 0.012,
        'lumi_fb_display': '36.3',
    }),
    ('UL16preVFP', {
        'short_name': 'UL16preVFP',
        'year': '2016 (early)',
        'long_name': 'Ultra Legacy 2016 (early)',
        # 'lumi_fb': 19.52,
        # 'lumi_pb': 19520.,
        'lumi_fb': 19.536, # brilcalc gives 19.536411965198516846 with PHYSICS normtag
        'lumi_pb': 19536.,
        'lumi_unc': 0.012,
        'lumi_fb_display': '19.5',
    }),
    ('UL16postVFP', {
        'short_name': 'UL16postVFP',
        'year': '2016 (late)',
        'long_name': 'Ultra Legacy 2016 (late)',
        # 'lumi_fb': 16.81,
        # 'lumi_pb': 16810.,
        'lumi_fb': 16.797, # brilcalc gives 16.796968109 with PHYSICS normtag
        'lumi_pb': 16797.,
        'lumi_unc': 0.012,
        'lumi_fb_display': '16.8',
    }),
    ('UL17', {
        'short_name': 'UL17',
        'year': '2017',
        'long_name': 'Ultra Legacy 2017',
        'lumi_fb': 41.480, # TopSystematics source and my own brilcalc give the same results (41.4796805287616 with PHYSICS normtag)
        'lumi_pb': 41480.,
        'lumi_unc': 0.023,
        'lumi_fb_display': '41.5',
    }),
    ('UL18', {
        'short_name': 'UL18',
        'year': '2018',
        'long_name': 'Ultra Legacy 2018',
        'lumi_fb': 59.832, # TopSystematics source and my own brilcalc give the same results (59.8324753390886 with PHYSICS normtag)
        'lumi_pb': 59832.,
        'lumi_unc': 0.025,
        'lumi_fb_display': '59.8',
    }),
    ('run2', {
        'short_name': 'ULRunII',
        'year': '2016 + 2017 + 2018',
        'long_name': 'Run II Ultra Legacy',
        'lumi_fb': 137.645, # sum from above brilcalc calculations: 137.645535941767
        'lumi_pb': 137645.,
        'lumi_unc': 0.016,
        'lumi_fb_display': '138',
    }),
])

_DEEPCSV_WPS = {
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation106XUL16preVFP
    'UL16preVFP': {
        'loose': 0.2027,
        'medium': 0.6001,
        'tight': 0.8819,
    },
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation106XUL16postVFP
    'UL16postVFP': {
        'loose': 0.1918,
        'medium': 0.5847,
        'tight': 0.8767,
    },
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation106XUL17
    'UL17': {
        'loose': 0.1355,
        'medium': 0.4506,
        'tight': 0.7738,
    },
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation106XUL18
    'UL18': {
        'loose': 0.1208,
        'medium': 0.4168,
        'tight': 0.7665,
    }
}

_DEEPJET_WPS = {
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation106XUL16preVFP
    'UL16preVFP': {
        'loose': 0.0508,
        'medium': 0.2598,
        'tight': 0.6502,
    },
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation106XUL16postVFP
    'UL16postVFP': {
        'loose': 0.0480,
        'medium': 0.2489,
        'tight': 0.6377,
    },
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation106XUL17
    'UL17': {
        'loose': 0.0532,
        'medium': 0.3040,
        'tight': 0.7476,
    },
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation106XUL18
    'UL18': {
        'loose': 0.0490,
        'medium': 0.2783,
        'tight': 0.7100,
    }
}

class VarInterval:

    def __init__(self, var, var_min=None, var_max=None, inversed=False):
        self.var = var
        self.var_min = var_min or 0
        self.var_max = var_max or 999999
        self.name = '_'.join([var, str(var_min)+'to'+str(var_max or 'Inf')])
        self.inversed = inversed # defines in which direction to scan; e.g. either 1 to 0 or 0 to 1

_PT_INTERVALS = [
    VarInterval('pt', 200),
    VarInterval('pt', 200, 250),
    VarInterval('pt', 250),
    VarInterval('pt', 250, 300),
    VarInterval('pt', 300),
    VarInterval('pt', 300, 350),
    VarInterval('pt', 300, 400),
    VarInterval('pt', 350),
    VarInterval('pt', 350, 400),
    VarInterval('pt', 400),
    VarInterval('pt', 400, 480),
    VarInterval('pt', 480),
    VarInterval('pt', 480, 600),
    VarInterval('pt', 600),
    VarInterval('pt', 600, 800),
    VarInterval('pt', 800),
    VarInterval('pt', 800, 1000),
    VarInterval('pt', 1000),
    VarInterval('pt', 1000, 1500),
    VarInterval('pt', 1500),
    # VarInterval('pt', 600, 620), # just for quick tests
]
_PT_INTERVALS = {pt_interval.name: pt_interval for pt_interval in _PT_INTERVALS}

class Tagger:

    def __init__(self, name, scan_var, tag_rule = []):
        self.name = name
        self.scan_var = scan_var
        self.tag_rule = tag_rule
        self.var_intervals = _PT_INTERVALS # just use all for now

    def set_tag_rule(self, tag_rule):
        self.tag_rule = tag_rule

    def get_tag_rule(self):
        if len(tag_rule):
            return ' & '.join(['('+x+')' for x in self.tag_rule])
        else:
            return 'True'

    def get_pt_min(self):
        self.pt_min = np.inf
        for pt_interval in self.var_intervals:
            if not pt_interval.var.startswith('pt'): continue
            self.pt_min = min(self.pt_min, pt_interval.var_min)
        return self.pt_min

# _TAGGERS = {
#     'hotvr_t_incl__tau': Tagger(
#         'tau32',
#         [],
#     ),
#     'hotvr_t__tau': Tagger(
#         'tau32',
#         ['mass > 140', 'mass < 220', 'fpt1 < 0.8', 'nsub > 2', 'mpair > 50'],
#     ),
#     'ak8_t_incl__tau': Tagger(
#         'tau32',
#         [],
#     ),
#     'ak8_t__tau': Tagger(
#         'tau32',
#         ['msd > 105', 'msd < 210'],
#     ),
#     'ak8_t_btagDJet__tau': Tagger(
#         'tau32',
#         ['msd > 105', 'msd < 210', 'subdeepjet > '+str(_DEEPJET_WPS[year]['loose'])],
#     ),
#     'ak8_t_btagDCSV__tau': Tagger(
#         'tau32',
#         ['msd > 105', 'msd < 210', 'subdeepcsv > '+str(_DEEPCSV_WPS[year]['loose'])],
#     ),
#     'ak8_t__partnet': Tagger(
#         'partnet_TvsQCD',
#         ['msd > 105', 'msd < 210'],
#     ),
#     'ak8_t__deepak8': Tagger(
#         'deepak8_TvsQCD',
#         ['msd > 105', 'msd < 210'],
#     ),
#     'ak8_w_incl__tau': Tagger(
#         'tau21',
#         [],
#     ),
#     'ak8_w__tau': Tagger(
#         'tau21',
#         ['msd > 65', 'msd < 105'],
#     ),
#     'ak8_w_incl__partnet': Tagger(
#         'partnet_WvsQCD',
#         [],
#     ),
#     'ak8_w__partnet': Tagger(
#         'partnet_WvsQCD',
#         ['msd > 65', 'msd < 105'],
#     ),
#     'ak8_w__deepak8': Tagger(
#         'deepak8_WvsQCD',
#         ['msd > 65', 'msd < 105'],
#     ),
# }

_TAGGERS = [
    Tagger('hotvr_t_incl__tau',
        'tau32',
        [],
    ),
    Tagger('hotvr_t__tau',
        'tau32',
        ['mass > 140', 'mass < 220', 'fpt1 < 0.8', 'nsub > 2', 'mpair > 50'],
    ),
    Tagger('ak8_t_incl__tau',
        'tau32',
        [],
    ),
    Tagger('ak8_t__tau',
        'tau32',
        ['msd > 105', 'msd < 210'],
    ),
    # Tagger('ak8_t_btagDJet__tau',
    #     'tau32',
    #     ['msd > 105', 'msd < 210', 'subdeepjet > '+str(_DEEPJET_WPS[year]['loose'])],
    # ),
    # Tagger('ak8_t_btagDCSV__tau',
    #     'tau32',
    #     ['msd > 105', 'msd < 210', 'subdeepcsv > '+str(_DEEPCSV_WPS[year]['loose'])],
    # ),
    Tagger('ak8_t__partnet',
        'partnet_TvsQCD',
        ['msd > 105', 'msd < 210'],
    ),
    Tagger('ak8_t__deepak8',
        'deepak8_TvsQCD',
        ['msd > 105', 'msd < 210'],
    ),
    Tagger('ak8_w_incl__tau',
        'tau21',
        [],
    ),
    Tagger('ak8_w__tau',
        'tau21',
        ['msd > 65', 'msd < 105'],
    ),
    Tagger('ak8_w_incl__partnet',
        'partnet_WvsQCD',
        [],
    ),
    Tagger('ak8_w__partnet',
        'partnet_WvsQCD',
        ['msd > 65', 'msd < 105'],
    ),
    Tagger('ak8_w__deepak8',
        'deepak8_WvsQCD',
        ['msd > 65', 'msd < 105'],
    ),
]
_TAGGERS = {tagger.name: tagger for tagger in _TAGGERS}

class Band():

    def __init__(self, name, index):
        self.name = name
        self.index = index

_BANDS = [
    Band('Main', 0),
    Band('QCD', 1),
]
_BANDS = {band.name: band for band in _BANDS}

#kWhite  = 0,   kBlack  = 1,   kGray    = 920,  kRed    = 632,  kGreen  = 416,
#kBlue   = 600, kYellow = 400, kMagenta = 616,  kCyan   = 432,  kOrange = 800,
#kSpring = 820, kTeal   = 840, kAzure   =  860, kViolet = 880,  kPink   = 900
_TCOLORS = {
    # Colors close to matplotlib's default color cycle:
    'pyplot_blue': 800 + 7,
    'pyplot_orange': 860 + 9,
    'pyplot_green': 416 + 2,
    'pyplot_red': 632 + 1,
}
