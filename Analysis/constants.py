from __future__ import print_function

import sys
from collections import OrderedDict
import numpy as np
import sympy
import os
import pandas as pd
import glob
import json


# Source: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopSystematics#Luminosity
# See also: https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM#SummaryTable (includes links to CMS publications which you can cite)
# For correlations, see: https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun2#Combination_and_correlations
_YEARS = OrderedDict([
    ('UL16', {
        'short_name': 'UL16',
        'year': '2016',
        'long_name': 'Ultra Legacy 2016',
        'lumi_fb': 36.333, # TopSystematics source and my own brilcalc give the same results (36.333380073.916976929 with PHYSICS normtag)
        'lumi_pb': 36333.,
        'lumi_unc': 0.012,
        'lumi_unc_detailed': {
            'Uncorrelated16': 0.01,
            'Correlated161718': 0.006,
        },
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
        'lumi_unc_detailed': {
            'Uncorrelated16': 0.01,
            'Correlated161718': 0.006,
        },
        'lumi_fb_display': '19.5',
        'index': 1, # same as in LegacyTopTagging/include/Constants.h
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
        'lumi_unc_detailed': {
            'Uncorrelated16': 0.01,
            'Correlated161718': 0.006,
        },
        'lumi_fb_display': '16.8',
        'index': 2, # same as in LegacyTopTagging/include/Constants.h
    }),
    ('UL17', {
        'short_name': 'UL17',
        'year': '2017',
        'long_name': 'Ultra Legacy 2017',
        'lumi_fb': 41.480, # TopSystematics source and my own brilcalc give the same results (41.4796805287616 with PHYSICS normtag)
        'lumi_pb': 41480.,
        'lumi_unc': 0.023,
        'lumi_unc_detailed': {
            'Uncorrelated17': 0.02,
            'Correlated161718': 0.009,
            'Correlated1718': 0.006,
        },
        'lumi_fb_display': '41.5',
        'index': 3, # same as in LegacyTopTagging/include/Constants.h
    }),
    ('UL18', {
        'short_name': 'UL18',
        'year': '2018',
        'long_name': 'Ultra Legacy 2018',
        'lumi_fb': 59.832, # TopSystematics source and my own brilcalc give the same results (59.8324753390886 with PHYSICS normtag)
        'lumi_pb': 59832.,
        'lumi_unc': 0.025,
        'lumi_unc_detailed': {
            'Uncorrelated18': 0.015,
            'Correlated161718': 0.02,
            'Correlated1718': 0.002,
        },
        'lumi_fb_display': '59.8',
        'index': 4, # same as in LegacyTopTagging/include/Constants.h
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

    def __init__(self, var, var_min=None, var_max=None, inversed=False, fit=False, total_range=False):
        self.var = var
        self.var_min = var_min or 0
        self.var_max = var_max or 999999
        self.max_set = True if var_max else False
        self.name = '_'.join([var, str(var_min)+'to'+str(var_max or 'Inf')])
        self.inversed = inversed # defines in which direction to scan; e.g. either 1 to 0 or 0 to 1
        self.fit = fit
        self.total_range = total_range

    def get_interval_tex(self):
        return '{}--{}'.format(self.var_min, self.var_max if self.max_set else r'$\infty$')

_PT_INTERVALS = [
    VarInterval('pt', 200),
    VarInterval('pt', 200, 250),
    VarInterval('pt', 200, 300),
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
    VarInterval('pt', 500, 600), # for DeepAK8
    VarInterval('pt', 600),
    VarInterval('pt', 600, 800),
    VarInterval('pt', 800),
    VarInterval('pt', 800, 1000),
    VarInterval('pt', 1000),
    VarInterval('pt', 1000, 1500),
    VarInterval('pt', 1500),
    VarInterval('pt', 600, 620), # just for quick tests
]
_PT_INTERVALS = {pt_interval.name: pt_interval for pt_interval in _PT_INTERVALS}

_PT_INTERVALS_TANDP_HOTVR = [
    VarInterval('pt', 200,
        # fit=True,
        total_range=True,
    ),
    VarInterval('pt', 200, 250,
        fit=True,
    ),
    VarInterval('pt', 250, 300,
        fit=True,
    ),
    VarInterval('pt', 300, 350,
        fit=True,
    ),
    VarInterval('pt', 350, 400,
        fit=True,
    ),
    VarInterval('pt', 300, 400,
        # fit=True,
    ),
    VarInterval('pt', 400,
        # fit=True,
    ),
    VarInterval('pt', 400, 480,
        fit=True,
    ),
    VarInterval('pt', 480, 600,
        fit=True,
    ),
    VarInterval('pt', 600,
        fit=True,
    ),
]
_PT_INTERVALS_TANDP_HOTVR = {pt_interval.name: pt_interval for pt_interval in _PT_INTERVALS_TANDP_HOTVR}

_PT_INTERVALS_TANDP_AK8_T = [
    VarInterval('pt', 300,
        # fit=True,
        total_range=True,
    ),
    VarInterval('pt', 300, 350,
        # fit=True,
    ),
    VarInterval('pt', 350, 400,
        # fit=True,
    ),
    VarInterval('pt', 300, 400,
        fit=True,
    ),
    VarInterval('pt', 400,
        # fit=True,
        # total_range=True,
    ),
    VarInterval('pt', 400, 480,
        fit=True,
    ),
    VarInterval('pt', 480, 600,
        fit=True,
    ),
    VarInterval('pt', 600,
        fit=True,
    ),
]
_PT_INTERVALS_TANDP_AK8_T = {pt_interval.name: pt_interval for pt_interval in _PT_INTERVALS_TANDP_AK8_T}

_PT_INTERVALS_TANDP_AK8_W = [
    VarInterval('pt', 200,
        # fit=True,
        total_range=True,
    ),
    VarInterval('pt', 200, 250,
        fit=True,
    ),
    # VarInterval('pt', 200, 300,
    #     fit=True,
    # ),
    VarInterval('pt', 250, 300,
        fit=True,
    ),
    VarInterval('pt', 300,
        # fit=True,
    ),
    VarInterval('pt', 300, 350,
        fit=True,
    ),
    VarInterval('pt', 350, 400,
        fit=True,
    ),
    VarInterval('pt', 300, 400,
        # fit=True,
    ),
    VarInterval('pt', 400,
        fit=True,
    ),
    VarInterval('pt', 400, 480,
        # fit=True,
    ),
    VarInterval('pt', 480, 600,
        # fit=True,
    ),
    VarInterval('pt', 600,
        # fit=True,
    ),
]
_PT_INTERVALS_TANDP_AK8_W = {pt_interval.name: pt_interval for pt_interval in _PT_INTERVALS_TANDP_AK8_W}


class WorkingPoint:

    def __init__(self, bkg_eff, cut_value, name=None, null=False, alt_name=None, alt_name_short=None):
        self.bkg_eff = bkg_eff
        self.cut_value = cut_value # can be float or dict (key=year, value=float)
        self.name = name or 'BkgEff' + '{:.3f}'.format(bkg_eff).replace('.', 'p')
        self.null = null
        self.alt_name = alt_name
        self.alt_name_short = alt_name_short

    def get_cut_value(self, year=None):
        '''
        Use this function, not self.cut_value directly!
        '''
        if isinstance(self.cut_value, dict):
            if year is not None:
                return self.cut_value[year]
            else:
                sys.exit('Cut values of working point {} are year-dependent. But year not provided'.format(self.name))
        else:
            return self.cut_value

_NULL_WP = WorkingPoint(1., np.inf, 'NullWP', null=True)


class Tagger:

    def __init__(self, name, scan_var, tag_rule = [], wps = [], tandp_rule='', label='', use_all_pt_intervals=False):
        self.name = name
        self.scan_var = scan_var
        self.tag_rule = tag_rule
        if self.name.startswith('ak8_t'):
            self.var_intervals = _PT_INTERVALS_TANDP_AK8_T
            self.fit_variable = 'output_probejet_AK8_mSD'
        elif self.name.startswith('ak8_w'):
            self.var_intervals = _PT_INTERVALS_TANDP_AK8_W
            self.fit_variable = 'output_probejet_AK8_mSD'
        elif self.name.startswith('hotvr_t'):
            self.var_intervals = _PT_INTERVALS_TANDP_HOTVR
            self.fit_variable = 'output_probejet_HOTVR_mass'
        else:
            self.var_intervals = _PT_INTERVALS # just use all for now
            self.fit_variable = 'VAR_NOT_DEFINED' # fixme
        if use_all_pt_intervals:
            self.var_intervals = _PT_INTERVALS
        self.wps = wps # cut values for scan_var; should be list of instances of class "WorkingPoint"
        self.tandp_rule = tandp_rule
        self.label = label

    def set_tag_rule(self, tag_rule):
        self.tag_rule = tag_rule

    def get_tag_rule(self, year=None):
        if len(self.tag_rule):
            tag_rule = ' & '.join(['('+x+')' for x in self.tag_rule])
        else:
            tag_rule = 'True'
        replaces = {}
        if year:
            replaces['DEEPCSV_SCORE'] = '{:.5f}'.format(_DEEPCSV_WPS[year]['loose'])
            replaces['DEEPJET_SCORE'] = '{:.5f}'.format(_DEEPJET_WPS[year]['loose'])
        return tag_rule.format(**replaces)

    def get_wp(self, wp_index=None, year=None, wp_name=None):
        # if wp_index not given, will return list of all wps for the given year
        # if year not given, will return
        if isinstance(self.wps, dict):
            if not year:
                sys.exit('Working points for tagger are year-dependent, but you have not specified a year')
            elif year not in self.wps:
                sys.exit('Chosen year "{year}" has no WP defined'.format(year=year))
            if wp_index != None:
                return self.wps[year][wp_index]
            else:
                return self.wps[year]
        elif isinstance(self.wps, list):
            if wp_index != None:
                return self.wps[wp_index]
            elif wp_name != None:
                result = None
                for wp in self.wps:
                    if wp_name == wp.name:
                        result = wp
                        break
                if result is not None:
                    return result
                else:
                    sys.exit('WP with name {} not found'.format(wp_name))
            else:
                return self.wps
        else:
            sys.exit('Format of self.wps not correct. Needs to be either dict of lists (keys=years) or plain list if no year dependency')

    def get_tandp_rule(self, wp_index, year=None):
        replaces = {}
        # replaces['WP_VALUE'] = self.get_wp(wp_index, year).cut_value
        replaces['WP_VALUE'] = self.get_wp(wp_index, year).get_cut_value(year)
        if year:
            replaces['DEEPCSV_SCORE'] = '{:.5f}'.format(_DEEPCSV_WPS[year]['loose'])
            replaces['DEEPJET_SCORE'] = '{:.5f}'.format(_DEEPJET_WPS[year]['loose'])
        return self.tandp_rule.format(**replaces)

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
        tandp_rule='(output_probejet_HOTVR_tau32 < {WP_VALUE})',
        wps=[
            WorkingPoint(1., 0., name='Inclusive')
        ],
        label='#tau_{3}/#tau_{2} < {WP_VALUE} + HOTVR cuts',
    ),
    Tagger('hotvr_t__tau',
        'tau32',
        ['mass > 140', 'mass < 220', 'fpt1 < 0.8', 'nsub > 2', 'mpair > 50'],
        tandp_rule='(output_probejet_HOTVR_fpt1 < 0.8) & (output_probejet_HOTVR_nsub > 2) & (output_probejet_HOTVR_mpair > 50) & (output_probejet_HOTVR_tau32 < {WP_VALUE})',
        wps=[
        # Varified on Nov 15, 2022 with my WorkingPointStudy code that for UL with pt > 400 GeV, |eta| < 2.5, a tau32 cut of < 0.550 or 0.551 (depending on year) would lead to 3% bkg. eff. in QCD HT, HT > 300 GeV
        # We stick to 0.56 for consistency with previous
            WorkingPoint(9.999, 0.56, name='Standard') # providing bkg_eff does not make sense here; I have not derived it on my own; should be 3% or so
        ],
        label='#tau_{3}/#tau_{2} < {WP_VALUE} + HOTVR cuts',
    ),
    Tagger('ak8_t_incl__tau',
        'tau32',
        [],
        tandp_rule='(output_probejet_AK8_tau32 < {WP_VALUE})',
        wps=[
            WorkingPoint(1., 0., name='Inclusive')
        ],
        label='#tau_{3}/#tau_{2} < {WP_VALUE}',
    ),
    Tagger('ak8_t__tau',
        'tau32',
        ['msd > 105', 'msd < 210'],
        tandp_rule='(output_probejet_AK8_tau32 < {WP_VALUE})',
        wps=[
            WorkingPoint(0.001, 0.38, alt_name='very tight', alt_name_short='vt'),
            WorkingPoint(0.005, 0.47, alt_name='tight', alt_name_short='t'),
            WorkingPoint(0.010, 0.52, alt_name='medium', alt_name_short='m'),
            WorkingPoint(0.025, 0.61, alt_name='loose', alt_name_short='l'),
            WorkingPoint(0.050, 0.69, alt_name='very loose', alt_name_short='vl'),
        ],
        label='#tau_{3}/#tau_{2} < {WP_VALUE}',
    ),
    Tagger('ak8_t_btagDJet__tau',
        'tau32',
        # ['msd > 105', 'msd < 210', 'subdeepjet > '+str(_DEEPJET_WPS[year]['loose'])],
        ['msd > 105', 'msd < 210', 'subdeepjet > {DEEPJET_SCORE}'],
        tandp_rule='(output_probejet_AK8_maxDeepJet > {DEEPJET_SCORE}) & (output_probejet_AK8_tau32 < {WP_VALUE})',
        wps=[
            WorkingPoint(0.001, 0.38, alt_name='very tight', alt_name_short='vt'),
            WorkingPoint(0.005, 0.47, alt_name='tight', alt_name_short='t'),
            WorkingPoint(0.010, 0.52, alt_name='medium', alt_name_short='m'),
            WorkingPoint(0.025, 0.61, alt_name='loose', alt_name_short='l'),
            WorkingPoint(0.050, 0.69, alt_name='very loose', alt_name_short='vl'),
        ],
        label='#tau_{3}/#tau_{2} < {WP_VALUE} + loose DeepJet subjet b tag',
    ),
    Tagger('ak8_t_btagDCSV__tau',
        'tau32',
        # ['msd > 105', 'msd < 210', 'subdeepcsv > '+str(_DEEPCSV_WPS[year]['loose'])],
        ['msd > 105', 'msd < 210', 'subdeepcsv > {DEEPCSV_SCORE}'],
        tandp_rule='(output_probejet_AK8_maxDeepCSV > {DEEPCSV_SCORE}) & (output_probejet_AK8_tau32 < {WP_VALUE})',
        wps=[
            WorkingPoint(0.001, 0.38, alt_name='very tight', alt_name_short='vt'),
            WorkingPoint(0.005, 0.47, alt_name='tight', alt_name_short='t'),
            WorkingPoint(0.010, 0.52, alt_name='medium', alt_name_short='m'),
            WorkingPoint(0.025, 0.61, alt_name='loose', alt_name_short='l'),
            WorkingPoint(0.050, 0.69, alt_name='very loose', alt_name_short='vl'),
        ],
        label='#tau_{3}/#tau_{2} < {WP_VALUE} + loose DeepCSV subjet b tag',
    ),
    Tagger('ak8_t_incl__partnet',
        'partnet_TvsQCD',
        [],
        tandp_rule='(output_probejet_AK8_ParticleNet_TvsQCD > {WP_VALUE})',
    ),
    Tagger('ak8_t__partnet',
        'partnet_TvsQCD',
        ['msd > 105', 'msd < 210'],
        tandp_rule='(output_probejet_AK8_ParticleNet_TvsQCD > {WP_VALUE})',
    ),
    Tagger('ak8_t__deepak8',
        'deepak8_TvsQCD',
        ['msd > 105', 'msd < 210'],
        tandp_rule='(output_probejet_AK8_DeepAK8_TvsQCD > {WP_VALUE})',
    ),
    Tagger('ak8_t_incl__deepak8',
        'deepak8_TvsQCD',
        [],
        tandp_rule='(output_probejet_AK8_DeepAK8_TvsQCD > {WP_VALUE})',
    ),
    Tagger('ak8_t__MDdeepak8',
        'MDdeepak8_TvsQCD',
        ['msd > 105', 'msd < 210'],
        tandp_rule='(output_probejet_AK8_MDDeepAK8_TvsQCD > {WP_VALUE})',
        wps=[
        # Re-varified on Nov 15, 2022 for pt > 200 (first entry) and pt > 400 (second entry)
            WorkingPoint(0.010,
            {
                'UL16preVFP': 0.485,
                'UL16postVFP': 0.475,
                'UL17': 0.487,
                'UL18': 0.477,
            },
            name='CustomPt400'
            ),
        ],
        # wps={
        # # Re-varified on Nov 15, 2022 for pt > 200 (first entry) and pt > 400 (second entry)
        #     'UL16preVFP': [
        #         WorkingPoint(0.010, 0.485, name='CustomPt400'),
        #     ],
        #     'UL16postVFP': [
        #         WorkingPoint(0.010, 0.475, name='CustomPt400'),
        #     ],
        #     'UL17': [
        #         WorkingPoint(0.010, 0.487, name='CustomPt400'),
        #     ],
        #     'UL18': [
        #         WorkingPoint(0.010, 0.477, name='CustomPt400'),
        #     ],
        # },
        label='MD-DeepAK8 TvsQCD > {WP_VALUE}',
    ),
    Tagger('ak8_t_incl__MDdeepak8',
        'MDdeepak8_TvsQCD',
        [],
        tandp_rule='(output_probejet_AK8_MDDeepAK8_TvsQCD > {WP_VALUE})',
    ),
    Tagger('ak8_w_incl__tau',
        'tau21',
        [],
        tandp_rule='(output_probejet_AK8_tau21 < {WP_VALUE})',
        wps=[
            WorkingPoint(1., 0., name='Inclusive')
        ],
        label='#tau_{2}/#tau_{1} < {WP_VALUE}',
    ),
    Tagger('ak8_w__tau',
        'tau21',
        ['msd > 65', 'msd < 105'],
        tandp_rule='(output_probejet_AK8_tau21 < {WP_VALUE})',
    ),
    Tagger('ak8_w_incl__partnet',
        'partnet_WvsQCD',
        [],
        tandp_rule='(output_probejet_AK8_ParticleNet_WvsQCD > {WP_VALUE})',
    ),
    Tagger('ak8_w__partnet',
        'partnet_WvsQCD',
        ['msd > 65', 'msd < 105'],
        tandp_rule='(output_probejet_AK8_ParticleNet_WvsQCD > {WP_VALUE})',
        wps=[
        # Re-varified on Nov 15, 2022 for pt > 200 (first entry) and pt > 400 (second entry)
            WorkingPoint(0.030,
            {
                'UL16preVFP': 0.737,
                'UL16postVFP': 0.729,
                'UL17': 0.746,
                'UL18': 0.734,
            },
            name='CustomPt400'
            ),
        ],
        # wps={
        # # Re-varified on Nov 15, 2022 for pt > 200 (first entry) and pt > 400 (second entry)
        #     'UL16preVFP': [
        #         # WorkingPoint(0.030, 0.871, name='CustomPt200'),
        #         WorkingPoint(0.030, 0.737, name='CustomPt400'),
        #     ],
        #     'UL16postVFP': [
        #         # WorkingPoint(0.030, 0.868, name='CustomPt200'),
        #         WorkingPoint(0.030, 0.729, name='CustomPt400'),
        #     ],
        #     'UL17': [
        #         # WorkingPoint(0.030, 0.868, name='CustomPt200'),
        #         WorkingPoint(0.030, 0.746, name='CustomPt400'),
        #     ],
        #     'UL18': [
        #         # WorkingPoint(0.030, 0.864, name='CustomPt200'),
        #         WorkingPoint(0.030, 0.734, name='CustomPt400'),
        #     ],
        # },
        label='ParticleNet WvsQCD > {WP_VALUE}',
    ),
    Tagger('ak8_w__deepak8',
        'deepak8_WvsQCD',
        ['msd > 65', 'msd < 105'],
        tandp_rule='(output_probejet_AK8_DeepAK8_WvsQCD > {WP_VALUE})',
    ),
    Tagger('ak8_w_incl__deepak8',
        'deepak8_WvsQCD',
        [],
        tandp_rule='(output_probejet_AK8_DeepAK8_WvsQCD > {WP_VALUE})',
    ),
    Tagger('ak8_w__MDdeepak8',
        'MDdeepak8_WvsQCD',
        ['msd > 65', 'msd < 105'],
        tandp_rule='(output_probejet_AK8_MDDeepAK8_WvsQCD > {WP_VALUE})',
    ),
    Tagger('ak8_w_incl__MDdeepak8',
        'MDdeepak8_WvsQCD',
        [],
        tandp_rule='(output_probejet_AK8_MDDeepAK8_WvsQCD > {WP_VALUE})',
    ),
]
_TAGGERS = {tagger.name: tagger for tagger in _TAGGERS}

class Band():

    def __init__(self, name, index):
        self.name = name
        self.index = index

_BANDS = [
    Band('Main', 0),
    # Band('QCD', 1),
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

# Correlations taken from: https://docs.google.com/spreadsheets/d/1JZfk78_9SD225bcUuTWVo4i02vwI5FfeVKH-dwzUdhM/edit#gid=1345121349
# This spreadsheet is linked here: https://twiki.cern.ch/twiki/bin/view/CMS/JECUncertaintySources#Run2_JEC_uncertainty_correlation
# keys are the source names, compatible with the names given in the JEC textFiles. Values are the correlations between RunII years
# 0.5 are the optimal values; can be simplified to 1.0 if necessary

class JECSmearSource():

    def __init__(self, name, era_correlation, is_split_full, is_split_reduced):
        self.name = name
        self.era_correlation = era_correlation
        self.is_split_full = is_split_full
        self.is_split_reduced = is_split_reduced

    def get_name(year=None):
        if year is None or not ('year' in self.name):
            return self.name
        else:
            return self.name.format(year=year)


class JECSmearSources():

    def __init__(self):
        self.base = [
            JECSmearSource('AbsoluteMPFBias', 1.0, True, False),
            JECSmearSource('AbsoluteScale', 1.0, True, False),
            JECSmearSource('AbsoluteStat', 0.0, True, False),
            JECSmearSource('FlavorQCD', 1.0, True, False),
            JECSmearSource('Fragmentation', 1.0, True, False),
            JECSmearSource('PileUpDataMC', 0.5, True, False),
            JECSmearSource('PileUpPtBB', 0.5, True, False),
            JECSmearSource('PileUpPtEC1', 0.5, True, False),
            JECSmearSource('PileUpPtEC2', 0.5, True, False),
            JECSmearSource('PileUpPtHF', 0.5, True, False),
            JECSmearSource('PileUpPtRef', 0.5, True, False),
            JECSmearSource('RelativeFSR', 0.5, True, False),
            JECSmearSource('RelativeJEREC1', 0.0, True, False),
            JECSmearSource('RelativeJEREC2', 0.0, True, False),
            JECSmearSource('RelativeJERHF', 0.5, True, False),
            JECSmearSource('RelativePtBB', 0.5, True, False),
            JECSmearSource('RelativePtEC1', 0.0, True, False),
            JECSmearSource('RelativePtEC2', 0.0, True, False),
            JECSmearSource('RelativePtHF', 0.5, True, False),
            JECSmearSource('RelativeBal', 0.5, True, False),
            JECSmearSource('RelativeSample', 0.0, True, False),
            JECSmearSource('RelativeStatEC', 0.0, True, False),
            JECSmearSource('RelativeStatFSR', 0.0, True, False),
            JECSmearSource('RelativeStatHF', 0.0, True, False),
            JECSmearSource('SinglePionECAL', 1.0, True, False),
            JECSmearSource('SinglePionHCAL', 1.0, True, False),
            JECSmearSource('TimePtEta', 0.0, True, False),
        ]




_JECSMEAR_SOURCES = {
    'AbsoluteMPFBias': 1.0,
    'AbsoluteScale': 1.0,
    'AbsoluteStat': 0.0,
    'FlavorQCD': 1.0,
    'Fragmentation': 1.0,
    'PileUpDataMC': 0.5,
    'PileUpPtBB': 0.5,
    'PileUpPtEC1': 0.5,
    'PileUpPtEC2': 0.5,
    'PileUpPtHF': 0.5,
    'PileUpPtRef': 0.5,
    'RelativeFSR': 0.5,
    'RelativeJEREC1': 0.0,
    'RelativeJEREC2': 0.0,
    'RelativeJERHF': 0.5,
    'RelativePtBB': 0.5,
    'RelativePtEC1': 0.0,
    'RelativePtEC2': 0.0,
    'RelativePtHF': 0.5,
    'RelativeBal': 0.5,
    'RelativeSample': 0.0,
    'RelativeStatEC': 0.0,
    'RelativeStatFSR': 0.0,
    'RelativeStatHF': 0.0,
    'SinglePionECAL': 1.0,
    'SinglePionHCAL': 1.0,
    'TimePtEta': 0.0,
}

# class Systematic:
#
#     def __init__(self, name, weight_alias=None):
#         self.name = name
#         self.weight_alias = weight_alias or 'weight'
#         self.weight_based = True if weight_alias else False
#         self.variations = [] # only used when compiling compressed list of systematics
#
# _SYSTEMATICS = [
#     Systematic('nominal',
#         None,
#     ),
#
#     # b tagging
#     Systematic('btag_bc_down',
#         'weight / weight_btag_central * weight_btag_bc_down',
#     ),
#     Systematic('btag_bc_correlated_down',
#         'weight / weight_btag_central * weight_btag_bc_down_correlated',
#     ),
#     Systematic('btag_bc_jes_down',
#         'weight / weight_btag_central * weight_btag_bc_down_jes',
#     ),
#     Systematic('btag_bc_pileup_down',
#         'weight / weight_btag_central * weight_btag_bc_down_pileup',
#     ),
#     Systematic('btag_bc_statistic_down',
#         'weight / weight_btag_central * weight_btag_bc_down_statistic',
#     ),
#     Systematic('btag_bc_type3_down',
#         'weight / weight_btag_central * weight_btag_bc_down_type3',
#     ),
#     Systematic('btag_bc_uncorrelated_down',
#         'weight / weight_btag_central * weight_btag_bc_down_uncorrelated',
#     ),
#
#     Systematic('btag_bc_up',
#         'weight / weight_btag_central * weight_btag_bc_up',
#     ),
#     Systematic('btag_bc_correlated_up',
#         'weight / weight_btag_central * weight_btag_bc_up_correlated',
#     ),
#     Systematic('btag_bc_jes_up',
#         'weight / weight_btag_central * weight_btag_bc_up_jes',
#     ),
#     Systematic('btag_bc_pileup_up',
#         'weight / weight_btag_central * weight_btag_bc_up_pileup',
#     ),
#     Systematic('btag_bc_statistic_up',
#         'weight / weight_btag_central * weight_btag_bc_up_statistic',
#     ),
#     Systematic('btag_bc_type3_up',
#         'weight / weight_btag_central * weight_btag_bc_up_type3',
#     ),
#     Systematic('btag_bc_uncorrelated_up',
#         'weight / weight_btag_central * weight_btag_bc_up_uncorrelated',
#     ),
#
#     Systematic('btag_light_down',
#         'weight / weight_btag_central * weight_btag_light_down',
#     ),
#     Systematic('btag_light_correlated_down',
#         'weight / weight_btag_central * weight_btag_light_down_correlated',
#     ),
#     Systematic('btag_light_uncorrelated_down',
#         'weight / weight_btag_central * weight_btag_light_down_uncorrelated',
#     ),
#
#     Systematic('btag_light_up',
#         'weight / weight_btag_central * weight_btag_light_up',
#     ),
#     Systematic('btag_light_correlated_up',
#         'weight / weight_btag_central * weight_btag_light_up_correlated',
#     ),
#     Systematic('btag_light_uncorrelated_up',
#         'weight / weight_btag_central * weight_btag_light_up_uncorrelated',
#     ),
#
#     # Scale
#     Systematic('murmuf_downdown',
#         'weight * weight_murmuf_downdown',
#     ),
#     Systematic('murmuf_downnone',
#         'weight * weight_murmuf_downnone',
#     ),
#     Systematic('murmuf_nonedown',
#         'weight * weight_murmuf_nonedown',
#     ),
#     Systematic('murmuf_noneup',
#         'weight * weight_murmuf_noneup',
#     ),
#     Systematic('murmuf_upnone',
#         'weight * weight_murmuf_upnone',
#     ),
#     Systematic('murmuf_upup',
#         'weight * weight_murmuf_upup',
#     ),
#
#     # FSR inclusive
#     Systematic('fsr_2_down',
#         'weight * weight_fsr_2_down',
#     ),
#     Systematic('fsr_2_up',
#         'weight * weight_fsr_2_up',
#     ),
#     # FSR splitted
#     Systematic('fsr_g2gg_mur_down',
#         'weight * weight_fsr_g2gg_mur_down',
#     ),
#     Systematic('fsr_g2gg_mur_up',
#         'weight * weight_fsr_g2gg_mur_up',
#     ),
#     Systematic('fsr_g2qq_mur_down',
#         'weight * weight_fsr_g2qq_mur_down',
#     ),
#     Systematic('fsr_g2qq_mur_up',
#         'weight * weight_fsr_g2qq_mur_up',
#     ),
#     Systematic('fsr_q2qg_mur_down',
#         'weight * weight_fsr_q2qg_mur_down',
#     ),
#     Systematic('fsr_q2qg_mur_up',
#         'weight * weight_fsr_q2qg_mur_up',
#     ),
#     Systematic('fsr_x2xg_mur_down',
#         'weight * weight_fsr_x2xg_mur_down',
#     ),
#     Systematic('fsr_x2xg_mur_up',
#         'weight * weight_fsr_x2xg_mur_up',
#     ),
#     Systematic('fsr_g2gg_cns_down',
#         'weight * weight_fsr_g2gg_cns_down',
#     ),
#     Systematic('fsr_g2gg_cns_up',
#         'weight * weight_fsr_g2gg_cns_up',
#     ),
#     Systematic('fsr_g2qq_cns_down',
#         'weight * weight_fsr_g2qq_cns_down',
#     ),
#     Systematic('fsr_g2qq_cns_up',
#         'weight * weight_fsr_g2qq_cns_up',
#     ),
#     Systematic('fsr_q2qg_cns_down',
#         'weight * weight_fsr_q2qg_cns_down',
#     ),
#     Systematic('fsr_q2qg_cns_up',
#         'weight * weight_fsr_q2qg_cns_up',
#     ),
#     Systematic('fsr_x2xg_cns_down',
#         'weight * weight_fsr_x2xg_cns_down',
#     ),
#     Systematic('fsr_x2xg_cns_up',
#         'weight * weight_fsr_x2xg_cns_up',
#     ),
#
#     # ISR inclusive
#     Systematic('isr_2_down',
#         'weight * weight_isr_2_down',
#     ),
#     Systematic('isr_2_up',
#         'weight * weight_isr_2_up',
#     ),
#     # ISR splitted
#     Systematic('isr_g2gg_mur_down',
#         'weight * weight_isr_g2gg_mur_down',
#     ),
#     Systematic('isr_g2gg_mur_up',
#         'weight * weight_isr_g2gg_mur_up',
#     ),
#     Systematic('isr_g2qq_mur_down',
#         'weight * weight_isr_g2qq_mur_down',
#     ),
#     Systematic('isr_g2qq_mur_up',
#         'weight * weight_isr_g2qq_mur_up',
#     ),
#     Systematic('isr_q2qg_mur_down',
#         'weight * weight_isr_q2qg_mur_down',
#     ),
#     Systematic('isr_q2qg_mur_up',
#         'weight * weight_isr_q2qg_mur_up',
#     ),
#     Systematic('isr_x2xg_mur_down',
#         'weight * weight_isr_x2xg_mur_down',
#     ),
#     Systematic('isr_x2xg_mur_up',
#         'weight * weight_isr_x2xg_mur_up',
#     ),
#     Systematic('isr_g2gg_cns_down',
#         'weight * weight_isr_g2gg_cns_down',
#     ),
#     Systematic('isr_g2gg_cns_up',
#         'weight * weight_isr_g2gg_cns_up',
#     ),
#     Systematic('isr_g2qq_cns_down',
#         'weight * weight_isr_g2qq_cns_down',
#     ),
#     Systematic('isr_g2qq_cns_up',
#         'weight * weight_isr_g2qq_cns_up',
#     ),
#     Systematic('isr_q2qg_cns_down',
#         'weight * weight_isr_q2qg_cns_down',
#     ),
#     Systematic('isr_q2qg_cns_up',
#         'weight * weight_isr_q2qg_cns_up',
#     ),
#     Systematic('isr_x2xg_cns_down',
#         'weight * weight_isr_x2xg_cns_down',
#     ),
#     Systematic('isr_x2xg_cns_up',
#         'weight * weight_isr_x2xg_cns_up',
#     ),
#
#     # Prefiring
#     Systematic('prefire_down',
#         'weight / weight_prefire * weight_prefire_down',
#     ),
#     Systematic('prefire_up',
#         'weight / weight_prefire * weight_prefire_up',
#     ),
#
#     # Pileup
#     Systematic('pu_down',
#         'weight / weight_pu * weight_pu_down',
#     ),
#     Systematic('pu_up',
#         'weight / weight_pu * weight_pu_up',
#     ),
#
#     # Systematic('sfelec_id_down',
#     #     'weight / weight_sfelec_id * weight_sfelec_id_down',
#     # ),
#     # Systematic('sfelec_id_up',
#     #     'weight / weight_sfelec_id * weight_sfelec_id_up',
#     # ),
#     #
#     # Systematic('sfelec_reco_down',
#     #     'weight / weight_sfelec_reco * weight_sfelec_reco_down',
#     # ),
#     # Systematic('sfelec_reco_up',
#     #     'weight / weight_sfelec_reco * weight_sfelec_reco_up',
#     # ),
#
#     # Muon ID
#     Systematic('sfmu_id_down',
#         'weight / weight_sfmu_id * weight_sfmu_id_down',
#     ),
#     Systematic('sfmu_id_up',
#         'weight / weight_sfmu_id * weight_sfmu_id_up',
#     ),
#
#     # Muon Trigger
#     Systematic('sfmu_trigger_down',
#         'weight / weight_sfmu_trigger * weight_sfmu_trigger_down',
#     ),
#     Systematic('sfmu_trigger_up',
#         'weight / weight_sfmu_trigger * weight_sfmu_trigger_up',
#     ),
#
#     # TopPt Parameter A
#     Systematic('toppt_a_down',
#         'weight / weight_toppt_applied * weight_toppt_a_down',
#     ),
#     Systematic('toppt_a_up',
#         'weight / weight_toppt_applied * weight_toppt_a_up',
#     ),
#
#     # TopPt Parameter B
#     Systematic('toppt_b_down',
#         'weight / weight_toppt_applied * weight_toppt_b_down',
#     ),
#     Systematic('toppt_b_up',
#         'weight / weight_toppt_applied * weight_toppt_b_up',
#     ),
#
#     # JES
#     Systematic('jesTotal_down',
#         None,
#     ),
#     Systematic('jesTotal_up',
#         None,
#     ),
#
#     # JER
#     Systematic('jer_down',
#         None,
#     ),
#     Systematic('jer_up',
#         None,
#     ),
#
#     # Top mass
#     Systematic('mtop_mtop171p5',
#         None,
#     ),
#     Systematic('mtop_mtop173p5',
#         None,
#     ),
#
#     # Event tune
#     Systematic('tune_down',
#         None,
#     ),
#     Systematic('tune_up',
#         None,
#     ),
#
#     # ME-PS matching
#     Systematic('hdamp_down',
#         None,
#     ),
#     Systematic('hdamp_up',
#         None,
#     ),
#
#     # Color reconnection
#     Systematic('cr_cr1',
#         None,
#     ),
#     Systematic('cr_cr2',
#         None,
#     ),
#     Systematic('cr_erdon',
#         None,
#     ),
# ]
# for k in _JECSMEAR_SOURCES.keys():
#     _SYSTEMATICS.append(
#         Systematic('jes'+k+'_down',
#         None,
#         )
#     )
#     _SYSTEMATICS.append(
#         Systematic('jes'+k+'_up',
#         None,
#         )
#     )
# _SYSTEMATICS = {syst.name: syst for syst in _SYSTEMATICS}

class Variation:

    def __init__(self, short_name, weight_alias=None):
        self.short_name = short_name
        self.name = None # holds full name inlcuding name of systematic
        self.weight_alias = weight_alias or 'weight'
        self.weight_based = True if weight_alias else False

class Systematic:

    def __init__(self, name, variations, correlation_eras=None, correlation_procs=None, combine_name=None, tandp=False, symmetrize=False):
        self.name = name
        self.variations = {k: Variation(k, v) for k, v in variations.items()}
        for k in self.variations.keys():
            self.variations[k].name = self.name+'_'+self.variations[k].short_name
        self.correlation_eras = 1.0 if correlation_eras is None else correlation_eras # correlation between years (1.0 = 100%)
        self.correlation_procs = 'shared' if correlation_procs is None else correlation_procs # options: None/'shared' (uncertainty is shared between all processes), 'qcdewk' (differentiate between QCD-induced and EWK-induced processes), 'allprocs' (differentiate between all base processes)
        self.combine_name = combine_name or self.name # in combine, systematics cannot have underscores in their names; in case self.name has an underscore, we need to provide a new name here to be used in the combine framework
        self.tandp = tandp # decide whether to use this systematic in the final combine fits (and for prefit plotting)
        self.symmetrize = symmetrize # whether or not to symmetrize the uncertainty # NOT implemented yet #FIXME

class Systematics:

    def __init__(self, include_jes_splits=True, blacklist=[], whitelist=[]):
        self.include_jes_splits = include_jes_splits
        self.nominal = Variation('nominal')
        self.nominal.name = self.nominal.short_name
        self.base = [
            #________________________________________
            Systematic('btag_bc', {
                'down': 'weight / weight_btag_central * weight_btag_bc_down',
                'up': 'weight / weight_btag_central * weight_btag_bc_up',
                },
                correlation_eras=0.0,
                combine_name='btagBC',
            ),
            Systematic('btag_bc_correlated', {
                'down': 'weight / weight_btag_central * weight_btag_bc_down_correlated',
                'up': 'weight / weight_btag_central * weight_btag_bc_up_correlated',
                },
                correlation_eras=1.0,
                correlation_procs=None,
                combine_name='btagBCCorrelated',
                # tandp=True,
            ),
            Systematic('btag_bc_jes', {
                'down': 'weight / weight_btag_central * weight_btag_bc_down_jes',
                'up': 'weight / weight_btag_central * weight_btag_bc_up_jes',
                },
                correlation_eras=1.0,
                correlation_procs=None,
                combine_name='btagBCJES',
                tandp=True,
            ),
            Systematic('btag_bc_pileup', {
                'down': 'weight / weight_btag_central * weight_btag_bc_down_pileup',
                'up': 'weight / weight_btag_central * weight_btag_bc_up_pileup',
                },
                correlation_eras=1.0,
                correlation_procs=None,
                combine_name='btagBCPileup',
                tandp=True,
            ),
            Systematic('btag_bc_statistic', {
                'down': 'weight / weight_btag_central * weight_btag_bc_down_statistic',
                'up': 'weight / weight_btag_central * weight_btag_bc_up_statistic',
                },
                correlation_eras=0.0,
                correlation_procs=None,
                combine_name='btagBCStatistic',
                tandp=True,
            ),
            Systematic('btag_bc_type3', {
                'down': 'weight / weight_btag_central * weight_btag_bc_down_type3',
                'up': 'weight / weight_btag_central * weight_btag_bc_up_type3',
                },
                correlation_eras=1.0,
                correlation_procs=None,
                combine_name='btagBCType3',
                tandp=True,
            ),
            Systematic('btag_bc_uncorrelated', {
                'down': 'weight / weight_btag_central * weight_btag_bc_down_uncorrelated',
                'up': 'weight / weight_btag_central * weight_btag_bc_up_uncorrelated',
                },
                correlation_eras=0.0,
                correlation_procs=None,
                combine_name='btagBCUncorrelated',
                # tandp=True,
            ),
            #________________________________________
            Systematic('btag_light', {
                'down': 'weight / weight_btag_central * weight_btag_light_down',
                'up': 'weight / weight_btag_central * weight_btag_light_up',
                },
                correlation_eras=0.0,
                combine_name='btagLight',
            ),
            Systematic('btag_light_correlated', {
                'down': 'weight / weight_btag_central * weight_btag_light_down_correlated',
                'up': 'weight / weight_btag_central * weight_btag_light_up_correlated',
                },
                correlation_eras=1.0,
                correlation_procs=None,
                combine_name='btagLightCorrelated',
                tandp=True,
            ),
            Systematic('btag_light_uncorrelated', {
                'down': 'weight / weight_btag_central * weight_btag_light_down_uncorrelated',
                'up': 'weight / weight_btag_central * weight_btag_light_up_uncorrelated',
                },
                correlation_eras=0.0,
                correlation_procs=None,
                combine_name='btagLightUncorrelated',
                tandp=True,
            ),
            #________________________________________
            Systematic('murmuf', {
                'downdown': 'weight * weight_murmuf_downdown',
                'downnone': 'weight * weight_murmuf_downnone',
                'nonedown': 'weight * weight_murmuf_nonedown',
                'noneup': 'weight * weight_murmuf_noneup',
                'upnone': 'weight * weight_murmuf_upnone',
                'upup': 'weight * weight_murmuf_upup',
                },
                correlation_eras=1.0,
                # correlation_procs='allprocs',
                correlation_procs='split',
                combine_name='Scale',
                tandp=True,
            ),
            #________________________________________
            Systematic('fsr_2', {
                'down': 'weight * weight_fsr_2_down',
                'up': 'weight * weight_fsr_2_up',
                },
                combine_name='FSR2',
                correlation_eras=1.0,
                # correlation_procs='qcdewk',
                tandp=True,
            ),
            Systematic('fsr_g2gg_mur', {
                'down': 'weight * weight_fsr_g2gg_mur_down',
                'up': 'weight * weight_fsr_g2gg_mur_up',
                },
                combine_name='FSRg2ggMuR',
                correlation_eras=1.0,
                # tandp=True,
            ),
            Systematic('fsr_g2qq_mur', {
                'down': 'weight * weight_fsr_g2qq_mur_down',
                'up': 'weight * weight_fsr_g2qq_mur_up',
                },
                combine_name='FSRg2qqMuR',
                correlation_eras=1.0,
                # tandp=True,
            ),
            Systematic('fsr_q2qg_mur', {
                'down': 'weight * weight_fsr_q2qg_mur_down',
                'up': 'weight * weight_fsr_q2qg_mur_up',
                },
                combine_name='FSRq2qgMuR',
                correlation_eras=1.0,
                # tandp=True,
            ),
            Systematic('fsr_x2xg_mur', {
                'down': 'weight * weight_fsr_x2xg_mur_down',
                'up': 'weight * weight_fsr_x2xg_mur_up',
                },
                combine_name='FSRx2xgMuR',
                correlation_eras=1.0,
                # tandp=True,
            ),
            Systematic('fsr_g2gg_cns', {
                'down': 'weight * weight_fsr_g2gg_cns_down',
                'up': 'weight * weight_fsr_g2gg_cns_up',
                },
                combine_name='FSRg2ggCNS',
                correlation_eras=1.0,
                # tandp=True,
            ),
            Systematic('fsr_g2qq_cns', {
                'down': 'weight * weight_fsr_g2qq_cns_down',
                'up': 'weight * weight_fsr_g2qq_cns_up',
                },
                combine_name='FSRg2qqCNS',
                correlation_eras=1.0,
                # tandp=True,
            ),
            Systematic('fsr_q2qg_cns', {
                'down': 'weight * weight_fsr_q2qg_cns_down',
                'up': 'weight * weight_fsr_q2qg_cns_up',
                },
                combine_name='FSRq2qgCNS',
                correlation_eras=1.0,
                # tandp=True,
            ),
            Systematic('fsr_x2xg_cns', {
                'down': 'weight * weight_fsr_x2xg_cns_down',
                'up': 'weight * weight_fsr_x2xg_cns_up',
                },
                combine_name='FSRx2xgCNS',
                correlation_eras=1.0,
                # tandp=True,
            ),
            #________________________________________
            Systematic('isr_2', {
                'down': 'weight * weight_isr_2_down',
                'up': 'weight * weight_isr_2_up',
                },
                combine_name='ISR2',
                correlation_eras=1.0,
                correlation_procs='split',
                # tandp=True,
            ),
            Systematic('isr_g2gg_mur', {
                'down': 'weight * weight_isr_g2gg_mur_down',
                'up': 'weight * weight_isr_g2gg_mur_up',
                },
                combine_name='ISRg2ggMuR',
                correlation_eras=1.0,
                correlation_procs='split',
                # tandp=True,
            ),
            Systematic('isr_g2qq_mur', {
                'down': 'weight * weight_isr_g2qq_mur_down',
                'up': 'weight * weight_isr_g2qq_mur_up',
                },
                combine_name='ISRg2qqMuR',
                correlation_eras=1.0,
                correlation_procs='split',
                # tandp=True,
            ),
            Systematic('isr_q2qg_mur', {
                'down': 'weight * weight_isr_q2qg_mur_down',
                'up': 'weight * weight_isr_q2qg_mur_up',
                },
                combine_name='ISRq2qgMuR',
                correlation_eras=1.0,
                correlation_procs='split',
                # tandp=True,
            ),
            Systematic('isr_x2xg_mur', {
                'down': 'weight * weight_isr_x2xg_mur_down',
                'up': 'weight * weight_isr_x2xg_mur_up',
                },
                combine_name='ISRx2xgMuR',
                correlation_eras=1.0,
                correlation_procs='split',
                # tandp=True,
            ),
            Systematic('isr_g2gg_cns', {
                'down': 'weight * weight_isr_g2gg_cns_down',
                'up': 'weight * weight_isr_g2gg_cns_up',
                },
                combine_name='ISRg2ggCNS',
                correlation_eras=1.0,
                correlation_procs='split',
                # tandp=True,
            ),
            Systematic('isr_g2qq_cns', {
                'down': 'weight * weight_isr_g2qq_cns_down',
                'up': 'weight * weight_isr_g2qq_cns_up',
                },
                combine_name='ISRg2qqCNS',
                correlation_eras=1.0,
                correlation_procs='split',
                # tandp=True,
            ),
            Systematic('isr_q2qg_cns', {
                'down': 'weight * weight_isr_q2qg_cns_down',
                'up': 'weight * weight_isr_q2qg_cns_up',
                },
                combine_name='ISRq2qgCNS',
                correlation_eras=1.0,
                correlation_procs='split',
                # tandp=True,
            ),
            Systematic('isr_x2xg_cns', {
                'down': 'weight * weight_isr_x2xg_cns_down',
                'up': 'weight * weight_isr_x2xg_cns_up',
                },
                combine_name='ISRx2xgCNS',
                correlation_eras=1.0,
                correlation_procs='split',
                # tandp=True,
            ),
            #________________________________________
            Systematic('prefire', {
                'down': 'weight / weight_prefire * weight_prefire_down',
                'up': 'weight / weight_prefire * weight_prefire_up',
                },
                # tandp=True,
            ),
            #________________________________________
            Systematic('pu', {
                'down': 'weight / weight_pu * weight_pu_down',
                'up': 'weight / weight_pu * weight_pu_up',
                },
                correlation_eras=1.0,
                tandp=True,
            ),
            #________________________________________
            Systematic('sfmu_id', {
                'down': 'weight / weight_sfmu_id * weight_sfmu_id_down',
                'up': 'weight / weight_sfmu_id * weight_sfmu_id_up',
                },
                combine_name='sfmuID',
                # tandp=True,
            ),
            Systematic('sfmu_iso', {
                'down': 'weight / weight_sfmu_iso * weight_sfmu_iso_down',
                'up': 'weight / weight_sfmu_iso * weight_sfmu_iso_up',
                },
                combine_name='sfmuISO',
            ),
            Systematic('sfmu_trigger', {
                'down': 'weight / weight_sfmu_trigger * weight_sfmu_trigger_down',
                'up': 'weight / weight_sfmu_trigger * weight_sfmu_trigger_up',
                },
                combine_name='sfmuTrigger',
                # tandp=True,
            ),
            Systematic('sfelec_id', {
                'down': 'weight / weight_sfelec_id * weight_sfelec_id_down',
                'up': 'weight / weight_sfelec_id * weight_sfelec_id_up',
                },
                combine_name='sfelecID',
                # tandp=True,
            ),
            Systematic('sfelec_reco', {
                'down': 'weight / weight_sfelec_reco * weight_sfelec_reco_down',
                'up': 'weight / weight_sfelec_reco * weight_sfelec_reco_up',
                },
                combine_name='sfelecReco',
                # tandp=True,
            ),
            Systematic('sfelec_trigger', {
                'down': 'weight / weight_sfelec_trigger * weight_sfelec_trigger_down',
                'up': 'weight / weight_sfelec_trigger * weight_sfelec_trigger_up',
                },
                combine_name='sfelecTrigger',
            ),
            #________________________________________
            Systematic('toppt_a', {
                'down': 'weight / weight_toppt_applied * weight_toppt_a_down',
                'up': 'weight / weight_toppt_applied * weight_toppt_a_up',
                },
                correlation_eras=1.0,
                combine_name='topptA',
                # tandp=True, # it's only a constant normalization uncertainty... should not be varied
            ),
            Systematic('toppt_b', {
                'down': 'weight / weight_toppt_applied * weight_toppt_b_down',
                'up': 'weight / weight_toppt_applied * weight_toppt_b_up',
                },
                correlation_eras=1.0,
                combine_name='topptB',
                # tandp=True,
            ),
            #________________________________________
            Systematic('jesTotal', {
                'down': None,
                'up': None,
                },
                correlation_eras=0.0,
                correlation_procs=None,
                tandp=True,
            ),
            #________________________________________
            Systematic('jer', {
                'down': None,
                'up': None,
                },
                correlation_eras=0.0,
                correlation_procs=None,
                tandp=True,
            ),
            #________________________________________
            Systematic('mtop', {
                'mtop171p5': None,
                'mtop173p5': None,
                },
                correlation_eras=1.0,
                correlation_procs=None,
                # tandp=True,
            ),
            #________________________________________
            Systematic('tune', {
                'down': None,
                'up': None,
                },
                correlation_eras=1.0,
                correlation_procs=None,
                # tandp=True,
            ),
            #________________________________________
            Systematic('hdamp', {
                'down': None,
                'up': None,
                },
                correlation_eras=1.0,
                correlation_procs=None,
                # tandp=True,
            ),
            #________________________________________
            Systematic('cr', {
                'cr1': None,
                'cr2': None,
                'erdon': None,
                },
                correlation_eras=1.0,
                correlation_procs=None,
                # tandp=True,
            ),
        ]
        if self.include_jes_splits:
            for name, corr_eras in _JECSMEAR_SOURCES.items():
                self.base.append(Systematic('jes'+name, {
                    'down': None,
                    'up': None,
                    },
                    correlation_eras=corr_eras,
                    correlation_procs=None,
                    # tandp=True,
                ))
        self.base = {syst.name: syst for syst in self.base}

        ### whitelist has higher priority than blacklist. Only works with "startswith"-matches
        keys = [k for k in self.base.keys()]
        for black in blacklist:
            for k in keys:
                if k.startswith(black):
                    if len(whitelist)!=0:
                        for white in whitelist:
                            if not k.startswith(white):
                                del self.base[k]
                    else:
                        del self.base[k]
        if len(blacklist)==0:
            for k in keys:
                for white in whitelist:
                    if not k.startswith(white):
                        del self.base[k]

    def get_all_variations(self, include_nominal=True):

        all_variations = {}
        all_variations[self.nominal.name] = self.nominal
        for syst in self.base.values():
            for variation in syst.variations.values():
                all_variations[variation.name] = variation
        return all_variations


def get_variable_binning_xlabel_xunit(variable_name, tagger_name=None, fit_variable=False):

    logy = False
    leg_offset_x = 0.
    leg_offset_y = 0.

    # binning = None # FIXME: other variables than mSD and mass use other binning than the following one...
    if variable_name=='mSD' or variable_name=='mass':

        if tagger_name==None:
            sys.exit('Tagger name not given. Cannot determine binning for mSD or mass')
        elif '_t_' in tagger_name:
            if fit_variable:
                binning = np.array([50, 70, 85, 105, 120, 140, 155, 170, 185, 200, 210, 220, 230, 250, 275, 300, 350, 400, 450, 500], dtype=float)
            else:
                binning = np.linspace(0, 500, num=51)
        elif '_w_' in tagger_name:
            if fit_variable:
                binning = np.array([50, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 130, 140, 155, 170, 185, 200, 250, 300], dtype=float)
            else:
                binning = np.linspace(0, 300, num=51)
        else:
            binning = None
            sys.exit('Not a W or top tagger: no defined binning')

        if variable_name=='mSD':
            xlabel = 'Probe jet #it{m}_{SD}'
        elif variable_name=='mass':
            xlabel = 'Probe jet #it{m}_{jet}'
        xunit = 'GeV'

    elif variable_name=='tau32':
        binning = np.linspace(0, 1, num=51)
        xlabel = 'Probe jet #tau_{3}/#tau_{2}'
        xunit = None
        leg_offset_x = -0.35

    elif variable_name=='maxDeepJet':#, 'maxDeepCSV', 'MDdeepak8_TvsQCD', 'ParticleNet_WvsQCD', 'fpt1']:
        binning = np.linspace(0, 1, num=51)
        xlabel = 'Max. DeepJet score of probe subjets'
        xunit = None
        logy = True
        leg_offset_x = -0.23
        leg_offset_y = 0.06

    elif variable_name=='maxDeepCSV':#, 'maxDeepCSV', 'MDdeepak8_TvsQCD', 'ParticleNet_WvsQCD', 'fpt1']:
        binning = np.linspace(0, 1, num=51)
        xlabel = 'Max. DeepCSV score of probe subjets'
        xunit = None
        logy = True
        leg_offset_x = -0.23
        leg_offset_y = 0.06

    elif variable_name=='MDDeepAK8_TvsQCD':#, 'maxDeepCSV', 'MDdeepak8_TvsQCD', 'ParticleNet_WvsQCD', 'fpt1']:
        binning = np.linspace(0, 1, num=51)
        xlabel = 'Probe jet MD-DeepAK8 TvsQCD score'
        xunit = None
        logy = True
        leg_offset_x = -0.23
        leg_offset_y = 0.06

    elif variable_name=='ParticleNet_WvsQCD':#, 'maxDeepCSV', 'MDdeepak8_TvsQCD', 'ParticleNet_WvsQCD', 'fpt1']:
        binning = np.linspace(0, 1, num=51)
        xlabel = 'Probe jet ParticleNet WvsQCD score'
        xunit = None
        logy = True
        leg_offset_x = -0.23
        leg_offset_y = 0.06

    elif variable_name=='fpt1':
        binning = np.linspace(0, 1, num=51)
        xlabel = '#it{p}_{T} fraction of leading probe subjet'
        xunit = None
        leg_offset_x = -0.41
        leg_offset_y = 0.06

    elif variable_name=='nsub':
        binning = np.linspace(-0.5, 10.5, num=12)
        xlabel = 'Probe subjet multiplicity'
        xunit = None

    elif variable_name=='mpair':
        binning = np.linspace(0, 250, num=51)
        xlabel = 'Min. #it{m}_{ij} of three leading probe subjets'
        xunit = 'GeV'

    elif variable_name=='pt':
        binning = np.linspace(0, 1000, num=51)
        xlabel = 'Probe jet #it{p}_{T}'
        xunit = 'GeV'

    else:
        binning = None
        sys.exit('Variable "{}" has no defined binning'.format(variable_name))

    return binning, xlabel, xunit, logy, leg_offset_x, leg_offset_y


class Primes():

    def __init__(self):
        self.calls_to_get_next_prime = 0

    def get_next_prime(self):
        self.calls_to_get_next_prime += 1
        return sympy.prime(self.calls_to_get_next_prime)


class NormFactsQCD():

    def __init__(self, years=None):

        possible_years = ['UL16preVFP', 'UL16postVFP', 'UL17', 'UL18']
        if years is None:
            self.years = possible_years
        elif isinstance(years, list):
            self.years = years
        elif isinstance(years, str) and years in possible_years:
            self.years = years
        else:
            sys.exit('Check years argument')

        variations = [
            'murmuf_upup',
            'murmuf_upnone',
            'murmuf_noneup',
            'murmuf_downdown',
            'murmuf_downnone',
            'murmuf_nonedown',
        ]

        self.base_path_to_json_files = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/HighPtSingleTop/output/Analysis/normFacts/')
        self.qcd_norm_dict = {
            'dataset': [],
        }
        for variation in variations:
            self.qcd_norm_dict[variation] = []

        for year in self.years:

            json_path = os.path.join(self.base_path_to_json_files, year, 'norm_facts/*_norm.json')
            json_list = glob.glob(json_path)

            for json_file_path in json_list:
                dataset_name = os.path.basename(json_file_path).replace('_norm.json', '')
                with open(json_file_path) as json_file:
                    json_dict = json.load(json_file)
                    self.qcd_norm_dict['dataset'].append(dataset_name)
                    for variation in variations:
                        self.qcd_norm_dict[variation].append(json_dict['QCD'][variation])
        self.qcd_norm_df = pd.DataFrame(self.qcd_norm_dict)
        # print(self.qcd_norm_df)

normFactsQCD = NormFactsQCD()

# if __name__ == '__main__':
#
#     n = NormFactsQCD()
