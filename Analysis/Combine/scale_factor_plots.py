import os
import sys
import subprocess
import ROOT as root
from copy import deepcopy
from collections import OrderedDict
from timeit import default_timer as timer
import numpy as np
import warnings
import re
from tqdm import tqdm

import sys
sys.path.append(os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/Analysis'))
from constants import Systematics, _BANDS, _TAGGERS, _PT_INTERVALS_TANDP_AK8_T, _PT_INTERVALS_TANDP_AK8_W, _PT_INTERVALS_TANDP_HOTVR, _YEARS, get_variable_binning_xlabel_xunit, Primes



all_years = [
'UL16preVFP',
'UL16postVFP',
'UL17',
'UL18',
]
regions = [
'Pass',
'Fail',
]
channels = [
'muo',
'ele',
]
# prepostfits = [
# 'prefitRaw',
# # 'postfitCombine',
# ]
taggers = [
'ak8_t__tau',
'ak8_t_btagDJet__tau',
'ak8_t_btagDCSV__tau',
'hotvr_t__tau',
'ak8_w__partnet',
'ak8_t__MDdeepak8',
]
taggers = {k: _TAGGERS[k] for k in taggers}

class ScaleFactorPlots():

    def __init__(self,
        tagger_name,
        wp,
        # mode='TagAndProbe',
        mscSplitting = 'mscTop3',
        include_total_range = False,
        # years = None, # either list of strings or single string (will be converted to list with one string element)
        # task_name_suffix = '',
    ):

        pass



        self.task_name_suffix += task_name_suffix # FIXME be more creative with this name
        self.task_name = '-'.join(['combineTask', self.tagger.name, self.wp.name]) + ('-'+self.task_name_suffix if len(self.task_name_suffix) else '')
        self.workdir = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel/combine', self.tagger.name, 'workdirs', self.task_name)


    def read_scale_factors(self):

        pass
