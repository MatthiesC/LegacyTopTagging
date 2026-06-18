#!/usr/bin/env python2

from __future__ import print_function

import os

import ROOT as root
# from copy import deepcopy
# from collections import OrderedDict
# import numpy as np

import sys
sys.path.append(os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/Analysis'))
from constants import _YEARS


all_years = [
'UL16preVFP',
'UL16postVFP',
'UL17',
'UL18',
]

all_channels = [
'muo',
'ele',
]

base_path = os.path.join(
    os.environ.get("CMSSW_BASE"),
    'src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel',
)


name_of_histo = sys.argv[1]
if len(sys.argv) > 1:
    binx1 = int(sys.argv[2])
    binx2 = int(sys.argv[3])


n_data = 0
n_data_muo = 0
n_data_ele = 0
for year in all_years:
    n_data_ele_tmp = 0
    n_data_muo_tmp = 0
    for channel in all_channels:
        infile_name = 'uhh2.AnalysisModuleRunner.DATA.DATA.root'
        infile_path = os.path.join(
            base_path,
            year,
            channel,
            'nominal/hadded',
            infile_name,
        )
        infile = root.TFile.Open(infile_path, 'READ')
        # inhist = infile.Get('Presel_Main_Common/count')
        # inhist = infile.Get('Presel_Main_HOTVR/hotvrjet1_eta')
        inhist = infile.Get(name_of_histo)
        # n_data_tmp = inhist.GetBinContent(1)
        # n_data_tmp = inhist.Integral(31, 101)
        n_data_tmp = inhist.Integral(binx1, binx2)
        print(year, channel, n_data_tmp)
        n_data += n_data_tmp
        if channel == 'ele':
            n_data_ele_tmp = n_data_tmp
            n_data_ele += n_data_tmp
        if channel == 'muo':
            n_data_muo_tmp = n_data_tmp
            n_data_muo += n_data_tmp
    print('muo/ele', n_data_muo_tmp / n_data_ele_tmp)

print(n_data)
print('muo/ele', n_data_muo / n_data_ele)
