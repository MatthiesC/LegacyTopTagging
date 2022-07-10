import os
import sys
import argparse
from tqdm import tqdm

import uproot
# import awkward as ak
import numpy as np
from hist import Hist
import hist

sys.path.append(os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/Analysis'))
from constants import _TAGGERS, _BANDS


all_years = ['UL16preVFP', 'UL16postVFP', 'UL17', 'UL18']
all_channels = ['muo']

taggers = ['ak8_t__tau', 'hotvr_t__tau', 'ak8_w__partnet']
taggers = {k: _TAGGERS[k] for k in taggers}

if not sys.argv[1:]: sys.exit('No arguments provided. Exit.')
parser = argparse.ArgumentParser()
parser.add_argument('-y', '--years', nargs='+', choices=all_years, default=all_years)
parser.add_argument('-c', '--channels', nargs='+', choices=all_channels, default=all_channels)
parser.add_argument('-t', '--tagger', default='ak8_t__tau', choices=taggers.keys())
args = parser.parse_args(sys.argv[1:])

processes = [
'DATA.DATA',
'MC.QCD',
'MC.nonQCD',
]

for year in args.years:

    for channel in args.channels:

        workDir = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel', year, channel)
        inputDir = os.path.join(workDir, 'nominal/hadded/')
        inputFilePathBase = os.path.join(inputDir, 'uhh2.AnalysisModuleRunner.{0}.root')

        if args.tagger.startswith('ak8_t'):
            tagger_base = 'ak8_t'
            base_cut = '(output_has_probejet_AK8 == True) & (output_probejet_AK8_pt > 300)'
        elif args.tagger.startswith('ak8_w'):
            tagger_base = 'ak8_w'
            base_cut = '(output_has_probejet_AK8 == True) & (output_probejet_AK8_pt > 200)'
        elif args.tagger.startswith('hotvr_t'):
            tagger_base = 'hotvr_t'
            base_cut = '(output_has_probejet_HOTVR == True) & (output_probejet_HOTVR_pt > 200)'
        # else:
        #     tagger_base = None
        #     base_cut = None

        batches = {}
        hists = {}
        mtw_binning = np.linspace(0, 200, 41) # 40 bins
        for process in processes:
            batches[process] = [inputFilePathBase.format(process)+':AnalysisTree']
            for band_k, band in _BANDS.items():
                hists.setdefault(band_k, {})[process] = Hist(
                    hist.axis.Variable(
                        mtw_binning, label='mtw', name='mtw', underflow=True, overflow=True
                    ),
                    storage=hist.storage.Weight(), # like ROOT's Sumw2()
                )

        for process in processes:

            for batch in tqdm(uproot.iterate(batches[process], expressions=['weight', 'output_mtw', 'band'], cut=base_cut, library='pd')):

                weight = batch['weight'].to_numpy()
                mtw = batch['output_mtw'].to_numpy()

                for band_k, band in _BANDS.items():
                    hists[band_k][process].fill(mtw[batch['band'] == band.index], weight=weight[batch['band'] == band.index])

        outputFileName = 'qcd_normalization_histograms__'+tagger_base+'.root'
        outputDir = os.path.join(workDir, 'qcd_normalization_study')
        os.system('mkdir -p '+outputDir)
        outputFilePath = os.path.join(outputDir, outputFileName)

        with uproot.recreate(outputFilePath) as outfile:
            for band_k, value in hists.items():
                for process, histogram in value.items():
                    out_name = band_k+'/'+process
                    outfile[out_name] = histogram
