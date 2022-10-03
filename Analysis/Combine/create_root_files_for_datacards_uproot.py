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
from constants import WorkingPoint, _TAGGERS, _DEEPCSV_WPS, _DEEPJET_WPS, _BANDS, _SYSTEMATICS, _PT_INTERVALS_TANDP_HOTVR, _PT_INTERVALS_TANDP_AK8_T, _PT_INTERVALS_TANDP_AK8_W


all_years = ['UL16preVFP', 'UL16postVFP', 'UL17', 'UL18']
all_channels = ['muo']

taggers = ['ak8_t__tau', 'ak8_t_btagDJet__tau', 'hotvr_t__tau']
taggers = {k: _TAGGERS[k] for k in taggers}

# set working points:
for tagger_k, tagger in taggers.items():
    if tagger_k == 'ak8_t__tau' or tagger_k == 'ak8_t_btagDJet__tau' or tagger_k == 'ak8_t_btagDCSV__tau':
        tagger.wps = [
            WorkingPoint(0.001, 0.38),
            WorkingPoint(0.005, 0.47),
            WorkingPoint(0.010, 0.52),
            WorkingPoint(0.025, 0.61),
            WorkingPoint(0.050, 0.69),
        ]
        tagger.var_intervals = _PT_INTERVALS_TANDP_AK8_T
        tagger.fit_variable = 'output_probejet_AK8_mSD'
    if tagger_k == 'hotvr_t__tau':
        hotvr_wp = WorkingPoint(9.999, 0.56) # providing bkg_eff does not make sense here; I have not derived it on my own; should be 3% or so
        hotvr_wp.name = 'Standard'
        tagger.wps = [hotvr_wp]
        tagger.var_intervals = _PT_INTERVALS_TANDP_HOTVR
        tagger.fit_variable = 'output_probejet_HOTVR_mass'

if not sys.argv[1:]: sys.exit('No arguments provided. Exit.')
parser = argparse.ArgumentParser()
parser.add_argument('-y', '--years', nargs='+', choices=all_years, default=all_years)
parser.add_argument('-c', '--channels', nargs='+', choices=all_channels, default=all_channels)
parser.add_argument('-t', '--tagger', default='ak8_t__tau', choices=taggers.keys())
args = parser.parse_args(sys.argv[1:])

years = args.years
channels = args.channels
the_tagger = taggers[args.tagger]

# MSc = Merge Scenario
processes = [
# 'DATA',
# 'QCD',
# 'nonQCD',
# 'TTbar',
'TTbar__MSc_FullyMerged',
# # 'TTbar__MSc_WMerged',
# # 'TTbar__MSc_QBMerged',
# # 'TTbar__MSc_NotMerged',
# 'TTbar__MSc_YllufMerged',
# # 'TTbar__MSc_SemiMerged',
# 'TTbar__MSc_Background',
# 'ST',
# 'ST__MSc_FullyMerged',
# # 'ST__MSc_WMerged',
# # 'ST__MSc_QBMerged',
# # 'ST__MSc_NotMerged',
# 'ST__MSc_YllufMerged',
# # 'ST__MSc_SemiMerged',
# 'ST__MSc_Background',
# 'VJetsAndVV',
]

regions = [
    'Pass',
    'Fail',
]

bands = _BANDS


def create_input_hists(tagger, year, variable):

    probejetalgo = ''
    if tagger.name.startswith('ak8'):
        probejetalgo = 'AK8'
    elif tagger.name.startswith('hotvr'):
        probejetalgo = 'HOTVR'

    # binning = None # FIXME: other variables than mSD and mass use other binning than the following one...
    if variable.endswith('_mSD') or variable.endswith('_mass'):
        binning = np.array([50, 70, 85, 105, 120, 140, 155, 170, 185, 200, 210, 220, 230, 250, 275, 300, 350, 400, 450, 500], dtype=float)
    else:
        binning = None
        sys.exit('Variable "{}" has no defined binning'.format(variable))

    inDirBase = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel', year)
    outDirBase = os.path.join(inDirBase, 'combine', tagger.name)

    # cut_rule = '(output_has_probejet_{} == True)'.format(probejetalgo)

    for wp in tagger.wps:

        tandp_rule = tagger.tandp_rule.replace('WP_VALUE', '{:.5f}'.format(wp.cut_value)).replace('DEEPCSV_SCORE', '{:.5f}'.format(_DEEPCSV_WPS[year]['loose'])).replace('DEEPJET_SCORE', '{:.5f}'.format(_DEEPJET_WPS[year]['loose']))

        for pt_bin in tagger.var_intervals.values():

            task_name = '-'.join([tagger.name, wp.name, pt_bin.name, year, variable])
            outDir = os.path.join(outDirBase, task_name)
            os.system('mkdir -p '+outDir)
            outFileName = 'BasicHists-'+task_name+'.root'
            outFilePath = os.path.join(outDir, outFileName)

            batches = {}
            hists = {}

            for band_k, band in bands.items():

                for region in regions: # Pass or Fail

                    for channel in channels:

                        for process in processes:

                            for syst, syst_v in _SYSTEMATICS.items():

                                # if process == 'DATA' and syst != 'nominal':
                                #     continue

                                inDir = os.path.join(inDirBase, channel, 'nominal', 'hadded')
                                if not syst_v.weight_based:
                                    inDir.replace('/nominal/', '/syst_'+syst+'/')
                                process_ = process.split('__MSc_')[0]
                                inFileName = 'uhh2.AnalysisModuleRunner.{MCDATA}.{PROCESS}.root'.format(MCDATA=('DATA' if process == 'DATA' else 'MC'), PROCESS=process_)
                                inFilePath = os.path.join(inDir, inFileName)
                                if not syst_v.weight_based and not os.path.isfile(inFilePath):
                                    # if for this process, the extra syst file does not exist (e.g. TTbar systs for DYJets sample...), then fallback to nominal
                                    inDir.replace('/syst_'+syst+'/', '/nominal/')
                                    inFilePath.replace('/syst_'+syst+'/', '/nominal/')

                                batches.setdefault(band_k, {}).setdefault(region, {}).setdefault(channel, {}).setdefault(process, {})[syst] = [
                                    inFilePath+':AnalysisTree'
                                ]

                                # declare Hist object:
                                hists.setdefault(band_k, {}).setdefault(region, {}).setdefault(channel, {}).setdefault(process, {})[syst] = Hist(
                                    hist.axis.Variable(
                                        binning, label=variable, name=variable, underflow=True, overflow=True
                                    ),
                                    storage=hist.storage.Weight(), # like ROOT's Sumw2()
                                )

            for band_k, band in bands.items():

                print('Working on band:', band_k)

                # cut_rule += '& (band == {})'.format(str(band.index))
                #
                # # write the pt cut rule here instead of in the previous for-loop step since it is probably computationally more effective to have the band cut rule first
                # cut_rule += '& (output_probejet_{}_pt > {:.5f})'.format(probejetalgo, pt_bin.var_min)
                # cut_rule += '& (output_probejet_{}_pt <= {:.5f})'.format(probejetalgo, pt_bin.var_max)

                for region in regions:

                    print('Working on region:', region)

                    # # mass cuts are not included as defined in constants.py
                    # if region == 'Pass':
                    #     cut_rule += '& ({})'.format(tandp_rule)
                    # elif region == 'Fail':
                    #     cut_rule += '& ~({})'.format(tandp_rule) # tilde is like a "not"

                    for channel in channels:

                        print('Working on channel:', channel)

                        for process in processes:

                            print('Working on process:', process)

                            cut_rule = '(output_has_probejet_{} == True)'.format(probejetalgo)

                            cut_rule += ' & (band == {})'.format(str(band.index))

                            # write the pt cut rule here instead of in the previous for-loop step since it is probably computationally more effective to have the band cut rule first
                            cut_rule += ' & (output_probejet_{}_pt > {:.5f})'.format(probejetalgo, pt_bin.var_min)
                            cut_rule += ' & (output_probejet_{}_pt <= {:.5f})'.format(probejetalgo, pt_bin.var_max)

                            # mass cuts are not included as defined in constants.py
                            if region == 'Pass':
                                cut_rule += ' & ({})'.format(tandp_rule)
                            elif region == 'Fail':
                                cut_rule += ' & ~({})'.format(tandp_rule) # tilde is like a "not"

                            msc = process.split('__MSc_')[-1]
                            if msc == 'FullyMerged':
                                cut_rule += ' & (output_merge_scenario_{0} == 1)'.format(probejetalgo)
                            elif msc == 'WMerged':
                                cut_rule += ' & (output_merge_scenario_{0} == 2)'.format(probejetalgo)
                            elif msc == 'QBMerged':
                                cut_rule += ' & (output_merge_scenario_{0} == 3)'.format(probejetalgo)
                            elif msc == 'NotMerged':
                                cut_rule += ' & (output_merge_scenario_{0} == 4)'.format(probejetalgo)
                            elif msc == 'YllufMerged':
                                cut_rule += ' & ((output_merge_scenario_{0} == 2) | (output_merge_scenario_{0} == 3) | (output_merge_scenario_{0} == 4))'.format(probejetalgo)
                            elif msc == 'SemiMerged':
                                cut_rule += ' & ((output_merge_scenario_{0} == 2) | (output_merge_scenario_{0} == 3))'.format(probejetalgo)
                            elif msc == 'Background':
                                cut_rule += ' & (output_merge_scenario_{0} == -1)'.format(probejetalgo)

                            print(cut_rule)

                            for syst, syst_v in _SYSTEMATICS.items():

                                print('Working on systematic:', syst)

                                weight_alias = syst_v.weight_alias
                                expressions = [
                                    variable,
                                    'the_weight',
                                    'output_has_probejet_{}'.format(probejetalgo),
                                    'band',
                                    'output_probejet_{}_pt'.format(probejetalgo),
                                    'output_merge_scenario_{}'.format(probejetalgo),
                                ]
                                # manually add the variables used in the tandp_rule to the expressions list (pretty hacky):
                                for var_string in tandp_rule.replace('(', ' ').replace(')', ' ').split(' '):
                                    if var_string.startswith('output'):
                                        expressions.append(var_string)

                                for batch in tqdm(uproot.iterate(batches[band_k][region][channel][process][syst], expressions=expressions, aliases={'the_weight': weight_alias}, cut=cut_rule, library='pd')):

                                    arr_weight = batch['the_weight'].to_numpy()
                                    arr_variable = batch[variable].to_numpy()

                                    hists[band_k][region][channel][process][syst].fill(arr_variable, weight=arr_weight)


            print('Writing output file...')
            with uproot.recreate(outFilePath) as outFile:

                for band_k, band in bands.items():

                    for region in regions:

                        for channel in channels:

                            for process in processes:

                                for syst in _SYSTEMATICS.keys():

                                    histName = '/'.join([band_k, region, channel, process + '__' + syst])
                                    outFile[histName] = hists[band_k][region][channel][process][syst]
            print('Wrote', outFilePath)

if __name__ == '__main__':

    for year in years:

        create_input_hists(the_tagger, year, the_tagger.fit_variable)
