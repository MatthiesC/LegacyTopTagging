#!/usr/bin/env /nfs/dust/cms/user/matthies/anaconda3/envs/ana3/bin/python

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
# from constants import WorkingPoint, _TAGGERS, _DEEPCSV_WPS, _DEEPJET_WPS, _BANDS, _SYSTEMATICS, _PT_INTERVALS_TANDP_HOTVR, _PT_INTERVALS_TANDP_AK8_T, _PT_INTERVALS_TANDP_AK8_W
from constants import WorkingPoint, _TAGGERS, _BANDS, _PT_INTERVALS_TANDP_HOTVR, _PT_INTERVALS_TANDP_AK8_T, _PT_INTERVALS_TANDP_AK8_W, Systematics, get_variable_binning_xlabel_xunit, _NULL_WP
from constants import normFactsQCD

# systematics = Systematics(blacklist=['sfelec', 'sfmu_iso'])
# systematics = Systematics(blacklist=['sfmu_iso'])
systematics = Systematics(blacklist=['sfelec_trigger', 'sfmu_iso'])
_SYSTEMATICS = systematics.get_all_variations()
# for k,v in _SYSTEMATICS.items():
#     print(k, v.weight_based)

from parallel_threading import run_with_pool


all_years = ['UL16preVFP', 'UL16postVFP', 'UL17', 'UL18']
all_channels = [
'muo',
'ele',
]

taggers = ['ak8_t__tau', 'ak8_t_btagDJet__tau', 'ak8_t_btagDCSV__tau', 'hotvr_t__tau', 'ak8_w__partnet', 'ak8_t__MDdeepak8']
taggers = {k: _TAGGERS[k] for k in taggers}

# # set working points:
# for tagger_k, tagger in taggers.items():
#     if tagger_k == 'ak8_t__tau' or tagger_k == 'ak8_t_btagDJet__tau' or tagger_k == 'ak8_t_btagDCSV__tau' or tagger_k == 'ak8_t__MDdeepak8':
#         # tagger.wps = [
#         #     WorkingPoint(0.001, 0.38),
#         #     WorkingPoint(0.005, 0.47),
#         #     WorkingPoint(0.010, 0.52),
#         #     WorkingPoint(0.025, 0.61),
#         #     WorkingPoint(0.050, 0.69),
#         # ]
#         tagger.var_intervals = _PT_INTERVALS_TANDP_AK8_T
#         tagger.fit_variable = 'output_probejet_AK8_mSD'
#     if tagger_k == 'hotvr_t__tau':
#         # hotvr_wp = WorkingPoint(9.999, 0.56) # providing bkg_eff does not make sense here; I have not derived it on my own; should be 3% or so
#         # hotvr_wp.name = 'Standard'
#         # tagger.wps = [hotvr_wp]
#         tagger.var_intervals = _PT_INTERVALS_TANDP_HOTVR
#         tagger.fit_variable = 'output_probejet_HOTVR_mass'
#     if tagger_k == 'ak8_w__partnet':
#         # tagger.wps = [
#         #     WorkingPoint(0.030, 0.00), # cut value year dependent! Need to hack it in the code below where we define the cut rules...
#         # ]
#         tagger.var_intervals = _PT_INTERVALS_TANDP_AK8_W
#         tagger.fit_variable = 'output_probejet_AK8_mSD'

# MSc = Merge Scenario
processes = [
'DATA',
'QCD',
# 'nonQCD', # I don't have those hadded
'TTbar',
'TTbar__MSc_FullyMerged',
'TTbar__MSc_WMerged',
'TTbar__MSc_QBMerged',
'TTbar__MSc_NotMerged',
'TTbar__MSc_YllufMerged',
'TTbar__MSc_SemiMerged',
'TTbar__MSc_NotTopOrWMerged',
'TTbar__MSc_Background',
'ST',
'ST__MSc_FullyMerged',
'ST__MSc_WMerged',
'ST__MSc_QBMerged',
'ST__MSc_NotMerged',
'ST__MSc_YllufMerged',
'ST__MSc_SemiMerged',
'ST__MSc_NotTopOrWMerged',
'ST__MSc_Background',
'VJetsAndVV',
]

regions = [
    'Pass',
    'Fail',
]

bands = _BANDS


if not sys.argv[1:]: sys.exit('No arguments provided. Exit.')
parser = argparse.ArgumentParser()
parser.add_argument('-v', '--vars', nargs='+') # Optional; by default, will use the_tagger.fit_variable
parser.add_argument('-t', '--tagger', default='ak8_t__tau', choices=taggers.keys())
parser.add_argument('-w', '--wps', nargs='+', default=['0']) # integer indices of wps, following the indexing given at the top of this python script
parser.add_argument('-y', '--years', nargs='+', choices=all_years, default=all_years)
parser.add_argument('-s', '--systs', nargs='+', default=['nominal'])
# parser.add_argument('-b', '--bands', nargs='+', choices=bands.keys())
# parser.add_argument('-r', '--regions', nargs='+', choices=regions)
parser.add_argument('-c', '--channels', nargs='+', choices=all_channels, default=all_channels)
args = parser.parse_args(sys.argv[1:])

years = args.years
channels = args.channels
the_tagger = taggers[args.tagger]

def create_input_hists(variable, tagger, year, arg_wp_index=None, arg_syst=None, arg_band=None, arg_region=None, arg_channel=None, arg_fit_variable=True):

    probejetalgo = ''
    if tagger.name.startswith('ak8'):
        probejetalgo = 'AK8'
    elif tagger.name.startswith('hotvr'):
        probejetalgo = 'HOTVR'

    binning, _, _, _, _, _ = get_variable_binning_xlabel_xunit(variable_name=variable.split('_'+probejetalgo+'_')[-1], tagger_name=tagger.name, fit_variable=arg_fit_variable)

    # # binning = None # FIXME: other variables than mSD and mass use other binning than the following one...
    # if variable.endswith('_mSD') or variable.endswith('_mass'):
    #     if '_t_' in tagger.name:
    #         binning = np.array([50, 70, 85, 105, 120, 140, 155, 170, 185, 200, 210, 220, 230, 250, 275, 300, 350, 400, 450, 500], dtype=float)
    #     elif '_w_' in tagger.name:
    #         binning = np.array([50, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 130, 140, 155, 170, 185, 200, 250, 300], dtype=float)
    #     else:
    #         binning = None
    #         sys.exit('Not a W or top tagger: no defined binning')
    # else:
    #     binning = None
    #     sys.exit('Variable "{}" has no defined binning'.format(variable))

    inDirBase = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel', year)
    outDirBase = os.path.join(inDirBase, 'combine', tagger.name)

    # cut_rule = '(output_has_probejet_{} == True)'.format(probejetalgo)

    # for wp in ([arg_wp] if arg_wp else tagger.wps):

        # #___________
        # # HACK since my custom PartNet WP cuts are year dependent:
        # if tagger.name == 'ak8_w__partnet':
        #     if year == 'UL16preVFP':
        #         wp.cut_value = 0.871
        #     elif year == 'UL16postVFP':
        #         wp.cut_value = 0.868
        #     elif year == 'UL17':
        #         wp.cut_value = 0.868
        #     elif year == 'UL18':
        #         wp.cut_value = 0.864
        # #___________ HACK END

        # tandp_rule = tagger.tandp_rule.replace('WP_VALUE', '{:.5f}'.format(wp.cut_value)).replace('DEEPCSV_SCORE', '{:.5f}'.format(_DEEPCSV_WPS[year]['loose'])).replace('DEEPJET_SCORE', '{:.5f}'.format(_DEEPJET_WPS[year]['loose']))

    if arg_wp_index != -1:
        tandp_rule = tagger.get_tandp_rule(arg_wp_index, year)
        wp = tagger.get_wp(wp_index, year)
    else:
        tandp_rule = 'True'
        wp = _NULL_WP

    for pt_bin in tagger.var_intervals.values():

        if arg_wp_index == -1 and pt_bin.total_range != True:
            continue

        batches = {}
        hists = {}

        for syst, syst_v in ({arg_syst.name: arg_syst}.items() if arg_syst else _SYSTEMATICS.items()):

            for band_k, band in ({arg_band.name: arg_band}.items() if arg_band else bands.items()):

                for region in ([arg_region] if arg_region else regions): # Pass or Fail
                    if region == 'Fail' and arg_wp_index == -1:
                        continue

                    for channel in ([arg_channel] if arg_channel else channels):

                        for process in processes:

                            # if process == 'DATA' and syst != 'nominal':
                            #     continue

                            inDir = os.path.join(inDirBase, channel, 'nominal', 'hadded')
                            if not syst_v.weight_based:
                                inDir = inDir.replace('/nominal/', '/syst_'+syst+'/')
                            process_ = process.split('__MSc_')[0]
                            inFileName = 'uhh2.AnalysisModuleRunner.{MCDATA}.{PROCESS}.root'.format(MCDATA=('DATA' if process == 'DATA' else 'MC'), PROCESS=process_)
                            # print(inFileName)
                            inFilePath = os.path.join(inDir, inFileName)
                            if not syst_v.weight_based and not os.path.isfile(inFilePath):
                                # if for this process, the extra syst file does not exist (e.g. TTbar systs for DYJets sample...), then fallback to nominal
                                inDir = inDir.replace('/syst_'+syst+'/', '/nominal/')
                                inFilePath = inFilePath.replace('/syst_'+syst+'/', '/nominal/')

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

        for syst, syst_v in ({arg_syst.name: arg_syst}.items() if arg_syst else _SYSTEMATICS.items()):

            is_murmuf = syst.startswith('murmuf')

            print('Working on systematic:', syst)

            for band_k, band in ({arg_band.name: arg_band}.items() if arg_band else bands.items()):

                print('Working on band:', band_k)

                # cut_rule += '& (band == {})'.format(str(band.index))
                #
                # # write the pt cut rule here instead of in the previous for-loop step since it is probably computationally more effective to have the band cut rule first
                # cut_rule += '& (output_probejet_{}_pt > {:.5f})'.format(probejetalgo, pt_bin.var_min)
                # cut_rule += '& (output_probejet_{}_pt <= {:.5f})'.format(probejetalgo, pt_bin.var_max)

                for region in ([arg_region] if arg_region else regions): # Pass or Fail
                    if region == 'Fail' and arg_wp_index == -1:
                        continue

                    print('Working on region:', region)

                    # # mass cuts are not included as defined in constants.py
                    # if region == 'Pass':
                    #     cut_rule += '& ({})'.format(tandp_rule)
                    # elif region == 'Fail':
                    #     cut_rule += '& ~({})'.format(tandp_rule) # tilde is like a "not"

                    for channel in ([arg_channel] if arg_channel else channels):

                        print('Working on channel:', channel)

                        for process in processes:

                            is_data = (process == 'DATA')

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
                            elif msc == 'NotTopOrWMerged':
                                cut_rule += ' & ((output_merge_scenario_{0} == 3) | (output_merge_scenario_{0} == 4))'.format(probejetalgo)
                            elif msc == 'Background':
                                cut_rule += ' & (output_merge_scenario_{0} == -1)'.format(probejetalgo)

                            print(cut_rule)

                            weight_alias = syst_v.weight_alias
                            expressions = [
                                # variable,
                                'the_weight',
                                # 'dataset',
                                'output_has_probejet_{}'.format(probejetalgo),
                                'band',
                                'output_probejet_{}_pt'.format(probejetalgo),
                                'output_merge_scenario_{}'.format(probejetalgo),
                            ]

                            if is_murmuf:
                                expressions.append('dataset')

                            if variable not in expressions: # avoid error due to having same variable twice in expressions list
                                expressions.append(variable)

                            # manually add the variables used in the tandp_rule to the expressions list (pretty hacky):
                            for var_string in tandp_rule.replace('(', ' ').replace(')', ' ').split(' '):
                                if var_string.startswith('output'):
                                    expressions.append(var_string)

                            for batch in tqdm(uproot.iterate(batches[band_k][region][channel][process][syst], expressions=expressions, aliases={'the_weight': weight_alias}, cut=cut_rule, library='pd')):

                                if is_murmuf and not is_data: # no murmuf norm factors available for data
                                    batch['dataset'] = batch['dataset'].apply(lambda dataset: dataset.replace('__AllMergeScenarios', '')) # get rid of __AllMergeScenarios in the dataset names
                                    batch = batch.merge(normFactsQCD.qcd_norm_df, on='dataset')
                                    norm_factor_name = syst
                                    batch['the_weight'] = batch['the_weight'] * batch[norm_factor_name]

                                arr_weight = batch['the_weight'].to_numpy()
                                arr_variable = batch[variable].to_numpy()

                                hists[band_k][region][channel][process][syst].fill(arr_variable, weight=arr_weight)



            for syst, syst_v in ({arg_syst.name: arg_syst}.items() if arg_syst else _SYSTEMATICS.items()):

                task_name = '-'.join([tagger.name, wp.name, pt_bin.name, year, syst, variable])
                # outDir = os.path.join(outDirBase, task_name)
                outDir = os.path.join(outDirBase, 'BasicHists')
                os.system('mkdir -p '+outDir)
                outFileName = 'BasicHists-'+task_name+'.root'
                outFilePath = os.path.join(outDir, outFileName)

                print('Writing output file...')
                with uproot.recreate(outFilePath) as outFile:

                    for band_k, band in ({arg_band.name: arg_band}.items() if arg_band else bands.items()):

                        for region in ([arg_region] if arg_region else regions): # Pass or Fail
                            if region == 'Fail' and arg_wp_index == -1:
                                continue

                            for channel in ([arg_channel] if arg_channel else channels):

                                for process in processes:

                                    histName = '/'.join([band_k, region, channel, process + '__' + syst])
                                    outFile[histName] = hists[band_k][region][channel][process][syst]
        print('Wrote', outFilePath)

if __name__ == '__main__':

    #_______________________________________
    ## Running ROOT in pool mode doesn't work (probably due to ROOT extensive use of global variables...)

    # sets_of_args = []
    # for year in years:
    #     for wp in the_tagger.wps:
    #         set_of_args = {}
    #         set_of_args['tagger'] = the_tagger
    #         set_of_args['year'] = year
    #         set_of_args['variable'] = the_tagger.fit_variable
    #         set_of_args['wp'] = wp
    #         sets_of_args.append(set_of_args)
    #
    # run_with_pool(function=create_input_hists, sets_of_args=sets_of_args)

    # create_input_hists(variable=the_tagger.fit_variable, tagger=the_tagger, arg_wp=the_tagger.wps[0], year='UL17', arg_syst=_SYSTEMATICS['nominal'], arg_band=_BANDS['Main'], arg_region='Fail', arg_channel='muo')


    # #_______________________________________
    # ## For running individual jobs manually
    #
    # # create_input_hists(variable=the_tagger.fit_variable, tagger=the_tagger, arg_wp=the_tagger.wps[0], year='UL17', arg_syst=_SYSTEMATICS['nominal'])
    #
    # the_vars = [the_tagger.fit_variable]
    # if args.vars:
    #     the_vars = args.vars
    #
    # # the_wps = the_tagger.wps
    # # if args.wps:
    # #     the_wps = [the_tagger.wps[int(x)] for x in args.wps]
    #
    # the_systs = args.systs
    #
    # for var in the_vars:
    #     for year in years:
    #         the_wp_indices = range(0, len(the_tagger.get_wp(year=year)))
    #         if args.wps:
    #             the_wp_indices = [int(x) for x in args.wps]
    #         for wp_index in the_wp_indices:
    #             for syst in the_systs:
    #                 create_input_hists(variable=var, tagger=the_tagger, arg_wp_index=wp_index, year=year, arg_syst=_SYSTEMATICS[syst])


    #_______________________________________
    ## For running individual jobs manually for non-WP plots

    # create_input_hists(variable=the_tagger.fit_variable, tagger=the_tagger, arg_wp=the_tagger.wps[0], year='UL17', arg_syst=_SYSTEMATICS['nominal'])

    # the_vars = [the_tagger.fit_variable]
    the_vars = []
    if the_tagger.name.startswith('ak8_t'):
        the_vars = [
            'output_probejet_AK8_tau32',
            'output_probejet_AK8_maxDeepJet',
            'output_probejet_AK8_maxDeepCSV',
            'output_probejet_AK8_MDDeepAK8_TvsQCD',
            'output_probejet_AK8_mSD',
            'output_probejet_AK8_mass',
            'output_probejet_AK8_pt',
        ]
    elif the_tagger.name.startswith('ak8_w'):
        the_vars = [
            'output_probejet_AK8_ParticleNet_WvsQCD',
            'output_probejet_AK8_mSD',
            'output_probejet_AK8_mass',
            'output_probejet_AK8_pt',
        ]
    elif the_tagger.name.startswith('hotvr_t'):
        the_vars = [
            'output_probejet_HOTVR_tau32',
            'output_probejet_HOTVR_nsub',
            'output_probejet_HOTVR_fpt1',
            'output_probejet_HOTVR_mpair',
            'output_probejet_HOTVR_mass',
            'output_probejet_HOTVR_pt',
        ]
    # if args.vars:
    #     the_vars = args.vars

    # the_wps = the_tagger.wps
    # if args.wps:
    #     the_wps = [the_tagger.wps[int(x)] for x in args.wps]

    the_systs = args.systs

    for var in the_vars:
        for year in years:
            for syst in the_systs:
                create_input_hists(variable=var, tagger=the_tagger, arg_wp_index=-1, year=year, arg_syst=_SYSTEMATICS[syst], arg_fit_variable=False)


    # #_______________________________________
    # ## code part for condor jobs
    #
    # the_vars = [the_tagger.fit_variable]
    # if args.vars:
    #     the_vars = args.vars
    #
    # the_systs = args.systs
    #
    # for var in the_vars:
    #     for year in years:
    #
    #         # the_wps = the_tagger.wps
    #         the_wp_indices = range(0, len(the_tagger.get_wp(year=year)))
    #         if args.wps:
    #             the_wp_indices = [int(x) for x in args.wps]
    #
    #         for wp_index in the_wp_indices:
    #             for syst in the_systs:
    #                 create_input_hists(variable=var, tagger=the_tagger, arg_wp_index=wp_index, year=year, arg_syst=_SYSTEMATICS[syst])
