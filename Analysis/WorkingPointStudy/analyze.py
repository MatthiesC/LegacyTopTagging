from __future__ import print_function

import os
import sys
import numpy as np
from tqdm import tqdm
# from array import array
import argparse
# import matplotlib
# matplotlib.use('Agg')
# from matplotlib import pyplot as plt
# from matplotlib import rc

# import uproot
# from ROOT import TFile, TGraph

sys.path.append(os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/Analysis'))
from constants import _TAGGERS, VarInterval
from parallel_threading import run_with_pool

options_year = ['UL16preVFP', 'UL16postVFP', 'UL17', 'UL18']

if not sys.argv[1:]: sys.exit('No arguments provided. Exit.')
parser = argparse.ArgumentParser()
parser.add_argument('-y', '--year', type=str, choices=options_year)
parser.add_argument('-r', '--recalculate', action='store_true', help='Recalculate the efficiencies and do not use already existing numpy outfiles for the TGraphs')
parser.add_argument('-p', '--pool', action='store_true', help='Use parallel threading (each analyzed tagger gets one thread)')
args = parser.parse_args(sys.argv[1:])

recalculate = args.recalculate
year = args.year
if args.pool:
    print('Running in Pool mode')

if recalculate:
    import uproot
else:
    from ROOT import TGraph, TFile

inputDir = os.environ.get("CMSSW_BASE")+'/src/UHH2/LegacyTopTagging/output/WorkingPointStudy/'+year+'/nominal/'
file_name_prefix = 'uhh2.AnalysisModuleRunner.MC.'

scan_vars = {
    'tau32': VarInterval('tau32', 0, 1, inversed=True),
    'tau21': VarInterval('tau21', 0, 1, inversed=True),
    'partnet_TvsQCD': VarInterval('partnet_TvsQCD', 0, 1),
    'partnet_WvsQCD': VarInterval('partnet_WvsQCD', 0, 1),
    'deepak8_TvsQCD': VarInterval('deepak8_TvsQCD', 0, 1),
    'deepak8_WvsQCD': VarInterval('deepak8_WvsQCD', 0, 1),
    'MDdeepak8_TvsQCD': VarInterval('MDdeepak8_TvsQCD', 0, 1),
    'MDdeepak8_WvsQCD': VarInterval('MDdeepak8_WvsQCD', 0, 1),
}

def get_file_names(key, year, sig_wjets = False):
    file_name_bkg, file_name_sig = file_name_prefix, file_name_prefix
    if sig_wjets:
        if key.startswith('ak8_w'):
            file_name_sig += 'WJetsToQQ_HT200toInf_'+year+'.root.restructured.AK8'
        else:
            sys.exit('Requesting WJetsToQQ sample for signal but this tagger is not a W tagger')
    else:
        if key.startswith('ak8_w'):
            file_name_sig += 'TTbarToHadronic_'+year+'.root.restructured_w.AK8'
        elif key.startswith('ak8_t'):
            file_name_sig += 'TTbarToHadronic_'+year+'.root.restructured_t.AK8'
        elif key.startswith('hotvr_t'):
            file_name_sig += 'TTbarToHadronic_'+year+'.root.restructured_t.HOTVR'
        else:
            sys.exit('Don\'t know what signal sample to use')
    if key.startswith('ak8'):
        file_name_bkg += 'QCD_HT200toInf_'+year+'.root.restructured.AK8'
    elif key.startswith('hotvr'):
        file_name_bkg += 'QCD_HT200toInf_'+year+'.root.restructured.HOTVR'
    else:
        sys.exit('Don\'t know what background sample to use')
    return file_name_bkg, file_name_sig

def get_matching_rule(key):
    if key.startswith('ak8'):
        return 'dr < 0.6' # using 0.6 instead of 0.8 here deliberately; 0.6 is more in line with what has been done in the Heavy Res Paper
    elif key.startswith('hotvr'):
        return 'dr < reff' # Heavy Res Paper actually also used a fixed radius for HOTVR; but using Reff here is more in line with what I did for the HOTVR jet pt response study
    else:
        sys.exit('Cannot determine matching rule based on tagger name')

def scan_eff_vs_variable(arr_weight_, arr_variable_, scan_var, norm=None, steps=1000, do_pbar=True):
    arr_weight = arr_weight_
    arr_variable = arr_variable_
    x = np.array([scan_var.var_max if scan_var.inversed else scan_var.var_min], dtype='f')
    y = np.array([np.sum(arr_weight)], dtype='f')
    the_range = reversed(range(steps)) if scan_var.inversed else range(1, steps+1)
    for i in (tqdm(the_range, desc='Scan eff. vs. {0}'.format(scan_var.var), total=steps, dynamic_ncols=True, leave=False) if do_pbar else the_range):
        cutoff = float(i)*(scan_var.var_max - scan_var.var_min)/steps + scan_var.var_min
        x = np.append(x, [cutoff])
        cutoff_rule = arr_variable < cutoff if scan_var.inversed else arr_variable > cutoff
        arr_weight = arr_weight[cutoff_rule]
        arr_variable = arr_variable[cutoff_rule]
        y = np.append(y, [np.sum(arr_weight)])
    # revert x and y again since we looped in reversed order of scan_var range
    x = x[::-1]
    y = y[::-1]
    if norm:
        return x, y/norm # norm should be sum of weights of empty tag_rule
    else:
        return x, y

# for tagger_k, tagger_v in _TAGGERS.items():
def analyze_tagger(tagger_k, pool=False):
    tagger_v = _TAGGERS[tagger_k]
    if not pool: print('Tagger:', tagger_k)
    if recalculate:
        file_path_bkg, file_path_sig = get_file_names(tagger_k, year)
        tree_bkg = uproot.open(os.path.join(inputDir, file_path_bkg))['my_tree']
        tree_sig = uproot.open(os.path.join(inputDir, file_path_sig))['my_tree']
        scan_var = scan_vars[tagger_v.scan_var]
    for i_var, var in enumerate(tagger_v.var_intervals.values()):
        print('{0}Interval:'.format('Tagger: {0}, '.format(tagger_k) if pool else ''), var.name, '({0}/{1})'.format(i_var+1, len(tagger_v.var_intervals)))
        outputDir = os.path.join(inputDir, tagger_k, var.name)
        # os.makedirs(outputDir, exist_ok=True)
        os.system('mkdir -p '+outputDir)
        if recalculate:
            if not pool: print('(Re-)calculate efficiencies')
            if not pool: print('Background')
            expr = '({tag_rule}) & ({var_rule})'.format(tag_rule=tagger_v.get_tag_rule(year=year), var_rule='({0} > {1}) & ({0} < {2})'.format(var.var, var.var_min, var.var_max))
            arrays_bkg = tree_bkg.arrays(['weight', scan_var.var], expr)
            expr_norm = '({var_rule})'.format(var_rule='({0} > {1}) & ({0} < {2})'.format(var.var, var.var_min, var.var_max))
            norm = np.sum(tree_bkg.arrays(['weight'], expr_norm)['weight'])
            var_cuts, eff_bkg = scan_eff_vs_variable(arrays_bkg['weight'], arrays_bkg[scan_var.var], scan_var, norm=norm, do_pbar=(not pool))
            np.save(os.path.join(outputDir, 'eff_bkg.npy'), eff_bkg)
            if not pool: print('Signal')
            expr = '(({expr}) & ({matching_rule}))'.format(expr=expr, matching_rule=get_matching_rule(tagger_k))
            arrays_sig = tree_sig.arrays(['weight', scan_var.var], expr)
            expr_norm = '(({expr_norm}) & ({matching_rule}))'.format(expr_norm=expr_norm, matching_rule=get_matching_rule(tagger_k))
            norm = np.sum(tree_sig.arrays(['weight'], expr_norm)['weight'])
            var_cuts, eff_sig = scan_eff_vs_variable(arrays_sig['weight'], arrays_sig[scan_var.var], scan_var, norm=norm, do_pbar=(not pool))
            np.save(os.path.join(outputDir, 'eff_sig.npy'), eff_sig)
            np.save(os.path.join(outputDir, 'var_cuts.npy'), var_cuts)
        else:
            if not pool: print('Using previously saved efficiencies')
            eff_bkg = np.load(os.path.join(outputDir, 'eff_bkg.npy'))
            eff_bkg = eff_bkg.astype(float)
            eff_sig = np.load(os.path.join(outputDir, 'eff_sig.npy'))
            eff_sig = eff_sig.astype(float)
            var_cuts = np.load(os.path.join(outputDir, 'var_cuts.npy'))
            var_cuts = var_cuts.astype(float)

            if not pool: print('Conversion to ROOT TGraph objects')
            root_file_path = os.path.join(outputDir, 'graphs.root')
            root_file = TFile.Open(root_file_path, 'RECREATE')
            root_file.cd()
            graph = TGraph(len(var_cuts), var_cuts, eff_bkg)
            graph.Write('eff_bkg')
            graph = TGraph(len(var_cuts), var_cuts, eff_sig)
            graph.Write('eff_sig')
            graph = TGraph(len(var_cuts), eff_sig, eff_bkg)
            graph.Write('roc')
            root_file.Close()

if __name__=='__main__':
    # _TAGGERS_to_analyze = _TAGGERS.keys()
    _TAGGERS_to_analyze = [
        # 'ak8_t_incl__MDdeepak8',
        # 'ak8_t__MDdeepak8',
        # 'ak8_w_incl__MDdeepak8',
        # 'ak8_w__MDdeepak8',
        # 'ak8_t_incl__deepak8',
        # 'ak8_t__deepak8',
        # 'ak8_w_incl__deepak8',
        # 'ak8_w__deepak8',
        # 'ak8_t_incl__partnet',
        # 'ak8_t__partnet',
        # 'ak8_w_incl__partnet',
        # 'ak8_w__partnet',
        # 'hotvr_t__tau',
        # 'hotvr_t_incl__tau',
        # 'ak8_t_incl__tau',
        # 'ak8_t__tau',
        # 'ak8_w_incl__tau',
        # 'ak8_w__tau',


        'ak8_t_btagDCSV__tau',
        'ak8_t_btagDJet__tau',
    ]
    sets_of_args = []
    for tagger in _TAGGERS_to_analyze:
        set_of_args = {}
        set_of_args['tagger_k'] = tagger
        set_of_args['pool'] = args.pool # in pool mode, we do not want the individual pbars that would spam the terminal endlessly
        sets_of_args.append(set_of_args)
        if not args.pool: analyze_tagger(**set_of_args)
    if args.pool: run_with_pool(function=analyze_tagger, sets_of_args=sets_of_args)
