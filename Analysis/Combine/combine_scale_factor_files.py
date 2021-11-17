import os
import sys
import ROOT as root
from collections import OrderedDict
import numpy as np

categories = OrderedDict([ # name as used in this framework: name as it should be in the final root file
    ('FullyMerged', 'FullyMerged'),
    ('YllufMerged', 'NotMerged'),
])

wps = {
    'AK8': {
        'BkgEff0p001': 'wp0p38_vt',
        'BkgEff0p005': 'wp0p47_t',
        'BkgEff0p010': 'wp0p52_m',
        'BkgEff0p025': 'wp0p61_l',
        'BkgEff0p050': 'wp0p69_vl',
    },
    'HOTVR': {
        'Standard': 'default',
    },
}

wp_to_MisTagEffString = {
    'BkgEff0p001': 'mis0p001',
    'BkgEff0p005': 'mis0p005',
    'BkgEff0p010': 'mis0p010',
    'BkgEff0p025': 'mis0p025',
    'BkgEff0p050': 'mis0p050',
}

def take_Yerrors_from_second(graph1, graph2):
    '''
    :param graph1/2: TGraphAsymmErrors
    Keep the x values and x erros from first.
    '''
    # sanity check:
    if graph1.GetN() != graph2.GetN():
        sys.exit('The two graphs do not have same number of points')
    x = []
    y = []
    exl = []
    exh = []
    eyl = []
    eyh = []
    for ibin in range(graph1.GetN()):
        x_ = root.Double()
        y_ = root.Double()
        graph1.GetPoint(ibin, x_, y_)
        x.append(float(x_))
        y.append(float(y_))
        exl.append(graph1.GetErrorXlow(ibin))
        exh.append(graph1.GetErrorXhigh(ibin))
        eyl.append(graph2.GetErrorYlow(ibin))
        eyh.append(graph2.GetErrorYhigh(ibin))
    new_graph = root.TGraphAsymmErrors(len(x), np.array(x), np.array(y), np.array(exl), np.array(exh), np.array(eyl), np.array(eyh))
    return new_graph

def combine_scale_factor_files(tasks):
    '''
    :param tasks: list of entities of ComplexCombineTask
    '''
    # sanity check:
    year = tasks[0].Year
    for task in tasks:
        if task.Year != year:
            sys.exit('Not all tasks have the same year!!! Abort.')

    outfile_name = 'TopTaggingScaleFactors_RunIISummer19'+year['short_name']+'_PUPPIv15.root'
    workdir = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/Analysis/Combine/workdir')
    outfile_path = os.path.join(workdir, outfile_name)
    outfile = root.TFile.Open(outfile_path, 'RECREATE')
    outfile.cd()

    for task in tasks:
        infile = root.TFile.Open(task.sf_file_path_obs, 'READ')
        sf_graphs = OrderedDict()
        folder_name = []
        folder_name.append(task.Algo)
        folder_name.append('PUPPI')
        if task.JetVersion == 'BTag':
            folder_name.append('subjetbtag')
        folder_name.append(wps[task.Algo][task.WP])
        if task.Algo == 'AK8' and task.JetVersion == 'All':
            folder_name.append(wp_to_MisTagEffString[task.WP])
        folder_name = '_'.join(folder_name)
        folder = outfile.mkdir(folder_name)
        folder.cd()
        for cat, translated_cat in categories.items():
            sf_graphs[cat] = OrderedDict()
            sf_graphs[cat]['tot'] = infile.Get(cat+'_tot')
            sf_graphs[cat]['stat'] = take_Yerrors_from_second(sf_graphs[cat]['tot'], infile.Get(cat+'_stat'))
            sf_graphs[cat]['syst'] = take_Yerrors_from_second(sf_graphs[cat]['tot'], infile.Get(cat+'_syst'))
            for unc in ['tot', 'stat', 'syst']:
                sf_graphs[cat][unc].Write(translated_cat+'_'+unc)
