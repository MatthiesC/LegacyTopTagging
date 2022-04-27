import os
import sys
import numpy as np
from tqdm import tqdm
import ROOT
from array import array
import argparse
# import matplotlib
# matplotlib.use('Agg')
# from matplotlib import pyplot as plt
# from matplotlib import rc

# true: re-analyze efficiency vs. tau32 cuts and store results in numpy files
# false: only load previously saved numpy files and store results in ROOT format
options_recalculate = ['True', 'False']
options_year = ['UL16preVFP', 'UL16postVFP', 'UL17', 'UL18']

if not sys.argv[1:]: sys.exit('No arguments provided. Exit.')
parser = argparse.ArgumentParser()
parser.add_argument('-y', '--year', type=str, choices=options_year)
parser.add_argument('-r', '--recalculate', action='store_true', help='Recalculate the efficiencies and do not use already existing numpy outfiles for the TGraphs')
parser.add_argument('-m', '--massCut', type=float, help='Entirely ignore jets with a SoftDrop mass smaller than the given value.')
args = parser.parse_args(sys.argv[1:])

recalculate = args.recalculate
year = args.year
massCut_given = type(args.massCut) == float
msd_threshold = None
if massCut_given:
    msd_threshold = float(args.massCut)

outputDirPath = os.environ.get("CMSSW_BASE")+'/src/UHH2/LegacyTopTagging/output/WorkingPointStudy/'+year+'/'
fileName_prefix = 'uhh2.AnalysisModuleRunner.MC.'
fileName_postfix = '.npy'

fileName_qcd = 'QCD_HT300toInf_'+year
filePath_qcd = outputDirPath+fileName_prefix+fileName_qcd+fileName_postfix
fileName_ttbar = 'TTbarToHadronic_'+year
filePath_ttbar = outputDirPath+fileName_prefix+fileName_ttbar+fileName_postfix

workdir = outputDirPath+'workdir_npy/'

array_ttbar_raw = None
array_qcd_raw = None
if recalculate:
    array_ttbar_raw = np.load(filePath_ttbar)
    array_qcd_raw = np.load(filePath_qcd)
    print 'Raw numpy files loaded into memory.'
    if massCut_given:
        print 'Rejecting jets with SoftDrop mass smaller than', msd_threshold
        n_ttbar_before = len(array_ttbar_raw)
        array_ttbar_raw = array_ttbar_raw[array_ttbar_raw[:,2] > msd_threshold]
        print 'Percentage of top jets fulfilling the mass threshold: %f' % (100.*len(array_ttbar_raw)/n_ttbar_before)
        n_qcd_before = len(array_qcd_raw)
        array_qcd_raw = array_qcd_raw[array_qcd_raw[:,2] > msd_threshold]
        print 'Percentage of QCD jets fulfilling the mass threshold: %f' % (100.*len(array_qcd_raw)/n_qcd_before)

# indices:
# 0: weight
# 1: pt
# 2: msd
# 3: subdeepcsv
# 4: subdeepjet
# 5: tau32
# 6: dr (only for TTbar)

# define here the softdrop mass window
msd_min = 105.0
msd_max = 210.0

deepcsv = {
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

deepjet = {
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


# define here which b-tagging selection you want to use
deepcsv_value = deepcsv[year]['loose']
deepjet_value = deepjet[year]['loose']

# define here the jet pt intervals you want to analyze
pt_intervals = {
    '300toInf': {
        'pt_min': 300.,
        'pt_max': np.inf,
    },
    '400toInf': {
        'pt_min': 400.,
        'pt_max': np.inf,
    },
    '300to400': {
        'pt_min': 300.,
        'pt_max': 400.,
    },
    '400to480': {
        'pt_min': 400.,
        'pt_max': 480.,
    },
    '480to600': {
        'pt_min': 480.,
        'pt_max': 600.,
    },
    '600toInf': {
        'pt_min': 600.,
        'pt_max': np.inf,
    },
    '1000toInf': {
    'pt_min': 1000.,
    'pt_max': np.inf,
    },
    # '600to620': { # for quick test runs
    #     'pt_min': 600.,
    #     'pt_max': 620.,
    # },
}


def get_eff_vs_tau32cut(ARG_array, ARG_sumofweights=1, ARG_nx=1000):
    weights = ARG_array[:,0].flatten()
    tau32 = ARG_array[:,5].flatten()
    x = np.array([1.], dtype='f') # tau32 cut value
    y = np.array([weights.sum()], dtype='f') # efficiency
    for i in tqdm(reversed(range(ARG_nx)), desc='Calculating efficiency vs. N-subjettiness cut', total=ARG_nx, dynamic_ncols=True, leave=False):
        tau32cut = float(i)/ARG_nx
        x = np.append(x, [tau32cut])
        weights = weights[tau32 < tau32cut]
        tau32 = tau32[tau32 < tau32cut]
        new_y = weights.sum()
        y = np.append(y, [new_y])
    # revert x and y
    x = x[::-1]
    y = y[::-1]
    return x, y/ARG_sumofweights


for pt in pt_intervals.keys():
    print 'Working on jet pt =', pt
    workdir_pt = workdir+'Pt'+pt+'/'
    os.system('mkdir -p '+workdir_pt)
    tau32cuts, eff_qcd, eff_ttbar = None, None, None
    tau32cuts_msd, eff_qcd_msd, eff_ttbar_msd = None, None, None
    tau32cuts_msd_deepcsv, eff_qcd_msd_deepcsv, eff_ttbar_msd_deepcsv = None, None, None
    tau32cuts_msd_deepjet, eff_qcd_msd_deepjet, eff_ttbar_msd_deepjet = None, None, None
    if recalculate:
        print '> inclusive'
        nx = 1000
        print ' >> QCD'
        array_qcd = array_qcd_raw[(array_qcd_raw[:,1] > pt_intervals[pt]['pt_min']) & (array_qcd_raw[:,1] < pt_intervals[pt]['pt_max'])]
        total_sum_qcd = array_qcd[:,0].sum()
        tau32cuts, eff_qcd = get_eff_vs_tau32cut(array_qcd, total_sum_qcd, nx)
        np.save(workdir_pt+'tau32cuts.npy', tau32cuts)
        np.save(workdir_pt+'eff_qcd.npy', eff_qcd)
        print ' >> TTbar'
        array_ttbar = array_ttbar_raw[(array_ttbar_raw[:,1] > pt_intervals[pt]['pt_min']) & (array_ttbar_raw[:,1] < pt_intervals[pt]['pt_max'])]
        total_sum_ttbar = array_ttbar[:,0].sum()
        tau32cuts, eff_ttbar = get_eff_vs_tau32cut(array_ttbar, total_sum_ttbar, nx)
        np.save(workdir_pt+'eff_ttbar.npy', eff_ttbar)
        print '> mSD window'
        nx_msd = 1000
        print ' >> QCD'
        array_qcd_msd = array_qcd[(array_qcd[:,2] > msd_min) & (array_qcd[:,2] < msd_max)]
        tau32cuts_msd, eff_qcd_msd = get_eff_vs_tau32cut(array_qcd_msd, total_sum_qcd, nx)
        np.save(workdir_pt+'tau32cuts_msd.npy', tau32cuts_msd)
        np.save(workdir_pt+'eff_qcd_msd.npy', eff_qcd_msd)
        print ' >> TTbar'
        array_ttbar_msd = array_ttbar[(array_ttbar[:,2] > msd_min) & (array_ttbar[:,2] < msd_max)]
        tau32cuts_msd, eff_ttbar_msd = get_eff_vs_tau32cut(array_ttbar_msd, total_sum_ttbar, nx)
        np.save(workdir_pt+'eff_ttbar_msd.npy', eff_ttbar_msd)
        print '> mSD window + DeepCSV b-tag'
        nx_msd_deepcsv = 1000
        print ' >> QCD'
        array_qcd_msd_deepcsv = array_qcd_msd[array_qcd_msd[:,3] > deepcsv_value]
        tau32cuts_msd_deepcsv, eff_qcd_msd_deepcsv = get_eff_vs_tau32cut(array_qcd_msd_deepcsv, total_sum_qcd, nx)
        np.save(workdir_pt+'tau32cuts_msd_deepcsv.npy', tau32cuts_msd_deepcsv)
        np.save(workdir_pt+'eff_qcd_msd_deepcsv.npy', eff_qcd_msd_deepcsv)
        print ' >> TTbar'
        array_ttbar_msd_deepcsv = array_ttbar_msd[array_ttbar_msd[:,3] > deepcsv_value]
        tau32cuts_msd_deepcsv, eff_ttbar_msd_deepcsv = get_eff_vs_tau32cut(array_ttbar_msd_deepcsv, total_sum_ttbar, nx)
        np.save(workdir_pt+'eff_ttbar_msd_deepcsv.npy', eff_ttbar_msd_deepcsv)
        print '> mSD window + DeepJet b-tag'
        nx_msd_deepjet = 1000
        print ' >> QCD'
        array_qcd_msd_deepjet = array_qcd_msd[array_qcd_msd[:,4] > deepjet_value]
        tau32cuts_msd_deepjet, eff_qcd_msd_deepjet = get_eff_vs_tau32cut(array_qcd_msd_deepjet, total_sum_qcd, nx)
        np.save(workdir_pt+'tau32cuts_msd_deepjet.npy', tau32cuts_msd_deepjet)
        np.save(workdir_pt+'eff_qcd_msd_deepjet.npy', eff_qcd_msd_deepjet)
        print ' >> TTbar'
        array_ttbar_msd_deepjet = array_ttbar_msd[array_ttbar_msd[:,4] > deepjet_value]
        tau32cuts_msd_deepjet, eff_ttbar_msd_deepjet = get_eff_vs_tau32cut(array_ttbar_msd_deepjet, total_sum_ttbar, nx)
        np.save(workdir_pt+'eff_ttbar_msd_deepjet.npy', eff_ttbar_msd_deepjet)
    else:
        print 'Using previously saved efficiencies'
        tau32cuts = np.load(workdir_pt+'tau32cuts.npy')
        eff_qcd = np.load(workdir_pt+'eff_qcd.npy')
        eff_ttbar = np.load(workdir_pt+'eff_ttbar.npy')
        tau32cuts_msd = np.load(workdir_pt+'tau32cuts_msd.npy')
        eff_qcd_msd = np.load(workdir_pt+'eff_qcd_msd.npy')
        eff_ttbar_msd = np.load(workdir_pt+'eff_ttbar_msd.npy')
        tau32cuts_msd_deepcsv = np.load(workdir_pt+'tau32cuts_msd_deepcsv.npy')
        eff_qcd_msd_deepcsv = np.load(workdir_pt+'eff_qcd_msd_deepcsv.npy')
        eff_ttbar_msd_deepcsv = np.load(workdir_pt+'eff_ttbar_msd_deepcsv.npy')
        tau32cuts_msd_deepjet = np.load(workdir_pt+'tau32cuts_msd_deepjet.npy')
        eff_qcd_msd_deepjet = np.load(workdir_pt+'eff_qcd_msd_deepjet.npy')
        eff_ttbar_msd_deepjet = np.load(workdir_pt+'eff_ttbar_msd_deepjet.npy')

        print 'Conversion to ROOT TGraph objects'
        root_file = ROOT.TFile.Open(workdir_pt+'root_Pt'+pt+'.root', 'RECREATE')
        root_file.cd()

        graph = ROOT.TGraph(len(eff_qcd), tau32cuts, eff_qcd)
        graph.Write("eff_qcd")
        graph = ROOT.TGraph(len(eff_ttbar), tau32cuts, eff_ttbar)
        graph.Write("eff_ttbar")
        graph = ROOT.TGraph(len(eff_qcd), eff_ttbar, eff_qcd)
        graph.Write("roc")

        graph = ROOT.TGraph(len(eff_qcd_msd), tau32cuts_msd, eff_qcd_msd)
        graph.Write("eff_qcd_msd")
        graph = ROOT.TGraph(len(eff_ttbar_msd), tau32cuts_msd, eff_ttbar_msd)
        graph.Write("eff_ttbar_msd")
        graph = ROOT.TGraph(len(eff_qcd_msd), eff_ttbar_msd, eff_qcd_msd)
        graph.Write("roc_msd")

        graph = ROOT.TGraph(len(eff_qcd_msd_deepcsv), tau32cuts_msd_deepcsv, eff_qcd_msd_deepcsv)
        graph.Write("eff_qcd_msd_deepcsv")
        graph = ROOT.TGraph(len(eff_ttbar_msd_deepcsv), tau32cuts_msd_deepcsv, eff_ttbar_msd_deepcsv)
        graph.Write("eff_ttbar_msd_deepcsv")
        graph = ROOT.TGraph(len(eff_qcd_msd_deepcsv), eff_ttbar_msd_deepcsv, eff_qcd_msd_deepcsv)
        graph.Write("roc_msd_deepcsv")

        graph = ROOT.TGraph(len(eff_qcd_msd_deepjet), tau32cuts_msd_deepjet, eff_qcd_msd_deepjet)
        graph.Write("eff_qcd_msd_deepjet")
        graph = ROOT.TGraph(len(eff_ttbar_msd_deepjet), tau32cuts_msd_deepjet, eff_ttbar_msd_deepjet)
        graph.Write("eff_ttbar_msd_deepjet")
        graph = ROOT.TGraph(len(eff_qcd_msd_deepjet), eff_ttbar_msd_deepjet, eff_qcd_msd_deepjet)
        graph.Write("roc_msd_deepjet")

        root_file.Close()
