import os
import numpy as np
from tqdm import tqdm
import ROOT
from array import array
# import matplotlib
# matplotlib.use('Agg')
# from matplotlib import pyplot as plt
# from matplotlib import rc


outputDirPath = os.environ.get("CMSSW_BASE")+'/src/UHH2/LegacyTopTagging/output/WorkingPointStudy/UL17/'
fileName_prefix = 'uhh2.AnalysisModuleRunner.MC.'
fileName_postfix = '.npy'

fileName_qcd = 'QCD_HT300toInf_UL17'
filePath_qcd = outputDirPath+fileName_prefix+fileName_qcd+fileName_postfix
fileName_ttbar = 'TTbarToHadronic_UL17'
filePath_ttbar = outputDirPath+fileName_prefix+fileName_ttbar+fileName_postfix

workdir = outputDirPath+'workdir_npy/'

array_ttbar_raw = np.load(filePath_ttbar)
array_qcd_raw = np.load(filePath_qcd)
print 'Raw numpy files loaded into memory.'

# indices:
# 0: weight
# 1: pt
# 2: msd
# 3: subdeepcsv
# 4: tau32
# 5: dr (only for TTbar)

# define here the softdrop mass window
msd_min = 105.0
msd_max = 210.0

deepcsv = {
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

# define here which b-tagging selection you want to use
deepcsv_value = deepcsv['UL17']['loose']

# define here the jet pt intervals you want to analyze
pt_intervals = {
    '300toInf': {
        'pt_min': 300.,
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
    # '600to620': { # for quick test runs
    #     'pt_min': 600.,
    #     'pt_max': 620.,
    # },
}

recalculate = True


def get_eff_vs_tau32cut(ARG_array, ARG_sumofweights=1, ARG_nx=1000):
    weights = ARG_array[:,0].flatten()
    tau32 = ARG_array[:,4].flatten()
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
    tau32cuts_msd_btag, eff_qcd_msd_btag, eff_ttbar_msd_btag = None, None, None
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
        print '> mSD window + b-tag'
        nx_msd_btag = 1000
        print ' >> QCD'
        array_qcd_msd_btag = array_qcd_msd[array_qcd_msd[:,3] > deepcsv_value]
        tau32cuts_msd_btag, eff_qcd_msd_btag = get_eff_vs_tau32cut(array_qcd_msd_btag, total_sum_qcd, nx)
        np.save(workdir_pt+'tau32cuts_msd_btag.npy', tau32cuts_msd_btag)
        np.save(workdir_pt+'eff_qcd_msd_btag.npy', eff_qcd_msd_btag)
        print ' >> TTbar'
        array_ttbar_msd_btag = array_ttbar_msd[array_ttbar_msd[:,3] > deepcsv_value]
        tau32cuts_msd_btag, eff_ttbar_msd_btag = get_eff_vs_tau32cut(array_ttbar_msd_btag, total_sum_ttbar, nx)
        np.save(workdir_pt+'eff_ttbar_msd_btag.npy', eff_ttbar_msd_btag)
    else:
        print 'Using previously saved efficiencies'
        tau32cuts = np.load(workdir_pt+'tau32cuts.npy')
        eff_qcd = np.load(workdir_pt+'eff_qcd.npy')
        eff_ttbar = np.load(workdir_pt+'eff_ttbar.npy')
        tau32cuts_msd = np.load(workdir_pt+'tau32cuts_msd.npy')
        eff_qcd_msd = np.load(workdir_pt+'eff_qcd_msd.npy')
        eff_ttbar_msd = np.load(workdir_pt+'eff_ttbar_msd.npy')
        tau32cuts_msd_btag = np.load(workdir_pt+'tau32cuts_msd_btag.npy')
        eff_qcd_msd_btag = np.load(workdir_pt+'eff_qcd_msd_btag.npy')
        eff_ttbar_msd_btag = np.load(workdir_pt+'eff_ttbar_msd_btag.npy')

    print 'Re-conversion to ROOT'
    # rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
    # # rc('text', usetex=True)
    # fig, ax = plt.subplots()
    # ax.plot(tau32cuts, eff_qcd, color='royalblue', linestyle='solid', linewidth=1)
    # ax.plot(tau32cuts, eff_qcd_msd, color='darkorange', linestyle='dashed', linewidth=1)
    # ax.plot(tau32cuts, eff_qcd_msd_btag, color='magenta', linestyle='dashdot', linewidth=1)
    # fig.gca().set_xlabel(r'$\tau_3/\tau_2$ upper limit')
    # fig.gca().set_ylabel(r'$\epsilon_\mathrm{B}$')
    # ax.set_yscale('log')
    # fig.savefig(workdir_pt+'plot_'+pt+'_eff_qcd.pdf')

    root_file = ROOT.TFile.Open(workdir_pt+'root_Pt'+pt+'.root', 'RECREATE')
    root_file.cd()

    # Do some really cheesy stuff to convert back to ROOT ...

    h_tau32cuts = ROOT.TH1F("tau32cuts", "tau32cuts", len(tau32cuts), 0, 1)
    for i in range(len(tau32cuts)): h_tau32cuts.SetBinContent(i, tau32cuts[i])
    h_tau32cuts.Write()

    h_eff_qcd = ROOT.TH1F("eff_qcd", "eff_qcd", len(eff_qcd), 0, 1)
    for i in range(len(eff_qcd)): h_eff_qcd.SetBinContent(i, eff_qcd[i])
    h_eff_qcd.Write()

    h_eff_ttbar = ROOT.TH1F("eff_ttbar", "eff_ttbar", len(eff_ttbar), 0, 1)
    for i in range(len(eff_ttbar)): h_eff_ttbar.SetBinContent(i, eff_ttbar[i])
    h_eff_ttbar.Write()


    h_tau32cuts_msd = ROOT.TH1F("tau32cuts_msd", "tau32cuts_msd", len(tau32cuts_msd), 0, 1)
    for i in range(len(tau32cuts_msd)): h_tau32cuts_msd.SetBinContent(i, tau32cuts_msd[i])
    h_tau32cuts_msd.Write()

    h_eff_qcd_msd = ROOT.TH1F("eff_qcd_msd", "eff_qcd_msd", len(eff_qcd_msd), 0, 1)
    for i in range(len(eff_qcd_msd)): h_eff_qcd_msd.SetBinContent(i, eff_qcd_msd[i])
    h_eff_qcd_msd.Write()

    h_eff_ttbar_msd = ROOT.TH1F("eff_ttbar_msd", "eff_ttbar_msd", len(eff_ttbar_msd), 0, 1)
    for i in range(len(eff_ttbar_msd)): h_eff_ttbar_msd.SetBinContent(i, eff_ttbar_msd[i])
    h_eff_ttbar_msd.Write()


    h_tau32cuts_msd_btag = ROOT.TH1F("tau32cuts_msd_btag", "tau32cuts_msd_btag", len(tau32cuts_msd_btag), 0, 1)
    for i in range(len(tau32cuts_msd_btag)): h_tau32cuts_msd_btag.SetBinContent(i, tau32cuts_msd_btag[i])
    h_tau32cuts_msd_btag.Write()

    h_eff_qcd_msd_btag = ROOT.TH1F("eff_qcd_msd_btag", "eff_qcd_msd_btag", len(eff_qcd_msd_btag), 0, 1)
    for i in range(len(eff_qcd_msd_btag)): h_eff_qcd_msd_btag.SetBinContent(i, eff_qcd_msd_btag[i])
    h_eff_qcd_msd_btag.Write()

    h_eff_ttbar_msd_btag = ROOT.TH1F("eff_ttbar_msd_btag", "eff_ttbar_msd_btag", len(eff_ttbar_msd_btag), 0, 1)
    for i in range(len(eff_ttbar_msd_btag)): h_eff_ttbar_msd_btag.SetBinContent(i, eff_ttbar_msd_btag[i])
    h_eff_ttbar_msd_btag.Write()

    root_file.Close()
