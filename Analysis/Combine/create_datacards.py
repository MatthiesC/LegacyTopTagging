#!/usr/bin/env python

import os
import sys
import argparse
from collections import OrderedDict
import itertools
import ROOT as root
import subprocess
from timeit import default_timer as timer
import numpy as np
import warnings
import re

from constants import years

# https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiLUM
# years = {
#     'UL17': {
#         'lumi': 41.53,
#         'lumi_unc': 0.023,
#     },
#     'UL18': {
#         'lumi': 59.74,
#         'lumi_unc': 0.025,
#     },
# }

parser = argparse.ArgumentParser()
parser.add_argument('-y', '--year', choices=years.keys(), nargs=1)
parser.add_argument('-p', '--reduced-pt-bins', action='store_true')
parser.add_argument('-i', '--run-impacts', action='store_true')
args = parser.parse_args(sys.argv[1:])

arg_year = years.get(args.year[0])

pt_intervals = {
    '200to250': {
        'pt_min': 200.,
        'pt_max': 250.,
    },
    '250to300': {
        'pt_min': 250.,
        'pt_max': 300.,
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
        # 'pt_max': np.inf,
        'pt_max': 1000.,
    },
    '400toInf': {
        'pt_min': 400.,
        # 'pt_max': np.inf,
        'pt_max': 1000.,
    },
}

if args.reduced_pt_bins:
    pt_bins = {
        'AK8': [
            'Pt300to400',
            'Pt400toInf',
        ],
        'HOTVR': [
            'Pt200to250',
            'Pt250to300',
            'Pt300to400',
            # 'Pt400toInf',
            'Pt400to480',
            'Pt480to600',
            'Pt600toInf',
        ],
    }
else:
    pt_bins = {
        'AK8': [
            'Pt300to400',
            'Pt400to480',
            'Pt480to600',
            'Pt600toInf',
        ],
        'HOTVR': [
            'Pt200to250',
            'Pt250to300',
            'Pt300to400',
            'Pt400to480',
            'Pt480to600',
            'Pt600toInf',
        ],
    }

jet_versions = {
    'AK8': [
        'All',
        # 'Mass',
        'BTag',
        # 'MassAndBTag',
    ],
    'HOTVR': [
        # 'All',
        # 'Mass',
        'HOTVRCuts',
        # 'HOTVRCutsAndMass',
    ],
}

wps = {
    'AK8': [
        'BkgEff0p001',
        'BkgEff0p005',
        'BkgEff0p010',
        'BkgEff0p025',
        'BkgEff0p050',
    ],
    'HOTVR': [
        'Standard',
    ],
}

regions = [
    'Pass', # passed tau32
    'Fail', # failed tau32
    # 'PassW', # failed tau32, passed tau21
    # 'FailW',  # failed tau32, failed tau21
]

base_processes = OrderedDict([ # processes not divided into merge categories, prior xsec uncertainty in percent
    ('TTbar', 0.05),
    ('ST', 0.20),
    ('WJets', 0.20),
    ('DYJetsAndDiboson', 0.50),
    ('QCD', 1.00),
])

processes = OrderedDict([ # ids <= 0 represent signals
    # ('TTbar', {
    #     'id': -100,
    #     'base_process': 'TTbar',
    # }),
    ('TTbar_FullyMerged', {
        'id': -101,
        'base_process': 'TTbar',
    }),
    # ('TTbar_SemiMerged', {
    #     'id': -102,
    #     'base_process': 'TTbar',
    # }),
    # ('TTbar_WMerged', {
    #     'id': -103,
    #     'base_process': 'TTbar',
    # }),
    # ('TTbar_QBMerged', {
    #     'id': -104,
    #     'base_process': 'TTbar',
    # }),
    # ('TTbar_BkgOrNotMerged', {
    #     'id': -105,
    #     'base_process': 'TTbar',
    # }),
    # ('TTbar_NotFullyOrWMerged', {
    #     'id': -106,
    #     'base_process': 'TTbar',
    # }),
    ('TTbar_YllufMerged', {
        'id': -107,
        'base_process': 'TTbar',
    }),
    ('TTbar_Background', {
        'id': 108,
        'base_process': 'TTbar',
    }),
    # ('ST', {
    #     'id': -200,
    #     'base_process': 'ST',
    # }),
    ('ST_FullyMerged', {
        'id': -201,
        'base_process': 'ST',
    }),
    # ('ST_SemiMerged', {
    #     'id': -202,
    #     'base_process': 'ST',
    # }),
    # ('ST_WMerged', {
    #     'id': -203,
    #     'base_process': 'ST',
    # }),
    # ('ST_QBMerged', {
    #     'id': -204,
    #     'base_process': 'ST',
    # }),
    # ('ST_BkgOrNotMerged', {
    #     'id': -205,
    #     'base_process': 'ST',
    # }),
    # ('ST_NotFullyOrWMerged', {
    #     'id': -206,
    #     'base_process': 'ST',
    # }),
    ('ST_YllufMerged', {
        'id': -207,
        'base_process': 'ST',
    }),
    ('ST_Background', {
        'id': 208,
        'base_process': 'ST',
    }),
    ('WJetsToLNu', {
        'id': 300,
        'base_process': 'WJets',
    }),
    ('DYJetsToLLAndDiboson', {
        'id': 400,
        'base_process': 'DYJetsAndDiboson',
    }),
    ('QCD_Mu', {
        'id': 500,
        'base_process': 'QCD',
    }),
])

systematics = [
    'btaggingbc',
    'btaggingudsg',
    'jec',
    'jer',
    'scale',
    'pileup',
    'fsr',
    # 'tageff3prong',
    # 'tageff2prong',
    # 'tageff1prong',
    # 'tageffFully',
    # 'tageffYlluf',
    # 'tageffBkgrd',
    # 'tau21',
    'topptA',
    'topptB',
]

global_categories = [
    'FullyMerged',
    # 'SemiMerged',
    # 'WMerged',
    # 'BkgOrNotMerged',
    # 'NotFullyOrWMerged',
    'YllufMerged',
]

def expobs(observed, short=True):
    if short:
        return ('obs' if observed else 'exp')
    else:
        return ('observed' if observed else 'expected')

def get_binning(th1):
    binning = np.array([th1.GetXaxis().GetBinLowEdge(1)], dtype=float)
    for ibin in range(1, th1.GetNbinsX()+1):
        binning = np.append(binning, th1.GetXaxis().GetBinUpEdge(ibin))
    return binning

def return_rebinned(th1, binning):
    if th1.GetNbinsX() != len(binning)-1:
        warnings.warn('Cannot rebin TH1: given binning scheme does not match number of bins of TH1')
        return th1
    th1_new = root.TH1F(th1.GetName(), th1.GetTitle(), len(binning)-1, binning)
    for ibin in range(len(th1)+2):
        th1_new.SetBinContent(ibin, th1.GetBinContent(ibin))
        th1_new.SetBinError(ibin, th1.GetBinError(ibin))
    return th1_new

def empty_histogram(binning, name="some_histogram", title="some_title"):
    th1_new = root.TH1F(name, title, len(binning)-1, binning)
    return th1_new

class CombineTask:

    def __init__(self, year, algo, pt_bin, jet_version, wp, is_child_of_complex_task=False):
        self.Year = year
        # self.PtBin = pt_bin
        self.Algo = algo
        self.JetVersion = jet_version
        self.WP = wp
        self.Name = '_'.join([algo, pt_bin, jet_version, wp])
        self.IsChild = is_child_of_complex_task
        print 'New CombineTask:', self.Name
        self.WorkDir = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel', year.get('short_name'), 'combine', algo, self.Name)
        print 'Work directory:', self.WorkDir
        if not os.path.isdir(self.WorkDir): sys.exit('Work directory does not exist. Abort.')
        self.LogDir = os.path.join(self.WorkDir, 'log')
        self.InputRootFilePath = os.path.join(self.WorkDir, '_'.join([self.Name, '_ProbeJetHists', algo])+'.root')
        if not os.path.isfile(self.InputRootFilePath): sys.exit('ROOT file with input histograms does not exist. Abort.')
        self.DataCardRootFilePath = self.InputRootFilePath.replace('.root', '_forCombine.root')
        self.DataCardPath = os.path.join(self.WorkDir, '_'.join([self.Name, 'datacard.txt']))
        self.WorkSpacePath = os.path.join(self.WorkDir, '_'.join([self.Name, 'workspace.root']))
        self.PrePostFitShapesPaths = OrderedDict()
        self.PrePostFitShapesForPlotsPaths = OrderedDict()
        for x in ['exp', 'obs']:
            self.PrePostFitShapesPaths[x] = os.path.join(self.WorkDir, '_'.join([self.Name, 'prepostfitshapes_'+x+'.root']))
            self.PrePostFitShapesForPlotsPaths[x] = OrderedDict()
            for y in ['prefit', 'postfit']:
                self.PrePostFitShapesForPlotsPaths[x][y] = self.PrePostFitShapesPaths.get(x).replace('.root', '_forPlots_'+y+'.root')
        self.ControlPlotsPath = os.path.join(self.WorkDir, 'plots')
        self.Variable = 'mSD' if algo=='AK8' else 'mass' if algo=='HOTVR' else None
        self.Channels = OrderedDict() # dict of dicts
        for region in regions:
            self.Channels['_'.join([self.Name, region, self.Variable])] = dict()
        self.Processes = OrderedDict()
        for process in processes.keys():
            process_dict = processes.get(process)
            process_dict['prior_name'] = process
            if processes.get(process).get('id') <= 0: # if process is a signal
                if self.IsChild:
                    new_process_name = '_'.join([process, self.Name]) # self.Name ensures that the TTbar signals are unique per algo/pt_bin/jet_version/wp combination
                else:
                    new_process_name = process
                is_signal = True
            else:
                new_process_name = process
                is_signal = False
            process_dict['is_signal'] = is_signal
            process_dict['rateParam'] = '_'.join(['rate', new_process_name])
            self.Processes[new_process_name] = process_dict
        self.SignalRateParams = OrderedDict()
        for process in self.Processes.keys():
            if self.Processes.get(process).get('is_signal'):
                self.SignalRateParams[process] = self.Processes.get(process).get('rateParam')
        self.POIs = OrderedDict([
            ('SF_'+cat+('_'+self.Name if self.IsChild else ''), self.POI_init()) for cat in global_categories
        ])
        # self.POIs = OrderedDict([
        #     # ('r_overall'+('_'+self.Name if self.IsChild else ''), self.POI_init()),
        #     # ('SF_TTbar'+('_'+self.Name if self.IsChild else ''), self.POI_init()),
        #     ('SF_FullyMerged'+('_'+self.Name if self.IsChild else ''), self.POI_init()),
        #     # ('SF_SemiMerged'+('_'+self.Name if self.IsChild else ''), self.POI_init()),
        #     # ('SF_WMerged'+('_'+self.Name if self.IsChild else ''), self.POI_init()),
        #     # ('SF_QBMerged'+('_'+self.Name if self.IsChild else ''), self.POI_init()),
        #     # ('SF_BkgOrNotMerged'+('_'+self.Name if self.IsChild else ''), self.POI_init()),
        #     # ('SF_NotFullyOrWMerged'+('_'+self.Name if self.IsChild else ''), self.POI_init()),
        #     ('SF_YllufMerged'+('_'+self.Name if self.IsChild else ''), self.POI_init()),
        #     # ('r_pass_FullyMerged', self.POI_init(False)),
        #     # ('r_pass_SemiMerged', self.POI_init(False)),
        #     # ('r_pass_BkgOrNotMerged', self.POI_init(False)),
        #     # ('AntiSF_FullyMerged', self.POI_init(False)),
        #     # ('AntiSF_SemiMerged', self.POI_init(False)),
        #     # ('AntiSF_BkgOrNotMerged', self.POI_init(False)),
        #     # ('r_fail_FullyMerged', self.POI_init(False)),
        #     # ('r_fail_SemiMerged', self.POI_init(False)),
        #     # ('r_fail_BkgOrNotMerged', self.POI_init(False)),
        # ])
        self.ImpactPlotPaths = list()

    def POI_init(self, ind=True): # boolean: independent POI or not
        return OrderedDict([('independent', ind)])

    def write_rootfile(self):
        infile = root.TFile.Open(self.InputRootFilePath, 'READ')
        outfile = root.TFile.Open(self.DataCardRootFilePath, 'RECREATE')
        for channel in self.Channels.keys():
            target_folder = outfile.mkdir(channel)
            target_folder.cd()
            #________________________________
            # Data
            data_obs = infile.Get(channel+'/data_obs')
            data_obs.Write('data_obs')
            self.Channels.get(channel)['binning'] = get_binning(data_obs)
            self.Channels.get(channel)['yields'] = OrderedDict()
            self.Channels.get(channel)['skip'] = OrderedDict()
            self.Channels.get(channel).get('yields')['data_obs'] = data_obs.Integral()
            #________________________________
            # Individual processes with individual uncertainties
            for process in self.Processes.keys():
                prior_process_name = self.Processes.get(process).get('prior_name')
                hist = infile.Get(channel+'/'+prior_process_name)
                self.Channels.get(channel).get('yields')[process] = OrderedDict()
                self.Channels.get(channel).get('yields').get(process)['nominal'] = hist.Integral()
                skip = False
                skip = (True if skip == True else hist.Integral() <= 0)
                hist.Write(process)
                for systematic in systematics:
                    hist_up = infile.Get(channel+'/'+prior_process_name+'_'+systematic+'Up')
                    self.Channels.get(channel).get('yields').get(process)[systematic+'Up'] = hist_up.Integral()
                    skip = (True if skip == True else hist_up.Integral() <= 0)
                    hist_up.Write(process+'_'+systematic+'Up')
                    hist_down = infile.Get(channel+'/'+prior_process_name+'_'+systematic+'Down')
                    self.Channels.get(channel).get('yields').get(process)[systematic+'Down'] = hist_down.Integral()
                    skip = (True if skip == True else hist_down.Integral() <= 0)
                    hist_down.Write(process+'_'+systematic+'Down')
                self.Channels.get(channel).get('skip')[process] = skip
        infile.Close()
        outfile.Close()
        print 'Wrote', self.DataCardRootFilePath
        return self.DataCardRootFilePath

    def write_datacard(self):
        with open(self.DataCardPath, 'w') as file:
            file.write('imax *\n')
            file.write('jmax *\n')
            file.write('kmax *\n')
            file.write('-'*42+'\n')
            #________________________________
            file.write('shapes * * '+os.path.basename(self.DataCardRootFilePath)+' $CHANNEL/$PROCESS $CHANNEL/$PROCESS_$SYSTEMATIC\n')
            file.write('-'*42+'\n')
            #________________________________
            file.write('bin '+' '.join(self.Channels.keys())+'\n')
            file.write('observation'+' -1'*len(self.Channels)+'\n')
            file.write('-'*42+'\n')
            #________________________________
            newline1 = 'bin'
            newline2 = 'process'
            newline3 = 'process'
            newline4 = 'rate'
            for channel in self.Channels.keys():
                for process in self.Processes.keys():
                    if self.Channels.get(channel).get('skip').get(process):
                        continue
                    newline1 += ' '+channel
                    newline2 += ' '+process
                    newline3 += ' '+str(self.Processes.get(process).get('id'))
                    newline4 += ' -1'
            file.write(newline1+'\n')
            file.write(newline2+'\n')
            file.write(newline3+'\n')
            file.write(newline4+'\n')
            file.write('-'*42+'\n')
            #________________________________
            # Luminosity
            file.write('_'.join(['lumi', self.Year.get('short_name')])+' lnN'+(' '+str(1+self.Year.get('lumi_unc')))*len(self.Channels)*len(processes.keys())+'\n')
            # Cross sections (duplicate/colliding with rateParams? Matthias: yes)
            # Only use xsec uncertainty for backgrounds. Do not introduce xsec uncertainty for signal processes but let them float via rateParams (self.SignalRateParams)
            for base_process in base_processes.keys():
                newline1 = '_'.join(['xsec', base_process])+' lnN'
                entries = 0 # number of columns without '-'
                for channel in self.Channels.keys():
                    for process, process_info in self.Processes.items():
                        if self.Channels.get(channel).get('skip').get(process):
                            continue
                        # if process_info.get('base_process') == base_process and not process_info.get('is_signal'):
                        if process_info.get('base_process') == base_process:
                            newline1 += ' '+str(1+base_processes.get(base_process))
                            entries += 1
                        else:
                            newline1 += ' -'
                if entries > 0:
                    file.write(newline1+'\n')
            # Other systematics
            for systematic in systematics:
                newline1 = systematic+' shape'
                for channel in self.Channels.keys():
                    for process, process_info in self.Processes.items():
                        if self.Channels.get(channel).get('skip').get(process):
                            continue
                        if 'tageff' in systematic:
                            if process_info.get('is_signal'):
                                newline1 += ' 1.0'
                            else:
                                newline1 += ' -'
                        else:
                            newline1 += ' 1.0'
                file.write(newline1+'\n')
            # Rate parameters only for signal (backgrounds will be constrained to their nominal yields):
            # for signal_process, signal_rate_param in self.SignalRateParams.items():
            #     file.write(signal_rate_param+' rateParam * '+signal_process+' 1.0 [0.01,5.0] \n')
            # MC statistics
            # https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part2/bin-wise-stats.html
            file.write('* autoMCStats 0 1 1\n') # event-threshold=0, include-signal=1, hist-mode=1
        print 'Wrote', self.DataCardPath
        return self.DataCardPath

    def combine_task(self, task_name, command, run):
        if run:
            time_start = timer()
            os.system('mkdir -p '+self.LogDir)
            logfile_out = open(os.path.join(self.LogDir, task_name+'.log.out'), 'w')
            logfile_err = open(os.path.join(self.LogDir, task_name+'.log.err'), 'w')
            p = subprocess.Popen((command), shell=True, cwd=self.WorkDir, stdout=logfile_out, stderr=logfile_err)
            p.wait()
            logfile_out.close()
            logfile_err.close()
            print '  Runtime: {:.1f} sec'.format(timer() - time_start) # in seconds
            return p.returncode # int
        else:
            return command # string

    def create_workspace(self, run=True):
        print 'Creating workspace'
        command = 'text2workspace.py \\'
        command += self.DataCardPath+' \\'
        command += '-o '+self.WorkSpacePath+' \\'
        command += '-P HiggsAnalysis.CombinedLimit.my_combine_physics_model:lttModel \\'
        # if self.IsChild:
        #     command += '--PO postfix='+self.Name+' \\'
        # command += '--PO verbose \\'
        # for signal_process, signal_rate_param in self.SignalRateParams.items():
        #     command += '''--PO 'map=.*/'''+signal_process+':'+signal_rate_param+'''[1,0,5]\' \\'''
        return self.combine_task('create_workspace', command, run)

    def generate_toys(self, run=True):
        print 'Creating toys'
        command = 'combine -M GenerateOnly \\'
        command += '-t -1 \\'
        # command += '--setParameters '+','.join([x+'=1.' for x in self.SignalRateParams.values()])+' \\'
        # command += '--freezeParameters '+','.join([x for x in self.SignalRateParams.values()])+' \\'
        command += '--setParameters '+','.join([x+'=1.' for x, ind in self.POIs.items() if ind.get('independent')])+' \\'
        command += '--freezeParameters '+','.join([x for x, ind in self.POIs.items() if ind.get('independent')])+' \\'
        command += '--saveToys \\'
        command += '-n _toy \\'
        command += self.WorkSpacePath+' \\'
        return self.combine_task('generate_toys', command, run)

    def multidimfit(self, observed=False, freezeSyst=False, run=True):
        print 'Performing maximum-likelihood fit (MultiDimFit):', expobs(observed, False), ('with frozen systematics' if freezeSyst else '')
        command = 'combine -M MultiDimFit \\'
        command += '-v 2 \\' # more verbosity
        command += '--cminSingleNuisFit \\'
        command += '--cminFallbackAlgo Minuit2,Simplex,0:0.1 \\' # fallback algorithm if default fails; see https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/runningthetool/#generic-minimizer-options
        if freezeSyst:
            datacard_path = os.path.join(self.WorkDir, 'higgsCombine_'+expobs(observed)+'.MultiDimFit.mH120.root')
        else:
            datacard_path = self.WorkSpacePath
        if not os.path.isfile(datacard_path):
            print 'Warning in CombineTask.multidimfit():', datacard_path, 'does not exist'
        command += '--datacard '+datacard_path+' \\'
        command += '-n _'+expobs(observed)+('_freezeSyst' if freezeSyst else '')+' \\'
        if not observed:
            command += '-t -1 \\'
            command += '--toysFile '+os.path.join(self.WorkDir, 'higgsCombine_toy.GenerateOnly.mH120.123456.root')+' \\'
        command += '--algo singles \\'
        if freezeSyst:
            command += '--freezeParameters allConstrainedNuisances --snapshotName MultiDimFit \\'
        command += '--saveFitResult \\'
        command += '--saveWorkspace \\'
        # command += '--fastScan \\'
        # command += '--robustFit=1 \\'
        # command += '--robustHesse=1 \\'
        task_name = '_'.join(['multidimfit', expobs(observed)])
        if freezeSyst:
            task_name += '_freezeSyst'
        return self.combine_task(task_name, command, run)

    def prepostfitshapes(self, observed=False, run=True):
        print 'Calculating prefit and postfit shapes:', expobs(observed, False)
        command = 'PostFitShapesFromWorkspace \\'
        command += '-w '+self.WorkSpacePath+' \\'
        command += '-o '+(self.PrePostFitShapesPaths.get(expobs(observed)))+' \\'
        command += '-f multidimfit_'+expobs(observed)+'.root:fit_mdf \\'
        command += '--postfit \\'
        command += '--sampling \\'
        command += '--print \\'
        command += '--total-shapes \\'
        command += '--covariance \\'
        task_name = '_'.join(['prepostfitshapes', expobs(observed)])
        return self.combine_task(task_name, command, run)

    def prepostfitshapes_for_plots(self, observed=False):
        print 'Writing correctly binned prefit and postfit shapes:', expobs(observed, False)
        infile = root.TFile.Open(self.PrePostFitShapesPaths.get(expobs(observed)), 'READ')
        outfile_paths = list()
        for p in ['prefit', 'postfit']:
            outfile_path = self.PrePostFitShapesForPlotsPaths.get(expobs(observed)).get(p)
            outfile_paths.append(outfile_path)
            outfile = root.TFile.Open(outfile_path, 'RECREATE')
            for channel, channel_info in self.Channels.items():
                target_folder = outfile.mkdir(channel[4:] if channel[:2] == 'ch' else channel)
                target_folder.cd()
                infile_folder = '_'.join([channel, p])
                #________________________________
                # Data
                data_obs_raw = infile.Get(infile_folder+'/data_obs')
                data_obs = return_rebinned(data_obs_raw, channel_info.get('binning'))
                data_obs.Write('data_obs')
                #________________________________
                # Total prediction stack with full uncertainty
                full_stack_raw = infile.Get(infile_folder+'/TotalProcs')
                full_stack = return_rebinned(full_stack_raw, channel_info.get('binning'))
                full_stack.Write('TotalProcs')
                #________________________________
                # Individual processes with full uncertainties
                for process, skip in self.Channels.get(channel).get('skip').items():
                    if skip:
                        # fill in empty histogram
                        process_hist = empty_histogram(channel_info.get('binning'))
                    else:
                        process_hist_raw = infile.Get(infile_folder+'/'+process)
                        process_hist = return_rebinned(process_hist_raw, channel_info.get('binning'))
                    # hist_name = '_'.join(process.split('_')[:2]) if 'TTbar' in process else process
                    hist_name = '_'.join(process.split('_')[:2]) if self.Processes.get(process).get('is_signal') else process
                    process_hist.Write(hist_name)
            outfile.Close()
            print 'Wrote', outfile_path
        infile.Close()
        return outfile_paths

    def plot_prepostfitshapes(self, observed=False, run=True):
        print 'Plotting pre- and post-fit distributions:', expobs(observed, False)
        file_path = os.path.join(self.WorkDir, 'plotting_commands_prepostfitshapes_'+expobs(observed)+'.txt')
        with open(file_path, 'w') as file:
            for p in ['prefit', 'postfit']:
                for channel in self.Channels.keys():
                    channel_name = channel[4:] if channel[:2] == 'ch' else channel
                    newline = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/Analysis/Combine/bin', 'plots')+' '
                    args = list()
                    args.append(p+'Comb')
                    args.append(self.PrePostFitShapesForPlotsPaths.get(expobs(observed)).get(p)) # rootFilePath
                    args.append(channel_name) # folderName / channelName
                    args.append(os.path.join(self.WorkDir, '..', 'plots'))
                    args.append('__'.join([self.Algo, channel_name, p+'Comb.pdf']))
                    args.append(self.Algo+' PUPPI')
                    # args.append('''UL17, 41.5 fb^{#minus1} (13 TeV)''')
                    args.append(self.Year.get('short_name')+''', '''+'''{:.1f}'''.format(self.Year.get('lumi_fb'))+''' fb^{#minus1} (13 TeV)''')
                    args.append('Work in Progress')
                    args.append(str(self.Year.get('lumi_unc'))) # not used for prefitComb and postfitComb plots, only for regular prefit plots; needs to be here so that plotting script has correct number of command line arguments
                    for arg in args:
                        arg = re.escape(arg)
                        newline += ' '+arg
                    newline += '\n'
                    file.write(newline)
        print 'Wrote', file_path
        task_name = '_'.join(['plot_prepostfitshapes', expobs(observed)])
        command = ' '.join(['python', os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/Analysis/Combine/parallel_threading.py'), file_path])
        return self.combine_task(task_name, command, run)

    def impacts(self, observed=False, run=True, parallel=True):
        print 'Calculating nuisance parameter impacts:', expobs(observed, False)
        results = OrderedDict()
        #________________________________
        print 'Performing an initial fit for the signal strength and its uncertainty'
        command1 = 'combineTool.py -M Impacts \\'
        command1 += '-m 0 \\' # required argument
        command1 += '-d '+self.WorkSpacePath+' \\'
        command1 += '--doInitialFit \\'
        command1 += '--robustFit 1 \\'
        if not observed:
            command1 += '-t -1 \\'
            command1 += '--toysFile '+os.path.join(self.WorkDir, 'higgsCombine_toy.GenerateOnly.mH120.123456.root')+' \\'
        task_name1 = '_'.join(['impacts', expobs(observed), 'command1'])
        results[task_name1] = self.combine_task(task_name1, command1, run)
        #________________________________
        print 'Calculating impacts for all nuisance parameters'
        command2 = command1.replace('--doInitialFit', '--doFits')
        if parallel:
            command2 += '--parallel 12 \\' # spawn 6 parallel jobs
        task_name2 = '_'.join(['impacts', expobs(observed), 'command2'])
        results[task_name2] = self.combine_task(task_name2, command2, run)
        #________________________________
        print 'Collecting output, converting to .json format'
        impactsJsonPath = os.path.join(self.WorkDir, 'impacts_'+expobs(observed)+'.json')
        command3 = 'combineTool.py -M Impacts \\'
        command3 += '-m 0 \\' # required argument
        command3 += '-d '+self.WorkSpacePath+' \\'
        command3 += '-o '+impactsJsonPath+' \\'
        task_name3 = '_'.join(['impacts', expobs(observed), 'command3'])
        results[task_name3] = self.combine_task(task_name3, command3, run)
        #________________________________
        os.system('mkdir -p '+self.ControlPlotsPath)
        plot_command_base = 'plotImpacts.py \\'
        plot_command_base += '-i '+impactsJsonPath+' \\'
        plot_command_base += '--POI {0} \\'
        plot_command_base += '-o {1} \\'
        # plot_command_base += '-t '+globalRenameJsonFile+' \\' # rename parameters
        for signal_rate_param in self.POIs.keys():
            print 'Plotting impacts for parameter:', signal_rate_param
            plot_path = os.path.join(self.WorkDir, impactsJsonPath.replace('.json', '_'+signal_rate_param))
            plot_command = plot_command_base.format(signal_rate_param, os.path.join(os.path.basename(self.ControlPlotsPath), os.path.basename(plot_path))) # need to use basename, else combine will try to save pdf to ./{plot_path} --> error
            self.ImpactPlotPaths.append(plot_path+'.pdf')
            plot_task_name = '_'.join(['plot_impacts', expobs(observed), signal_rate_param])
            results[plot_task_name] = self.combine_task(plot_task_name, plot_command, run)
        return results

    def write_scale_factor_file(self, observed=False):
        print 'Writing scale factors to disk'
        # Check out this link to learn how to extract post-fit nuisance parameter values from RooFitResult:
        # https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/102x/src/MultiDimFit.cc#L434
        mdf_file = root.TFile.Open(os.path.join(self.WorkDir, 'multidimfit_'+expobs(observed)+'.root'), 'READ')
        mdf_file_stat = root.TFile.Open(os.path.join(self.WorkDir, 'multidimfit_'+expobs(observed)+'_freezeSyst.root'), 'READ')
        fit_result = mdf_file.Get('fit_mdf') # RooFitResult
        fit_result_stat = mdf_file_stat.Get('fit_mdf') # RooFitResult
        for poi_name in self.POIs.keys():
            rf = fit_result.floatParsFinal().find(poi_name)
            rf_stat = fit_result_stat.floatParsFinal().find(poi_name)
            bestFitVal = rf.getVal()
            if abs(bestFitVal - rf_stat.getVal()) > 0.01: # check if first fit and syst-frozen fit lead to approximately same central value
                print 'bestFitVal regular:', bestFitVal
                print 'bestFitVal frozenSyst:', rf_stat.getVal()
                print 'difference:', bestFitVal - rf_stat.getVal()
                print 'Warning: best-fit scale factor central values of regular fit and fit with frozen systematics differ from each other'
            #________________________________
            hiErr = +(rf.getMax('err68') - bestFitVal if rf.hasRange('err68') else rf.getAsymErrorHi())
            loErr = -(rf.getMin('err68') - bestFitVal if rf.hasRange('err68') else rf.getAsymErrorLo())
            maxError = max(max(hiErr, loErr), rf.getError())
            if abs(hiErr) < 0.001*maxError:
                # print " Warning - No valid high-error found, will report difference to maximum of range for : ", rf.GetName()
                # hiErr = -bestFitVal + rf.getMax()
                print 'Warning: no valid asymmetric high-error found, will use getError() result'
                hiErr = rf.getError()
            if abs(loErr) < 0.001*maxError:
                # print " Warning - No valid low-error found, will report difference to minimum of range for : ", rf.GetName()
                # loErr = +bestFitVal - rf.getMin()
                print 'Warning: no valid asymmetric low-error found, will use getError() result'
                loErr = rf.getError()
            #________________________________
            hiErr_stat = +(rf_stat.getMax('err68') - bestFitVal if rf_stat.hasRange('err68') else rf_stat.getAsymErrorHi())
            loErr_stat = -(rf_stat.getMin('err68') - bestFitVal if rf_stat.hasRange('err68') else rf_stat.getAsymErrorLo())
            maxError_stat = max(max(hiErr_stat, loErr_stat), rf_stat.getError())
            if abs(hiErr_stat) < 0.001*maxError_stat:
                # print " Warning - No valid high-error found, will report difference to maximum of range for : ", rf.GetName()
                # hiErr_stat = -bestFitVal + rf_stat.getMax()
                print 'Warning: no valid asymmetric high-error found, will use getError() result'
                hiErr_stat = rf.getError()
            if abs(loErr_stat) < 0.001*maxError_stat:
                # print " Warning - No valid low-error found, will report difference to minimum of range for : ", rf.GetName()
                # loErr_stat = +bestFitVal - rf_stat.getMin()
                print 'Warning: no valid asymmetric low-error found, will use getError() result'
                loErr_stat = rf.getError()
            #________________________________
            self.POIs.get(poi_name)[expobs(observed)] = OrderedDict()
            self.POIs.get(poi_name).get(expobs(observed))['bestFitVal'] = bestFitVal
            self.POIs.get(poi_name).get(expobs(observed))['loErr'] = loErr
            self.POIs.get(poi_name).get(expobs(observed))['hiErr'] = hiErr
            self.POIs.get(poi_name).get(expobs(observed))['loErr_stat'] = loErr_stat
            self.POIs.get(poi_name).get(expobs(observed))['hiErr_stat'] = hiErr_stat
            loErr_syst = loErr*loErr - loErr_stat*loErr_stat
            if loErr_syst >= 0:
                loErr_syst = root.TMath.Sqrt(loErr_syst)
            else:
                print 'Warning: lower total error smaller than lower statistical component; something went wrong! Lower systematic component set to zero'
                loErr_syst = 0
            hiErr_syst = hiErr*hiErr - hiErr_stat*hiErr_stat
            if hiErr_syst >= 0:
                hiErr_syst = root.TMath.Sqrt(hiErr_syst)
            else:
                print 'Warning: upper total error smaller than upper statistical component; something went wrong! Upper systematic component set to zero'
                hiErr_syst = 0
            self.POIs.get(poi_name).get(expobs(observed))['loErr_syst'] = loErr_syst
            self.POIs.get(poi_name).get(expobs(observed))['hiErr_syst'] = hiErr_syst
            print "  {:55s}:  {:.3f}  -{:.3f} / +{:.3f} (tot.)  [ -{:.3f} / +{:.3f} (stat.)  |  -{:.3f} / +{:.3f} (syst.) ]".format(poi_name, bestFitVal, loErr, hiErr, loErr_stat, hiErr_stat, loErr_syst, hiErr_syst)
        sf_file_path = os.path.join(self.WorkDir, 'scale_factors_'+expobs(observed)+'.root')
        if observed:
            self.sf_file_path_obs = sf_file_path
        else:
            self.sf_file_path_exp = sf_file_path
        sf_file = root.TFile.Open(sf_file_path, 'RECREATE')
        sf_file.cd()
        for cat in global_categories:
            x = []
            y = []
            exl = []
            exh = []
            eyl = []
            eyh = []
            eyl_stat = []
            eyh_stat = []
            eyl_syst = []
            eyh_syst = []
            for poi_name, poi_info in self.POIs.items():
                if cat not in poi_name:
                    continue
                pt_interval_identified = False
                pt_max = 0.
                pt_min = 0.
                for pt_interval_key, pt_interval_info in pt_intervals.items():
                    if pt_interval_key in poi_name:
                        pt_interval_identified = True
                        pt_max = pt_interval_info.get('pt_max')
                        pt_min = pt_interval_info.get('pt_min')
                        break
                if not pt_interval_identified:
                    for pt_interval_key, pt_interval_info in pt_intervals.items():
                        if pt_interval_key in self.Name:
                            pt_interval_identified = True
                            pt_max = pt_interval_info.get('pt_max')
                            pt_min = pt_interval_info.get('pt_min')
                            break
                if not pt_interval_identified:
                    print 'Warning: could not convert POI name into correct pt range'
                half_x_width = (pt_max - pt_min)/2.
                x.append(half_x_width + pt_min)
                exl.append(half_x_width)
                exh.append(half_x_width)
                y.append(poi_info.get(expobs(observed))['bestFitVal'])
                eyl.append(poi_info.get(expobs(observed))['loErr'])
                eyh.append(poi_info.get(expobs(observed))['hiErr'])
                eyl_stat.append(poi_info.get(expobs(observed))['loErr_stat'])
                eyh_stat.append(poi_info.get(expobs(observed))['hiErr_stat'])
                eyl_syst.append(poi_info.get(expobs(observed))['loErr_syst'])
                eyh_syst.append(poi_info.get(expobs(observed))['hiErr_syst'])
            sf_graph = root.TGraphAsymmErrors(len(x), np.array(x), np.array(y), np.array(exl), np.array(exh), np.array(eyl), np.array(eyh))
            sf_graph_stat = root.TGraphAsymmErrors(len(x), np.array(x), np.array(y), np.array(exl), np.array(exh), np.array(eyl_stat), np.array(eyh_stat))
            sf_graph_syst = root.TGraphAsymmErrors(len(x), np.array(x), np.array(y), np.array(exl), np.array(exh), np.array(eyl_syst), np.array(eyh_syst))
            sf_graph.Write(cat+'_tot')
            sf_graph_stat.Write(cat+'_stat')
            sf_graph_syst.Write(cat+'_syst')
        print 'Wrote', sf_file_path

    def run_all(self):
        # self.write_rootfile()
        # self.write_datacard()
        # self.create_workspace()

        # Expected:
        # self.generate_toys()
        # self.multidimfit()
        # self.multidimfit(freezeSyst=True)
        # self.write_scale_factor_file()
        # self.prepostfitshapes()
        # self.prepostfitshapes_for_plots()
        if args.run_impacts:
            self.impacts()

        # Observed:
        # self.multidimfit(True)
        # self.multidimfit(True, freezeSyst=True)
        self.write_scale_factor_file(True)
        # self.prepostfitshapes(True)
        # self.prepostfitshapes_for_plots(True)
        # self.plot_prepostfitshapes(True)
        if args.run_impacts:
            self.impacts(True)


class ComplexCombineTask(CombineTask):

    def __init__(self, year, algo, jet_version, wp):
        self.Year = year
        self.Algo = algo
        self.JetVersion = jet_version
        self.WP = wp
        self.Name = '_'.join([algo, 'PtAll', jet_version, wp])
        print 'New ComplexCombineTask:', self.Name
        self.WorkDir = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel', year.get('short_name'), 'combine', algo, self.Name)
        os.system('mkdir -p '+self.WorkDir)
        print 'Work directory:', self.WorkDir
        self.LogDir = os.path.join(self.WorkDir, 'log')
        self.DataCardPath = os.path.join(self.WorkDir, '_'.join([self.Name, 'datacard.txt']))
        self.WorkSpacePath = os.path.join(self.WorkDir, '_'.join([self.Name, 'workspace.root']))
        self.PrePostFitShapesPaths = OrderedDict()
        self.PrePostFitShapesForPlotsPaths = OrderedDict()
        for x in ['exp', 'obs']:
            self.PrePostFitShapesPaths[x] = os.path.join(self.WorkDir, '_'.join([self.Name, 'prepostfitshapes_'+x+'.root']))
            self.PrePostFitShapesForPlotsPaths[x] = OrderedDict()
            for y in ['prefit', 'postfit']:
                self.PrePostFitShapesForPlotsPaths[x][y] = self.PrePostFitShapesPaths.get(x).replace('.root', '_forPlots_'+y+'.root')
        self.ControlPlotsPath = os.path.join(self.WorkDir, 'plots')
        self.Variable = 'mSD' if algo=='AK8' else 'mass' if algo=='HOTVR' else None
        # self.Channels = list()
        self.Channels = OrderedDict()
        self.Processes = OrderedDict()
        self.SignalRateParams = OrderedDict()
        self.POIs = OrderedDict()
        self.ChildCombineTasks = list()
        for pt_bin in pt_bins.get(algo): # for now, a ComplexCombineTask is a simultaneous CombineTask in all pt_bins
            self.ChildCombineTasks.append(CombineTask(year, algo, pt_bin, jet_version, wp, True))
        number_of_ccts = 0
        for cct in self.ChildCombineTasks:
            number_of_ccts += 1
            cct.write_rootfile()
            cct.write_datacard()
            for channel, channel_info in cct.Channels.items():
                channel_name = 'ch'+str(number_of_ccts)+'_'+channel # after datacard combination, channels have prefixes like 'ch1_', 'ch1_', 'ch2_', 'ch2_',  etc.
                self.Channels[channel_name] = channel_info
            for process in cct.Processes.keys():
                 # Processes like WJets will be called multiple times here, but in the end this dict will contain
                 # it only once. On the other hand, the TTbar signal processes will be truly unique per pt_bin
                self.Processes[process] = cct.Processes.get(process)
            for signal_process, signal_rate_param in cct.SignalRateParams.items():
                self.SignalRateParams[signal_process] = signal_rate_param
            for poi, ind in cct.POIs.items():
                self.POIs[poi] = self.POI_init(ind.get('independent'))
        self.ImpactPlotPaths = list()

    def write_rootfile(self):
        '''Disabled method of superclass. ComplexCombineTasks inherit DataCardRootFilePaths from their ChildCombineTasks.'''
        pass

    def write_datacard(self):
        print 'Combining datacards'
        command = 'combineCards.py '
        command += ' '.join([cct.DataCardPath for cct in self.ChildCombineTasks])
        command += ' &> '+self.DataCardPath
        return self.combine_task('combineCards', command, True)

    def create_workspace(self, run=True):
        print 'Creating workspace'
        command = 'text2workspace.py \\'
        command += self.DataCardPath+' \\'
        command += '-o '+self.WorkSpacePath+' \\'
        command += '-P HiggsAnalysis.CombinedLimit.my_combine_physics_model:lttModel \\'
        command += '--PO tasks='+','.join([cct.Name for cct in self.ChildCombineTasks])+' \\'
        return self.combine_task('create_workspace', command, run)




# c = CombineTask('UL17', 'HOTVR', 'Pt250to300', 'HOTVRCuts', 'Standard')
# c = CombineTask('UL17', 'HOTVR', 'Pt200to250', 'HOTVRCuts', 'Standard')
# c.write_rootfile()
# c.write_datacard()
# c.create_workspace()
# c.run_all()

# c.generate_toys()
# c.multidimfit()
# c.prepostfitshapes()
# c.impacts()



all_tasks = list()
all_tasks = list(itertools.product(['AK8'], jet_versions.get('AK8'), wps.get('AK8')))
all_tasks.extend(list((itertools.product(['HOTVR'], jet_versions.get('HOTVR'), wps.get('HOTVR')))))

all_tasks = [OrderedDict([
    ('probejet_coll', x[0]),
    ('jet_version', x[1]),
    ('wp', x[2]),
    # ('variables', set_variables(x)), ### Schreibe wirklich nur die Variablen aus, die Du brauchst!!!!! (entferne den obigen for-loop in create_input_histograms)
]) for x in all_tasks]

# print all_tasks
#

all_tasks = [
OrderedDict([('probejet_coll', 'AK8'), ('jet_version', 'All'), ('wp', 'BkgEff0p001')]),
OrderedDict([('probejet_coll', 'AK8'), ('jet_version', 'All'), ('wp', 'BkgEff0p005')]),
OrderedDict([('probejet_coll', 'AK8'), ('jet_version', 'All'), ('wp', 'BkgEff0p010')]),
OrderedDict([('probejet_coll', 'AK8'), ('jet_version', 'All'), ('wp', 'BkgEff0p025')]),
OrderedDict([('probejet_coll', 'AK8'), ('jet_version', 'All'), ('wp', 'BkgEff0p050')]),
OrderedDict([('probejet_coll', 'AK8'), ('jet_version', 'BTag'), ('wp', 'BkgEff0p001')]),
OrderedDict([('probejet_coll', 'AK8'), ('jet_version', 'BTag'), ('wp', 'BkgEff0p005')]),
OrderedDict([('probejet_coll', 'AK8'), ('jet_version', 'BTag'), ('wp', 'BkgEff0p010')]),
OrderedDict([('probejet_coll', 'AK8'), ('jet_version', 'BTag'), ('wp', 'BkgEff0p025')]),
OrderedDict([('probejet_coll', 'AK8'), ('jet_version', 'BTag'), ('wp', 'BkgEff0p050')]),
OrderedDict([('probejet_coll', 'HOTVR'), ('jet_version', 'HOTVRCuts'), ('wp', 'Standard')])
]

tasks = []

for task in all_tasks:
    c = ComplexCombineTask(arg_year, task.get('probejet_coll'), task.get('jet_version'), task.get('wp'))
    c.write_datacard()
    c.create_workspace()
    c.run_all()
    # if args.run_impacts:
    #     c.impacts()
    #     c.impacts(True)
    tasks.append(c)

from combine_scale_factor_files import combine_scale_factor_files
combine_scale_factor_files(tasks)

# c = ComplexCombineTask('UL17', 'HOTVR', 'HOTVRCuts', 'Standard')
# c.write_datacard()
# c.create_workspace()
# c.run_all()



# # d = ComplexCombineTask('UL17', 'HOTVR', 'HOTVRCuts', 'Standard')
# d = ComplexCombineTask('UL17', 'AK8', 'All', 'BkgEff0p050')
# # print d.Processes
# # # print d.SignalRateParams
# # print d.POIs
# # print d.WorkDir
# d.write_datacard()
# d.create_workspace()
# d.run_all()
