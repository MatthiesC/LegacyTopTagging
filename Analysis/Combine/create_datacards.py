#!/usr/bin/env python

import os
import sys
from collections import OrderedDict
import itertools
import ROOT as root
import subprocess
from timeit import default_timer as timer
import numpy as np
import warnings

# https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiLUM
years = {
    'UL17': {
        'lumi': 41.53,
        'lumi_unc': 0.023,
    },
    'UL18': {
        'lumi': 59.74,
        'lumi_unc': 0.025,
    },
}

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
    'Pass',
    'Fail',
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
    # ('ST', {
    #     'id': 200,
    #     'base_process': 'ST',
    # }),
    ('ST_FullyMerged', {
        'id': 201,
        'base_process': 'ST',
    }),
    # ('ST_SemiMerged', {
    #     'id': 202,
    #     'base_process': 'ST',
    # }),
    # ('ST_WMerged', {
    #     'id': 203,
    #     'base_process': 'ST',
    # }),
    # ('ST_QBMerged', {
    #     'id': 204,
    #     'base_process': 'ST',
    # }),
    # ('ST_BkgOrNotMerged', {
    #     'id': 205,
    #     'base_process': 'ST',
    # }),
    # ('ST_NotFullyOrWMerged', {
    #     'id': 206,
    #     'base_process': 'ST',
    # }),
    ('ST_YllufMerged', {
        'id': 207,
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
    # 'topptA',
    # 'topptB',
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
        self.Name = '_'.join([algo, pt_bin, jet_version, wp])
        self.IsChild = is_child_of_complex_task
        print 'New CombineTask:', self.Name
        self.WorkDir = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel', year, 'combine', algo, self.Name)
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
                self.PrePostFitShapesForPlotsPaths[x][y] = self.PrePostFitShapesPaths.get(x).replace('.root', '_forPlots.root')
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
        self.POIs = OrderedDict([ # boolean: independent POI or not
            # ('r_overall'+('_'+self.Name if self.IsChild else ''), True),
            # ('SF_TTbar'+('_'+self.Name if self.IsChild else ''), True),
            ('SF_FullyMerged'+('_'+self.Name if self.IsChild else ''), True),
            # ('SF_SemiMerged'+('_'+self.Name if self.IsChild else ''), True),
            # ('SF_WMerged'+('_'+self.Name if self.IsChild else ''), True),
            # ('SF_QBMerged'+('_'+self.Name if self.IsChild else ''), True),
            # ('SF_BkgOrNotMerged'+('_'+self.Name if self.IsChild else ''), True),
            # ('SF_NotFullyOrWMerged'+('_'+self.Name if self.IsChild else ''), True),
            ('SF_YllufMerged'+('_'+self.Name if self.IsChild else ''), True),
            # ('r_pass_FullyMerged', False),
            # ('r_pass_SemiMerged', False),
            # ('r_pass_BkgOrNotMerged', False),
            # ('AntiSF_FullyMerged', False),
            # ('AntiSF_SemiMerged', False),
            # ('AntiSF_BkgOrNotMerged', False),
            # ('r_fail_FullyMerged', False),
            # ('r_fail_SemiMerged', False),
            # ('r_fail_BkgOrNotMerged', False),
        ])
        self.ImpactPlotPaths = list()

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
            file.write('_'.join(['lumi', self.Year])+' lnN'+(' '+str(1+years.get(self.Year).get('lumi_unc')))*len(self.Channels)*len(processes.keys())+'\n')
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
        command += '--setParameters '+','.join([x+'=1.' for x, ind in self.POIs.items() if ind])+' \\'
        command += '--freezeParameters '+','.join([x for x, ind in self.POIs.items() if ind])+' \\'
        command += '--saveToys \\'
        command += '-n _toy \\'
        command += self.WorkSpacePath+' \\'
        return self.combine_task('generate_toys', command, run)

    def multidimfit(self, observed=False, run=True):
        print 'Performing maximum-likelihood fit (MultiDimFit):', expobs(observed, False)
        command = 'combine -M MultiDimFit \\'
        if observed:
            command += '-n _obs \\'
        else:
            command += '-n _exp \\'
            command += '-t -1 \\'
            command += '--toysFile '+os.path.join(self.WorkDir, 'higgsCombine_toy.GenerateOnly.mH120.123456.root')+' \\'
        command += '--algo singles \\'
        command += '--datacard '+self.WorkSpacePath+' \\'
        command += '--saveFitResult \\'
        # command += '--saveWorkspace \\'
        # command += '--fastScan \\'
        # command += '--robustFit=1 \\'
        # command += '--robustHesse=1 \\'
        task_name = '_'.join(['multidimfit', expobs(observed)])
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
                target_folder = outfile.mkdir(channel)
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
                for process in self.Processes.keys():
                    if self.Channels.get(channel).get('skip').get(process):
                        # fill in empty histogram
                        process_hist = empty_histogram(channel_info.get('binning'))
                    else:
                        process_hist_raw = infile.Get(infile_folder+'/'+process)
                        process_hist = return_rebinned(process_hist_raw, channel_info.get('binning'))
                    process_hist.Write(process)
            outfile.Close()
            print 'Wrote', outfile_path
        infile.Close()
        return outfile_paths

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


    def run_all(self):
        # self.write_rootfile()
        # self.write_datacard()
        # self.create_workspace()
        # Expected:
        self.generate_toys()
        self.multidimfit()
        # self.prepostfitshapes()
        # self.prepostfitshapes_for_plots()
        # self.impacts()
        # Observed:
        self.multidimfit(True)
        # self.prepostfitshapes(True)
        # self.prepostfitshapes_for_plots(True)
        # self.impacts(True)


class ComplexCombineTask(CombineTask):

    def __init__(self, year, algo, jet_version, wp):
        self.Year = year
        self.Name = '_'.join([algo, 'PtAll', jet_version, wp])
        print 'New ComplexCombineTask:', self.Name
        self.WorkDir = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel', year, 'combine', algo, self.Name)
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
                self.PrePostFitShapesForPlotsPaths[x][y] = self.PrePostFitShapesPaths.get(x).replace('.root', '_forPlots.root')
        self.ControlPlotsPath = os.path.join(self.WorkDir, 'plots')
        self.Variable = 'mSD' if algo=='AK8' else 'mass' if algo=='HOTVR' else None
        # self.Channels = list()
        self.Processes = OrderedDict()
        self.SignalRateParams = OrderedDict()
        self.POIs = OrderedDict()
        self.ChildCombineTasks = list()
        for pt_bin in pt_bins.get(algo): # for now, a ComplexCombineTask is a simultaneous CombineTask in all pt_bins
            self.ChildCombineTasks.append(CombineTask(year, algo, pt_bin, jet_version, wp, True))
        for cct in self.ChildCombineTasks:
            cct.write_rootfile()
            cct.write_datacard()
            for process in cct.Processes.keys():
                 # Processes like WJets will be called multiple times here, but in the end this dict will contain
                 # it only once. On the other hand, the TTbar signal processes will be truly unique per pt_bin
                self.Processes[process] = cct.Processes.get(process)
            for signal_process, signal_rate_param in cct.SignalRateParams.items():
                self.SignalRateParams[signal_process] = signal_rate_param
            for poi, ind in cct.POIs.items():
                self.POIs[poi] = ind
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
# for task in all_tasks:
#     c = ComplexCombineTask('UL17', task.get('probejet_coll'), task.get('jet_version'), task.get('wp'))
#     c.write_datacard()
#     c.create_workspace()
#     c.run_all()

c = ComplexCombineTask('UL17', 'HOTVR', 'HOTVRCuts', 'Standard')
c.write_datacard()
c.create_workspace()
c.run_all()



# # d = ComplexCombineTask('UL17', 'HOTVR', 'HOTVRCuts', 'Standard')
# d = ComplexCombineTask('UL17', 'AK8', 'All', 'BkgEff0p050')
# # print d.Processes
# # # print d.SignalRateParams
# # print d.POIs
# # print d.WorkDir
# d.write_datacard()
# d.create_workspace()
# d.run_all()
