#!/usr/bin/env python2

from __future__ import print_function

import os
import sys
import subprocess
import ROOT as root
from copy import deepcopy
from collections import OrderedDict
from timeit import default_timer as timer
import numpy as np
import warnings
import re
from tqdm import tqdm

import sys
sys.path.append(os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/Analysis'))
from constants import Systematics, _BANDS, _TAGGERS, _PT_INTERVALS_TANDP_AK8_T, _PT_INTERVALS_TANDP_AK8_W, _PT_INTERVALS_TANDP_HOTVR, _YEARS, get_variable_binning_xlabel_xunit, Primes

systematics = Systematics(include_jes_splits=False, blacklist=['sfelec_trigger', 'sfmu_iso'])
systs = systematics.base

sys.path.append(os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/NicePlots/python'))
from plotter import NiceStackWithRatio, Process, human_format


years = [
'UL16preVFP',
'UL16postVFP',
'UL17',
'UL18',
]
regions = [
'Pass',
'Fail',
]
channels = [
'muo',
'ele',
]
# prepostfits = [
# 'prefitRaw',
# # 'postfitCombine',
# ]
taggers = [
'ak8_t__tau',
'ak8_t_btagDJet__tau',
'ak8_t_btagDCSV__tau',
'hotvr_t__tau',
'ak8_w__partnet',
'ak8_t__MDdeepak8',
]
taggers = {k: _TAGGERS[k] for k in taggers}


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


class ScaleFactorFits():

    def __init__(self,
        tagger_name,
        wp,
        mode='TagAndProbe',
        mscSplitting = 'mscTop3',
        total_range = False,
        task_name_suffix = '',
    ):

        print('Hello World from ScaleFactorFits! :-)')

        if mode not in ['TagAndProbe', 'Hybrid', 'ThetaLike']:
            sys.exit('Mode not supported:', mode)
        self.mode = mode

        self.primes = Primes()

        self.years = years
        self.regions = regions
        self.channels = channels
        self.task_name_suffix = ''

        self.mscSplitting = mscSplitting
        if self.mscSplitting == 'mscTop2':
            self.merge_scenarios = OrderedDict([
                ('FullyMerged', 'ThetaLike' if self.mode == 'ThetaLike' else 'TagAndProbe'),
                ('YllufMerged', 'TagAndProbe' if self.mode == 'TagAndProbe' else 'ThetaLike'),
            ])
            self.processes_Plotter = [
                Process('VJetsAndVV', 'V+jets, VV', root.kSpring-3),
                Process('ST__MSc_Background', 'Background', root.kOrange+4),
                Process('ST__MSc_YllufMerged', 'Not merged', root.kOrange-3),
                Process('ST__MSc_FullyMerged', 'Fully merged', root.kOrange),
                Process('TTbar__MSc_Background', 'Background', root.kPink-7),
                Process('TTbar__MSc_YllufMerged', 'Not merged', root.kPink+4),
                Process('TTbar__MSc_FullyMerged', 'Fully merged', root.kPink-3),
                Process('QCD', 'Multijet', root.kAzure+10),
            ]
        elif self.mscSplitting == 'mscTop3':
            self.merge_scenarios = OrderedDict([
                ('FullyMerged', 'ThetaLike' if self.mode == 'ThetaLike' else 'TagAndProbe'),
                ('SemiMerged', 'TagAndProbe' if self.mode == 'TagAndProbe' else 'ThetaLike'),
                ('NotMerged', 'TagAndProbe' if self.mode == 'TagAndProbe' else 'ThetaLike'),
            ])
            self.processes_Plotter = [
                Process('VJetsAndVV', 'V+jets, VV', root.kSpring-3),
                Process('ST__MSc_Background', 'Background', root.kOrange+4),
                Process('ST__MSc_NotMerged', 'Not merged', root.kOrange-3),
                Process('ST__MSc_SemiMerged', 'Semi-merged', root.kYellow),
                Process('ST__MSc_FullyMerged', 'Fully merged', root.kOrange),
                Process('TTbar__MSc_Background', 'Background', root.kPink-7),
                Process('TTbar__MSc_NotMerged', 'Not merged', root.kPink+4),
                Process('TTbar__MSc_SemiMerged', 'Semi-merged', root.kPink+6),
                Process('TTbar__MSc_FullyMerged', 'Fully merged', root.kPink-3),
                Process('QCD', 'Multijet', root.kAzure+10),
            ]
        elif self.mscSplitting == 'mscW3':
            self.merge_scenarios = OrderedDict([
                ('FullyMerged', 'TagAndProbe' if self.mode == 'TagAndProbe' else 'ThetaLike'),
                ('WMerged', 'ThetaLike' if self.mode == 'ThetaLike' else 'TagAndProbe'),
                ('NotTopOrWMerged', 'TagAndProbe' if self.mode == 'TagAndProbe' else 'ThetaLike'),
            ])
            self.processes_Plotter = [
                Process('VJetsAndVV', 'V+jets, VV', root.kSpring-3),
                Process('ST__MSc_Background', 'Background', root.kOrange+4),
                Process('ST__MSc_NotTopOrWMerged', 'Not merged', root.kOrange-3),
                Process('ST__MSc_WMerged', 'W merged', root.kYellow),
                Process('ST__MSc_FullyMerged', 't merged', root.kOrange),
                Process('TTbar__MSc_Background', 'Background', root.kPink-7),
                Process('TTbar__MSc_NotTopOrWMerged', 'Not merged', root.kPink+4),
                Process('TTbar__MSc_WMerged', 'W merged', root.kPink+6),
                Process('TTbar__MSc_FullyMerged', 't merged', root.kPink-3),
                Process('QCD', 'Multijet', root.kAzure+10),
            ]
        else:
            sys.exit('Invalid merge scenario splitting:', mscSplitting)
        processes_Plotter_temp = OrderedDict()
        for index, proc in enumerate(self.processes_Plotter):
            proc.index = index
            processes_Plotter_temp[proc.name] = proc
        self.processes_Plotter = processes_Plotter_temp

        self.tagger = taggers[tagger_name]
        if isinstance(wp, str):
            self.wp = self.tagger.get_wp(wp_name=wp)
        elif isinstance(wp, (int, long,)):
            self.wp = self.tagger.get_wp(wp_index=wp)
        else:
            sys.exit('Cannot determine WP')

        sorted_pt_bins = []
        for pt_bin in self.tagger.var_intervals.values():
            if pt_bin.fit == True:
                sorted_pt_bins.append(pt_bin)
        sorted_pt_bins = sorted(sorted_pt_bins, key=lambda x: x.var_min, reverse=False)
        self.pt_bins = OrderedDict()
        for pt_bin in sorted_pt_bins:
            self.pt_bins[pt_bin.name] = pt_bin

        self.pt_bin_total_range = None
        for pt_bin in self.tagger.var_intervals.values():
            if pt_bin.total_range == True:
                self.pt_bin_total_range = pt_bin
                break
        self.total_range = total_range or False
        if self.total_range:
            self.pt_bins = OrderedDict([(self.pt_bin_total_range.name, self.pt_bin_total_range)])
            self.task_name_suffix += 'PtTotal'
        else:
            self.task_name_suffix += 'PtSplit'

        if len(self.years) == 4:
            self.task_name_suffix += '-Run2'
        elif len(self.years) == 1:
            self.task_name_suffix += '-'+self.years[0]
        else:
            self.task_name_suffix += '-'
            for year in self.years:
                self.task_name_suffix += year

        self.task_name_suffix += task_name_suffix # FIXME be more creative with this name
        self.task_name = '-'.join(['combineTask', self.tagger.name, self.wp.name]) + ('-'+self.task_name_suffix if len(self.task_name_suffix) else '')
        self.workdir = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel/combine', self.tagger.name, 'workdirs', self.task_name)
        print('Working directory:', self.workdir)
        self.logdir = os.path.join(self.workdir, 'log')
        command = "mkdir -p {}".format(self.logdir)
        p = subprocess.Popen((command), shell=True)
        p.wait()

        self.rootfile_path = os.path.join(self.workdir, self.task_name+'-DataCardTemplates.root')
        self.datacard_path = os.path.join(self.workdir, self.task_name+'-DataCard.dat')
        self.workspace_path = os.path.join(self.workdir, self.task_name+'-WorkSpace.root')

        self.prepostfitshapes_paths = OrderedDict()
        self.prepostfitshapes_for_plots_paths = OrderedDict()
        for x in ['exp', 'obs']:
            self.prepostfitshapes_paths[x] = os.path.join(self.workdir, self.task_name+'-PrePostFitShapes_'+x+'.root')
            self.prepostfitshapes_for_plots_paths[x] = OrderedDict()
            for y in ['prefit', 'postfit']:
                self.prepostfitshapes_for_plots_paths[x][y] = self.prepostfitshapes_paths.get(x).replace('.root', '_forPlots_'+y+'.root')
        self.control_plots_dir = os.path.join(self.workdir, 'plots_templates')
        command = "mkdir -p {}".format(self.control_plots_dir)
        p = subprocess.Popen((command), shell=True)
        p.wait()

        self.combine_channels = OrderedDict()
        for year in self.years:
            for region in self.regions:
                for channel in self.channels:
                    for pt_bin in self.pt_bins.values():
                        combine_channel_name = '-'.join(['CombBin', 'Main', region, channel, year, pt_bin.name])
                        combine_channel_dict = OrderedDict([
                            ('name', combine_channel_name),
                            ('year', year),
                            ('region', region),
                            ('channel', channel),
                            ('pt_bin', pt_bin),
                            ('id_factor', self.primes.get_next_prime()),
                        ])
                        self.combine_channels.setdefault(combine_channel_name, combine_channel_dict)

        self.sf_naming_scheme = '__'.join(['SF', r'{msc}', r'{year}', r'{pt_bin_name}'])
        self.pois = OrderedDict()
        for msc in self.merge_scenarios.keys():
            for year in self.years:
                for pt_bin in self.pt_bins.values():
                    # poi_name = '__'.join(['SF', msc, year, pt_bin.name])
                    poi_name = self.sf_naming_scheme.format(msc=msc, year=year, pt_bin_name=pt_bin.name)
                    poi_dict = OrderedDict([
                        ('name', poi_name),
                        ('msc', msc),
                        ('year', year),
                        ('pt_bin', pt_bin),
                    ])
                    self.pois.setdefault(poi_name, poi_dict)

        # numbers are the assumed xsec uncertainty
        self.base_processes = OrderedDict([
            ('TTbar', 0.05),
            ('ST', 0.5),
            ('VJetsAndVV', 0.5),
            ('QCD', 1.00),
        ])

        self.processes = OrderedDict()
        for tt_or_st in ['TTbar', 'ST']:
            for msc in self.merge_scenarios.keys():
                process_name = tt_or_st+'__MSc_'+msc
                fit_name_base = tt_or_st+'__MSc_'+msc+r'__{year}__{pt_bin}'
                process_dict = OrderedDict([
                    ('name', process_name),
                    ('is_signal', True),
                    ('fit_name_base', fit_name_base),
                    ('id_factor', self.primes.get_next_prime()),
                    ('base_process', tt_or_st),
                ])
                self.processes.setdefault(process_name, process_dict)
            process_name = tt_or_st+'__MSc_Background'
            process_dict = OrderedDict([
                ('name', process_name),
                ('is_signal', False),
                ('fit_name_base', process_name),
                ('id_factor', self.primes.get_next_prime()),
                ('base_process', tt_or_st),
            ])
            self.processes.setdefault(process_name, process_dict)
        for proc in ['VJetsAndVV', 'QCD']:
            process_name = proc
            process_dict = OrderedDict([
                ('name', process_name),
                ('is_signal', False),
                ('fit_name_base', process_name),
                ('id_factor', self.primes.get_next_prime()),
                ('base_process', proc),
            ])
            self.processes.setdefault(process_name, process_dict)

        self.systs = OrderedDict()
        for syst in systs.values():
            if syst.tandp:
                self.systs[syst.name] = syst


    def get_fit_name(self, process, combine_channel):

        '''helper function'''

        return process['fit_name_base'].format(year=combine_channel['year'], pt_bin=combine_channel['pt_bin'].name)


    def get_process_id(self, process, combine_channel):

        '''helper function'''

        # id = -1 if process['is_signal'] else 1
        # id *= process['id_factor']
        # id *= combine_channel['id_factor']

        id = process['id_factor']
        if process['is_signal']:
            id *= -1
            id *= combine_channel['id_factor'] # signals are treated as independent processes per combine_channel; backgrounds need to have same id across combine_channels

        return id


    def write_rootfile(self):

        '''
        Needs to be run before `self.write_datacard()`
        '''

        print('Writing ROOT file for datacard...')

        outfile = root.TFile.Open(self.rootfile_path, 'RECREATE')
        n_histograms = len(self.combine_channels)*(1+len(self.processes)*len(self.systs)*2)
        pbar = tqdm(total=n_histograms, desc='Histograms read', dynamic_ncols=True, leave=False)

        for combine_channel in self.combine_channels.values():
            substring = '-'.join([self.tagger.name, self.wp.name, combine_channel['pt_bin'].name, combine_channel['year'], self.tagger.fit_variable])
            infile_path = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel', combine_channel['year'], 'combine', self.tagger.name, substring, 'Templates-'+substring+'.root')
            infile = root.TFile.Open(infile_path, 'READ')
            infile_folder = '_'.join(['Main', combine_channel['region'], combine_channel['channel']])
            target_folder = outfile.mkdir(combine_channel['name'])
            target_folder.cd()
            #________________________________
            # Data
            data_obs = infile.Get(infile_folder+'/data_obs')
            data_obs.Write('data_obs')
            pbar.update(1)
            combine_channel['binning'] = get_binning(data_obs)
            combine_channel.setdefault('yields', OrderedDict())['data_obs'] = data_obs.Integral()
            combine_channel.setdefault('skip', OrderedDict())['data_obs'] = False
            #________________________________
            # Individual processes with individual uncertainties
            for process in self.processes.values():
                hist = infile.Get(infile_folder+'/'+process['name'])
                combine_channel.setdefault('yields', OrderedDict()).setdefault(process['name'], OrderedDict())['nominal'] = hist.Integral()
                skip = False
                skip = (skip or hist.Integral() <= 0)
                fit_name = self.get_fit_name(process, combine_channel)
                hist.Write(fit_name)
                pbar.update(1)
                for syst in self.systs.values():
                    hist_up = infile.Get(infile_folder+'/'+process['name']+'_'+syst.combine_name+'Up')
                    combine_channel.setdefault('yields', OrderedDict()).setdefault(process['name'], OrderedDict())[syst.combine_name+'Up'] = hist_up.Integral()
                    skip = (skip or hist_up.Integral() <= 0)
                    hist_up.Write(fit_name+'_'+syst.combine_name+'Up')
                    pbar.update(1)
                    hist_down = infile.Get(infile_folder+'/'+process['name']+'_'+syst.combine_name+'Down')
                    combine_channel.setdefault('yields', OrderedDict()).setdefault(process['name'], OrderedDict())[syst.combine_name+'Down'] = hist_down.Integral()
                    skip = (skip or hist_down.Integral() <= 0)
                    hist_down.Write(fit_name+'_'+syst.combine_name+'Down')
                    pbar.update(1)
                combine_channel.setdefault('skip', OrderedDict())[process['name']] = skip
            infile.Close()
        pbar.close()
        outfile.Close()
        print('Wrote', self.rootfile_path)
        return self.rootfile_path



    def write_datacard(self):

        '''
        Needs to be run after `self.write_rootfile()`
        '''

        print('Writing datacard...')

        with open(self.datacard_path, 'w') as file:
            file.write('imax *\n')
            file.write('jmax *\n')
            file.write('kmax *\n')
            file.write('-'*42+'\n')
            #________________________________
            file.write('shapes * * '+os.path.basename(self.rootfile_path)+' $CHANNEL/$PROCESS $CHANNEL/$PROCESS_$SYSTEMATIC\n')
            file.write('-'*42+'\n')
            #________________________________
            file.write('bin '+' '.join(self.combine_channels.keys())+'\n')
            file.write('observation'+' -1'*len(self.combine_channels)+'\n')
            file.write('-'*42+'\n')
            #________________________________
            newline1 = 'bin'
            newline2 = 'process'
            newline3 = 'process'
            newline4 = 'rate'
            for combine_channel in self.combine_channels.values():
                for process in self.processes.values():
                    if combine_channel['skip'][process['name']] == True:
                        continue
                    newline1 += ' '+combine_channel['name']
                    newline2 += ' '+self.get_fit_name(process, combine_channel)
                    newline3 += ' '+str(self.get_process_id(process, combine_channel))
                    newline4 += ' -1'
            file.write(newline1+'\n')
            file.write(newline2+'\n')
            file.write(newline3+'\n')
            file.write(newline4+'\n')
            file.write('-'*42+'\n')
            #________________________________
            # Luminosity
            for year in self.years:
                newline_lumi = 'lumi_'+year+' lnN'
                for combine_channel in self.combine_channels.values():
                    for process in self.processes.values():
                        if combine_channel['skip'][process['name']] == True:
                            continue
                        newline_lumi += ' '
                        if year == combine_channel['year']:
                            newline_lumi += str(1+_YEARS.get(year).get('lumi_unc'))
                        else:
                            newline_lumi += '-'
                file.write(newline_lumi+'\n')
            #________________________________
            # Cross sections / normalization uncertainty
            for base_process, xsec_uncert in self.base_processes.items():
                newline_xsec = 'xsec_'+base_process+' lnN'
                for combine_channel in self.combine_channels.values():
                    for process in self.processes.values():
                        if combine_channel['skip'][process['name']] == True:
                            continue
                        newline_xsec += ' '
                        if base_process == process['base_process']:
                            newline_xsec += str(1+xsec_uncert)
                        else:
                            newline_xsec += '-'
                file.write(newline_xsec+'\n')
            #________________________________
            # Other systematics

            #________________________________
            # MC statistics
            # https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part2/bin-wise-stats.html
            file.write('* autoMCStats 0 1 1\n') # event-threshold=0, include-signal=1, hist-mode=1
        print('Wrote', self.datacard_path)
        return self.datacard_path



    def combine_task(self, task_name, command, run):
        if run:
            time_start = timer()
            os.system('mkdir -p '+self.logdir)
            logfile_out = open(os.path.join(self.logdir, task_name+'.log.out'), 'w')
            logfile_err = open(os.path.join(self.logdir, task_name+'.log.err'), 'w')
            p = subprocess.Popen((command), shell=True, cwd=self.workdir, stdout=logfile_out, stderr=logfile_err)
            p.wait()
            logfile_out.close()
            logfile_err.close()
            print('  Runtime: {:.1f} sec'.format(timer() - time_start)) # in seconds
            return p.returncode # int
        else:
            return command # string



    def create_workspace(self, run=True):
        print('Creating workspace')
        command = 'text2workspace.py \\'
        command += self.datacard_path+' \\'
        command += '-o '+self.workspace_path+' \\'
        command += '-P HiggsAnalysis.CombinedLimit.my_combine_physics_model:lttModel \\'
        command += '--PO sf_naming_scheme='+self.sf_naming_scheme+' \\'
        command += '--PO sf_range=[1,0.2,2.0] \\'
        # command += '--PO merge_scenarios='+','.join(self.merge_scenarios.keys())+' \\'
        command += '--PO merge_scenarios='+','.join([k+':'+v for k,v in self.merge_scenarios.items()])+' \\'
        command += '--PO years='+','.join(self.years)+' \\'
        command += '--PO pt_bin_names='+','.join([pt_bin.name for pt_bin in self.pt_bins.values()])+' \\'
        return self.combine_task('create_workspace', command, run)



    def generate_toys(self, run=True):
        print('Creating toys')
        command = 'combine -M GenerateOnly \\'
        command += '-t -1 \\'
        command += '--setParameters '+','.join([x+'=1.' for x in self.pois.keys()])+' \\'
        command += '--freezeParameters '+','.join([x for x in self.pois.keys()])+' \\'
        command += '--saveToys \\'
        command += '-n _toy \\'
        command += self.workspace_path+' \\'
        return self.combine_task('generate_toys', command, run)



    def multidimfit(self, observed=False, freezeSyst=False, run=True):
        print('Performing maximum-likelihood fit (MultiDimFit):', expobs(observed, False), ('with frozen systematics' if freezeSyst else ''))
        command = 'combine -M MultiDimFit \\'
        command += '-v 2 \\' # more verbosity
        command += '--cminSingleNuisFit \\'
        command += '--cminFallbackAlgo Minuit2,Simplex,0:0.1 \\' # fallback algorithm if default fails; see https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/runningthetool/#generic-minimizer-options
        if freezeSyst:
            datacard_path = os.path.join(self.workdir, 'higgsCombine_'+expobs(observed)+'.MultiDimFit.mH120.root')
        else:
            datacard_path = self.workspace_path
        if not os.path.isfile(datacard_path):
            print('Warning in CombineTask.multidimfit():', datacard_path, 'does not exist')
        command += '--datacard '+datacard_path+' \\'
        command += '-n _'+expobs(observed)+('_freezeSyst' if freezeSyst else '')+' \\'
        if not observed:
            command += '-t -1 \\'
            command += '--toysFile '+os.path.join(self.workdir, 'higgsCombine_toy.GenerateOnly.mH120.123456.root')+' \\'
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
        print('Calculating prefit and postfit shapes:', expobs(observed, False))
        command = 'PostFitShapesFromWorkspace \\'
        command += '-w '+self.workspace_path+' \\'
        command += '-o '+(self.prepostfitshapes_paths.get(expobs(observed)))+' \\'
        command += '-f multidimfit_'+expobs(observed)+'.root:fit_mdf \\'
        command += '--postfit \\'
        command += '--sampling \\'
        command += '--print \\'
        command += '--total-shapes \\'
        command += '--covariance \\'
        task_name = '_'.join(['prepostfitshapes', expobs(observed)])
        return self.combine_task(task_name, command, run)



    def prepostfitshapes_for_plots(self, observed=False):
        print('Writing correctly binned prefit and postfit shapes:', expobs(observed, False))
        infile = root.TFile.Open(self.prepostfitshapes_paths.get(expobs(observed)), 'READ')
        outfile_paths = list()
        for prepost in ['prefit', 'postfit']:
            outfile_path = self.prepostfitshapes_for_plots_paths.get(expobs(observed)).get(prepost)
            outfile_paths.append(outfile_path)
            outfile = root.TFile.Open(outfile_path, 'RECREATE')
            for combine_channel_name, combine_channel in self.combine_channels.items():
                target_folder = outfile.mkdir(combine_channel_name)
                target_folder.cd()
                infile_folder = '_'.join([combine_channel_name, prepost])
                #________________________________
                # Data
                data_obs_raw = infile.Get(infile_folder+'/data_obs')
                data_obs = return_rebinned(data_obs_raw, combine_channel.get('binning'))
                data_obs.Write('data_obs')
                #________________________________
                # Total prediction stack with full uncertainty
                full_stack_raw = infile.Get(infile_folder+'/TotalProcs')
                full_stack = return_rebinned(full_stack_raw, combine_channel.get('binning'))
                full_stack.Write('TotalProcs')
                #________________________________
                # Individual processes with full uncertainties
                # for combine_channel in self.combine_channels.values():
                for process in self.processes.values():
                    if combine_channel['skip'][process['name']] == True:
                        # fill in empty histogram
                        process_hist = empty_histogram(combine_channel.get('binning'))
                    else:
                        process_hist_raw = infile.Get(infile_folder+'/'+self.get_fit_name(process, combine_channel))
                        process_hist = return_rebinned(process_hist_raw, combine_channel.get('binning'))
                    hist_name = process['name']
                    process_hist.Write(hist_name)
            outfile.Close()
            print('Wrote', outfile_path)
        infile.Close()
        return outfile_paths



    def plot_prepostfitshapes(self, observed=False, prepostfit='prefitCombine', do_legend=True):

        '''prepostfit is either "prefitCombine" or "postfitCombine"'''

        print('Plotting pre- and post-fit distributions:', expobs(observed, False))

        if self.tagger.name.startswith('ak8_'):
            probejetalgo = 'AK8'
        elif self.tagger.name.startswith('hotvr_'):
            probejetalgo = 'HOTVR'
        _, x_axis_title, x_axis_unit, logy, leg_offset_x, leg_offset_y = get_variable_binning_xlabel_xunit(variable_name=self.tagger.fit_variable.split('_'+probejetalgo+'_')[-1], tagger_name=self.tagger.name, fit_variable=True)

        mscSplitting = self.mscSplitting
        processes_Plotter = self.processes_Plotter

        for combine_channel in self.combine_channels.values():

            year = combine_channel['year']
            region = combine_channel['region']
            channel = combine_channel['channel']
            pt_bin = combine_channel['pt_bin']

            nice = NiceStackWithRatio(
                infile_path = self.prepostfitshapes_for_plots_paths.get(expobs(observed)).get(prepostfit.replace('Combine', '')),
                infile_directory = combine_channel['name'], # the directory within the ROOT file
                x_axis_title = x_axis_title,
                x_axis_unit = x_axis_unit,
                prepostfit = prepostfit,
                processes = processes_Plotter.values(),
                # syst_names = systs_Plotter,
                lumi_unc = _YEARS.get(year).get('lumi_unc'),
                # divide_by_bin_width = False,
                data_name = 'data_obs',
                text_prelim = 'Private Work',
                # text_top_left = _YEARS.get(year).get('long_name'),
                text_top_left = 'T&P '+('e' if channel == 'ele' else '#mu')+'+jets, UL '+_YEARS.get(year).get('year'),
                text_top_right = _YEARS.get(year).get('lumi_fb_display')+' fb^{#minus1} (13 TeV)',
                # nostack = True,
                logy = logy,
            )

            task_name = '-'.join([self.tagger.name, self.wp.name, pt_bin.name, year, self.tagger.fit_variable])
            region_substring = '_'.join(['Main', region, channel])
            plotName = 'Plot-'+prepostfit+'-'+task_name+'-'+region_substring+'-'+mscSplitting+('-withLeg' if do_legend else '-noLeg')+'.pdf'
            plotDir = self.control_plots_dir

            nice.plot()

            #________________
            # Some text labels

            nice.canvas.cd()

            pt_text_string = self.tagger.label.replace('{WP_VALUE}', '{}'.format(self.wp.get_cut_value(year)))+' [#bf{'+region+'}]'
            pt_text_string2 = ''
            if pt_bin.max_set:
                pt_text_string2 += '#it{p}_{T} #in ({pt_min}, {pt_max}] GeV, |#eta| < 2.5'.replace('{pt_min}', '{}'.format(pt_bin.var_min)).replace('{pt_max}', '{}'.format(pt_bin.var_max))
            else:
                pt_text_string2 += '#it{p}_{T} > {pt_min} GeV, |#eta| < 2.5'.replace('{pt_min}', '{}'.format(pt_bin.var_min))

            tlatex_pt = root.TLatex(nice.coord.graph_to_pad_x(0.95), nice.coord.graph_to_pad_y(0.84), pt_text_string)
            tlatex_pt.SetTextAlign(31) # left top
            tlatex_pt.SetTextFont(42)
            tlatex_pt.SetTextSize(0.025)
            tlatex_pt.SetNDC()
            tlatex_pt.Draw()

            tlatex_pt2 = root.TLatex(nice.coord.graph_to_pad_x(0.95), nice.coord.graph_to_pad_y(0.78), pt_text_string2)
            tlatex_pt2.SetTextAlign(31) # left top
            tlatex_pt2.SetTextFont(42)
            tlatex_pt2.SetTextSize(0.025)
            tlatex_pt2.SetNDC()
            tlatex_pt2.Draw()

            tlatex_probejetalgo = root.TLatex(nice.coord.graph_to_pad_x(0.95), nice.coord.graph_to_pad_y(0.95), probejetalgo+' PUPPI')
            tlatex_probejetalgo.SetTextAlign(33) # right top
            tlatex_probejetalgo.SetTextFont(62)
            tlatex_probejetalgo.SetTextSize(0.05)
            tlatex_probejetalgo.SetTextColor(root.kGray+1)
            tlatex_probejetalgo.SetNDC()
            tlatex_probejetalgo.Draw()

            #________________
            # Legend

            if do_legend:

                legend = root.TLegend(nice.coord.graph_to_pad_x(0.45+leg_offset_x), nice.coord.graph_to_pad_y(0.5+leg_offset_y), nice.coord.graph_to_pad_x(0.7+leg_offset_x), nice.coord.graph_to_pad_y(0.73+leg_offset_y))
                if mscSplitting == 'mscTop3' or mscSplitting == 'mscW3':
                    legend.SetHeader('')
                legend.AddEntry(nice.data_hist, 'Data', 'ep')
                legend.AddEntry(nice.stack.GetStack().At(processes_Plotter['VJetsAndVV'].index), processes_Plotter['VJetsAndVV'].legend, 'f')
                legend.AddEntry(nice.stack.GetStack().At(processes_Plotter['QCD'].index), processes_Plotter['QCD'].legend, 'f')
                legend.AddEntry(nice.stack_unc, 'Uncertainty', 'f')
                legend.SetTextSize(0.025)
                legend.SetBorderSize(0)
                legend.SetFillStyle(0)
                legend.Draw()

                legend2 = root.TLegend(nice.coord.graph_to_pad_x(0.45+0.22+leg_offset_x), nice.coord.graph_to_pad_y(0.5+leg_offset_y), nice.coord.graph_to_pad_x(0.7+0.22+leg_offset_x), nice.coord.graph_to_pad_y(0.73+leg_offset_y))
                legend2.SetHeader('t#bar{t} / single t categories:')
                if mscSplitting == 'mscTop2':
                    legend2.AddEntry(nice.stack.GetStack().At(processes_Plotter['TTbar__MSc_FullyMerged'].index), '/', 'f')
                    legend2.AddEntry(nice.stack.GetStack().At(processes_Plotter['TTbar__MSc_YllufMerged'].index), '/', 'f')
                    legend2.AddEntry(nice.stack.GetStack().At(processes_Plotter['TTbar__MSc_Background'].index), '/', 'f')
                elif mscSplitting == 'mscTop3':
                    legend2.AddEntry(nice.stack.GetStack().At(processes_Plotter['TTbar__MSc_FullyMerged'].index), '/', 'f')
                    legend2.AddEntry(nice.stack.GetStack().At(processes_Plotter['TTbar__MSc_SemiMerged'].index), '/', 'f')
                    legend2.AddEntry(nice.stack.GetStack().At(processes_Plotter['TTbar__MSc_NotMerged'].index), '/', 'f')
                    legend2.AddEntry(nice.stack.GetStack().At(processes_Plotter['TTbar__MSc_Background'].index), '/', 'f')
                elif mscSplitting == 'mscW3':
                    legend2.AddEntry(nice.stack.GetStack().At(processes_Plotter['TTbar__MSc_FullyMerged'].index), '/', 'f')
                    legend2.AddEntry(nice.stack.GetStack().At(processes_Plotter['TTbar__MSc_WMerged'].index), '/', 'f')
                    legend2.AddEntry(nice.stack.GetStack().At(processes_Plotter['TTbar__MSc_NotTopOrWMerged'].index), '/', 'f')
                    legend2.AddEntry(nice.stack.GetStack().At(processes_Plotter['TTbar__MSc_Background'].index), '/', 'f')
                legend2.SetTextSize(0.025)
                legend2.SetBorderSize(0)
                legend2.SetFillStyle(0)
                legend2.Draw()

                legend3 = root.TLegend(nice.coord.graph_to_pad_x(0.45+0.22+0.067+leg_offset_x), nice.coord.graph_to_pad_y(0.5+leg_offset_y), nice.coord.graph_to_pad_x(0.7+0.22+0.067+leg_offset_x), nice.coord.graph_to_pad_y(0.73+leg_offset_y))
                legend3.SetHeader('')
                if mscSplitting == 'mscTop2':
                    legend3.AddEntry(nice.stack.GetStack().At(processes_Plotter['ST__MSc_FullyMerged'].index), processes_Plotter['ST__MSc_FullyMerged'].legend, 'f')
                    legend3.AddEntry(nice.stack.GetStack().At(processes_Plotter['ST__MSc_YllufMerged'].index), processes_Plotter['ST__MSc_YllufMerged'].legend, 'f')
                    legend3.AddEntry(nice.stack.GetStack().At(processes_Plotter['ST__MSc_Background'].index), processes_Plotter['ST__MSc_Background'].legend, 'f')
                elif mscSplitting == 'mscTop3':
                    legend3.AddEntry(nice.stack.GetStack().At(processes_Plotter['ST__MSc_FullyMerged'].index), processes_Plotter['ST__MSc_FullyMerged'].legend, 'f')
                    legend3.AddEntry(nice.stack.GetStack().At(processes_Plotter['ST__MSc_SemiMerged'].index), processes_Plotter['ST__MSc_SemiMerged'].legend, 'f')
                    legend3.AddEntry(nice.stack.GetStack().At(processes_Plotter['ST__MSc_NotMerged'].index), processes_Plotter['ST__MSc_NotMerged'].legend, 'f')
                    legend3.AddEntry(nice.stack.GetStack().At(processes_Plotter['ST__MSc_Background'].index), processes_Plotter['ST__MSc_Background'].legend, 'f')
                elif mscSplitting == 'mscW3':
                    legend3.AddEntry(nice.stack.GetStack().At(processes_Plotter['ST__MSc_FullyMerged'].index), processes_Plotter['ST__MSc_FullyMerged'].legend, 'f')
                    legend3.AddEntry(nice.stack.GetStack().At(processes_Plotter['ST__MSc_WMerged'].index), processes_Plotter['ST__MSc_WMerged'].legend, 'f')
                    legend3.AddEntry(nice.stack.GetStack().At(processes_Plotter['ST__MSc_NotTopOrWMerged'].index), processes_Plotter['ST__MSc_NotTopOrWMerged'].legend, 'f')
                    legend3.AddEntry(nice.stack.GetStack().At(processes_Plotter['ST__MSc_Background'].index), processes_Plotter['ST__MSc_Background'].legend, 'f')
                legend3.SetTextSize(0.025)
                legend3.SetBorderSize(0)
                legend3.SetFillStyle(0)
                legend3.Draw()

            pt_text_string3 = '#minus postfit #minus' if prepostfit=='postfitCombine' else '#minus prefit #minus'
            tlatex_pt3 = root.TLatex(nice.coord.graph_to_pad_x(0.95), nice.coord.graph_to_pad_y(0.4), pt_text_string3)
            tlatex_pt3.SetTextAlign(31) # left top
            tlatex_pt3.SetTextFont(72)
            tlatex_pt3.SetTextSize(nice.text_size)
            tlatex_pt3.SetNDC()
            tlatex_pt3.Draw()

            nice.save_plot(plotName, plotDir)



    def write_scale_factor_file(self, observed=False):
        print('Writing scale factors to disk')
        # Check out this link to learn how to extract post-fit nuisance parameter values from RooFitResult:
        # https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/102x/src/MultiDimFit.cc#L434
        mdf_file = root.TFile.Open(os.path.join(self.workdir, 'multidimfit_'+expobs(observed)+'.root'), 'READ')
        mdf_file_stat = root.TFile.Open(os.path.join(self.workdir, 'multidimfit_'+expobs(observed)+'_freezeSyst.root'), 'READ')
        fit_result = mdf_file.Get('fit_mdf') # RooFitResult
        fit_result_stat = mdf_file_stat.Get('fit_mdf') # RooFitResult
        for poi_name in self.pois.keys():
            rf = fit_result.floatParsFinal().find(poi_name)
            rf_stat = fit_result_stat.floatParsFinal().find(poi_name)
            bestFitVal = rf.getVal()
            if abs(bestFitVal - rf_stat.getVal()) > 0.01: # check if first fit and syst-frozen fit lead to approximately same central value
                print('bestFitVal regular:', bestFitVal)
                print('bestFitVal frozenSyst:', rf_stat.getVal())
                print('difference:', bestFitVal - rf_stat.getVal())
                print('Warning: best-fit scale factor central values of regular fit and fit with frozen systematics differ from each other')
            #________________________________
            hiErr = +(rf.getMax('err68') - bestFitVal if rf.hasRange('err68') else rf.getAsymErrorHi())
            loErr = -(rf.getMin('err68') - bestFitVal if rf.hasRange('err68') else rf.getAsymErrorLo())
            maxError = max(max(hiErr, loErr), rf.getError())
            if abs(hiErr) < 0.001*maxError:
                # print " Warning - No valid high-error found, will report difference to maximum of range for : ", rf.GetName()
                # hiErr = -bestFitVal + rf.getMax()
                print('Warning: no valid asymmetric high-error found, will use getError() result')
                hiErr = rf.getError()
            if abs(loErr) < 0.001*maxError:
                # print " Warning - No valid low-error found, will report difference to minimum of range for : ", rf.GetName()
                # loErr = +bestFitVal - rf.getMin()
                print('Warning: no valid asymmetric low-error found, will use getError() result')
                loErr = rf.getError()
            #________________________________
            hiErr_stat = +(rf_stat.getMax('err68') - bestFitVal if rf_stat.hasRange('err68') else rf_stat.getAsymErrorHi())
            loErr_stat = -(rf_stat.getMin('err68') - bestFitVal if rf_stat.hasRange('err68') else rf_stat.getAsymErrorLo())
            maxError_stat = max(max(hiErr_stat, loErr_stat), rf_stat.getError())
            if abs(hiErr_stat) < 0.001*maxError_stat:
                # print " Warning - No valid high-error found, will report difference to maximum of range for : ", rf.GetName()
                # hiErr_stat = -bestFitVal + rf_stat.getMax()
                print('Warning: no valid asymmetric high-error found, will use getError() result')
                hiErr_stat = rf.getError()
            if abs(loErr_stat) < 0.001*maxError_stat:
                # print " Warning - No valid low-error found, will report difference to minimum of range for : ", rf.GetName()
                # loErr_stat = +bestFitVal - rf_stat.getMin()
                print('Warning: no valid asymmetric low-error found, will use getError() result')
                loErr_stat = rf.getError()
            #________________________________
            self.pois.get(poi_name)[expobs(observed)] = OrderedDict()
            self.pois.get(poi_name).get(expobs(observed))['bestFitVal'] = bestFitVal
            self.pois.get(poi_name).get(expobs(observed))['loErr'] = loErr
            self.pois.get(poi_name).get(expobs(observed))['hiErr'] = hiErr
            self.pois.get(poi_name).get(expobs(observed))['loErr_stat'] = loErr_stat
            self.pois.get(poi_name).get(expobs(observed))['hiErr_stat'] = hiErr_stat
            loErr_syst = loErr*loErr - loErr_stat*loErr_stat
            if loErr_syst >= 0:
                loErr_syst = root.TMath.Sqrt(loErr_syst)
            else:
                print('Warning: lower total error smaller than lower statistical component; something went wrong! Lower systematic component set to zero')
                loErr_syst = 0
            hiErr_syst = hiErr*hiErr - hiErr_stat*hiErr_stat
            if hiErr_syst >= 0:
                hiErr_syst = root.TMath.Sqrt(hiErr_syst)
            else:
                print('Warning: upper total error smaller than upper statistical component; something went wrong! Upper systematic component set to zero')
                hiErr_syst = 0
            self.pois.get(poi_name).get(expobs(observed))['loErr_syst'] = loErr_syst
            self.pois.get(poi_name).get(expobs(observed))['hiErr_syst'] = hiErr_syst
            print("  {:55s}:  {:.3f}  -{:.3f} / +{:.3f} (tot.)  [ -{:.3f} / +{:.3f} (stat.)  |  -{:.3f} / +{:.3f} (syst.) ]".format(poi_name, bestFitVal, loErr, hiErr, loErr_stat, hiErr_stat, loErr_syst, hiErr_syst))
        sf_file_path = os.path.join(self.workdir, 'scale_factors_'+expobs(observed)+'.root')
        if observed:
            self.sf_file_path_obs = sf_file_path
        else:
            self.sf_file_path_exp = sf_file_path
        sf_file = root.TFile.Open(sf_file_path, 'RECREATE')
        sf_file.cd()
        for year in self.years:
            for msc in self.merge_scenarios:
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
                for poi_name, poi_info in self.pois.items():
                    # if msc not in poi_name:
                    if msc != poi_info['msc']:
                        continue
                    if year != poi_info['year']:
                        continue
                    pt_bin = poi_info['pt_bin']
                    pt_min = pt_bin.var_min
                    pt_max = pt_bin.var_max
                    pt_max = min(pt_max, 1000.)
                    if pt_min >= pt_max:
                        sys.exit('pt_min > pt_max. Check definition of this pt interval')
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
                sf_graph.Write(msc+'_'+year+'_tot')
                sf_graph_stat.Write(msc+'_'+year+'_stat')
                sf_graph_syst.Write(msc+'_'+year+'_syst')
        print('Wrote', sf_file_path)






if __name__=='__main__':

    x = ScaleFactorFits(
        tagger_name = 'ak8_t__tau',
        wp = 0,
        total_range = True,
        # mode='Hybrid',
    )

    # Need to run this always:
    x.write_rootfile()
    x.write_datacard()
    x.create_workspace()

    # Expected:
    # x.generate_toys()
    # x.multidimfit(observed=False)
    # x.multidimfit(observed=False, freezeSyst=True)
    # x.write_scale_factor_file(observed=False)
    # x.prepostfitshapes(observed=False)
    # x.prepostfitshapes_for_plots(observed=False)

    # Observed:
    # x.multidimfit(observed=True)
    # x.multidimfit(observed=True, freezeSyst=True)
    # x.write_scale_factor_file(observed=True)
    # x.prepostfitshapes(observed=True)
    # x.prepostfitshapes_for_plots(observed=True)
    x.plot_prepostfitshapes(observed=True, prepostfit='prefitCombine')
    x.plot_prepostfitshapes(observed=True, prepostfit='postfitCombine')
