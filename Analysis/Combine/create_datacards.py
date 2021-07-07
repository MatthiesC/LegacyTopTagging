import os
import sys
from collections import OrderedDict
import itertools
import ROOT as root
import subprocess

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

base_processes = OrderedDict([
    ('TTbar', 0.10),
    ('ST', 0.20),
    ('WJets', 0.20),
    ('DYJetsAndDiboson', 0.20),
    ('QCD', 1.00),
])

processes = OrderedDict([
    ('TTbar_FullyMerged', {
        'id': -101,
        'base_process': 'TTbar',
    }),
    ('TTbar_SemiMerged', {
        'id': -102,
        'base_process': 'TTbar',
    }),
    # ('TTbar_WMerged', {
    #     'id': -103,
    #     'base_process': 'TTbar',
    # }),
    # ('TTbar_QBMerged', {
    #     'id': -104,
    #     'base_process': 'TTbar',
    # }),
    ('TTbar_BkgOrNotMerged', {
        'id': -105,
        'base_process': 'TTbar',
    }),
    ('ST', {
        'id': 200,
        'base_process': 'ST',
    }),
    # ('ST_FullyMerged', {
    #     'id': 201,
    #     'base_process': 'ST',
    # }),
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
    'tageff3prong',
    'tageff2prong',
    'tageff1prong',
]

def expobs(observed, short=True):
    if short:
        return ('obs' if observed else 'exp')
    else:
        return ('observed' if observed else 'expected')

class CombineTask:

    def __init__(self, year, algo, pt_bin, jet_version, wp):
        self.Year = year
        self.Name = '_'.join([algo, pt_bin, jet_version, wp])
        print 'New CombineTask:', self.Name
        self.WorkDir = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel', year, 'combine', algo, self.Name)
        self.LogDir = os.path.join(self.WorkDir, 'log')
        if not os.path.isdir(self.WorkDir): sys.exit('Work directory does not exist. Abort.')
        self.InputRootFilePath = os.path.join(self.WorkDir, '_'.join([self.Name, '_ProbeJetHists', algo])+'.root')
        if not os.path.isfile(self.InputRootFilePath): sys.exit('ROOT file with input histograms does not exist. Abort.')
        self.DataCardRootFilePath = self.InputRootFilePath.replace('.root', '_forCombine.root')
        self.DataCardPath = os.path.join(self.WorkDir, '_'.join([self.Name, 'datacard.txt']))
        self.WorkSpacePath = os.path.join(self.WorkDir, '_'.join([self.Name, 'workspace.root']))
        self.PrePostFitShapesExpPath = os.path.join(self.WorkDir, '_'.join([self.Name, 'prepostfitshapes_exp.root']))
        self.PrePostFitShapesObsPath = os.path.join(self.WorkDir, '_'.join([self.Name, 'prepostfitshapes_obs.root']))
        self.Variable = 'mSD' if algo=='AK8' else 'mass' if algo=='HOTVR' else None
        self.Channels = ['_'.join([pt_bin, jet_version, wp, x, self.Variable]) for x in regions]
        self.Processes = OrderedDict()
        for process in processes.keys():
            process_dict = processes.get(process)
            process_dict['prior_name'] = process
            if processes.get(process).get('id') <= 0: # if process is a signal
                new_process_name = '_'.join([process, self.Name])
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
        self.ImpactPlotPaths = list()

    def write_rootfile(self):
        infile = root.TFile.Open(self.InputRootFilePath, 'READ')
        outfile = root.TFile.Open(self.DataCardRootFilePath, 'RECREATE')
        for channel in self.Channels:
            target_folder = outfile.mkdir(channel)
            target_folder.cd()
            infile.Get(channel+'/data_obs').Write('data_obs')
            for process in self.Processes.keys():
                prior_process_name = self.Processes.get(process).get('prior_name')
                infile.Get(channel+'/'+prior_process_name).Write(process)
                for systematic in systematics:
                    infile.Get(channel+'/'+prior_process_name+'_'+systematic+'Up').Write(process+'_'+systematic+'Up')
                    infile.Get(channel+'/'+prior_process_name+'_'+systematic+'Down').Write(process+'_'+systematic+'Down')
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
            file.write('bin '+' '.join(self.Channels)+'\n')
            file.write('observation'+' -1'*len(self.Channels)+'\n')
            file.write('-'*42+'\n')
            #________________________________
            newline1 = 'bin'
            newline2 = 'process'
            newline3 = 'process'
            newline4 = 'rate'
            for channel in self.Channels:
                for process in self.Processes.keys():
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
            # Cross sections
            for base_process in base_processes.keys():
                newline1 = '_'.join(['xsec', base_process])+' lnN'
                for channel in self.Channels:
                    for process in self.Processes.keys():
                        if self.Processes.get(process).get('base_process') == base_process:
                            newline1 += ' '+str(1+base_processes.get(base_process))
                        else:
                            newline1 += ' -'
                file.write(newline1+'\n')
            # Other systematics
            for systematic in systematics:
                newline1 = systematic+' shape'
                for channel in self.Channels:
                    for process in self.Processes.keys():
                        newline1 += ' 1.0'
                file.write(newline1+'\n')
            # Rate parameters
            for process in self.Processes.keys():
                file.write(self.Processes.get(process).get('rateParam')+' rateParam * '+process+' 1.0 [0.2,5.0] \n')
            # MC statistics
            file.write('* autoMCStats 0 1 1\n') # event-threshold=0, include-signal=1, hist-mode=1
        print 'Wrote', self.DataCardPath
        return self.DataCardPath

    def combine_task(self, task_name, command, run):
        if run:
            os.system('mkdir -p '+self.LogDir)
            logfile_out = open(os.path.join(self.LogDir, task_name+'.log.out'), 'w')
            logfile_err = open(os.path.join(self.LogDir, task_name+'.log.err'), 'w')
            p = subprocess.Popen((command), shell=True, cwd=self.WorkDir, stdout=logfile_out, stderr=logfile_err)
            p.wait()
            logfile_out.close()
            logfile_err.close()
            return p.returncode # int
        else:
            return command # string

    def create_workspace(self, run=True):
        print 'Creating workspace'
        command = 'text2workspace.py \\'
        command += self.DataCardPath+' \\'
        command += '-o '+self.WorkSpacePath+' \\'
        command += '-P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel \\'
        command += '--PO verbose \\'
        for process in self.Processes.keys():
            if self.Processes.get(process).get('is_signal'):
                command += '''--PO 'map=.*/'''+process+':'+self.Processes.get(process).get('rateParam')+'''[1,0,5]\' \\'''
        return self.combine_task('create_workspace', command, run)

    def generate_toys(self, run=True):
        print 'Creating toys'
        command = 'combine -M GenerateOnly \\'
        command += '-t -1 \\'
        command += '--setParameters '+','.join([x+'=1.' for x in self.SignalRateParams.values()])+' \\'
        command += '--freezeParameters '+','.join([x for x in self.SignalRateParams.values()])+' \\'
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
        task_name = '_'.join(['multidimfit', expobs(observed)])
        return self.combine_task(task_name, command, run)

    def prepostfitshapes(self, observed=False, run=True):
        print 'Calculating prefit and postfit shapes:', expobs(observed, False)
        command = 'PostFitShapesFromWorkspace \\'
        command += '-w '+self.WorkSpacePath+' \\'
        command += '-o '+(self.PrePostFitShapesObsPath if observed else self.PrePostFitShapesExpPath)+' \\'
        command += '-f multidimfit_'+expobs(observed)+'.root:fit_mdf \\'
        command += '--postfit \\'
        command += '--sampling \\'
        command += '--print \\'
        command += '--total-shapes \\'
        command += '--covariance \\'
        task_name = '_'.join(['prepostfitshapes', expobs(observed)])
        return self.combine_task(task_name, command, run)

    def impacts(self, observed=False, run=True):
        print 'Calculating nuisance parameter impacts:', expobs(observed, False)
        results = OrderedDict()
        #________________________________
        # Do an initial fit
        command1 = 'combineTool.py -M Impacts \\'
        command1 += '-m 0 \\' # required argument
        command1 += '-d '+self.WorkSpacePath+' \\'
        command1 += '--doInitialFit \\'
        command1 += '--robustFit 1 \\'
        command1 += '-t -1 \\'
        command1 += '--toysFile '+os.path.join(self.WorkDir, 'higgsCombine_toy.GenerateOnly.mH120.123456.root')+' \\'
        task_name1 = '_'.join(['impacts', expobs(observed), 'command1'])
        results[task_name1] = self.combine_task(task_name1, command1, run)
        #________________________________
        #
        command2 = command1.replace('--doInitialFit', '--doFits')
        task_name2 = '_'.join(['impacts', expobs(observed), 'command2'])
        results[task_name2] = self.combine_task(task_name2, command2, run)
        #________________________________
        #
        impactsJsonPath = os.path.join(self.WorkDir, 'impacts_'+expobs(observed)+'.json')
        command3 = 'combineTool.py -M Impacts \\'
        command3 += '-m 0 \\' # required argument
        command3 += '-d '+self.WorkSpacePath+' \\'
        command3 += '-o '+impactsJsonPath+' \\'
        task_name3 = '_'.join(['impacts', expobs(observed), 'command3'])
        results[task_name3] = self.combine_task(task_name3, command3, run)
        #________________________________
        #
        plot_command_base = 'plotImpacts.py \\'
        plot_command_base += '-i '+impactsJsonPath+' \\'
        plot_command_base += '--POI {0} \\'
        plot_command_base += '-o {1} \\'
        # plot_command_base += '-t '+globalRenameJsonFile+' \\' # rename parameters
        for signal_rate_param in self.SignalRateParams.values():
            print 'Plotting impacts for parameter:', signal_rate_param
            plot_path = os.path.join(self.WorkDir, impactsJsonPath.replace('.json', '_'+signal_rate_param))
            plot_command = plot_command_base.format(signal_rate_param, os.path.basename(plot_path))
            self.ImpactPlotPaths.append(plot_path+'.pdf')
            plot_task_name = '_'.join(['plot_impacts', expobs(observed), signal_rate_param])
            results[plot_task_name] = self.combine_task(plot_task_name, plot_command, run)
        return results

    def run_all(self):
        self.write_rootfile()
        self.write_datacard()
        self.create_workspace()
        # Expected:
        self.generate_toys()
        self.multidimfit()
        self.prepostfitshapes()
        self.impacts()
        # Observed:
        self.multidimfit(True)
        self.prepostfitshapes(True)
        self.impacts(True)



    # def impacts_observed():
    #     pass

# PostFitShapesFromWorkspace \
#     -w HOTVR_Pt400to480_HOTVRCuts_Standard_datacard.root \
#     -o HOTVR_Pt400to480_HOTVRCuts_Standard_shapes.root \
#     -f multidimfit_nominal_obs.root:fit_mdf \
#     --postfit \
#     --sampling \
#     --print \
#     --total-shapes \
#     --covariance

# def write_datacard(algo, pt_bin, jet_version, wp, region):
#     task_name = '_'.join([pt_bin, jet_version, wp])
#     workdir = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel/UL17/combine', algo, task_name)
#     input_root_file_name =

all_tasks = list(itertools.product(['AK8'], pt_bins.get('AK8'), jet_versions.get('AK8'), wps.get('AK8')))
all_tasks.extend(list((itertools.product(['HOTVR'], pt_bins.get('HOTVR'), jet_versions.get('HOTVR'), wps.get('HOTVR')))))

all_tasks = [OrderedDict([
    ('probejet_coll', x[0]),
    ('pt_bin', x[1]),
    ('jet_version', x[2]),
    ('wp', x[3]),
    # ('variables', set_variables(x)), ### Schreibe wirklich nur die Variablen aus, die Du brauchst!!!!! (entferne den obigen for-loop in create_input_histograms)
]) for x in all_tasks]

print all_tasks

# for task in all_tasks:
#     c = CombineTask('UL17', task.get('probejet_coll'), task.get('pt_bin'), task.get('jet_version'), task.get('wp'))
#     c.write_rootfile()
#     c.write_datacard()





c = CombineTask('UL17', 'HOTVR', 'Pt300to400', 'HOTVRCuts', 'Standard')
c.run_all()
# c.write_rootfile()
# c.write_datacard()
# c.create_workspace()
# c.generate_toys()
# c.multidimfit()
# c.prepostfitshapes()
# c.impacts()

# text2workspace.py \
#    HOTVR_Pt400to480_HOTVRCuts_Standard_datacard.txt \
#    -o HOTVR_Pt400to480_HOTVRCuts_Standard_datacard.root \
#    -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel \
#    --PO verbose \
#    --PO 'map=.*/TTbar_FullyMerged:rate_TTbar_FullyMerged[1,0,5]' \
#    --PO 'map=.*/TTbar_SemiMerged:rate_TTbar_SemiMerged[1,0,5]' \
#    --PO 'map=.*/TTbar_BkgOrNotMerged:rate_TTbar_BkgOrNotMerged[1,0,5]'
# combine -M MultiDimFit \
#     -n _nominal_obs \
#     --algo singles \
#     --datacard HOTVR_Pt400to480_HOTVRCuts_Standard_datacard.root \
#     --saveFitResult
# combine -M GenerateOnly \
#     -t -1 \
#     --setParameters rate_TTbar_FullyMerged=1.,rate_TTbar_SemiMerged=1.,rate_TTbar_BkgOrNotMerged=1. \
#     --freezeParameters rate_TTbar_FullyMerged,rate_TTbar_SemiMerged,rate_TTbar_BkgOrNotMerged \
#     --saveToys \
#     -n _toy HOTVR_Pt400to480_HOTVRCuts_Standard_datacard.root
# combine -M MultiDimFit \
#     -n _nominal_exp \
#     --algo singles \
#     -t -1 \
#     --toysFile higgsCombine_toy.GenerateOnly.mH120.123456.root \
#     --datacard HOTVR_Pt400to480_HOTVRCuts_Standard_datacard.root \
#     --saveFitResult
# PostFitShapesFromWorkspace \
#     -w HOTVR_Pt400to480_HOTVRCuts_Standard_datacard.root \
#     -o HOTVR_Pt400to480_HOTVRCuts_Standard_shapes.root \
#     -f multidimfit_nominal_obs.root:fit_mdf \
#     --postfit \
#     --sampling \
#     --print \
#     --total-shapes \
#     --covariance



# text2workspace.py \
#    HOTVR_PtCombined_HOTVRCuts_Standard_datacard.txt \
#    -o HOTVR_PtCombined_HOTVRCuts_Standard_datacard.root \
#    -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel \
#    --PO verbose \
#    --PO 'map=.*/TTbar_FullyMerged_HOTVR_Pt200to250_HOTVRCuts_Standard:rate_TTbar_FullyMerged_HOTVR_Pt200to250_HOTVRCuts_Standard[1,0,5]' \
#    --PO 'map=.*/TTbar_SemiMerged_HOTVR_Pt200to250_HOTVRCuts_Standard:rate_TTbar_SemiMerged_HOTVR_Pt200to250_HOTVRCuts_Standard[1,0,5]' \
#    --PO 'map=.*/TTbar_BkgOrNotMerged_HOTVR_Pt200to250_HOTVRCuts_Standard:rate_TTbar_BkgOrNotMerged_HOTVR_Pt200to250_HOTVRCuts_Standard[1,0,5]' \
#    --PO 'map=.*/TTbar_FullyMerged_HOTVR_Pt250to300_HOTVRCuts_Standard:rate_TTbar_FullyMerged_HOTVR_Pt250to300_HOTVRCuts_Standard[1,0,5]' \
#    --PO 'map=.*/TTbar_SemiMerged_HOTVR_Pt250to300_HOTVRCuts_Standard:rate_TTbar_SemiMerged_HOTVR_Pt250to300_HOTVRCuts_Standard[1,0,5]' \
#    --PO 'map=.*/TTbar_BkgOrNotMerged_HOTVR_Pt250to300_HOTVRCuts_Standard:rate_TTbar_BkgOrNotMerged_HOTVR_Pt250to300_HOTVRCuts_Standard[1,0,5]' \
#    --PO 'map=.*/TTbar_FullyMerged_HOTVR_Pt300to400_HOTVRCuts_Standard:rate_TTbar_FullyMerged_HOTVR_Pt300to400_HOTVRCuts_Standard[1,0,5]' \
#    --PO 'map=.*/TTbar_SemiMerged_HOTVR_Pt300to400_HOTVRCuts_Standard:rate_TTbar_SemiMerged_HOTVR_Pt300to400_HOTVRCuts_Standard[1,0,5]' \
#    --PO 'map=.*/TTbar_BkgOrNotMerged_HOTVR_Pt300to400_HOTVRCuts_Standard:rate_TTbar_BkgOrNotMerged_HOTVR_Pt300to400_HOTVRCuts_Standard[1,0,5]' \
#    --PO 'map=.*/TTbar_FullyMerged_HOTVR_Pt400to480_HOTVRCuts_Standard:rate_TTbar_FullyMerged_HOTVR_Pt400to480_HOTVRCuts_Standard[1,0,5]' \
#    --PO 'map=.*/TTbar_SemiMerged_HOTVR_Pt400to480_HOTVRCuts_Standard:rate_TTbar_SemiMerged_HOTVR_Pt400to480_HOTVRCuts_Standard[1,0,5]' \
#    --PO 'map=.*/TTbar_BkgOrNotMerged_HOTVR_Pt400to480_HOTVRCuts_Standard:rate_TTbar_BkgOrNotMerged_HOTVR_Pt400to480_HOTVRCuts_Standard[1,0,5]' \
#    --PO 'map=.*/TTbar_FullyMerged_HOTVR_Pt480to600_HOTVRCuts_Standard:rate_TTbar_FullyMerged_HOTVR_Pt480to600_HOTVRCuts_Standard[1,0,5]' \
#    --PO 'map=.*/TTbar_SemiMerged_HOTVR_Pt480to600_HOTVRCuts_Standard:rate_TTbar_SemiMerged_HOTVR_Pt480to600_HOTVRCuts_Standard[1,0,5]' \
#    --PO 'map=.*/TTbar_BkgOrNotMerged_HOTVR_Pt480to600_HOTVRCuts_Standard:rate_TTbar_BkgOrNotMerged_HOTVR_Pt480to600_HOTVRCuts_Standard[1,0,5]' \
#    --PO 'map=.*/TTbar_FullyMerged_HOTVR_Pt600toInf_HOTVRCuts_Standard:rate_TTbar_FullyMerged_HOTVR_Pt600toInf_HOTVRCuts_Standard[1,0,5]' \
#    --PO 'map=.*/TTbar_SemiMerged_HOTVR_Pt600toInf_HOTVRCuts_Standard:rate_TTbar_SemiMerged_HOTVR_Pt600toInf_HOTVRCuts_Standard[1,0,5]' \
#    --PO 'map=.*/TTbar_BkgOrNotMerged_HOTVR_Pt600toInf_HOTVRCuts_Standard:rate_TTbar_BkgOrNotMerged_HOTVR_Pt600toInf_HOTVRCuts_Standard[1,0,5]'













if(True):
    pass
