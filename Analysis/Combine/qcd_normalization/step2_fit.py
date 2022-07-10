import os
import sys
import argparse
# from tqdm import tqdm
import subprocess
from timeit import default_timer as timer
import warnings

import ROOT as root

from copy import deepcopy
import numpy as np

sys.path.append(os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/Analysis'))
from constants import _TAGGERS, _BANDS, _YEARS


def combine_task(task_name, command, logDir, workDir, run=True):
    if run:
        time_start = timer()
        os.system('mkdir -p '+logDir)
        logfile_out = open(os.path.join(logDir, task_name+'.log.out'), 'w')
        logfile_err = open(os.path.join(logDir, task_name+'.log.err'), 'w')
        p = subprocess.Popen((command), shell=True, cwd=workDir, stdout=logfile_out, stderr=logfile_err)
        p.wait()
        logfile_out.close()
        logfile_err.close()
        print '  Runtime: {:.1f} sec'.format(timer() - time_start) # in seconds
        return p.returncode # int
    else:
        return command # string

def return_rebinned(th1, binning):
    if th1.GetNbinsX() != len(binning)-1:
        warnings.warn('Cannot rebin TH1: given binning scheme does not match number of bins of TH1')
        return th1
    th1_new = root.TH1F(th1.GetName(), th1.GetTitle(), len(binning)-1, binning)
    for ibin in range(len(th1)+2):
        th1_new.SetBinContent(ibin, th1.GetBinContent(ibin))
        th1_new.SetBinError(ibin, th1.GetBinError(ibin))
    return th1_new


all_years = ['UL16preVFP', 'UL16postVFP', 'UL17', 'UL18']
all_channels = ['muo']

taggers = ['ak8_t__tau', 'hotvr_t__tau', 'ak8_w__partnet']
taggers = {k: _TAGGERS[k] for k in taggers}

if not sys.argv[1:]: sys.exit('No arguments provided. Exit.')
parser = argparse.ArgumentParser()
parser.add_argument('-y', '--years', nargs='+', choices=all_years, default=all_years)
parser.add_argument('-c', '--channels', nargs='+', choices=all_channels, default=all_channels)
parser.add_argument('-t', '--tagger', default='ak8_t__tau', choices=taggers.keys())
parser.add_argument('-o', '--use-overflow', action='store_true')
args = parser.parse_args(sys.argv[1:])

processes = [
'DATA.DATA',
'MC.QCD',
'MC.nonQCD',
]

for year in args.years:

    for channel in args.channels:

        if args.tagger.startswith('ak8_t'):
            tagger_base = 'ak8_t'
        elif args.tagger.startswith('ak8_w'):
            tagger_base = 'ak8_w'
        elif args.tagger.startswith('hotvr_t'):
            tagger_base = 'hotvr_t'

        workDir = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel', year, channel)
        inputFileName = 'qcd_normalization_histograms__'+tagger_base+'.root'
        inputDir = os.path.join(workDir, 'qcd_normalization_study')
        inputFilePath = os.path.join(inputDir, inputFileName)

        infile = root.TFile.Open(inputFilePath, 'READ')

        hists = {}
        binning = np.array([], dtype=np.float32)
        binning_set = False
        for band_k, band in _BANDS.items():
            for process in processes:
                hist = infile.Get(band_k+'/'+process)
                # hist.Rebin(2) # rebin from 40 to 20 bins
                n_bins = hist.GetNbinsX()
                if not binning_set:
                    for i_bin in range(1, n_bins+1):
                        binning = np.append(binning, [hist.GetXaxis().GetBinLowEdge(i_bin)])
                    binning = np.append(binning, [hist.GetXaxis().GetBinUpEdge(n_bins)])
                    binning_set = True
                if args.use_overflow:
                    # merge overflow into last bin:
                    hist.SetBinContent(n_bins, hist.GetBinContent(n_bins) + hist.GetBinContent(n_bins + 1))
                    hist.SetBinError(n_bins, np.sqrt(hist.GetBinError(n_bins)**2 + hist.GetBinError(n_bins + 1)**2))
                    # "delete" overflow bin:
                    hist.SetBinContent(n_bins + 1, 0)
                    hist.SetBinError(n_bins + 1, 0)
                # write to dict:
                hists.setdefault(band_k, {})[process] = hist

        # subtract non-QCD from Data in QCD sideband region:
        hist_QCD_DD = deepcopy(hists['QCD']['DATA.DATA'])
        hist_QCD_DD.Add(hists['QCD']['MC.nonQCD'], -1.)

        # write for combine:
        outFileName = inputFileName.replace('.root', '_forCombine.root')
        outputDir = inputDir
        outFilePath = os.path.join(outputDir, outFileName)
        outfile = root.TFile.Open(outFilePath, 'RECREATE')

        folder = outfile.mkdir('qcd_fit')
        folder.cd()
        hists['Main']['DATA.DATA'].Write('data_obs')
        hists['Main']['MC.nonQCD'].Write('nonQCD')
        hist_QCD_DD.Write('QCD_DD')

        outfile.Close() # need to close this before trying to create workspace, else text2workspace fails because it cannot look up still open ROOT file

        # do combine stuff:
        dcName = 'datacard_'+tagger_base+'.txt'
        dcPath = os.path.join(outputDir, dcName)
        with open(dcPath, 'w') as dc:
            dc.write('''imax *\n''')
            dc.write('''jmax *\n''')
            dc.write('''kmax *\n''')
            dc.write('''----------\n''')
            dc.write('''shapes * * {root_file} $CHANNEL/$PROCESS $CHANNEL/$PROCESS_$SYSTEMATIC\n'''.format(root_file=outFilePath))
            dc.write('''----------\n''')
            dc.write('''bin qcd_fit\n''')
            dc.write('''observation -1\n''')
            dc.write('''----------\n''')
            dc.write('''bin qcd_fit qcd_fit\n''')
            dc.write('''process nonQCD QCD_DD\n''')
            # dc.write('''process -2 -1\n''')
            dc.write('''process 2 -1\n''')  # FIXME
            dc.write('''rate -1 -1\n''')
            dc.write('''----------\n''')
            dc.write('''lumi lnN {0} {0}\n'''.format(str(_YEARS.get(year)['lumi_unc'] + 1)))
            dc.write('''norm_nonQCD lnN 1.2 -\n''') # 20% should be okay, considering that 80% of nonQCD is TTbar
            # dc.write('''norm_QCD_DD lnN - 2.0\n''') # 100% normalization uncertainty on DD template
            dc.write('''* autoMCStats 0 1 1\n''')

        # Create workspace
        print 'Create Workspace'
        wsName = 'workspace_'+tagger_base+'.root'
        wsPath = os.path.join(outputDir, wsName)
        logDir = os.path.join(outputDir, 'log_'+tagger_base)
        print logDir
        command = 'text2workspace.py \\'
        command += dcPath+' \\'
        command += '-o '+wsPath+' \\'
        command += '-P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel \\'
        # command += """--PO 'map=qcd_fit/nonQCD:r_nonQCD[1,0,5]' --PO 'map=qcd_fit/QCD_DD:r_QCD_DD[1,0,100]' \\"""
        command += """--PO 'map=qcd_fit/QCD_DD:r_QCD_DD[1,0,100]' \\"""  # FIXME
        command += '--X-allow-no-background \\'
        combine_task('create_workspace', command, logDir, outputDir)

        # Perform fit
        print 'MultiDimFit'
        command = 'combine -M MultiDimFit \\'
        command += '-v 2 \\' # more verbosity
        # command += '--cminSingleNuisFit \\'
        # command += '--cminFallbackAlgo Minuit2,Simplex,0:0.1 \\' # fallback algorithm if default fails; see https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/runningthetool/#generic-minimizer-options
        command += '--datacard '+wsPath+' \\'
        command += '--algo singles \\'
        command += '--saveFitResult \\'
        command += '--saveWorkspace \\'
        command += '-n _'+tagger_base+' \\'
        # command += '--fastScan \\'
        # command += '--robustFit=1 \\'
        # command += '--robustHesse=1 \\'
        command += '--cminDefaultMinimizerStrategy=0 \\'
        combine_task('multidimfit', command, logDir, outputDir)

        # pre/psost-fit shapes
        outfileName_prepostfitshapes = 'prepostfitshapes_'+tagger_base+'.root'
        outfilePath_prepostfitshapes = os.path.join(outputDir, outfileName_prepostfitshapes)
        print 'PostFitShapesFromWorkspace'
        command = 'PostFitShapesFromWorkspace \\'
        command += '-w '+wsPath+' \\'
        command += '-o '+outfilePath_prepostfitshapes+' \\'
        command += '-f multidimfit_'+tagger_base+'.root:fit_mdf \\'
        command += '--postfit \\'
        command += '--sampling \\'
        command += '--print \\'
        command += '--total-shapes \\'
        command += '--covariance \\'
        combine_task('prepostfitshapes', command, logDir, outputDir)

        # write the found QCD normalization factor to disk
        print 'Writing scale factors to disk'
        # Check out this link to learn how to extract post-fit nuisance parameter values from RooFitResult:
        # https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/102x/src/MultiDimFit.cc#L434
        mdf_file = root.TFile.Open(os.path.join(outputDir, 'multidimfit_'+tagger_base+'.root'), 'READ')
        fit_result = mdf_file.Get('fit_mdf') # RooFitResult

        # r_nonQCD = fit_result.floatParsFinal().find('r_nonQCD')
        # r_nonQCD_central = r_nonQCD.getVal()
        # r_nonQCD_up = +(r_nonQCD.getMax('err68') - r_nonQCD_central if r_nonQCD.hasRange('err68') else r_nonQCD.getAsymErrorHi())
        # r_nonQCD_down = -(r_nonQCD.getMin('err68') - r_nonQCD_central if r_nonQCD.hasRange('err68') else r_nonQCD.getAsymErrorLo())
        # string_nonQCD = 'nonQCD norm factor: {0}  + {1} / - {2}'.format(r_nonQCD_central, r_nonQCD_up, r_nonQCD_down)
        # print string_nonQCD

        r_QCD_DD = fit_result.floatParsFinal().find('r_QCD_DD')
        r_QCD_DD_central = r_QCD_DD.getVal()
        r_QCD_DD_up = +(r_QCD_DD.getMax('err68') - r_QCD_DD_central if r_QCD_DD.hasRange('err68') else r_QCD_DD.getAsymErrorHi())
        r_QCD_DD_down = -(r_QCD_DD.getMin('err68') - r_QCD_DD_central if r_QCD_DD.hasRange('err68') else r_QCD_DD.getAsymErrorLo())
        string_QCD_DD = 'QCD_DD norm factor: {0}  + {1} / - {2}'.format(r_QCD_DD_central, r_QCD_DD_up, r_QCD_DD_down)
        print string_QCD_DD

        outfileName_NormFactors = 'final_norm_factors__'+tagger_base+'.root'
        outfilePath_NormFactors = os.path.join(outputDir, outfileName_NormFactors)
        with open(outfilePath_NormFactors.replace('.root', '.txt'), 'w') as outfile_NormFactors_txt:
            # outfile_NormFactors_txt.write(string_nonQCD+'\n')
            outfile_NormFactors_txt.write(string_QCD_DD+'\n')
        outfile_NormFactors = root.TFile.Open(outfilePath_NormFactors, 'RECREATE')
        outfile_NormFactors.cd()

        r_vec = root.TVectorF(1) # float vector with one element
        # r_vec[0] = r_nonQCD_central
        # outfile_NormFactors.WriteObject(r_vec, 'r_nonQCD_central')
        # r_vec[0] = r_nonQCD_up
        # outfile_NormFactors.WriteObject(r_vec, 'r_nonQCD_up')
        # r_vec[0] = r_nonQCD_down
        # outfile_NormFactors.WriteObject(r_vec, 'r_nonQCD_down')
        r_vec[0] = r_QCD_DD_central
        outfile_NormFactors.WriteObject(r_vec, 'r_QCD_DD_central')
        r_vec[0] = r_QCD_DD_up
        outfile_NormFactors.WriteObject(r_vec, 'r_QCD_DD_up')
        r_vec[0] = r_QCD_DD_down
        outfile_NormFactors.WriteObject(r_vec, 'r_QCD_DD_down')

        outfile_NormFactors.Close()


        print 'Writing correctly binned prefit and postfit shapes'
        file_PrePostFitShapes = root.TFile.Open(outfilePath_prepostfitshapes, 'READ')
        outfilePath_prepostfitshapes_rebinned = outfilePath_prepostfitshapes.replace('.root', '_rebinned.root')
        file_PrePostFitShapes_rebinned = root.TFile.Open(outfilePath_prepostfitshapes_rebinned, 'RECREATE')
        for p in ['prefit', 'postfit']:
            infile_folder = 'qcd_fit_'+p
            folder = file_PrePostFitShapes_rebinned.mkdir(p)
            folder.cd()
            data_obs_raw = file_PrePostFitShapes.Get(infile_folder+'/data_obs')
            data_obs = return_rebinned(data_obs_raw, binning)
            data_obs.Write('data_obs')
            full_stack_raw = file_PrePostFitShapes.Get(infile_folder+'/TotalProcs')
            full_stack = return_rebinned(full_stack_raw, binning)
            full_stack.Write('TotalProcs')
            stack_nonQCD_raw = file_PrePostFitShapes.Get(infile_folder+'/nonQCD')
            stack_nonQCD = return_rebinned(stack_nonQCD_raw, binning)
            stack_nonQCD.Write('nonQCD')
            stack_QCD_DD_raw = file_PrePostFitShapes.Get(infile_folder+'/QCD_DD')
            stack_QCD_DD = return_rebinned(stack_QCD_DD_raw, binning)
            stack_QCD_DD.Write('QCD_DD')
        file_PrePostFitShapes_rebinned.Close()
        file_PrePostFitShapes.Close()
