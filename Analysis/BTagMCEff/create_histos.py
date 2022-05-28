#!/usr/bin/env python

import os
import sys
import argparse
import subprocess

#____________________________________________________________________________________________________
# We need the varial library for this script
# Clone it if it doesn't already exist
try:
    import varial
except ImportError:
    varial_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'Varial')

    if not os.path.isdir(varial_dir):
        print 'Cannot find Varial. Trying to clone it...'
        git_url = 'https://github.com/UHH2/Varial'
        cmd = 'git clone ' + git_url
        subprocess.check_output(cmd, shell=True)
        # Check again
        if not os.path.isdir(varial_dir):
            print 'Error cloning Varial from', git_url
            print 'Please try this manually'
            exit(-1)

    if os.path.isdir(varial_dir):
        sys.path.append(varial_dir)
        import varial
    else:
        print 'ImportError. You need Varial for this script:'
        print 'https://github.com/HeinerTholen/Varial'
        exit(-1)

#____________________________________________________________________________________________________
uhh2Dir = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2')
outputDir = os.path.join(uhh2Dir, 'LegacyTopTagging/output/TagAndProbe/mainsel')
histDir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'files')
print histDir

year_choices = ['UL16preVFP', 'UL16postVFP', 'UL17', 'UL18']
channel_choices = ['muo']

if not sys.argv[1:]: sys.exit('No arguments provided. Exit.')
parser = argparse.ArgumentParser()
parser.add_argument('-y', '--years', choices=year_choices, nargs='*', default=year_choices)
parser.add_argument('-c', '--channels', choices=channel_choices, nargs='*', default=channel_choices)
parser.add_argument('-f', '--folder', type=str, default='')
args = parser.parse_args(sys.argv[1:])

#____________________________________________________________________________________________________
folder = args.folder
btag_histo_names = (
    'BTagMCEffFlavBPassing',
    'BTagMCEffFlavBTotal',
    'BTagMCEffFlavCPassing',
    'BTagMCEffFlavCTotal',
    'BTagMCEffFlavUDSGPassing',
    'BTagMCEffFlavUDSGTotal',
)
if folder == '':
    filter_func = lambda w: w.name in btag_histo_names
else:
    filter_func = lambda w: w.name in btag_histo_names and folder in w.in_file_path

for year in args.years:
    for channel in args.channels:
        print 'Working on:', year
        input_pattern = os.path.join(outputDir, year, channel, 'nominal', 'uhh2.AnalysisModuleRunner.MC.*AllMergeScenarios*.root')
        print 'Using these MC files as input:', input_pattern
        # for file in varial.gen.resolve_file_pattern(input_pattern):
        #     print file

        # load and merge histograms
        pipe = varial.gen.open_filter_load(input_pattern, filter_func)
        pipe = varial.gen.sort_group_merge(pipe, lambda w: w.name)
        pipe = list(pipe)  # pull everything through the generators

        # do some checks
        if len(pipe) != 6:
            print 'Error. I need histograms for all of these names: ', btag_histo_names
            print 'I only found these: ', list(w.in_file_path for w in pipe)
            exit(-1)

        empty_hists = list(w.name for w in pipe if not w.obj.Integral())
        if empty_hists:
            print 'Error: integral of histograms "%s" is 0.0! Exit.' % empty_hists
            exit(-1)

        # eff histos
        pipe = varial.gen.gen_make_eff_graphs(
            pipe,
            postfix_sub='Passing',
            postfix_tot='Total',
            new_postfix='Eff',
            yield_everything=True,
            eff_func=varial.op.div,
        )

        # write out
        os.system('mkdir -p '+histDir)
        outfile = os.path.join(histDir, 'BTagMCEfficiencyHists'+folder+'_'+year+'_'+channel+'.root')
        print 'Saving histograms to file:', outfile
        f = varial.ROOT.TFile(outfile, 'RECREATE')
        for w in pipe:
            print '  writing histogram:', w.name
            w.obj.SetName(w.name)
            w.obj.Write()
        f.Close()
