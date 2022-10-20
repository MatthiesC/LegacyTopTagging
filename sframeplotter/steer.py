#!/usr/bin/env python

import sys
import argparse
import os
import subprocess


#------------------#
# GLOBAL VARIABLES #
#------------------#

sys.path.append(os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/Analysis'))
from constants import _YEARS

# lumis = {
#     '2016': '35.92',
#     '2017': '41.53',
#     '2018': '59.74',
#     'run2': '137.19',
# }

# years = ['2016', '2017', '2018']
all_years = ['UL16preVFP', 'UL16postVFP', 'UL17', 'UL18']
all_channels = ['ele', 'muo']
# all_channels = ['muo']
combos = list()

#--------#
# PARSER #
#--------#

if not sys.argv[1:]: sys.exit('No arguments provided. Exit.')
parser = argparse.ArgumentParser()
parser.add_argument('--all', action='store_true')
parser.add_argument('--runii', action='store_true')
parser.add_argument('-c', '--channels', choices=all_channels, nargs='*', default=[])
parser.add_argument('-y', '--years', choices=all_years, nargs='*', default=[])
parser.add_argument('-s', '--singleeps', action='store_true')
parser.add_argument('-l', '--legend', action='store_true')
parser.add_argument('-r', '--shapenorm', action='store_true')
parser.add_argument('-n', '--no', choices=['data', 'qcd'], nargs='*', default=[])
parser.add_argument('-t', '--tar', action='store_true')
parser.add_argument('-p', '--plot', action='store_true')
parser.add_argument('-q', '--qcd', action='store_true') # QCD sideband region
args = parser.parse_args(sys.argv[1:])

if args.all == True:
    if args.runii == True:
        sys.exit('Do not use "--runii" option jointly with "--all" option. Abort.')
    if len(args.channels) or len(args.years):
        sys.exit('Do not use "--all" option jointly with --channels and/or --years options. Abort.')
    for year in all_years:
        for channel in all_channels:
            combos.append((year, channel))

elif args.runii == True:
    if not len(args.channels):
        combos.append(('run2', 'both'))
    else:
        for channel in args.channels:
            combos.append(('run2', channel))

else:
    if not len(args.channels) and not len(args.years):
        sys.exit('No years or channels given. Abort.')
    for year in (args.years if len(args.years) else all_years):
        for channel in (args.channels if len(args.channels) else all_channels):
            combos.append((year, channel))

if not len(combos):
    sys.exit('Nothing to plot. Exit.')

if args.tar and args.plot:
    sys.exit('Cannot plot and tar at once.')

print 'Working on:'
print combos
# sys.exit()

#---------#
# PROGRAM #
#---------#

print 'Using SFramePlotter at:', subprocess.check_output('which Plots', shell=True).strip('\n')

uhh2Dir = os.environ.get('CMSSW_BASE')+'/src/UHH2/'
mainselDir = uhh2Dir+'LegacyTopTagging/output/TagAndProbe/mainsel/'
templateDir = uhh2Dir+'LegacyTopTagging/sframeplotter/'
workDir = templateDir+'workdir/'
if not os.path.exists(workDir): os.mkdir(workDir)
sframeplotterBase = os.environ.get('CMSSW_BASE')+'/../SFramePlotter/'

FNULL = open(os.devnull, 'w')

for year, channel in combos:
    template_name = 'template.steer'
    if len(args.no):
        if 'data' in args.no: template_name = template_name.replace('.', '_woData.')
        if 'qcd' in args.no: template_name = template_name.replace('.', '_woQCD.')
    print 'Using template:', template_name
    template_file = open(templateDir+template_name, 'r')
    steerFilePath = workDir+'_'.join(['mainsel', year, channel])+'.steer'
    fCycleName = mainselDir+year+'/'+channel+'/nominal/hadded/uhh2.AnalysisModuleRunner'
    if args.qcd: fCycleName = fCycleName.replace('/nominal/', '/QCDsideband/')
    outputDirName = ('plots_single' if args.singleeps else 'plots')+'_'+'_'.join([year, channel])
    if args.qcd: outputDirName += '_QCDsideband'
    thisMainselDir = mainselDir+year+'/'+channel+'/'
    outputDir = thisMainselDir+outputDirName
    fOutputPsFile = outputDir+'/'+'_'.join(['mainsel', year, channel])+'.ps'
    if not os.path.exists(outputDir):
        print 'Create new directory:', outputDir
        os.mkdir(outputDir)
    with open(steerFilePath, 'w') as steerFile:
        print 'Create new steer file:', steerFilePath
        for line in template_file:
            newline = line
            newline = newline.replace('<<<fCycleName>>>', fCycleName)
            newline = newline.replace('<<<fOutputPsFile>>>', fOutputPsFile)
            newline = newline.replace('<<<fLumi>>>', str(_YEARS[year]['lumi_fb']))
            newline = newline.replace('<<<bSingleEPS>>>', 'true' if args.singleeps else 'false')
            newline = newline.replace('<<<bDrawLegend>>>', 'true' if args.legend else 'false')
            newline = newline.replace('<<<bShapeNorm>>>', 'true' if args.shapenorm else 'false')
            steerFile.write(newline)
    template_file.close()
    # Plots does not except absolute file paths to steer file, thus get relative path of steer file with working directory = sframeplotterBase
    steerFilePathRelative = steerFilePath.replace(uhh2Dir, '../'+os.environ.get('CMSSW_BASE').split('/')[-1]+'/src/UHH2/')
    command_Plots = 'Plots -f '+steerFilePathRelative
    command_targz = 'tar czf '+outputDirName+'.tar.gz '+outputDirName
    if not (args.plot or args.tar):
        print 'Command for SFramePlotter:', command_Plots
        print 'Command for tar/gzip:', command_targz
        print 'If you wish to run the plotter (compression), use option "-p" ("-t").'
    if args.plot:
        subprocess.Popen((command_Plots), shell=True, cwd=sframeplotterBase, stdout=FNULL, stderr=FNULL)
        # subprocess.Popen((command_Plots), shell=True, stdout=FNULL, stderr=FNULL)
    if args.tar:
        subprocess.Popen((command_targz), shell=True, cwd=thisMainselDir, stdout=FNULL, stderr=FNULL)
