#!/usr/bin/env python3

import os
import sys
from collections import OrderedDict
import argparse
from termcolor import colored
import subprocess

# to include CrossSectionHelper:
# sys.path.append(os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/common/UHH2-datasets'))
# from CrossSectionHelper import MCSampleValuesHelper
# helper = MCSampleValuesHelper()

sys.path.append(os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/Analysis'))
from constants import _YEARS

from database import samplesDict
from xmlCreator import sampleEntity


class wpXmlCreator:

    def __init__(self, years):

        self.years = years
        self.workdirNames = OrderedDict()
        for year in self.years:
            self.workdirNames[year] = 'workdir_WorkingPointStudy_'+year
        self.configDir = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/config/config_WorkingPointStudy')
        os.makedirs(self.configDir, exist_ok=True)
        self.xmlFilePaths = OrderedDict()

        self.samples = OrderedDict()
        for year in self.years:
            for k, v in samplesDict.items():
                use_me = v.get('years') == None or year in v.get('years', [])
                use_me = use_me and ('wp' in v.get('analysis', []))
                if use_me:
                    sample_entity = sampleEntity('106X_v2', (k, v, year,))
                    if not os.path.isfile(sample_entity.xmlPath):
                        print(colored('XML for sample  '+sample_entity.nickName+' ('+year+')  does not exist. Skipping this sample', 'red'))
                        continue
                    self.samples.setdefault(year, []).append(sample_entity)

    def write_xml(self, year):

        xml_name = 'parsedConfigFile_WorkingPointStudy_'+year+'.xml'
        xml_path = os.path.join(self.configDir, xml_name)
        self.xmlFilePaths[year] = xml_path

        with open(xml_path, 'w') as file:
            file.write('''<?xml version="1.0" encoding="UTF-8"?>\n''')
            file.write('''\n''')
            file.write('''<!DOCTYPE JobConfiguration PUBLIC "" "JobConfig.dtd"[\n''')
            file.write('''\n''')
            file.write('''<!ENTITY TargetLumi "'''+str(_YEARS[year].get('lumi_pb'))+'''">\n''')
            file.write('''<!ENTITY OUTPUTdir "'''+os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/output/WorkingPointStudy', year)+'''">\n''')
            file.write('''<!ENTITY b_Cacheable "False">\n''')
            file.write('''<!ENTITY NEVT "-1">\n''')
            file.write('''<!ENTITY YEARsuffix "_'''+year+'''">\n''')
            file.write('''<!ENTITY PROOFdir "'''+os.path.join('/nfs/dust/cms/user', os.environ.get('USER'), '.proof2')+'''">\n''')
            file.write('''\n''')
            for s in self.samples[year]:
                file.write('''<!ENTITY '''+s.nickName+''' SYSTEM "'''+s.xmlPath+'''">\n''')
            file.write('''\n''')
            file.write(''']>\n''')
            file.write('''\n''')
            file.write('''<!--\n''')
            file.write('''<ConfigParse NEventsBreak="0" FileSplit="25" AutoResubmit="5"/>\n''')
            file.write('''<ConfigSGE RAM="4" DISK="3" Mail="'''+os.environ.get('USER')+'''@mail.desy.de" Notification="as" Workdir="'''+self.workdirNames[year]+'''"/>\n''')
            file.write('''-->\n''')
            file.write('''\n''')
            file.write('''<!-- OutputLevel controls which messages are printed; set to VERBOSE or DEBUG for more verbosity, to WARNING or ERROR for less -->\n''')
            file.write('''<JobConfiguration JobName="ExampleCycleJob" OutputLevel="INFO">\n''')
            file.write('''<Library Name="libSUHH2LegacyTopTagging"/>\n''')
            file.write('''<Package Name="SUHH2LegacyTopTagging.par"/>\n''')
            file.write('''<Cycle Name="uhh2::AnalysisModuleRunner" OutputDirectory="&OUTPUTdir;/" PostFix="" TargetLumi="&TargetLumi;">\n''')
            file.write('''\n''')
            for s in self.samples[year]:
                file.write('''<InputData Lumi="'''+str(s.lumi)+'''" NEventsMax="&NEVT;" Type="'''+('DATA' if s.is_data else 'MC')+'''" Version="'''+s.nickName+'''&YEARsuffix;" Cacheable="&b_Cacheable;"> &'''+s.nickName+'''; <InputTree Name="AnalysisTree"/> <OutputTree Name="AnalysisTree"/> </InputData>\n''')
            file.write('''\n''')
            file.write('''<UserConfig>\n''')
            file.write('''\n''')
            file.write('''<Item Name="TopJetCollection" Value="jetsAk8PuppiSubstructure_SoftDropPuppi"/>\n''')
            file.write('''<Item Name="GenTopJetCollection" Value="genjetsAk8SubstructureSoftDrop"/>\n''')
            file.write('''<Item Name="GenParticleCollection" Value="GenParticles"/>\n''')
            file.write('''<Item Name="GenInfoName" Value="genInfo"/>\n''')
            file.write('''\n''')
            file.write('''<Item Name="pileup_directory" Value="'''+os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/common/UHH2-data', year, 'MyMCPileupHistogram_'+year+'.root')+'''"/>\n''')
            file.write('''<Item Name="pileup_directory_data" Value="'''+os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/common/UHH2-data', year, 'MyDataPileupHistogram_'+year+'.root')+'''"/>\n''')
            file.write('''\n''')
            file.write('''<!-- Keys for systematic uncertainties -->\n''')
            file.write('''<Item Name="jecsmear_direction" Value="nominal"/>\n''')
            file.write('''<Item Name="jersmear_direction" Value="nominal"/>\n''')
            file.write('''\n''')
            file.write('''<!-- Tell AnalysisModuleRunner NOT to use the MC event weight from SFrame; rather let MCLumiWeight (called via CommonModules) calculate the MC event weight. The MC event weight assigned by MCLumiWeight is InputData.Lumi / Cycle.TargetLumi. -->\n''')
            file.write('''<Item Name="use_sframe_weight" Value="false"/>\n''')
            file.write('''<Item Name="AnalysisModule" Value="WorkingPointModule"/>\n''')
            file.write('''\n''')
            file.write('''<!-- Switch for debugging of the central AnalysisModule -->\n''')
            file.write('''<Item Name="debug" Value="false"/>\n''')
            file.write('''\n''')
            file.write('''</UserConfig>\n''')
            file.write('''\n''')
            file.write('''</Cycle>\n''')
            file.write('''</JobConfiguration>\n''')

        print('Created '+xml_path)
        return xml_path

    def write_bash_scripts(self):

        scriptFilePath_sframe_batch = os.path.join(self.configDir, 'run_all_sframe_batch.sh')
        with open(scriptFilePath_sframe_batch, 'w') as outfile:
            outfile.write('#!/bin/bash\n')
            newline_base = 'sframe_batch.py $1 '
            for year in self.years:
                outfile.write(newline_base+self.xmlFilePaths[year]+'\n')
        p = subprocess.Popen('chmod +x '+scriptFilePath_sframe_batch, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        p.wait()
        print('Created '+scriptFilePath_sframe_batch)

        scriptFilePath_run_local = os.path.join(self.configDir, 'run_all_local.sh')
        with open(scriptFilePath_run_local, 'w') as outfile:
            outfile.write('#!/bin/bash\n')
            newline_base = 'python run_local.py $1 '
            for year in self.years:
                outfile.write(newline_base+self.workdirNames[year]+'\n')
        p = subprocess.Popen('chmod +x '+scriptFilePath_run_local, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        p.wait()
        print('Created '+scriptFilePath_run_local)
        link_run_local_command = 'ln -s '+os.path.abspath(os.path.join(self.configDir, '..', 'run_local.py'))+' '+os.path.join(self.configDir, 'run_local.py')
        p = subprocess.Popen(link_run_local_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        p.wait()

    def run(self):

        for year in self.years:
            self.write_xml(year)
        self.write_bash_scripts()

if __name__=='__main__':

    years = ['UL16preVFP', 'UL16postVFP', 'UL17', 'UL18']

    if not sys.argv[1:]: sys.exit('No arguments provided. Exit.')
    parser = argparse.ArgumentParser()
    parser.add_argument('--all', action='store_true', help='Create XML files for all years.')
    parser.add_argument('-y', '--years', choices=years, nargs='*', default=[])
    args = parser.parse_args(sys.argv[1:])

    if(args.all == True):
        if len(args.years) != 0:
            sys.exit('Not allowed to use "--all" option jointly with manually given year argument. Exit.')
        args.years = years

    print('Going to create XML files for:')
    print('  Years: '+', '.join(str(x) for x in args.years))

    x = wpXmlCreator(args.years)
    x.run()
