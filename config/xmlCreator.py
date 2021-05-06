#!/usr/bin/env python3

import os
import sys
import csv
from collections import OrderedDict
import argparse
from itertools import permutations


class configContainer:

   '''Container class for configurating XML files'''

   userName = str()
   uhh2Dir = str()
   userMail = str()
   yearVars = dict()
   used_samples = OrderedDict()


   def __init__(self):

      self.userName = os.environ.get('USER')
      self.uhh2Dir = os.environ.get('CMSSW_BASE')+'/src/UHH2/'
      self.userMail = '@'.join([self.userName, 'mail.desy.de']) # Avoid spam due to public code on GitHub

      self.yearVars['yearVersions'] = { # e.g. 'v3' in case  of 2016v3
         'UL17': '',
         'UL18': '',
      }

      # Set these values such that there are no more than 2,500 jobs per preselection. This way, you can submit two preselections in parallel to avoid going over 5,000 jobs (current user limit for NAF)
      self.yearVars['preselFileSplit'] = {
         'UL17': '50',
         'UL18': '50',
      }

      self.yearVars['targetLumis'] = {
         'UL17': 41480.,
         'UL18': 59830.,
      }

      self.yearVars['lumiFiles'] = {
         'UL17': self.uhh2Dir+'common/data/UL17/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON_normtag.root',
         'UL18': self.uhh2Dir+'common/data/UL18/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON_normtag.root',
      }

      self.yearVars['pileupFiles'] = {
         'mc': {
            'UL17': self.uhh2Dir+'common/data/UL17/MyMCPileupHistogram_UL17.root',
            'UL18': self.uhh2Dir+'common/data/UL18/MyMCPileupHistogram_UL18.root',
         },
         'data': {
            'UL17': self.uhh2Dir+'common/data/UL17/MyDataPileupHistogram_UL17.root',
            'UL18': self.uhh2Dir+'common/data/UL18/MyDataPileupHistogram_UL18.root',
         },
         'dataUp': {
            'UL17': self.uhh2Dir+'common/data/UL17/MyDataPileupHistogram_UL17_72383.root',
            'UL18': self.uhh2Dir+'common/data/UL18/MyDataPileupHistogram_UL18_72383.root',
         },
         'dataDown': {
            'UL17': self.uhh2Dir+'common/data/UL17/MyDataPileupHistogram_UL17_66017.root',
            'UL18': self.uhh2Dir+'common/data/UL18/MyDataPileupHistogram_UL18_66017.root',
         },
      }

      self.yearVars['triggerSFFiles'] = {
         'UL17': self.uhh2Dir+'common/data/UL17/Efficiencies_muon_generalTracks_Z_Run2017_UL_SingleMuonTriggers.root',
         'UL18': self.uhh2Dir+'common/data/UL18/Efficiencies_muon_generalTracks_Z_Run2018_UL_SingleMuonTriggers.root',
      }

      self.yearVars['muonidSFFiles'] = {
         'UL17': self.uhh2Dir+'common/data/UL17/Efficiencies_muon_generalTracks_Z_Run2017_UL_ID.root',
         'UL18': self.uhh2Dir+'common/data/UL18/Efficiencies_muon_generalTracks_Z_Run2018_UL_ID.root',
      }

      self.yearVars['deepjetMCEffFiles'] = {
         'UL17': self.uhh2Dir+'LegacyTopTagging/data/',
         'UL18': self.uhh2Dir+'LegacyTopTagging/data/',
      }

      self.yearVars['deepjetSFFiles'] = {
         'UL17': self.uhh2Dir+'common/data/UL17/DeepJet_106XUL17SF_WPonly.csv',
         'UL18': self.uhh2Dir+'common/data/UL18/DeepJet_106XUL18SF_WPonly.csv',
      }

      self.systematics = list()


   @staticmethod
   def read_database(years: list):

      used_samples = OrderedDict()
      for year in years:
         used_samples[year] = list()
         with open('database.csv', 'r') as file:
            reader = csv.DictReader(file)
            for row in reader:
               use_me = row['use_for_sf']=='True' and row['year']==year
               if use_me:
                  used_samples[year].append(sampleEntity(row))
      configContainer.used_samples = used_samples


   def setup_systematics(self, selection: str, year: str):

      self.systematics.append(systEntity('jec', 'jecsmear_direction'))
      self.systematics.append(systEntity('jer', 'jersmear_direction'))
      if selection=='mainsel':
         self.systematics.append(systEntity('mur', 'ScaleVariationMuR'))
         self.systematics.append(systEntity('muf', 'ScaleVariationMuF'))
         self.systematics.append(systEntity('murmuf', 'N/A')) # Need to access ScaleVariationMuR and ScaleVariationMuF in another way
         self.systematics.append(systEntity('pileup', 'SystDirection_PS'), directions=['FSRup_2', 'FSRdown_2', 'ISRup_2', 'ISRdown_2'])
         self.systematics.append(systEntity('pileup', 'SystDirection_Pileup'))
         # if year in ['2016', '2017']:
         #    self.systematics.append(systEntity('prefiring', 'SystDirection_Prefiring'))
         self.systematics.append(systEntity('muontrigger', 'SystDirection_MuonTrigger'))
         self.systematics.append(systEntity('muonid', 'SystDirection_MuonID'))
         # self.systematics.append(systEntity('muoniso', 'SystDirection_MuonIso'))
         self.systematics.append(systEntity('btagging', 'SystDirection_BTag'))


class sampleEntity:

   '''Container to hold information about a data or MC sample, as read from CSV database'''

   def __init__(self, csvRow: OrderedDict):

      self.year = csvRow['year']
      self.is_data = True if csvRow['is_data']=='True' else False
      self.nickName = csvRow['nick_name']
      self.n_das = int(csvRow['n_das'])
      self.n_pnfs = int(csvRow['n_pnfs'])
      self.pnfs_sum_of_weights = float(self.n_pnfs) if csvRow['pnfs_sum_of_weights']=='-' else float(csvRow['pnfs_sum_of_weights'])
      self.xs_gen = None if csvRow['xs_gen']=='-' else float(csvRow['xs_gen'])
      self.xs_theo = None if csvRow['xs_theo']=='-' else float(csvRow['xs_theo'])
      if not self.is_data and self.xs_gen==None and self.xs_theo==None:
         sys.exit('No cross section given for sample labelled as simulation: '+self.nickName+' ('+self.year+')')
      self.xsection = self.xs_gen if self.xs_theo==None else self.xs_theo
      self.xmlPath = csvRow['xml_path']
      self.lumi = 1. if self.is_data else self.pnfs_sum_of_weights/self.xsection
      self.mainsel_versions = list()

      # init mainsel_versions
      merge_scenarios = ['FullyMerged', 'WMerged', 'QBMerged', 'NotMerged', 'Background']
      if self.nickName.startswith('ST_tW') or self.nickName.startswith('ST_tChannel') or self.nickName.startswith('ST_sChannel_had'): # in case we ever use the hadronicDecays s-channel sample
         for m in merge_scenarios:
            self.mainsel_versions.append('_'.join([self.nickName, m]))
      elif self.nickName.startswith('TTbarToSemiLeptonic'):
         for m in merge_scenarios:
            if m=='Background': continue
            self.mainsel_versions.append('_'.join([self.nickName, m]))
      # else:
      self.mainsel_versions.append(self.nickName)


class systEntity:

   '''Container to hold information about a systematic uncertainty'''

   def __init__(self, shortName: str, ctxName: str, defaultValue='nominal', directions=['up', 'down']):

      self.shortName = shortName
      self.ctxName = ctxName
      self.defaultValue = defaultValue
      self.directions = directions


class xmlCreator:

   '''Creates XML files for SFrame'''

   def __init__(self, selection: str, year: str):

      confCon = configContainer()
      self.uhh2Dir = confCon.uhh2Dir
      self.userName = confCon.userName
      self.userMail = confCon.userMail
      self.yearVars = confCon.yearVars
      self.sample_list = confCon.used_samples[year]
      confCon.setup_systematics(selection, year)
      self.systematics = confCon.systematics

      if selection not in ['presel', 'mainsel']:
         sys.exit('Given value of argument "selection" not valid. Abort.')
      self.selection = selection
      self.is_mainsel = True if selection=='mainsel' else False

      if year not in ['UL17', 'UL18']:
         sys.exit('Given value of argument "year" not valid. Abort.')
      self.year = year
      self.yearVersion = self.yearVars['yearVersions'][year]

      self.outputDirBase = self.uhh2Dir+'LegacyTopTagging/output/'
      if not os.path.isdir(self.outputDirBase):
         sys.exit('Warning: Make sure to create output directory via "ln -s". Abort.')
      self.outputDirBase += 'TagAndProbe/'

      self.xmlFileName = '_'.join(['parsedConfigFile', self.selection, self.year])+'.xml'
      self.xmlFilePathBase = self.uhh2Dir+'LegacyTopTagging/config/'+'_'.join(['config', self.selection, self.year])+'/'
      os.makedirs(self.xmlFilePathBase, exist_ok=True)
      self.xmlFilePath = self.xmlFilePathBase+self.xmlFileName

      self.write_xml_successful = False
      self.systXmlFilePaths = list()


   def write_xml(self):

      with open(self.xmlFilePath, 'w') as file:
         file.write('''<?xml version="1.0" encoding="UTF-8"?>\n''')
         file.write('''\n''')
         file.write('''<!DOCTYPE JobConfiguration PUBLIC "" "JobConfig.dtd"[\n''')
         file.write('''\n''')
         file.write('''<!ENTITY TargetLumi "'''+str(self.yearVars['targetLumis'][self.year])+'''">\n''')
         if self.is_mainsel:
            file.write('''<!ENTITY PRESELdir "'''+(self.outputDirBase+'presel/'+self.year+'/nominal/')+'''">\n''')
            file.write('''<!ENTITY PRESELfilename "uhh2.AnalysisModuleRunner">\n''')
         file.write('''<!ENTITY OUTPUTdir "'''+(self.outputDirBase+self.selection+'/'+self.year+'/nominal/')+'''">\n''')
         file.write('''<!ENTITY b_Cacheable "False">\n''')
         file.write('''<!ENTITY NEVT "-1">\n''')
         file.write('''<!ENTITY YEARsuffix "_'''+self.year+self.yearVersion+'''">\n''')
         file.write('''<!ENTITY PROOFdir "/nfs/dust/cms/user/'''+self.userName+'''/.proof2">\n''')
         file.write('''\n''')
         for s in self.sample_list:
            if self.is_mainsel:
               file.write('''<!ENTITY '''+s.nickName+''' "&PRESELdir;/&PRESELfilename;'''+('.DATA.' if s.is_data else '.MC.')+s.nickName+'''&YEARsuffix;.root">\n''')
            else:
               file.write('''<!ENTITY '''+s.nickName+''' SYSTEM "'''+self.uhh2Dir+'LegacyTopTagging/config/datasets/'+s.xmlPath+'''">\n''')
         file.write('''\n''')
         file.write(''']>\n''')
         file.write('''\n''')
         file.write('''<!--\n''')
         file.write('''<ConfigParse NEventsBreak="'''+('200000' if self.is_mainsel else '0')+'''" FileSplit="'''+('0' if self.is_mainsel else self.yearVars['preselFileSplit'][self.year])+'''" AutoResubmit="5"/>\n''')
         file.write('''<ConfigSGE RAM="4" DISK="3" Mail="'''+self.userMail+'''" Notification="as" Workdir="'''+'''_'''.join(['workdir', self.selection, self.year])+'''"/>\n''')
         file.write('''-->\n''')
         file.write('''\n''')
         file.write('''<!-- OutputLevel controls which messages are printed; set to VERBOSE or DEBUG for more verbosity, to WARNING or ERROR for less -->\n''')
         file.write('''<JobConfiguration JobName="ExampleCycleJob" OutputLevel="INFO">\n''')
         file.write('''<Library Name="libSUHH2LegacyTopTagging"/>\n''')
         file.write('''<Package Name="SUHH2LegacyTopTagging.par"/>\n''')
         file.write('''<Cycle Name="uhh2::AnalysisModuleRunner" OutputDirectory="&OUTPUTdir;/" PostFix="" TargetLumi="&TargetLumi;">\n''')
         file.write('''\n''')
         for s in self.sample_list:
            if self.is_mainsel:
               for v in s.mainsel_versions:
                  file.write('''<InputData Lumi="'''+str(s.lumi)+'''" NEventsMax="&NEVT;" Type="'''+('DATA' if s.is_data else 'MC')+'''" Version="'''+v+'''&YEARsuffix;" Cacheable="&b_Cacheable;"> <In FileName="&'''+s.nickName+''';" Lumi="0.0"/> <InputTree Name="AnalysisTree"/> <OutputTree Name="AnalysisTree"/> </InputData>\n''')
            else:
               file.write('''<InputData Lumi="'''+str(s.lumi)+'''" NEventsMax="&NEVT;" Type="'''+('DATA' if s.is_data else 'MC')+'''" Version="'''+s.nickName+'''&YEARsuffix;" Cacheable="&b_Cacheable;"> &'''+s.nickName+'''; <InputTree Name="AnalysisTree"/> <OutputTree Name="AnalysisTree"/> </InputData>\n''')
         file.write('''\n''')
         file.write('''<UserConfig>\n''')
         file.write('''\n''')
         file.write('''<Item Name="PrimaryVertexCollection" Value="offlineSlimmedPrimaryVertices"/>\n''')
         file.write('''<Item Name="METName" Value="slimmedMETs"/>\n''')
         file.write('''<Item Name="ElectronCollection" Value="slimmedElectronsUSER"/>\n''')
         file.write('''<Item Name="MuonCollection" Value="slimmedMuonsUSER"/>\n''')
         file.write('''<Item Name="JetCollection" Value="jetsAk4CHS"/>\n''')
         file.write('''<Item Name="GenJetCollection" Value="slimmedGenJets"/>\n''')
         file.write('''<Item Name="TopJetCollection" Value="hotvrPuppi"/>\n''')
         file.write('''<Item Name="GenTopJetCollection" Value="hotvrGen"/>\n''')
         file.write('''<Item Name="GenParticleCollection" Value="GenParticles"/>\n''')
         file.write('''<Item Name="GenInfoName" Value="genInfo"/>\n''')
         file.write('''<Item Name="additionalBranches" Value="jetsAk8PuppiSubstructure_SoftDropPuppi genjetsAk8SubstructureSoftDrop"/>\n''')
         if self.is_mainsel:
            file.write('''<Item Name="AK8Collection_rec" Value="jetsAk8PuppiSubstructure_SoftDropPuppi"/>\n''')
            file.write('''<Item Name="AK8Collection_gen" Value="genjetsAk8SubstructureSoftDrop"/>\n''')
         file.write('''\n''')
         file.write('''<Item Name="lumi_file" Value="'''+self.yearVars['lumiFiles'][self.year]+'''"/>\n''')
         file.write('''<Item Name="lumihists_lumi_per_bin" Value="500."/>\n''')
         file.write('''\n''')
         file.write('''<Item Name="MuonIDScaleFactorFile" Value="'''+self.yearVars['muonidSFFiles'][self.year]+'''"/>\n''')
         file.write('''<Item Name="TriggerScaleFactorFile" Value="'''+self.yearVars['triggerSFFiles'][self.year]+'''"/>\n''')
         file.write('''\n''')
         file.write('''<Item Name="pileup_directory" Value="'''+self.yearVars['pileupFiles']['mc'][self.year]+'''"/>\n''')
         file.write('''<Item Name="pileup_directory_data" Value="'''+self.yearVars['pileupFiles']['data'][self.year]+'''"/>\n''')
         if self.is_mainsel:
            file.write('''<Item Name="pileup_directory_data_up" Value="'''+self.yearVars['pileupFiles']['dataUp'][self.year]+'''"/>\n''')
            file.write('''<Item Name="pileup_directory_data_down" Value="'''+self.yearVars['pileupFiles']['dataDown'][self.year]+'''"/>\n''')
            file.write('''\n''')
            # file.write('''<Item Name="BTagMCEffFile" Value="'''+self.yearVars['deepjetMCEffFiles'][self.year]+'''"/>\n''')
            file.write('''<Item Name="BTagScaleFactorFile" Value="'''+self.yearVars['deepjetSFFiles'][self.year]+'''"/>\n''')
         file.write('''\n''')
         file.write('''<!-- Keys for systematic uncertainties -->\n''')
         for syst in self.systematics:
            if syst.shortName == 'murmuf': continue
            file.write('''<Item Name="'''+syst.ctxName+'''" Value="'''+syst.defaultValue+'''"/>\n''')
         file.write('''\n''')
         file.write('''<!-- Tell AnalysisModuleRunner NOT to use the MC event weight from SFrame; rather let MCLumiWeight (called via CommonModules) calculate the MC event weight. The MC event weight assigned by MCLumiWeight is InputData.Lumi / Cycle.TargetLumi. -->\n''')
         file.write('''<Item Name="use_sframe_weight" Value="false"/>\n''')
         file.write('''<Item Name="AnalysisModule" Value="'''+('TagAndProbeMainSelectionModule' if self.is_mainsel else 'TagAndProbePreSelectionModule')+'''"/>\n''')
         file.write('''<Item Name="uhh2Dir" Value="'''+self.uhh2Dir+'''"/>\n''')
         file.write('''\n''')
         file.write('''<!-- Switch for debugging of the central AnalysisModule -->\n''')
         file.write('''<Item Name="debug" Value="false"/>\n''')
         file.write('''\n''')
         file.write('''</UserConfig>\n''')
         file.write('''\n''')
         file.write('''</Cycle>\n''')
         file.write('''</JobConfiguration>\n''')

      self.write_xml_successful = True

      print('Created '+self.xmlFilePath)

      return self.xmlFilePath


   def write_systematics_xml(self, syst: systEntity):

      if not self.write_xml_successful:
         sys.exit('xmlCreator::write_xml() not called. Danger of parsing potentially outdated XML file. Exit.')

      for direction in syst.directions:
         systXmlFilePath = self.xmlFilePath.replace('.xml', '_')+'_'.join(['syst', syst.shortName, direction])+'.xml'
         infile = open(self.xmlFilePath, 'r')
         with open(systXmlFilePath, 'w') as outfile:
            for line in infile:
               newline = line
               # if newline.startswith('<!ENTITY PRESELdir') or newline.startswith('<!ENTITY OUTPUTdir'):
                  # newline = newline.replace('/nominal/', '/'+'_'.join(['syst', syst.shortName, direction])+'/')
               if newline.startswith('<!ENTITY PRESELdir') and (syst.shortName == 'jec' or syst.shortName == 'jer'):
                  newline = newline.replace('/nominal/', '/'+'_'.join(['syst', syst.shortName, direction])+'/')
               if newline.startswith('<!ENTITY OUTPUTdir'):
                  newline = newline.replace('/nominal/', '/'+'_'.join(['syst', syst.shortName, direction])+'/')
               if newline.startswith('<ConfigSGE'):
                  newline = newline.replace('"/>', '_'+'_'.join(['syst', syst.shortName, direction])+'"/>')
               if syst.shortName != 'murmuf' and newline.startswith('<Item Name="'+syst.ctxName):
                  newline = newline.replace(syst.defaultValue, direction)
               if syst.shortName == 'murmuf' and newline.startswith('<Item Name="ScaleVariationMu'):
                  newline = newline.replace(syst.defaultValue, direction)
               outfile.write(newline)
         infile.close()

         print('Created '+systXmlFilePath)

         self.systXmlFilePaths.append(systXmlFilePath)


   def write_all_systematics_xmls(self):

      for syst in self.systematics:

         self.write_systematics_xml(syst)

      return self.systXmlFilePaths


if __name__=='__main__':

   selections = ['presel', 'mainsel']
   years = ['UL17', 'UL18']

   if not sys.argv[1:]: sys.exit('No arguments provided. Exit.')
   parser = argparse.ArgumentParser()
   parser.add_argument('--all', action='store_true', help='Create XML files for all selections and years.')
   parser.add_argument('--syst', action='store_true', help='Create XML files for systematic uncertainties.')
   parser.add_argument('-s', '--selections', choices=selections, nargs='*', default=[])
   parser.add_argument('-y', '--years', choices=years, nargs='*', default=[])
   parser.add_argument('-a', '--auto-complete', action='store_true', help='Auto-complete arguments if not all arguments for selections and years are given.')
   args = parser.parse_args(sys.argv[1:])

   if(args.all == True):
      if(len(args.selections) + len(args.years) != 0):
         sys.exit('Not allowed to use "--all" option jointly with manually given selection or year argument. Exit.')
      args.selections = selections
      args.years = years
      if(args.auto_complete):
         print('Warning: You already specified "--all". Therefore, "--auto-complete" will not have any effect.')
   else:
      if args.auto_complete:
         print('Auto-completing arguments.')
         if not args.selections: args.selections = selections
         if not args.years: args.years = years
      else:
         for p in permutations([args.selections, args.years]):
            if p[0] and not p[1]:
               sys.exit('You specified arguments for at least one of the two options: "--selections" or "--years", but not for both of them. Also, you did not specify "--auto-complete" to compensate for this. Exit.')

   print('Going to create XML files for:')
   print('  Selections: '+', '.join(str(x) for x in args.selections))
   print('  Years: '+', '.join(str(x) for x in args.years))

   configContainer.read_database(args.years)

   for selection in args.selections:
      for year in args.years:
         x = xmlCreator(selection, year)
         x.write_xml()
         if args.syst:
            x.write_all_systematics_xmls()
