from __future__ import print_function

import sys

from HiggsAnalysis.CombinedLimit.PhysicsModel import PhysicsModel

# https://hypernews.cern.ch/HyperNews/CMS/get/higgs-combination/1653.html
# https://gitlab.cern.ch/gouskos/boostedjetcalibration/-/blob/master/TagAndProbeExtended.py

class LegacyTopTaggingModel(PhysicsModel):

    def __init__(self):
        super(LegacyTopTaggingModel, self).__init__()
        self.merge_scenarios = {}
        self.sf_naming_scheme = None
        self.sf_range = None
        self.pt_bin_names = None
        self.years = None
        self.regions = ['Pass', 'Fail']

    def setPhysicsOptions(self, physOptions):
        for po in physOptions:
            if po.startswith('merge_scenarios='):
                # self.merge_scenarios = po.replace('merge_scenarios=', '').split(',')
                split = po.replace('merge_scenarios=', '').split(',')
                for s in split:
                    merge_scenario = s.split(':')[0]
                    use_TagAndProbe = s.split(':')[1] == 'TagAndProbe'
                    self.merge_scenarios[merge_scenario] = use_TagAndProbe
            if po.startswith('sf_naming_scheme='):
                self.sf_naming_scheme = po.replace('sf_naming_scheme=', '')
            if po.startswith('sf_range='):
                self.sf_range = po.replace('sf_range=', '')
            if po.startswith('pt_bin_names='):
                self.pt_bin_names = po.replace('pt_bin_names=', '').split(',')
            if po.startswith('years='):
                self.years = po.replace('years=', '').split(',')

    def getMergeScenarioFromProcess(self, process):
        split = process.split('__')
        for s in split:
            for msc in self.merge_scenarios.keys():
                msc_string = 'MSc_'+msc
                if s == msc_string:
                    return msc
        return 'dummy'
        sys.exit('Merge scenario not identified')

    def getYearFromProcess(self, process):
        split = process.split('__')
        for s in split:
            if s in self.years:
                return s
        return 'dummy'
        sys.exit('Year not identified')

    def getPtBinFromProcess(self, process):
        split = process.split('__')
        for s in split:
            if s in self.pt_bin_names:
                return s
        return 'dummy'
        sys.exit('PtBin not identified')

    def getRegionFromBin(self, bin):
        split = bin.split('-')
        for s in split:
            if s in self.regions:
                return s
        return 'dummy'
        sys.exit('Region not identified')

    def getYearFromBin(self, bin):
        split = bin.split('-')
        for s in split:
            if s in self.years:
                return s
        return 'dummy'
        sys.exit('Year not identified')

    def getPtBinFromBin(self, bin):
        split = bin.split('-')
        for s in split:
            if s in self.pt_bin_names:
                return s
        return 'dummy'
        sys.exit('PtBin not identified')

    def doParametersOfInterest(self):
        '''Create POI and other parameters, and define the POI set.'''
        pois = []

        expected_yields = {}
        for bin in self.DC.bins:
            for process in self.DC.exp.get(bin).keys():
                expected_yields.setdefault(bin, {})[process] = self.DC.exp.get(bin).get(process)

        # get sum of pass/fail expected yields for TTbar/ST for each combine_channel/msc combo
        expected_yields_top = {}
        for msc in self.merge_scenarios.keys():
            for year in self.years:
                for pt_bin_name in self.pt_bin_names:
                    for region in self.regions:
                        for bin in self.DC.bins:
                            if not year == self.getYearFromBin(bin): continue
                            if not pt_bin_name == self.getPtBinFromBin(bin): continue
                            if not region == self.getRegionFromBin(bin): continue
                            for process in self.DC.exp.get(bin).keys():
                                if not msc == self.getMergeScenarioFromProcess(process): continue
                                expected_yields_top.setdefault(msc, {}).setdefault(year, {}).setdefault(pt_bin_name, {}).setdefault(region, 0)
                                expected_yields_top[msc][year][pt_bin_name][region] += self.DC.exp.get(bin).get(process)

        print(expected_yields_top)

        self.sf_naming_scheme = '__'.join(['SF', r'{msc}', r'{year}', r'{pt_bin_name}'])

        for msc, use_TagAndProbe in self.merge_scenarios.items():
            for year in self.years:
                for pt_bin_name in self.pt_bin_names:
                    # Tagging scale factors (pass region)
                    sf_name = self.sf_naming_scheme.format(msc=msc, year=year, pt_bin_name=pt_bin_name)
                    sf_range = self.sf_range
                    self.modelBuilder.doVar(sf_name+sf_range)
                    pois.append(sf_name)
                    if use_TagAndProbe: # do not create the antisf variable if not needed
                        # Anti-tagging scale factors (fail region)
                        antisf_name = sf_name.replace('SF_', 'AntiSF_')
                        self.modelBuilder.factory_(
                            # 'expr::{antisf_name}("max(0.,1.+(1.-@0)*{N_pass}/{N_fail})", {sf_name})'.format(antisf_name=antisf_name, N_pass=exp_pass.get(task).get(cat), N_fail=exp_fail.get(task).get(cat), sf_name=sf_name)
                            'expr::{antisf_name}("max(0.,1.+(1.-@0)*{N_pass}/{N_fail})", {sf_name})'.format(antisf_name=antisf_name, N_pass=expected_yields_top[msc][year][pt_bin_name]['Fail'], N_fail=expected_yields_top[msc][year][pt_bin_name]['Fail'], sf_name=sf_name)
                        ) # FIXME: This could probably crash in the unlikely case that N_fail = 0 (division by zero)

        #________________________________
        self.modelBuilder.doSet("POI", ','.join(pois))

    def getYieldScale(self, bin, process):
        '''Return the name of a RooAbsReal to scale this yield by, or the two special values 1 and 0 (don't scale, and set to zero).'''
        if self.DC.isSignal.get(process):
            msc = self.getMergeScenarioFromProcess(process)
            use_TagAndProbe = self.merge_scenarios.get(msc)
            year = self.getYearFromBin(bin)
            pt_bin_name = self.getPtBinFromBin(bin)
            region = self.getRegionFromBin(bin)
            sf_name = self.sf_naming_scheme.format(msc=msc, year=year, pt_bin_name=pt_bin_name)
            if use_TagAndProbe:
                if region == 'Pass':
                    return sf_name
                elif region == 'Fail':
                    return sf_name.replace('SF_', 'AntiSF_')
                else:
                    return 1
            else:
                return sf_name
        else:
            return 1

lttModel = LegacyTopTaggingModel()
