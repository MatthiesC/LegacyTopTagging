from HiggsAnalysis.CombinedLimit.PhysicsModel import PhysicsModel

# https://hypernews.cern.ch/HyperNews/CMS/get/higgs-combination/1653.html
# https://gitlab.cern.ch/gouskos/boostedjetcalibration/-/blob/master/TagAndProbeExtended.py

class LegacyTopTaggingModel(PhysicsModel):

    def __init__(self):
        super(LegacyTopTaggingModel, self).__init__()
        self.ttbar_categories = None
        self.sf_name_prefix = 'SF_'
        self.antisf_name_prefix = 'AntiSF_'
        # self.postfix = None
        self.tasks = None

    def setPhysicsOptions(self, physOptions):
        for po in physOptions:
            if po.startswith('categories='):
                # physOptions.remove(po)
                self.ttbar_categories = po.replace('categories=', '').split(',')
            # if po.startswith('postfix='):
            #     self.postfix = po.replace('postfix=', '')
            if po.startswith('tasks='):
                self.tasks = po.replace('tasks=', '').split(',')
        if not self.ttbar_categories:
            self.ttbar_categories = [
                # 'TTbar',
                'FullyMerged',
                # 'SemiMerged',
                # 'WMerged',
                # 'QBMerged',
                # 'BkgOrNotMerged',
                # 'NotFullyOrWMerged',
                'YllufMerged',
            ]
            print 'No ttbar categories provided via "--PO categories=catA,catB,...". Assuming these:'
            print self.ttbar_categories

    def getProcessCategory(self, process):
        result = None
        if self.DC.isSignal.get(process):
            for cat in self.ttbar_categories:
                if cat in process.split('_'):
                # if cat in process:
                    result = cat
                    break
        return result

    def getTaskFromBin(self, bin):
        if self.tasks:
            result = None
            for task in self.tasks:
                if task in bin:
                    result = task
                    break
            return result
        else:
            return 'dummy'

    def getRegionFromBin(self, bin):
        split = bin.split('_')
        for s in split:
            if 'fail' in s.lower() or 'pass' in s.lower():
                return s
        return 'Region not identified'

    def doParametersOfInterest(self):
        '''Create POI and other parameters, and define the POI set.'''
        pois = []
        #________________________________
        # Get the expected ('prefit') yields of signal processes
        exp_pass = dict()
        exp_fail = dict()
        exp_passW = dict()
        exp_failW = dict()
        # expected_yields = dict()
        print self.DC.bins
        for bin in self.DC.bins:
            task = self.getTaskFromBin(bin)
            if self.getRegionFromBin(bin) == 'Pass':
                exp_pass[task] = dict()
            elif self.getRegionFromBin(bin) == 'Fail':
                exp_fail[task] = dict()
            elif self.getRegionFromBin(bin) == 'PassW':
                exp_passW[task] = dict()
            elif self.getRegionFromBin(bin) == 'FailW':
                exp_failW[task] = dict()
            # expected_yields[bin] = dict()
            for process in self.DC.exp.get(bin).keys():
                # expected_yields[bin][process] = self.DC.exp.get(bin).get(process)
                cat = self.getProcessCategory(process)
                if self.getRegionFromBin(bin) == 'Pass':
                    exp_pass[task][cat] = self.DC.exp.get(bin).get(process)
                elif self.getRegionFromBin(bin) == 'Fail':
                    exp_fail[task][cat] = self.DC.exp.get(bin).get(process)
                elif self.getRegionFromBin(bin) == 'PassW':
                    exp_passW[task][cat] = self.DC.exp.get(bin).get(process)
                elif self.getRegionFromBin(bin) == 'FailW':
                    exp_failW[task][cat] = self.DC.exp.get(bin).get(process)
        print 'exp_pass:', exp_pass
        print 'exp_fail:', exp_fail
        print 'exp_passW:', exp_passW
        print 'exp_failW:', exp_failW
        if self.tasks:
            for task in self.tasks:
                #________________________________
                # Overall ttbar rate
                # r_name = 'r_overall_'+task
                # r_range = '[1,0,5]'
                # self.modelBuilder.doVar(r_name+r_range)
                # pois.append(r_name)
                #________________________________
                for cat in self.ttbar_categories:
                    #________________________________
                    # Tagging scale factor (pass region)
                    sf_name = self.sf_name_prefix+cat+'_'+task
                    sf_range = '[1,0.2,2.0]'
                    self.modelBuilder.doVar(sf_name+sf_range)
                    pois.append(sf_name)
                    # r_pass_name = 'r_pass_%s_%s' % (task, cat)
                    # self.modelBuilder.factory_(
                    #     'expr::{r_pass_name}("@0*@1", {sf_name}, {r_name})'.format(r_pass_name=r_pass_name, sf_name=sf_name, r_name=r_name)
                    # )
                    # pois.append(r_pass_name)
                    #________________________________
                    # Anti-tagging scale factor (fail region)
                    antisf_name = self.antisf_name_prefix+cat+'_'+task
                    self.modelBuilder.factory_(
                        # 'expr::{antisf_name}("max(0.,1.+(1.-@0)*{N_pass}/{N_fail})", {sf_name})'.format(antisf_name=antisf_name, N_pass=exp_pass.get(cat), N_fail=exp_fail.get(cat), sf_name=sf_name)
                        'expr::{antisf_name}("max(0.,1.+(1.-@0)*{N_pass}/{N_fail})", {sf_name})'.format(antisf_name=antisf_name, N_pass=exp_pass.get(task).get(cat), N_fail=exp_fail.get(task).get(cat), sf_name=sf_name)
                        # 'expr::{antisf_name}("max(0.,1.+(1.-@0)*{N_pass}/{N_fail})", {sf_name})'.format(antisf_name=antisf_name, N_pass=exp_pass.get(task).get(cat), N_fail=exp_passW.get(task).get(cat)+exp_failW.get(task).get(cat), sf_name=sf_name)
                    ) # FIXME: This could probably crash in the unlikely case that N_fail = 0 (division by zero)
                    # pois.append(antisf_name)
                    # r_fail_name = 'r_fail_%s_%s' % (task, cat)
                    # self.modelBuilder.factory_(
                    #     'expr::{r_fail_name}("@0*@1", {antisf_name}, {r_name})'.format(r_fail_name=r_fail_name, antisf_name=antisf_name, r_name=r_name)
                    # )
                    # pois.append(r_fail_name)
        else:
            #________________________________
            # Overall ttbar rate
            # r_name = 'r_overall'
            # r_range = '[1,0,5]'
            # self.modelBuilder.doVar(r_name+r_range)
            # pois.append(r_name)
            #________________________________
            for cat in self.ttbar_categories:
                #________________________________
                # Tagging scale factor (pass region)
                sf_name = self.sf_name_prefix+cat
                sf_range = '[1,0.2,2.0]'
                self.modelBuilder.doVar(sf_name+sf_range)
                pois.append(sf_name)
                # r_pass_name = 'r_pass_%s' % cat
                # self.modelBuilder.factory_(
                #     'expr::{r_pass_name}("@0*@1", {sf_name}, {r_name})'.format(r_pass_name=r_pass_name, sf_name=sf_name, r_name=r_name)
                # )
                # pois.append(r_pass_name)
                #________________________________
                # Anti-tagging scale factor (fail region)
                antisf_name = self.antisf_name_prefix+cat
                self.modelBuilder.factory_(
                    # 'expr::{antisf_name}("max(0.,1.+(1.-@0)*{N_pass}/{N_fail})", {sf_name})'.format(antisf_name=antisf_name, N_pass=exp_pass.get('dummy').get(cat), N_fail=exp_fail.get('dummy').get(cat), sf_name=sf_name)
                    'expr::{antisf_name}("max(0.,1.+(1.-@0)*{N_pass}/{N_fail})", {sf_name})'.format(antisf_name=antisf_name, N_pass=exp_pass.get('dummy').get(cat), N_fail=exp_passW.get('dummy').get(cat)+exp_failW.get('dummy').get(cat), sf_name=sf_name)
                ) # FIXME: This could probably crash in the unlikely case that N_fail = 0 (division by zero)
                # pois.append(antisf_name)
                # r_fail_name = 'r_fail_%s' % cat
                # self.modelBuilder.factory_(
                #     'expr::{r_fail_name}("@0*@1", {antisf_name}, {r_name})'.format(r_fail_name=r_fail_name, antisf_name=antisf_name, r_name=r_name)
                # )
                # pois.append(r_fail_name)
        #________________________________
        self.modelBuilder.doSet("POI", ','.join(pois))

    def getYieldScale(self, bin, process):
        '''Return the name of a RooAbsReal to scale this yield by, or the two special values 1 and 0 (don't scale, and set to zero).'''
        task = self.getTaskFromBin(bin)
        if self.DC.isSignal.get(process):
            cat = self.getProcessCategory(process)
            if self.getRegionFromBin(bin) == 'Pass':
                if self.tasks:
                    # return 'r_pass_%s_%s' % (task, cat)
                    return self.sf_name_prefix+cat+'_'+task
                else:
                    # return 'r_pass_%s' % cat
                    return self.sf_name_prefix+cat
            # elif 'fail' in bin.lower():
            elif self.getRegionFromBin(bin) == 'Fail' or self.getRegionFromBin(bin) == 'PassW' or self.getRegionFromBin(bin) == 'FailW':
                if self.tasks:
                    # return 'r_fail_%s_%s' % (task, cat)
                    return self.antisf_name_prefix+cat+'_'+task
                else:
                    # return 'r_fail_%s' % cat
                    return self.antisf_name_prefix+cat
            else:
                return 1
        else:
            return 1

lttModel = LegacyTopTaggingModel()
