import os
from subprocess import call

import sys
sys.path.append(os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/Analysis'))
from constants import Systematics

# systematics = Systematics(blacklist=['sfelec', 'sfmu_iso'])
# systematics = Systematics(blacklist=['sfmu_iso'])
systematics = Systematics(blacklist=['sfelec_trigger', 'sfmu_iso'])
# systematics = Systematics(whitelist=['mtop'])
# systematics = Systematics(whitelist=['murmuf'])
_SYSTEMATICS = systematics.get_all_variations()
# del _SYSTEMATICS['nominal']
print(_SYSTEMATICS.keys())
# sys.exit()

class JobSubmitter():

    def __init__(self, tagger_name, wp_index, year, only_nominal=False):

        self.workdir_name = 'workdir_uproot'
        self.workdir = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/Analysis/Combine', self.workdir_name)
        os.system('mkdir -p '+self.workdir)
        self.tagger_name = tagger_name
        self.wp_index = wp_index
        self.year = year
        # self.systematics_queue = '\n'.join(['nominal', 'pu_down']) # FIXME
        if only_nominal:
            self.systematics_queue = 'nominal'
        else:
            self.systematics_queue = '\n'.join(_SYSTEMATICS.keys()) # FIXME
        self.batch_name = '-'.join(['UpRoot', self.tagger_name, 'wp'+str(self.wp_index), self.year])
        self.submit_file_name = 'Condor.'+self.batch_name+'.submit'
        self.submit_file_path = os.path.join(self.workdir, self.submit_file_name)

    def write_submit_file(self):

        with open(self.submit_file_path, 'w') as submit_file:
            submit_file.write("""#HTC Submission File
# +MyProject        =  "af-cms"
Requirements = ( OpSysAndVer == "CentOS7" )
universe          = vanilla
notification      = Error
notify_user       = """+os.environ.get('USER')+"""@mail.desy.de
output            = """+self.workdir+"""/"""+self.batch_name+""".$(Systematic).$(ClusterId).$(Process).out
error             = """+self.workdir+"""/"""+self.batch_name+""".$(Systematic).$(ClusterId).$(Process).err
log               = """+self.workdir+"""/"""+self.batch_name+""".$(Cluster).log
RequestCpus       = 1
RequestMemory     = 4G
RequestDisk       = 4G
getenv            = True
JobBatchName      = """+self.batch_name+"""
executable        = create_root_files_for_datacards_uproot.py
arguments         = -t """+self.tagger_name+""" -w """+str(self.wp_index)+""" -y """+self.year+""" -s $(Systematic)
queue Systematic from (
"""+self.systematics_queue+"""
)
            """)

    def submit(self, dryrun=False):

        call("condor_submit {} {}".format(self.submit_file_path, '--dry-run {}/dryrun.log'.format(self.workdir) if dryrun else ''), shell=True)


if __name__=='__main__':

    years = [
    'UL16preVFP',
    'UL16postVFP',
    'UL17',
    'UL18',
    ]

    ## ADJUST HERE WHICH TAGGER/WP YOU WANT TO SUBMIT
    ## Then just do `python submit_uproot_jobs.py`
    # key = tagger_name, value = number of wps
    # don't submit all at once! Will go over 5,000 jobs limit! (With all systs incl. all JES splits, FSR/ISR splits, one WP creates 146 jobs/year)
    taggers = {
        'ak8_t__tau': 5,
        # 'ak8_t_btagDJet__tau': 5,
        # 'ak8_t_btagDCSV__tau': 5,
        # 'hotvr_t__tau': 1,
        # 'ak8_w__partnet': 1,
        # 'ak8_t__MDdeepak8': 1,
    }

    null_wp = False
    # null_wp = True

    for tagger_name, n_wps in taggers.items():
        if null_wp:
            for year in years:
                js = JobSubmitter(tagger_name, 999, year)
                js.write_submit_file()
                # js.submit(dryrun=True)
                js.submit()
        else:
            for wp_index in range(n_wps):
                for year in years:
                    js = JobSubmitter(tagger_name, wp_index, year)
                    js.write_submit_file()
                    # js.submit(dryrun=True)
                    js.submit()
                # if null_wp==True:
                #     break

    # #____________
    # ## For submitting individual job batches:
    # js = JobSubmitter('ak8_t_btagDJet__tau', 1, 'UL16postVFP')
    # js.write_submit_file()
    # # js.submit(dryrun=True)
    # js.submit()
