# Example usage:
# python get_NoTopPtRw_numpys.py /data/dust/user/matthies/uhh2-106X-v2/CMSSW_10_6_28/src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel/run2/combine/ak8_w__partnet/ak8_w__partnet-NullWP-pt_200toInf-run2-output_probejet_AK8_pt-NoTopPtRw/BasicHists-ak8_w__partnet-NullWP-pt_200toInf-run2-nominal-output_probejet_AK8_pt.root

import numpy as np
import ROOT as root
import sys
from copy import deepcopy


infile_name = sys.argv[-1]
infile = root.TFile.Open(infile_name, 'READ')

processes = [
    'TTbar__nominal',
    'ST__nominal',
    'VJetsAndVV__nominal',
    'QCD__nominal',
]

hist_sum = None

for process in processes:
    for channel in ['muo', 'ele']:
        inhist = infile.Get('Main/Pass/'+channel+'/'+process)
        if hist_sum is None:
            hist_sum = deepcopy(inhist)
        else:
            hist_sum.Add(inhist)

contents = []
for i_bin in range(1, hist_sum.GetNbinsX()+1):
    contents.append(hist_sum.GetBinContent(i_bin))

print(contents)
contents = np.array(contents, dtype=np.float32)
print(contents)

tagger = infile_name.split('/')[-1].split('-')[1]
outfile_name = 'get_NoTopPtRw_numpys_'+tagger+'.npy'
with open(outfile_name, 'wb') as outfile:
    np.save(outfile, contents)
