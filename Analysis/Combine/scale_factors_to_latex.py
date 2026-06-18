from __future__ import print_function

import os
import sys
import ROOT as root

# argv 1: ak8_t__tau or other
# argv 2: BkgEff0p001 or other

tagger = sys.argv[1]
wp = sys.argv[2]
cat = sys.argv[3]

cats = {
    "fm": "FullyMerged",
    "sm": "SemiMerged",
    "nm": "NotMerged",
    "wm": "WMerged",
    "xm": "NotTopOrWMerged",
}

infile = os.path.join(
    "scaleFactors_2023-08-10",
    "scaleFactors-"+tagger+"-"+wp,
    "scaleFactors-"+tagger+"-"+wp+".root",
)

years = [
    "UL16preVFP",
    "UL16postVFP",
    "UL17",
    "UL18",
]

gtypes = [
    "tot",
    "stat",
    "syst",
    "uncorr",
    "corr",
]

class PtBin:
    def __init__(self, pt_low, pt_high=None):
        self.pt_low = pt_low
        self.pt_high = pt_high or None
    @property
    def pt_central(self):
        return (self.pt_high + self.pt_low) * 0.5 if self.pt_high is not None else self.pt_low + 1.

if tagger.startswith("ak8_t"):
    pt_bins = [
        PtBin(300,400),
        PtBin(400,480),
        PtBin(480,600),
        PtBin(600),
    ]
elif tagger.startswith("ak8_w"):
    pt_bins = [
        PtBin(200,250),
        PtBin(250,300),
        PtBin(300,350),
        PtBin(350,400),
        PtBin(400),
    ]
elif tagger.startswith("hotvr_t"):
    pt_bins = [
        PtBin(200,250),
        PtBin(250,300),
        PtBin(300,350),
        PtBin(350,400),
        PtBin(400,480),
        PtBin(480,600),
        PtBin(600),
    ]

if not os.path.isfile(infile):
    sys.exit("File does not exist")

print("Using file:", infile)

infile = root.TFile.Open(infile, "READ")
graphs = {}
for year in years:
    for gt in gtypes:
        graphs.setdefault(year, {})[gt] = infile.Get(tagger+"-"+wp+"-"+year+"/"+cats[cat]+"_"+year+"_"+gt)

for comb in [["stat", "syst"], ["uncorr", "corr"]]:
    print()
    for i_bin, pt_bin in enumerate(pt_bins):
        line = "\quad{}"
        if pt_bin.pt_high is not None:
            line += str(int(pt_bin.pt_low))+"--"+str(int(pt_bin.pt_high))
        else:
            line += "$>"+str(int(pt_bin.pt_low))+"$"
        for year in years:
            x_value = root.Double()
            y_value = root.Double()
            graphs[year]["tot"].GetPoint(i_bin, x_value, y_value)
            v1 = float(y_value)
            v2 = float(graphs[year][comb[0]].GetErrorYhigh(i_bin))
            v3 = float(graphs[year][comb[0]].GetErrorYlow(i_bin))
            v4 = float(graphs[year][comb[1]].GetErrorYhigh(i_bin))
            v5 = float(graphs[year][comb[1]].GetErrorYlow(i_bin))
            line += " & \scaleFactor{"
            line += "{:.3f}".format(v1)
            line += "}"
            for v in [v2,v3,v4,v5]:
                line += "{"
                line += "{:.3f}".format(v)[1:] if v < 1. else "{:.3f}".format(v)
                line += "}"
        line += " \\\\"
        print(line)
        # line = """\quad{}#1 & \scaleFactor{1.000}{.002}{.004}{.002}{.004}"""
