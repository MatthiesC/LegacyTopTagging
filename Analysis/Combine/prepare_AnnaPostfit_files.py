from __future__ import print_function

import ROOT as root

import os
import sys
from copy import deepcopy
import math

tagger_algo = 'HOTVR'
#tagger_algo = 'AK8'

tagger_name = '{}_t__tau'.format(tagger_algo.lower())

print(tagger_name)

tau32_limit = 0.56 if tagger_algo == 'HOTVR' else 0.68


path_base = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel/run2/combine/{}/AnnaPostfit-{}-NullWP-pt_400toInf-run2-output_probejet_{}_tau32'.format(tagger_name, tagger_name, tagger_algo))
path_infile = os.path.join(path_base, 'Templates-Anna-{}-NullWP-pt_400toInf-run2-output_probejet_{}_tau32-both.root'.format(tagger_name, tagger_algo))
path_outfile = path_infile.replace('.root', '-postfit.root')

#print(path_infile)
#print(path_outfile)

channels = ['ele', 'muo']
years = ['UL16preVFP', 'UL16postVFP', 'UL17', 'UL18']
sf_directions = ['sf_central', 'sf_up', 'sf_down']
sf_types = ['corr', 'uncorr']
infiles_sf = {}
for year in years:
    for sf_direction in sf_directions:
        for sf_type in sf_types:
            path_infile_sf = os.path.join(
                path_base,
                'AnnaPostfit_BasicHists_'+year,
                'AnnaPostfit_BasicHists-'+sf_type+'-'+sf_direction,
                'BasicHists-{}-NullWP-pt_400toInf-{}-nominal-output_probejet_{}_tau32.root'.format(
                    tagger_name, year, tagger_algo
                )
            )
            asdf = root.TFile.Open(path_infile_sf, 'READ')
            root.SetOwnership(asdf, False)  # <-- Prevent ROOT from deleting the object
            infiles_sf.setdefault(year, {}).setdefault(sf_direction, {})[sf_type] = asdf

#print( infiles_sf['UL16preVFP']['sf_central']['corr'].Get('Main/Pass/ele/DATA__nominal').GetBinContent(20) )
#exit()

infile = root.TFile.Open(path_infile, "READ")

# Check if the file is open
if not infile or infile.IsZombie():
    print("Error: Could not open file")
    exit()

infile_hist_names = []
# Function to iterate over keys in the file and find TH1F histograms
def iterate_histograms(directory, recursive_path=None):
    for key in directory.GetListOfKeys():
        obj = key.ReadObj()
        if isinstance(obj, root.TH1):
            #print(f"Found TH1F: {obj.GetName()}")
            if recursive_path is not None:
                infile_hist_names.append('/'.join([recursive_path, key.GetName()]))
            else:
                infile_hist_names.append(key.GetName())
        elif isinstance(obj, root.TDirectory):
            if recursive_path is not None:
                recursive_path = '/'.join([recursive_path, key.GetName()])
            else:
                recursive_path = key.GetName()
            print("Entering directory:", recursive_path)
            iterate_histograms(obj, recursive_path=recursive_path)  # Recursive call for subdirectories

iterate_histograms(infile)
#print(infile_hist_names)

mscs = ["FullyMerged", "SemiMerged", "NotMerged"]

outfile = root.TFile.Open(path_outfile, 'RECREATE')
print('created', path_outfile)
outdir = outfile.mkdir('Main_Pass_both')

for infile_hist_name in infile_hist_names:
    outdir.cd()
    inhist = infile.Get(infile_hist_name)

    is_msc_hist = False
    for msc in mscs:
        #if "__MSc_"+msc in infile_hist_name and not (infile_hist_name.endswith(msc)):
        if "__MSc_"+msc in infile_hist_name:
            is_msc_hist = True

    is_jer = '_jer' in infile_hist_name

    if is_msc_hist and not is_jer:
        #inhist_nominal_name = '_'.join(infile_hist_name.split('_')[:-1])  # FIXME
        #inhist_nominal = infile.Get(inhist_nominal_name)
        inhist_nominal_name = infile_hist_name.split('_')[:-1] if infile_hist_name.endswith('Up') or infile_hist_name.endswith('Down') else infile_hist_name.split('_')
        inhist_nominal_name = '_'.join(inhist_nominal_name) + '__nominal'
        inhist_nominal_name = inhist_nominal_name.split('/')[-1]
        #print(inhist_nominal_name)
        #exit()
        sum_years = None
        for year in years:
            for channel in channels:
                if sum_years is None:
                    sum_years = deepcopy(infiles_sf[year]['sf_central']['uncorr'].Get('Main/Pass/'+channel+'/'+inhist_nominal_name))
                else:
                    sum_years.Add(
                        infiles_sf[year]['sf_central']['uncorr'].Get('Main/Pass/'+channel+'/'+inhist_nominal_name)
                    )
        for ibin in range(inhist.GetNbinsX()+2):
            if inhist.GetBinLowEdge(ibin) >= tau32_limit:
                continue
            inhist.SetBinContent(ibin, sum_years.GetBinContent(ibin))
            inhist.SetBinError(ibin, sum_years.GetBinError(ibin))

    if is_jer:
        inhist_nominal_name = '_'.join(infile_hist_name.split('_')[:-1]) + '__nominal'
        inhist_nominal_name = inhist_nominal_name.split('/')[-1]
        sf_direction = 'sf_up' if infile_hist_name.endswith('Up') else 'sf_down'
        for ibin in range(inhist.GetNbinsX()+2):
            if inhist.GetBinLowEdge(ibin) >= tau32_limit:
                continue
            binc = 0
            bine_corr = 0
            bine_uncorr = 0
            for year in years:
                sum_year = 0
                for channel in channels:
                    inhist_nominal_path = 'Main/Pass/'+channel+'/'+inhist_nominal_name
                    binc_temp = infiles_sf[year]['sf_central']['corr'].Get(inhist_nominal_path).GetBinContent(ibin)
                    bine_corr += 0 if binc_temp == 0 else abs(
                        binc_temp \
                        - infiles_sf[year][sf_direction]['corr'].Get(inhist_nominal_path).GetBinContent(ibin)
                    ) / binc_temp
                    binc += binc_temp  # corr/uncorr is irrelevant for sf_central
                    sum_year += binc_temp
                bine_uncorr_temp = math.sqrt(
                    (
                        abs(
                            infiles_sf[year]['sf_central']['uncorr'].Get('Main/Pass/ele/'+inhist_nominal_name).GetBinContent(ibin) \
                            - infiles_sf[year][sf_direction]['uncorr'].Get('Main/Pass/ele/'+inhist_nominal_name).GetBinContent(ibin)
                        ) \
                        + abs(
                            infiles_sf[year]['sf_central']['uncorr'].Get('Main/Pass/muo/'+inhist_nominal_name).GetBinContent(ibin) \
                            - infiles_sf[year][sf_direction]['uncorr'].Get('Main/Pass/muo/'+inhist_nominal_name).GetBinContent(ibin)
                        ) \
                    )**2 \
                    + (
                        abs(
                            infiles_sf[year]['sf_central']['corr'].Get('Main/Pass/ele/'+inhist_nominal_name).GetBinContent(ibin) \
                            - infiles_sf[year][sf_direction]['corr'].Get('Main/Pass/ele/'+inhist_nominal_name).GetBinContent(ibin)
                        ) \
                        + abs(
                            infiles_sf[year]['sf_central']['corr'].Get('Main/Pass/muo/'+inhist_nominal_name).GetBinContent(ibin) \
                            - infiles_sf[year][sf_direction]['corr'].Get('Main/Pass/muo/'+inhist_nominal_name).GetBinContent(ibin)
                        ) \
                    )**2
                )
                bine_uncorr += (bine_uncorr_temp / sum_year)**2 if sum_year != 0 else 0
            bine_corr = 0
            if sf_direction == 'sf_up':
                binc = binc * ( 1 + math.sqrt(bine_corr**2 + bine_uncorr) )
            else:
                binc = binc * ( 1 - math.sqrt(bine_corr**2 + bine_uncorr) )
            inhist.SetBinContent(ibin, binc)
            inhist.SetBinError(ibin, 0)

    inhist.Write(infile_hist_name.split('/')[-1])
    #print(inhist.GetName())
    #inhist.Write(inhist.GetName())
