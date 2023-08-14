import ROOT as root
import os

mscSplits = [
'FullyMerged',
'SemiMerged',
'NotMerged',
]

graph_types = [
'tot',
'stat',
'syst',
'uncorr',
'corr',
]

wps_AK8 = [
'BkgEff0p001',
'BkgEff0p005',
'BkgEff0p010',
'BkgEff0p025',
'BkgEff0p050',
]

wps = {
    'AK8': {
        'BkgEff0p001': 'wp0p38_vt',
        'BkgEff0p005': 'wp0p47_t',
        'BkgEff0p010': 'wp0p52_m',
        'BkgEff0p025': 'wp0p61_l',
        'BkgEff0p050': 'wp0p69_vl',
    },
    'HOTVR': {
        'Standard': 'default',
    },
}

wp_to_MisTagEffString = {
    'BkgEff0p001': 'mis0p001',
    'BkgEff0p005': 'mis0p005',
    'BkgEff0p010': 'mis0p010',
    'BkgEff0p025': 'mis0p025',
    'BkgEff0p050': 'mis0p050',
}

workdir_base = 'scaleFactors_2023-08-10/'
outdir = os.path.join(workdir_base, 'scaleFactors_for_TopTaggingScaleFactors_repo')
os.system('mkdir -p '+outdir)

years = [
'UL16preVFP',
'UL16postVFP',
'UL17',
'UL18',
]

for year in years:

    outname = 'TopTaggingScaleFactors_RunIISummer20{year}_PUPPIv15.root'.format(year=year)
    outpath = os.path.join(outdir, outname)
    outfile = root.TFile.Open(outpath, 'RECREATE')

    for wp in wps_AK8:
        foldername_out = 'AK8_PUPPI_'+wps['AK8'][wp]+'_'+wp_to_MisTagEffString[wp]
        outfolder = outfile.mkdir(foldername_out)
        inpath = os.path.join(workdir_base, 'scaleFactors-ak8_t__tau-'+wp, 'scaleFactors-ak8_t__tau-'+wp+'.root')
        foldername_in = 'ak8_t__tau-'+wp+'-'+year
        infile = root.TFile.Open(inpath, 'READ')
        for msc in mscSplits:
            for graph_type in graph_types:
                graph_name_in = msc+'_'+year+'_'+graph_type
                graph = infile.Get(foldername_in+'/'+graph_name_in)
                graph_name_out = msc+'_'+graph_type
                outfolder.cd()
                graph.Write(graph_name_out)

    for wp in wps_AK8:
        foldername_out = 'AK8_PUPPI_subjetbtag_'+wps['AK8'][wp]
        outfolder = outfile.mkdir(foldername_out)
        inpath = os.path.join(workdir_base, 'scaleFactors-ak8_t_btagDJet__tau-'+wp, 'scaleFactors-ak8_t_btagDJet__tau-'+wp+'.root')
        foldername_in = 'ak8_t_btagDJet__tau-'+wp+'-'+year
        infile = root.TFile.Open(inpath, 'READ')
        for msc in mscSplits:
            for graph_type in graph_types:
                graph_name_in = msc+'_'+year+'_'+graph_type
                graph = infile.Get(foldername_in+'/'+graph_name_in)
                graph_name_out = msc+'_'+graph_type
                outfolder.cd()
                graph.Write(graph_name_out)

    for wp in ['Standard']:
        foldername_out = 'HOTVR_PUPPI_'+wps['HOTVR'][wp]
        outfolder = outfile.mkdir(foldername_out)
        inpath = os.path.join(workdir_base, 'scaleFactors-hotvr_t__tau-'+wp, 'scaleFactors-hotvr_t__tau-'+wp+'.root')
        foldername_in = 'hotvr_t__tau-'+wp+'-'+year
        infile = root.TFile.Open(inpath, 'READ')
        for msc in mscSplits:
            for graph_type in graph_types:
                graph_name_in = msc+'_'+year+'_'+graph_type
                graph = infile.Get(foldername_in+'/'+graph_name_in)
                graph_name_out = msc+'_'+graph_type
                outfolder.cd()
                graph.Write(graph_name_out)

    outfile.Close()
