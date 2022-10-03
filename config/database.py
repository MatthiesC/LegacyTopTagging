samplesDict = {
    #_________________________________________________
    'DATA_SingleMuon_A': {
        'db_name': 'SingleMuon_RunA',
        'years': ['UL18'],
        'channel': ['muo'],
    },
    'DATA_SingleMuon_B': {
        'db_name': 'SingleMuon_RunB',
        'years': ['UL16preVFP', 'UL17', 'UL18'],
        'channel': ['muo'],
    },
    'DATA_SingleMuon_C': {
        'db_name': 'SingleMuon_RunC',
        'years': ['UL16preVFP', 'UL17', 'UL18'],
        'channel': ['muo'],
    },
    'DATA_SingleMuon_D': {
        'db_name': 'SingleMuon_RunD',
        'years': ['UL16preVFP', 'UL17', 'UL18'],
        'channel': ['muo'],
    },
    'DATA_SingleMuon_E': {
        'db_name': 'SingleMuon_RunE',
        'years': ['UL16preVFP', 'UL17'],
        'channel': ['muo'],
    },
    'DATA_SingleMuon_F': {
        'db_name': 'SingleMuon_RunF',
        'years': ['UL16preVFP', 'UL16postVFP', 'UL17'],
        'channel': ['muo'],
    },
    'DATA_SingleMuon_G': {
        'db_name': 'SingleMuon_RunG',
        'years': ['UL16postVFP'],
        'channel': ['muo'],
    },
    'DATA_SingleMuon_H': {
        'db_name': 'SingleMuon_RunH',
        'years': ['UL16postVFP'],
        'channel': ['muo'],
    },
    #_________________________________________________
    # SingleElectron
    'DATA_SingleElectron_B': {
        'db_name': 'SingleElectron_RunB',
        'years': ['UL16preVFP', 'UL17'],
        'channel': ['ele'],
    },
    'DATA_SingleElectron_C': {
        'db_name': 'SingleElectron_RunC',
        'years': ['UL16preVFP', 'UL17'],
        'channel': ['ele'],
    },
    'DATA_SingleElectron_D': {
        'db_name': 'SingleElectron_RunD',
        'years': ['UL16preVFP', 'UL17'],
        'channel': ['ele'],
    },
    'DATA_SingleElectron_E': {
        'db_name': 'SingleElectron_RunE',
        'years': ['UL16preVFP', 'UL17'],
        'channel': ['ele'],
    },
    'DATA_SingleElectron_F': {
        'db_name': 'SingleElectron_RunF',
        'years': ['UL16preVFP', 'UL16postVFP', 'UL17'],
        'channel': ['ele'],
    },
    'DATA_SingleElectron_G': {
        'db_name': 'SingleElectron_RunG',
        'years': ['UL16postVFP'],
        'channel': ['ele'],
    },
    'DATA_SingleElectron_H': {
        'db_name': 'SingleElectron_RunH',
        'years': ['UL16postVFP'],
        'channel': ['ele'],
    },
    #_________________________________________________
    # SinglePhoton
    'DATA_SinglePhoton_B': {
        'db_name': 'SinglePhoton_RunB',
        'years': ['UL16preVFP', 'UL17'],
        'channel': ['ele'],
    },
    'DATA_SinglePhoton_C': {
        'db_name': 'SinglePhoton_RunC',
        'years': ['UL16preVFP', 'UL17'],
        'channel': ['ele'],
    },
    'DATA_SinglePhoton_D': {
        'db_name': 'SinglePhoton_RunD',
        'years': ['UL16preVFP', 'UL17'],
        'channel': ['ele'],
    },
    'DATA_SinglePhoton_E': {
        'db_name': 'SinglePhoton_RunE',
        'years': ['UL16preVFP', 'UL17'],
        'channel': ['ele'],
    },
    'DATA_SinglePhoton_F': {
        'db_name': 'SinglePhoton_RunF',
        'years': ['UL16preVFP', 'UL16postVFP', 'UL17'],
        'channel': ['ele'],
    },
    'DATA_SinglePhoton_G': {
        'db_name': 'SinglePhoton_RunG',
        'years': ['UL16postVFP'],
        'channel': ['ele'],
    },
    'DATA_SinglePhoton_H': {
        'db_name': 'SinglePhoton_RunH',
        'years': ['UL16postVFP'],
        'channel': ['ele'],
    },
    #_________________________________________________
    # EGamma
    'DATA_EGamma_A': {
        'db_name': 'EGamma_RunA',
        'years': ['UL18'],
        'channel': ['ele'],
    },
    'DATA_EGamma_B': {
        'db_name': 'EGamma_RunB',
        'years': ['UL18'],
        'channel': ['ele'],
    },
    'DATA_EGamma_C': {
        'db_name': 'EGamma_RunC',
        'years': ['UL18'],
        'channel': ['ele'],
    },
    'DATA_EGamma_D': {
        'db_name': 'EGamma_RunD',
        'years': ['UL18'],
        'channel': ['ele'],
    },
    #_________________________________________________
    'TTbarTo2L2Nu': {
        'db_name': 'TTTo2L2Nu',
        'extra_systs': {
            'hdamp_down': 'TTTo2L2Nu_hdampDOWN',
            'hdamp_up': 'TTTo2L2Nu_hdampUP',
            'tune_down': 'TTTo2L2Nu_TuneCP5down',
            'tune_up': 'TTTo2L2Nu_TuneCP5up',
            'mtop_mtop171p5': 'TTTo2L2Nu_mtop171p5',
            'mtop_mtop173p5': 'TTTo2L2Nu_mtop173p5',
            'cr_cr1': 'TTTo2L2Nu_CR1',
            'cr_cr2': 'TTTo2L2Nu_CR2',
            'cr_erdon': 'TTTo2L2Nu_erdON',
        },
    },
    'TTbarToSemiLeptonic': {
        'db_name': 'TTToSemiLeptonic',
        'extra_systs': {
            'hdamp_down': 'TTToSemiLeptonic_hdampDOWN',
            'hdamp_up': 'TTToSemiLeptonic_hdampUP',
            'tune_down': 'TTToSemiLeptonic_TuneCP5down',
            'tune_up': 'TTToSemiLeptonic_TuneCP5up',
            'mtop_mtop171p5': 'TTToSemiLeptonic_mtop171p5',
            'mtop_mtop173p5': 'TTToSemiLeptonic_mtop173p5',
            'cr_cr1': 'TTToSemiLeptonic_CR1',
            'cr_cr2': 'TTToSemiLeptonic_CR2',
            'cr_erdon': 'TTToSemiLeptonic_erdON',
        },
    },
    'TTbarToHadronic': {
        'db_name': 'TTToHadronic',
        'analysis': ['wp', 'sf', 'tw'],
        'extra_systs': {
            'hdamp_down': 'TTToHadronic_hdampDOWN',
            'hdamp_up': 'TTToHadronic_hdampUP',
            'tune_down': 'TTToHadronic_TuneCP5down',
            'tune_up': 'TTToHadronic_TuneCP5up',
            'mtop_mtop171p5': 'TTToHadronic_mtop171p5',
            'mtop_mtop173p5': 'TTToHadronic_mtop173p5',
            'cr_cr1': 'TTToHadronic_CR1',
            'cr_cr2': 'TTToHadronic_CR2',
            'cr_erdon': 'TTToHadronic_erdON',
        },
    },
    #_________________________________________________
    'TTbarMtt700to1000': {
        'db_name': 'TT_Mtt-700to1000',
        # 'analysis': ['tw'],
        'corr': True,
    },
    'TTbarMtt1000toInf': {
        'db_name': 'TT_Mtt-1000toInf',
        # 'analysis': ['tw'],
        'corr': True,
    },
    #_________________________________________________
    'ST_tW_DR_inclusiveDecays_T': {
        'db_name': 'ST_tW_top_5f_inclusiveDecays',
    },
    'ST_tW_DR_inclusiveDecays_Tbar': {
        'db_name': 'ST_tW_antitop_5f_inclusiveDecays',
    },
    'ST_tW_DR_NoFullyHadronic_T': {
        'db_name': 'ST_tW_top_5f_NoFullyHadronicDecays',
    },
    'ST_tW_DR_NoFullyHadronic_Tbar': {
        'db_name': 'ST_tW_antitop_5f_NoFullyHadronicDecays',
    },
    'ST_tW_DR_NoFullyHadronic_PDFWeights_T': {
        'db_name': 'ST_tW_top_5f_NoFullyHadronicDecays_PDFWeights',
    },
    'ST_tW_DR_NoFullyHadronic_PDFWeights_Tbar': {
        'db_name': 'ST_tW_antitop_5f_NoFullyHadronicDecays_PDFWeights',
    },
    'ST_tW_DS_NoFullyHadronic_T': {
        'db_name': 'ST_tW_top_5f_DS_NoFullyHadronicDecays',
        'analysis': ['tw'],
    },
    'ST_tW_DS_NoFullyHadronic_Tbar': {
        'db_name': 'ST_tW_antitop_5f_DS_NoFullyHadronicDecays',
        'analysis': ['tw'],
    },
    #_________________________________________________
    'ST_tChannel_T': {
        'db_name': 'ST_t-channel_top_4f_InclusiveDecays',
    },
    'ST_tChannel_Tbar': {
        'db_name': 'ST_t-channel_antitop_4f_InclusiveDecays',
    },
    #_________________________________________________
    'ST_sChannel_leptonDecays': {
        'db_name': 'ST_s-channel_4f_leptonDecays',
    },
    #_________________________________________________
    'WJetsToLNu_HT70to100': {
        'db_name': 'WJetsToLNu_HT-70to100',
    },
    'WJetsToLNu_HT100to200': {
        'db_name': 'WJetsToLNu_HT-100to200',
    },
    'WJetsToLNu_HT200to400': {
        'db_name': 'WJetsToLNu_HT-200to400',
    },
    'WJetsToLNu_HT400to600': {
        'db_name': 'WJetsToLNu_HT-400to600',
    },
    'WJetsToLNu_HT600to800': {
        'db_name': 'WJetsToLNu_HT-600to800',
    },
    'WJetsToLNu_HT800to1200': {
        'db_name': 'WJetsToLNu_HT-800to1200',
    },
    'WJetsToLNu_HT1200to2500': {
        'db_name': 'WJetsToLNu_HT-1200to2500',
    },
    'WJetsToLNu_HT2500toInf': {
        'db_name': 'WJetsToLNu_HT-2500toInf',
    },
    #_________________________________________________
    'DYJetsToLL_HT70to100': {
        'db_name': 'DYJetsToLL_M-50_HT-70to100',
    },
    'DYJetsToLL_HT100to200': {
        'db_name': 'DYJetsToLL_M-50_HT-100to200',
    },
    'DYJetsToLL_HT200to400': {
        'db_name': 'DYJetsToLL_M-50_HT-200to400',
    },
    'DYJetsToLL_HT400to600': {
        'db_name': 'DYJetsToLL_M-50_HT-400to600',
    },
    'DYJetsToLL_HT600to800': {
        'db_name': 'DYJetsToLL_M-50_HT-600to800',
    },
    'DYJetsToLL_HT800to1200': {
        'db_name': 'DYJetsToLL_M-50_HT-800to1200',
    },
    'DYJetsToLL_HT1200to2500': {
        'db_name': 'DYJetsToLL_M-50_HT-1200to2500',
    },
    'DYJetsToLL_HT2500toInf': {
        'db_name': 'DYJetsToLL_M-50_HT-2500toInf',
    },
    #_________________________________________________
    'Diboson_WW': {
        'db_name': 'WW',
    },
    'Diboson_WZ': {
        'db_name': 'WZ',
    },
    'Diboson_ZZ': {
        'db_name': 'ZZ',
    },
    #_________________________________________________
    'QCD_Mu_Pt15to20': {
        'db_name': 'QCD_Pt-15To20_MuEnrichedPt5',
        'channel': ['muo'],
    },
    'QCD_Mu_Pt20to30': {
        'db_name': 'QCD_Pt-20To30_MuEnrichedPt5',
        'channel': ['muo'],
    },
    'QCD_Mu_Pt30to50': {
        'db_name': 'QCD_Pt-30To50_MuEnrichedPt5',
        'channel': ['muo'],
    },
    'QCD_Mu_Pt50to80': {
        'db_name': 'QCD_Pt-50To80_MuEnrichedPt5',
        'channel': ['muo'],
    },
    'QCD_Mu_Pt80to120': {
        'db_name': 'QCD_Pt-80To120_MuEnrichedPt5',
        'channel': ['muo'],
    },
    'QCD_Mu_Pt120to170': {
        'db_name': 'QCD_Pt-120To170_MuEnrichedPt5',
        'channel': ['muo'],
    },
    'QCD_Mu_Pt170to300': {
        'db_name': 'QCD_Pt-170To300_MuEnrichedPt5',
        'channel': ['muo'],
    },
    'QCD_Mu_Pt300to470': {
        'db_name': 'QCD_Pt-300To470_MuEnrichedPt5',
        'channel': ['muo'],
    },
    'QCD_Mu_Pt470to600': {
        'db_name': 'QCD_Pt-470To600_MuEnrichedPt5',
        'channel': ['muo'],
    },
    'QCD_Mu_Pt600to800': {
        'db_name': 'QCD_Pt-600To800_MuEnrichedPt5',
        'channel': ['muo'],
    },
    'QCD_Mu_Pt800to1000': {
        'db_name': 'QCD_Pt-800To1000_MuEnrichedPt5',
        'channel': ['muo'],
    },
    'QCD_Mu_Pt1000toInf': {
        'db_name': 'QCD_Pt-1000_MuEnrichedPt5',
        'channel': ['muo'],
    },
    #_________________________________________________
    'QCD_EM_Pt15to20': {
        'db_name': 'QCD_Pt-15to20_EMEnriched',
        # 'analysis': ['tw'],
        'channel': ['ele'],
    },
    'QCD_EM_Pt20to30': {
        'db_name': 'QCD_Pt-20to30_EMEnriched',
        # 'analysis': ['tw'],
        'channel': ['ele'],
    },
    'QCD_EM_Pt30to50': {
        'db_name': 'QCD_Pt-30to50_EMEnriched',
        # 'analysis': ['tw'],
        'channel': ['ele'],
    },
    'QCD_EM_Pt50to80': {
        'db_name': 'QCD_Pt-50to80_EMEnriched',
        # 'analysis': ['tw'],
        'channel': ['ele'],
    },
    'QCD_EM_Pt80to120': {
        'db_name': 'QCD_Pt-80to120_EMEnriched',
        # 'analysis': ['tw'],
        'channel': ['ele'],
    },
    'QCD_EM_Pt120to170': {
        'db_name': 'QCD_Pt-120to170_EMEnriched',
        # 'analysis': ['tw'],
        'channel': ['ele'],
    },
    'QCD_EM_Pt170to300': {
        'db_name': 'QCD_Pt-170to300_EMEnriched',
        # 'analysis': ['tw'],
        'channel': ['ele'],
    },
    'QCD_EM_Pt300toInf': {
        'db_name': 'QCD_Pt-300toInf_EMEnriched',
        # 'analysis': ['tw'],
        'channel': ['ele'],
    },
    #_________________________________________________
    'QCD_bcToE_Pt15to20': {
        'db_name': 'QCD_Pt_15to20_bcToE',
        # 'analysis': ['tw'],
        'channel': ['ele'],
    },
    'QCD_bcToE_Pt20to30': {
        'db_name': 'QCD_Pt_20to30_bcToE',
        # 'analysis': ['tw'],
        'channel': ['ele'],
    },
    'QCD_bcToE_Pt30to80': {
        'db_name': 'QCD_Pt_30to80_bcToE',
        # 'analysis': ['tw'],
        'channel': ['ele'],
    },
    'QCD_bcToE_Pt80to170': {
        'db_name': 'QCD_Pt_80to170_bcToE',
        # 'analysis': ['tw'],
        'channel': ['ele'],
    },
    'QCD_bcToE_Pt170to250': {
        'db_name': 'QCD_Pt_170to250_bcToE',
        # 'analysis': ['tw'],
        'channel': ['ele'],
    },
    'QCD_bcToE_Pt250toInf': {
        'db_name': 'QCD_Pt_250toInf_bcToE',
        # 'analysis': ['tw'],
        'channel': ['ele'],
    },
    #_________________________________________________
    'QCD_HT200to300': {
        'db_name': 'QCD_HT200to300',
        'analysis': ['wp'],
    },
    'QCD_HT300to500': {
        'db_name': 'QCD_HT300to500',
        'analysis': ['wp'],
    },
    'QCD_HT500to700': {
        'db_name': 'QCD_HT500to700',
        'analysis': ['wp'],
    },
    'QCD_HT700to1000': {
        'db_name': 'QCD_HT700to1000',
        'analysis': ['wp'],
    },
    'QCD_HT1000to1500': {
        'db_name': 'QCD_HT1000to1500',
        'analysis': ['wp'],
    },
    'QCD_HT1500to2000': {
        'db_name': 'QCD_HT1500to2000',
        'analysis': ['wp'],
    },
    'QCD_HT2000toInf': {
        'db_name': 'QCD_HT2000toInf',
        'analysis': ['wp'],
    },
    #_________________________________________________
    'WJetsToQQ_HT200to400': {
        'db_name': 'WJetsToQQ_HT200to400',
        'analysis': ['wp'],
    },
    'WJetsToQQ_HT400to600': {
        'db_name': 'WJetsToQQ_HT400to600',
        'analysis': ['wp'],
    },
    'WJetsToQQ_HT600to800': {
        'db_name': 'WJetsToQQ_HT600to800',
        'analysis': ['wp'],
    },
    'WJetsToQQ_HT800toInf': {
        'db_name': 'WJetsToQQ_HT800toInf',
        'analysis': ['wp'],
    },
}
