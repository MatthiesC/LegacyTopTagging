from collections import OrderedDict

# Source: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopSystematics#Luminosity
_YEARS = OrderedDict([
    ('UL16', {
        'short_name': 'UL16',
        'year': '2016',
        'long_name': 'Ultra Legacy 2016',
        'lumi_fb': 36.33,
        'lumi_pb': 36330.,
        'lumi_unc': 0.012,
    }),
    ('UL16preVFP', {
        'short_name': 'UL16preVFP',
        'year': '2016 (early)',
        'long_name': 'Ultra Legacy 2016 (early)',
        'lumi_fb': 19.52,
        'lumi_pb': 19520.,
        'lumi_unc': 0.012,
    }),
    ('UL16postVFP', {
        'short_name': 'UL16postVFP',
        'year': '2016 (late)',
        'long_name': 'Ultra Legacy 2016 (late)',
        'lumi_fb': 16.81,
        'lumi_pb': 16810.,
        'lumi_unc': 0.012,
    }),
    ('UL17', {
        'short_name': 'UL17',
        'year': '2017',
        'long_name': 'Ultra Legacy 2017',
        'lumi_fb': 41.48,
        'lumi_pb': 41480.,
        'lumi_unc': 0.023,
    }),
    ('UL18', {
        'short_name': 'UL18',
        'year': '2018',
        'long_name': 'Ultra Legacy 2018',
        'lumi_fb': 59.83,
        'lumi_pb': 59830.,
        'lumi_unc': 0.025,
    }),
    ('ULRunII', {
        'short_name': 'ULRunII',
        'year': '2016 + 2017 + 2018',
        'long_name': 'Run II Ultra Legacy',
        'lumi_fb': 137.65,
        'lumi_pb': 137650.,
        'lumi_unc': 0.016,
    }),
])
