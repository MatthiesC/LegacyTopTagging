years = ['UL18']
# years = ['UL16preVFP', 'UL16postVFP', 'UL17', 'UL18']

commands = list()

for year in years:
    # commands.append('''root -l -q -b 'restructure_root_trees.cxx("'''+year+'''", "qcd_ak8_200")' ''')
    # commands.append('''root -l -q -b 'restructure_root_trees.cxx("'''+year+'''", "qcd_ak8_300")' ''')
    # commands.append('''root -l -q -b 'restructure_root_trees.cxx("'''+year+'''", "qcd_hotvr_200")' ''')
    commands.append('''root -l -q -b 'restructure_root_trees.cxx("'''+year+'''", "wjets_ak8_200")' ''')
    commands.append('''root -l -q -b 'restructure_root_trees.cxx("'''+year+'''", "ttbar")' ''')

for c in commands:
    c = c.replace("\'", r"'")
    print c

# exit()

import sys
sys.path.append('..')
from parallel_threading import run_with_pool

# run_with_pool(commands=commands)
run_with_pool(commands=commands, max_workers=len(commands)) # careful with this!
