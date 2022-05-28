import sys

filename = sys.argv[-1]

systs = {}
with open(filename, 'r') as f:
    first_line = True
    for l in f.readlines():
        if first_line:
            first_line = False
            continue
        l = l.strip().split(',')
        # print(l)
        if l[1] != "mujets": continue
        # systs.setdefault(l[2].replace('up_', '').replace('down_', ''), {}).setdefault(l[3], {})
        systs.setdefault(l[2], {}).setdefault(l[3], {})

sorted_systs = {}
for i in range(len(systs)):
    key = sorted(systs.keys())[i]
    sorted_systs[key] = systs[key]

systs_ = {}
for k, v in sorted_systs.items():
    systs_[k] = list(v.keys())
    # if 'jes' in k: continue
    print(k+' :  '+', '.join([str(x) for x in systs_[k]]))
# systs = systs_
# print(systs)
