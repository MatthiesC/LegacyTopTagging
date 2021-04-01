# WorkingPointStudy

- Run SFrame on `config/WorkingPointStudyUL17.xml`
- Run `root -l -q -b restructure_root_trees.cxx` to produce ROOT files containing flat TTrees with all the jets needed for the WP study (may take some minutes)
- Run `python root_to_numpy.py` to convert the TTrees into numpy format
- Run `python analyze.py` to calculate efficiencies vs. tau32 cuts and to ultimately get TGraphs stored in yet another set of ROOT files which can then be used for plotting
- Run `root -l -q -b plots.cxx` to produce efficiency and ROC plots
