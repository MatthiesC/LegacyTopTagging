# Working point study

Usage example for UL17:
- Go to `.../LegacyTopTagging/config`, run `./wpXmlCreator.py -y UL17`, go to the new subdirectory `config_WorkingPointStudy` and run SFrame on the newly created config XML
- After SFrame hadding is complete, manually hadd all QCD and WJetsToQQ HT bin outputs respectively into a single file called `<your/path/to/LegacyTopTaggingOutput>/WorkingPointStudy/UL17/nominal/uhh2.AnalysisModuleRunner.MC.<QCD,WJetsToQQ>_HT200toInf.root` <br />
  Being in `.../LegacyTopTagging/output/WorkingPointStudy` you can use <br />
  `for i in UL*; do hadd ${i}/nominal/uhh2.AnalysisModuleRunner.MC.WJetsToQQ_HT200toInf_${i}.root ${i}/nominal/uhh2.AnalysisModuleRunner.MC.WJetsToQQ_HT*_${i}.root; done;` <br />
  `for i in UL*; do hadd ${i}/nominal/uhh2.AnalysisModuleRunner.MC.QCD_HT200toInf_${i}.root ${i}/nominal/uhh2.AnalysisModuleRunner.MC.QCD_HT*_${i}.root; done;` <br />
  `shopt -s extglob; for i in UL*; do hadd -f ${i}/nominal/uhh2.AnalysisModuleRunner.MC.QCD_HT300toInf_${i}.root ${i}/nominal/uhh2.AnalysisModuleRunner.MC.QCD_HT!(200to*)_${i}.root; done; shopt -u extglob;`
- Run `root -l -q -b 'restructure_root_trees.cxx("UL17")'` to produce ROOT files containing flat TTrees with all the jets needed for the WP study (may take some minutes)
- [DEPRECATED: now using `uproot`!] ~~Run `python root_to_numpy.py -y UL17` to convert the TTrees into numpy format~~
- Run `pyconda3 analyze.py -y UL17 -r` to calculate efficiencies vs. tau32 cuts (or vs. other variables). Uses O(50) GB RAM (depending on the size of the TTrees) and runs for ca. an hour; `pyconda3` in an alias to a python executable that does support the `uproot` package (you probably need to install Anaconda3 for this first)
- To ultimately get TGraphs stored in yet another set of ROOT files which can then be used for plotting, run `python analyze.py -y UL17` (same command as before, but without further arguments); here `python` is again just the python executable that comes with CMSSW

- Run `root -l -q -b 'plots.cxx("UL17", "ak8_t__tau")'` to produce efficiency and ROC plots. The actual working point analysis happens on the go within this script. The working points defined via targeted background efficiencies can be customized within the main function of this script
- Final plots are located at `<your/path/to/LegacyTopTaggingOutput>/WorkingPointStudy/UL17/nominal/ak8_t__tau/{**/,}plot*.pdf`
