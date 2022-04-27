# Working point study

Usage example for UL17:
- Go to `.../LegacyTopTagging/config`, run `./wpXmlCreator.py -y UL17`, go to the new subdirectory `config_WorkingPointStudy` and run SFrame on the newly created config XML
- After SFrame hadding is complete, manually hadd all QCD output into a file called `<your/path/to/LegacyTopTaggingOutput>/WorkingPointStudy/UL17/uhh2.AnalysisModuleRunner.MC.QCD_HT300toInf.root`
- Run `root -l -q -b 'restructure_root_trees.cxx("UL17")'` to produce ROOT files containing flat TTrees with all the jets needed for the WP study (may take some minutes)
- Run `python root_to_numpy.py -y UL17` to convert the TTrees into numpy format
- Run `python analyze.py -y UL17 -r` to calculate efficiencies vs. tau32 cuts. Uses ca. 15 GB RAM (depending on the size of the TTrees) and runs for ca. half an hour. There is the option to specify a threshold on the SoftDrop mass of studied jets: adding the parameter `-m 10` will reject jets with mSD < 10 GeV
- To ultimately get TGraphs stored in yet another set of ROOT files which can then be used for plotting, run `python analyze.py -y UL17` (same command as before, but without further arguments)
- Run `root -l -q -b 'plots.cxx("UL17")'` to produce efficiency and ROC plots. The actual working point analysis happens on the go within this script. The working points defined via targeted background efficiencies can be customized within the main function of this script
- Final plots are located at `<your/path/to/LegacyTopTaggingOutput>/WorkingPointStudy/UL17/workdir_npy/{**/,}plot*.pdf`
