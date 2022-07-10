# LegacyTopTagging

UHH2 sub-framework for measuring top quark tagging efficiencies in AK8 and HOTVR jets for the Ultra Legacy datasets 2016, 2017, and 2018.

## Setup

- Please install UHH2 version `RunII_106X_v1`
- Once the UHH2 installation is done, clone this repository into `<your-path-to>/UHH2/`: <br />
`git clone https://github.com/MatthiesC/LegacyTopTagging.git` <br />
NB: This repo is not dependent on other UHH2 sub-frameworks like https://github.com/UHH2/HOTVR
- Please create a directory into which all output files should be stored, e.g. `/nfs/dust/cms/user/<username>/LegacyTopTaggingOutput/`; then create a symbolic link to that output directory within the framework: `ln -s /nfs/dust/cms/user/<username>/LegacyTopTaggingOutput/ <your-path-to>/UHH2/LegacyTopTagging/output`. This step is **obligatory**; this framework considers `<your-path-to>/UHH2/LegacyTopTagging/output` to be the base output directory -- since you do not want to spam your actual code environment with ultra-large output files, you should use this symbolic link
- Install the combine framework for the template fits:
  ```
  cd $CMSSW_BASE/src
  git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
  cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit
  git fetch origin
  git checkout v8.2.0
  cd $CMSSW_BASE/src
  git clone https://github.com/cms-analysis/CombineHarvester.git CombineHarvester
  scramv1 b clean; scramv1 b
  ```
- The CMSSW versions 10.6.X do not yet include the `uproot` package in their python3. Thus, you need to install Anaconda3 and `conda install -c conda-forge uproot`. (**FIXME**: Which python3 version should you install?) Declare a bash alias like `alias pyconda3=/path/to/your/ana3env/bin/python` and use this on python scripts that need `uproot`. **FIXME**: Right now, there seems to be a conflict between the CMSSW ROOT version and the ana3 ROOT version. Thus, you cannot do `import ROOT` with `pyconda3`

## Code conventions / tips

- Keep the code user-independent! Only declare relative paths, using e.g. `os.environ.get('CMSSW_BASE')` (python) or `std::getenv('CMSSW_BASE')` (C++). Note that these methods will be able to find the `CMSSW_BASE` only after you set up the environment (`cmsenv`)
- All SFrame-related C++ code (i.e. everything in `src/` and `include/`) is to be nested in the `uhh2::ltt` namespace
- All types of constants (working point cuts, mass cuts, etc.) are to be declared in `include/Constants.h`
- String comparisons within `uhh2::AnalysisModule::process()`, `uhh2::Selection::passes()`, and `uhh2::Hists::fill()` methods are to be avoided! Rather use enumerators and, if needed, map them to strings (cf. `include/Constants.h` for examples)
