# LegacyTopTagging

UHH2 sub-framework for measuring top quark tagging efficiencies in AK8 and HOTVR jets for the Ultra Legacy datasets 2016, 2017, and 2018.

## Setup

- Please install UHH2 version `RunII_106X_v1`
- Once the UHH2 installation is done, clone this repository into `<your-path-to>/UHH2/`: <br />
`git clone https://github.com/MatthiesC/LegacyTopTagging.git` <br />
NB: This repository is not dependent on other UHH2 sub-frameworks like `UHH2/HOTVR` (https://github.com/UHH2/HOTVR.git)!
- Please create a directory into which all output files should be stored, e.g. `/nfs/dust/cms/user/<username>/LegacyTopTaggingOutput/`; then create a symbolic link to that output directory within the framework: `ln -s /nfs/dust/cms/user/<username>/LegacyTopTaggingOutput/ <your-path-to>/UHH2/LegacyTopTagging/output`. This step is **obligatory**; this framework considers `<your-path-to>/UHH2/LegacyTopTagging/output` to be the base output directory -- since you do not want to spam your framework with ultra-large output files, you should use this symbolic link
