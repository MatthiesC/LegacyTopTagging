# Template Fits

All output goes to `<your/path/to/LegacyTopTaggingOutput>/TagAndProbe/mainsel/UL17/Combine/`

Usage example:
- Run `python create_root_files_for_datacards.py` to read-in and re-organize histograms of the ProbeJetHists classes. It also creates a text file (see next step) with commands to plot prefit histograms
- Run `make` to compile all ROOT macros placed in `src/`. This will give a new binary `bin/plots`, so far used to plot prefit histograms
- Run `python plots.py workdir/workfile_plotting_commands_prefit.txt` to get prefit histograms


**NEW usage example:**

All output goes to `<your/path/to/LegacyTopTaggingOutput>/TagAndProbe/mainsel/UL17/combine/` and similar

- Adjust settings in `submit_uproot_jobs.py` and then run `python submit_uproot_jobs.py` (choose which tagger to run, how many WPs, which years...)
- If you want to run specific tagger/year/wp etc. locally, you can just run `pyconda3 create_root_files_for_datacards_uproot.py` with appropriate parameters (you can also choose different variable than jet mass)
- Rearrange histograms in combine-friendly format: `python rearrange_basic_hists_from_uproot.py` (need to adjust settings in the file, no argparse implemented)
