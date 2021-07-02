# Template Fits

All output goes to `<your/path/to/LegacyTopTaggingOutput>/TagAndProbe/mainsel/UL17/Combine/`

Usage example:
- Run `python create_root_files_for_datacards.py` to read-in and re-organize histograms of the ProbeJetHists classes. It also creates a text file (see next step) with commands to plot prefit histograms
- Run `make` to compile all ROOT macros placed in `src/`. This will give a new binary `bin/plots`, so far used to plot prefit histograms
- Run `python plots.py workdir/workfile_plotting_commands_prefit.txt` to get prefit histograms
