- `pyconda3 step1_write_hists.py -y UL18 -c muo -t ak8_t__tau` or: `pyconda3 step1_write_hists.py -t ak8_t__tau` to iterate over all years/channels
- `python step2_fit.py -y UL18 -c muo -t ak8_t__tau` to fit with Combine and get the QCD transfer normalization factor from QCD sideband to main band
- `root -l -q -b plot_CompareShapes.cxx`
- `python plot_Stacks.py -t ak8_t` (without `__tau` now)
<!-- - `python latex_slides.py` (output goes here in this directory) -->

Other taggers you must run: `hotvr_t__tau`, `ak8_w__partnet` (the `__tau`/`__partnet` actually does not matter; we will get transfer factors for `ak8_t`, `ak8_w`, `hotvr_t` in the end; it is just important to trim down the T&P selection to having an AK8/HOTVR probejet with the appropriate pT cut)

All output goes to `.../src/UHH2/LegacyTopTagging/output/TagAndProbe/mainsel/UL18/muo/qcd_normalization_study/` and similar for other years/channels
