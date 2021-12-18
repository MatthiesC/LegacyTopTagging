import ROOT as root
import os
from copy import deepcopy

basepath = os.path.join(os.environ.get('CMSSW_BASE'), 'src/UHH2/LegacyTopTagging/output/TT_Mtt_study')

file_TT_inclusive = root.TFile.Open(os.path.join(basepath, 'uhh2.AnalysisModuleRunner.MC.TT_inclusive_UL16postVFP.root'), 'READ')
file_TT_Mtt700to1000 = root.TFile.Open(os.path.join(basepath, 'uhh2.AnalysisModuleRunner.MC.TT_Mtt700to1000_UL16postVFP.root'), 'READ')
file_TT_Mtt1000toInf = root.TFile.Open(os.path.join(basepath, 'uhh2.AnalysisModuleRunner.MC.TT_Mtt1000toInf_UL16postVFP.root'), 'READ')

hist_TT_inclusive = file_TT_inclusive.Get('TT_Mtt_Hists/mtt')
hist_TT_Mtt700to1000 = file_TT_Mtt700to1000.Get('TT_Mtt_Hists/mtt')
hist_TT_Mtt1000toInf = file_TT_Mtt1000toInf.Get('TT_Mtt_Hists/mtt')

int_TT_inclusive = hist_TT_inclusive.Integral(0, 2001)
int_TT_inclusive_Mtt700to1000 = hist_TT_inclusive.Integral(701, 1000)
int_TT_inclusive_Mtt1000toInf = hist_TT_inclusive.Integral(1001, 2001)
int_TT_Mtt700to1000 = hist_TT_Mtt700to1000.Integral(0, 2001)
int_TT_Mtt1000toInf = hist_TT_Mtt1000toInf.Integral(0, 2001)

# hist_TT_inclusive.Rebin(10)
# hist_TT_Mtt700to1000.Rebin(10)
# hist_TT_Mtt1000toInf.Rebin(10)

print str(int_TT_inclusive)
print str(int_TT_inclusive_Mtt700to1000)
print str(int_TT_inclusive_Mtt1000toInf)
print str(int_TT_Mtt700to1000)
print str(int_TT_Mtt1000toInf)

corr_Mtt700to1000 = int_TT_inclusive_Mtt700to1000/int_TT_Mtt700to1000
corr_Mtt1000toInf = int_TT_inclusive_Mtt1000toInf/int_TT_Mtt1000toInf

print 'Correction factor Mtt700to1000:', corr_Mtt700to1000
print 'Correction factor Mtt1000toInf:', corr_Mtt1000toInf

c = root.TCanvas('canvas', '', 800, 800)
root.gStyle.SetOptStat(0)
c.cd()
p = root.TPad('pad', '', 0, 0, 1, 1)
p.SetMargin(0.15, 0.05, 0.15, 0.05)
p.Draw()
p.cd()

for hist in [hist_TT_inclusive, hist_TT_Mtt700to1000, hist_TT_Mtt1000toInf]:
    hist.SetMarkerStyle(0)
    hist.SetMarkerSize(0)

leg = root.TLegend(0.25, 0.5, 0.92, 0.92)

draw_option = 'hist l'

# hist_TT_inclusive.SetTitle('inclusive')
hist_TT_inclusive.SetTitle('')
hist_TT_inclusive.Draw(draw_option)
leg.AddEntry(hist_TT_inclusive, 'inclusive (had+semilep+dilep), 831.76 pb (NNLO+NNLL)')
hist_TT_inclusive.SetLineColor(1)

hist_TT_Mtt700to1000.SetLineColor(2)
# hist_TT_Mtt700to1000.SetTitle('700to1000')
hist_TT_Mtt700to1000.Draw('same '+draw_option)
leg.AddEntry(hist_TT_Mtt700to1000, 'Mtt-700to1000, 64.72 pb (NLO from GenXSecAnalyzer)')

hist_TT_Mtt1000toInf.SetLineColor(4)
# hist_TT_Mtt1000toInf.SetTitle('1000toInf')
hist_TT_Mtt1000toInf.Draw('same '+draw_option)
leg.AddEntry(hist_TT_Mtt1000toInf, 'Mtt-1000toInf, 16.44 pb (NLO from GenXSecAnalyzer)')

hist_TT_Mtt700to1000_corr = deepcopy(hist_TT_Mtt700to1000)
hist_TT_Mtt700to1000_corr.Scale(corr_Mtt700to1000)
hist_TT_Mtt700to1000_corr.SetLineStyle(2)
hist_TT_Mtt700to1000_corr.SetLineColor(3)
# hist_TT_Mtt700to1000_corr.SetTitle('700to1000 (corrected)')
hist_TT_Mtt700to1000_corr.Draw('same '+draw_option)
leg.AddEntry(hist_TT_Mtt700to1000_corr, 'Mtt-700to1000 (x%.2f)' % corr_Mtt700to1000)

hist_TT_Mtt1000toInf_corr = deepcopy(hist_TT_Mtt1000toInf)
hist_TT_Mtt1000toInf_corr.Scale(corr_Mtt1000toInf)
hist_TT_Mtt1000toInf_corr.SetLineStyle(2)
hist_TT_Mtt1000toInf_corr.SetLineColor(7)
# hist_TT_Mtt1000toInf_corr.SetTitle('1000toInf (corrected)')
hist_TT_Mtt1000toInf_corr.Draw('same '+draw_option)
leg.AddEntry(hist_TT_Mtt1000toInf_corr, 'Mtt-1000toInf (x%.2f)' % corr_Mtt1000toInf)



hist_TT_inclusive.GetXaxis().SetRangeUser(600, 1600)
hist_TT_inclusive.GetXaxis().SetTitle('#it{m}_{t#bar{t}} [GeV]')
hist_TT_inclusive.GetXaxis().SetTitleOffset(1.5)
# hist_TT_inclusive.GetYaxis().SetTitle('Expected events for #it{L}_{int} = 10000 pb^{#minus1}')
hist_TT_inclusive.GetYaxis().SetTitle('Expected events [arb. units]')

leg.Draw()

c.SaveAs('mtt_plot.pdf')
