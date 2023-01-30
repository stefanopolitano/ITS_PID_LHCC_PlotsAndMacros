'''
Script for the plot of big cluster origini
'''
from os import system
import yaml
import sys
from ROOT import TFile, TLatex, TCanvas, gStyle, TGraphAsymmErrors, TLegend, TH1F, TF1, kGray # pylint: disable=import-error,no-name-in-module
sys.path.append('..')
from utils.AnalysisUtils import fit_reso
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle, SetXsystForLogScale, LatLabel, kDrays, kHInelastic, kPrimary, LineAtOne, SetLegendStyle, kAzureCool, kDplusPrompt, kRed

SetGlobalStyle(padbottommargin=0.14, padleftmargin=0.15,
               padrightmargin=0.1, padtopmargin=0.05,
               titleoffsety=1.4, maxdigits=3,
               titlesizex=1200/25550,
               titlesizey=1200/25550,
               labelsizex=900/25545,
               labelsizey=900/25545)

with open('input_config.yml', 'r') as ymlCfgFile:
    inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)

inFile = TFile.Open(inputCfg['clsizecoslam'])

hP = inFile.Get('0')
hP.SetName('hP')
hK = inFile.Get('1')
hK.SetName('hK')
hPi = inFile.Get('2')
hPi.SetName('hPi')
SetObjectStyle(hPi, markercolor=kAzureCool,
               fillcolor=kAzureCool, fillstyle=1001,
               fillalpha=0.3, linewidth=1,
               linecolor=kAzureCool, markerstyle=20,
               markersize=1)
SetObjectStyle(hK, markercolor=kDplusPrompt,
               fillcolor=kDplusPrompt,
               fillstyle=1001, fillalpha=0.3,
               linewidth=1, linecolor=kDplusPrompt,
               markerstyle=20, markersize=1)
SetObjectStyle(hP, markercolor=kRed+1,
               fillcolor=kRed+1, fillstyle=1001,
               fillalpha=0.3, linewidth=1,
               linecolor=kRed+1, markerstyle=20,
               markersize=1)

legend = TLegend(0.63, 0.4, 0.9, 0.7)
SetLegendStyle(legend,
               ncolumns=1,
               textsize=900/25545)
legend.AddEntry(hPi, '#pi', 'flp')
legend.AddEntry(hK, 'K', 'flp')
legend.AddEntry(hP, 'p', 'flp')

# resolution fit (TO DO: move to flarefly)
outfile_reso = TFile('../plots/Resolution.root', 'recreate')
reso1, reso_unc1, gaus1 = fit_reso(hP, nsigma=1)
gaus1.SetName('fPreso1')
reso2, reso_unc2, gaus2 = fit_reso(hP, nsigma=2)
gaus2.SetName('fPreso2')
reso3, reso_unc3, gaus3 = fit_reso(hP, nsigma=3)
gaus3.SetName('fPreso3')
print('\033\[1mResolution Proton\033\[0m')
print(f'Resolution: {reso1} +- {reso_unc1}')
print(f'Resolution: {reso2} +- {reso_unc2}')
print(f'Resolution: {reso3} +- {reso_unc3}')
print(f'RMS: {hP.GetRMS()/hP.GetMean()} +- {hP.GetRMSError()/hP.GetMean()}')
canv = TCanvas('canvP', 'canv', 1200, 1200)
hFrame = canv.cd().DrawFrame(hP.GetMean()-3.5*hP.GetRMS(),
                             0.0001,
                             hP.GetMean()+3.5*hP.GetRMS(),
                             500,
                             "; #LT cluster size #GT #times cos#lambda; counts")
hP.Draw('hist same')
LatLabel('Proton resolution', 0.2, 0.85, 0.05)
leghp = TLegend(0.55, 0.6, 0.9, 0.85)
gaus1.Draw('same')
gaus2.Draw('same')
gaus3.Draw('same')
leghp.AddEntry(gaus1, 'reso (1#sigma) = %.3f #pm %.3f' % (reso1, reso_unc1), 'l')
leghp.AddEntry(gaus2, 'reso (2#sigma) = %.3f #pm %.3f' % (reso2, reso_unc2), 'l')
leghp.AddEntry(gaus3, 'reso (3#sigma) = %.3f #pm %.3f' % (reso3, reso_unc3), 'l')
leghp.Draw('same')
canv.SaveAs('../plots/ResolutionProton.pdf')
canv.Write()
hP.Write()
gaus1.Write()
gaus2.Write()
gaus3.Write()
input('Press enter to continue')

reso1, reso_unc1, gaus1 = fit_reso(hK, nsigma=1)
gaus1.SetName('fKreso1')
reso2, reso_unc2, gaus2 = fit_reso(hK, nsigma=2)
gaus2.SetName('fKreso2')
reso3, reso_unc3, gaus3 = fit_reso(hK, nsigma=3)
gaus3.SetName('fKreso3')
print('\033\[1mResolution Kaon\033\[0m')
print(f'Resolution: {reso1} +- {reso_unc1}')
print(f'Resolution: {reso2} +- {reso_unc2}')
print(f'Resolution: {reso3} +- {reso_unc3}')
print(f'RMS: {hK.GetRMS()/hK.GetMean()} +- {hK.GetRMSError()/hK.GetMean()}')
canv = TCanvas('canvK', 'canv', 1200, 1200)
ymax = hK.GetMaximum()*1.2
hFrame = canv.cd().DrawFrame(hK.GetMean()-3.5*hK.GetRMS(),
                             0.0001,
                             hK.GetMean()+3.5*hK.GetRMS(),
                             ymax,
                             "; #LT cluster size #GT #times cos#lambda; counts")
hK.Draw('hist same')
LatLabel('Kaon resolution', 0.2, 0.85, 0.05)
leghk = TLegend(0.6, 0.6, 0.9, 0.85)
gaus1.Draw('same')
gaus2.Draw('same')
gaus3.Draw('same')
leghk.AddEntry(gaus1, 'reso (1#sigma) = %.3f #pm %.3f' % (reso1, reso_unc1), 'l')
leghk.AddEntry(gaus2, 'reso (2#sigma) = %.3f #pm %.3f' % (reso2, reso_unc2), 'l')
leghk.AddEntry(gaus3, 'reso (3#sigma) = %.3f #pm %.3f' % (reso3, reso_unc3), 'l')
leghk.Draw('same')
canv.SaveAs('../plots/ResolutionKaon.pdf')
canv.Write()
hK.Write()
gaus1.Write()
gaus2.Write()
gaus3.Write()
input('Press enter to continue')

reso1, reso_unc1, gaus1 = fit_reso(hPi, nsigma=1)
gaus1.SetName('hPireso1')
reso2, reso_unc2, gaus2 = fit_reso(hPi, nsigma=2)
gaus2.SetName('hPireso2')
reso3, reso_unc3, gaus3 = fit_reso(hPi, nsigma=3)
gaus3.SetName('hPireso3')
print('\033Resolution Pion\033\[0m')
print(f'Resolution: {reso1} +- {reso_unc1}')
print(f'Resolution: {reso2} +- {reso_unc2}')
print(f'Resolution: {reso3} +- {reso_unc3}')
print(f'RMS: {hPi.GetRMS()/hPi.GetMean()} +- {hPi.GetRMSError()/hPi.GetMean()}')
input('Press enter to continue')
canv = TCanvas('canvPi', 'canv', 1200, 1200)
ymax = hPi.GetMaximum()*1.2
hFrame = canv.cd().DrawFrame(hPi.GetMean()-3.5*hPi.GetRMS(),
                             0.0001,
                             hPi.GetMean()+3.5*hPi.GetRMS(),
                             ymax,
                             "; #LT cluster size #GT #times cos#lambda; counts")
hPi.Draw('hist same')
LatLabel('Pion resolution', 0.2, 0.85, 0.05)
leghpi = TLegend(0.6, 0.6, 0.9, 0.85)
gaus1.Draw('same')
gaus2.Draw('same')
gaus3.Draw('same')
leghpi.AddEntry(gaus1, 'reso (1#sigma) = %.3f #pm %.3f' % (reso1, reso_unc1), 'l')
leghpi.AddEntry(gaus2, 'reso (2#sigma) = %.3f #pm %.3f' % (reso2, reso_unc2), 'l')
leghpi.AddEntry(gaus3, 'reso (3#sigma) = %.3f #pm %.3f' % (reso3, reso_unc3), 'l')
leghpi.Draw('same')
canv.SaveAs('../plots/ResolutionPi.pdf')
canv.Write()
hPi.Write()
gaus1.Write()
gaus2.Write()
gaus3.Write()
outfile_reso.Close()

# cluster size vs cos lambda
cCLSize = TCanvas('cCLSize', '', 800, 800)
#cCLSize.cd().SetLogy()
hFrame = cCLSize.cd().DrawFrame(0.,
                                0.0001,
                                12.,
                                0.1,
                                "; #LT cluster size #GT #times cos#lambda; norm. counts")
hFrame.GetYaxis().SetDecimals()
LatLabel('This Analysis', 0.58, 0.88, 0.05)
LatLabel('Run 3', 0.58, 0.83, 0.05)
LatLabel('pp, #sqrt{#it{s}} = 13 TeV', 0.58, 0.78, 0.03)
LatLabel('0.3 < #it{p}^{ITS-TPC} < 0.4 GeV/#it{c}', 0.55, 0.73, 0.03)
hPi.DrawNormalized('same')
hK.DrawNormalized('same')
hP.DrawNormalized('same')
legend.Draw()
cCLSize.Update()

for outformat in inputCfg["outformat"]:
    cCLSize.SaveAs(f'{inputCfg["outdir"]}/ClusterSizeCosLam.{outformat}')

input('Press enter to continue')
