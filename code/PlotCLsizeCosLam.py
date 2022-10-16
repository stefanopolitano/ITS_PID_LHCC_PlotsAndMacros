'''
Script for the plot of big cluster origini
'''
from os import system
import yaml
import sys
from ROOT import TFile, TLatex, TCanvas, gStyle, TGraphAsymmErrors, TLegend, TH1F, TF1, kGray # pylint: disable=import-error,no-name-in-module
sys.path.append('..')
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
hK = inFile.Get('1')
hPi = inFile.Get('2')
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

cCLSize = TCanvas('cCLSize', '', 800, 800)
#cCLSize.cd().SetLogy()
hFrame = cCLSize.cd().DrawFrame(0.,
                                0.0001,
                                12.,
                                0.1,
                                "; cluster size #times #LT cos#lambda #GT; norm. counts")
hFrame.GetYaxis().SetDecimals()
LatLabel('This Analyses', 0.58, 0.88, 0.05)
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
