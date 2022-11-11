'''
Script for the plot of Ln big clusters
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

layer = 'L0'
if inputCfg['bigcluster']['layer']: layer = inputCfg['bigcluster']['layer']
else : print('No layer specified, using L0')


inFileData = TFile.Open(inputCfg['bigcluster']['data'])
hData = inFileData.Get(f'hClusterSize{layer}Oct')
normData = hData.Integral(1, hData.GetNbinsX()) # scaling suggest by Felix to avoid the first bin (which depends on masking)
SetObjectStyle(hData, markercolor=kAzureCool,
               fillcolor=kAzureCool, fillstyle=1001,
               fillalpha=0.3, linewidth=1,
               linecolor=kAzureCool, markerstyle=20,
               markersize=1)
hData.Scale(1./normData)
if layer != 'L0' and layer != 'L6':
    input('L0 or L6 only in MC, please check input_config.yml')
    sys.exit(1)
infileMC = TFile.Open(inputCfg['bigcluster']['mc'])
hMC = infileMC.Get(f'hC{layer}')
hMCshift = TH1F(f'hMCshift{layer}', '', hMC.GetNbinsX(), hMC.GetXaxis().GetXmin()-0.5, hMC.GetXaxis().GetXmax()-0.5)
for ibin in range(1, hMC.GetNbinsX()+1):
    hMCshift.SetBinContent(ibin, hMC.GetBinContent(ibin))
hMC = hMCshift
normMC = hMC.Integral(1, hMC.GetNbinsX())
SetObjectStyle(hMC, markercolor=kRed+1,
               fillcolor=kRed+1, fillstyle=1001,
               fillalpha=0.3, linewidth=1,
               linecolor=kRed+1, markerstyle=20,
               markersize=1)
hMC.Scale(1./normMC)

legend = TLegend(0.51, 0.6, 0.9, 0.75)
SetLegendStyle(legend,
               ncolumns=1,
               textsize=800/25545,
               header='Pilot beam data')
legend.AddEntry(hData, 'Data - run (505658)', 'flp')
legend.AddEntry(hMC, 'Anchored MC', 'flp')

cCLSize = TCanvas('cCLSize', '', 800, 800)
cCLSize.cd().SetLogy()
hFrame = cCLSize.cd().DrawFrame(0.,
                                1.e-8,
                                100.,
                                10,
                                f";cluster size on {layer}; norm. counts")
hFrame.GetYaxis().SetDecimals()
LatLabel('This Analysis', 0.52, 0.88, 0.05)
LatLabel('Run 3', 0.52, 0.83, 0.05)
LatLabel('pp, #sqrt{#it{s}} = 900 GeV', 0.52, 0.78, 0.03)
hData.DrawCopy('samehiste')
hMC.DrawCopy('samehiste')
legend.Draw()
cCLSize.Update()

for outformat in inputCfg["outformat"]:
    cCLSize.SaveAs(f'{inputCfg["outdir"]}/BigClustersDatavsMC_{layer}.{outformat}')

input('Press enter to continue')
