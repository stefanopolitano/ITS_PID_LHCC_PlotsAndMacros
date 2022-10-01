'''
Script for the plot of MC beta
'''
from os import system
import yaml
import sys
from ROOT import TFile, TLatex, TCanvas, gStyle, TGraphAsymmErrors, TLegend, TH1F, TDirectoryFile, TMath, TF1 # pylint: disable=import-error,no-name-in-module
sys.path.append('..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle, SetXsystForLogScale, LatLabel

SetGlobalStyle(padbottommargin=0.14, padleftmargin=0.15, padrightmargin=0.15, padtopmargin=0.05, titleoffsety=1.4, maxdigits=2, palette=53)

gStyle.SetPalette(53)

with open('input_config.yml', 'r') as ymlCfgFile:
    inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)

inFile = TFile.Open(inputCfg['beta_MC'])
hbeta_vs_mc = inFile.Get('scatter_plot')

cBeta = TCanvas('cBeta', '', 900, 800)
cBeta.cd().SetLogz()
hFrame = cBeta.cd().DrawFrame(0, 0., 1., 1.1, ';#it{p}_{T} (GeV/#it{c}); #beta (ML)')
hFrame.GetYaxis().SetDecimals()
hFrame.GetXaxis().SetDecimals()
hbeta_vs_mc.Draw('samecolz')
#latALICE = LatLabel('ALICE', 0.6, 0.24, 0.04) # wait for approval
latSystem = LatLabel('Run 3 MC ', 0.55, 0.24, 0.04)
latSystem = LatLabel('pp, #sqrt{#it{s}} = 13.6 TeV', 0.55, 0.18, 0.04)
cBeta.Update()
for outformat in inputCfg["outformat"]:
    cBeta.SaveAs(f'{inputCfg["outdir"]}/beta_vs_p_MC_011022.{outformat}')

input('Press enter to exit')
