'''
Script for the plot of MC beta, Armenteros, ITS-TPC data from first plenary 
'''
from os import system
import yaml
import sys
from ROOT import TFile, TLatex, TCanvas, gStyle, TGraphAsymmErrors, TLegend, TH1F, TDirectoryFile, TMath, TF1 # pylint: disable=import-error,no-name-in-module
sys.path.append('..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle, SetXsystForLogScale, LatLabel

SetGlobalStyle(padbottommargin=0.14, padleftmargin=0.15, padrightmargin=0.15,
               padtopmargin=0.05, titleoffsety=1.4, maxdigits=2, palette=53,
               titlesizex=0.05, titlesizey=0.05, labelsize=0.03)

gStyle.SetPalette(53)

with open('input_config.yml', 'r') as ymlCfgFile:
    inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)

inFiles = ['beta_MC',
           'armenteros',
           'beta_data_plenary'] #'beta_MC', 'armenteros', 'beta_data_plenary'
histoName = {'beta_MC': 'scatter_plot',
             'armenteros': 'h_th2_arm',
             'beta_data_plenary': 'hBgVsPReg'}
axisTitles = {'beta_MC': ';#it{p} (GeV/#it{c}); #beta (ML)',
              'armenteros': ';#alpha^{Arm.}; q_{T}^{Arm.}',
              'beta_data_plenary': ';#it{p} (GeV/#it{c}); #beta (ML)'}
axisLimits = {'beta_MC': [0, 0, 1., 1.1],
              'armenteros': [-1, 0., 1, 0.3],
              'beta_data_plenary': [0., 0, 1., 1.1]}
latexLabels = {'beta_MC': [0.55, 0.3, 0.24, 0.18], # x is fixed for all, y is variable
               'armenteros': [0.55, 0.88, 0.84, 0.8],
               'beta_data_plenary': [0.55, 0.3, 0.24, 0.18]}

for infile in inFiles:
    inFile = TFile.Open(inputCfg[infile])
    hbeta_vs_mc = inFile.Get(f'{histoName[infile]}')

    cBeta = TCanvas('cBeta', '', 900, 800)
    cBeta.cd().SetLogz()
    hFrame = cBeta.cd().DrawFrame(axisLimits[infile][0], axisLimits[infile][1],
                                  axisLimits[infile][2], axisLimits[infile][3],
                                  axisTitles[infile])
    hFrame.GetYaxis().SetDecimals()
    hFrame.GetXaxis().SetDecimals()
    hFrame.GetZaxis().SetDecimals()
    hbeta_vs_mc.Draw('samecolz')
    #latALICE = LatLabel('ALICE', 0.6, 0.24, 0.04) # wait for approval
    latALICE = LatLabel('This Analyses', latexLabels[infile][0],
                        latexLabels[infile][1], 0.05)
    run_label = 'Run 3 MC' if infile == 'beta_MC' else 'Run 3'
    latSystem = LatLabel(f'{run_label}', latexLabels[infile][0],
                             latexLabels[infile][2], 0.04)
    latSystem = LatLabel('pp, #sqrt{#it{s}} = 900 GeV', latexLabels[infile][0],
                         latexLabels[infile][3], 0.04)
    cBeta.Update()
    for outformat in inputCfg["outformat"]:
        cBeta.SaveAs(f'{inputCfg["outdir"]}/{infile}_011022.{outformat}')

    input('Press enter to continue')
