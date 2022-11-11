'''
Script for the plot of big cluster origini
'''
from os import system
import yaml
import sys
from ROOT import TFile, TLatex, TCanvas, gStyle, TGraphAsymmErrors, TLegend, TH1F, TF1, kGray, kRed, kAzure, kGreen # pylint: disable=import-error,no-name-in-module
sys.path.append('..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle, SetXsystForLogScale, LatLabel, kDrays, kHInelastic, kPrimary, LineAtOne, SetLegendStyle, kAzureCool, kDplusPrompt

SetGlobalStyle(padbottommargin=0.14, padleftmargin=0.15,
               padrightmargin=0.1, padtopmargin=0.05,
               titleoffsety=1.4, maxdigits=3,
               titlesizex=1200/25550,
               titlesizey=1200/25550,
               labelsizex=900/25545,
               labelsizey=900/25545)

with open('input_config.yml', 'r') as ymlCfgFile:
    inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)

histoName = ['h_proton_clus_slice',
             'h_proton_clus_slice_0',
             'h_proton_clus_slice_1',
             'h_proton_clus_slice_2']

colors = {'h_proton_clus_slice': kGray+1,
          'h_proton_clus_slice_0': kRed+1,
          'h_proton_clus_slice_1': kAzure+4,
          'h_proton_clus_slice_2': kGreen+2}

labels = {'h_proton_clus_slice': 'Total',
          'h_proton_clus_slice_0': '1^{st} BC',
          'h_proton_clus_slice_1': '2^{nd} BC',
          'h_proton_clus_slice_2': '3^{rd} BC'}

inFile = TFile.Open(inputCfg['cluster_splitting'])
axisLimits = {'pad1': {'x': [0.01, 2000], 'y': [1.e-1, 1.e+9]},
              'pad2': {'x': [0.01, 20], 'y': [0, 1000]}}
axisTitles = {'pad1': ';BC; Counts',
              'pad2': ';Cluster Size; Counts'}

histos, histos_BC, hFrame = ([] for i in range(3))
legend = TLegend(0.5, 0.4, 0.9, 0.7)
SetLegendStyle(legend,
               header='BC Origin',
               ncolumns=2,
               textsize=900/25545)
for _, histo_name in enumerate(histoName):
    print(histo_name)
    histos.append(inFile.Get(f'{histo_name}'))
    SetObjectStyle(histos[-1], markercolor=colors[histo_name],
                   fillcolor=colors[histo_name],
                   fillstyle=1001,
                   fillalpha=0.3,
                   linewidth=1,
                   linecolor=colors[histo_name],
                   markerstyle=20, markersize=0.5)
    legend.AddEntry(histos[-1], labels[histo_name], 'fpl')

pads = 2
cBeta = TCanvas('cBeta', '', 800, 400)
cBeta.Divide(pads, 1)

for pad in range(1, pads+1):
    cBeta.cd(pad)
    if pad != 2:
        cBeta.cd(pad).SetLogy()
    hFrame.append(cBeta.cd(pad).DrawFrame(axisLimits[f'pad{pad}']['x'][0],
                                          axisLimits[f'pad{pad}']['y'][0],
                                          axisLimits[f'pad{pad}']['x'][1],
                                          axisLimits[f'pad{pad}']['y'][1],
                                          axisTitles[f'pad{pad}']))
    hFrame[-1].GetYaxis().SetDecimals()
    if pad == 2:
        for ihisto, histo in enumerate(histos):
            histo.Draw('histesame')
        LatLabel('This Analysis', 0.25, 0.85, 0.05)
        LatLabel('Run 3', 0.25, 0.8, 0.05)
        LatLabel('pp, #sqrt{#it{s}} = 13.6 GeV', 0.25, 0.75, 0.05)
        legend.Draw()
    if pad == 1:
        for iBC in range(3):
            histos_BC.append(inFile.Get(f'h_rof_bc{iBC}'))
            SetObjectStyle(histos_BC[iBC], markercolor=colors[histoName[iBC]],
                           fillcolor=colors[histoName[iBC]],
                           fillsyle=1001,
                           linewidth=1,
                           linecolor=colors[histoName[iBC]],
                           markerstyle=20, markersize=0.5)
            histos_BC[iBC].Draw('histsame')
    cBeta.Update()

for outformat in inputCfg["outformat"]:
    cBeta.SaveAs(f'{inputCfg["outdir"]}/ClusterSplitting_111122.{outformat}')

input('Press enter to continue')
