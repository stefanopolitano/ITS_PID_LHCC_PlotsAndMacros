'''
Script for the plot of big cluster origini
'''
from os import system
import yaml
import sys
from ROOT import TFile, TLatex, TCanvas, gStyle, TGraphAsymmErrors, TLegend, TH1F, TF1, kGray # pylint: disable=import-error,no-name-in-module
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

histoName = ['hCLsize',
             'hCLsizeSel_procPrimary',
             'hCLsizeSel_procd-rays',
             'hCLsizeSel_procHInhelastic'] # keep '' as first element

colors = {'hCLsize': kGray+1,
          'hCLsizeSel_procPrimary': kHInelastic,
          'hCLsizeSel_procd-rays': kDrays,
          'hCLsizeSel_procHInhelastic': kPrimary}

labels = {'hCLsize': 'Total',
          'hCLsizeSel_procPrimary': 'Primary',
          'hCLsizeSel_procd-rays': 'd-rays',
          'hCLsizeSel_procHInhelastic': 'HInelastic'}

inFile = TFile.Open(inputCfg['big_cluster_origin_mc']['input'])
axisLimits = {'pad1': {'x': [0.01, 100], 'y': [1.e-1, 1.e+9]},
              'pad2': {'x': [0.01, 100], 'y': [0, 1.1]},
              'pad3': {'x': [0, 1.], 'y': [1, 1.e+6]}}
axisTitles = {'pad1': ';Cluster Size; Counts',
              'pad2': ';Cluster Size; Ratio to Total',
              'pad3': ';Energy (MeV); Counts'}

histos, histos_ratio, hFrame = ([] for i in range(3))
legend = TLegend(0.5, 0.4, 0.9, 0.7)
SetLegendStyle(legend,
               header='Big Cluster Origin',
               ncolumns=2,
               textsize=900/25545)
for _, histo_name in enumerate(histoName):
    histos.append(inFile.Get(f'{histo_name}'))
    SetObjectStyle(histos[-1], markercolor=colors[histo_name],
                   fillcolor=colors[histo_name],
                   fillstyle=1001,
                   fillalpha=0.3,
                   linewidth=1,
                   linecolor=colors[histo_name],
                   markerstyle=20, markersize=0.5)
    legend.AddEntry(histos[-1], labels[histo_name], 'fpl')
    histos_ratio.append(histos[-1].Clone(f'{histo_name}_ratio'))
    histos_ratio[-1].Divide(histos[0])

cBeta = TCanvas('cBeta', '', 1200, 400)
pads = 3 if inputCfg['big_cluster_origin_mc']['drays_energy'] else 2
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
    if pad == 1:
        for ihisto, histo in enumerate(histos):
            if ihisto == 0:
                histo.Draw('histesame')
            else:
                histo.Draw('histesame')
        LatLabel('This Analyses', 0.25, 0.85, 0.05)
        LatLabel('Run 3 MC', 0.25, 0.8, 0.05)
        LatLabel('pp, #sqrt{#it{s}} = 900 GeV', 0.25, 0.75, 0.05)
        legend.Draw()
    if pad == 2:
        for iratio, historatio in enumerate(histos_ratio):
            if iratio == 0:
                continue
            historatio.Draw('histsame')
        line = LineAtOne(histos_ratio[0].GetXaxis().GetXmin(),
                         histos_ratio[0].GetXaxis().GetXmax())
        line.Draw('same')
    if pad == 3 and inputCfg['big_cluster_origin_mc']['drays_energy']:
        print('Plotting d-rays energy')
        energy_file = TFile.Open(inputCfg['big_cluster_origin_mc']['drays_input'])
        hEkin = energy_file.Get('hEkin')
        hEtotal = energy_file.Get('hEtotal')
        SetObjectStyle(hEkin, markercolor=kAzureCool,
                       linecolor=kAzureCool,
                       fillcolor=kAzureCool,
                       fillalpha=0.3,
                       linewidth=1)
        SetObjectStyle(hEtotal, markercolor=kDplusPrompt,
                       linecolor=kDplusPrompt,
                       fillcolor=kDplusPrompt,
                       fillalpha=0.3,
                       linewidth=1)
        hEkin.Draw('histsame')
        hEtotal.Draw('histsame')
        legDelta = TLegend(0.2, 0.7, 0.6, 0.9)
        legDelta.SetHeader('d-rays Energy')
        legDelta.AddEntry(hEkin, 'Kinetic Energy', 'fpl')
        legDelta.AddEntry(hEtotal, 'Total Energy', 'fpl')
        legDelta.Draw()
    cBeta.Update()

for outformat in inputCfg["outformat"]:
    cBeta.SaveAs(f'{inputCfg["outdir"]}/BigClusterOriginMC_{pads}pads_011022.{outformat}')

input('Press enter to continue')
