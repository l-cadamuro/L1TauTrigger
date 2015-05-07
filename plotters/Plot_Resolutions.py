from ROOT import *
from LucaStyle import *

LS = SetLucaStyle(0)
gROOT.ForceStyle(1)

fileStage2 = TFile.Open ("../Common/resolutions_histo_gg_Stage2_WithHPSIso_eta2p1_ptMin20.00.root")
fileRunI = TFile.Open ("../Common/resolutions_histo_gg_RunI_ptMin20.00.root")

# speed up with structure 0: STage2 -- 1: Run1
fileList = [fileStage2, fileRunI]

pt_EB = []
pt_EE = []
pt_All = []
phi = []
eta = []
grResolEta = []
grResolPt = []

for i in [0,1]:
    pt_EB.append(fileList[i].Get("pt_resol_EB"))
    pt_EE.append(fileList[i].Get("pt_resol_EE"))
    pt_All.append(fileList[i].Get("pt_resol_All"))
    phi.append(fileList[i].Get("phi_resol_All"))
    eta.append(fileList[i].Get("eta_resol_All"))
    grResolEta.append(fileList[i].Get("resolution_vs_eta_plot"))
    grResolPt.append(fileList[i].Get("resolution_vs_pT_plot"))
# set all npx to 1000
for i in [0,1]:
    pt_EB[i].GetFunction("CBFuncAsymm").SetNpx(1000) #non eccede mai la risoluzione dei 100 px... ?
    pt_EE[i].GetFunction("CBFuncAsymm").SetNpx(1000) # ma se scendo sotto (e.g. 10) si vede la differenza!
    pt_All[i].GetFunction("CBFuncAsymm").SetNpx(1000)

eta[0].GetFunction("CBFunc").SetNpx(1000)
phi[0].GetFunction("CBFunc").SetNpx(1000)


###########################

canv1 = TCanvas ("canv1", "canv1")

headLeg = TLegend (0.05, 0.926, 0.887, 0.996)
headLeg.SetFillStyle(0)
headLeg.SetTextSize(0.030)
#headLeg.SetTextFont(62)
headLeg.AddEntry("NULL","CMS Simulation 2015: gg #rightarrow H #rightarrow #tau #tau - #sqrt{s}=13 TeV, bx=25ns, <PU>=40","h")

#############
# barrel / endcap

SetColors(pt_EB[0], kRed)
SetColors(pt_EE[0], kBlue)
pt_EE[0].GetFunction("CBFuncAsymm").SetLineColor(kAzure+5)
pt_EB[0].GetFunction("CBFuncAsymm").SetLineColor(kOrange+8)

leg = TLegend (0.50, 0.66, 0.88, 0.797)
leg.SetFillStyle(0)
#leg.SetTextFont(72)
leg.AddEntry(pt_EB[0], "barrel: #eta < 1.305", "lp")
leg.AddEntry(pt_EE[0], "endcap: |#eta| > 1.479", "lp")

pt_EB[0].Draw()
pt_EE[0].Draw("same")
headLeg.Draw()
leg.Draw()

canv1.Print ("resolution_barr_endc.pdf")


############
# eta Stage2 / run I

canv2 = TCanvas ("canv2", "canv2")

legEta = TLegend (0.18, 0.5857, 0.56, 0.7255)
legEta.SetFillStyle(0)
legEta.SetTextSize(0.030)
legEta.AddEntry (eta[0], "Upgrade 2016", "lp")
legEta.AddEntry (eta[1], "Run I", "lp")

legSigmaEta = TLegend (0.58, 0.60, 0.95, 0.78)
legSigmaEta.SetFillStyle(0)
legSigmaEta.AddEntry ("NULL", "#sigma = (3.21 #pm 0.01) #upoint 10^{-2}", "h")

SetColors(eta[0], kBlue)
SetColors(eta[1], kRed)
eta[0].GetFunction("CBFunc").SetLineColor(kAzure+5)

eta[0].SetMaximum(0.08)
eta[0].Draw()
eta[1].Draw("same")
headLeg.Draw()
legEta.Draw()
legSigmaEta.Draw()
canv2.Print("eta_resol_vs_RunI.pdf", "pdf")


############
# phi Stage2 / run I


canv3 = TCanvas ("canv3", "canv3")
legPhi = TLegend (00.18, 0.5857, 0.56, 0.7255)
legPhi.SetFillStyle(0)
legPhi.SetTextSize(0.030)
legPhi.AddEntry (phi[0], "Upgrade 2016", "lp")
legPhi.AddEntry (phi[1], "Run I", "lp")

legSigmaPhi = TLegend (0.58, 0.60, 0.95, 0.78)
legSigmaPhi.SetFillStyle(0)
legSigmaPhi.AddEntry ("NULL", "#sigma = (3.47 #pm 0.01) #upoint 10^{-2}", "h")

SetColors(phi[0], kBlue)
SetColors(phi[1], kRed)
phi[0].GetFunction("CBFunc").SetLineColor(kAzure+5)

phi[0].SetMaximum(0.08)
phi[0].Draw()
phi[1].Draw("same")
headLeg.Draw()
legPhi.Draw()
legSigmaPhi.Draw()
canv3.Print("phi_resol_vs_RunI.pdf", "pdf")

###########



############
# compare pt resolution (inclusive) for Run I and now
canv4 = TCanvas ("canv4", "canv4")

SetColors(pt_All[0], kBlue)
SetColors(pt_All[1], kRed)
pt_All[0].GetFunction("CBFuncAsymm").SetLineColor(kAzure+5)
pt_All[1].GetFunction("CBFuncAsymm").SetLineColor(kOrange+8)
legRun = TLegend (0.50, 0.66, 0.88, 0.797)
legRun.SetFillStyle(0)
#leg.SetTextFont(72)
legRun.AddEntry(pt_All[0], "Upgrade 2016", "lp")
legRun.AddEntry(pt_All[1], "Run I (rescaled)", "lp")

pt_All[0].SetMaximum(0.08)
pt_All[0].Draw()
pt_All[1].Draw("same")
#pt_EB[0].SetMaximum(0.08) # tried as a test
#pt_EB[0].Draw()
#pt_EB[1].Draw("same")

headLeg.Draw()
legRun.Draw()

canv4.Print ("resolution_Stage2_RunI.pdf")

############
# compare resolution plots - eta
canv5 = TCanvas ("canv5", "canv5")
SetColors(grResolEta[0], kBlue)
SetColors(grResolEta[1], kRed)

legRunResEta = TLegend (0.19, 0.18, 0.57, 0.31)
#legRunResEta.SetFillStyle(0)
legRunResEta.SetFillColor(kWhite)
legRunResEta.SetBorderSize(1)
#legRunResEta.SetLineColor(kBlack)
#leg.SetTextFont(72)
legRunResEta.AddEntry(grResolEta[0], "Upgrade 2016", "lp")
legRunResEta.AddEntry(grResolEta[1], "Run I", "lp")

grResolEta[0].GetYaxis().SetTitle ("Energy resolution")
grResolEta[0].SetMinimum(0.0)
grResolEta[0].SetMaximum(0.30)
grResolEta[0].Draw()
grResolEta[1].Draw("same")
headLeg.Draw()
legRunResEta.Draw()
canv5.Print("resolution_vs_eta.pdf", "pdf")


###############
# compare resolution plots - pt
canv6 = TCanvas ("canv6", "canv6")
SetColors(grResolPt[0], kBlue)
SetColors(grResolPt[1], kRed)

legRunResPt = TLegend (0.19, 0.18, 0.57, 0.31)
#legRunResEta.SetFillStyle(0)
legRunResPt.SetFillColor(kWhite)
legRunResPt.SetBorderSize(1)
#legRunResEta.SetLineColor(kBlack)
#leg.SetTextFont(72)
legRunResPt.AddEntry(grResolPt[0], "Upgrade 2016", "lp")
legRunResPt.AddEntry(grResolPt[1], "Run I", "lp")

grResolPt[0].GetYaxis().SetTitle ("Energy resolution")
grResolPt[0].SetMinimum(0.)
grResolPt[0].SetMaximum(0.35)
grResolPt[0].GetXaxis().SetRangeUser(20., 150.)
grResolPt[0].Draw()
grResolPt[1].Draw("same")
headLeg.Draw()
legRunResPt.Draw()
canv6.Print("resolution_vs_pt.pdf", "pdf")

raw_input()

