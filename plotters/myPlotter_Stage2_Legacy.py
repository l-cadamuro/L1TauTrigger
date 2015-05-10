import ROOT

ROOT.ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetTickLength(0.02,"X");
ROOT.gStyle.SetTickLength(0.02,"Y");
ROOT.gStyle.SetPadTickY(1);
ROOT.gStyle.SetPadTickX(1);
#/ROOT.gStyle.SetPadTopMargin(0.05);
ROOT.gStyle.SetPadBottomMargin(0.11);
ROOT.gStyle.SetPadLeftMargin(0.13);
ROOT.gStyle.SetTitleYOffset(1.1);
ROOT.gStyle.SetTitleXOffset(0.9);
ROOT.gStyle.SetLabelOffset(0.009, "XYZ");
ROOT.gStyle.SetPadGridX(1)
ROOT.gStyle.SetPadGridY(1)

#ROOT.gStyle.SetTitleSize (0.050)
#ROOT.gStyle.SetTitleSize (0.050)

# DUMMY HISTO con tutti i range etc..
hDummy = ROOT.TH1F("hDummy_", "", 1, 0, 110)
hDummy.SetAxisRange(0, 1.25, "Y")
hDummy.GetXaxis().SetTitle ("p_{T}^{offl} [GeV]")
hDummy.GetYaxis().SetTitle ("Efficiency")
hDummy.GetXaxis().SetTitleSize(0.050)
hDummy.GetYaxis().SetTitleSize(0.050)


# Legacy files
legacyFile = ROOT.TFile.Open ("/home/llr/cms/cadamuro/Level1_Stage2/TurnOnFitter/UMF/FitEfficiency/results/TurnOn_RunI_gg_NoHPSIso_bx25Pu40E13TeV_rescaledL1/Tau_stage1_stage2_EB_EE_All_Iso.root")

legacy_B = legacyFile.Get("histo_Stage1_Barrel_vs_Pt")
legacy_B.__class__ = ROOT.RooHist
legacy_E = legacyFile.Get("histo_Stage1_Endcaps_vs_Pt")
legacy_E.__class__ = ROOT.RooHist
legacy_All = legacyFile.Get("histo_Stage1_All_vs_Pt")
legacy_All.__class__ = ROOT.RooHist
legacy_fitB = legacyFile.Get("fit_Stage1_All_vs_Pt")
legacy_fitB.__class__ = ROOT.RooCurve
legacy_fitE = legacyFile.Get("fit_Stage1_Endcaps_vs_Pt")
legacy_fitE.__class__ = ROOT.RooCurve
legacy_fitAll = legacyFile.Get("fit_Stage1_All_vs_Pt")
legacy_fitAll.__class__ = ROOT.RooCurve
legacy_fitB.SetMinimum(20)
legacy_fitE.SetMinimum(20)

legacy_B.SetMarkerColor(ROOT.kOrange + 5)
legacy_E.SetMarkerColor(ROOT.kOrange + 5)
legacy_All.SetMarkerColor(ROOT.kOrange + 5)

legacy_B.SetLineColor(ROOT.kOrange + 5)
legacy_E.SetLineColor(ROOT.kOrange + 5)
legacy_All.SetLineColor(ROOT.kOrange + 5)

legacy_fitB.SetLineColor(ROOT.kOrange + 5)
legacy_fitE.SetLineColor(ROOT.kOrange + 5)
legacy_fitAll.SetLineColor(ROOT.kOrange + 5)

legacy_B.SetMarkerStyle(32)
legacy_E.SetMarkerStyle(32)
legacy_All.SetMarkerStyle(32)

inputFile = ROOT.TFile.Open("/home/llr/cms/cadamuro/Level1_Stage2/TurnOnFitter/UMF/FitEfficiency/results/TurnOn_LU70_NoHPSIso_AlwaysCBM_NoL1Iso/Tau_stage1_stage2_EB_EE_All_Iso.root") 
    
histo_noiso_EB = inputFile.Get("histo_Stage2_Barrel_vs_Pt")
histo_noiso_EB.__class__ = ROOT.RooHist
histo_noiso_EE = inputFile.Get("histo_Stage2_Endcaps_vs_Pt")
histo_noiso_EE.__class__ = ROOT.RooHist
histo_noiso_All = inputFile.Get("histo_Stage2_All_vs_Pt")
histo_noiso_All.__class__ = ROOT.RooHist
    
fit_noiso_EB   = inputFile.Get("fit_Stage2_Barrel_vs_Pt")
fit_noiso_EB.__class__ = ROOT.RooCurve
fit_noiso_EE   = inputFile.Get("fit_Stage2_Endcaps_vs_Pt")
fit_noiso_EE.__class__ = ROOT.RooCurve
fit_noiso_All   = inputFile.Get("fit_Stage2_All_vs_Pt")
fit_noiso_All.__class__ = ROOT.RooCurve

fit_noiso_EB.SetLineWidth(2)
fit_noiso_EE.SetLineWidth(2)
fit_noiso_All.SetLineWidth(2)
fit_noiso_EB.SetMarkerSize(2)
fit_noiso_EE.SetMarkerSize(2) #non so perche ma e la marjer size che controla la linea del fit
fit_noiso_All.SetMarkerSize(2)
fit_noiso_EB.SetLineColor(ROOT.kBlack)
fit_noiso_EE.SetLineColor(ROOT.kBlack)
fit_noiso_All.SetLineColor(ROOT.kBlack)

# Stage 2 files --> fill lists with fit and histo for each WP
WPs = [90,80,70]
colors = [ROOT.kGreen-3, ROOT.kBlue, ROOT.kRed]

histo_B = []
histo_E = []
histo_All = []
fit_B = []
fit_E = []
fit_All = []
histo_shape_B = []
histo_shape_E = []
histo_shape_All = []
fit_shape_B = []
fit_shape_E = []
fit_shape_All = []
# only iso
for index, WP in enumerate(WPs):
    inputFile = ROOT.TFile.Open("/home/llr/cms/cadamuro/Level1_Stage2/TurnOnFitter/UMF/FitEfficiency/results/TurnOn_LU{}_NoHPSIso_AlwaysCBM/Tau_stage1_stage2_EB_EE_All_Iso.root".format(WP)) # gg+VBF training
    #inputFile = ROOT.TFile.Open("./TurnOn_LU{}_NoHPSIso_AlwaysCBM_ShapeVeto/Tau_stage1_stage2_EB_EE_All_Iso.root".format(WP)) # gg+VBF training + shape veto
    
    histo_EB = inputFile.Get("histo_Stage2_Barrel_vs_Pt")
    histo_EB.__class__ = ROOT.RooHist
    histo_EE = inputFile.Get("histo_Stage2_Endcaps_vs_Pt")
    histo_EE.__class__ = ROOT.RooHist
    histo_AllThis = inputFile.Get("histo_Stage2_All_vs_Pt")
    histo_AllThis.__class__ = ROOT.RooHist

    fit_EB   = inputFile.Get("fit_Stage2_Barrel_vs_Pt")
    fit_EB.__class__ = ROOT.RooCurve
    fit_EE   = inputFile.Get("fit_Stage2_Endcaps_vs_Pt")
    fit_EE.__class__ = ROOT.RooCurve
    fit_AllThis   = inputFile.Get("fit_Stage2_All_vs_Pt")
    fit_AllThis.__class__ = ROOT.RooCurve

    histo_B.append(histo_EB)
    histo_E.append(histo_EE)
    histo_All.append(histo_AllThis)


    fit_B.append(fit_EB)
    fit_E.append(fit_EE)
    fit_All.append(fit_AllThis)

    histo_B[-1].SetName ("histo_B_%i" % index)
    histo_E[-1].SetName ("histo_E_%i" % index)
    histo_All[-1].SetName ("histo_All_%i" % index)
    fit_B[-1].SetName ("fit_B_%i" % index)
    fit_E[-1].SetName ("fit_E_%i" % index)
    fit_All[-1].SetName ("fit_All_%i" % index)


    colore = colors[index]
    histo_B[-1].SetLineColor(colore)
    histo_B[-1].SetMarkerColor(colore)

    histo_E[-1].SetLineColor(colore)
    histo_E[-1].SetMarkerColor(colore)

    histo_All[-1].SetLineColor(colore)
    histo_All[-1].SetMarkerColor(colore)
    
    fit_B[-1].SetLineColor(colore)
    fit_B[-1].SetMarkerColor(colore)

    fit_E[-1].SetLineColor(colore)
    fit_E[-1].SetMarkerColor(colore)

    fit_All[-1].SetLineColor(colore)
    fit_All[-1].SetMarkerColor(colore)

#also shape
for index, WP in enumerate(WPs):
    #inputFile = ROOT.TFile.Open("./TurnOn_gg_AllSampleTraining_{}/Tau_stage1_stage2_EB_EE_All_Iso.root".format(WP)) # all trainings!
    #inputFile = ROOT.TFile.Open("./TurnOn_LU{}_NoHPSIso_AlwaysCBM/Tau_stage1_stage2_EB_EE_All_Iso.root".format(WP)) # gg+VBF training
    inputFile = ROOT.TFile.Open("/home/llr/cms/cadamuro/Level1_Stage2/TurnOnFitter/UMF/FitEfficiency/results/TurnOn_LU{}_NoHPSIso_AlwaysCBM_ShapeVeto/Tau_stage1_stage2_EB_EE_All_Iso.root".format(WP)) # gg+VBF training + shape veto
    
    histo_EB = inputFile.Get("histo_Stage2_Barrel_vs_Pt")
    histo_EB.__class__ = ROOT.RooHist
    histo_EE = inputFile.Get("histo_Stage2_Endcaps_vs_Pt")
    histo_EE.__class__ = ROOT.RooHist
    histo_AllThis = inputFile.Get("histo_Stage2_All_vs_Pt")
    histo_AllThis.__class__ = ROOT.RooHist
    
    fit_EB   = inputFile.Get("fit_Stage2_Barrel_vs_Pt")
    fit_EB.__class__ = ROOT.RooCurve
    fit_EE   = inputFile.Get("fit_Stage2_Endcaps_vs_Pt")
    fit_EE.__class__ = ROOT.RooCurve
    fit_AllThis   = inputFile.Get("fit_Stage2_All_vs_Pt")
    fit_AllThis.__class__ = ROOT.RooCurve

    histo_shape_B.append(histo_EB)
    histo_shape_E.append(histo_EE)
    histo_shape_All.append(histo_AllThis)

    fit_shape_B.append(fit_EB)
    fit_shape_E.append(fit_EE)
    fit_shape_All.append(fit_AllThis)

    histo_shape_B[-1].SetName ("histo_shape_B_%i" % index)
    histo_shape_E[-1].SetName ("histo_shape_E_%i" % index)
    histo_shape_All[-1].SetName ("histo_shape_All_%i" % index)
    fit_shape_B[-1].SetName ("fit_shape_B_%i" % index)
    fit_shape_E[-1].SetName ("fit_shape_E_%i" % index)
    fit_shape_All[-1].SetName ("fit_shape_All_%i" % index)

    colore = colors[index]
    histo_shape_B[-1].SetLineColor(colore)
    histo_shape_B[-1].SetMarkerColor(colore)

    histo_shape_E[-1].SetLineColor(colore)
    histo_shape_E[-1].SetMarkerColor(colore)
    
    histo_shape_All[-1].SetLineColor(colore)
    histo_shape_All[-1].SetMarkerColor(colore)

    histo_shape_B[-1].SetLineStyle(7)
    histo_shape_B[-1].SetMarkerStyle(24)
    histo_shape_E[-1].SetLineStyle(7)
    histo_shape_E[-1].SetMarkerStyle(24)
    histo_shape_All[-1].SetLineStyle(7)
    histo_shape_All[-1].SetMarkerStyle(24)

    fit_shape_B[-1].SetLineColor(colore)
    fit_shape_B[-1].SetLineStyle(7)
    fit_shape_B[-1].SetMarkerColor(colore)


    fit_shape_E[-1].SetLineColor(colore)
    fit_shape_E[-1].SetLineStyle(7)
    fit_shape_E[-1].SetMarkerColor(colore)

    fit_shape_All[-1].SetLineColor(colore)
    fit_shape_All[-1].SetLineStyle(7)
    fit_shape_All[-1].SetMarkerColor(colore)

# ================================================================================
# ================================================================================
# ================================================================================

# legacy_B.Draw()
canvBarrel = ROOT.TCanvas ("canvBarrel", "canvBarrel", 600, 600)

# plot barrel only

# plot header: CMS simulation...
legend1 = ROOT.TLegend (0.05, 0.926, 0.887, 0.996)
legend1.AddEntry("NULL","CMS Simulation 2015: gg #rightarrow H #rightarrow #tau #tau - #sqrt{s}=13 TeV, bx=25ns, <PU>=40","h")
legend1.SetFillStyle(0)
legend1.SetBorderSize(0)
legend1.SetTextSize(0.030)
legend1.SetTextFont(62)

# legend for threshold
legendTh = ROOT.TLegend (0.58, 0.17, 0.91, 0.23)
legendTh.AddEntry("NULL","L1 threshold: 30 GeV","h")
legendTh.SetFillStyle(0)
legendTh.SetBorderSize(0)
legendTh.SetTextSize(0.030)

#legend for barrel / endcap
legendBar = ROOT.TLegend (0.58, 0.13, 0.91, 0.17)
legendBar.AddEntry("NULL","Barrel","h")
legendBar.SetFillStyle(0)
legendBar.SetBorderSize(0)
legendBar.SetTextSize(0.030)

#legend for barrel / endcap
legendEnd = ROOT.TLegend (0.58, 0.13, 0.91, 0.17)
legendEnd.AddEntry("NULL","Endcap","h")
legendEnd.SetFillStyle(0)
legendEnd.SetBorderSize(0)
legendEnd.SetTextSize(0.030)

#legend for barrel / endcap / all
legendAll = ROOT.TLegend (0.58, 0.13, 0.91, 0.17)
legendAll.AddEntry("NULL","Inclusive","h")
legendAll.SetFillStyle(0)
legendAll.SetBorderSize(0)
legendAll.SetTextSize(0.030)


# legend for curve meaning -- split in two columns, it is filled by line --> weird order in add entries...

fakeHistoShape = ROOT.TH1D ("fakeHistoShape", "", 1, 0, 1)
fakeHistoShape.SetLineColor(ROOT.kBlack)
fakeHistoShape.SetMarkerColor(ROOT.kBlack)
fakeHistoShape.SetMarkerStyle(24)
fakeHistoShape.SetLineStyle(7)

legendCu = ROOT.TLegend (0.132, 0.758, 0.886, 0.869)
legendCu.SetNColumns(2)

legendCu.AddEntry(histo_noiso_EB,"Upg. 2016 no iso","pl")
legendCu.AddEntry(histo_B[0],"Upg. 2016 iso WP 90%","pl")

legendCu.AddEntry(fakeHistoShape,"Upg. 2016 iso + shape veto","pl")
legendCu.AddEntry(histo_B[1],"Upg. 2016 iso WP 80%","pl")

legendCu.AddEntry(legacy_B,"Run I rescaled","pl")
legendCu.AddEntry(histo_B[2],"Upg. 2016 iso WP 70%","pl")

legendCu.SetFillStyle(0)
legendCu.SetBorderSize(0)
legendCu.SetTextSize(0.028)

#####
## DRAW!
hDummy.Draw()
histo_noiso_EB.Draw("p same")
fit_noiso_EB.Draw("l same")
for index in range(0, len(histo_B)):
    histo_B[index].Draw("p same")
    fit_B[index].Draw("l same")
    histo_shape_B[index].Draw("p same")
    fit_shape_B[index].Draw("l same")
legacy_B.Draw("p same")
#legacy_fitB.Draw("l same")
legend1.Draw()
legendTh.Draw()
legendCu.Draw()
legendBar.Draw()
canvBarrel.Print ("TurnOnPlots/barrel.pdf", "pdf")


# plot endcap only
canvEndcap = ROOT.TCanvas ("canvEndcap", "canvEndcap", 600, 600)
hDummy.Draw()
histo_noiso_EE.Draw("p same")
fit_noiso_EE.Draw("l same")
for index in range(0, len(histo_B)):
    histo_E[index].Draw("p same")
    fit_E[index].Draw("l same")
    histo_shape_E[index].Draw("p same")
    fit_shape_E[index].Draw("l same")
legacy_E.Draw("p same")
#legacy_fitE.Draw("l same")
legend1.Draw()
legendTh.Draw()
legendCu.Draw()
legendEnd.Draw()
canvEndcap.Print ("TurnOnPlots/endcap.pdf", "pdf")


#plot barrel + endcao
# plot endcap only
canvAll = ROOT.TCanvas ("canvAll", "canvAll", 600, 600)
hDummy.Draw()
histo_noiso_All.Draw("p same")
fit_noiso_All.Draw("l same")
for index in range(0, len(histo_All)):
    histo_All[index].Draw("p same")
    fit_All[index].Draw("l same")
    histo_shape_All[index].Draw("p same")
    fit_shape_All[index].Draw("l same")
legacy_All.Draw("p same")
#legacy_fitE.Draw("l same")
legend1.Draw()
legendTh.Draw()
legendCu.Draw()
legendAll.Draw()
canvAll.Print ("TurnOnPlots/All.pdf", "pdf")



raw_input()



