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
hDummy.SetAxisRange(0, 1.05, "Y")
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
legacy_fitB = legacyFile.Get("fit_Stage1_Barrel_vs_Pt")
legacy_fitB.__class__ = ROOT.RooCurve
legacy_fitE = legacyFile.Get("fit_Stage1_Endcaps_vs_Pt")
legacy_fitE.__class__ = ROOT.RooCurve

legacy_B.SetMarkerColor(ROOT.kOrange + 2)
legacy_E.SetMarkerColor(ROOT.kOrange + 2)
legacy_B.SetLineColor(ROOT.kOrange + 2)
legacy_E.SetLineColor(ROOT.kOrange + 2)
legacy_fitB.SetLineColor(ROOT.kOrange + 2)
legacy_fitE.SetLineColor(ROOT.kOrange + 2)
legacy_B.SetMarkerStyle(32)
legacy_E.SetMarkerStyle(32)

# non isolated Stage 2
inputFile = ROOT.TFile.Open("/home/llr/cms/cadamuro/Level1_Stage2/TurnOnFitter/UMF/FitEfficiency/results/TurnOn_LU70_NoHPSIso_AlwaysCBM_NoL1Iso/Tau_stage1_stage2_EB_EE_All_Iso.root") 
    
histo_noiso_EB = inputFile.Get("histo_Stage2_Barrel_vs_Pt")
histo_noiso_EB.__class__ = ROOT.RooHist
histo_noiso_EE = inputFile.Get("histo_Stage2_Endcaps_vs_Pt")
histo_noiso_EE.__class__ = ROOT.RooHist

    
fit_noiso_EB   = inputFile.Get("fit_Stage2_Barrel_vs_Pt")
fit_noiso_EB.__class__ = ROOT.RooCurve
fit_noiso_EE   = inputFile.Get("fit_Stage2_Endcaps_vs_Pt")
fit_noiso_EE.__class__ = ROOT.RooCurve

fit_noiso_EB.SetLineWidth(2)
fit_noiso_EE.SetLineWidth(2)
fit_noiso_EB.SetMarkerSize(2)
fit_noiso_EE.SetMarkerSize(2) #non so perche ma e la marjer size che controla la linea del fit

# Stage 2 files --> fill lists with fit and histo for each WP
WPs = [90,80,70]
colors = [ROOT.kGreen-3, ROOT.kBlue, ROOT.kRed]

histo_B = []
histo_E = []
fit_B = []
fit_E = []
histo_shape_B = []
histo_shape_E = []
fit_shape_B = []
fit_shape_E = []

# only iso
for index, WP in enumerate(WPs):
    inputFile = ROOT.TFile.Open("/home/llr/cms/cadamuro/Level1_Stage2/TurnOnFitter/UMF/FitEfficiency/results/TurnOn_LU{}_NoHPSIso_AlwaysCBM/Tau_stage1_stage2_EB_EE_All_Iso.root".format(WP)) # gg+VBF training
    #inputFile = ROOT.TFile.Open("./TurnOn_LU{}_NoHPSIso_AlwaysCBM_ShapeVeto/Tau_stage1_stage2_EB_EE_All_Iso.root".format(WP)) # gg+VBF training + shape veto
    
    histo_EB = inputFile.Get("histo_Stage2_Barrel_vs_Pt")
    histo_EB.__class__ = ROOT.RooHist
    histo_EE = inputFile.Get("histo_Stage2_Endcaps_vs_Pt")
    histo_EE.__class__ = ROOT.RooHist
    
    fit_EB   = inputFile.Get("fit_Stage2_Barrel_vs_Pt")
    fit_EB.__class__ = ROOT.RooCurve
    fit_EE   = inputFile.Get("fit_Stage2_Endcaps_vs_Pt")
    fit_EE.__class__ = ROOT.RooCurve

    histo_B.append(histo_EB)
    histo_E.append(histo_EE)

    fit_B.append(fit_EB)
    fit_E.append(fit_EE)

    histo_B[-1].SetName ("histo_B_%i" % index)
    histo_E[-1].SetName ("histo_E_%i" % index)
    fit_B[-1].SetName ("fit_B_%i" % index)
    fit_E[-1].SetName ("fit_E_%i" % index)


    colore = colors[index]
    histo_B[-1].SetLineColor(colore)
    histo_B[-1].SetMarkerColor(colore)

    histo_E[-1].SetLineColor(colore)
    histo_E[-1].SetMarkerColor(colore)

    fit_B[-1].SetLineColor(colore)
    fit_B[-1].SetMarkerColor(colore)

    fit_E[-1].SetLineColor(colore)
    fit_E[-1].SetMarkerColor(colore)

#also shape
for index, WP in enumerate(WPs):
    #inputFile = ROOT.TFile.Open("./TurnOn_gg_AllSampleTraining_{}/Tau_stage1_stage2_EB_EE_All_Iso.root".format(WP)) # all trainings!
    #inputFile = ROOT.TFile.Open("./TurnOn_LU{}_NoHPSIso_AlwaysCBM/Tau_stage1_stage2_EB_EE_All_Iso.root".format(WP)) # gg+VBF training
    inputFile = ROOT.TFile.Open("/home/llr/cms/cadamuro/Level1_Stage2/TurnOnFitter/UMF/FitEfficiency/results/TurnOn_LU{}_NoHPSIso_AlwaysCBM_ShapeVeto/Tau_stage1_stage2_EB_EE_All_Iso.root".format(WP)) # gg+VBF training + shape veto
    
    histo_EB = inputFile.Get("histo_Stage2_Barrel_vs_Pt")
    histo_EB.__class__ = ROOT.RooHist
    histo_EE = inputFile.Get("histo_Stage2_Endcaps_vs_Pt")
    histo_EE.__class__ = ROOT.RooHist

    
    fit_EB   = inputFile.Get("fit_Stage2_Barrel_vs_Pt")
    fit_EB.__class__ = ROOT.RooCurve
    fit_EE   = inputFile.Get("fit_Stage2_Endcaps_vs_Pt")
    fit_EE.__class__ = ROOT.RooCurve

    histo_shape_B.append(histo_EB)
    histo_shape_E.append(histo_EE)

    fit_shape_B.append(fit_EB)
    fit_shape_E.append(fit_EE)

    histo_shape_B[-1].SetName ("histo_shape_B_%i" % index)
    histo_shape_E[-1].SetName ("histo_shape_E_%i" % index)
    fit_shape_B[-1].SetName ("fit_shape_B_%i" % index)
    fit_shape_E[-1].SetName ("fit_shape_E_%i" % index)

    colore = colors[index]
    histo_shape_B[-1].SetLineColor(colore)
    histo_shape_B[-1].SetMarkerColor(colore)

    histo_shape_E[-1].SetLineColor(colore)
    histo_shape_E[-1].SetMarkerColor(colore)

    histo_shape_B[-1].SetLineStyle(7)
    histo_shape_B[-1].SetMarkerStyle(24)
    histo_shape_E[-1].SetLineStyle(7)
    histo_shape_E[-1].SetMarkerStyle(24)

    fit_shape_B[-1].SetLineColor(colore)
    fit_shape_B[-1].SetLineStyle(7)
    fit_shape_B[-1].SetMarkerColor(colore)


    fit_shape_E[-1].SetLineColor(colore)
    fit_shape_E[-1].SetLineStyle(7)
    fit_shape_E[-1].SetMarkerColor(colore)

# ================================================================================
# ================================================================================
# ================================================================================

############
#compare iso e single WP for barrel/endcap (no shape, too messy!)
canvCompare = ROOT.TCanvas ("canvCompare", "canvCompare", 600, 600)

# plot header: CMS simulation...
legend1 = ROOT.TLegend (0.05, 0.926, 0.887, 0.996)
legend1.AddEntry("NULL","CMS Simulation 2015: gg #rightarrow H #rightarrow #tau #tau - #sqrt{s}=13 TeV, bx=25ns, <PU>=40","h")
legend1.SetFillStyle(0)
legend1.SetBorderSize(0)
legend1.SetTextSize(0.030)
legend1.SetTextFont(62)

# legend for threshold
legendTh = ROOT.TLegend (0.14, 0.82, 0.47, 0.88)
legendTh.AddEntry("NULL","L1 threshold: 30 GeV","h")
legendTh.SetFillStyle(0)
legendTh.SetBorderSize(0)
legendTh.SetTextSize(0.030)

hDummy.Draw()
histo_noiso_EB.SetLineColor(ROOT.kBlack)
histo_noiso_EB.SetMarkerColor(ROOT.kBlack)
fit_noiso_EB.SetLineColor(ROOT.kBlack)
fit_noiso_EB.SetMarkerColor(ROOT.kBlack)

histo_noiso_EB.Draw("p same")
fit_noiso_EB.Draw("l same")

#change endcap style
fit_noiso_EE.SetLineStyle(7)
histo_noiso_EE.SetLineStyle(7)
histo_noiso_EE.SetMarkerStyle(25)
fit_noiso_EE.SetLineColor(ROOT.kAzure)
histo_noiso_EE.SetLineColor(ROOT.kAzure)
histo_noiso_EE.SetMarkerColor(ROOT.kAzure)
fit_noiso_EE.SetMarkerColor(ROOT.kAzure)

histo_noiso_EE.Draw("p same")
fit_noiso_EE.Draw("l same")


histo_B[0].SetLineColor(ROOT.kRed)
histo_B[0].SetMarkerColor(ROOT.kRed)
fit_B[0].SetLineColor(ROOT.kRed)
fit_B[0].SetMarkerColor(ROOT.kRed)
histo_B[0].Draw("p same")
fit_B[0].Draw("l same")
fit_E[0].SetLineStyle(7)
histo_E[0].SetLineStyle(7)
histo_E[0].SetMarkerStyle(25)
histo_E[0].SetLineColor(ROOT.kGreen-3)
histo_E[0].SetMarkerColor(ROOT.kGreen-3)
fit_E[0].SetLineColor(ROOT.kGreen-3)
fit_E[0].SetMarkerColor(ROOT.kGreen-3)
histo_E[0].Draw("p same")
fit_E[0].Draw("l same")
    
legendComp = ROOT.TLegend (0.43, 0.30, 0.90, 0.51)
legendComp.AddEntry("NULL","L1 threshold: 30 GeV","h")
legendComp.AddEntry(histo_noiso_EB,"No isolation, barrel","pl")
legendComp.AddEntry(histo_noiso_EE,"No isolation, endcap","pl")

legendComp.AddEntry(histo_B[0],"Isolation WP %i%%, barrel" % WPs[0],"pl")
legendComp.AddEntry(histo_E[0],"Isolation WP %i%%, endcap" % WPs[0],"pl")

legendComp.SetFillStyle(0)
legendComp.SetBorderSize(0)
legendComp.SetTextSize(0.030)

legend1.Draw()
#legendTh.Draw()
legendComp.Draw()
canvCompare.Print ("TurnOnPlots/end_bar_compare.pdf", "pdf")


raw_input()



