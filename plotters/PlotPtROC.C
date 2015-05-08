#include "LucaStyle.h"

TGraph* DrawROC (TFile* fback, TFile* fsig)
{
	TH1D* hback = (TH1D*) fback->Get("RelativeRate_diTau");
	TH1D* hsig = (TH1D*) fsig->Get("h_Efficiency");
	
	TGraph* gr = new TGraph;
	gr->SetLineWidth(3);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(1);

	int counter = 0;
	for (int i = 5; i < 71; i++) //do in range 5-70 GeV
	{
		gr -> SetPoint (counter, hsig->GetBinContent (i), 1.-hback->GetBinContent(i));
		counter++;
	}
	return gr;
	
}

TGraph* DrawROCNoIso (TFile* fback, TFile* fsig)
{
	TH1D* hback = (TH1D*) fback->Get("RelativeRate_diTau_NoIso");
	TH1D* hsig = (TH1D*) fsig->Get("h_Efficiency");
	
	TGraph* gr = new TGraph;
	gr->SetLineWidth(2);
	//gr->SetLineStyle(2);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(1);
	
	int counter = 0;
	for (int i = 5; i < 71; i++) //do in range 5-70 GeV
	{
		gr -> SetPoint (counter, hsig->GetBinContent (i), 1.-hback->GetBinContent(i));
		counter++;
	}
	return gr;
	
}

TMarker* DrawMarkerWP (int WP, TFile* fback, TFile* fsig, int markerStyle, int color, int isIso = 0)
{
    TH1D* hback;
    if (isIso == 0) hback = (TH1D*) fback->Get("RelativeRate_diTau");
	else hback = (TH1D*) fback->Get("RelativeRate_diTau_NoIso");
	TH1D* hsig = (TH1D*) fsig->Get("h_Efficiency");

    TMarker* m = new TMarker (hsig->GetBinContent(WP), 1.-hback->GetBinContent(WP), 0);
    m->SetMarkerStyle(markerStyle);
    m->SetMarkerSize(1.5);
    m->SetMarkerColor(color);
    return m;
}

void PlotPtROC()
{
	TGaxis::SetMaxDigits(2);
	//gStyle->SetOptStat(0);
	SetLucaStyle(0);


	// files di rate
	TFile* f70 = new TFile ("/home/llr/cms/cadamuro/Level1_Stage2/CMSSW_7_1_0_pre8/src/macros/Htautau_fullyHad/ggFusionSamples/RateEvaluator/rateL1Tau_70eff_alsoShape_singleShapesVeto_eta2p1_ForEPS.root");
	TFile* f80 = new TFile ("/home/llr/cms/cadamuro/Level1_Stage2/CMSSW_7_1_0_pre8/src/macros/Htautau_fullyHad/ggFusionSamples/RateEvaluator/rateL1Tau_80eff_alsoShape_singleShapesVeto_eta2p1_ForEPS.root");
	TFile* f90 = new TFile ("/home/llr/cms/cadamuro/Level1_Stage2/CMSSW_7_1_0_pre8/src/macros/Htautau_fullyHad/ggFusionSamples/RateEvaluator/rateL1Tau_90eff_alsoShape_singleShapesVeto_eta2p1_ForEPS.root");
	
	// files di efficienza
	TFile* f70sig = new TFile ("../Stage2Trigger/eff_L1Tau_Stage2_gg_LUTProcessed70_AlwaysCBMWithHpsIso.root");
	TFile* f80sig = new TFile ("../Stage2Trigger/eff_L1Tau_Stage2_gg_LUTProcessed80_AlwaysCBMWithHpsIso.root");
	TFile* f90sig = new TFile ("../Stage2Trigger/eff_L1Tau_Stage2_gg_LUTProcessed90_AlwaysCBMWithHpsIso.root");
	
	TFile* fNoIsoSig = new TFile ("../Stage2Trigger/eff_L1Tau_Stage2_gg_NoL1Iso_AlwaysCBMWithHpsIso.root");

	TFile* fRunI = new TFile ("../RunILegacyTrigger/rateL1Tau_L1Legacy_ZeroBiasPU50bx25E13TeV_forEPSRunIScaled.root");
	TFile* fRunISig = new TFile ("../RunILegacyTrigger/eff_L1Tau_RunI_gg_WithHpsIso.root");
	
	TGraph* gr70 = DrawROC(f70, f70sig);
	TGraph* gr80 = DrawROC(f80, f80sig);
	TGraph* gr90 = DrawROC(f90, f90sig);
	TGraph* grNoIso = DrawROCNoIso (f70, fNoIsoSig);
	TGraph* grRunI = DrawROCNoIso (fRunI, fRunISig);
	
	gr70->SetLineColor(kRed);
	gr80->SetLineColor(kBlue);
	gr90->SetLineColor(kGreen-3);
	grNoIso->SetLineColor(kBlack);
	grRunI->SetLineColor(kOrange+5);
	
	TMarker* cut25_eff70 = DrawMarkerWP (25, f70, f70sig, 24, kRed);
	TMarker* cut25_eff80 = DrawMarkerWP (25, f80, f80sig, 24, kBlue);
	TMarker* cut25_eff90 = DrawMarkerWP (25, f90, f90sig, 24, kGreen-3);
	TMarker* cut25_RunI = DrawMarkerWP (25, fRunI, fRunISig, 24, kOrange+5);
	TMarker* cut25_NoIso = DrawMarkerWP (25, f70, fNoIsoSig, 24, kBlack, 1 );

	TMarker* cut30_eff70 = DrawMarkerWP (30, f70, f70sig, 25, kRed);
	TMarker* cut30_eff80 = DrawMarkerWP (30, f80, f80sig, 25, kBlue);
	TMarker* cut30_eff90 = DrawMarkerWP (30, f90, f90sig, 25, kGreen-3);
	TMarker* cut30_RunI = DrawMarkerWP (30, fRunI, fRunISig, 25, kOrange+5);
	TMarker* cut30_NoIso = DrawMarkerWP (30, f70, fNoIsoSig, 25, kBlack, 1);

	TMarker* cut35_eff70 = DrawMarkerWP (35, f70, f70sig, 26, kRed);
	TMarker* cut35_eff80 = DrawMarkerWP (35, f80, f80sig, 26, kBlue);
	TMarker* cut35_eff90 = DrawMarkerWP (35, f90, f90sig, 26, kGreen-3);
	TMarker* cut35_RunI = DrawMarkerWP (35, fRunI, fRunISig, 26, kOrange+5);
	TMarker* cut35_NoIso = DrawMarkerWP (35, f70, fNoIsoSig, 26, kBlack, 1);

	TMarker* cut40_eff70 = DrawMarkerWP (40, f70, f70sig, 30, kRed);
	TMarker* cut40_eff80 = DrawMarkerWP (40, f80, f80sig, 30, kBlue);
	TMarker* cut40_eff90 = DrawMarkerWP (40, f90, f90sig, 30, kGreen-3);
	TMarker* cut40_RunI = DrawMarkerWP (40, fRunI, fRunISig, 30, kOrange+5);
	TMarker* cut40_NoIso = DrawMarkerWP (40, f70, fNoIsoSig, 30, kBlack, 1);

	// legenda
	TLegend* leg = new TLegend (0.597, 0.38, 0.88, 0.89);
	//leg ->SetFillColor(kWhite);
	//leg->SetBorderSize(1);
	leg->SetFillStyle(0);
	leg->SetTextSize(0.027);
	leg->SetTextFont(62);

	leg->AddEntry (grNoIso, "Upgr. 2016, no iso", "l");
	leg->AddEntry (gr90, "Upgr. 2016, WP 90%", "l");
	leg->AddEntry (gr80, "Upgr. 2016, WP 80%", "l");
	leg->AddEntry (gr70, "Upgr. 2016, WP 70%", "l");
	leg->AddEntry (grRunI, "Run I", "l");
	leg->AddEntry (cut25_NoIso, "L1 E_{t} > 25 GeV", "p");
	leg->AddEntry (cut30_NoIso, "L1 E_{t} > 30 GeV", "p");
	leg->AddEntry (cut35_NoIso, "L1 E_{t} > 35 GeV", "p");
	leg->AddEntry (cut40_NoIso, "L1 E_{t} > 40 GeV", "p");

/*
	leg->SetNColumns(2);
	leg->AddEntry (grNoIso, "Upgr. 2016, no iso", "l");
	leg->AddEntry (cut25_eff70, "L1 pt > 25 GeV", "p");
	leg->AddEntry (gr90, "Upgr. 2016, iso WP 90%", "l");
	leg->AddEntry (cut30_eff70, "L1 pt > 30 GeV", "p");
	leg->AddEntry (gr80, "Upgr. 2016, iso WP 80%", "l");
	leg->AddEntry (cut35_eff70, "L1 pt > 35 GeV", "p");
	leg->AddEntry (gr70, "Upgr. 2016, iso WP 70%", "l");
	leg->AddEntry (cut40_eff70, "L1 pt > 40 GeV", "p");
	leg->AddEntry (grRunI, "Run I", "l");
*/	

	TCanvas* c1 = new TCanvas ("c", "c", 600, 600);
	c1->SetGridx();
	c1->SetGridy();
	gr70->GetHistogram()->SetMinimum(0.999);
	gr70->GetHistogram()->SetMaximum(1.);
	gr70->GetXaxis()->SetLimits(0.0,1);
	gr70->GetYaxis()->SetTitleOffset (1.7);
	gr70->GetXaxis()->SetTitle ("Signal efficiency");
	gr70->GetYaxis()->SetTitle ("Backgr. rejection");
	//gr70->SetTitle ("Double #tau ROC curves for L1 threshold values");
	gr70->SetTitle ("");
	
	gr70->Draw("AC");
	gr80->Draw("C same");
	gr90->Draw("C same");
	grNoIso->Draw("C same");
	grRunI->Draw("C same");

    // draw markers
   	cut25_eff70->Draw("same");
	cut25_eff80->Draw("same");
	cut25_eff90->Draw("same");
	cut25_RunI->Draw("same");
	cut25_NoIso->Draw("same");

	cut30_eff70->Draw("same");
	cut30_eff80->Draw("same");
	cut30_eff90->Draw("same");
	cut30_RunI->Draw("same");
	cut30_NoIso->Draw("same");

	cut35_eff70->Draw("same");
	cut35_eff80->Draw("same");
	cut35_eff90->Draw("same");
	cut35_RunI->Draw("same");
	cut35_NoIso->Draw("same");
	
	cut40_eff70->Draw("same");
	cut40_eff80->Draw("same");
	cut40_eff90->Draw("same");
	cut40_RunI->Draw("same");
	cut40_NoIso->Draw("same");
	
	leg->Draw();


	// header
	TLegend* legend1;
    legend1 = new TLegend (0.05, 0.926, 0.887, 0.996);
    legend1->AddEntry("NULL","CMS Simulation 2015: gg #rightarrow H #rightarrow #tau#tau   #sqrt{s}=13 TeV, bx=25ns, <PU>=40","h");
    legend1->AddEntry("NULL","CMS Simulation 2015: Minimum Bias  #sqrt{s}=13 TeV, bx=25ns, <PU>=40","h");
    legend1->SetTextFont(62);
    //legend1->AddEntry("NULL","L1 Threshold : 30 GeV","h");
    legend1->SetLineColor(0);
    legend1->SetBorderSize(0);
    legend1->SetTextSize(0.030);
    legend1->SetFillColor(0);
    legend1->SetFillStyle(0);

    legend1->Draw();

    c1->Print ("ROC_gg.pdf", "pdf");

}
