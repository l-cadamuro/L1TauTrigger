// get the closer bin that alow a certain rates
int findCloser (TH1D* histo, float rate)
{
    int bin = -1;
    for (int i = 1; i <= histo->GetNbinsX(); i++)
    {
        // empty bins don't count
        if (histo->GetBinContent(i) == 0) continue;
        float content = histo->GetBinContent(i);
        if (content < rate)
        {
            bin = i;
            break;
        }
    }
    
    return bin;
}

void PlotDiTauRates_ForEPS()
{
	float absRateScale = 0.001*(2448. * 299792. / 27.); // [khz] * nbunch * ligh speed / LHC lenght
    //float absRateScale  = 1.;
    
    bool PlotRunI = true;

    gStyle->SetOptStat(0);
	gStyle ->SetOptStat(0);
    gStyle->SetTickLength(0.02,"X");
    gStyle->SetTickLength(0.02,"Y");
	gStyle->SetPadTickY(1);
    gStyle->SetPadTickX(1);
    gStyle->SetPadGridY(1);
    gStyle->SetPadGridX(1);
    //gStyle->SetPadTopMargin(0.05);
    //gStyle->SetPadBottomMargin(0.13);
    gStyle->SetPadLeftMargin(0.16);
    gStyle->SetTitleYOffset(1.4);
    gStyle->SetTitleXOffset(0.8);
    gStyle->SetLabelOffset(0.009, "XYZ");

    // for the moment, use the files already produced for EPS in old folder
	TFile* f70 = new TFile ("/home/llr/cms/cadamuro/Level1_Stage2/CMSSW_7_1_0_pre8/src/macros/Htautau_fullyHad/ggFusionSamples/RateEvaluator/rateL1Tau_70eff_alsoShape_singleShapesVeto_eta2p1_ForEPS.root");
	TFile* f80 = new TFile ("/home/llr/cms/cadamuro/Level1_Stage2/CMSSW_7_1_0_pre8/src/macros/Htautau_fullyHad/ggFusionSamples/RateEvaluator/rateL1Tau_80eff_alsoShape_singleShapesVeto_eta2p1_ForEPS.root");
	TFile* f90 = new TFile ("/home/llr/cms/cadamuro/Level1_Stage2/CMSSW_7_1_0_pre8/src/macros/Htautau_fullyHad/ggFusionSamples/RateEvaluator/rateL1Tau_90eff_alsoShape_singleShapesVeto_eta2p1_ForEPS.root");
	TFile* fStage1Iso = new TFile ("..//RunILegacyTrigger/rateL1Tau_L1Legacy_ZeroBiasPU50bx25E13TeV_forEPSRunIScaled.root");

	TH1D* rate_preIso = (TH1D*) f70->Get("RelativeRate_diTau_NoIso");
	TH1D* r70 = (TH1D*) f70 -> Get ("RelativeRate_diTau");
	TH1D* r80 = (TH1D*) f80 -> Get ("RelativeRate_diTau");
	TH1D* r90 = (TH1D*) f90 -> Get ("RelativeRate_diTau");
	TH1D* rStage1Iso = (TH1D*) fStage1Iso -> Get ("RelativeRate_diTau");

	TH1D* r70shape = (TH1D*) f70 -> Get ("RelativeRate_diTau_IsoAndShape");
	TH1D* r80shape = (TH1D*) f80 -> Get ("RelativeRate_diTau_IsoAndShape");
	TH1D* r90shape = (TH1D*) f90 -> Get ("RelativeRate_diTau_IsoAndShape");
    
    
	// set to 0 all bins under 20 in Run I trigger
	for (int i = 1; i <= 15; i++) rStage1Iso->SetBinContent(i, 0);
	// avoid replicated points --> set to zero if previous bin is the same
	/*
    for (int i = 2; i <= rStage1Iso->GetNbinsX(); i++)
		if (rStage1Iso->GetBinContent(i) == rStage1Iso->GetBinContent(i-1)
			|| rStage1Iso->GetBinContent(i) == rStage1Iso->GetBinContent(i-2)
			|| rStage1Iso->GetBinContent(i) == rStage1Iso->GetBinContent(i-3)) rStage1Iso->SetBinContent(i, 0);
    */

	rStage1Iso -> SetMarkerStyle(8);
	rStage1Iso -> SetMarkerSize(0.8);
	rStage1Iso -> SetMarkerColor(kOrange+2);


	r70->GetXaxis()->SetRange(5, 100);
	r70->GetYaxis()->SetRange(0, 1);

	TCanvas* c1 = new TCanvas ("c", "c", 600, 600);
	rate_preIso->SetLineColor(kBlack);
	r70->SetLineColor(kRed);
	r80->SetLineColor(kBlue);
	r90->SetLineColor(kGreen);
    r70shape->SetLineColor(kRed);
    r80shape->SetLineColor(kBlue);
    r90shape->SetLineColor(kGreen);
	rStage1Iso->SetLineColor(kOrange+2);
	
	rate_preIso->SetLineWidth(2);
	rate_preIso->SetLineStyle(1);
	r70->SetLineWidth(2);
	r80->SetLineWidth(2);
	r90->SetLineWidth(2);
    rStage1Iso->SetLineStyle(2);
	rStage1Iso->SetLineWidth(2);
	r70shape->SetLineWidth(2);
	r80shape->SetLineWidth(2);
	r90shape->SetLineWidth(2);
	r70shape->SetLineStyle(7);
	r80shape->SetLineStyle(7);
	r90shape->SetLineStyle(7);


	r70->SetTitle ("Di-#tau rate reduction");
	r70->GetXaxis()->SetTitle ("L1 threshold [GeV]");
	r70->GetYaxis()->SetTitle ("rate reduction");

	c1->SetLogy();
	c1->SetGridx();
	c1->SetGridy();

	TH1D* fakeHisto = new TH1D (*r90shape);
	fakeHisto -> SetLineColor(kBlack);

	TLegend* l = new TLegend (0.51, 0.59, 0.85, 0.89);
	l->SetFillStyle(0);
	l->SetBorderSize(0);
	//l->SetTextFont(62);
	l->SetTextSize(0.030);
	l->AddEntry (rate_preIso, "No isolation", "l");
	l->AddEntry (r90, "Eff. WP 90%", "l");
	l->AddEntry (r80, "Eff. WP 80%", "l");
	l->AddEntry (r70, "Eff. WP 70%", "l");
	//l->AddEntry (r90shape, "Eff. WP 90% + shape veto", "l");
	//l->AddEntry (r80shape, "Eff. WP 80% + shape veto", "l");
	//l->AddEntry (r70shape, "Eff. WP 70% + shape veto", "l");
	l->AddEntry (fakeHisto, "With shape veto", "l");
	if (PlotRunI) l->AddEntry (rStage1Iso, "Run I scaled (k = 0.628)", "p");
	

	r70->GetXaxis()->SetTitleSize(0.050);
    r70->GetYaxis()->SetTitleSize(0.050);
    r70->GetYaxis()->SetTitleOffset(1.4);
    r70->GetXaxis()->SetTitleOffset(0.9);
	r70->SetTitle("");
	
    // scale all to abs rate
    r70->Scale(absRateScale);
    r80->Scale(absRateScale);
    r90->Scale(absRateScale);
    r70shape->Scale(absRateScale);
    r80shape->Scale(absRateScale);
    r90shape->Scale(absRateScale);
    rate_preIso->Scale(absRateScale);
    rStage1Iso->Scale(absRateScale);
    r70->GetYaxis()->SetTitle("Rate [kHz]");
    r70->GetXaxis()->SetRange(15., 100.);
    r70->SetMinimum(1.e-3);
    r70->SetMaximum(1.e3);

    r70->Draw();
	r80->Draw("same");
	r90->Draw("same");
	r70shape->Draw("same");
	r80shape->Draw("same");
	r90shape->Draw("same");
	rate_preIso->Draw("same");
	if (PlotRunI) rStage1Iso->Draw("Psame");
	l->Draw();

	TLegend* legend1;
    legend1 = new TLegend (0.05, 0.926, 0.887, 0.996);
    legend1->AddEntry("NULL","CMS Simulation: Minimum Bias         #sqrt{s}=13 TeV, bx=25ns, PU=40","h");
    legend1->SetTextFont(62);
    //legend1->AddEntry("NULL","L1 Threshold : 30 GeV","h");
    legend1->SetLineColor(0);
    legend1->SetBorderSize(0);
    legend1->SetTextSize(0.030);
    legend1->SetFillColor(0);
    legend1->SetFillStyle(0);

    legend1->Draw();

    if (PlotRunI)
        c->Print ("rate_WithLegacy.pdf", "pdf");
    else
        c->Print ("rate_NOLegacy.pdf", "pdf");

    
    // get thresholds coreesponding to certain ratres
    float rateThis = 1.;
    cout << " Rate   Run I    NoIso    90   90sh   80    80sh   70s   70sh" << endl;
    cout << rateThis << " " << findCloser (rStage1Iso, rateThis) << " " << findCloser (rate_preIso, rateThis)
         << " " << findCloser (r90, rateThis) << " " << findCloser (r90shape, rateThis)
         << " " << findCloser (r80, rateThis) << " " << findCloser (r80shape, rateThis)
         << " " << findCloser (r70, rateThis) << " " << findCloser (r70shape, rateThis) << endl;

/*
    // print the rate for some interesting thresholds
    int thresholds [] = {20, 25, 30, 35, 40, 45, 50};
    int nThresh = sizeof(thresholds)/sizeof(int);
    
    cout << " Thr   Run I    NoIso    90   80   70" << endl;
    for (int i = 0; i < nThresh; i++)
    {
        cout << i << " " << 
    }
*/    
    
}
