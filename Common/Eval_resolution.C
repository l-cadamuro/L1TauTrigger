// use filtered tau to obtain separately barrel/endcap resolution
//
// all aestethics settings are removed, and put in a dedicated plotter macro

//#include "shape_functions.h" // read flags saved as hwQual
//#include "tdrStyle.h"

// single CB
/*
double CrystalBall
(double* x, double* par){ 
  //http://en.wikipedia.org/wiki/Crystal_Ball_function 
  double xcur = x[0]; 
  double alpha = par[0]; 
  double n = par[1]; 
  double mu = par[2]; 
  double sigma = par[3]; 
  double N = par[4]; 
  TF1* exp = new TF1("exp","exp(x)",1e-20,1e20); 
  double A; double B; 
  if (alpha < 0){ 
    A = pow((n/(-1*alpha)),n)*exp->Eval((-1)*alpha*alpha/2); 
    B = n/(-1*alpha) + alpha;} 
  else { 
    A = pow((n/alpha),n)*exp->Eval((-1)*alpha*alpha/2); 
    B = n/alpha - alpha;} 
double f; 
  if ((xcur-mu)/sigma > (-1)*alpha) 
    f = N*exp->Eval((-1)*(xcur-mu)*(xcur-mu)/(2*sigma*sigma)); 
  else 
    f = N*A*pow((B- (xcur-mu)/sigma),(-1*n)); 
  delete exp; 
  return f; 
} 
*/

// double CB symmetric (two-sided tails)
double CrystalBall
(double* x, double* par){ 
  //http://en.wikipedia.org/wiki/Crystal_Ball_function 
  double xcur = x[0]; 
  double alpha = par[0]; 
  double n = par[1]; 
  double mu = par[2]; 
  double sigma = par[3]; 
  double N = par[4]; 
  TF1* exp = new TF1("exp","exp(x)",1e-20,1e20); 
  double A; double B; 
  if (alpha < 0){ 
    A = pow((n/(-1*alpha)),n)*exp->Eval((-1)*alpha*alpha/2); 
    B = n/(-1*alpha) + alpha;} 
  else { 
    A = pow((n/alpha),n)*exp->Eval((-1)*alpha*alpha/2); 
    B = n/alpha - alpha;} 
    double f; 
  if (TMath::Abs((xcur-mu)/sigma) < TMath::Abs(alpha) ) 
    f = N*exp->Eval((-1)*(xcur-mu)*(xcur-mu)/(2*sigma*sigma)); 
  else if (((xcur-mu)/sigma) < (-1.)*alpha )
    f = N*A*pow((B- (xcur-mu)/sigma),(-1*n)); // left tail
  else
    f = N*A*pow( (B- (mu-xcur)/sigma),(-1*n)); // right tail
  delete exp; 
  return f; 
} 

// double CB asymmetric
double DoubleCrystalBall
(double* x, double* par){ 
  //http://en.wikipedia.org/wiki/Crystal_Ball_function 
  double xcur = x[0]; 
  double alphaL = par[0]; 
  double nL = par[1]; 
  double alphaR = par[2]; 
  double nR = par[3]; 

  double mu = par[4]; 
  double sigma = par[5]; 
  double N = par[6]; 
 
  TF1* exp = new TF1("exp","exp(x)",1e-20,1e20); 
  double AL; double BL; double AR; double BR; 
 
  if (alphaL < 0){ 
    AL = pow((nL/(-1*alphaL)),nL)*exp->Eval((-1)*alphaL*alphaL/2); 
    BL = nL/(-1*alphaL) + alphaL;} 
  else { 
    AL = pow((nL/alphaL),nL)*exp->Eval((-1)*alphaL*alphaL/2); 
    BL = nL/alphaL - alphaL;} 

  if (alphaR < 0){ 
    AR = pow((nR/(-1*alphaR)),nR)*exp->Eval((-1)*alphaR*alphaR/2); 
    BR = nR/(-1*alphaR) + alphaR;} 
  else { 
    AR = pow((nR/alphaR),nR)*exp->Eval((-1)*alphaR*alphaR/2); 
    BR = nR/alphaR - alphaR;} 
   

    double f; 
  if ( ((xcur-mu)/sigma) > (-1.)*alphaL  && ((xcur-mu)/sigma) < (1.)*alphaR) 
    f = N*exp->Eval((-1)*(xcur-mu)*(xcur-mu)/(2*sigma*sigma)); 
  
  // left
  else if ( ((xcur-mu)/sigma) <= (-1.)*alphaL )
    f = N*AL*pow((BL- (xcur-mu)/sigma),(-1*nL)); // left tail
  //right
  else
    f = N*AR*pow( (BR- (mu-xcur)/sigma),(-1*nR)); // right tail
  delete exp; 
  return f; 
} 

// double CB asymmetric
double DoubleCrystalBallDoubleGaus
(double* x, double* par){ 
  //http://en.wikipedia.org/wiki/Crystal_Ball_function 
  double xcur = x[0]; 
  double alphaL = par[0]; 
  double nL = par[1]; 
  double alphaR = par[2]; 
  double nR = par[3]; 

  double mu = par[4]; 
  double sigmaL = par[5]; 
  double sigmaR = par[6]; 

  double N = par[7]; 
 
  TF1* exp = new TF1("exp","exp(x)",1e-20,1e20); 
  double AL; double BL; double AR; double BR; 
 
  if (alphaL < 0){ 
    AL = pow((nL/(-1*alphaL)),nL)*exp->Eval((-1)*alphaL*alphaL/2); 
    BL = nL/(-1*alphaL) + alphaL;} 
  else { 
    AL = pow((nL/alphaL),nL)*exp->Eval((-1)*alphaL*alphaL/2); 
    BL = nL/alphaL - alphaL;} 

  if (alphaR < 0){ 
    AR = pow((nR/(-1*alphaR)),nR)*exp->Eval((-1)*alphaR*alphaR/2); 
    BR = nR/(-1*alphaR) + alphaR;} 
  else { 
    AR = pow((nR/alphaR),nR)*exp->Eval((-1)*alphaR*alphaR/2); 
    BR = nR/alphaR - alphaR;} 
   

  double f; 
  if ( ((xcur-mu)/sigmaL) > (-1.)*alphaL  && ((xcur-mu)/sigmaR) < (1.)*alphaR) 
  {
    // gaussian left side
    if ((xcur-mu) < 0)
        f = N*exp->Eval((-1)*(xcur-mu)*(xcur-mu)/(2*sigmaL*sigmaL)); 
    else
        f = N*exp->Eval((-1)*(xcur-mu)*(xcur-mu)/(2*sigmaR*sigmaR)); 
  }
  // left
  else if ( ((xcur-mu)/sigmaL) <= (-1.)*alphaL )
    f = N*AL*pow((BL- (xcur-mu)/sigmaL),(-1*nL)); // left tail
  //right
  else
    f = N*AR*pow( (BR- (mu-xcur)/sigmaR),(-1*nR)); // right tail
  delete exp; 
  return f; 
} 

double DoubleGaussian (double* x, double* par)
{
    double xcur = x[0];
    double mean = par[0];
    double sigmaL = par[1];
    double sigmaR = par[2];
    double N = par[3];

    TF1* exp = new TF1("exp","exp(x)",1e-20,1e20);
    
    double f;
    if ( (xcur-mean) >= 0 )
        f = N*exp->Eval ((-1)*(xcur-mean)*(xcur-mean)/(2*sigmaL*sigmaL));
    else
        f = N*exp->Eval ((-1)*(xcur-mean)*(xcur-mean)/(2*sigmaR*sigmaR));
    delete exp;
    return f;
}

// get the smallest interval containing 68% of events
double EffectiveRMS (TH1D* h)
{
    int nbins = h->GetNbinsX();
    double effRMS = h->GetBinLowEdge (nbins+1) - h->GetBinLowEdge(1);
    double totEvts = 1.*h->Integral();
    double bufRMS = 0;
    cout << effRMS << endl;

    // loop on bins
    for (int ibin = 1; ibin <= nbins; ibin++)
    {
        //bool RMSfound = false;
        for (int iUp = ibin+1; iUp <= nbins; iUp++)
        {
            double currInt = h->Integral (ibin, iUp);
            if (currInt/totEvts >= 0.68)
            {
                //RMSfound = true;
                bufRMS = h->GetBinLowEdge (iUp) - h->GetBinLowEdge(ibin);
                if (bufRMS < effRMS) effRMS = bufRMS;
                break;
            }
        }
    }

    // convert in real units ??
    return effRMS;

}

void Eval_resolution()
{

    TH1::SetDefaultSumw2();
	//setTDRStyle();
    TF1* CBFunc = new TF1("CBFunc",&CrystalBall,-0.3,0.3,5);
    TF1* CBFuncAsymm = new TF1("CBFuncAsymm",&DoubleCrystalBall,0.,3.,7);
    TF1* CBFuncAsymmDoubleGaus = new TF1("CBFuncAsymmDoubleGaus",&DoubleCrystalBallDoubleGaus,0.,3.,8);

    bool doFits = true;
    bool doFitsAngular = false; // set to false for Run I, clearly not a CB!!
    bool scaleL1Pt = true; // true to scale run I
    double scaleFactor = 1.075/1.713; // mean of the Run I resolution distribution taken from the histogram, num is mean of Stage 2 to put both in the same spot!
    //double fitLowPt = 0.6; double fitHighPt = 1.6; // "different fit range"
    double fitLowPt = 0.8; double fitHighPt = 1.5; // normally used


    //TF1* CBFuncAsymm_clone = new TF1("CBFuncAsymm_clone",&DoubleCrystalBall,0.,3.,7);
    //TF1* CBFuncAsymmDoubleGaus_clone = new TF1("CBFuncAsymmDoubleGaus_clone",&DoubleCrystalBallDoubleGaus,0.,3.,8);
    //CBFuncAsymmDoubleGaus_clone->SetLineColor(kBlue);
    //CBFuncAsymm_clone->SetLineColor(kBlue);
    //TF1* DoubleGaus = new TF1("DoubleGaus",&DoubleGaussian,0.,3.,4);
    
    /*
    CBFunc->SetNpx(1000);
    CBFuncAsymm->SetNpx(1000);
    CBFuncAsymm_clone->SetNpx(1000);
    CBFuncAsymmDoubleGaus -> SetNpx(1000);
    CBFuncAsymmDoubleGaus_clone -> SetNpx(1000);
    

	  gStyle ->SetOptStat(0);
    gStyle->SetTickLength(0.02,"X");
    gStyle->SetTickLength(0.02,"Y");
	  gStyle->SetPadTickY(1);
    gStyle->SetPadTickX(1);
    //gStyle->SetPadTopMargin(0.05);
    //gStyle->SetPadBottomMargin(0.13);
    gStyle->SetPadLeftMargin(0.16);
    gStyle->SetTitleYOffset(1.4);
    gStyle->SetTitleXOffset(0.9);
    gStyle->SetLabelOffset(0.009, "XYZ");

    //gStyle->SetPadRightMargin(0.02);
    //TApplication theApp("App",&argc,argv);
	  */
	// open ROOT file -> filtered paricles!
	
	// uncalib pT
	//TFile* fInput = new TFile ("filtered_taus_WithHPSIsorequirement.root"); // with hps iso
	//TFile* fInput = new TFile ("filtered_taus_NoHPSIsorequirement.root"); // without hps iso
  //TFile* fInput = new TFile ("filtered_taus_WithHPSIsorequirement_ReverseMatch_v7.root");
  //TFile* fInput = new TFile ("filtered_taus_NoHPSIsorequirement_ReverseMatch_v7.root"); // min dR
  //TFile* fInput = new TFile ("filtered_taus_tree_reverse_noHPSIso_maxPt.root"); // max pT --> no change wrt mindR
  //TFile* fInput = new TFile ("filtered_taus_tree_reverse_noHPSIso_noAmbigResol.root"); // no amvbiguities resolution as Luca M was doing
  //TFile* fInput = new TFile ("filtered_taus_gg_noHPSIso_etaOptim_PUS_forEPS.root");
  //TFile* fInput = new TFile ("filtered_taus_VBF_noHPSIso_etaOptim_PUS_forEPS.root");
 
  // firmware optim, eta 2.1, gg, With Hps iso
  //TFile* fInput = new TFile ("/home/llr/cms/cadamuro/Level1_Stage2_PUSubDevel/PUSMacros/filtered_taus_gg_WithHPSIso_eta2p1_PUS_forEPS.root"); 
  //TString fOutName = "resolutions_histo_gg_Stage2_WithHPSIso_eta2p1.root";
  
  // Run I, gg, With Hps Iso
  TFile* fInput = new TFile ("../RunILegacyTrigger/filtered_taus_L1LegacyRunI_gg_bx25Pu40E13TeVWithHpsIso.root"); 
  TString fOutName = "resolutions_histo_gg_RunI_Scaled.root";

	TTree* tInput = (TTree*) fInput->Get("filtered_taus_tree");
	
	// set branch and variables
	float L1_pt, L1_eta, L1_phi;
	int EvtNtt, L1_IsoEt, L1_flags;
	float hps_pt, hps_eta, hps_phi;
	int gen_decaymode;
		
	// branch	
	tInput ->SetBranchAddress ("L1_pt", &L1_pt);    // nTT is stored for each tau, will only take it for one tau in the vector
	tInput ->SetBranchAddress ("L1_eta", &L1_eta);
	tInput ->SetBranchAddress ("L1_phi", &L1_phi);  
	//tInput ->SetBranchAddress ("EvtNtt", &EvtNtt);
	//tInput ->SetBranchAddress ("L1_IsoEt", &L1_IsoEt);
	//tInput ->SetBranchAddress ("L1_flags", &L1_flags);

	tInput ->SetBranchAddress ("hps_pt", &hps_pt); 
	tInput ->SetBranchAddress ("hps_eta", &hps_eta); 
	tInput ->SetBranchAddress ("hps_phi", &hps_phi);   

	tInput->SetBranchAddress ("gen_decaymode", &gen_decaymode);


  ///////////////////////////////////////
  /////// HISTOGRAMS ////////////////////
  ///////////////////////////////////////
  // save all oobjects to ROOT file, it's getting too messy here!!
  TFile* fOut = new TFile (fOutName, "recreate");
  
  // istotgrams
	TH1D* pt_resol_EB =   new TH1D ("pt_resol_EB", "; E_{T}(L1)/p_{T}(offline); a.u.", 100, 0, 3);
	TH1D* pt_resol_EE =	  new TH1D ("pt_resol_EE", "; E_{T}(L1)/p_{T}(offline); a.u.", 100, 0, 3);
	TH1D* pt_resol_OneBin =	  new TH1D ("pt_resol_OneBin", "Energy response; E_{T}(L1)/p_{T}(offline); a.u.", 100, 0, 3);
  TH1D* pt_resol_All =   new TH1D ("pt_resol_All", "; E_{T}(L1)/p_{T}(offline); a.u.", 100, 0, 3);

	TH1D* eta_resol_EB =  new TH1D ("eta_resol_EB", "#eta resolution; #eta(L1) - #eta(offline); a.u.", 100, -0.3, 0.3);
	TH1D* eta_resol_EE =  new TH1D ("eta_resol_EE", "#eta resolution; #eta(L1) - #eta(offline); a.u.", 100, -0.3, 0.3);
	TH1D* eta_resol_All =  new TH1D ("eta_resol_All", "; #eta(L1) - #eta(offline); a.u.", 100, -0.3, 0.3);
	
	TH1D* phi_resol_EB =  new TH1D ("phi_resol_EB", "#phi resolution; #phi(L1) - #phi(offline); a.u.", 100, -0.3, 0.3);
	TH1D* phi_resol_EE =  new TH1D ("phi_resol_EE", "#phi resolution; #phi(L1) - #phi(offline); a.u.", 100, -0.3, 0.3);
	TH1D* phi_resol_All =  new TH1D ("phi_resol_All", "; #phi(L1) - #phi(offline); a.u.", 100, -0.3, 0.3);

	TH1D* pt_resol_EB_diff =   new TH1D ("pt_resol_EB_diff", "pT resolution; pT(L1) - pT(offline); a.u.", 100, -50, 50);
	TH1D* pt_resol_EE_diff =	  new TH1D ("pt_resol_EE_diff", "pT resolution; pT(L1) - pT(offline); a.u.", 100, -50, 50);

  // block for different eta regions, to check if the effects comes for the beamspot positions --> eta resolution is better at high eta
  // probably an effect of beamspot position (there is no difference as a function of phi)
  TH2D* eta_resol_binned = new TH2D ("eta_resol_binned", "#eta bin vs #eta resolution; #eta; resol", 10, 0, 2.1, 100, -0.3, 0.3);

  double lowEta = 0.; double highEta = 2.1;
  //double lowEta = 1.2; double highEta = 1.8;
  double lowPt = 0.; double highPt = 100.;

  TH2D* pT_resol_binned = new TH2D ("pT_resol_binned_hpseta", "eta bin vs pT resolution; #eta; resol", 10, lowEta, highEta, 100, 0, 3);
  TH2D* pT_resol_binned_hpsPt = new TH2D ("pT_resol_binned_hpsPt", "pt bin vs pT resolution; pt offline; resol", 10, lowPt, highPt, 100, 0, 3);

  TH1D* pt_resol_decMode[4];
	TH1D* eta_resol_decMode[4];
	TH1D* phi_resol_decMode[4];



	//int colori [4] = {kRed, kBlue, kGreen, kOrange+2}; 
	for (int i = 0; i < 4; i++)
	{
		pt_resol_decMode[i] = new TH1D (Form("hpt_decaymode_%i", i), "Energy response; pT(L1)/pT(reco); a.u.", 100, 0, 3);
		eta_resol_decMode[i] = new TH1D (Form("heta_decaymode_%i", i), "#eta resolution; #eta(L1) - #eta(reco); a.u.", 100, -0.3, 0.3);
		phi_resol_decMode[i] = new TH1D (Form("hphi_decaymode_%i", i), "#phi resolution; #phi(L1) - #phi(reco); a.u.", 100, -0.3, 0.3);
    // pt_resol_decMode[i]->SetLineColor(colori[i]);
		// eta_resol_decMode[i]->SetLineColor(colori[i]);
		// phi_resol_decMode[i]->SetLineColor(colori[i]);
		// pt_resol_decMode[i]->SetLineWidth(2);
		// eta_resol_decMode[i]->SetLineWidth(2);
		// phi_resol_decMode[i]->SetLineWidth(2);

	}


  /*
	pt_resol_EB ->SetLineColor(kRed);
	eta_resol_EB -> SetLineColor(kRed);
	pt_resol_EB_diff ->SetLineColor(kRed);
  phi_resol_EB ->SetLineColor(kRed);
  */

	int nEvents = tInput->GetEntries(); 


	// loop on all events
	for (int iEv = 0; iEv < nEvents; iEv++)
	{
		tInput->GetEntry(iEv);
		if (iEv%10000 == 0) cout << iEv << endl;

		float absEta = TMath::Abs(L1_eta);
    //int prodMode = GetTauClusterProductionMode (L1_flags); // 0 : can Be Merged, 3: not merged (no neighbours)

		// discard L1 taus with pt < 20 GeV (all , pre rescaling)
		if (L1_pt < 20) continue;
    if (scaleL1Pt) L1_pt *= scaleFactor;
    //if (prodMode != 3) continue;

		// fill istos
		if (TMath::Abs(hps_eta) < 1.305)
		{
			pt_resol_EB -> Fill (L1_pt / hps_pt);
			eta_resol_EB -> Fill (L1_eta - hps_eta);
			pt_resol_EB_diff -> Fill (L1_pt - hps_pt);
      phi_resol_EB ->Fill (L1_phi - hps_phi);
		}		
		else if (TMath::Abs(hps_eta) > 1.479 && TMath::Abs(hps_eta) < 2.5)
		{
			pt_resol_EE -> Fill (L1_pt / hps_pt);
			eta_resol_EE -> Fill (L1_eta - hps_eta);
			pt_resol_EE_diff -> Fill (L1_pt - hps_pt);
      phi_resol_EE ->Fill (L1_phi - hps_phi);
		}
        
		eta_resol_All -> Fill (L1_eta - hps_eta);
		phi_resol_All -> Fill (L1_phi - hps_phi);
    pt_resol_All -> Fill (L1_pt / hps_pt);
		
    if (TMath::Abs(hps_eta) > 0.7 &&  TMath::Abs(hps_eta) < 0.77 ) pt_resol_OneBin->Fill(L1_pt / hps_pt);

    // binned eta
    eta_resol_binned->Fill (TMath::Abs(hps_eta), L1_eta - hps_eta);	
    // don't want outliers
    if ( L1_pt / hps_pt > 0 &&  L1_pt / hps_pt < 3)
    {
        pT_resol_binned ->Fill (TMath::Abs(hps_eta), L1_pt / hps_pt);	
        pT_resol_binned_hpsPt ->Fill (hps_pt, L1_pt / hps_pt);
    }

    // gen decay mode
    int decMode = -1; // switch for binning
    if (gen_decaymode == 0 ) 						decMode = 0;
    if (gen_decaymode == 1 || gen_decaymode == 2)   decMode = 1;
    if (gen_decaymode == 4 ) 						decMode = 2;
    if (gen_decaymode == 5 || gen_decaymode == 6 )  decMode = 3;

    if (decMode != -1)
    {
    	pt_resol_decMode[decMode]->Fill(L1_pt / hps_pt);
      eta_resol_decMode[decMode]->Fill(L1_eta - hps_eta);
      phi_resol_decMode[decMode]->Fill(L1_phi - hps_phi);
    }
	}

//////////////////////////////////////////////////////////////
///////////////////////// NORMALIZE  /////////////////////////
//////////////////////////////////////////////////////////////


  pt_resol_EE->Scale(1./pt_resol_EE->Integral());
  pt_resol_EB->Scale(1./pt_resol_EB->Integral());
  pt_resol_All->Scale(1./pt_resol_All->Integral());

  eta_resol_EE->Scale(1./eta_resol_EE->Integral());
  eta_resol_EB->Scale(1./eta_resol_EB->Integral());
  eta_resol_All->Scale(1./eta_resol_All->Integral());

  phi_resol_EE->Scale(1./phi_resol_EE->Integral());
  phi_resol_EB->Scale(1./phi_resol_EB->Integral());
  phi_resol_All->Scale(1./phi_resol_All->Integral());

  pt_resol_EE_diff->Scale(1./pt_resol_EE_diff->Integral());
  pt_resol_EB_diff->Scale(1./pt_resol_EB_diff->Integral());

//////////////////////////////////////////////////////////////
///////////////////////// MAKE FITS  /////////////////////////
//////////////////////////////////////////////////////////////

/*
	TLegend* l = new TLegend (0.50, 0.63, 0.88, 0.82);
	l->AddEntry (pt_resol_EB, "barrel: |#eta| < 1.305", "lep");
	l->AddEntry (pt_resol_EE, "endcap: |#eta| > 1.479", "lep");
	l->SetFillStyle(0);
    l->SetBorderSize(0);
	
    pt_resol_EE->SetLineWidth(2);
    pt_resol_EB->SetLineWidth(2);
*/

/*
	TCanvas * c1 = new TCanvas ("c1", "c1", 600, 600);
    TLegend* legend1 = new TLegend (0.07, 0.85, 0.9, 1.12);
    legend1->AddEntry("NULL","CMS Simulation: gg #rightarrow H #rightarrow #tau #tau         #sqrt{s}=13 TeV, bx=25ns, PU=40","h");
    //legend1->AddEntry("NULL","L1 Threshold : 30 GeV","h");
    legend1->SetLineColor(0);
    legend1->SetBorderSize(0);
    legend1->SetTextSize(0.030);
    legend1->SetFillColor(0);
    legend1->SetFillStyle(0);
*/
  // prelimianry gaussian + later CB fit
    
  TFitResultPtr gausPrelim = pt_resol_EB ->Fit("gaus", "S0", "", 0.9, 1.3);
  TF1 *myfit = (TF1*) pt_resol_EB->GetFunction("gaus");
  double constant = myfit->GetParameter(0);
  double mean = myfit->GetParameter(1);
  double sigma = myfit->GetParameter(2);

  cout << myfit->GetParameter(0) << " " << myfit->GetParameter(1)<< " " << myfit->GetParameter(2) << endl;
  
  /*
  CBFuncAsymm->SetParLimits(4, mean-0.2*sigma, mean+0.2*sigma);
  CBFuncAsymm->SetParLimits(5, 0.5*sigma, 2.*sigma);
  CBFuncAsymm->SetParLimits(0, 0.1, 1.5.);
  CBFuncAsymm->SetParLimits(1, 0., 30.);
  CBFuncAsymm->SetParLimits(2, 0.1, 3.);
  CBFuncAsymm->SetParLimits(3, 0., 30.);
  */

  CBFuncAsymm->SetParameters(0.9, 4.3, 1.3, 4.3, mean, sigma, constant);
  if (doFits) pt_resol_EB ->Fit("CBFuncAsymm", "", "", 0.6, 1.8);

  //CBFuncAsymmDoubleGaus->SetParameters(1.5, 4, 1, 4, 1., 0.25, 0.25, 0.06);
  //CBFuncAsymmDoubleGaus_clone->SetParameters(1.5, 4, 1, 4, 1., 0.25, 0.25, 0.06);
  //DoubleGaus->SetParameters (1, 1, 1, 0.06);
  //legend1->SetTextFont(62);
    //pt_resol_EB ->Fit("CBFuncAsymm");
    //pt_resol_EB ->Fit("DoubleGaus");
    //TFitResultPtr fitresEB = pt_resol_EB ->Fit("CBFuncAsymmDoubleGaus", "S");
    //TFitResultPtr fitresEE = pt_resol_EE ->Fit("CBFuncAsymmDoubleGaus_clone", "S");
  CBFuncAsymm->SetParameters(1.5, 4, 1, 4, 1, 0.25, 0.06);
  if (doFits) pt_resol_EE ->Fit("CBFuncAsymm", "", "", 0.6, 1.6);    
  if (doFits) pt_resol_All ->Fit("CBFuncAsymm", "", "", 0.6, 1.6);
  // draw this to check ft quality  
  //TCanvas* c1 = new TCanvas ("cc", "cc", 600, 600);
  //pt_resol_EE->Draw();
  //pt_resol_EB->Draw();

/*
    pt_resol_EB ->SetMarkerStyle(20);
    pt_resol_EB ->SetMarkerSize(0.8);
    pt_resol_EB ->SetMarkerColor(kRed+1);

    pt_resol_EE ->SetMarkerStyle(8);
    pt_resol_EE ->SetMarkerSize(0.8);
    pt_resol_EE ->SetMarkerColor(kBlue+1);

    pt_resol_EE->GetXaxis()->SetTitleSize(0.050);
    pt_resol_EE->GetYaxis()->SetTitleSize(0.050);
    pt_resol_EB->GetXaxis()->SetTitleSize(0.050);
    pt_resol_EB->GetYaxis()->SetTitleSize(0.050);

    pt_resol_EB ->Draw();
	pt_resol_EE ->Draw("same");
*/
/*	
cout << "BARREL/ENDCAP RMS/mean: " <<  pt_resol_EB->GetRMS()/pt_resol_EB->GetMean() << " " << pt_resol_EE->GetRMS()/pt_resol_EE->GetMean() << endl; 
    c1->SetGridx();
    c1->SetGridy();
    l->Draw();
    //CBFuncAsymm->Draw();
	legend1->Draw();

    c1->Print ("pt_resol_rel.pdf", "pdf");
*/
	

  //TCanvas* c3 = new TCanvas ("c3", "c3", 600, 600);
	//phi_resol_EE ->Draw();
	//phi_resol_EB ->Draw("same");
	//l->Draw();
    
  //phi_resol_All->Fit("gaus", "N");
  CBFunc->SetParameters(3, 1, 0, 0.05, 0.06);
  if (doFitsAngular)
  {
    phi_resol_All->Fit("CBFunc");
    eta_resol_All->Fit("CBFunc");
  }

  /*  
  //cout << CBFunc->GetParameter (3) << " " << CBFunc->GetParError (3) << phi_resol_All->GetRMS() << endl; //gaus->getParameter (0) <<  endl;
  TLegend* lFitPhi = new TLegend (0.58, 0.60, 0.95, 0.78);
  lFitPhi->SetFillStyle(0);
  lFitPhi->SetBorderSize(0);
  TString textPhi = "#sigma = (3.47 #pm 0.01) #upoint 10^{-2}";
  lFitPhi->AddEntry ("NULL", textPhi, "h");
*/

/*
  phi_resol_All->SetMaximum (0.08);
  phi_resol_All->SetLineWidth(2);
  phi_resol_All->GetXaxis()->SetTitleSize(0.050);
  //phi_resol_All->GetXaxis()->SetTitleOffset(0.8);
  phi_resol_All->GetYaxis()->SetTitleSize(0.050);
  //phi_resol_All->GetYaxis()->SetTitleOffset(0.8);
  //c3->SetLeftMargin(0.90);
  double w = 600;
  double h = 600;
  //c3->SetWindowSize(w + (w - c3->GetWw()), h + (h - c3->GetWh()));

  c3->SetGridx();
  c3->SetGridy();
  phi_resol_All->SetMarkerStyle(8);
  phi_resol_All->SetMarkerSize(0.8);
  phi_resol_All->SetMarkerColor(kBlue);
  phi_resol_All->Draw();
  legend1->Draw();
  lFitPhi->Draw();
c3->Print("phi_resol.pdf", "pdf");
*/
    //TCanvas * c2 = new TCanvas ("c2", "c2", 600, 600);
    //eta_resol_EE ->Draw();
    //eta_resol_EB ->Draw("same");
    //eta_resol_All->Scale(1./eta_resol_All->Integral());

    
    /*
    eta_resol_All->SetMaximum (0.08);
    eta_resol_All->SetLineWidth(2);
    eta_resol_All->GetXaxis()->SetTitleSize(0.050);
    //eta_resol_All->GetXaxis()->SetTitleOffset(0.8);
    eta_resol_All->GetYaxis()->SetTitleSize(0.050);
    //eta_resol_All->GetYaxis()->SetTitleOffset(0.8);
    eta_resol_All->SetMarkerStyle(8);
    eta_resol_All->SetMarkerSize(0.8);
    eta_resol_All->SetMarkerColor(kBlue);

    c2->SetGridx();
    c2->SetGridy();
    //eta_resol_All->Draw();
    
    TLegend* lFitEta = new TLegend (0.58, 0.60, 0.95, 0.78);
    lFitEta->SetFillStyle(0);
    lFitEta->SetBorderSize(0);
    TString textEta = "#sigma = (3.21 #pm 0.01) #upoint 10^{-2}";
    lFitEta->AddEntry ("NULL", textEta, "h");

    legend1->Draw();
    lFitEta->Draw();
    //l->Draw();
    c2->Print ("eta_resol.pdf", "pdf");
    */
// ======================================
/*
	TCanvas * c4 = new TCanvas;
	pt_resol_EB_diff ->Draw();
	pt_resol_EE_diff ->Draw("same");
	l->Draw();
	//c3->Print ("pt_resol_diff.pdf", "pdf");

    TCanvas* c5 = new TCanvas;
    //eta_resol_binned->Draw("COLZ");
    pT_resol_binned->Draw("COLZ");
*/
    
    //TCanvas* c6 = new TCanvas;
  
  // Resolution as a function of eta
  const int nEtaBinsBuf = pT_resol_binned->GetNbinsX();
  const int nEtaBins = nEtaBinsBuf;
  TH1D* projs [nEtaBins];
  
  for (int i = 0; i < nEtaBins; i++) projs[i] = pT_resol_binned->ProjectionY(Form("h_%i", i) , i+1, i+1);
  float RMS [nEtaBins];
  float RMSErr [nEtaBins];

  //TF1* fGaus = new TF1 ("fGaus", "gaus", -20, -10);
  for (int i = 0; i < nEtaBins; i++)
  {
  	float thisMean = projs[i]->GetMean();
    float thisRMS= projs[i]->GetRMS();
    float thisIntegral= projs[i]->Integral();
    //fGaus->SetParameters (thisIntegral, thisMean, thisRMS);
    //projs[i] -> Fit ("fGaus");
    //RMS[i] = fGaus->GetParameter (2);

    RMS[i] = thisRMS/thisMean;
  	RMSErr[i] = projs[i]->GetRMSError()/thisMean; // supposing error on mean is neglibible
  	//cout << i << " " << RMS[i] << " +/- " << RMSErr[i] <<  " ( " << i*2.1/nEtaBins << " - ) " << (i+1)*2.1/nEtaBins << endl;
  }

  TH1D* gr = new TH1D ("resolution_vs_eta_plot", "; |#eta(offline)|; Energy response resolution", nEtaBins, lowEta, highEta);
  for (int i = 0; i < nEtaBins; i++)
  {
    gr->SetBinContent (i+1, RMS[i]);
    //gr->SetBinContent (i+1, EffectiveRMS(projs[i]));
    gr->SetBinError (i+1, RMSErr[i]);
  }


  //projs[0]->DrawNormalized();
  //for (int i = 1; i < nEtaBins; i++) projs[i] ->DrawNormalized("same");

  //////////////////////////////////////////////
  // now do Resolution as a function of pT
  const int nPtBinsBuf = pT_resol_binned_hpsPt->GetNbinsX();
  const int nPtBins = nPtBinsBuf;
  TH1D* projs_pt [nPtBins];
  for (int i = 0; i < nPtBins; i++) projs_pt[i] = pT_resol_binned_hpsPt->ProjectionY(Form("h_pt_%i", i) , i+1, i+1);
  float RMS_pt [nPtBins];
  float RMSErr_pt [nPtBins];
  for (int i = 0; i < nPtBins; i++)
  {
    float thisMean = projs_pt[i]->GetMean();
    float thisRMS= projs_pt[i]->GetRMS();
    float thisIntegral= projs_pt[i]->Integral();
    //fGaus->SetParameters (thisIntegral, thisMean, thisRMS);
    //projs[i] -> Fit ("fGaus");
    //RMS[i] = fGaus->GetParameter (2);

    RMS_pt[i] = thisRMS/thisMean;
    RMSErr_pt[i] = projs_pt[i]->GetRMSError()/thisMean;
    cout << i << " " << RMS_pt[i] << " +/- " << RMSErr_pt[i] << endl;
  }
  

  TH1D* gr_pt = new TH1D ("resolution_vs_pT_plot", ";pT(offline); Energy response RMS", nPtBins, lowPt, highPt);
  for (int i = 0; i < nPtBins; i++)
  {
    if (RMS_pt[i] > 0) // if zero will diverge
    {
      gr_pt->SetBinContent (i+1, RMS_pt[i]);
      //gr_pt->SetBinContent (i+1, EffectiveRMS(projs[i]));
      gr_pt->SetBinError (i+1, RMSErr_pt[i]);
    }
    else
       gr_pt->SetBinContent (i+1, 0);
  } 

  /*
  gr_pt->SetMarkerStyle(20);
  gr_pt->SetMarkerSize(0.7);
  //gr_pt->SetMinimum (0.20);
  //gr_pt->SetMaximum (0.38);
  gr_pt->GetXaxis()->SetTitleSize(0.050);
  gr_pt->GetYaxis()->SetTitleSize(0.050);
  */





/*
    TCanvas* c7 = new TCanvas ("c7", "c7", 600, 600);

  TLegend* legend2 = new TLegend (0.07, 0.85, 0.9, 1.12);
    legend2->AddEntry("NULL","CMS Simulation: gg #rightarrow H #rightarrow #tau #tau         #sqrt{s}=13 TeV, bx=25ns, PU=40","h");
    //legend1->AddEntry("NULL","L1 Threshold : 30 GeV","h");
    legend2->SetLineColor(0);
    legend2->SetBorderSize(0);
    legend2->SetTextSize(0.030);
    legend2->SetFillColor(0);
    legend2->SetFillStyle(0);

	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.7);
	gr->SetMinimum (0.20);
	gr->SetMaximum (0.38);
	gr->GetXaxis()->SetTitleSize(0.050);
  gr->GetYaxis()->SetTitleSize(0.050);
  //gr->GetXaxis()->SetTitleOffset(0.9);
  //gr->GetYaxis()->SetTitleOffset(1.4);
  gr->Draw("PE");
  legend2->Draw();
	c7->Print ("RMS_vs_eta.pdf", "pdf");
 */   
/*
	TCanvas* c8 = new TCanvas ("c8", "c8", 600, 600);
	pt_resol_OneBin->DrawNormalized();
	c8->Print ("res_oneBin.pdf", "pdf");

	TCanvas* c9 = new TCanvas ("c9", "c9", 600, 600);
	TLegend* legDec = new TLegend (0.6, 0.6, 0.88, 0.88);
	legDec->SetFillColor(kWhite);
	legDec->SetLineColor(kWhite);
	legDec->AddEntry (pt_resol_decMode[0], "1 prong 0#pi0");
	legDec->AddEntry (pt_resol_decMode[1], "1 prong >0#pi0");
	legDec->AddEntry (pt_resol_decMode[2], "3 prong 0#pi0");
	legDec->AddEntry (pt_resol_decMode[3], "3 prong >0#pi0");

	pt_resol_decMode[0]->DrawNormalized();
	for (int i = 0; i < 4; i++) pt_resol_decMode[i]->DrawNormalized("same");
	legDec->Draw();
	c9->Print ("pt_res_decMode.pdf", "pdf");

	TCanvas* c10 = new TCanvas ("c10", "c10", 600, 600);
	eta_resol_decMode[0]->DrawNormalized();
	for (int i = 0; i < 4; i++) eta_resol_decMode[i]->DrawNormalized("same");
	legDec->Draw();
	c10->Print ("eta_res_decMode.pdf", "pdf");

	TCanvas* c11 = new TCanvas ("c11", "c11", 600, 600);
	phi_resol_decMode[0]->Scale(1./phi_resol_decMode[0]->Integral());
	phi_resol_decMode[0]->SetMaximum(0.08);
	phi_resol_decMode[0]->DrawNormalized();
	for (int i = 0; i < 4; i++) phi_resol_decMode[i]->DrawNormalized("same");
	legDec->Draw();
	c11->Print ("phi_res_decMode.pdf", "pdf");
*/

  fOut->Write();

}
