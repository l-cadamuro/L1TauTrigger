// evaluates rate of SINGLE TAU
// same as EvalRate.cpp, but runs on a signal that is processed with the LUT in the emulator
// so isolation is checked using only iso == 1

#include <iostream>
#include <vector>
#include <algorithm>
#include "TCanvas.h"
#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"
//#include "/home/llr/cms/cadamuro/Level1_Stage2_PUSubDevel/PUSMacros/LUTReader.h"

//#include "/home/llr/cms/cadamuro/Level1_Stage2_PUSubDevel/PUSMacros/Soglie_PerL1Talk/BinningFile.h"
//#include "/home/llr/cms/cadamuro/Level1_Stage2_PUSubDevel/PUSMacros/IsoThrCalculator.h"

using namespace std;

// checks if a bit in the flags variable is set or not
bool CheckBit (int number, int bitpos)
{
    bool res = number & (1 << bitpos);
    return res;
}


// decodes the two bit that are related to tau production (merging, etc...)
int GetTauClusterProductionMode (int flags)
{
	// check bits 30 and 31
	bool bit30 = CheckBit (flags, 30);
	bool bit31 = CheckBit (flags, 31);
	
	if (!bit31 && !bit30) return 0; // 00
	if (!bit31 && bit30)  return 1; // 01
	if (bit31 && !bit30)  return 2; // 10
	if (bit31 && bit30)   return 3; // 11
}




// this function reduces the integer containing tau clusters flags into a smaller flag
// containing only shape information of main cluster (for the moment, to be used only on non-merged clusters)
//
// Cluster with symmetrical eta shapes (one becomes the other under eta reflection) are given the same code,
// the side can always be retrieved from the TRIM_LEFT flag
//
// The output shape flag uses 7 bits.
//
// The shape computation follows the same conventions as done in the gamma gamma emulator


int ClusterShape (int flags)
{
    int shape = 0; // initialize all flags to 0
    
    bool TRIM_LEFT = CheckBit (flags, 11);
    
    if ( CheckBit (flags, 2) )                shape |= (0x1 << 0);   
    if ( CheckBit (flags, 6) )                shape |= (0x1 << 1);
    
    if ( TRIM_LEFT && CheckBit (flags, 4) )   shape |= (0x1 << 2);   
    if ( !TRIM_LEFT && CheckBit (flags, 8) )  shape |= (0x1 << 2);   
   
    if ( TRIM_LEFT && CheckBit (flags, 3) )   shape |= (0x1 << 3);   
    if ( !TRIM_LEFT && CheckBit (flags, 1) )  shape |= (0x1 << 3);   

    if ( TRIM_LEFT && CheckBit (flags, 5) )   shape |= (0x1 << 4);   
    if ( !TRIM_LEFT && CheckBit (flags, 7) )  shape |= (0x1 << 4);   

    if ( CheckBit (flags, 9) )                shape |= (0x1 << 5);   
    if ( CheckBit (flags, 10) )               shape |= (0x1 << 6);

    return shape;
    
}



// compute number of towers in the event -- uses the towers input collection
// towers are std::vector<int> objects
int computeTowers(std::vector<int>* tower_eta)
{
    int counter = 0;
    for (int tow = 0; tow < tower_eta->size(); tow++)
    {
        if ( abs(tower_eta->at(tow)) <= 4 ) // only first 4 slices in eta {-4,-3,-2,-1,+1,+2,+3,+4}
            counter++;
        
        // the check on E is useless as the tower vector is filled only if EtEm > 0 || EtHad > 0
        //if ( ( tower_hwEtEm->at(tow) + tower_hwEtHad->at(tow) ) > 0) counter++; // towers must have E > 0 
    }

    return counter;
}


void CopyArray (const int from [], int to[], int entries)
{
    for (int i = 0; i < entries; i++)
        to[i] = from[i];
        
    return;
}


bool PassVeto (int shapeNumber, const int shapeVetoes [], const int shapeVetoesNr)
{
    bool isContained = false;
    
    for (int i = 0; i < shapeVetoesNr; i++)
    {
        isContained = ( isContained || (shapeNumber == shapeVetoes[i]));
    }
    
    return (!isContained);
}


int main(int argc, char** argv)
{

    // read IsoThresholds value from command line
    if (argc < 2)
    {
        cout << "Efficiency WP name not set; please use value (e.g. 90 = 90% eff)" << endl;
        return -1;
    }

    int effWPname = atoi (argv[1]);    
    //int etaLimit = atoi (argv[2]); cout << "I'm vetoing clusters with ieta > " << etaLimit << endl;
    
    // ZeroBias sample
    //TFile* fInput = new TFile ("/data_CMS/cms/cadamuro/L1_trigger_data/MergeTrees_ZeroBias/ZeroBias_EtaOptim_PUSDevel.root"); // ZeroBias with Eta Optim   
    //TTree* tInput = (TTree*) fInput->Get("l1TriggerNtuplizer/l1TriggerTree");
    //TString fOutName = Form("rateL1Tau_%ieff_alsoShape_singleShapesVeto_eta2p1_ForEPS.root", effWPname);

    // gg fusion sample for ROC drawing
    TFile* fInput = new TFile (Form("/data_CMS/cms/cadamuro/L1_trigger_data/MergeTrees_GluGlu/MergeTree_gg_LUT%iStage2_AlwaysCBM_0.root", effWPname));
    TTree* tInput = (TTree*) fInput->Get("produceNtuple_recoTree");
    TString fOutName = Form("eff_L1Tau_Stage2_gg_LUTProcessed%i_AlwaysCMS.root", effWPname);

    // set branch and variables
    
    // TT
    std::vector<int> *TT_hweta = 0;

    // stage 2 tau
    int l1t_stage2_tau_n;
    std::vector<double>  *TauCalibPt = 0;  // calibrated tau energy
    std::vector<double>  *TauEta     = 0;  // tau eta
    std::vector<double>  *TauPhi     = 0;  // tau phi
    std::vector<int>     *TauIsoEt   = 0;  // E(5x9) - E(tau)
    std::vector<int>     *flagsVec = 0;  // flags of tau
    std::vector<int>     *TauHwPt = 0;
    std::vector<int>     *TauHwEta = 0;
    std::vector<int>     *TauHwPhi = 0;
    
           
    // branch   
    tInput ->SetBranchAddress ("l1t_stage2_tower_hwEta", &TT_hweta);    

    tInput ->SetBranchAddress ("l1t_stage2_tau_hwQual", &flagsVec);    
    tInput ->SetBranchAddress ("l1t_stage2_tau_n", &l1t_stage2_tau_n);
    tInput ->SetBranchAddress ("l1t_stage2_tau_pt", &TauCalibPt);
    tInput ->SetBranchAddress ("l1t_stage2_tau_eta", &TauEta);
    tInput ->SetBranchAddress ("l1t_stage2_tau_phi", &TauPhi);
    tInput ->SetBranchAddress ("l1t_stage2_tau_hwIso", &TauIsoEt);
    tInput ->SetBranchAddress ("l1t_stage2_tau_hwPt", &TauHwPt);
    tInput ->SetBranchAddress ("l1t_stage2_tau_hwEta", &TauHwEta);
    tInput ->SetBranchAddress ("l1t_stage2_tau_hwPhi", &TauHwPhi);

    
    // speed up with branch status
    tInput->SetBranchStatus("*", 0);
    
    tInput ->SetBranchStatus ("l1t_stage2_tower_hwEta", 1);
    
    tInput->SetBranchStatus("l1t_stage2_tau_hwQual", 1);
    tInput->SetBranchStatus("l1t_stage2_tau_n", 1);
    tInput->SetBranchStatus("l1t_stage2_tau_pt", 1);
    tInput->SetBranchStatus("l1t_stage2_tau_eta", 1);
    tInput->SetBranchStatus("l1t_stage2_tau_phi", 1);
    tInput->SetBranchStatus("l1t_stage2_tau_hwIso", 1);
    tInput->SetBranchStatus("l1t_stage2_tau_hwPt", 1);
    tInput->SetBranchStatus("l1t_stage2_tau_hwEta", 1);
    tInput->SetBranchStatus("l1t_stage2_tau_hwPhi", 1);
    
    
    ///////////////
    //// HISTO ////
    ///////////////

    TFile* fOut = new TFile (fOutName, "recreate");
   
    TH1D* LeadPt_passIso = new TH1D ("LeadPt_passIso", "LeadPt_passIso", 5000, 0, 5000);
    TH1D* SubleadPt_passIso = new TH1D ("SubleadPt_passIso", "SubleadPt_passIso", 5000, 0, 5000);
    TH1D* RelativeRate_singleTau = new TH1D ("RelativeRate_singleTau", "Relative rate - single tau", 5000, 0, 5000);
    TH1D* RelativeRate_diTau = new TH1D ("RelativeRate_diTau", "Relative rate - double tau", 5000, 0, 5000);

    TH1D* LeadPt_NoIso = new TH1D ("LeadPt_NoIso", "LeadPt_NoIso", 5000, 0, 5000);
    TH1D* SubleadPt_NoIso = new TH1D ("SubleadPt_NoIso", "SubleadPt_NoIso", 5000, 0, 5000);
    TH1D* RelativeRate_singleTau_NoIso = new TH1D ("RelativeRate_singleTau_NoIso", "Relative rate - single tau", 5000, 0, 5000);
    TH1D* RelativeRate_diTau_NoIso = new TH1D ("RelativeRate_diTau_NoIso", "Relative rate - double tau", 5000, 0, 5000);

    TH1D* LeadPt_IsoAndShape = new TH1D ("LeadPt_IsoAndShape", "LeadPt_IsoAndShape", 5000, 0, 5000);
    TH1D* SubleadPt_IsoAndShape = new TH1D ("SubleadPt_IsoAndShape", "SubleadPt_IsoAndShape", 5000, 0, 5000);
    TH1D* RelativeRate_singleTau_IsoAndShape = new TH1D ("RelativeRate_singleTau_IsoAndShape", "Relative rate - single tau", 5000, 0, 5000);
    TH1D* RelativeRate_diTau_IsoAndShape = new TH1D ("RelativeRate_diTau_IsoAndShape", "Relative rate - double tau", 5000, 0, 5000);

    
    /*
    // initialize iso calculator
  	// intialize class to compute Isolation thershold - need binning information!
	int _NetaBins = NetaBins; // these variables are globally defined in the BinningFile.h
	int _NnTTBins = NnTTBins;
    int _NEtBins = NEtBins;
    
    int _etaBins [_NetaBins+1]; 
    int _nTTBins [_NnTTBins+1];
	int _EtBins  [_NEtBins+1];
	CopyArray (etaBins, _etaBins, _NetaBins+1);
	CopyArray (nTTBins, _nTTBins, _NnTTBins+1);
	CopyArray (EtBins, _EtBins, _NEtBins+1);

	IsoThrCalculator IsoCalc (_etaBins, _nTTBins, _EtBins, _NetaBins, _NnTTBins, _NEtBins, Form("/home/llr/cms/cadamuro/Level1_Stage2_PUSubDevel/PUSMacros/Soglie_PerL1Talk/IsoEt_fits_results_%i.txt", effWPname));
    */
    
    //TString LUTfile = "/home/llr/cms/cadamuro/Level1_Stage2_PUSubDevel/PUSMacros/tau_isolation_LUT_WP"; // "new" LUT
    //TString LUTfile = "/home/llr/cms/cadamuro/Stage2_730_emulator_finalPlots_PUSdevel/CMSSW_7_3_0_pre1/src/L1Trigger/L1TCalorimeter/data/tau_isolation_LUT_WP"; //"old" LUT (from old training)
    //LUTfile += effWPname;
    //LUTfile += ".txt";
    //LUTReader lread (LUTfile.Data());
    

	// declare array of cluster shapes to be vetoed
	const int shapesVeto [] = {24, 127, 115, 107, 123, 119, 111, 49, 90, 74, 57, 16, 8};
	const int shapesVetoNr = sizeof(shapesVeto)/sizeof(int);
	
	// pritn
	cout <<  "Veto: ";
	for (int i = 0; i < shapesVetoNr; i++)
		cout << shapesVeto[i] << " ";
	cout << endl;
	
	
    // analyze data    
    long int nEvents = tInput->GetEntries(); 
    //long int nEvents = 1000000; 
    
    std::vector<double> pt_pass; // just pass iso
    std::vector<double> pt_passAndShape; // also pass veto shape
    std::vector<double> pt_noiso;
    
    // loop on all events
    for (long int iEv = 0; iEv < nEvents; iEv++)
    {
        tInput->GetEntry(iEv);
        if (iEv%100000 == 0) cout << iEv << " / " << nEvents << endl;
        
        // clear pt vector for this event and compute number of TT
        pt_pass.clear();
        pt_noiso.clear();
        pt_passAndShape.clear();
        
        //int nTT = computeTowers (TT_hweta);
        
        // loop on all L1 taus --> save all taus passing Iso requirement + other selections
        for (int iL1 = 0; iL1 < l1t_stage2_tau_n; iL1++)
        {
            // selections
            double abseta = TMath::Abs( (*TauEta)[iL1] );
            int  abshweta = abs ( (*TauHwEta)[iL1] );
            double tauPt  = (*TauCalibPt)[iL1];
            int tauHwIso  = (*TauIsoEt)[iL1];
            int tauHwPt   = (*TauHwPt)[iL1];
            
            if (abseta < 2.2 && abshweta < 28)
            //if (abshweta <= etaLimit) 
            //if (abshweta <= 26)            
            {
                pt_noiso.push_back (tauPt);
            
                // apply iso
                //int thresh = lread.Threshold(abshweta, tauHwPt, nTT);
                //bool PassIso = ( tauHwIso < thresh );
                bool PassIso = ( tauHwIso  == 1 );
                //cout << "hweta - nTT - hwPt: " << abshweta << " - " << nTT << " - " << tauHwPt << " | IsoEt - thresh: " << tauHwIso << " - " << thresh << " | Pass: " << PassIso << endl;
                if (PassIso)
                {
                    pt_pass.push_back (tauPt);
                    
                    // now test shape!
                    int flags = (*flagsVec)[iL1];
                    int clusterShape = ClusterShape(flags);                    
                    // selection on shapes --  still to be optimized!!
                    //if (clusterShape < 99)
                    if (PassVeto(clusterShape, shapesVeto, shapesVetoNr))
					{
                        pt_passAndShape.push_back (tauPt);
                    }
                }
            }
        }

    
        // now that all taus in the event are analyzed, fill the histograms of lead and sublead pt (if missing, fil with a -1)
        std::sort (pt_pass.begin(), pt_pass.end());
        std::sort (pt_noiso.begin(), pt_noiso.end());
        std::sort (pt_passAndShape.begin(), pt_passAndShape.end());
        
        // no iso operations
        if (pt_noiso.size() >= 2 )
        {
            LeadPt_NoIso -> Fill ( *(pt_noiso.rbegin()) ); 
            SubleadPt_NoIso -> Fill ( *(pt_noiso.rbegin() + 1) );
        } 
        
        else if (pt_noiso.size() == 1) 
        {
            LeadPt_NoIso -> Fill ( *(pt_noiso.rbegin()) ); 
            SubleadPt_NoIso -> Fill ( -1 );
        } 
        
        else
        {
            LeadPt_NoIso -> Fill ( -1 ); 
            SubleadPt_NoIso -> Fill ( -1 );        
        }       


        
        // pass iso operations
        if (pt_pass.size() >= 2 )
        {
            LeadPt_passIso -> Fill ( *(pt_pass.rbegin()) ); 
            SubleadPt_passIso -> Fill ( *(pt_pass.rbegin() + 1) );
        } 
        
        else if (pt_pass.size() == 1) 
        {
            LeadPt_passIso -> Fill ( *(pt_pass.rbegin()) ); 
            SubleadPt_passIso -> Fill ( -1 );
        } 
        
        else
        {
            LeadPt_passIso -> Fill ( -1 ); 
            SubleadPt_passIso -> Fill ( -1 );        
        }      
        
        
        // also shape operations
        if (pt_passAndShape.size() >= 2 )
        {
            LeadPt_IsoAndShape -> Fill ( *(pt_passAndShape.rbegin()) ); 
            SubleadPt_IsoAndShape -> Fill ( *(pt_passAndShape.rbegin() + 1) );
        } 
        
        else if (pt_passAndShape.size() == 1) 
        {
            LeadPt_IsoAndShape -> Fill ( *(pt_passAndShape.rbegin()) ); 
            SubleadPt_IsoAndShape -> Fill ( -1 );
        } 
        
        else
        {
            LeadPt_IsoAndShape -> Fill ( -1 ); 
            SubleadPt_IsoAndShape -> Fill ( -1 );        
        }     
    }
    
    // compute rate plots
   cout << "Computing rates..." << endl; 
    
    for (int i = 1; i <= 5000; i++)
    {
        // no Iso
        int Tot_singleTau_noiso = LeadPt_NoIso -> Integral (0, 5001);
        int Tot_diTau_noiso = SubleadPt_NoIso -> Integral (0, 5001);
        
        double relRateSingle_noiso = 1.*(LeadPt_NoIso->Integral(i, 5001))/Tot_singleTau_noiso;
        double relRateDouble_noiso = 1.*(SubleadPt_NoIso->Integral(i, 5001))/Tot_diTau_noiso;
        
        RelativeRate_singleTau_NoIso -> SetBinContent (i, relRateSingle_noiso);
        RelativeRate_diTau_NoIso -> SetBinContent (i, relRateDouble_noiso);

        // with Iso
        int Tot_singleTau = LeadPt_passIso -> Integral (0, 5001);
        int Tot_diTau = SubleadPt_passIso -> Integral (0, 5001);
        
        double relRateSingle = 1.*(LeadPt_passIso->Integral(i, 5001))/Tot_singleTau;
        double relRateDouble = 1.*(SubleadPt_passIso->Integral(i, 5001))/Tot_diTau;
        
        RelativeRate_singleTau -> SetBinContent (i, relRateSingle);
        RelativeRate_diTau -> SetBinContent (i, relRateDouble);        

        // with Iso and shape
        int Tot_singleTau_shape = LeadPt_IsoAndShape -> Integral (0, 5001);
        int Tot_diTau_shape = SubleadPt_IsoAndShape -> Integral (0, 5001);
        
        double relRateSingle_shape = 1.*(LeadPt_IsoAndShape->Integral(i, 5001))/Tot_singleTau_shape;
        double relRateDouble_shape = 1.*(SubleadPt_IsoAndShape->Integral(i, 5001))/Tot_diTau_shape;
        
        RelativeRate_singleTau_IsoAndShape -> SetBinContent (i, relRateSingle_shape);
        RelativeRate_diTau_IsoAndShape -> SetBinContent (i, relRateDouble_shape);        

    }
    
    fOut -> Write();
}

