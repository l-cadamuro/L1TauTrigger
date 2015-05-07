// evaluates rate  - compile with c++ -lm -o EvalRate_OnL1Legacy EvalRate_OnL1Legacy.cpp `root-config --glibs --cflags`
// and launch the executable 

#include <iostream>
#include <vector>
#include <algorithm>
#include "TCanvas.h"
#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"

using namespace std;


int main(int argc, char** argv)
{
    //double scaleRunI = 1./1.713;
    double scaleRunI = 1.075/1.713;
    bool doScaleRunI = true;

    cout << "Do scale: " << doScaleRunI << endl;
    cout << "Scale factor: " << scaleRunI << endl;

    // ZeroBias sample L1 (processed by Luca M)
    TChain * tInput = new TChain ("recoTree");

    tInput -> Add ("/data_CMS/cms/cadamuro/L1_trigger_data/RunI_Legacy_trigger_trees/ZeroBias/ZeroBias_RunI_190Files_Pruned.root");
    TString fOutName = "rateL1Tau_L1Legacy_ZeroBiasPU50bx25E13TeV_forEPS";
    if (doScaleRunI) fOutName += "RunIScaled.root";
    else fOutName += "RunINotScaled.root";


    // set branch and variables
 
    int _trig_L1tau_N;
    std::vector<double>  * trig_L1tau_eta  = 0;
    std::vector<double>  * trig_L1tau_phi  = 0;
    //std::vector<double>  * trig_L1tau_energy = 0;
    std::vector<double>  * trig_L1tau_et   = 0;
 
    tInput ->SetBranchAddress ("trig_L1tau_N", &_trig_L1tau_N);
    tInput ->SetBranchAddress ("trig_L1tau_eta", &trig_L1tau_eta);
    tInput ->SetBranchAddress ("trig_L1tau_phi", &trig_L1tau_phi);
    //tInput ->SetBranchAddress ("trig_L1tau_energy", &trig_L1tau_energy);
    tInput ->SetBranchAddress ("trig_L1tau_et", &trig_L1tau_et);

    // speed up with branch status
    tInput->SetBranchStatus("*", 0);

    tInput ->SetBranchStatus ("trig_L1tau_N", 1);
    tInput ->SetBranchStatus ("trig_L1tau_eta", 1);
    tInput ->SetBranchStatus ("trig_L1tau_phi", 1);
    //tInput ->SetBranchStatus ("trig_L1tau_energy", 1);
    tInput ->SetBranchStatus ("trig_L1tau_et", 1);
    
    ///////////////
    //// HISTO ////
    ///////////////

    TFile* fOut = new TFile (fOutName, "recreate");
   
    // the tree histograms are filed in identical way, structure is simply copied from the Stage 2 rate Evaluator
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

    
    // analyze data    
    long int nEvents = tInput->GetEntries(); 
    //long int nEvents = 100000; 
    
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
        
        
        // loop on all L1 taus --> save all taus passing Iso requirement + other selections
        //cout << iEv << " " << _trig_L1tau_N << " " << trig_L1tau_eta->size() << endl;
        for (int iL1 = 0; iL1 < _trig_L1tau_N; iL1++)
        {
            // selections
            double abseta = TMath::Abs( (*trig_L1tau_eta)[iL1] );
            double tauPt  = (*trig_L1tau_et)[iL1];
            
            if (doScaleRunI) tauPt = scaleRunI*tauPt;

            if (abseta < 2.2)
            {
                pt_noiso.push_back (tauPt);
            
                // apply iso
                bool PassIso = true; // iso is applied by default in the collection
                if (PassIso)
                {
                    pt_pass.push_back (tauPt);
                    
                    
                    if (true) // there is no shape to test in Legacy trigger
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

