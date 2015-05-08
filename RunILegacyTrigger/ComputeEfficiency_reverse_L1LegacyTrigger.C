// this program does the same as ProduceFilteredTaus.C, but the logic of the matching is reversed:
//
// 1) the main loop is on all gen taus
// 2) for each gen, find the closest hps in dR (hps must also pass analysis selections)
// 3) once the hps has been found, save the closest matching L1
//
// in this way, a triplet is formed with a (gen, hps, L1) and it is saved as a "filtered tau"
//
// v2: uses as input the new AOD files with my corrections (max 200 taus pt ordered, was not like this before)
// v3: added tau charge
// v4: input mergeTree has also gen tau energy, and it is saved in the output, plus it has SkipIso flag from command line
// v5: input mergeTree has been created with new hps version that contains hps mass
// v6: nTT variable is computed using the TT collection available; the tau_hwQual vector now contains the flags related to tau algo
// v7: added hardware quantities in the output (to be used in PUS development)
// v8: added HPSIsoflag value (both 0 and 1 values will only be saved if skipHPSIso == true, else will only save those with 1)
// cmdLine: take infile and output path from cmd line 
// v9: add possibility to resolve ambiguities by highest pT and not closest dR (try to explain endcap / barrel angular resolution difference)

//
// computes the efficiency as a function of the threshold for ROC curve (rate vs efficiency)
//

#include "TString.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TTree.h"
#include "TFile.h"
#include "TProfile.h"
#include "TMath.h"
#include "TStyle.h"
#include <iostream>
#include <algorithm>

#define maxNHPS 200 // set it to 5 if using OLD hps file --> should not really matter if is bigger...

using namespace std;

// deltaR
float deltaR(float eta1, float phi1, float eta2, float phi2)
{
  float deta = eta1 - eta2 ;
  float dphi = phi1 - phi2 ;
  while (dphi > TMath::Pi() ) 
    dphi -= 2*TMath::Pi();
  while (dphi <= -TMath::Pi() ) 
    dphi += 2*TMath::Pi();
  return sqrt(deta*deta + dphi*dphi); 
}



int main(int argc, char** argv)
{
    bool SkipHPSIso = false;
    //int effWPName = 90;
    //bool AskL1Iso = true;
 
    if (argc < 2)
    {
        cout << "Usage: [skipHPSIso]" << endl;
        cout << "SkipHPSIso not inserted, using : " << SkipHPSIso << endl;
        //cout << "effWPName = " << effWPName << endl;
    }
    
    else
    {
        if (atof(argv[1]) != 0)
        {
            cout << "SkipHPSIso: true" << endl;
            SkipHPSIso = true;
        }
        else
        {
            cout << "SkipHPSIso: false" << endl;
            SkipHPSIso = false;
        }
    }
    //////////////////////
    // input settings
    //////////////////////
    
    TString inFile = "/data_CMS/cms/cadamuro/L1_trigger_data/RunI_Legacy_trigger_trees/ggFusion/RunI_LegacyL1_ggFusion_bx25Pu40E13TeV.root";
    TString outFile = "eff_L1Tau_RunI_gg_";
    if (SkipHPSIso) outFile += "NoHpsIso.root";
    else outFile += "WithHpsIso.root";
   
    TFile* fInput = new TFile (inFile);
    TTree* tInput = (TTree*) fInput->Get("produceNtuple/recoTree");
    
    // configuration - ambiguities handling
    const bool resolveL1Ambiguities = true; // if false: take the last "good" L1 (like Luca M was doing)
    const bool resolveL1bydR = false; // true: min dR -- false: max pT
    
    cout << "Resolve L1-HPS match ambiguities: " << resolveL1Ambiguities << endl;
    cout << "Resolve L1-HPS ambiguities by dR: " <<  resolveL1bydR << endl;
    
    double scaleRunI = 1.075/1.713;
    cout << "SCALING RUN I" << endl;

    // set branch and variables
    
   
    // stage 2 tau
    int trig_L1tau_N;
    std::vector<double>  *TauCalibPt = 0;  // calibrated tau energy
    std::vector<double>  *TauEta     = 0;  // tau eta
    std::vector<double>  *TauPhi     = 0;  // tau phi
   
    
    
    // hps tau
    float HpsPt[maxNHPS];   // hps pt
    float HpsEta[maxNHPS];  // hps eta
    float HpsPhi[maxNHPS];  // hps phi
    float HpsMass[maxNHPS];  // hps mass
    float HpsTauDiscrByDecMode[maxNHPS];
    float HpsTauDiscrByMediumIso[maxNHPS];
    int   HpsCharge[maxNHPS]; // tau charge
    int   nHps;
        
    // gen
    float GenPt[5];   // gen pt
    float GenEta[5];  // gen eta
    float GenPhi[5];  // gen phi
    float GenEnergy[5]; // gen energy
    int   GenDecayMode[5]; // gen decay mode
    int   nGen;
        
    // branch   

    tInput ->SetBranchAddress ("trig_L1tau_N", &trig_L1tau_N);
    tInput ->SetBranchAddress ("trig_L1tau_et", &TauCalibPt);
    tInput ->SetBranchAddress ("trig_L1tau_eta", &TauEta);
    tInput ->SetBranchAddress ("trig_L1tau_phi", &TauPhi);
  

    tInput ->SetBranchAddress ("hpsTau_pt", &HpsPt);
    tInput ->SetBranchAddress ("hpsTau_eta", &HpsEta);
    tInput ->SetBranchAddress ("hpsTau_phi", &HpsPhi);
    tInput ->SetBranchAddress ("hpsTau_mass", &HpsMass);
    tInput ->SetBranchAddress ("hpsTau_decayMode", &HpsTauDiscrByDecMode);
    tInput ->SetBranchAddress ("hpsTau_isoM", &HpsTauDiscrByMediumIso);
    tInput ->SetBranchAddress ("hpsTau_N", &nHps);
    tInput ->SetBranchAddress ("hpsTau_charge", &HpsCharge);

    tInput ->SetBranchAddress ("genTau_pt", &GenPt);
    tInput ->SetBranchAddress ("genTau_eta", &GenEta);
    tInput ->SetBranchAddress ("genTau_phi", &GenPhi);
    tInput ->SetBranchAddress ("genTau_energy", &GenEnergy);
    tInput ->SetBranchAddress ("genTau_DecayMode", &GenDecayMode);
    tInput ->SetBranchAddress ("genTau_N", &nGen);


    // speed up with branch status
    tInput->SetBranchStatus("*", 0);
    
    tInput->SetBranchStatus("trig_L1tau_N", 1);
    tInput->SetBranchStatus("trig_L1tau_et", 1);
    tInput->SetBranchStatus("trig_L1tau_eta", 1);
    tInput->SetBranchStatus("trig_L1tau_phi", 1);
     
    tInput->SetBranchStatus("hpsTau_pt", 1);
    tInput->SetBranchStatus("hpsTau_eta", 1);
    tInput->SetBranchStatus("hpsTau_phi", 1);
    tInput->SetBranchStatus("hpsTau_mass", 1);
    tInput->SetBranchStatus("hpsTau_decayMode", 1);
    tInput->SetBranchStatus("hpsTau_isoM", 1);
    tInput->SetBranchStatus("hpsTau_N", 1);
    tInput->SetBranchStatus("hpsTau_charge", 1);
    
    tInput->SetBranchStatus("genTau_pt", 1);
    tInput->SetBranchStatus("genTau_eta", 1);
    tInput->SetBranchStatus("genTau_phi", 1);
    tInput->SetBranchStatus("genTau_energy", 1);
    tInput->SetBranchStatus("genTau_DecayMode", 1);
    tInput->SetBranchStatus("genTau_N", 1);



    ///////////////////////
    // output settings
    ///////////////////////
    
    TFile* fOutput = new TFile (outFile, "recreate");
    //histo definition
    TH1D*h_LeadGenPt = new TH1D("h_LeadGenPt","h_LeadGenPt",300,0,300);
    TH1D*h_LeadGenEta = new TH1D("h_LeadGenEta","h_LeadGenEta",100,-10,10);
    TH1D*h_SubLeadGenPt = new TH1D("h_SubLeadGenPt","h_SubLeadGenPt",300,0,300);
    TH1D*h_SubLeadGenEta = new TH1D("h_SubLeadGenEta","h_SubLeadGenEta",100,-10,10);
   
    TH1D*h_LeadTauPt = new TH1D("h_LeadTauPt","h_LeadTauPt",300,0,300);
    TH1D*h_SubLeadTauPt = new TH1D("h_SubLeadTauPt","h_SubLeadTauPt",300,0,300);
    TH1D*h_IsoEt = new TH1D("h_IsoEt","h_IsoEt",100,0,100);
    TH2D*h_Scatter_IsoEt_Vs_Pt = new TH2D("h_Scatter_IsoEt_Vs_Pt","h_Scatter_IsoEt_Vs_Pt",300,0,300,100,0,100);
    TH1D*h_Efficiency = new TH1D("h_Efficiency","h_Efficiency",300,0,300);
    TH1D*hTauInEventAfterSelection = new TH1D("hTauInEventAfterSelection","hTauInEventAfterSelection",10,0,10);

    int numEntries = tInput->GetEntries() ;
    int nProcess = numEntries;      
 /*
    TTree* tOutput = new TTree ("filtered_taus_tree", "filtered_taus_tree");
    
    // variables
    float L1_pt;
    float L1_eta;
    float L1_phi;
    int   L1_flags;
    int EvtNtt, nTT;
    int   L1_IsoEt;
    int L1_hwpt;
    int L1_hweta;
    int L1_hwphi;
    float hps_pt;
    float hps_eta;
    float hps_phi;
    float hps_mass;
    int hps_charge;
    int hps_iso;
    float gen_pt;
    float gen_eta;
    float gen_phi;
    float gen_energy;
    int gen_decaymode;
    
    tOutput -> Branch ("L1_pt", &L1_pt, "L1_pt/F");
    tOutput -> Branch ("L1_eta", &L1_eta, "L1_eta/F");
    tOutput -> Branch ("L1_phi", &L1_phi, "L1_phi/F");
    tOutput -> Branch ("L1_flags", &L1_flags, "L1_flags/I");
    tOutput -> Branch ("L1_hwpt", &L1_hwpt, "L1_hwpt/I");
    tOutput -> Branch ("L1_hweta", &L1_hweta, "L1_hweta/I");
    tOutput -> Branch ("L1_hwphi", &L1_hwphi, "L1_hwphi/I");
    tOutput -> Branch ("EvtNtt", &EvtNtt, "EvtNtt/I");
    tOutput -> Branch ("L1_IsoEt", &L1_IsoEt, "L1_IsoEt/I");
    tOutput -> Branch ("hps_pt", &hps_pt, "hps_pt/F");
    tOutput -> Branch ("hps_eta", &hps_eta, "hps_eta/F");
    tOutput -> Branch ("hps_phi", &hps_phi, "hps_phi/F");
    tOutput -> Branch ("hps_mass", &hps_mass, "hps_mass/F");
    tOutput -> Branch ("hps_charge", &hps_charge, "hps_charge/I");
    tOutput -> Branch ("hps_iso", &hps_iso, "hps_iso/I");
    tOutput -> Branch ("gen_pt", &gen_pt, "gen_pt/F");
    tOutput -> Branch ("gen_eta", &gen_eta, "gen_eta/F");
    tOutput -> Branch ("gen_phi", &gen_phi, "gen_phi/F");
    tOutput -> Branch ("gen_energy", &gen_energy, "gen_energy/F");
    tOutput -> Branch ("gen_decaymode", &gen_decaymode, "gen_decaymode/I");
*/

    /////////////////////
    // matching phase
    /////////////////////

    int nEvents = tInput->GetEntries(); 
    //int nEvents = 1000; 
    
    int nTausIn = 0;  // total number of gen Taus processed (main loop)
    int nTausOut = 0; // total number of triplets (gen, hps, L1) saved

    // these 3 variables will save the address of the triplet elements
    int genIndex = -1;
    int hpsIndex = -1;
    int L1Index = -1;
    
    // these vector will contain particle candidates indexes and dRs
    vector<int> v_hpsIndex;
    vector<float> v_hpsdR;
    vector<int> v_L1Index;
    vector<float> v_L1dR; // use to resolve ambiguities with smallest dR
    vector<float> v_L1pT; // use to resolve ambiguities with highest pT
    
    // for histos setting
    double LeadingTau_Pt = 0;
    double SubLeadingTau_Pt = 0;
    int EffDenom = 0;
    int NhpsGoodTaus = 0;
    std::vector<double> _GoodL1TauPt;

    // loop on all events
    for (int iEv = 0; iEv < nEvents; iEv++)
    {
        tInput->GetEntry(iEv);
        if (iEv%10000 == 0) cout << iEv << endl;
        
        // compute evt nTT (global property of the event)
        //nTT = computeTowers(TT_hweta);
        
        // reset counters for this event
        NhpsGoodTaus = 0;
        _GoodL1TauPt.clear();

        // loop on all gen taus
        for (int iGen = 0; iGen < nGen; iGen++)
        {
            // reset triplet indexes, initializing genIndex to current gen particle index
            genIndex = iGen;
            hpsIndex = -1;
            L1Index = -1;
            
            v_hpsIndex.clear();
            v_hpsdR.clear();
            v_L1Index.clear();
            v_L1dR.clear();
            v_L1pT.clear();
            
            nTausIn++;
            
            /////// HPS match ////////
            // search for match with the hps candidate 
            for (int iHps = 0; iHps < nHps; iHps++)
            {
                bool analysisSel = false;
                analysisSel = (HpsPt[iHps] > 20. && TMath::Abs(HpsEta[iHps]) < 2.5 && HpsTauDiscrByDecMode[iHps] ==1 && (HpsTauDiscrByMediumIso[iHps] ==1 || SkipHPSIso) );
                float dRbuf = deltaR (GenEta[genIndex], GenPhi[genIndex], HpsEta[iHps], HpsPhi[iHps]);
                if (analysisSel && dRbuf < 0.1)
                {
                    v_hpsIndex.push_back(iHps);
                    v_hpsdR.push_back(dRbuf);
                }
            }
            
            // now find hps index corresponding to minimum dR
            vector<float>::iterator mindr_hps_ptr = min_element (begin(v_hpsdR), end(v_hpsdR));
            if (mindr_hps_ptr != end (v_hpsdR)) // check if vector is non empty
            {
                int hps_index_position = distance (begin (v_hpsdR), mindr_hps_ptr);
                hpsIndex = v_hpsIndex.at(hps_index_position);
            }
            
            /////// L1 match ////////
            // continue with L1 matching only if hps match was positive
            if (hpsIndex > -1)
            {
                // loop on all L1 candidates
                for (int iL1 = 0; iL1 < trig_L1tau_N; iL1++)
                {
                    //bool L1Isolato = ( (*TauIsoEt)[iL1] == 1); // run I always isolated...
                    float dRbuf2 = deltaR (HpsEta[hpsIndex], HpsPhi[hpsIndex], (*TauEta)[iL1], (*TauPhi)[iL1]);
                    if (dRbuf2 < 0.5 )
                    {
                        v_L1Index.push_back(iL1);
                        v_L1dR.push_back(dRbuf2);
                        v_L1pT.push_back( (*TauCalibPt)[iL1] );
                    } 
                }
                   
                // resolve ambiguities                
                if (resolveL1Ambiguities)
                {            
                    if (resolveL1bydR) // use dR
                    {   
                        // search for the L1 index corresponding to minimum dR
                        vector<float>::iterator mindr_L1_ptr = min_element(begin(v_L1dR), end(v_L1dR));
                        if (mindr_L1_ptr != end(v_L1dR)) // check if vector is empty
                        {
                            int L1_index_position = distance (begin(v_L1dR), mindr_L1_ptr);
                            L1Index = v_L1Index.at(L1_index_position);
                        }
                    }
                    
                    else // use max pT
                    {
                        vector<float>::iterator maxpT_L1_ptr = max_element(begin(v_L1pT), end(v_L1pT));
                        if (maxpT_L1_ptr != end(v_L1pT)) // check if vector is empty
                        {
                            int L1_index_position = distance (begin(v_L1pT), maxpT_L1_ptr);
                            L1Index = v_L1Index.at(L1_index_position);
                        }

                    
                    }
                    
                }
                
                else // do not resolve ambiguities, take last
                {
                    L1Index = v_L1Index.back();
                }                
                
            }

            if (hpsIndex > -1) NhpsGoodTaus++; 
            if (L1Index > -1) _GoodL1TauPt.push_back(scaleRunI * ((*TauCalibPt)[L1Index])); // NNED TO SCALE RUN I

            /*
            // save output only if triplet is complete
            if (genIndex > -1 && hpsIndex > -1 && L1Index > -1)
            {
                gen_pt = GenPt [genIndex];
                gen_eta = GenEta [genIndex];
                gen_phi = GenPhi [genIndex];
                gen_energy = GenEnergy [genIndex];
                gen_decaymode = GenDecayMode [genIndex];
            
                hps_pt = HpsPt[hpsIndex];
                hps_eta = HpsEta[hpsIndex];
                hps_phi = HpsPhi[hpsIndex];
                hps_mass = HpsMass[hpsIndex];
                hps_charge = HpsCharge[hpsIndex];
                hps_iso = HpsTauDiscrByMediumIso[hpsIndex];
                
                L1_pt = (*TauCalibPt)[L1Index];
                L1_eta = (*TauEta)[L1Index];
                L1_phi = (*TauPhi)[L1Index];
                L1_flags = (*flagsVec)[L1Index];
                L1_IsoEt = (*TauIsoEt)[L1Index];
                L1_hwpt = (*TauHwPt)[L1Index];
                L1_hweta = (*TauHwEta)[L1Index];
                L1_hwphi = (*TauHwPhi)[L1Index];
                
                EvtNtt = nTT;

                tOutput->Fill();
                nTausOut++;
            }
            */
        }

        // event finished, compute efficiencies
        if (NhpsGoodTaus >= 2)
        {
            EffDenom++;
                
            if (_GoodL1TauPt.size() < 2)
            SubLeadingTau_Pt = -1.;
        
            else
            {
                std::sort(_GoodL1TauPt.begin(),_GoodL1TauPt.end());
                //SubLeadingTau_Pt = (*(_GoodL1TauPt.rbegin()+1))*k_scale;
                SubLeadingTau_Pt = (*(_GoodL1TauPt.rbegin()+1));
            }
            
            h_SubLeadTauPt->Fill(SubLeadingTau_Pt);
        } 
        if (NhpsGoodTaus > 2) cout << "Warning: there are " << NhpsGoodTaus << " good hps taus! (more than 2)" << endl;
    }
    
    // loop on events finished, write tree and print statistics
    cout << "Finished, selected " << nTausOut << " triplets on " << nTausIn << " input gen level taus" << endl;
    cout << "SkipHPSIso value for this selection was: " << SkipHPSIso << endl;
    
    // compute effs
    for(int i=0;i<300;i++){
        //h_Efficiency->SetBinContent(i+1,(float)(h_SubLeadTauPt->Integral(i+1,300))/(1.*EffDenom));
    h_Efficiency->SetBinContent(i+1,(float)(h_SubLeadTauPt->Integral(i+1,300))/(h_SubLeadTauPt->Integral(0,300)));
        std::cout << "numerator: " << (h_SubLeadTauPt->Integral(i+1,300)) << "   denominator: " << EffDenom << std::endl;
    }

    // Save histo
    h_LeadGenPt->Write();
    h_LeadGenEta->Write();
    h_SubLeadGenPt->Write(); 
    h_SubLeadGenEta->Write();
    h_LeadTauPt->Write();
    h_SubLeadTauPt->Write();
    h_Efficiency->Write();
    h_IsoEt->Write();
    h_Scatter_IsoEt_Vs_Pt->Write();  
    hTauInEventAfterSelection->Write(); 
    
    //myChain->Delete();


    //tOutput->Write();
    //fOutput->Close();
    return 1;
}
