// output is a root file containing for each entry a "triplet" of L1 - hps - gen candidates
//
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
// v10: adding npu vertexes and npv

// produce the filtered taus trees on Legacy RunI trigger
// 

#include "TString.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TTree.h"
#include "TFile.h"
#include "TProfile.h"
#include "TMath.h"
#include "TStyle.h"
#include"TGraphAsymmErrors.h"
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
    
    if (argc < 2)
        cout << "SkipHPSIso not inserted from cmd line ==> using " << SkipHPSIso << endl;
    
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
    
    // gg fusion
    TFile* fInput = new TFile ("/data_CMS/cms/cadamuro/L1_trigger_data/RunI_Legacy_trigger_trees/ggFusion/RunI_LegacyL1_ggFusion_bx25Pu40E13TeV.root");
    TString outFile = "filtered_taus_L1LegacyRunI_gg_bx25Pu40E13TeV";
    //TString outFile = "MyTestOnlyPi0Dr1p0SolveByPt";
    if (SkipHPSIso) outFile += "NoHpsIso.root";
    else outFile += "WithHpsIso.root";

    TTree* tInput = (TTree*) fInput->Get("produceNtuple/recoTree");
    
    // configuration - ambiguities handling
    const bool resolveL1Ambiguities = true; // if false: take the last "good" L1 (like Luca M was doing)
    const bool resolveL1bydR = false; // true: min dR -- false: max pT
    
    cout << "Resolve L1-HPS match ambiguities: " << resolveL1Ambiguities << endl;
    cout << "Resolve L1-HPS ambiguities by dR: " <<  resolveL1bydR << endl;
    
    // set branch and variables
    
    // TT
    int in_npu;
    int in_npv;

    int trig_L1tau_N;
    std::vector<double>  *TauCalibPt = 0;  // calibrated tau energy (is an Et, not really a Pt)
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
    tInput ->SetBranchAddress ("nPU", &in_npu);
    tInput ->SetBranchAddress ("nVtx", &in_npv);

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

    tInput->SetBranchStatus("nPU", 1);
    tInput->SetBranchStatus("nVtx", 1);


    ///////////////////////
    // output settings
    ///////////////////////
    
    TFile* fOutput = new TFile (outFile, "recreate");
/*
    if (SkipHPSIso)
        fOutput = new TFile ("VBF_filtered_taus_NoHPSIsorequirement_ReverseMatch_v8.root", "recreate");
    else
        fOutput = new TFile ("VBF_filtered_taus_WithHPSIsorequirement_ReverseMatch_v8.root", "recreate");
*/        
    TTree* tOutput = new TTree ("filtered_taus_tree", "filtered_taus_tree");
    
    // variables
    float L1_pt;
    float L1_eta;
    float L1_phi;

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
    int npu;
    int npv;

    tOutput -> Branch ("L1_pt", &L1_pt, "L1_pt/F");
    tOutput -> Branch ("L1_eta", &L1_eta, "L1_eta/F");
    tOutput -> Branch ("L1_phi", &L1_phi, "L1_phi/F");

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
    tOutput -> Branch ("npu", &npu, "npu/I");
    tOutput -> Branch ("npv", &npv, "npv/I");

    /////////////////////
    // matching phase
    /////////////////////

    // histos for the matching efficiency as a function of pt, eta, phi
    // 0: hps-gen PRE MATCH  
    // 1: hps-gen POST MATCH  
    // 2: hps-L1 PRE MATCH  
    // 3: hps-L1 PRE MATCH  
    TH1D* matchEffEta [4];
    TH1D* matchEffPhi [4];
    TH1D* matchEffPt [4];
    for (int i = 0; i < 4; i++)
    {
        matchEffEta[i] = new TH1D (Form("match_eff_eta_%i", i), Form("match_eff_eta_%i;hps or gen eta; N", i), 150, -3., 3.);
        matchEffPhi[i] = new TH1D (Form("match_eff_phi_%i", i), Form("match_eff_phi_%i;hps or gen phi; N", i), 150, -1.*TMath::Pi(), TMath::Pi());
        matchEffPt[i] = new TH1D (Form("match_eff_pt_%i", i), Form("match_eff_pt_%i;hps or gen Pt; N", i), 150, 0, 150);
    }

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
    
    int multHpsGen = 0;
    int multHpsL1 = 0;

    // loop on all events
    for (int iEv = 0; iEv < nEvents; iEv++)
    {
        tInput->GetEntry(iEv);
        if (iEv%10000 == 0) cout << iEv << endl;
        
        // loop on all gen taus
        for (int iGen = 0; iGen < nGen; iGen++)
        {
            //if (GenDecayMode[iGen] != 0) continue; // JUST FOR TEST, filter on only one progs

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
            
            matchEffEta[0]->Fill (GenEta[genIndex]);
            matchEffPhi[0]->Fill (GenPhi[genIndex]);
            matchEffPt[0]->Fill (GenPt[genIndex]);
            /////// HPS match ////////
            // search for match with the hps candidate 
            for (int iHps = 0; iHps < nHps; iHps++)
            {
                bool analysisSel = false;
                analysisSel = (HpsPt[iHps] > 20. && TMath::Abs(HpsEta[iHps]) < 2.1 && HpsTauDiscrByDecMode[iHps] ==1 && (HpsTauDiscrByMediumIso[iHps] ==1 || SkipHPSIso) );
                float dRbuf = deltaR (GenEta[genIndex], GenPhi[genIndex], HpsEta[iHps], HpsPhi[iHps]);
                if (analysisSel && dRbuf < 0.1)
                {
                    v_hpsIndex.push_back(iHps);
                    v_hpsdR.push_back(dRbuf);
                }
            }
            
            if (v_hpsdR.size()>1) multHpsGen++;

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
                // update eff counters
                matchEffEta[1]->Fill (GenEta[genIndex]);
                matchEffPhi[1]->Fill (GenPhi[genIndex]);
                matchEffPt[1]->Fill (GenPt[genIndex]);

                matchEffEta[2]->Fill (HpsEta[hpsIndex]);
                matchEffPhi[2]->Fill (HpsPhi[hpsIndex]);
                matchEffPt[2]->Fill (HpsPt[hpsIndex]);

                // loop on all L1 candidates
                for (int iL1 = 0; iL1 < trig_L1tau_N; iL1++)
                {
                   
                    float dRbuf2 = deltaR (HpsEta[hpsIndex], HpsPhi[hpsIndex], (*TauEta)[iL1], (*TauPhi)[iL1]);
                    if (dRbuf2 < 0.5)
                    {
                        v_L1Index.push_back(iL1);
                        v_L1dR.push_back(dRbuf2);
                        v_L1pT.push_back( (*TauCalibPt)[iL1] );
                    } 
                }
                
                if (v_L1dR.size()>1) multHpsL1++;


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
            
            // save output only if triplet is complete
            if (genIndex > -1 && hpsIndex > -1 && L1Index > -1)
            {
                // update eff counters
                matchEffEta[3]->Fill (HpsEta[hpsIndex]);
                matchEffPhi[3]->Fill (HpsPhi[hpsIndex]);
                matchEffPt[3]->Fill (HpsPt[hpsIndex]);


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

                npu = in_npu;
                npv = in_npv;

                tOutput->Fill();
                nTausOut++;
            }
        }
    }
    
    // loop on events finished, write tree and print statistics
    cout << "Finished, selected " << nTausOut << " triplets on " << nTausIn << " input gen level taus" << endl;
    cout << "SkipHPSIso value for this selection was: " << SkipHPSIso << endl;
    cout << "Multiples gen - hps match resolved in " << multHpsGen << " cases" << endl;
    cout << "Multiples hps - L1 match resolved in " << multHpsL1 << " cases" << endl;

    // now compute matching efficiency with rations
    TGraphAsymmErrors* eta_HPS_GEN = new TGraphAsymmErrors();
    TGraphAsymmErrors* phi_HPS_GEN = new TGraphAsymmErrors();
    TGraphAsymmErrors* pt_HPS_GEN = new TGraphAsymmErrors();

    TGraphAsymmErrors* eta_HPS_L1 = new TGraphAsymmErrors();
    TGraphAsymmErrors* phi_HPS_L1 = new TGraphAsymmErrors();
    TGraphAsymmErrors* pt_HPS_L1 = new TGraphAsymmErrors();

    eta_HPS_GEN->SetName ("eff_eta_HPS_GEN");
    phi_HPS_GEN->SetName ("eff_phi_HPS_GEN");
    pt_HPS_GEN->SetName ("eff_pt_HPS_GEN");

    eta_HPS_L1->SetName ("eff_eta_HPS_L1");
    phi_HPS_L1->SetName ("eff_phi_HPS_L1");
    pt_HPS_L1->SetName ("eff_pt_HPS_L1");

    eta_HPS_GEN -> BayesDivide (matchEffEta[1], matchEffEta[0]);
    phi_HPS_GEN -> BayesDivide (matchEffPhi[1], matchEffPhi[0]);
    pt_HPS_GEN -> BayesDivide (matchEffPt[1], matchEffPt[0]);

    eta_HPS_L1 -> BayesDivide (matchEffEta[3], matchEffEta[2]);
    phi_HPS_L1 -> BayesDivide (matchEffPhi[3], matchEffPhi[2]);
    pt_HPS_L1 -> BayesDivide (matchEffPt[3], matchEffPt[2]);

    tOutput->Write();
    eta_HPS_GEN->Write();
    phi_HPS_GEN->Write();
    pt_HPS_GEN->Write();
    eta_HPS_L1->Write();
    phi_HPS_L1->Write();
    pt_HPS_L1->Write();
    fOutput->Close();

    return 1;
}
