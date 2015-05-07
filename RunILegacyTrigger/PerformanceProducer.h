#include "TH1F.h" 
#include "TMath.h" 
#include "TChain.h" 

#include <map>
#include <sstream>
#include "TString.h"
#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <TString.h>

// bool IsoValue_Stage2(int IsoEt, double L1_Et){
// 
//     bool passIso=false;
//     double isoRel = IsoEt/L1_Et;
// 
//     if(L1_Et<20){
//     	if( (isoRel)<(0.33) ){
//         //if( (IsoEt)<1) ){
//         
// 			passIso = true;
// 			//std::cout << "ascissa: " << L1_Et<< "il valore della funzione --> " << isoFunction(L1_Et) << " valore dell'isolamento: "<< isoFunction(L1_Et)<<  std::endl;
// 	    	//std::cout << "-------> good version - isoEt: " << (IsoEt/L1_Et) << " pass iso? " << passIso <<  std::endl; 
// 		}
//     }
// 	else if((L1_Et>=20 && L1_Et<30)){
// 	
//     	if( (isoRel)<(0.23) ){
// 			passIso = true;
// 			//std::cout << "ascissa: " << L1_Et<< "il valore della funzione --> " << isoFunction(L1_Et) << " valore dell'isolamento: "<< isoFunction(L1_Et)<<  std::endl;
// 	    	//std::cout << "-------> good version - isoEt: " << (IsoEt/L1_Et) << " pass iso? " << passIso <<  std::endl; 
// 		}
// 	}
// 	else if((L1_Et>=30 && L1_Et<34)){
// 	
//     	if( (isoRel)<(0.19) ){
// 			passIso = true;
// 			//std::cout << "ascissa: " << L1_Et<< "il valore della funzione --> " << isoFunction(L1_Et) << " valore dell'isolamento: "<< isoFunction(L1_Et)<<  std::endl;
// 	    	//std::cout << "-------> good version - isoEt: " << (IsoEt/L1_Et) << " pass iso? " << passIso <<  std::endl; 
// 		}
// 	}
// 	else if((L1_Et>=34 && L1_Et<40)){
// 	
//     	if( (isoRel)<(0.16) ){
// 			passIso = true;
// 			//std::cout << "ascissa: " << L1_Et<< "il valore della funzione --> " << isoFunction(L1_Et) << " valore dell'isolamento: "<< isoFunction(L1_Et)<<  std::endl;
// 	    	//std::cout << "-------> good version - isoEt: " << (IsoEt/L1_Et) << " pass iso? " << passIso <<  std::endl; 
// 		}
// 	}
// 	else if((L1_Et>=40 && L1_Et<50)){
// 	
//     	if( (isoRel)<(0.155) ){
// 			passIso = true;
// 			//std::cout << "ascissa: " << L1_Et<< "il valore della funzione --> " << isoFunction(L1_Et) << " valore dell'isolamento: "<< isoFunction(L1_Et)<<  std::endl;
// 	    	//std::cout << "-------> good version - isoEt: " << (IsoEt/L1_Et) << " pass iso? " << passIso <<  std::endl; 
// 		}
// 	}
// 	//else if((L1_Et>=50 && L1_Et<100)){
// 	else if(L1_Et>=50){
// 
//     	if( (isoRel)<(0.16) ){
// 			passIso = true;
// 			//std::cout << "ascissa: " << L1_Et<< "il valore della funzione --> " << isoFunction(L1_Et) << " valore dell'isolamento: "<< isoFunction(L1_Et)<<  std::endl;
// 	    	//std::cout << "-------> good version - isoEt: " << (IsoEt/L1_Et) << " pass iso? " << passIso <<  std::endl; 
// 		}
// 	}
// //     else if(L1_Et>=100){
// // 		passIso = true;
// // 		//std::cout << "Ultimo intervallo dell´isolamento " << " pass iso? " << passIso << std::endl;
// //     }    
// 
// 	return passIso;
// 
// } 

double deltaR(double eta1, double phi1, double eta2, double phi2)
{
  double deta = eta1 - eta2 ;
  double dphi = phi1 - phi2 ;
  while (dphi > TMath::Pi() ) 
    dphi -= 2*TMath::Pi();
  while (dphi <= -TMath::Pi() ) 
    dphi += 2*TMath::Pi();
  return sqrt(deta*deta + dphi*dphi); 
}

//============ VARIABLE INPUT TREE =================

	//stage-2 variables
	std::vector<int> *_l1t_stage2_tower_hwEta  = 0;

	int _l1t_stage2_tau_n;

	std::vector<double> *_l1t_stage2_tau_pt  = 0;
    std::vector<double> *_l1t_stage2_tau_eta = 0;	
    std::vector<double> *_l1t_stage2_tau_phi = 0;
	std::vector<int> 	*_l1t_stage2_tau_hwIso = 0;
	std::vector<int>        *_l1t_stage2_tau_hwQual = 0;
	std::vector<int>        *_l1t_stage2_tau_hwEta = 0;
	std::vector<int>        *_l1t_stage2_tau_hwPt = 0;
	std::vector<int>        *_l1t_stage2_tau_hwPhi = 0;

    //stage-1 variables
	int _l1t_stage1_tau_n;
	std::vector<double> *_l1t_stage1_tau_pt  = 0;
    std::vector<double> *_l1t_stage1_tau_eta = 0;	
    std::vector<double> *_l1t_stage1_tau_phi = 0;
	std::vector<int> 	*_l1t_stage1_tau_hwIso = 0;

	Int_t _ngenTau;
    float _genTauPt[5];
    float _genTauEta[5];
    float _genTauPhi[5];
    int _genTauDecayMode[5];

    Int_t _nhpsTau;
    int _hpsTau_DecayType[100];
    float _hpsTau_TOTEnergy[100];
    float _hpsTau_ecalEnergy[100];
    float _hpsTau_hcalEnergy[100];
    
    float _hpsTauEta[100];
    float _hpsTauPhi[100];
    float _hpsTauPt[100];
    
    float _hpsTauJetPt[100];
    float _hpsTauLeadPionPt[100];
    float _hpsTauLeadTrackPt[100];
    Int_t _hpsTauCharge[100];
    float _hpsTauChargedIso[100];
    float _hpsTauPhotonsIso[100];
    float _hpsTauDiscrByDecMode[100];
    float _hpsTauDiscrByLooseIso[100];
    float _hpsTauDiscrByMediumIso[100];
    float _hpsTauDiscrByLooseIsoMVA[100];
    float _hpsTauDiscrByMediumIsoMVA[100];
    float _hpsTauDiscrAgainstMuonLoose[100];
    float _hpsTauDiscrAgainstMuonTight[100];
    float _hpsTauDiscrAgainstElecLoose[100];
    float _hpsTauDiscrAgainstElecMedium[100];
    float _hpsTauDiscrAgainstElecTight[100];
    float _hpsTauDiscrAgainstElecMVA[100];

//============ VARIABLE OUTPUT TREE =================

	int _eventN;
    
	double _recoTau_pt;
	double _recoTau_et;
	double _recoTau_eta;
    double _recoTau_phi;

    double _L1tau_stage2_et;
    double _L1tau_stage2_eta;
	double _L1tau_stage2_EtRes;
    double _L1tau_stage1_et;
    double _L1tau_stage1_eta;
	double _L1tau_stage1_EtRes;
    
    int _Stage2_match_20;
    int _Stage2_match_25;
    int _Stage2_match_30;
    int _Stage2_match_35;
    int _Stage2_match_40;
    int _Stage2_match_50;

    
    int _Stage1_match_20;
    int _Stage1_match_25;
    int _Stage1_match_30;
    int _Stage1_match_35;
    int _Stage1_match_40;
    int _Stage1_match_50;

