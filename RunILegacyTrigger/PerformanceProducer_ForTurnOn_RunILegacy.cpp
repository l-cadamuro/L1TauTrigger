////////////////////////////////////////////////
/// HTauTau TriggerStudies : MuTau & MuMu T&P //
////////////////////////////////////////////////

#include "TH1F.h" 
#include "TH2D.h" 

#include "TMath.h" 
#include "TChain.h" 
#include "TFile.h"
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
#include "PerformanceProducer.h"
using namespace std;


int main (int argc, char** argv)
//void PerformanceProducer()
{	

	bool debug=false;
    
    // INPUT TREES //

   	vector<TString> list;
	
    // gg Run I Pu 40 bx 25 E 13 TeV
    list.push_back("/data_CMS/cms/cadamuro/L1_trigger_data/RunI_Legacy_trigger_trees/ggFusion/RunI_LegacyL1_ggFusion_bx25Pu40E13TeV.root");
    
    TChain * myChain = new TChain("produceNtuple/recoTree");
    int nFiles = list.size();
    
    for(int i=0;i<nFiles;i++)
    {
        myChain->Add(list[i]);
    }


	// stage 2 - tower colection
	myChain->SetBranchAddress ("l1t_stage2_tower_hwEta", &_l1t_stage2_tower_hwEta);

    // uso i nomi di stage 1, ma sto girando sul Run I!!
	//stage-1 branch
    myChain->SetBranchAddress("trig_L1tau_N",    &_l1t_stage1_tau_n);
    myChain->SetBranchAddress("trig_L1tau_et",   &_l1t_stage1_tau_pt);
    myChain->SetBranchAddress("trig_L1tau_eta",  &_l1t_stage1_tau_eta);
    myChain->SetBranchAddress("trig_L1tau_phi",  &_l1t_stage1_tau_phi);
    //myChain->SetBranchAddress("l1t_stage1_tau_hwIso",&_l1t_stage1_tau_hwIso);


    myChain->SetBranchAddress("genTau_N", &_ngenTau);
    myChain->SetBranchAddress("genTau_pt", &_genTauPt);
    myChain->SetBranchAddress("genTau_eta", &_genTauEta);
    myChain->SetBranchAddress("genTau_phi", &_genTauPhi);
    myChain->SetBranchAddress("genTau_DecayMode", &_genTauDecayMode);

    myChain->SetBranchAddress("hpsTau_N", &_nhpsTau);
    myChain->SetBranchAddress("hpsTau_eta", &_hpsTauEta);
    myChain->SetBranchAddress("hpsTau_phi", &_hpsTauPhi);
    myChain->SetBranchAddress("hpsTau_pt", &_hpsTauPt);
    myChain->SetBranchAddress("hpsTau_charge", &_hpsTauCharge);
    myChain->SetBranchAddress("hpsTau_chIso", &_hpsTauChargedIso);
    myChain->SetBranchAddress("hpsTau_phIso", &_hpsTauPhotonsIso);
    myChain->SetBranchAddress("hpsTau_decayMode", &_hpsTauDiscrByDecMode);
    myChain->SetBranchAddress("hpsTau_isoL", &_hpsTauDiscrByLooseIso);
    myChain->SetBranchAddress("hpsTau_isoM", &_hpsTauDiscrByMediumIso);
    myChain->SetBranchAddress("hpsTau_isoMVAL", &_hpsTauDiscrByLooseIsoMVA);
    myChain->SetBranchAddress("hpsTau_isoMVAM", &_hpsTauDiscrByMediumIsoMVA);
    myChain->SetBranchAddress("hpsTau_antiMuL", &_hpsTauDiscrAgainstMuonLoose);
    myChain->SetBranchAddress("hpsTau_antiMuT", &_hpsTauDiscrAgainstMuonTight);
    myChain->SetBranchAddress("hpsTau_antiElL", &_hpsTauDiscrAgainstElecLoose);
    myChain->SetBranchAddress("hpsTau_anitElM", &_hpsTauDiscrAgainstElecMedium);
    myChain->SetBranchAddress("hpsTau_antiElT", &_hpsTauDiscrAgainstElecTight);
    myChain->SetBranchAddress("hpsTau_antiElMVA", &_hpsTauDiscrAgainstElecMVA);
    myChain->SetBranchAddress("hpsTau_TOTEnergy", &_hpsTau_TOTEnergy);
    myChain->SetBranchAddress("hpsTau_ecalEnergy", &_hpsTau_ecalEnergy);
    myChain->SetBranchAddress("hpsTau_hcalEnergy", &_hpsTau_hcalEnergy);
    
    // OUTPUT FILE //
 
    TString name = "1";

	double L1Thr = 0.0;
    bool ApplyRescale = true;
    double RunIScale = 1.;
    if (ApplyRescale) RunIScale = 1.075/1.713;

	name = "./TurnOnTrees/Performance_RunI_ggToHToTauTau_PU40bx25E13TeV";
    if (ApplyRescale) name += "_Rescaled.root";
    else name += "_NotRescaled.root";
	
	TFile *outfile = new TFile(name,"RECREATE");
    
    
   
    //histo definition
    //TH1D*h_IsoEt = new TH1D("h_IsoEt","h_IsoEt",100,0,100);
    //TH2D*h_Scatter_IsoEt_Vs_Pt = new TH2D("h_Scatter_IsoEt_Vs_Pt","h_Scatter_IsoEt_Vs_Pt",300,0,300,100,0,100);
    TH1D*h_EtRes_stage1 = new TH1D("h_EtRes_stage1","h_EtRes_stage1",100,-1,3);
	//TH1D*h_EtRes_stage2 = new TH1D("h_EtRes_stage2","h_EtRes_stage2",100,-1,3);
    TH1D*h_EtaRes_stage1 = new TH1D("h_EtaRes_stage1","h_EtaRes_stage1",100,-0.5,0.5);
	//TH1D*h_EtaRes_stage2 = new TH1D("h_EtaRes_stage2","h_EtaRes_stage2",100,-0.5,0.5);
    TH2D*h_2DforTurnOn_stage1 = new TH2D("h_2DforTurnOn_stage1","h_2DforTurnOn_stage1",500,0,500,500,0,500);
    //TH2D*h_2DforTurnOn_stage2 = new TH2D("h_2DforTurnOn_stage2","h_2DforTurnOn_stage2",500,0,500,500,0,500);

    /////////////////
    // OUTPUT TREE //
    ////////////////
                   
    TTree* treeTnP = new TTree("treeTurnOn", "treeTurnOn");
	//treeTnP->Branch("EventN",	&_EventN);
    treeTnP->Branch("eventN",		&_eventN);  
    
    treeTnP->Branch("recoTau_pt",		&_recoTau_pt);  
    treeTnP->Branch("recoTau_et",		&_recoTau_et); 
    treeTnP->Branch("recoTau_eta",		&_recoTau_eta);  
    treeTnP->Branch("recoTau_phi",		&_recoTau_phi);  
        
    treeTnP->Branch("L1tau_stage1_et",		&_L1tau_stage1_et);  
    treeTnP->Branch("L1tau_stage1_eta",		&_L1tau_stage1_eta);  
    treeTnP->Branch("L1tau_stage1_EtRes",	&_L1tau_stage1_EtRes);  
    treeTnP->Branch("L1tau_stage2_et",		&_L1tau_stage2_et);  
    treeTnP->Branch("L1tau_stage2_eta",		&_L1tau_stage2_eta);  
    treeTnP->Branch("L1tau_stage2_EtRes",	&_L1tau_stage2_EtRes); 
    
    /*
	treeTnP->Branch("Stage2_match_20", &_Stage2_match_20, "Stage2_match_20/I");
    treeTnP->Branch("Stage2_match_25", &_Stage2_match_25, "Stage2_match_25/I");
    treeTnP->Branch("Stage2_match_30", &_Stage2_match_30, "Stage2_match_30/I");
    treeTnP->Branch("Stage2_match_35", &_Stage2_match_35, "Stage2_match_35/I");
	treeTnP->Branch("Stage2_match_40", &_Stage2_match_40, "Stage2_match_40/I");
    treeTnP->Branch("Stage2_match_50", &_Stage2_match_50, "Stage2_match_50/I");
    */

	treeTnP->Branch("Stage1_match_20", &_Stage1_match_20, "Stage1_match_20/I");
    treeTnP->Branch("Stage1_match_25", &_Stage1_match_25, "Stage1_match_25/I");
    treeTnP->Branch("Stage1_match_30", &_Stage1_match_30, "Stage1_match_30/I");
    treeTnP->Branch("Stage1_match_35", &_Stage1_match_35, "Stage1_match_35/I");
	treeTnP->Branch("Stage1_match_40", &_Stage1_match_40, "Stage1_match_40/I");
    treeTnP->Branch("Stage1_match_50", &_Stage1_match_50, "Stage1_match_50/I");
         
    int numEntries = myChain->GetEntries() ;
    int nProcess = numEntries;
    
    // Variables //
	
    //rescale response to 1.0
	
    /*
	double k,h;
	if (ApplyRescale)
	{
		k=(1/1.01);
		h=(1/0.874);
	}
	else
	{
		k=1.0;
		h=1.0;
	}
	*/
	
	// print messages
	//std::cout << "NoApplyIso: " << NoApplyL1Iso <<std::endl;
	std::cout << "ApplyRescale: " << ApplyRescale <<std::endl;
	std::cout << "Using file: " << std::endl;
	for (int i = 0; i < list.size(); i++)  std::cout << list[i] << std::endl;
	std::cout << "Output file is: " << name << std::endl;
	
	//check RECO TAU 
    int Event=0;
    _eventN=0;
    for (int iEvent = 0 ; iEvent < nProcess ; ++iEvent )
    {	
    _eventN=iEvent;
		//cout<<iEvent<< endl;
	
        if (iEvent%1000 == 0) cout<<iEvent<<"/"<<nProcess<<" processed"<<endl;
        myChain->GetEntry(iEvent);
		
        //inizialyzing turn-on variables
		_recoTau_pt	= -999;
		_recoTau_et	= -999;
		_recoTau_eta= -999;
   		_recoTau_phi= -999;

		_L1tau_stage2_eta 	= -999;
		_L1tau_stage2_et	= -999;


        /*
        _Stage2_match_20	 	= -999;
		_Stage2_match_25     	= -999;
		_Stage2_match_30     	= -999;
		_Stage2_match_35     	= -999;
		_Stage2_match_40     	= -999;
		_Stage2_match_50     	= -999;
        */

		_L1tau_stage1_eta	= -999;
		_L1tau_stage1_et	= -999;

        _Stage1_match_20	 	= -999;
		_Stage1_match_25     	= -999;
		_Stage1_match_30     	= -999;
		_Stage1_match_35     	= -999;
		_Stage1_match_40     	= -999;
		_Stage1_match_50     	= -999;

        if(debug) std::cout << "=========== Start loop over gen Tau =========" << std::endl;

		for(int iGen=0; iGen<_ngenTau; iGen++){
			for(int iTau=0; iTau<_nhpsTau; iTau++){
				double dR_GenReco = deltaR(_hpsTauEta[iTau],_hpsTauPhi[iTau],_genTauEta[iGen],_genTauPhi[iGen]);
            	if(dR_GenReco<0.1){ // matching reco-gen
					//if(fabs(_hpsTauEta[iTau])<2.1 && _hpsTauPt[iTau]>20 && _hpsTauDiscrByDecMode[iTau] == 1 &&  _hpsTauDiscrByMediumIso[iTau] ==1){
					if(fabs(_hpsTauEta[iTau])<2.1 && _hpsTauPt[iTau]>20 && _hpsTauDiscrByDecMode[iTau] == 1 ){
						if(debug) std::cout << "hpsTau Pt: " << _hpsTauPt[iTau] << " eta: " << _hpsTauEta[iTau]<< " phi:" << _hpsTauPhi[iTau] << std::endl;
 
                   		_recoTau_pt = _hpsTauPt[iTau];
                   		_recoTau_et = _hpsTau_TOTEnergy[iTau];
                   		_recoTau_eta = _hpsTauEta[iTau];
                   		_recoTau_phi = _hpsTauPhi[iTau];
                   		
						//////////////
						/// STAGE 2 //
						//////////////
						
                        bool matched_stage2=false;
                   		double dRmin=0.5;
                   		int idx2=0;
                   		if (debug) std::cout << "------------- Start loop over L1candidate Tau -------------" << std::endl;
/*
                   		for(int il1tau=0; il1tau<_l1t_stage2_tau_n; il1tau++){
                  			//bool IsoValue_Stage2(int IsoEt, double L1_Et);
                			//bool TauPassIso = IsoValue_Stage2((*_l1t_stage2_tau_hwIso)[il1tau], k*(*_l1t_stage2_tau_pt)[il1tau]);
							
							int L1_hweta = (*_l1t_stage2_tau_hwEta)[il1tau];
							int L1_hwpt = (*_l1t_stage2_tau_hwPt)[il1tau];
							bool TauPassIso = ( (*_l1t_stage2_tau_hwIso)[il1tau] == 1);
							int flags = (*_l1t_stage2_tau_hwQual)[il1tau];
                            int clusterShape = ClusterShape(flags);
                            bool shVeto = PassVeto(clusterShape, shapesVeto, shapesVetoNr);              
                            bool passTot = (shVeto && TauPassIso);
							
							// !!!!!!! HERE COMPUTE ISO

                       		if(fabs((*_l1t_stage2_tau_eta)[il1tau])<2.2 && (*_l1t_stage2_tau_pt)[il1tau] > L1Thr && (passTot==true || NoApplyL1Iso ) ){ // iso
							//if(fabs((*_l1t_stage2_tau_eta)[il1tau])<2.2 && (*_l1t_stage2_tau_pt)[il1tau] > 0 ){ // non-iso
 
                           		if (debug)std::cout << "--->stage-2 Pt: " << (*_l1t_stage2_tau_pt)[il1tau] << " eta: " << (*_l1t_stage2_tau_eta)[il1tau]  << " phi: " << (*_l1t_stage2_tau_phi)[il1tau] <<  std::endl;
                           		double dR = deltaR((*_l1t_stage2_tau_eta)[il1tau], (*_l1t_stage2_tau_phi)[il1tau], _hpsTauEta[iTau], _hpsTauPhi[iTau] );
 
                           		if(dR<dRmin){
                               		matched_stage2=true;

                               		if (debug) std::cout << "dRmin: " << dRmin << " matched with dR: " << dR<< " cluster idx: " << il1tau<<  std::endl;
                               		dRmin=dR;
                               		idx2=il1tau;
                               		if (debug)std::cout << "new dRmin: " << dRmin << std::endl;

                           		}
                       		}
                   		}
                   		if(debug) std::cout << "stage-2: " << matched_stage2 << std::endl;
						_L1tau_stage2_et = k*(*_l1t_stage2_tau_pt)[idx2];
						_L1tau_stage2_eta = (*_l1t_stage2_tau_eta)[idx2];
                   		if(matched_stage2==true){
                       		if (debug)std::cout << "------------- MATCHED CLOSEST TAU WITH DR: " << dRmin << " AND INDEX: " << idx2 << std::endl;
                       		if (debug)std::cout << "----------------------------------------------> good candidate Pt: "<< (*_l1t_stage2_tau_pt)[idx2] << " eta: " << (*_l1t_stage2_tau_eta)[idx2] <<" phi: " << (*_l1t_stage2_tau_phi)[idx2] <<  std::endl;
 
                       		double EtRes = ((k*(*_l1t_stage2_tau_pt)[idx2])/_hpsTauPt[iTau]);
                       		double EtaRes = (((*_l1t_stage2_tau_eta)[idx2])-(_hpsTauEta[iTau]));

                       		h_EtRes_stage2->Fill(EtRes);
                       		h_EtaRes_stage2->Fill(EtaRes);
                       		_L1tau_stage2_EtRes = EtRes;
                       		h_2DforTurnOn_stage2->Fill(_hpsTauPt[iTau],k*(*_l1t_stage2_tau_pt)[idx2]);
 							
                       		if(k*(*_l1t_stage2_tau_pt)[idx2]>=20)	  _Stage2_match_20 = 1;
                       		else if(k*(*_l1t_stage2_tau_pt)[idx2]<20) _Stage2_match_20 = 0;
                       		if(k*(*_l1t_stage2_tau_pt)[idx2]>=25)     _Stage2_match_25 = 1;
                       		else if(k*(*_l1t_stage2_tau_pt)[idx2]<25) _Stage2_match_25 = 0;
                       		if(k*(*_l1t_stage2_tau_pt)[idx2]>=30)     _Stage2_match_30 = 1;
                       		else if(k*(*_l1t_stage2_tau_pt)[idx2]<30) _Stage2_match_30 = 0;
                       		if(k*(*_l1t_stage2_tau_pt)[idx2]>=35)     _Stage2_match_35 = 1;
                       		else if(k*(*_l1t_stage2_tau_pt)[idx2]<35) _Stage2_match_35 = 0;
                       		if(k*(*_l1t_stage2_tau_pt)[idx2]>=40)     _Stage2_match_40 = 1;
                       		else if(k*(*_l1t_stage2_tau_pt)[idx2]<40) _Stage2_match_40 = 0;
                       		if(k*(*_l1t_stage2_tau_pt)[idx2]>=50)     _Stage2_match_50 = 1;
                       		else if(k*(*_l1t_stage2_tau_pt)[idx2]<50) _Stage2_match_50 = 0;	


                   		}
                   		else if(matched_stage2==false){
                       		h_2DforTurnOn_stage2->Fill(_hpsTauPt[iTau],0.0);
 
                       		//_L1tau_stage2_et = (*_l1t_stage2_tau_pt)[idx2];
                       		//_L1tau_stage2_eta = (*_l1t_stage2_tau_eta)[idx2];

                       		_Stage2_match_20        = 0;
                       		_Stage2_match_25        = 0;
                       		_Stage2_match_30        = 0;
                       		_Stage2_match_35        = 0;
                       		_Stage2_match_40        = 0;
                       		_Stage2_match_50        = 0;

                   		}
                   		if(debug)std::cout << "" << std::endl;
*/ 
 						
						//////////////
						/// STAGE 1 //
						//////////////

                   		bool matched_stage1=false;
                   		//double dRmin1=1.0;
						double dRmin1=0.5;
                   		int idx1=0;
                   		if (debug) std::cout << "------------- Start loop over L1candidate Tau -------------" << std::endl;

                   		for(int il1tau=0; il1tau<_l1t_stage1_tau_n; il1tau++){
 
                    		if(fabs((*_l1t_stage1_tau_eta)[il1tau])<2.2 && (*_l1t_stage1_tau_pt)[il1tau] > L1Thr){  //iso
							//if(fabs((*_l1t_stage1_tau_eta)[il1tau])<2.2 && (*_l1t_stage1_tau_pt)[il1tau] > 0 ){	//non-iso

                        		if (debug)std::cout << "--->stage-1 Pt: " << (*_l1t_stage1_tau_pt)[il1tau] << " eta: " << (*_l1t_stage1_tau_eta)[il1tau]  << " phi: " << (*_l1t_stage1_tau_phi)[il1tau] <<  std::endl;
                           		double dR = deltaR((*_l1t_stage1_tau_eta)[il1tau], (*_l1t_stage1_tau_phi)[il1tau], _hpsTauEta[iTau], _hpsTauPhi[iTau] );
 
                           		if(dR<dRmin1){
                            		matched_stage1=true;

                               		if (debug) std::cout << "dRmin1: " << dRmin1 << " matched with dR: " << dR<< " cluster idx: " << il1tau<<  std::endl;
                               		dRmin1=dR;
                               		idx1=il1tau;
                               		if (debug)std::cout << "new dRmin1: " << dRmin1 << std::endl;

                           		}
                       		}
                   		}
                   		if(debug) std::cout << "stage-1 :" << matched_stage1 << std::endl;
						_L1tau_stage1_et = RunIScale*(*_l1t_stage1_tau_pt)[idx1];
						_L1tau_stage1_eta = (*_l1t_stage1_tau_eta)[idx1];
                    	if(matched_stage1==true){
                        	if (debug)std::cout << "------------- MATCHED CLOSEST TAU WITH DR: " << dRmin1 << " AND INDEX: " << idx1 << std::endl;
                        	if (debug)std::cout << "----------------------------------------------> good candidate Pt: "<< (*_l1t_stage1_tau_pt)[idx1] << " eta: " << (*_l1t_stage1_tau_eta)[idx1] <<" phi: " << (*_l1t_stage1_tau_phi)[idx1] <<  std::endl;
  
                        	double EtRes = (RunIScale*((*_l1t_stage1_tau_pt)[idx1])/_hpsTauPt[iTau]);
                        	double EtaRes = (((*_l1t_stage1_tau_eta)[idx1])-(_hpsTauEta[iTau]));	
 
                        	h_EtRes_stage1->Fill(EtRes);
                        	h_EtaRes_stage1->Fill(EtaRes);
                        	_L1tau_stage1_EtRes = EtRes;
                        	h_2DforTurnOn_stage1->Fill(_hpsTauPt[iTau],RunIScale*(*_l1t_stage1_tau_pt)[idx1]);
  							
                            //std::cout << (*_l1t_stage1_tau_pt)[idx1] << std::endl;
                            
                       		if(RunIScale*(*_l1t_stage1_tau_pt)[idx1]>=20)     _Stage1_match_20 = 1;
                       		else if(RunIScale*(*_l1t_stage1_tau_pt)[idx1]<20) _Stage1_match_20 = 0;
							if(RunIScale*(*_l1t_stage1_tau_pt)[idx1]>=25)     _Stage1_match_25 = 1;
                       		else if(RunIScale*(*_l1t_stage1_tau_pt)[idx1]<25) _Stage1_match_25 = 0;
                       		if(RunIScale*(*_l1t_stage1_tau_pt)[idx1]>=30)     _Stage1_match_30 = 1;
                       		else if(RunIScale*(*_l1t_stage1_tau_pt)[idx1]<30) _Stage1_match_30 = 0;
                       		if(RunIScale*(*_l1t_stage1_tau_pt)[idx1]>=35)     _Stage1_match_35 = 1;
                       		else if(RunIScale*(*_l1t_stage1_tau_pt)[idx1]<35) _Stage1_match_35 = 0;
                       		if(RunIScale*(*_l1t_stage1_tau_pt)[idx1]>=40)     _Stage1_match_40 = 1;
                       		else if(RunIScale*(*_l1t_stage1_tau_pt)[idx1]<40) _Stage1_match_40 = 0;
                       		if(RunIScale*(*_l1t_stage1_tau_pt)[idx1]>=50)     _Stage1_match_50 = 1;
                       		else if(RunIScale*(*_l1t_stage1_tau_pt)[idx1]<50) _Stage1_match_50 = 0;
   
                    	}
                    	else if(matched_stage1==false){
                    		h_2DforTurnOn_stage1->Fill(_hpsTauPt[iTau],0.0);
                          
                         //_L1tau_stage1_et = (*_l1t_stage1_tau_pt)[idx1];
                          //_L1tau_stage1_eta = (*_l1t_stage1_tau_eta)[idx1];
 
                          _Stage1_match_20        = 0;
                          _Stage1_match_25        = 0;
                          _Stage1_match_30        = 0;
                          _Stage1_match_35        = 0;
                          _Stage1_match_40        = 0;
                          _Stage1_match_50        = 0;
 
                    	}
                   		if(debug)std::cout << "" << std::endl;
               		}   // selection on the hpsTau
           		}   //end if matching condition between gen-reco tau
			} 	//end loop over hpsTau
 		}	//end loop over genTau	  
    	
        if(_recoTau_pt>-999)treeTnP->Fill() ; 

    	Event++;
    } // end loop entries
 
    
    // Save tree
	treeTnP->Write();
    	    
	//h_IsoEt->Write();
    //h_Scatter_IsoEt_Vs_Pt->Write(); 
	h_EtRes_stage1->Scale(1/h_EtRes_stage1->Integral());
    h_EtRes_stage1->Write();
    //h_EtRes_stage2->Scale(1/h_EtRes_stage2->Integral());
    //h_EtRes_stage2->Write();
	h_EtaRes_stage1->Scale(1/h_EtaRes_stage1->Integral());
    h_EtaRes_stage1->Write();
    //h_EtaRes_stage2->Scale(1/h_EtaRes_stage2->Integral());
    //h_EtaRes_stage2->Write();    
    h_2DforTurnOn_stage1->Write();
	//h_2DforTurnOn_stage2->Write();
    myChain->Delete();
        
    return 1;
        
}
 
 
// int main()
//{PerformanceProducer();}
