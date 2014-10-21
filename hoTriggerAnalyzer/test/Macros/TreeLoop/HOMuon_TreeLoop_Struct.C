#define HOMuon_TreeLoop_Struct_cxx
// The class definition in HOMuon_TreeLoop_Struct.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("HOMuon_TreeLoop_Struct.C")
// Root > T->Process("HOMuon_TreeLoop_Struct.C","some options")
// Root > T->Process("HOMuon_TreeLoop_Struct.C+")
//

#include "HOMuon_TreeLoop_Struct.h"
#include <TH2.h>
#include <TStyle.h>
#include <TH1F.h>
#include <sstream>
#include <string>
#include <map>
#include "iostream"
#include "math.h"
#include <vector>

//#include "HoMuonTrigger/hoTriggerAnalyzer/test/Macros/HistogramBuilderTwo.h"

//#include "HoMuonTrigger/hoTriggerAnalyzer/interface/HistogramBuilder.h"

void HOMuon_TreeLoop_Struct::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   outRootFile = TFile::Open("L1MuonHistogram.root", "RECREATE");

}

void HOMuon_TreeLoop_Struct::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t HOMuon_TreeLoop_Struct::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either HOMuon_TreeLoop_Struct::GetEntry() or TBranch::GetEntry()
  // to read either all or the required parts of the data. When processing
  // keyed objects with PROOF, the object is already loaded and is available
  // via the fObject pointer.
  //
  // This function should contain the "body" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.


  fChain->GetEntry(entry);
  //if(entry > 100000) return kFALSE;  //To speed up testing only run over part of dataset.
  if(entry%1000==0) std::cout << entry << std::endl;
  
  double weight = Generator_Weights;

  //Must be in the Barrel
  for(unsigned int l1MuonB_index = 0; l1MuonB_index < L1Muon_Etas->size(); l1MuonB_index++){
    if(fabs(L1Muon_Etas->at(l1MuonB_index)) > barrel_eta) continue;
    if(weight < 1.0) continue;  //Checking if Large Weights are the Issue
    L1MuonCand l1MuonCand;
    //L1MuonCand * l1MuonCand = {0,0,0, false, false, 0, 0, 0, false, false, 0, 0, 0};
    setL1Info(l1MuonB_index, &l1MuonCand);
    setHOInfo(&l1MuonCand);
    setHLTInfo(&l1MuonCand);
    fillHistograms(&l1MuonCand, weight);
    //free(&l1MuonCand); //Clean Up Struct
  }
 return kTRUE;
}
 

void HOMuon_TreeLoop_Struct::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void HOMuon_TreeLoop_Struct::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
  
  cout << "working1" << endl;
  outRootFile->cd();
  outRootFile->Write();
  cout << "working2" << endl;  

}


/*
 * Analyzes Data stored in HO Tree to Make Histograms
 * The Key names are based on the conditions the events must pass
 */
void HOMuon_TreeLoop_Struct::fillHistograms(L1MuonCand * l1MuonCand, double weight){
  
  //L1 Muon in the Barrel
  string key = "L1MuonB";
  histogramBuilder.fillWeightHistograms(weight, key);
  fillL1Muon_HO_Histograms(l1MuonCand, weight, key);
  
  
  //L1 Muon in the Barrel with a MIP
  if(l1MuonCand->hasMIP){
    key = "L1MuonB_MIP";
    fillL1Muon_HO_Histograms(l1MuonCand, weight, key);
    fillDeltaL1Muon_HO_Histograms(l1MuonCand, weight, key);
  } else { //No MIP
    key = "L1MuonB_NoMIP";
    fillL1Muon_HO_Histograms(l1MuonCand, weight, key);
  }
  //L1 Muon in the Barrel, Event has HLT
  if(l1MuonCand->HLT_InEvent){
    key = "L1MuonB_HLTInEvent";
    fillL1Muon_HO_Histograms(l1MuonCand, weight, key);
    fillHLT_Histograms(l1MuonCand, weight, key);
    fillDeltaL1Muon_HLT_Histograms(l1MuonCand, weight, key);
  } else { //Event has No HLT
    key = "L1MuonB_NoHLTInEvent";
    fillL1Muon_HO_Histograms(l1MuonCand, weight,key);
  }
  
  //if(l1MuonCand->HLT_InEvent){ 
  //L1 Muon in the Barrel, Muon R Matched to HLT
  if(l1MuonCand->hasHLTRMatch) {
    key = "L1MuonB_HLTRMatch";
    fillL1Muon_HO_Histograms(l1MuonCand, weight, key);
    fillHLT_Histograms(l1MuonCand, weight, key);
    fillDeltaL1Muon_HLT_Histograms(l1MuonCand, weight, key);
  } else { // HLT is not Match to L1 Muon (but it is an HLT Event
    key = "L1MuonB_NoHLTRMatch";
    fillL1Muon_HO_Histograms(l1MuonCand, weight, key);
    //fillHLT_Histograms(l1MuonCand, weight, key);
    //fillDeltaL1Muon_HLT_Histograms(l1MuonCand, weight, key);
  }  

  if(l1MuonCand->hasMIP){
    if(l1MuonCand->hasHLTRMatch){
      key = "L1MuonB_MIP_HLTRMatch";
      fillL1Muon_HO_Histograms(l1MuonCand, weight, key);
      fillHLT_Histograms(l1MuonCand, weight, key);
      fillDeltaL1Muon_HLT_Histograms(l1MuonCand, weight, key);
    } else { // HLT is not Match to L1 Muon (but it is an HLT Event
      key = "L1MuonB_MIP_NoHLTRMatch";
      fillL1Muon_HO_Histograms(l1MuonCand, weight, key);
      //fillHLT_Histograms(l1MuonCand, weight, key);
      //fillDeltaL1Muon_HLT_Histograms(l1MuonCand, weight, key);
    }
  }
}
    
// Fill the L1Muon and HO Histograms
void HOMuon_TreeLoop_Struct::fillL1Muon_HO_Histograms( L1MuonCand * l1MuonCand, double weight, string key){
  string l1muon_key = key + "_l1muon";
  
  histogramBuilder.fillCountHistogram(key,weight);
  
  histogramBuilder.fillL1MuonPtHistograms(l1MuonCand->l1Muon_pt, l1muon_key, weight);
  
  histogramBuilder.fillEtaPhiHistograms(l1MuonCand->l1Muon_eta, 
					l1MuonCand->l1Muon_phi,
					l1muon_key, weight);
  //Information for Max HO Rec Hits Matched to a L1 Muon Barrel
  string horeco_key = key + "_horeco";
  histogramBuilder.fillEnergyHistograms(l1MuonCand->horeco_maxE_energy,
					horeco_key, weight);
  
}

  
void HOMuon_TreeLoop_Struct::fillDeltaL1Muon_HO_Histograms( L1MuonCand * l1MuonCand, double weight, string key){
  string l1muon_horeco_key = key + "_l1muon_horeco";
  histogramBuilder.fillDeltaEtaDeltaPhiHistograms(l1MuonCand->l1Muon_eta, l1MuonCand->horeco_maxE_eta,
						  l1MuonCand->l1Muon_phi, l1MuonCand->horeco_maxE_phi,
						  l1muon_horeco_key, weight);
}
 

// L1 Muons in the Barrel  Matched to HLT Props
void HOMuon_TreeLoop_Struct::fillHLT_Histograms( L1MuonCand * l1MuonCand, double weight, string key){
  string hlt_key = key +"_hlt";
  histogramBuilder.fillPtHistograms(l1MuonCand->hlt_minR_pt, hlt_key, weight);
  histogramBuilder.fillEtaPhiHistograms(l1MuonCand->hlt_minR_eta, 
					l1MuonCand->hlt_minR_phi,
					hlt_key, weight);
}

void HOMuon_TreeLoop_Struct::fillDeltaL1Muon_HLT_Histograms( L1MuonCand * l1MuonCand, double weight, string key){
  string l1muon_hlt_key = key + "l1muon_hlt";
  histogramBuilder.fillDeltaEtaDeltaPhiHistograms(l1MuonCand->l1Muon_eta, l1MuonCand->hlt_minR_eta,
						  l1MuonCand->l1Muon_phi, l1MuonCand->hlt_minR_phi,
						  l1muon_hlt_key, weight);
}


void HOMuon_TreeLoop_Struct::setL1Info(unsigned int l1Muon_index,  L1MuonCand * l1MuonCand){
  //l1MuonCand->l1Muon_index = l1Muon_index;
  //cout << L1Muon_Etas->at(l1Muon_index) << endl;
  l1MuonCand->l1Muon_eta = L1Muon_Etas->at(l1Muon_index);
  l1MuonCand->l1Muon_phi = L1Muon_Phis->at(l1Muon_index);
  l1MuonCand->l1Muon_pt = L1Muon_Pts->at(l1Muon_index);
}

void HOMuon_TreeLoop_Struct::setHOInfo( L1MuonCand * l1MuonCand){
  int maxE_horeco_index=-1;
  float maxE_horeco_energy =0; //Maximum HOReco Energy within DeltaR of L1Muon
  bool hasMIP = false;
  bool hasHO = false;
  for(unsigned int horeco_index = 0; horeco_index < HOReco_Etas->size(); horeco_index++){
    float horeco_eta, horeco_phi, horeco_energy;
    horeco_eta = HOReco_Etas->at(horeco_index);
    horeco_phi = HOReco_Phis->at(horeco_index);
    horeco_energy = HOReco_Energies->at(horeco_index);
    //Interested in Collection of HO Rec Hits, inside Delta R of L1Muon
    if(calcDeltaR2(l1MuonCand->l1Muon_eta, horeco_eta, 
		   l1MuonCand->l1Muon_phi, horeco_phi)< (RMip_Max*RMip_Max)){
      if(horeco_energy > maxE_horeco_energy){
	hasHO=true;
	maxE_horeco_index = horeco_index;
	maxE_horeco_energy=horeco_energy; 
      }
    }
  }
  l1MuonCand->hasHO = hasHO;
  if(hasHO){
    l1MuonCand->horeco_maxE_eta = HOReco_Etas->at(maxE_horeco_index);
    l1MuonCand->horeco_maxE_phi = HOReco_Phis->at(maxE_horeco_index);
    l1MuonCand->horeco_maxE_energy = HOReco_Energies->at(maxE_horeco_index);
  } else { //Call the Energy 0 and use the positon of the L1Muon
    l1MuonCand->horeco_maxE_eta = l1MuonCand->l1Muon_eta;
    l1MuonCand->horeco_maxE_phi = l1MuonCand->l1Muon_phi;
    l1MuonCand->horeco_maxE_energy = 0;
  }
  if(maxE_horeco_energy > threshold) hasMIP=true;
  l1MuonCand->hasMIP=hasMIP;
}

void HOMuon_TreeLoop_Struct::setHLTInfo( L1MuonCand * l1MuonCand){
  int hlt_index = -1;
  float min_hlt_R2 = 100; //Matches To Nearest HLT
  bool HLT_InEvent = false;
  bool hasHLTRMatch = false;
  float hlt_eta, hlt_phi;
  for(unsigned int hltProp_index = 0; hltProp_index < hltMu5PropToRPC1_Etas->size(); hltProp_index++){
    HLT_InEvent = true;
    hlt_eta = hltMu5PropToRPC1_Etas->at(hltProp_index);
    hlt_phi = hltMu5PropToRPC1_Phis->at(hltProp_index);
    float hlt_l1_DeltaR2 = calcDeltaR2(l1MuonCand->l1Muon_eta, hlt_eta, 
				       l1MuonCand->l1Muon_phi, hlt_phi);
   
    if(hlt_l1_DeltaR2 < min_hlt_R2){ 
      min_hlt_R2 = hlt_l1_DeltaR2; 
      hlt_index = hltProp_index;
    }
  }
  if( min_hlt_R2 < RHlt_Max*RHlt_Max) hasHLTRMatch = true;

  l1MuonCand->HLT_InEvent = HLT_InEvent;
  l1MuonCand->hasHLTRMatch = hasHLTRMatch;
  if(HLT_InEvent){
    //Should only be assigned, if a match exists
    l1MuonCand->hlt_minR_eta = hltMu5PropToRPC1_Etas->at(hlt_index);
    l1MuonCand->hlt_minR_phi = hltMu5PropToRPC1_Phis->at(hlt_index);
    l1MuonCand->hlt_minR_pt = hltMu5PropToRPC1_Pts->at(hlt_index);
  }
}

float HOMuon_TreeLoop_Struct::calcDeltaR2(float eta1, float eta2, float phi1, float phi2){
  
  float delta_eta, delta_phi, delta_R2;

  delta_eta = eta1 - eta2;
  delta_phi = WrapCheck(phi1,phi2); //Finds difference in phi                              
  delta_R2 = delta_eta*delta_eta+delta_phi*delta_phi;
  //The L1 Muon is compared with all HO Rec Hits above Threshold.                                          
  return delta_R2;
}

float HOMuon_TreeLoop_Struct::WrapCheck(float phi1, float phi2){
  static const float PI = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348;

  float delta_phi = phi1 - phi2;
  if (delta_phi > PI) delta_phi-= 2*PI;
  if(delta_phi < -1*PI) delta_phi += 2*PI;
  
  return delta_phi;
}                 


 
   // Compare L1 Muon Trigger to Propagated Generator Muons
   
  /*
  for(unsigned int l1Muon_index = 0; l1Muon_index < L1Muon_Etas->size(); l1Muon_index++){
    
    float l1Muon_eta, l1Muon_phi;
    l1Muon_eta = L1Muon_Etas->at(l1Muon_index);
    l1Muon_phi = L1Muon_Phis->at(l1Muon_index);

    for(unsigned int genMuonProp_index = 0; genMuonProp_index < GenMuonPropToRPC1_Etas->size(); 
	genMuonProp_index++){
   
      float genMuonProp_eta, genMuonProp_phi;
      genMuonProp_eta = GenMuonPropToRPC1_Etas->at(genMuonProp_index);
      genMuonProp_phi = GenMuonPropToRPC1_Phis->at(genMuonProp_index);
      

      std::string l1MuonMatchGenMuonProp_key = "L1MuonMatchGenMuonProp";
      histogramBuilder.fillDeltaEtaDeltaPhiHistograms(l1Muon_eta, genMuonProp_eta,
						      l1Muon_phi, genMuonProp_phi,
						      l1MuonMatchGenMuonProp_key);

      if(GenMuonPropToRPC1_pdgIds->at(genMuonProp_index) ==13){
	std::string l1MuonMatchGenMuonPropMuMinus_key = "L1MuonMatchGenMuonPropMuMinus";
	histogramBuilder.fillDeltaEtaDeltaPhiHistograms(l1Muon_eta, genMuonProp_eta,
							l1Muon_phi, genMuonProp_phi,
							l1MuonMatchGenMuonPropMuMinus_key);
      }

      if(GenMuonPropToRPC1_pdgIds->at(genMuonProp_index) ==-13){
	std::string l1MuonMatchGenMuonPropMuPlus_key = "L1MuonMatchGenMuonPropMuPlus";
	histogramBuilder.fillDeltaEtaDeltaPhiHistograms(l1Muon_eta, genMuonProp_eta,
							l1Muon_phi, genMuonProp_phi,
							l1MuonMatchGenMuonPropMuPlus_key);
      }
  */

      /*      
      //Check for Particle AntiParticle Pairs
      //cout << "Wrap Check " << histogramBuilder.WrapCheck(l1Muon_phi, genMuonProp_phi) << endl;
      //Demand there is a trigger matched to a generator on the Opposite end of the detector.
      if(fabs(histogramBuilder.WrapCheck(l1Muon_phi, genMuonProp_phi)) > 2){
	// Loop Over All Generator Particles Again
	for(unsigned int i = 0; i < GenMuonPropToRPC1_Etas->size(); i++){
	  //std::cout << GenMuonPropToRPC1_pdgIds->at(i) << ", "; 
	  float genMuonPropAP_eta, genMuonPropAP_phi;
	  genMuonPropAP_eta = GenMuonPropToRPC1_Etas->at(i);
	  genMuonPropAP_phi = GenMuonPropToRPC1_Phis->at(i);

	  std::string l1MuonMatchGenMuonPropAP_key = "L1MuonMatchGenMuonPropAP";
	  histogramBuilder.fillDeltaEtaDeltaPhiHistograms(l1Muon_eta, genMuonPropAP_eta,
							  l1Muon_phi, genMuonPropAP_phi,
							  l1MuonMatchGenMuonPropAP_key);
	}
      }
      */


