#define HOMuon_TreeLoop_Plotter_cxx
// The class definition in HOMuon_TreeLoop_Plotter.h has been generated automatically
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
// Root > T->Process("HOMuon_TreeLoop_Plotter.C")
// Root > T->Process("HOMuon_TreeLoop_Plotter.C","some options")
// Root > T->Process("HOMuon_TreeLoop_Plotter.C+")
//

#include "HOMuon_TreeLoop_Plotter.h"
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



void HOMuon_TreeLoop_Plotter::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   outRootFile = TFile::Open("L1MuonHistogram.root", "RECREATE");

}

void HOMuon_TreeLoop_Plotter::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t HOMuon_TreeLoop_Plotter::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either HOMuon_TreeLoop_Plotter::GetEntry() or TBranch::GetEntry()
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
  std::cout << entry << std::endl;
  
  /*
   * Weighting Information
   */  
  double weight = Generator_Weights;

  std::string countAllEvents_key = "Events_All";
  std::string countAllEventsW_key = "Weighted Events";
  histogramBuilder.fillCountHistogram(countAllEvents_key);
  histogramBuilder.fillCountHistogram(countAllEventsW_key, weight);

  /*
   * L1 Muons
   */
  std::string l1Muon_key = "l1Muon";
  for(unsigned int l1Muon_index = 0; l1Muon_index < L1Muon_Etas->size(); l1Muon_index++){
    histogramBuilder.fillEtaPhiHistograms(L1Muon_Etas->at(l1Muon_index),
                                          L1Muon_Phis->at(l1Muon_index),
                                          l1Muon_key, weight);
    histogramBuilder.fillL1MuonPtHistograms(L1Muon_Pts->at(l1Muon_index), l1Muon_key, weight);
    histogramBuilder.fillCountHistogram(l1Muon_key, weight);
  }

  /*
   * L1 Muon Barrel
   */
  
  string l1MuonB_key = "L1Muon_Barrel";
  vector<unsigned int> l1MuonBIndices;
  l1MuonBIndices.clear();
  for(unsigned int l1MuonB_index = 0; l1MuonB_index < L1Muon_Etas->size(); l1MuonB_index++){
    if(fabs(L1Muon_Etas->at(l1MuonB_index)) > barrel_eta) continue;
    l1MuonBIndices.push_back(l1MuonB_index);

    histogramBuilder.fillCountHistogram(l1MuonB_key,weight);
    histogramBuilder.fillL1MuonPtHistograms(L1Muon_Pts->at(l1MuonB_index), l1MuonB_key, weight);
    histogramBuilder.fillEtaPhiHistograms(L1Muon_Etas->at(l1MuonB_index), 
					  L1Muon_Phis->at(l1MuonB_index),
					  l1MuonB_key, weight);
    
  }

  /*
   * HO Rec Hit
   */
  for(unsigned int hoReco_index = 0; hoReco_index < HOReco_Etas->size(); hoReco_index++){
    std::string hoReco_key = "HO_Reco";
    histogramBuilder.fillEtaPhiHistograms(HOReco_Etas->at(hoReco_index),
                                          HOReco_Phis->at(hoReco_index),
                                          hoReco_key, weight);
    histogramBuilder.fillEnergyHistograms(HOReco_Pts->at(hoReco_index),
					  hoReco_key, weight); //Typo in Tree Builder should have been labeled HOeco_Energies
  }

  /*
   * Generator Muon
   */

  for(unsigned int genMuon_index = 0; genMuon_index < Generator_Etas->size(); genMuon_index++){
    std::string genMuon_key = "Generator Muons";
    histogramBuilder.fillCountHistogram(genMuon_key,weight);
    histogramBuilder.fillEtaPhiHistograms(Generator_Etas->at(genMuon_index),
                                          Generator_Phis->at(genMuon_index),
                                          genMuon_key);
  }

  /*
   * Generator Muon Propagated
   */

  for(unsigned int genMuonProp_index = 0; genMuonProp_index < GenMuonPropToHO_Etas->size(); genMuonProp_index++){
    std::string genMuon_Prop_key = "genMuon_Prop";
    histogramBuilder.fillCountHistogram(genMuon_Prop_key,weight);
    histogramBuilder.fillEtaPhiHistograms(GenMuonPropToHO_Etas->at(genMuonProp_index),
                                          GenMuonPropToHO_Phis->at(genMuonProp_index),
                                          genMuon_Prop_key);
  }


  /*
   * Higher Level Trigger Single Mu 5
   */

  for(unsigned int hltMu5_index = 0; hltMu5_index < hltMu5_Etas->size(); hltMu5_index++){
    std::string hltMu5_key = "hltMu5";
    histogramBuilder.fillCountHistogram(hltMu5_key,weight);
    histogramBuilder.fillEtaPhiHistograms(hltMu5_Etas->at(hltMu5_index),
					  hltMu5_Phis->at(hltMu5_index),
					  hltMu5_key);
  }

  /*
   * Higher Level Trigger Single Mu 5 Propagated
   */

  for(unsigned int hltProp_index = 0; hltProp_index < hltMu5PropToRPC1_Etas->size(); hltProp_index++){
    std::string hltMu5_Prop_key = "hltMu5_Prop";
    histogramBuilder.fillCountHistogram(hltMu5_Prop_key,weight);
    histogramBuilder.fillEtaPhiHistograms(hltMu5PropToRPC1_Etas->at(hltProp_index),
					  hltMu5PropToRPC1_Phis->at(hltProp_index),
					  hltMu5_Prop_key);
  }


  /*
   * HO Rec Hits Matched to a L1 Muon Barrel, Makes L1MuonBwithMips
   *Select events using indices saved in the L1MuonBIndices Vector
   */
  vector<unsigned int> l1MuonBMipIndices;
  l1MuonBMipIndices.clear();
  for(unsigned int i =0; i< l1MuonBIndices.size(); i++){
    float l1MuonB_eta, l1MuonB_phi, l1MuonB_pt;
    l1MuonB_eta = L1Muon_Etas->at(l1MuonBIndices[i]);
    l1MuonB_phi = L1Muon_Phis->at(l1MuonBIndices[i]);
    l1MuonB_pt = L1Muon_Pts->at(l1MuonBIndices[i]);
    bool hasMip = false; //Only takes one match to be true
    for(unsigned int hoReco_index = 0; hoReco_index < HOReco_Etas->size(); hoReco_index++){
      float horeco_eta, horeco_phi, horeco_energy;
      std::string hoRecoL1MuonBEvent_key = "HO_Reco_L1MuonBEvent";
      horeco_eta = HOReco_Etas->at(hoReco_index);
      horeco_phi = HOReco_Phis->at(hoReco_index);
      horeco_energy = HOReco_Pts->at(hoReco_index);
      
      //histogramBuilder.fillEtaPhiHistograms(horeco_eta,horeco_phi,
      //			    hoRecoL1MuonBEvent_key, weight);
      histogramBuilder.fillEnergyHistograms(horeco_energy,
					    hoRecoL1MuonBEvent_key, weight);
      histogramBuilder.fillDeltaEtaDeltaPhiHistograms(l1MuonB_eta, horeco_eta,
						      l1MuonB_phi, horeco_phi,
						      hoRecoL1MuonBEvent_key, weight);
      
      // Position Match HO Rec Hit to L1MuonB
      
       if(isInsideRCut(RMip_Max, l1MuonB_eta, horeco_eta, l1MuonB_phi, horeco_phi)){
	std::string hoRecoMatchL1MuonB_key = "HOReco_Match_L1MuonB";
	histogramBuilder.fillCountHistogram(hoRecoMatchL1MuonB_key,weight);
	histogramBuilder.fillEnergyHistograms(horeco_energy,
					      hoRecoMatchL1MuonB_key, weight);
	histogramBuilder.fillDeltaEtaDeltaPhiHistograms(l1MuonB_eta, horeco_eta,
							l1MuonB_phi, horeco_phi,
							hoRecoMatchL1MuonB_key, weight);
	
	//Select Rec Hits Above Threshold set Mips
	
	if(horeco_energy > threshold){
	  hasMip = true;
	  std::string hoReco_L1MuonBMip_key = "HOReco_for_L1MuonBMip";
	  histogramBuilder.fillCountHistogram(hoReco_L1MuonBMip_key,weight);
	  histogramBuilder.fillEnergyHistograms(horeco_energy,
						hoReco_L1MuonBMip_key, weight);
	  histogramBuilder.fillDeltaEtaDeltaPhiHistograms(l1MuonB_eta, horeco_eta,
							  l1MuonB_phi, horeco_phi,
							  hoReco_L1MuonBMip_key, weight);
	}
       }
    }
    //Fill histograms for L1MuonB with Mips
    if(hasMip){
      l1MuonBMipIndices.push_back(l1MuonBIndices[i]);
      std::string l1MuonBMip_key = "L1MuonBwithMip";
      histogramBuilder.fillCountHistogram(l1MuonBMip_key,weight);
      histogramBuilder.fillL1MuonPtHistograms(l1MuonB_pt, l1MuonBMip_key, weight);
      histogramBuilder.fillEtaPhiHistograms(l1MuonB_eta,l1MuonB_phi,
					    l1MuonBMip_key, weight);
    }
  }

  /*
   * L1 Muons in the Barrel  Matched to HLT Props
   */

  for(unsigned int i =0; i< l1MuonBIndices.size(); i++){
    float l1MuonB_eta, l1MuonB_phi, l1MuonB_pt;
    l1MuonB_eta = L1Muon_Etas->at(l1MuonBIndices[i]);
    l1MuonB_phi = L1Muon_Phis->at(l1MuonBIndices[i]);
    l1MuonB_pt = L1Muon_Pts->at(l1MuonBIndices[i]);
    bool hasHlt = false; //Only takes one match to be true
   
    for(unsigned int hltProp_index = 0; hltProp_index < hltMu5PropToRPC1_Etas->size(); hltProp_index++){
      float hlt_eta, hlt_phi; 
      hlt_eta = hltMu5PropToRPC1_Etas->at(hltProp_index);
      hlt_phi = hltMu5PropToRPC1_Phis->at(hltProp_index);
      //hlt_pt = hltMu5PropToRPC1_Pts->at(hltProp_index);
      std::string hltPropAndL1MuonB_key = "HLTProp_and_L1MuonBarrel";
      histogramBuilder.fillDeltaEtaDeltaPhiHistograms(l1MuonB_eta, hlt_eta,
						      l1MuonB_phi, hlt_phi,
						      hltPropAndL1MuonB_key, weight);
      if(isInsideRCut(RHlt_Max, l1MuonB_eta, hlt_eta, l1MuonB_phi, hlt_phi)){
	hasHlt = true; //Only takes one match to be true
	std::string hltPropMatchL1MuonB_key = "HLTProp_match_L1MuonBarrel";
	histogramBuilder.fillDeltaEtaDeltaPhiHistograms(l1MuonB_eta, hlt_eta,
							l1MuonB_phi, hlt_phi,
							hltPropMatchL1MuonB_key, weight);
      }
    }
    if(hasHlt){
      std::string l1MuonBMatchHLT_key = "L1MuonBarrelMatchHLT";
      histogramBuilder.fillCountHistogram(l1MuonBMatchHLT_key,weight);
      histogramBuilder.fillL1MuonPtHistograms(l1MuonB_pt, l1MuonBMatchHLT_key, weight);
      histogramBuilder.fillEtaPhiHistograms(l1MuonB_eta,l1MuonB_phi,
					    l1MuonBMatchHLT_key, weight);
    } else {
      std::string l1MuonBNoHLT_key = "L1MuonBarrel_NoHLT";
      histogramBuilder.fillCountHistogram(l1MuonBNoHLT_key,weight);
      histogramBuilder.fillL1MuonPtHistograms(l1MuonB_pt, l1MuonBNoHLT_key, weight);
      histogramBuilder.fillEtaPhiHistograms(l1MuonB_eta,l1MuonB_phi,
                                            l1MuonBNoHLT_key, weight);
    }
  }

  /*
   * L1 Muons in the Barrel With MIP  Matched to HLT Props
   */

  for(unsigned int i =0; i< l1MuonBMipIndices.size(); i++){
    float l1MuonBMip_eta, l1MuonBMip_phi, l1MuonBMip_pt;
    l1MuonBMip_eta = L1Muon_Etas->at(l1MuonBMipIndices[i]);
    l1MuonBMip_phi = L1Muon_Phis->at(l1MuonBMipIndices[i]);
    l1MuonBMip_pt = L1Muon_Pts->at(l1MuonBMipIndices[i]);
    bool hasHlt = false; //Only takes one match to be true
   
    for(unsigned int hltProp_index = 0; hltProp_index < hltMu5PropToRPC1_Etas->size(); hltProp_index++){
      float hlt_eta, hlt_phi; 
      hlt_eta = hltMu5PropToRPC1_Etas->at(hltProp_index);
      hlt_phi = hltMu5PropToRPC1_Phis->at(hltProp_index);
      std::string hltPropAndL1MuonBMip_key = "HLTProp_and_L1MuonwithMipBarrel";
      histogramBuilder.fillDeltaEtaDeltaPhiHistograms(l1MuonBMip_eta, hlt_eta,
						      l1MuonBMip_phi, hlt_phi,
						      hltPropAndL1MuonBMip_key, weight);
      if(isInsideRCut(RHlt_Max, l1MuonBMip_eta, hlt_eta, l1MuonBMip_phi, hlt_phi)){
	hasHlt = true; //Only takes one match to be true
	std::string hltPropMatchL1MuonBMip_key = "HLTProp_match_L1MuonMipBarrel";
	histogramBuilder.fillDeltaEtaDeltaPhiHistograms(l1MuonBMip_eta, hlt_eta,
							l1MuonBMip_phi, hlt_phi,
							hltPropMatchL1MuonBMip_key, weight);
      }
    }
    if(hasHlt){
      std::string l1MuonBMipMatchHLT_key = "L1MuonBarrelwithMipMatchHLT";
      histogramBuilder.fillCountHistogram(l1MuonBMipMatchHLT_key,weight);
      histogramBuilder.fillL1MuonPtHistograms(l1MuonBMip_pt, l1MuonBMipMatchHLT_key, weight);
      histogramBuilder.fillEtaPhiHistograms(l1MuonBMip_eta,l1MuonBMip_phi,
					    l1MuonBMipMatchHLT_key, weight);
    } else {
      std::string l1MuonBMipNoHLT_key = "L1MuonBarrelwithMip_NoHLT";
      histogramBuilder.fillCountHistogram(l1MuonBMipNoHLT_key,weight);
      histogramBuilder.fillL1MuonPtHistograms(l1MuonBMip_pt, l1MuonBMipNoHLT_key, weight);
      histogramBuilder.fillEtaPhiHistograms(l1MuonBMip_eta,l1MuonBMip_phi,
                                            l1MuonBMipNoHLT_key, weight);
    }
  }

  /*
   * Compare L1 Muon Trigger to Propagated Generator Muons
   */

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
    }    
  }
  return kTRUE;
}

void HOMuon_TreeLoop_Plotter::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void HOMuon_TreeLoop_Plotter::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
  
  cout << "working1" << endl;
  outRootFile->cd();
  outRootFile->Write();
  cout << "working2" << endl;
  

}


bool HOMuon_TreeLoop_Plotter::isInsideRCut(float deltaR_Max, float eta1, float eta2, float phi1, float phi2){

  float delta_eta, delta_phi;

  delta_eta = eta1 - eta2;
  delta_phi = histogramBuilder.WrapCheck(phi1,phi2); //Finds difference in phi                              

  //The L1 Muon is compared with all HO Rec Hits above Threshold.                                          
  if(pow(delta_eta,2)+pow(delta_phi,2) <= pow(deltaR_Max,2)) return true;
  return false;
}
