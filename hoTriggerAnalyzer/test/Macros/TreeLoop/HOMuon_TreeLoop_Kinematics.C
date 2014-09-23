#define HOMuon_TreeLoop_Kinematics_cxx
// The class definition in HOMuon_TreeLoop_Kinematics.h has been generated automatically
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
// Root > T->Process("HOMuon_TreeLoop_Kinematics.C")
// Root > T->Process("HOMuon_TreeLoop_Kinematics.C","some options")
// Root > T->Process("HOMuon_TreeLoop_Kinematics.C+")
//

#include "HOMuon_TreeLoop_Kinematics.h"
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



void HOMuon_TreeLoop_Kinematics::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   outRootFile = TFile::Open("L1Muon_KinematicHistograms.root", "RECREATE");

}

void HOMuon_TreeLoop_Kinematics::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t HOMuon_TreeLoop_Kinematics::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either HOMuon_TreeLoop_Kinematics::GetEntry() or TBranch::GetEntry()
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
  
  /*
   * Weighting Information
   */  
  double weight = Generator_Weights;

  std::string countAllEvents_key = "All_Events";
  std::string countAllEventsW_key = "All_Events_Weighted";
  histogramBuilder.fillCountHistogram(countAllEvents_key);
  histogramBuilder.fillCountHistogram(countAllEventsW_key, weight);
  histogramBuilder.fillWeightHistograms(weight, countAllEvents_key);
  histogramBuilder.fillWeightHistograms(weight, countAllEventsW_key, weight);

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
   * HO Rec Hit
   */
  for(unsigned int hoReco_index = 0; hoReco_index < HOReco_Etas->size(); hoReco_index++){
    std::string hoReco_key = "HO_Reco";
    histogramBuilder.fillEtaPhiHistograms(HOReco_Etas->at(hoReco_index),
                                          HOReco_Phis->at(hoReco_index),
                                          hoReco_key, weight);
    histogramBuilder.fillEnergyHistograms(HOReco_Energies->at(hoReco_index),
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
    histogramBuilder.fillPtHistograms(hltMu5_Pts->at(hltMu5_index),
				      hltMu5_key, weight);
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

  return kTRUE;
}

void HOMuon_TreeLoop_Kinematics::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void HOMuon_TreeLoop_Kinematics::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
  
  cout << "working1" << endl;
  outRootFile->cd();
  outRootFile->Write();
  cout << "working2" << endl;
  

}


bool HOMuon_TreeLoop_Kinematics::isInsideRCut(float deltaR_Max, float eta1, float eta2, float phi1, float phi2){

  float delta_eta, delta_phi;

  delta_eta = eta1 - eta2;
  delta_phi = histogramBuilder.WrapCheck(phi1,phi2); //Finds difference in phi                              

  //The L1 Muon is compared with all HO Rec Hits above Threshold.                                          
  if(pow(delta_eta,2)+pow(delta_phi,2) <= pow(deltaR_Max,2)) return true;
  return false;
}
