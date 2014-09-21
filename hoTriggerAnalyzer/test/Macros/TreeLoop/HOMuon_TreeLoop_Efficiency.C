#define HOMuon_TreeLoop_Efficiency_cxx
// The class definition in HOMuon_TreeLoop_Efficiency.h has been generated automatically
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
// Root > T->Process("HOMuon_TreeLoop_Efficiency.C")
// Root > T->Process("HOMuon_TreeLoop_Efficiency.C","some options")
// Root > T->Process("HOMuon_TreeLoop_Efficiency.C+")
//

#include "HOMuon_TreeLoop_Efficiency.h"
#include <TH2.h>
#include <TStyle.h>
#include "TLorentzVector.h"
#include "HistogramBuilderTwo.h"

void HOMuon_TreeLoop_Efficiency::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   
   for(int i_Tcut=0; i_Tcut < num_Tcuts; i_Tcut++){
     efficiencyNum_count[i_Tcut]=0;
     efficiencyDen_count[i_Tcut]=0;
     //accepted_count[i_Tcut]=0;
     fake_count[i_Tcut]=0;
   }
}

void HOMuon_TreeLoop_Efficiency::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t HOMuon_TreeLoop_Efficiency::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either HOMuon_TreeLoop_Efficiency::GetEntry() or TBranch::GetEntry()
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
  if(entry%100 == 0) cout << entry << endl;

  
  /*  
  // Select Only L1 Muons inside the Barrel
  vector<unsigned int> l1MuonBIndices;
  l1MuonBIndices.clear();
  for(unsigned int l1MuonB_index = 0; l1MuonB_index < L1Muon_Etas->size(); l1MuonB_index++){
    if(fabs(L1Muon_Etas->at(l1MuonB_index)) > barrel_eta) continue; //Must be in barrel
    l1MuonBIndices.push_back(l1MuonB_index);
  }

  for(unsigned int i =0; i< l1MuonBIndices.size(); i++){
    TLoretnzVector l1MuonB;  //Kinematic Information is Stored 
    l1MuonB_eta = L1Muon_Etas->at(l1MuonBIndices[i]);
    l1MuonB_phi = L1Muon_Phis->at(l1MuonBIndices[i]);
    l1MuonB_pt = L1Muon_Pts->at(l1MuonBIndices[i]);
    for(unsigned int hoReco_index = 0; hoReco_index < HOReco_Etas->size(); hoReco_index++){
      //Position Match HO to L1MuonB
      horeco_eta = HOReco_Etas->at(hoReco_index);
      horeco_phi = HOReco_Phis->at(hoReco_index);
      horeco_energy = HOReco_Energies->at(hoReco_index);
      if(isInsideRCut(RMip_Max, l1MuonB_eta, horeco_eta, l1MuonB_phi, horeco_phi)){
      }
    }
  }
    
  */

  return kTRUE;
}

void HOMuon_TreeLoop_Efficiency::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void HOMuon_TreeLoop_Efficiency::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

}

bool HOMuon_TreeLoop_Efficiency::isInsideRCut(float deltaR_Max, float eta1, float eta2, float phi1, float phi2){
  
  float delta_eta, delta_phi;

  delta_eta = eta1 - eta2;
  delta_phi = histogramBuilder.WrapCheck(phi1,phi2); //Finds difference in phi                                        

  //The L1 Muon is compared with all HO Rec Hits above Threshold.                                                     
  if(pow(delta_eta,2)+pow(delta_phi,2) <= pow(deltaR_Max,2)) return true;
  return false;
}
