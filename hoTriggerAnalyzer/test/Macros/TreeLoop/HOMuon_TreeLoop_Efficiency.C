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
   
   efficiencyDen_count=0;
   for(int i_Tcut=0; i_Tcut < num_Tcuts; i_Tcut++){
     thresholds[i_Tcut]=(double)0.05+0.05*i_Tcut;
     efficiencyNum_count[i_Tcut]=0;
     //accepted_count[i_Tcut]=0;
     //fake_count[i_Tcut]=0;
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
   * Weighting Information                                                      
   */
  double weight = Generator_Weights;
  
  /*  
   * Select L1 Muons inside the Barrel
   */
  vector<unsigned int> l1MuonBIndices;
  l1MuonBIndices.clear();
  for(unsigned int l1MuonB_index = 0; l1MuonB_index < L1Muon_Etas->size(); l1MuonB_index++){
    if(fabs(L1Muon_Etas->at(l1MuonB_index)) > barrel_eta) continue; //Must be in barrel
    l1MuonBIndices.push_back(l1MuonB_index);
  }

  /*
   * L1 Muons in the Barrel Matched to an HLT
   */

  for(unsigned int i =0; i< l1MuonBIndices.size(); i++){
    float l1MuonB_eta, l1MuonB_phi;
    l1MuonB_eta = L1Muon_Etas->at(l1MuonBIndices[i]);
    l1MuonB_phi = L1Muon_Phis->at(l1MuonBIndices[i]);
    bool hasHlt = false; //Only takes one match to be true
    for(unsigned int hltProp_index = 0; hltProp_index < hltMu5PropToRPC1_Etas->size(); hltProp_index++){
      float hlt_eta, hlt_phi;
      hlt_eta = hltMu5PropToRPC1_Etas->at(hltProp_index);
      hlt_phi = hltMu5PropToRPC1_Phis->at(hltProp_index);
      if(isInsideRCut(RHlt_Max, l1MuonB_eta, hlt_eta, l1MuonB_phi, hlt_phi)){
        hasHlt = true; //Only takes one match to be true                                                    
      }
    }
    if(hasHlt){
      efficiencyDen_count+=weight;
    }
  }

  /*
   * Select L1 Barrel Muons, with a MIP
   */

  vector<unsigned int> l1MuonBMipIndices[num_Tcuts];
  for(unsigned int i =0; i< l1MuonBIndices.size(); i++){
    float l1MuonB_eta, l1MuonB_phi;
    l1MuonB_eta = L1Muon_Etas->at(l1MuonBIndices[i]);
    l1MuonB_phi = L1Muon_Phis->at(l1MuonBIndices[i]);
    for(int i_Tcut=0; i_Tcut< num_Tcuts; i_Tcut++){
      l1MuonBMipIndices[i_Tcut].clear();
      bool hasMip = false; //Only takes one match to be true
      for(unsigned int hoReco_index = 0; hoReco_index < HOReco_Etas->size(); hoReco_index++){
	float horeco_eta, horeco_phi, horeco_energy;
	std::string hoRecoL1MuonBEvent_key = "HO_Reco_L1MuonBEvent";
	horeco_eta = HOReco_Etas->at(hoReco_index);
	horeco_phi = HOReco_Phis->at(hoReco_index);
	horeco_energy = HOreco_Energies->at(hoReco_index);
	//Position Match HO Rec Hit to L1MuonB 
	if(isInsideRCut(RMip_Max, l1MuonB_eta, horeco_eta, l1MuonB_phi, horeco_phi)){
	  //Select Rec Hits Above Threshold set Mips
	  if(horeco_energy > thresholds[i_Tcut]){
	    hasMip = true;
	  }
	}
      }
      if(hasMip){
	l1MuonBMipIndices[i_Tcut].push_back(l1MuonBIndices[i]);
      }
    }
  }

  /*                                                        
   * L1 Muons in the Barrel With MIP  Matched to HLT Props                                                   
   */
  
  for(int i_Tcut=0; i_Tcut < num_Tcuts; i_Tcut++){
    // if(l1MuonBMipIndices[i_Tcut].size()!=0) cout << "Number of Non Zero Mip Indices" << l1MuonBMipIndices[i_Tcut].size() << endl;
    for(unsigned int i =0; i< l1MuonBMipIndices[i_Tcut].size(); i++){
      float l1MuonBMip_eta, l1MuonBMip_phi;
      //cout << i_Tcut << ":" << l1MuonBMipIndices[i_Tcut][0] << endl;
      l1MuonBMip_eta = L1Muon_Etas->at(l1MuonBMipIndices[i_Tcut][i]);
      l1MuonBMip_phi = L1Muon_Phis->at(l1MuonBMipIndices[i_Tcut][i]);
      
      bool hasHlt = false; //Only takes one match to be true  
      for(unsigned int hltProp_index = 0; hltProp_index < hltMu5PropToRPC1_Etas->size(); hltProp_index++){
	float hlt_eta, hlt_phi;
	hlt_eta = hltMu5PropToRPC1_Etas->at(hltProp_index);
	hlt_phi = hltMu5PropToRPC1_Phis->at(hltProp_index);
	if(isInsideRCut(RHlt_Max, l1MuonBMip_eta, hlt_eta, l1MuonBMip_phi, hlt_phi)){
	  hasHlt = true; //Only takes one match to be true 
	}
	}
      if(hasHlt){
	efficiencyNum_count[i_Tcut]+=weight;
      }
      
    }
  }

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
  cout << "Efficiency Denominator Count:" << efficiencyDen_count << endl;
  cout << "Efficiency Numerator Count:" << endl;
  for(int i_Tcut; i_Tcut < num_Tcuts; i_Tcut++){
    cout << "     Threshold " << thresholds[i_Tcut] << ": " << efficiencyNum_count[i_Tcut] << endl;
  }
}

bool HOMuon_TreeLoop_Efficiency::isInsideRCut(float deltaR_Max, float eta1, float eta2, float phi1, float phi2){
  
  float delta_eta, delta_phi;

  delta_eta = eta1 - eta2;
  delta_phi = histogramBuilder.WrapCheck(phi1,phi2); //Finds difference in phi                                        

  //The L1 Muon is compared with all HO Rec Hits above Threshold.                                                     
  if(pow(delta_eta,2)+pow(delta_phi,2) <= pow(deltaR_Max,2)) return true;
  return false;
}
