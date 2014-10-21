//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Jul  6 14:52:42 2014 by ROOT version 5.34/10
// from TTree ho_muon_tree/Generator, L1, L1Muon,HLT, HLT Trigger Object, and HO Reco Data
// found on file: ../HOMuonTree.root
//////////////////////////////////////////////////////////

#ifndef HOMuon_TreeLoop_Plotter_h
#define HOMuon_TreeLoop_Plotter_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <set>
#include "HistogramBuilderTwo.cc"

// Header file for the classes stored in the TTree if any.
#include <vector>
//#include <map>

//Define Constants for the Tree Analyzer
static const float barrel_eta = 1.3;
static const float threshold = 0.3;
static const float RMip_Max = 0.2;
static const float RHlt_Max = 0.4;

// Fixed size dimensions of array or collections stored in the TTree if any.

class HOMuon_TreeLoop_Plotter : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Double_t        Generator_Weights;
   vector<int>     *Generator_pdgIds;
   vector<float>   *Generator_Etas;
   vector<float>   *Generator_Phis;
   vector<float>   *Generator_Pts;
   vector<int>     *GenMuonPropToHO_pdgIds;
   vector<float>   *GenMuonPropToHO_Etas;
   vector<float>   *GenMuonPropToHO_Phis;
   vector<float>   *GenMuonPropToHO_Pts;
   vector<int>     *GenMuonPropToRPC1_pdgIds;
   vector<float>   *GenMuonPropToRPC1_Etas;
   vector<float>   *GenMuonPropToRPC1_Phis;
   vector<float>   *GenMuonPropToRPC1_Pts;
   vector<float>   *L1Muon_Etas;
   vector<float>   *L1Muon_Phis;
   vector<float>   *L1Muon_Pts;
   vector<float>   *HOReco_Etas;
   vector<float>   *HOReco_Phis;
   vector<float>   *HOReco_Pts;
   vector<float>   *hltMu5_Etas;
   vector<float>   *hltMu5_Phis;
   vector<float>   *hltMu5_Pts;
   vector<float>   *hltMu5PropToRPC1_Etas;
   vector<float>   *hltMu5PropToRPC1_Phis;
   vector<float>   *hltMu5PropToRPC1_Pts;

   // List of branches
   TBranch        *b_Generator_Weights;   //!
   TBranch        *b_Generator_pdgIds;   //!
   TBranch        *b_Generator_Etas;   //!
   TBranch        *b_Generator_Phis;   //!
   TBranch        *b_Generator_Pts;   //!
   TBranch        *b_GenMuonPropToHO_pdgIds;   //!
   TBranch        *b_GenMuonPropToHO_Etas;   //!
   TBranch        *b_GenMuonPropToHO_Phis;   //!
   TBranch        *b_GenMuonPropToHO_Pts;   //!
   TBranch        *b_GenMuonPropToRPC1_pdgIds;   //!
   TBranch        *b_GenMuonPropToRPC1_Etas;   //!
   TBranch        *b_GenMuonPropToRPC1_Phis;   //!
   TBranch        *b_GenMuonPropToRPC1_Pts;   //!
   TBranch        *b_L1Muon_Etas;   //!
   TBranch        *b_L1Muon_Phis;   //!
   TBranch        *b_L1Muon_Pts;   //!
   TBranch        *b_HOReco_Etas;   //!
   TBranch        *b_HOReco_Phis;   //!
   TBranch        *b_HOReco_Pts;   //!
   TBranch        *b_hltMu5_Etas;   //!
   TBranch        *b_hltMu5_Phis;   //!
   TBranch        *b_hltMu5_Pts;   //!
   TBranch        *b_hltMu5PropToRPC1_Etas;   //!
   TBranch        *b_hltMu5PropToRPC1_Phis;   //!
   TBranch        *b_hltMu5PropToRPC1_Pts;   //!

   HOMuon_TreeLoop_Plotter(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~HOMuon_TreeLoop_Plotter() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   //My Stuff
   TFile * outRootFile;
   HistogramBuilderTwo histogramBuilder;

   /*
    * Helper Functions
    */
   bool isInsideRCut(float deltaR_Max, float eta1, float eta2, float phi1, float phi2);

   ClassDef(HOMuon_TreeLoop_Plotter,0);
};

#endif

#ifdef HOMuon_TreeLoop_Plotter_cxx
void HOMuon_TreeLoop_Plotter::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Generator_pdgIds = 0;
   Generator_Etas = 0;
   Generator_Phis = 0;
   Generator_Pts = 0;
   GenMuonPropToHO_pdgIds = 0;
   GenMuonPropToHO_Etas = 0;
   GenMuonPropToHO_Phis = 0;
   GenMuonPropToHO_Pts = 0;
   GenMuonPropToRPC1_pdgIds = 0;
   GenMuonPropToRPC1_Etas = 0;
   GenMuonPropToRPC1_Phis = 0;
   GenMuonPropToRPC1_Pts = 0;
   L1Muon_Etas = 0;
   L1Muon_Phis = 0;
   L1Muon_Pts = 0;
   HOReco_Etas = 0;
   HOReco_Phis = 0;
   HOReco_Pts = 0;
   hltMu5_Etas = 0;
   hltMu5_Phis = 0;
   hltMu5_Pts = 0;
   hltMu5PropToRPC1_Etas = 0;
   hltMu5PropToRPC1_Phis = 0;
   hltMu5PropToRPC1_Pts = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Generator_Weights", &Generator_Weights, &b_Generator_Weights);
   fChain->SetBranchAddress("Generator_pdgIds", &Generator_pdgIds, &b_Generator_pdgIds);
   fChain->SetBranchAddress("Generator_Etas", &Generator_Etas, &b_Generator_Etas);
   fChain->SetBranchAddress("Generator_Phis", &Generator_Phis, &b_Generator_Phis);
   fChain->SetBranchAddress("Generator_Pts", &Generator_Pts, &b_Generator_Pts);
   fChain->SetBranchAddress("GenMuonPropToHO_pdgIds", &GenMuonPropToHO_pdgIds, &b_GenMuonPropToHO_pdgIds);
   fChain->SetBranchAddress("GenMuonPropToHO_Etas", &GenMuonPropToHO_Etas, &b_GenMuonPropToHO_Etas);
   fChain->SetBranchAddress("GenMuonPropToHO_Phis", &GenMuonPropToHO_Phis, &b_GenMuonPropToHO_Phis);
   fChain->SetBranchAddress("GenMuonPropToHO_Pts", &GenMuonPropToHO_Pts, &b_GenMuonPropToHO_Pts);
   fChain->SetBranchAddress("GenMuonPropToRPC1_pdgIds", &GenMuonPropToRPC1_pdgIds, &b_GenMuonPropToRPC1_pdgIds);
   fChain->SetBranchAddress("GenMuonPropToRPC1_Etas", &GenMuonPropToRPC1_Etas, &b_GenMuonPropToRPC1_Etas);
   fChain->SetBranchAddress("GenMuonPropToRPC1_Phis", &GenMuonPropToRPC1_Phis, &b_GenMuonPropToRPC1_Phis);
   fChain->SetBranchAddress("GenMuonPropToRPC1_Pts", &GenMuonPropToRPC1_Pts, &b_GenMuonPropToRPC1_Pts);
   fChain->SetBranchAddress("L1Muon_Etas", &L1Muon_Etas, &b_L1Muon_Etas);
   fChain->SetBranchAddress("L1Muon_Phis", &L1Muon_Phis, &b_L1Muon_Phis);
   fChain->SetBranchAddress("L1Muon_Pts", &L1Muon_Pts, &b_L1Muon_Pts);
   fChain->SetBranchAddress("HOReco_Etas", &HOReco_Etas, &b_HOReco_Etas);
   fChain->SetBranchAddress("HOReco_Phis", &HOReco_Phis, &b_HOReco_Phis);
   fChain->SetBranchAddress("HOReco_Pts", &HOReco_Pts, &b_HOReco_Pts);
   fChain->SetBranchAddress("hltMu5_Etas", &hltMu5_Etas, &b_hltMu5_Etas);
   fChain->SetBranchAddress("hltMu5_Phis", &hltMu5_Phis, &b_hltMu5_Phis);
   fChain->SetBranchAddress("hltMu5_Pts", &hltMu5_Pts, &b_hltMu5_Pts);
   fChain->SetBranchAddress("hltMu5PropToRPC1_Etas", &hltMu5PropToRPC1_Etas, &b_hltMu5PropToRPC1_Etas);
   fChain->SetBranchAddress("hltMu5PropToRPC1_Phis", &hltMu5PropToRPC1_Phis, &b_hltMu5PropToRPC1_Phis);
   fChain->SetBranchAddress("hltMu5PropToRPC1_Pts", &hltMu5PropToRPC1_Pts, &b_hltMu5PropToRPC1_Pts);
}

Bool_t HOMuon_TreeLoop_Plotter::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef HOMuon_TreeLoop_Plotter_cxx
