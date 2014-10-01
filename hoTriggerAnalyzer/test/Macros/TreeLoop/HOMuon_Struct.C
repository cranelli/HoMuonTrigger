/*
 * Created  by Christopher Anelli 3.10.2014
 * Continue using an All SiPM Detector.
 * Back an Front Surfaces set by the caloGeometry
 * are now propagated to.
 */

#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "iostream"
#include <vector>
#include <TROOT.h>

//MakeSelector
//#include "HoMuonTrigger/hoTriggerAnalyzer/test/Macros/HOMuon_TreeLoop_Struct.h"
//MakeClass
//#include "HoMuonTrigger/hoTriggerAnalyzer/test/Macros/HOMuon_Loop_Struct.h"

/*
 *Definitions (easy to modify)
 */

#define ROOTFILE_DIR "/data/users/cranelli/HOL1Muon/Trees/"
#define VERSION "Version_5_1/"
#define ROOTFILE_NAME "HOMuonTree_Test.root"
#define TREE_LOC "demo/ho_muon_tree"
#define MAX_EVENTS 20000000  //Default Max

void HOMuon_Struct(){  
  
  std::stringstream rootFileLoc;
  rootFileLoc << ROOTFILE_DIR << VERSION << ROOTFILE_NAME;
  TFile* inFile = TFile::Open(rootFileLoc.str().c_str());
  
  //TFile* inFile = TFile::Open(ROOTFILE_LOC);
  //std::cout << inFile << std::endl;
  TTree* T = (TTree*)inFile->Get(TREE_LOC);
  cout<< T << endl;
  //gRoot->ProcessLine(".L HOMuon_TreeLoop_FrontBack_Struct.C+");
  T->Process("HOMuon_TreeLoop_Struct.C+", "",MAX_EVENTS);



  //For MakeClass
  //gROOT->ProcessLine(".L HOMuon_Loop_Struct.C");
  //gROOT->ProcessLine(".L ../../plugins/HistogramBuilder.cc");
  //HOMuon_Loop_Struct looper;
  //looper.Loop();
  //looper.Show(16);
  //cout << "working" << endl;
};
  
