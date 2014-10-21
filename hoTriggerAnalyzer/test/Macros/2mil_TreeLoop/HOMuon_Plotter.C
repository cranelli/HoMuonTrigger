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

/*
 *Definitions (easy to modify)
 */

#define ROOTFILE_DIR "/data/users/cranelli/HOL1Muon/Trees/"
#define VERSION "Version_4_2/"
#define ROOTFILE_NAME "HOMuonTree.root"
#define TREE_LOC "demo/ho_muon_tree"
#define MAX_EVENTS 1000000000  //Default Max

void HOMuon_Plotter(){  
  
  std::stringstream rootFileLoc;
  rootFileLoc << ROOTFILE_DIR << VERSION << ROOTFILE_NAME;
  TFile* inFile = TFile::Open(rootFileLoc.str().c_str());

  TTree* T = (TTree*)inFile->Get(TREE_LOC);
  cout<< T << endl;
  //gRoot->ProcessLine(".L HOMuon_TreeLoop_FrontBack_Plotter.C+");
  T->Process("HOMuon_TreeLoop_Plotter.C+", "",MAX_EVENTS);

};
  
