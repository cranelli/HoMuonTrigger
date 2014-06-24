#include "TFile.h"
//#include "TDirectoryFile.h"
#include "TH1.h"
#include "iostream.h"
//#include "TObject.h"
//#include "TList.h"

const char* fileDir = "/home/cranelli/HO_Muon/Fall13/RAW/CMSSW_6_2_9/src/Analysis/hoTriggerAnalyzer/test/";
const char* fileName="L1MuonHistogram.root";
const char* newFileName = "L1MuonHistogram_Plus.root";
const char* histDir="demo/";
const char* histNumName="L1MuonBarrelwithMipHLTMatch_Pt";
const char* histDenName="L1MuonBarrelwithHLTMatch_Pt";



using namespace std;

void divide()
{
   
  std::stringstream inFileLoc, newFileLoc;
  inFileLoc << fileDir << fileName;
  TFile* inFile = TFile::Open(inFileLoc.str().c_str());
  newFileLoc << fileDir << newFileName;
  TFile* newFile = TFile::Open(newFileLoc.str().c_str(),"RECREATE");
  std::stringstream histNumLoc, histDenLoc;
  histNumLoc <<histDir << histNumName;
  histDenLoc <<histDir << histDenName;
  
  TH1F* h1Num = inFile->Get(histNumLoc.str().c_str());
  TH1F* h1Den = inFile->Get(histDenLoc.str().c_str());
  
  //cout <<h1Num->GetNbinsX() << endl;
  //cout <<h1Den->GetNbinsX() << endl;

  //h1Num->Draw();
  //h1Den->Draw();
  //inFile->Draw("h1Num");

  //hNum->Divide(hDen);

  float variableBinArray[] = {0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,12,14,16,18,20,25,30,35,40,45,50,60,70,80,100,120,140,180};

  TH1F* h1Div = new TH1F("L1MuonBarrel_Pt_Signal", "L1MuonBarrel_Pt_Signal",33,variableBinArray);

  h1Div->Divide(h1Num, h1Den);
  
  //int num_bins = (int) h1->GetNbinX();
  //std::cout << h1->GetNbinsX() << std::endl;
  

  newFile->Write();

  inFile->Close();

}
