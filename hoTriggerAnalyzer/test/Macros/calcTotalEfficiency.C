#include "TFile.h"
#include "TH1F.h"
#include "iostream.h"
#include "TCanvas.h"
#include <string>

//const string newFileName = "L1MuonHistogram_Plus.root";




using namespace std;

void calcTotalEfficiency(string fileLoc = "/home/cranelli/HO_Muon/Spring14/FLAT_QCD/CMSSW_7_0_0/src/"
			 "HoMuonTrigger/hoTriggerAnalyzer/test/Macros/2mil_TreeLoop/L1MuonHistogram.root",
			 string numName = "L1MuonBwithMip_Count", string denName="L1Muon_Barrel_Count"){

  cout << "Numerator Name: " << numName << endl;;
  cout << "Denominator Name: " << denName << endl;
  
  TFile * inFile = TFile::Open(fileLoc.c_str(), "READ");
  TH1F* h1Num = inFile->Get(numName.c_str());
  TH1F* h1Den = inFile->Get(denName.c_str());
  
  cout << "Numerator Count: " << h1Num->GetBinContent(2) << endl;
  cout << "Denominator Count: " << h1Den->GetBinContent(2) << endl;
  TH1F * h1Div = h1Den->Clone();
  h1Div->Reset();
  h1Div->Divide(h1Num, h1Den, 1, 1, "B");
  cout << "Value: " << h1Div->GetBinContent(2) << endl;
  cout << "Upper Error: " << h1Div->GetBinErrorUp(2) << endl;
  cout << "Lower Error: " << h1Div->GetBinErrorLow(2) << endl;
}
