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
const char* histName="L1MuonBarrel_Pt";
const int num_bins=33;

using namespace std;

void divideByBinWidth()
{
   
  std::stringstream inFileLoc, newFileLoc;
  inFileLoc << fileDir << fileName;
  TFile* inFile = TFile::Open(inFileLoc.str().c_str(),"READ");
  newFileLoc << fileDir << newFileName;
  TFile* newFile = TFile::Open(newFileLoc.str().c_str(),"RECREATE");
  std::stringstream histLoc;
  histLoc <<histDir << histName;
  TH1F* h1 = inFile->Get(histLoc.str().c_str());
  //int num_bins = (int) h1->GetNbinX();
  //std::cout << h1->GetNbinsX() << std::endl;


  float variableBinArray[] = {0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,12,14,16,18,20,25,30,35,40,45,50,60,70,80,100,120,140,180};
  
  TH1F* h1Norm = new TH1F("L1MuonBarrel_Pt_Norm", "L1MuonBarrel_Pt_Norm",33,variableBinArray);
  
  for(int i=0; i<=num_bins; i++){
    float binContent = h1->GetBinContent(i);
    float binWidth = h1->GetBinWidth(i);
    float binCenter = h1->GetBinCenter(i);
    cout << i << " , " << binCenter << " , " << binContent << " ," << binWidth << endl; 
    //cout << binContent << endl;
    //cout << binWidth << endl;
    //cout <<binCenter << endl;
    h1Norm->Fill(binCenter, binContent/binWidth); 
  }
  newFile->Write();

  //std::cout << h1Norm->GetNbinsX() << std::endl;

  //TDirectoryFile* directory = inFile->Get(histDir);
  // TList* list = directory->GetListOfKeys();

  /*
  TObjLink *lnk = list->FirstLink();
  while(lnk) {
    lnk->GetObject();
    lnk = lnk->Next();
  }
  */

  /*
  TObject* key;
  TObject* obj;
  
  TIter next(directory->GetListOfKeys());
  while(key = next()){
    std::cout << key->GetName() << std::endl;
    obj = directory->Get(key->GetName());
    std::stringstream outputName;
    outputName << key->GetName() <<".png"; 
    obj->SaveAs(outputName.str().c_str());
  }
  */
  //
  
  /*
  auto bObject = list->First();
  auto eObject = list->Last();

  for(; bObject!=eObject; ++bObject){
    std::cout<< bObject.GetName() << std::endl;
  }
  */

  inFile->Close();
  
  

}
