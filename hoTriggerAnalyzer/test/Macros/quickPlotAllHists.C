#include "TFile.h"
#include "TDirectoryFile.h"
#include "iostream.h"
#include "TObject.h"
#include "TList.h"

const char* PATH = "/home/cranelli/HO_Muon/My_AODSIM/CMSSW_6_2_0/src/L1TriggerDPGUpgrade/caloInspector/test/";
const char* fileName="MuonPtHistogram_CrabTest.root";
const char* histDir="histogramBuilder";


void quickPlotAllHists()
{
   
  std::stringstream fileLoc;
  fileLoc << PATH << fileName;
  TFile* inFile = TFile::Open(fileLoc.str().c_str());

  TDirectoryFile* directory = inFile->Get(histDir);
  // TList* list = directory->GetListOfKeys();

  /*
  TObjLink *lnk = list->FirstLink();
  while(lnk) {
    lnk->GetObject();
    lnk = lnk->Next();
  }
  */

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
