#include "TFile.h"
//#include "TDirectoryFile.h"
#include "TH1.h"
#include "iostream.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
//#include "TObject.h"
//#include "TList.h"


const char* fileDir = "/home/cranelli/HO_Muon/Spring14/FLAT_QCD/CMSSW_7_0_0/src/HoMuonTrigger/hoTriggerAnalyzer/test/Macros/TreeLoop";
const char* ptHat = "";

//const char* fileDir = "/data/users/cranelli/HOL1Muon/HOL1Muon_Histograms/QCD/Version_1_3/";
//#const char* ptHat = "WeightedMerger"; 
//"/home/cranelli/HO_Muon/Fall13/RAW/CMSSW_6_2_9/src/Analysis/hoTriggerAnalyzer/test/";
const char* fileName="L1MuonHistogram.root";
const char* newFileName = "L1MuonHistogram_Plus.root";
const char* histDir="";

using namespace std;

TFile * inFile;
TFile * newFile;

/*
 * Using a series of helper Macros,
 * This produces the secondary Plots of interests.
 * Made using the primary histograms.
 */

void secondaryPlots(){
  std::stringstream inFileLoc, newFileLoc;
  inFileLoc << fileDir << ptHat<< "/"<< fileName;
  inFile = TFile::Open(inFileLoc.str().c_str(),"READ");
  newFileLoc << fileDir << ptHat << "/" << newFileName;
  newFile = TFile::Open(newFileLoc.str().c_str(),"RECREATE");


  // Normal L1 Muon Barrel_Pt
  TH1F* h1L1Norm = divideByBinWidth("L1Muon_Barrel_Pt");
  SetAxises(h1L1Norm, "Pt", "dN/dx");
  //h1Norm->Draw("HIST");
  h1L1Norm->Write();
  
  TH1F* h1L1MipNorm = divideByBinWidth("L1MuonBwithMip_Pt");
  SetAxises(h1L1MipNorm, "Pt", "dN/dx");
  h1L1MipNorm->Write();

  TH1F* h1L1HLTNorm = divideByBinWidth("L1MuonBarrelwithMipMatchHLT_Pt");
  SetAxises(h1L1HLTNorm, "Pt", "dN/dx");
  h1L1HLTNorm->Write();
  
  
  /*
  //Signal 
  std::stringstream Signal;
  Signal  <<"Signal" << endl;
  //Signal << ptHat <<"_Signal" << endl;
  divide("L1MuonBarrelwithMipMatchHLT_Pt", "L1MuonBarrelMatchHLT_Pt", Signal.str().c_str());
  */

  //Background
  TH1F* back_num = subtract("L1MuonBwithMip_Pt", "L1MuonBarrelwithMipMatchHLT_Pt", "L1MuonBarrelwithMipNOHLTMatch_Pt");
  TH1F* back_den = subtract("L1Muon_Barrel_Pt", "L1MuonBarrelMatchHLT_Pt", "L1MuonBarrelwithoutHLTMatch_Pt");
  
  std::stringstream Background;
  Background <<"Background" << endl;
  divide(back_num, back_den, Background.str().c_str());
  

  //newFile->Write();
  inFile->Close();
}

/*
 * Creates a new Histogram.  By taking the entries of the orginal histogram,
 * and dividing them by their bin width.
 */

TH1F* divideByBinWidth(char * histName){  
  cout << histName << endl;
  std::stringstream histLoc, normName;
  histLoc <<histDir << histName;
  TH1F* h1 = inFile->Get(histLoc.str().c_str());
  
  normName << histName <<"_Norm";
  
  TH1F* h1Norm = (TH1F*) h1->Clone(normName.str().c_str());
  h1Norm->Reset();
  h1Norm->SetTitle(normName.str().c_str());
  h1Norm->GetYaxis()->SetTitle("#frac{dN}{dx}");

  h1->Sumw2();
  int num_bins = (int) h1Norm->GetNbinsX();
  for(int i=0; i<=num_bins; i++){
    float binContent = h1->GetBinContent(i);
    float binWidth = h1->GetBinWidth(i);
    float binCenter = h1->GetBinCenter(i);
    //cout << i << " , " << binCenter << " , " << binContent << " ," << binWidth << endl; 
    h1Norm->Fill(binCenter, binContent/binWidth);
    h1Norm->SetBinError(i, h1->GetBinError(i)/binWidth);  //Manually set
    //Check Errors
    /*
    cout << i << " " << h1->GetBinError(i) << endl;
    cout << h1Norm->GetBinError(i) << " " << sqrt(h1Norm.GetBinContent(i))<< " " 
	 << h1->GetBinError(i)/binWidth << endl; 
    */
    //cout << i << ", " << h1Norm_test->GetBinCenter(i) << endl;
  }
  return h1Norm;
  //h1Norm->Write();
  
}


/*
 * Wrapper For divide, that takes the names as strings
 */

void divide(char * histNumName, char *histDenName, char *quotientName)
{
  
  std::stringstream histNumLoc, histDenLoc;
  histNumLoc <<histDir << histNumName;
  histDenLoc <<histDir << histDenName;

  TH1F* h1Num = inFile->Get(histNumLoc.str().c_str());
  TH1F* h1Den = inFile->Get(histDenLoc.str().c_str());

  divide(h1Num, h1Den, quotientName);
}

/*
 * Divides the first histogram by the second histogram.
 */

void divide(TH1F* h1Num, TH1F* h1Den, char *quotientName)
{
  examineBinErrors(h1Num);
  examineBinErrors(h1Den);
  TH1F* h1Div = (TH1F*) h1Num->Clone(quotientName);
  h1Div->Reset();
  h1Div->SetTitle(quotientName);
  h1Div->GetYaxis()->SetTitle(quotientName);
  
  //h1Div->Sumw2()
  h1Num->Sumw2();
  h1Den->Sumw2();
  h1Div->Divide(h1Num, h1Den);
    
  cout << "Binomial Errors" << endl;
  binomialErrors = new TGraphAsymmErrors(h1Num, h1Den);
  binomialErrors.Print();

  c1 = new TCanvas("c1","A Simple Graph with assymetric error bars",200,10,700,500);  
  binomialErrors->SetTitle(quotientName);
  binomialErrors->GetXaxis()->SetTitle("Pt (GeV)");
  binomialErrors->GetYaxis()->SetTitle("Ratio");
  binomialErrors->SetMarkerColor(4);
  binomialErrors->SetMarkerStyle(21);
  binomialErrors->Draw("AP");
  //binomialErrors->Write("AP");
  
  //cout << "Quotient Content and Errors" << endl;
  //examineBinErrors(h1Div);
  
  h1Div->Write();
  
}


TH1F* subtract(char * h1FirstName, char *h1SecondName, char* difName)
{

  
  std::stringstream h1FirstLoc, h1SecondLoc;
  h1FirstLoc <<histDir << h1FirstName;
  h1SecondLoc <<histDir << h1SecondName;

  TH1F* h1First = inFile->Get(h1FirstLoc.str().c_str());
  TH1F* h1Second = inFile->Get(h1SecondLoc.str().c_str());

  TH1F* h1Dif = (TH1F*) h1First->Clone(difName);
  h1Dif->Reset();
  h1Dif->SetTitle(difName);
 
  h1Dif->Add(h1First, h1Second, 1.0, -1.0); //Subtraction

  h1Dif->Write();

  return h1Dif;

}

//Examine The Errors
void examineBinErrors(TH1F* h1){
  
  int num_bins = (int) h1->GetNbinsX();
  for(int i=0; i<=num_bins; i++){
    float binContent = h1->GetBinContent(i);
    float binError = h1->GetBinError(i);
    cout << binContent << " " << binError << " " << sqrt(binContent) << endl;
  }
}

void SetAxises(TH1F * h1, std::string xTitle, std::string yTitle){
  // X Axis
  h1->GetXaxis()->SetTitle(xTitle.c_str());
  h1->GetXaxis()->SetTitleFont(42);
  h1->GetXaxis()->SetTitleSize(0.05);
  h1->GetXaxis()->SetTitleOffset(0.9);

  // Y Axis 
  h1->GetYaxis()->SetTitle(yTitle.c_str());
  h1->GetYaxis()->SetTitleFont(42);
  h1->GetYaxis()->SetTitleSize(0.05);
  h1->GetYaxis()->SetTitleOffset(1.0);

  h1->SetLineWidth(2);
  h1->SetLineColor(kBlue);
}

