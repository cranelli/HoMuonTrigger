#include "TFile.h"
#include "TH1F.h"
#include "iostream.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include <string>
#include <vector>


const string fileDir = "/home/cranelli/HO_Muon/Spring14/FLAT_QCD/CMSSW_7_0_0/src/HoMuonTrigger/hoTriggerAnalyzer/test/Macros/TreeLoop/";
const string fileName="new_L1MuonHistogram.root";
const string newFileName = "L1MuonHistogram_Plus.root";
//const string histDir="";

using namespace std;

TFile * inFile = TFile::Open((fileDir+fileName).c_str(), "READ");
TFile * newFile =TFile::Open((fileDir+newFileName).c_str(),"RECREATE");

/*
 * Get Histograms
 */

//L1 Muon Pt Histograms
const int NUM_L1MUONPT_HISTS = 6;
TH1F* h1L1MuonB_l1muon_Pt = inFile->Get("L1MuonB_l1muon_Pt");
TH1F* h1L1MuonB_MIP_l1muon_Pt = inFile->Get("L1MuonB_MIP_l1muon_Pt");
TH1F* h1L1MuonB_HLT_l1muon_Pt = inFile->Get("L1MuonB_HLTRMatch_l1muon_Pt");
TH1F* h1L1MuonB_MIP_HLT_l1muon_Pt = inFile->Get("L1MuonB_MIP_HLTRMatch_l1muon_Pt");
TH1F* h1L1MuonB_NoHLT_l1muon_Pt = inFile->Get("L1MuonB_NoHLTRMatch_l1muon_Pt");
TH1F* h1L1MuonB_MIP_NoHLT_l1muon_Pt = inFile->Get("L1MuonB_MIP_NoHLTRMatch_l1muon_Pt");
TH1F* l1muon_ptHists[NUM_L1MUONPT_HISTS] = { h1L1MuonB_l1muon_Pt, h1L1MuonB_MIP_l1muon_Pt, 
					     h1L1MuonB_HLT_l1muon_Pt,h1L1MuonB_MIP_HLT_l1muon_Pt, 
					     h1L1MuonB_NoHLT_l1muon_Pt, h1L1MuonB_MIP_NoHLT_l1muon_Pt };

//L1 Muon Eta Histograms
const int NUM_ETA_HISTS = 6;
TH1F* h1L1MuonB_l1muon_Eta = inFile->Get("L1MuonB_l1muon_Eta");
TH1F* h1L1MuonB_MIP_l1muon_Eta = inFile->Get("L1MuonB_MIP_l1muon_Eta");
TH1F* h1L1MuonB_HLT_l1muon_Eta = inFile->Get("L1MuonB_HLTRMatch_l1muon_Eta");
TH1F* h1L1MuonB_MIP_HLT_l1muon_Eta = inFile->Get("L1MuonB_MIP_HLTRMatch_l1muon_Eta");
TH1F* h1L1MuonB_NoHLT_l1muon_Eta = inFile->Get("L1MuonB_NoHLTRMatch_l1muon_Eta");
TH1F* h1L1MuonB_MIP_NoHLT_l1muon_Eta = inFile->Get("L1MuonB_MIP_NoHLTRMatch_l1muon_Eta");
TH1F* etaHists[NUM_ETA_HISTS] = { h1L1MuonB_l1muon_Eta, h1L1MuonB_MIP_l1muon_Eta, 
				  h1L1MuonB_HLT_l1muon_Eta, h1L1MuonB_MIP_HLT_l1muon_Eta, 
				  h1L1MuonB_NoHLT_l1muon_Eta, h1L1MuonB_MIP_NoHLT_l1muon_Eta };

//HLT Pt Histograms
const int NUM_HLTPT_HISTS = 2;
TH1F* h1L1MuonB_HLT_hlt_Pt = inFile->Get("L1MuonB_HLTRMatch_hlt_Pt");
TH1F* h1L1MuonB_MIP_HLT_hlt_Pt = inFile->Get("L1MuonB_MIP_HLTRMatch_hlt_Pt");
TH1F* hlt_ptHists[NUM_HLTPT_HISTS] = { h1L1MuonB_HLT_hlt_Pt, h1L1MuonB_MIP_HLT_hlt_Pt};


//HLT Pt Histograms
const int NUM_HLTETA_HISTS = 2;
TH1F* h1L1MuonB_HLT_hlt_Eta = inFile->Get("L1MuonB_HLTRMatch_hlt_Eta");
TH1F* h1L1MuonB_MIP_HLT_hlt_Eta = inFile->Get("L1MuonB_MIP_HLTRMatch_hlt_Eta");
TH1F* hlt_EtaHists[NUM_HLTETA_HISTS] = { h1L1MuonB_HLT_hlt_Eta, h1L1MuonB_MIP_HLT_hlt_Eta};

//All Histograms
const int NUM_HISTS = 16;
TH1F* hists[NUM_HISTS] = { h1L1MuonB_l1muon_Pt, h1L1MuonB_MIP_l1muon_Pt, 
			   h1L1MuonB_HLT_l1muon_Pt, h1L1MuonB_MIP_HLT_l1muon_Pt, 
			   h1L1MuonB_NoHLT_l1muon_Pt, h1L1MuonB_MIP_NoHLT_l1muon_Pt, 
			   h1L1MuonB_l1muon_Eta, h1L1MuonB_MIP_l1muon_Eta, 
			   h1L1MuonB_HLT_l1muon_Eta, h1L1MuonB_MIP_HLT_l1muon_Eta,
			   h1L1MuonB_NoHLT_l1muon_Eta, h1L1MuonB_MIP_NoHLT_l1muon_Eta,
			   h1L1MuonB_HLT_hlt_Pt, h1L1MuonB_MIP_HLT_hlt_Pt,
			   h1L1MuonB_HLT_hlt_Eta, h1L1MuonB_MIP_HLT_hlt_Eta};
	    
			   
			   

/*
 * Using a series of helper functions, this 
 * produces the secondary Plots of interests.
 * Made using arrays of the primary histograms.
 */

void secondaryPlots(){
  //Sumw2
  //HistsSumw2();
  //Divide By Binwidth and Merge
  L1MuonPtHistsDivideByBinWidth();
  L1MuonPtHistsMergeDivideByBinWidth();  
  EtaHistsMerge();

  //Efficiency and Background
  //divideBinomialErrors(h1L1MuonB_MIP_NoHLT_l1muon_Eta, h1L1MuonB_NoHLT_l1muon_Eta, "Background_l1muon_Eta");
  
  divideBinomialErrors(h1L1MuonB_MIP_HLT_l1muon_Pt, h1L1MuonB_HLT_l1muon_Pt,"Efficiency_l1muon_Pt");
  divideBinomialErrors(h1L1MuonB_MIP_HLT_l1muon_Eta, h1L1MuonB_HLT_l1muon_Eta,"Efficiency_l1muon_Eta");

  L1MuonPtMerge_divideBinomialErrors(h1L1MuonB_MIP_HLT_l1muon_Pt, h1L1MuonB_HLT_l1muon_Pt,"Efficiency_Merge_l1muon_Pt");
  L1MuonPtMerge_divideBinomialErrors(h1L1MuonB_MIP_NoHLT_l1muon_Pt, h1L1MuonB_NoHLT_l1muon_Pt,"Background_Merge_l1muon_Pt");
  EtaMerge_divideBinomialErrors(h1L1MuonB_MIP_HLT_l1muon_Eta, h1L1MuonB_HLT_l1muon_Eta,"Efficiency_Merge_l1muon_Eta");
  EtaMerge_divideBinomialErrors(h1L1MuonB_MIP_NoHLT_l1muon_Eta, h1L1MuonB_NoHLT_l1muon_Eta,"Background_Merge_El1muon_ta");

  //Specical Case HLT
  L1MuonPtMerge_divideBinomialErrors(h1L1MuonB_MIP_HLT_hlt_Pt, h1L1MuonB_HLT_hlt_Pt, "Efficiency_hlt_Pt");
  EtaMerge_divideBinomialErrors(h1L1MuonB_MIP_HLT_hlt_Eta, h1L1MuonB_HLT_hlt_Eta, "Efficiency_hlt_Eta");



  //assymetricDivide(h1L1MuonB_MIP_HLT_l1muon_Pt, h1L1MuonB_HLT_l1muon_Pt,"Efficiency_Pt"); 
  //efficiency(h1L1BMipHlt_L1MuonPt, h1L1BHlt_L1MuonPt, "Efficiency_L1MuonPt");
  


  //Signal 
  std::stringstream Signal;
  Signal  <<"Signal" << endl;
  //Signal << ptHat <<"_Signal" << endl;
  //divide("L1MuonBarrelwithMipMatchHLT_Pt", "L1MuonBarrelMatchHLT_Pt", Signal.str().c_str());
  //divide(h1L1MuonBMipHltMerge, h1L1MuonBHltMerge, "Efficiency_Merge");
  //divide(h1L1BMipHltEtaRebin, h1L1BHltEtaRebin, "Eta_Efficiency_Rebin");

  
  //Background
  
  //TH1F* back_num = subtract("L1MuonBwithMip_Pt", "L1MuonBarrelwithMipMatchHLT_Pt", "L1MuonBarrelwithMipNOHLTMatch_Pt");
  //TH1F* back_den = subtract("L1Muon_Barrel_Pt", "L1MuonBarrelMatchHLT_Pt", "L1MuonBarrelwithoutHLTMatch_Pt");

  // TH1F* back_num_Merge = subtract("L1MuonBwithMip_Pt_Merge", "L1MuonBarrelwithMipMatchHLT_Pt_Merge", "L1MuonBarrelwithMipNOHLTMatch_Pt_Merge");
  //TH1F* back_den_Merge = subtract("L1Muon_Barrel_Pt_Merge", "L1MuonBarrelMatchHLT_Pt_Merge", "L1MuonBarrelwithoutHLTMatch_Pt_Merge");

  //std::stringstream Background;
  //Background <<"Background" << endl;
  
  //divide(back_num, back_den, "Background");
  //divide(back_num_Merge, back_den_Merge, "Background_Merge");

  //newFile->Write();
  inFile->Close();
}



/*
 * Loop over array of Hists, call Sumw2 on each
 */
/*
void HistsSumw2(){
  cout << "Call Sumw2 on:" << endl;
  for(int i=0; i< NUM_HISTS; i++){
    cout << "   " << hists[i]->GetName() <<endl;
    hists[i]->Sumw2();
  }
}
*/

/*
 * Divides Pt Histograms in Array, by their bin width,
 * labels axises, and writes to file.
 */
void L1MuonPtHistsDivideByBinWidth(){
  cout << "Divide By Bin Width:" << endl;
  for(int i=0; i<NUM_L1MUONPT_HISTS; i++){
  //for(int i =0; i< sizeof(l1muon_ptHists[])/sizeof(TH1F*); ++i){
    //for(auto ptHist = l1muon_ptHists->begin(); ptHist != l1muon_ptHists->end(); ++ptHist){
    cout << "   " << l1muon_ptHists[i]->GetName() <<endl;
    TH1F* h1PtNorm = divideByBinWidth(l1muon_ptHists[i]);
    SetAxises(h1PtNorm, "Pt (GeV)", "dN/dx");
    h1PtNorm->Write();
    }
}

/*
 * Define Merge Bin Array
 * Merge Bins in Histograms (Mainly to deal with low Statistics
 * at high and low Pt).
 * Divide by Bin Width
 */

void L1MuonPtHistsMergeDivideByBinWidth(){
  const int numBins_Pt=20;
  Double_t mergeVariableBinArray[numBins_Pt+1] = {0, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 10, 12, 14, 16, 18, 20, 25, 30, 40, 60, 100, 140};

  cout << "Merge and Divide By Bin Width:" << endl;
  for(int i=0; i<NUM_L1MUONPT_HISTS; i++){
    cout << "   " << l1muon_ptHists[i]->GetName() <<endl;
    TH1F* histMerged = merge(l1muon_ptHists[i], numBins_Pt, mergeVariableBinArray);
    
    TH1F* histMergedNorm = divideByBinWidth(histMerged);
    SetAxises(histMergedNorm, "Pt (GeV)", "dN/dx");
    histMergedNorm->Write();
  }
}

L1MuonPtMerge_divideBinomialErrors(TH1F* h1Num, TH1F* h1Den, char *quotientName){
  const int numBins_Pt=20;
  Double_t mergeVariableBinArray[numBins_Pt+1] = {0, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 10, 12, 14, 16, 18, 20, 25, 30, 40, 60, 100, 140};
  TH1F* h1NumMerged = merge(h1Num, numBins_Pt, mergeVariableBinArray);
  TH1F* h1DenMerged = merge(h1Den, numBins_Pt, mergeVariableBinArray);
  divideBinomialErrors(h1NumMerged, h1DenMerged,quotientName);
}
  
    
void EtaHistsMerge(){
  const int numBins_Eta=26;
  Double_t etaBinLowEdges[numBins_Eta+1];
  for(int i_bin=0; i_bin <= numBins_Eta; i_bin++) etaBinLowEdges[i_bin]= -1.3+0.1*i_bin;
  cout << "Merge:" << endl;
  for(int i=0; i<NUM_ETA_HISTS; i++){
    cout << "   " << etaHists[i]->GetName() <<endl;
    TH1F* histMerged =merge(etaHists[i], numBins_Eta, etaBinLowEdges);
    SetAxises(histMerged, "Eta", "Weighted Counts");
    histMerged->Write();
  }
}

EtaMerge_divideBinomialErrors(TH1F* h1Num, TH1F* h1Den, char *quotientName){
  const int numBins_Eta=26;
  Double_t etaBinLowEdges[numBins_Eta+1];
  for(int i_bin=0; i_bin <= numBins_Eta; i_bin++) etaBinLowEdges[i_bin]= -1.3+0.1*i_bin;
  TH1F* h1NumMerged = merge(h1Num, numBins_Eta, etaBinLowEdges);
  TH1F* h1DenMerged = merge(h1Den, numBins_Eta, etaBinLowEdges);
  divideBinomialErrors(h1NumMerged, h1DenMerged,quotientName);
  
}

/*
 * My wrapper for the ROOT Rebin function
 */
TH1F* merge(TH1F* hist, int numBins, Double_t* xBinLowEdges){
  
  std::stringstream mergeName;
  mergeName << hist->GetName() <<"_Merged";
  //std::cout << histLoc.str().c_str() << endl;
  //hist->Sumw2();  
  TH1F* histRebin = hist->Rebin(numBins, mergeName.str().c_str(), xBinLowEdges); //Creates a new histrogram called mergeName
  //TH1F* histRebin = (TH1F*)inFile->Get(rebinName.str().c_str());  //Take Histogram from Memory
  return histRebin;
}

/*
 * Creates a new Histogram.  By taking the entries of the orginal histogram,
 * and dividing them by their bin width.
 */

TH1F* divideByBinWidth(TH1F* h1){  
  std::stringstream  normName;
  normName << h1->GetName() <<"_Norm";
  TH1F* h1Norm = (TH1F*) h1->Clone(normName.str().c_str());
  h1Norm->Reset();
  h1Norm->SetTitle(normName.str().c_str());
  h1Norm->GetYaxis()->SetTitle("#frac{dN}{dx}");
  // h1->Sumw2();
  int num_bins = (int) h1Norm->GetNbinsX();
  for(int i=0; i<=num_bins; i++){
    float binContent = h1->GetBinContent(i);
    float binWidth = h1->GetBinWidth(i);
    float binCenter = h1->GetBinCenter(i);
    h1Norm->Fill(binCenter, binContent/binWidth);
    h1Norm->SetBinError(i, h1->GetBinError(i)/binWidth);  //Manually set
    /*
     * Check Errors
    cout << i << " " << h1->GetBinError(i) << endl;
    cout << h1Norm->GetBinError(i) << " " << sqrt(h1Norm.GetBinContent(i))<< " " 
    << h1->GetBinError(i)/binWidth << endl;
    cout << i << ", " << h1Norm_test->GetBinCenter(i) << endl;
    */    
  }
  return h1Norm;
}


/*
 * Divides the first histogram by the second histogram.
 */

void assymetricDivide(TH1F* h1Num, TH1F* h1Den, char *quotientName)
{
 //TGraphAsymmErrors * binomialErrors(h1Num, h1Den);
  binomialErrors = new TGraphAsymmErrors(h1Num, h1Den, "n");
  // binomialErrors->Print();
  //c1 = new TCanvas("c1","A Simple Graph with assymetric error bars",200,10,700,500);  
  binomialErrors->SetTitle(quotientName);
  //binomialErrors->GetXaxis()->SetTitle("Pt (GeV)");
  
  //binomialErrors->GetXaxis()->SetTitle("Pt");
  //binomialErrors->GetYaxis()->SetTitle("Ratio");
  binomialErrors->SetMarkerColor(4);
  binomialErrors->GetYaxis()->SetRange(0.0,1.0);
  //binomialErrors->SetMarkerStyle(21);
  binomialErrors->Draw("AP");
  binomialErrors->Write();
  //binomialErrors->Write("AP");
  examineBinErrors(h1Num);
  examineBinErrors(h1Den);
  
  //TH1F* h1Div = (TH1F*) h1Num->Clone(quotientName);
  //h1Div->Reset();
  //h1Div->SetTitle(quotientName);
  //h1Div->GetYaxis()->SetTitle(quotientName);
  
  //h1Div->Sumw2()
  //h1Num->Sumw2();
  //h1Den->Sumw2();
  //h1Div->Divide(h1Num, h1Den);
    
  //cout << "Binomial Errors" << endl; 
  //cout << "Quotient Content and Errors" << endl;
  //examineBinErrors(h1Div);
 
  //h1Div->Write();
  
}

//"Divide two histograms with setting at binomial errors.

void divideBinomialErrors(TH1F* h1Num, TH1F* h1Den, char *quotientName){
  TH1F* h1Quotient = (TH1F*) h1Num->Clone(quotientName);
  h1Quotient->Reset();
  h1Quotient->Divide(h1Num, h1Den, 1, 1, "B"); //b specifies binomial errors, handles weights properly.
  h1Quotient->SetMaximum(1.04);
  h1Quotient->SetMinimum(0);
  gStyle->SetOptStat(0);
  h1Quotient->Write();
  h1Quotient->Draw();
}



/*
void efficiency(TH1F* h1Num, TH1F* h1Den, char *quotientName)
{
  examineBinErrors(h1Num);
  examineBinErrors(h1Den);
    
  //cout << "Binomial Errors" << endl;
  ///TGraphAsymmErrors * binomialErrors(h1Num, h1Den);
  binomialErrors = new TGraphAsymmErrors(h1Num, h1Den, "n");
  
  //binomialErrors->Print();
  
  c1 = new TCanvas("c1","A Simple Graph with assymetric error bars",200,10,700,500);  
  binomialErrors->Draw("AP");
  
  
  //TEfficiency* efficiency = new TEfficiency(*h1Num, *h1Den);
  //efficiency->Draw("AP");

  //binomialErrors->SetTitle(quotientName);
  //binomialErrors->GetXaxis()->SetTitle("Pt (GeV)");
  
  //binomialErrors->GetXaxis()->SetTitle("Pt");
  //binomialErrors->GetYaxis()->SetTitle("Ratio");
  //binomialErrors->SetMarkerColor(4);
  //binomialErrors->SetMarkerStyle(21);
  //binomialErrors->Draw("AP");
  //binomialErrors->Write("AP");
}
*/

/*
 * Wrapper For divide, that takes the names as strings
 */

/*
void divide(char * histNumName, char *histDenName, char *quotientName)
{
  
  std::stringstream histNumLoc, histDenLoc;
  histNumLoc <<histDir << histNumName;
  histDenLoc <<histDir << histDenName;

  TH1F* h1Num = inFile->Get(histNumLoc.str().c_str());
  TH1F* h1Den = inFile->Get(histDenLoc.str().c_str());

  divide(h1Num, h1Den, quotientName);
}
*/	       


TH1F* subtract(char * h1FirstName, char *h1SecondName, char* difName)
{

  
  std::stringstream h1FirstLoc, h1SecondLoc;
  h1FirstLoc << h1FirstName;
  h1SecondLoc << h1SecondName;

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
