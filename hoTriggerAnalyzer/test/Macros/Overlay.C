#include "TFile.h"
#include "TH1F"
#include <iostream>

/*
 * Quick Macro to Overlay Two Histograms
 * First One is Plotted in Blue
 * Second One is Plotted in Red
 */


const char* rootName = "TreeLoop/L1MuonHistogram_Plus.root";
const char* hist1Name = "L1Muon_Barrel_Pt_Norm";
const char* hist2Name = "L1MuonBwithMip_Pt_Norm";
const char* hist3Name = "L1MuonBarrelMatchHLT_Pt_Norm";
const char* hist4Name = "L1MuonBarrelwithMipMatchHLT_Pt_Norm";
const char* hist5Name = "L1MuonBarrel_NoHLT_Pt_Norm";
const char* hist6Name = "L1MuonBarrelwithMip_NoHLT_Pt_Norm";

void Overlay(){
  TFile * rootFile = TFile::Open(rootName, "READ");
  TH1F * hist5 = rootFile->Get(hist5Name);
  hist5->SetLineColor(kBlack);
  hist5->Draw();
  TH1F * hist6 = rootFile->Get(hist6Name);
  hist6->SetLineColor(kRed);
  hist6->Draw("same");
  //TH1F * hist3 = rootFile->Get(hist3Name);
  //hist3->SetLineColor(kBlack); //Black
  //hist3->Draw("same");

  //rootFile->Close();
}
