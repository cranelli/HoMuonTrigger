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
const char* hist3Name = "L1MuonBarrelwithMipMatchHLT_Pt_Norm";

void Overlay(){
  TFile * rootFile = TFile::Open(rootName, "READ");
  TH1F * hist1 = rootFile->Get(hist1Name);
  hist1->Draw();
  TH1F * hist2 = rootFile->Get(hist2Name);
  hist2->SetLineColor(kRed);
  hist2->Draw("same");
  TH1F * hist3 = rootFile->Get(hist3Name);
  hist3->SetLineColor(kBlack); //Black
  hist3->Draw("same");

  //rootFile->Close();
}
