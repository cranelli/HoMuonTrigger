#include "HoMuonTrigger/hoTriggerAnalyzer/interface/HistogramBuilder.h"

/*
 * The HistogramBuilder Class contains  
 * functions to build generic histograms of                           
 * different types.  The type of object (or a specific                      
 * selection cut) can be sepecified by the key.
 * All histograms are saved through the TFileService 
 *                           
 * Created by Christopher Anelli
 * On 6.15.2014
*/

#include "TH1F.h"
#include "TH2F.h"
#include "HoMuonTrigger/hoTriggerAnalyzer/interface/CommonFunctions.h"

// Constructor
/*
HistogramBuilder::HistogramBuilder(){
};
*/

/*
 *Counting Histograms
 *Fills the 1 bin.
 */
void HistogramBuilder::fillCountHistogram(std::string key){                     
  if(!_h1Counter.count(key)){                                                   
    _h1Counter[key] = _fileService->make<TH1F>(Form("%s_Count",key.c_str()),    
                                               Form("%s Count",key.c_str()),    
                                               2, 0, 2);
    _h1Counter[key]->GetYaxis()->SetTitle("Counts");
  }                                                                             
  _h1Counter[key]->Fill(1);                                                     
}                                                                               

/*      
 *Trigger Histograms
 */
void HistogramBuilder::fillTrigHistograms(bool trigDecision,std::string key){ 
  if(!_h1Trig.count(key)){                                                   
    _h1Trig[key] = _fileService->make<TH1F>(Form("%s_Trig",key.c_str()),     
                                            Form("%s Trigger",key.c_str()),  
                                            2, 0, 2);
    SetAxises(_h1Trig[key],"Trigger Decision", "Counts");
  }                                                      
                       
  _h1Trig[key]->Fill(trigDecision);                                        
     
}

/*      
 *Weight Histograms
 */
void HistogramBuilder::fillWeightHistograms(float weight,std::string key){ 
  if(!_h1Weight.count(key)){                                                   
    _h1Weight[key] = _fileService->make<TH1F>(Form("%s_Weight",key.c_str()),     
                                            Form("%s Weight",key.c_str()),  
                                            2, 0, 2);
    SetAxises(_h1Weight[key],"Weight", "Counts");
  }                                                      
                       
  _h1Weight[key]->Fill(weight);                                        
     
}
      

/*                                          
 *Energy Histograms
 */

void HistogramBuilder::fillEnergyHistograms(float energy, std::string key){  
  if(!_h1Energy.count(key)){                                             
    _h1Energy[key] = _fileService->make<TH1F>(Form("%s_Energy",key.c_str()),
                                            Form("%s Energy",key.c_str()), 
                                            2100, -5.0, 100.0);
    SetAxises(_h1Energy[key],"Energy (GeV)", "Counts");
  }                                                                          
  _h1Energy[key]->Fill(energy);                                            
}  

/*                                                                             
 *Eta Phi Histograms                                                         
 */

void HistogramBuilder::fillEtaPhiHistograms(float eta, float phi, std::string key){
  if(!_h1Eta.count(key)){
    _h1Eta[key] = _fileService->make<TH1F>(Form("%s_Eta",key.c_str()),
                                           Form("%s Eta",key.c_str()),
                                           500, -1.5, 1.5);  //HO has 72 iphis and 30 ietas
    SetAxises(_h1Eta[key], "#eta", "Counts");
    
  }
  _h1Eta[key]->Fill(eta);

  if(!_h1Phi.count(key)){
    _h1Phi[key] = _fileService->make<TH1F>(Form("%s_Phi",key.c_str()),
                                           Form("%s Phi",key.c_str()),
                                           500, -3.14, 3.14);  //HO has 72 iphis and 30 ietas                         
    SetAxises(_h1Phi[key], "#phi", "Counts");
  }
  _h1Phi[key]->Fill(phi);

  if(!_h2EtaPhiMap.count(key)){
    _h2EtaPhiMap[key] = _fileService->make<TH2F>(Form("%s_EtaPhi",key.c_str()),
                                                 Form("%s_EtaPhi",key.c_str()),
						 500, -1.5, 1.5, 500, -3.14, 3.14);
    SetAxises(_h2EtaPhiMap[key], "#eta","#phi");
  }
  _h2EtaPhiMap[key]->Fill(eta, phi);

};

/*
 *Delta Eta Delta Phi Histograms
 */

void HistogramBuilder::fillDeltaEtaDeltaPhiHistograms(float eta1, float eta2, 
						      float phi1, float phi2, 
						      std::string key){
  float deltaEta, deltaPhi;
  CommonFunctions commonFunctions;
  deltaEta = eta1 - eta2;
  deltaPhi = commonFunctions.WrapCheck(phi1, phi2);

  //Delta Eta Histograms Fill
  if(!_h1DeltaEta.count(key)){
        _h1DeltaEta[key] =  _fileService->make<TH1F>(Form("%s_DeltaEta",
							  key.c_str()), 
						     Form("#Delta #eta %s",
							  key.c_str()),   
						     2000, -2.6, 2.6);
	SetAxises(_h1DeltaEta[key],"#Delta #eta", "Counts");
  }
  _h1DeltaEta[key]->Fill(deltaEta);
    
  //Delta Phi Histograms Fill
  if(!_h1DeltaPhi.count(key)){
    _h1DeltaPhi[key] = _fileService->make<TH1F>(Form("%s_DeltaPhi",
						     key.c_str()), 
						Form("%s #Delta #Phi",
						     key.c_str()),    
                                                2000, -3.14, 3.14);
    SetAxises(_h1DeltaPhi[key],"#Delta #phi", "Counts");
  }
  _h1DeltaPhi[key]->Fill(deltaPhi);
  
  //DeltaEta Delta Phi Histograms Fill
  if(!_h2DeltaEtaDeltaPhi.count(key)){
    _h2DeltaEtaDeltaPhi[key] = _fileService->make<TH2F>(Form("%s_DeltaEtaDeltaPhi",key.c_str()),Form("%s #Delta#eta #Delta#Phi",key.c_str()), 2000, -2.6, 2.6, 2000, -3.14, 3.14);
    SetAxises(_h2DeltaEtaDeltaPhi[key],"#Delta #eta", "#Delta #phi");
  }
  _h2DeltaEtaDeltaPhi[key]->Fill(deltaEta, deltaPhi);
} 


/*
 *Pt Histograms
 *has variable binning
 */
void HistogramBuilder::fillPtHistograms(float pt, std::string key){
  if(!_h1Pt.count(key)){
        
    _h1Pt[key] = _fileService->make<TH1F>(Form("%s_Pt",key.c_str()),
					  Form("%s Pt",key.c_str()),
					  800,0,400);
    SetAxises(_h1Pt[key],"Pt (GeV)", "Counts");
  } 
  _h1Pt[key]->Fill(pt);

if(!_h1InvPt.count(key)){
        
    _h1InvPt[key] = _fileService->make<TH1F>(Form("%s_InvPt",key.c_str()),
					  Form("%s Inv Pt",key.c_str()),
					  800,-1,1);
    SetAxises(_h1InvPt[key], "#frac{1}{P_t}", "Counts");
  } 
 _h1InvPt[key]->Fill(1.0/pt);
}


/*
 *L1Muon Pt Histograms
 *has variable binning
 */
void HistogramBuilder::fillL1MuonPtHistograms(float pt, std::string key){
  if(!_h1L1MuonPt.count(key)){
    
    float variableBinArray[] = {0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,10,12,14,16,18,20,25,30,35,40,45,50,60,70,80,100,120,140,180};
    
    _h1L1MuonPt[key] = _fileService->make<TH1F>(Form("%s_Pt",key.c_str()),
                                                Form("%s Pt",key.c_str()),
                                                33,
                                                variableBinArray);
    SetAxises(_h1L1MuonPt[key], "Pt", "Counts");
  } 
  _h1L1MuonPt[key]->Fill(pt);
}         

/*
int HistogramBuilder::getMyNumber(){
  return mynumber;
};
*/


/*
 *
 * Helper Functions
 *
 */


void HistogramBuilder::SetAxises(TH1F * h1, std::string xTitle, std::string yTitle){
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
void HistogramBuilder::SetAxises(TH2F * h2, std::string xTitle, std::string yTitle){
  // X Axis
  h2->GetXaxis()->SetTitle(xTitle.c_str());
  h2->GetXaxis()->SetTitleFont(42);
  h2->GetXaxis()->SetTitleSize(0.05);
  h2->GetXaxis()->SetTitleOffset(0.9);

  // Y Axis
  h2->GetYaxis()->SetTitle(yTitle.c_str());
  h2->GetYaxis()->SetTitleFont(42);
  h2->GetYaxis()->SetTitleSize(0.05);
  h2->GetYaxis()->SetTitleOffset(1.0);

  h2->SetLineWidth(2);
  h2->SetLineColor(kBlue);
}
