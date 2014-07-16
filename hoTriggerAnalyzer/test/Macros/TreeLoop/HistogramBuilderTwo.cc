#include "HistogramBuilderTwo.h"

/*
 * The HistogramBuilderTwo Class contains  
 * functions to build generic histograms of                           
 * different types.  The type of object (or a specific                      
 * selection cut) can be sepecified by the key.
 * Same as Histogram Builder, except files are no longer
 * held through TFileService.
 *                           
 * Created by Christopher Anelli
 * On 6.15.2014
*/

#include "TH1F.h"
#include "TH2F.h"
#include "math.h"
#include <iostream>
//#include "HoMuonTrigger/hoTriggerAnalyzer/interface/CommonFunctions.h"




// Constructor
/*
HistogramBuilderTwo::HistogramBuilderTwo(){
};
*/

/*
 *Counting Histograms
 *Fills the 1 bin.
 */
void HistogramBuilderTwo::fillCountHistogram(std::string key, double weight){
  if(!_h1Counter.count(key)){                                                   
    _h1Counter[key] = new TH1F(Form("%s_Count",key.c_str()),    
                                               Form("%s Count",key.c_str()),    
                                               2, 0, 2);
    _h1Counter[key]->GetYaxis()->SetTitle("Counts");
  }                                                                             
  _h1Counter[key]->Fill(1, weight);
}                                                                               

/*      
 *Trigger Histograms
 */
void HistogramBuilderTwo::fillTrigHistograms(bool trigDecision,std::string key, double weight){ 
  if(!_h1Trig.count(key)){                                                   
    _h1Trig[key] = new TH1F(Form("%s_Trig",key.c_str()),     
                                            Form("%s Trigger",key.c_str()),  
                                            2, 0, 2);
    SetAxises(_h1Trig[key],"Trigger Decision", "Counts");
  }                                                      
                       
  _h1Trig[key]->Fill(trigDecision, weight);                                        
     
}

/*      
 *Weight Histograms
 */
void HistogramBuilderTwo::fillWeightHistograms(float weight_val,std::string key, double weight){ 
  if(!_h1Weight.count(key)){                                                   
    _h1Weight[key] = new TH1F(Form("%s_Weight",key.c_str()),     
					      Form("%s Weight",key.c_str()),  
					      1000, 0, 1);
    SetAxises(_h1Weight[key],"Weight", "Counts");
  }                                                      
                       
  _h1Weight[key]->Fill(weight_val, weight);                                        
     
}
      

/*                                          
 *Energy Histograms
 */

void HistogramBuilderTwo::fillEnergyHistograms(float energy, std::string key, double weight){  
  if(!_h1Energy.count(key)){                                             
    _h1Energy[key] = new TH1F(Form("%s_Energy",key.c_str()),
                                            Form("%s Energy",key.c_str()), 
                                            2100, -5.0, 100.0);
    SetAxises(_h1Energy[key],"Energy (GeV)", "Counts");
  }                                                                          
  _h1Energy[key]->Fill(energy, weight);                                            
}  

/*                                                                             
 *Eta Phi Histograms                                                         
 */

void HistogramBuilderTwo::fillEtaPhiHistograms(float eta, float phi, std::string key, double weight){
  if(!_h1Eta.count(key)){
    _h1Eta[key] = new TH1F(Form("%s_Eta",key.c_str()),
                                           Form("%s Eta",key.c_str()),
                                           500, -1.5, 1.5);  //HO has 72 iphis and 30 ietas
    SetAxises(_h1Eta[key], "#eta", "Counts");
    
  }
  _h1Eta[key]->Fill(eta, weight);

  if(!_h1Phi.count(key)){
    _h1Phi[key] = new TH1F(Form("%s_Phi",key.c_str()),
                                           Form("%s Phi",key.c_str()),
                                           500, -3.2, 3.2);  //HO has 72 iphis and 30 ietas                         
    SetAxises(_h1Phi[key], "#phi", "Counts");
  }
  _h1Phi[key]->Fill(phi, weight);

  if(!_h2EtaPhiMap.count(key)){
    _h2EtaPhiMap[key] = new TH2F(Form("%s_EtaPhi",key.c_str()),
                                                 Form("%s_EtaPhi",key.c_str()),
						 500, -1.5, 1.5, 500, -3.14, 3.14);
    SetAxises(_h2EtaPhiMap[key], "#eta","#phi");
  }
  _h2EtaPhiMap[key]->Fill(eta, phi, weight);

};

/*
 *Delta Eta Delta Phi Histograms
 */

void HistogramBuilderTwo::fillDeltaEtaDeltaPhiHistograms(float eta1, float eta2, 
						      float phi1, float phi2, 
						      std::string key, double weight){
  float deltaEta, deltaPhi;
  //CommonFunctions commonFunctions;
  deltaEta = eta1 - eta2;
  //deltaPhi = commonFunctions.WrapCheck(Phi1, phi2);
  deltaPhi = WrapCheck(phi1, phi2);

  //Delta Eta Histograms Fill
  if(!_h1DeltaEta.count(key)){
        _h1DeltaEta[key] =  new TH1F(Form("%s_DeltaEta",
							  key.c_str()), 
						     Form("#Delta #eta %s",
							  key.c_str()),   
						     2000, -2.6, 2.6);
	SetAxises(_h1DeltaEta[key],"#Delta #eta", "Counts");
  }
  _h1DeltaEta[key]->Fill(deltaEta, weight);
    
  //Delta Phi Histograms Fill
  if(!_h1DeltaPhi.count(key)){
    _h1DeltaPhi[key] = new TH1F(Form("%s_DeltaPhi",
						     key.c_str()), 
						Form("%s #Delta #Phi",
						     key.c_str()),    
                                                2000, -3.2, 3.2);
    SetAxises(_h1DeltaPhi[key],"#Delta #phi", "Counts");
  }
  _h1DeltaPhi[key]->Fill(deltaPhi, weight);
  
  //DeltaEta Delta Phi Histograms Fill
  if(!_h2DeltaEtaDeltaPhi.count(key)){
    _h2DeltaEtaDeltaPhi[key] = new TH2F(Form("%s_DeltaEtaDeltaPhi",key.c_str()),Form("%s #Delta#eta #Delta#Phi",key.c_str()), 2000, -2.6, 2.6, 2000, -3.14, 3.14);
    SetAxises(_h2DeltaEtaDeltaPhi[key],"#Delta #eta", "#Delta #phi");
  }
  _h2DeltaEtaDeltaPhi[key]->Fill(deltaEta, deltaPhi, weight);
} 


/*
 *Pt Histograms
 *has variable binning
 */
void HistogramBuilderTwo::fillPtHistograms(float pt, std::string key, double weight){
  if(!_h1Pt.count(key)){
        
    _h1Pt[key] = new TH1F(Form("%s_Pt",key.c_str()),
					  Form("%s Pt",key.c_str()),
					  6000,0,3000);
    SetAxises(_h1Pt[key],"Pt (GeV)", "Counts");
  } 
  _h1Pt[key]->Fill(pt, weight);

if(!_h1InvPt.count(key)){
        
    _h1InvPt[key] = new TH1F(Form("%s_InvPt",key.c_str()),
					  Form("%s Inv Pt",key.c_str()),
					  800,-1,1);
    SetAxises(_h1InvPt[key], "#frac{1}{P_t}", "Counts");
  } 
 _h1InvPt[key]->Fill(1.0/pt, weight);
}


/*
 *L1Muon Pt Histograms
 *has variable binning
 */
void HistogramBuilderTwo::fillL1MuonPtHistograms(float pt, std::string key, double weight){
  if(!_h1L1MuonPt.count(key)){
    
    float variableBinArray[] = {0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,10,12,14,16,18,20,25,30,35,40,45,50,60,70,80,100,120,140,180};
    
    _h1L1MuonPt[key] = new TH1F(Form("%s_Pt",key.c_str()),
                                                Form("%s Pt",key.c_str()),
                                                32,
                                                variableBinArray);
    SetAxises(_h1L1MuonPt[key], "Pt", "Counts");
  } 
  _h1L1MuonPt[key]->Fill(pt, weight);
}         

/*
int HistogramBuilderTwo::getMyNumber(){
  return mynumber;
};
*/


/*
 *
 * Helper Functions
 *
 */


void HistogramBuilderTwo::SetAxises(TH1F * h1, std::string xTitle, std::string yTitle){
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
void HistogramBuilderTwo::SetAxises(TH2F * h2, std::string xTitle, std::string yTitle){
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

/*
 * Handles wrapping between two angles.
 * So they are never more than 2 PI radians separated.
 * Returns values between -PI and PI.
 */



float HistogramBuilderTwo::WrapCheck(float phi1, float phi2){
  static const float PI = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348;
  
  float delta_phi = phi1 - phi2;
  if (delta_phi > PI) delta_phi-= 2*PI;
  if(delta_phi < -1*PI) delta_phi += 2*PI;
    
  return delta_phi;
  /*
  float delta_case;
   delta_case = delta_phi;
  if(delta_phi < -PI){
    delta_case = (2*PI + delta_phi);
  }
  if(delta_phi > PI){
    delta_case = (delta_phi - 2*PI);
  }
    
  //float delta_mod = fmod(pos_angle + PI, 2*PI) - PI;
 
  // Use the float modulus operator
  */
};
