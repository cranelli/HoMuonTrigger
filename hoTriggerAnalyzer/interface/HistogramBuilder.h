#ifndef __HOMUON__HISTOGRAMBUILDER_H__
#define __HOMUON__HISTOGRAMBUILDER_H__

/*
 * Class Test Class
 * Author Chris Anelli
 * 6.13.2013
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "TH1F.h"
#include "TH2.h"
#include <map>



class HistogramBuilder {
  
  //int mynumber = 13;
  
 public:

  //HistogramBuilder(); //Constructor
 
  edm::Service<TFileService> _fileService;
  
  /*
   * Functions for HistogramBuilder                                       
   */
  
  void fillL1MuonPtHistograms(float pt, std::string key);
  void fillPtHistograms(float pt, std::string key);
  void fillEnergyHistograms(float energy, std::string key);                   
  void fillEtaPhiHistograms(float eta, float phi, std::string key);
  void fillDeltaEtaDeltaPhiHistograms(float eta1, float eta2, 
				      float phi1, float phi2,  std::string key);
  void fillCountHistogram(std::string key);                                   
  void fillTrigHistograms(bool trigDecision,std::string key);

  //int getMyNumber(); //{return mynumber;}

 private:

  std::map<std::string,TH1F*> _h1L1MuonPt;
  std::map<std::string,TH1F*> _h1Pt;
  std::map<std::string,TH1F*> _h1InvPt;
  std::map<std::string,TH1F*> _h1Energy;                                 
  std::map<std::string,TH1F*> _h1Eta;
  std::map<std::string,TH1F*> _h1Phi;
  std::map<std::string,TH2F*> _h2EtaPhiMap;
  std::map<std::string,TH1F*> _h1DeltaEta;                                   
  std::map<std::string,TH1F*> _h1DeltaPhi;
  std::map<std::string,TH2F*> _h2DeltaEtaDeltaPhi;
  std::map<std::string,TH1F*> _h1Trig;
  std::map<std::string,TH1F*> _h1Counter; 
  
  /*
   * Helper Functions
   */

  void SetAxises(TH1F * h1, std::string xtitle, std::string ytitle);
  void SetAxises(TH2F * h2, std::string xtitle, std::string ytitle);
  //float WrapCheck(float phi1, float phi2);
  
};

#endif
