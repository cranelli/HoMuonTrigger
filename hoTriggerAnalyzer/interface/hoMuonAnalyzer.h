#ifndef __HOMUON_HOMUONANALYZER_H__
#define __HOMUON_HOMUONANALYZER_H__
//
// Class:hoMuonAnalyzer
// 
/*
 Description: Currently the header file for the entire HOMuon Trigger Analysis

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Christopher Anelli  
//         Created:  Fri, 16 May 2014 04:20:05 GMT


// system include files
#include <memory>

// user include files


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"
//#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/HcalRecHit/interface/HORecHit.h"

#include "Analysis/hoTriggerAnalyzer/interface/HistogramBuilder.h"

#include <vector>
#include <iostream>
#include <map>
#include <list>

using namespace::std;

static const float threshold = 0.2;
static const float deltaR_Max = 0.3;

//
// class declaration
//

class hoMuonAnalyzer : public edm::EDAnalyzer {
public:
  explicit hoMuonAnalyzer(const edm::ParameterSet&);
  ~hoMuonAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& 
				   descriptions);
  

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, 
		       const edm::EventSetup&) override;
  virtual void endJob() override;
  virtual void beginRun (const edm::Run& iRun, 
			 const edm::EventSetup& evSetup);
  virtual void endRun(const edm::Run& iRun, 
		      const edm::EventSetup& evSetup);
  
  edm::Service<TFileService> _fileService;

  edm::InputTag _genInput;
  edm::InputTag _l1MuonInput;
  //edm::InputTag _stdMuInput;
  edm::InputTag _horecoInput;
  //edm::InputTag _l1GtTmLInputTag;
  edm::InputTag _hltSumAODInput;
  HistogramBuilder histogramBuilder;

  
  // I would prefer to run without an InputTag, the L1GtUtility should be     
  // able to find it automatically from the Providence information.        
  //edm::InputTag m_l1GtTmLInputTag;
  
  L1GtUtils m_l1GtUtils;
  //HLTConfigProvider hltConfig_;  //Copied from HLT TrigReport_h
  string m_nameAlgTechTrig;
  bool trigDecision;

  //toFigureOutL1VariableBinning
  std::list<float> listL1MuonPt;

  /*
   * Maps of selected hlt triggers to get the trigger decisions,
   * and hlt filters to get the trigger objects.
   */
  void defineTriggersOfInterest();
  map<string, string> hltNamesOfInterest;
  std::map<std::string, edm::InputTag> hltFiltersOfInterest;
  
  //For Filtering
  // bool hoBelowThreshold(HORecHit horeco);


  //Helper Functions for hoMuonAnalyzer Initilization
  //void initializeHistograms();


  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
};

#endif
