// -*- C++ -*-
//
// Package:    caloInspector
// Class:      caloInspector
// 
/**\class caloInspector caloInspector.cc L1TriggerDPGUpgrade/caloInspector/plugins/caloInspector.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Christopher Anelli  
//         Created:  Fri, 16 May 2014 04:20:05 GMT
// $Id$
//
//

//hoMuonAnalyzer header file
#include "Analysis/hoTriggerAnalyzer/interface/hoMuonAnalyzer.h"

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

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/HcalRecHit/interface/HORecHit.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
//#include "DataFormats/HcalDetId/interface/HcalDetId.h"

#include "Analysis/hoTriggerAnalyzer/interface/HistogramBuilder.h"
#include "Analysis/hoTriggerAnalyzer/interface/CommonFunctions.h"

#include <vector>
#include <iostream>
#include "math.h"

using namespace::std;

bool hoBelowThreshold(HORecHit horeco);
//bool isMipMatch(l1extra::L1MuonParticle l1muon, vector<HORecHit> & hoRecoHitsAboveThreshold);
bool isInsideRCut(float eta1, float eta2, float phi1, float phi2);
//float WrapCheck(float phi1, float phi2);

hoMuonAnalyzer::hoMuonAnalyzer(const edm::ParameterSet& iConfig){
  
  //now do what ever initialization is needed
 
  //Get Input Tags from the Configuration
  
  _genInput = iConfig.getParameter<edm::InputTag>("genSrc");
  _l1MuonInput = iConfig.getParameter<edm::InputTag>("l1MuonSrc");
  _horecoInput = iConfig.getParameter<edm::InputTag>("horecoSrc");
  //_stdMuInput = iConfig.getParameter<edm::InputTag>("stdMuSrc");
  _hltSumAODInput = iConfig.getParameter<edm::InputTag>("hltSumAODSrc");
  //_hltTriggerResults = iConfig.getParameter<edm::InputTag>("hltResultsSrc");
  // m_l1GtTmLInputTag = iConfig.getParameter<edm::InputTag> ("L1GtTmLInputTag");
  
  m_nameAlgTechTrig="L1_SingleMu7";
  //m_nameAlgTechTrig="L1_AlwaysTrue";
  
  //edm::Service<TFileService> _fileService;
    

  //initializeHistograms();

  defineTriggersOfInterest();
	
  /*
   * Vector to store hold unique L1MuonPt values, for variable binning.
   */
  //listL1MuonPt = new vector();

  
}


hoMuonAnalyzer::~hoMuonAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  
  //f->Close();
  
}


//
// member functions
//

// ------------ method called for each event  ------------
void
hoMuonAnalyzer::analyze(const edm::Event& iEvent, 
		       const edm::EventSetup& iSetup)
{
   using namespace edm;
   
   /*
    * Get Event Data and Event Setup
    */

   Handle<reco::GenParticleCollection> truthParticles;
   iEvent.getByLabel(_genInput,truthParticles);

   Handle<l1extra::L1MuonParticleCollection> l1Muons;
   iEvent.getByLabel(_l1MuonInput, l1Muons);

   //Handle<reco::TrackCollection> standAloneMuons;
   //iEvent.getByLabel(_stdMuInput,standAloneMuons);

   Handle<HORecHitCollection> hoRecoHits;
   iEvent.getByLabel(_horecoInput, hoRecoHits);

   Handle<trigger::TriggerEvent> aodTriggerEvent;
   iEvent.getByLabel(_hltSumAODInput, aodTriggerEvent);

   InputTag _hltTrigResultsLabel = edm::InputTag("TriggerResults", "", "HLT");
   Handle<TriggerResults> hltTriggerResults;
   iEvent.getByLabel(_hltTrigResultsLabel, hltTriggerResults);
   const TriggerNames & names = iEvent.triggerNames(*hltTriggerResults);
   //cout << "Number of names " << names.size() << endl;
   
   //For first Event List Triggers
   /*
      for(int i =0; i < (int) names.size(); i++){
     cout << names.triggerName(i) << endl;
   }
   */
   
   ESHandle<CaloGeometry> caloGeo;
   iSetup.get<CaloGeometryRecord>().get(caloGeo);


   /*
    * Set Up Level 1 Global Trigger Utility
    */
         
   bool useL1EventSetup = true;
   bool useL1GtTriggerMenuLite = true;

   m_l1GtUtils.getL1GtRunCache(iEvent, iSetup, useL1EventSetup, 
			       useL1GtTriggerMenuLite);
   //cout << "L1 Gt Readout Record From Provenance: " 
   //<< m_l1GtUtils.provL1GtReadoutRecordInputTag() << endl;
   //cout << "L1 Gt Record From Provenance: " 
   //<< m_l1GtUtils.provL1GtRecordInputTag() << endl;

   //m_l1GtUtils.getL1GtRunCache(iEvent, iSetup, useL1EventSetup, 
   //		       useL1GtTriggerMenuLite, m_l1GtTmLInputTag);
   
   //bool l1GtUtilsIsValid = m_l1GtUtils.isValid();
   
   //m_l1GtUtils.retrieveL1GtTriggerMenuLite(iEvent, m_l1GtTmLinputTag);
   

   /*
    *
    *  Start of Analysis
    *
    */

   histogramBuilder.fillCountHistogram("Events");
   
   

   /*
    * Level 1 Muons
    */

   string l1muon_key = "L1Muon";

   auto bl1Muon = l1Muons->cbegin();
   auto el1Muon = l1Muons->cend();

   for( ; bl1Muon != el1Muon; ++bl1Muon ) {
     histogramBuilder.fillCountHistogram(l1muon_key);
     histogramBuilder.fillL1MuonPtHistograms(bl1Muon->pt(), l1muon_key);
     histogramBuilder.fillEtaPhiHistograms(bl1Muon->eta(), bl1Muon->phi(), 
					   l1muon_key);
     //fillEtaPhiHistograms(bl1Muon->eta(), bl1Muon->phi(), l1muon_key);
     //For variable binning
     listL1MuonPt.push_back(bl1Muon->pt());
   }

   /*
    * HLT Triggers
    */

   
   //map<string, string>::iterator itHltNames;
   
   // Loop Over all HLT Trigger Names of Interst
   auto bHLTNames = hltNamesOfInterest.begin();
   for(; bHLTNames != hltNamesOfInterest.end(); ++bHLTNames){
     string hlt_key = bHLTNames->first;
     string hltName = bHLTNames->second;
     for(int hltNameIndex =0; hltNameIndex < (int) names.size(); hltNameIndex++){
       if(names.triggerName(hltNameIndex)==hltName){
	 histogramBuilder.fillTrigHistograms(hltTriggerResults->accept(hltNameIndex), hlt_key);
	 if(hltTriggerResults->accept(hltNameIndex)){
	   histogramBuilder.fillCountHistogram(hlt_key);
	 }
       }
     }
   }

   /*
   for(int i =0; i< (int) hltTriggerResults->size(); i++) {
     if(!(hltTriggerResults->accept(i))) continue;
     std::string name(names.triggerName(i));
     if(name == hltIsoMu24Name){
       histogramBuilder.fillCountHistogram(hltIsoMu24_key);
       cout << name << endl;
     }
   */
	
	

   // cout << "Trigger Results size " << hltTriggerResults->size() << endl;
   //std::vector<std::string> hltNames = hltTriggerResults->getTriggerNames();
   //cout << " HLT names size " << hltNames.size() << endl;
   
   

   /*
    * HLT Filters
    */
   
   //Define TriggerObject Collection
   trigger::TriggerObjectCollection hltAllObjects = aodTriggerEvent->getObjects();
    
   auto bHLTFilters = hltFiltersOfInterest.begin();
   for(; bHLTFilters != hltFiltersOfInterest.end(); ++bHLTFilters){
     string hlt_key = bHLTFilters->first;
     InputTag hltLastFilterTag = bHLTFilters->second;
     //find the filter index
     size_t filterIndex = aodTriggerEvent->filterIndex(hltLastFilterTag);  //(hltMu5LastFilterTag);
     if(filterIndex < aodTriggerEvent->sizeFilters()){
       //cout << "Trigger Object should be Present" << endl;
       const trigger::Keys &keys = aodTriggerEvent->filterKeys(filterIndex);
       for(size_t j=0; j <keys.size(); j++){
	 //cout << "Size of Keys"<< keys.size()  << endl;
	 histogramBuilder.fillEtaPhiHistograms(hltAllObjects[keys[j]].eta(),
					       hltAllObjects[keys[j]].phi(),
					       hlt_key);
	 histogramBuilder.fillPtHistograms(hltAllObjects[keys[j]].pt(),
					   hlt_key);
       }
     }
     
     /*
      * Sub Code to see which triggers were in an event,
      * with a selected filter.
     for(int i =0; i< (int) hltTriggerResults->size(); i++) {
       if(!(hltTriggerResults->accept(i))) continue;
       std::string name(names.triggerName(i));
       //cout << name << endl;
     }
     */
   }
   

   /*
   //For each HLT Filter in the Map, compare with all Filters in the HLT Trigger Event.
   for(itFilters=filtersMap.begin(); itFilters!=filtersMap.end(); ++itFilters){
     for(int i =0; i < aodTriggerEvent->sizeFilters(); i++){
       //Loop over all the HLT Triggers
       if(aodTriggerEvent->filterTag(i).label()!= itFilters->second) continue;
       
       //histogramBuilder.fillEtaPhiHistograms(
       cout << itFilters->first << endl; 
       cout << aodTriggerEvent->filterTag(i).label() << endl;
       trigger::Keys triggerKeys = aodTriggerEvent->filterKeys(i);
       //auto btriggerKeys = triggerKeys.begin();
       cout << triggerKeys.size() << endl;
       for(unsigned int j =0; j < triggerKeys.size(); j++){
	 cout << "The Trigger Key at " << j << " " << triggerKeys[j] << endl;
       }
     }
   } 
   */


   /*
    * HO Reco Hits
    */

   string horeco_key = "horeco";

   //cout << hoRecoHits.size() << endl;
   auto bho_reco = hoRecoHits->begin();
   auto eho_reco = hoRecoHits->end();
   histogramBuilder.fillCountHistogram(horeco_key);
   for(; bho_reco != eho_reco; ++bho_reco){
     //cout << caloGeo->getPosition(bho_reco->id()).eta() << endl;
     
     //h1HORecoEnergy->HistogramBuilder.Fill( bho_reco->energy());

     histogramBuilder.fillEnergyHistograms(bho_reco->energy(), horeco_key);

     float ho_eta, ho_phi;
     ho_eta = caloGeo->getPosition(bho_reco->id()).eta();
     ho_phi = caloGeo->getPosition(bho_reco->id()).phi();
     histogramBuilder.fillEtaPhiHistograms(ho_eta, ho_phi, horeco_key);
   }


   /*
    * L1 Trigger Decisions
    */
   
   // Select on events that pass a specific L1Trigger Decision
   int iErrorCode = -1; 
   bool trigDecision = m_l1GtUtils.decision(iEvent, m_nameAlgTechTrig, iErrorCode);
   //cout << "Error Code: " << iErrorCode << " and Trigger Decision: " << trigDecision << endl;
   
   if(iErrorCode == 0){
     histogramBuilder.fillTrigHistograms(trigDecision,m_nameAlgTechTrig);
     if(trigDecision){
       histogramBuilder.fillCountHistogram(m_nameAlgTechTrig);
     }  
     
   } else if (iErrorCode == 1) {
     cout<< "trigger does not exist in the L1 menu" << endl;
   } else {
     // error - see error code
     cout << "Error Code " << iErrorCode;
   }
   
  /*
     bool decisionBeforeMaskAlgTechTrig = m_l1GtUtils.decisionBeforeMask(iEvent, 
     m_nameAlgTechTrig, 
     iErrorCode);
   */

   /*
    *Stand Alone Muons
    */

   /*
   auto bstandAloneMuon = standAloneMuons->cbegin();
   auto estandAloneMuon = standAloneMuons->cend();
   
   for( ; bstandAloneMuon != estandAloneMuon; ++bstandAloneMuon ) {
     h1SAMuonPt->HistogramBuilder.Fill(bstandAloneMuon->pt());
      }
   */


   /*
    * HO Rec Hits Above Threshold
    */
   
   string horecoT_key ="horecoAboveThreshold";

   //Filter out HO Rec Hits below Threshold.
    
   std::vector<HORecHit> hoRecoHitsAboveThreshold (hoRecoHits->size());
   
   auto it =std::remove_copy_if (hoRecoHits->begin(), hoRecoHits->end(),
				 hoRecoHitsAboveThreshold.begin(), 
				 hoBelowThreshold);
   hoRecoHitsAboveThreshold.resize(std::distance(hoRecoHitsAboveThreshold.begin(),it));
   
   //Handle<HORecHitCollection> hoRecHitsAboveThreshold;
   //HORecHitCollection* hoRecHitsAboveThreshold = new HORecHitCollection(hoRecoHits->size());
   
   auto bho_recoT = hoRecoHitsAboveThreshold.begin();
   auto eho_recoT = hoRecoHitsAboveThreshold.end();
   
   for(; bho_recoT != eho_recoT; ++bho_recoT){
     histogramBuilder.fillCountHistogram(horecoT_key);
     histogramBuilder.fillEnergyHistograms(bho_recoT->energy(), horecoT_key);

     float hoT_eta, hoT_phi;
     hoT_eta = caloGeo->getPosition(bho_recoT->id()).eta();
     hoT_phi = caloGeo->getPosition(bho_recoT->id()).phi();
     histogramBuilder.fillEtaPhiHistograms(hoT_eta, hoT_phi, horecoT_key);
   }

   /*
    * L1 Muons Matched to a MIP
    */

   string l1MuonMipMatch_key = "L1MuonwithMipMatch";
   
   bl1Muon = l1Muons->cbegin();
   el1Muon = l1Muons->cend();
   
   for( ; bl1Muon != el1Muon; ++bl1Muon ){

     //bool isMipMatch=checkMipMatch();
     bho_recoT = hoRecoHitsAboveThreshold.begin();
     eho_recoT = hoRecoHitsAboveThreshold.end();
   
     bool mipMatch = false;
     for(; bho_recoT != eho_recoT; ++bho_recoT){
       float l1Muon_eta, horeco_eta, l1Muon_phi, horeco_phi;
       l1Muon_eta = bl1Muon->eta();
       l1Muon_phi = bl1Muon->phi();
       horeco_eta = caloGeo->getPosition(bho_recoT->id()).eta();
       horeco_phi = caloGeo->getPosition(bho_recoT->id()).phi();

       string l1MuonhoReco_key = "L1MuonandHOReco";
       histogramBuilder.fillDeltaEtaDeltaPhiHistograms(l1Muon_eta, horeco_eta,
						       l1Muon_phi, horeco_phi,
						       l1MuonhoReco_key);
       if(isInsideRCut(l1Muon_eta, horeco_eta, l1Muon_phi, horeco_phi)){
	 mipMatch=true; //Only need a single match
	 //NB It is possible for there to be more than one matched Mip.
	 string hoRecoMipMatch_key = "HORecowithMipMatch";
	 histogramBuilder.fillCountHistogram(hoRecoMipMatch_key);
	 histogramBuilder.fillEtaPhiHistograms(caloGeo->getPosition(bho_recoT->id()).eta(),
					       caloGeo->getPosition(bho_recoT->id()).phi(), 
					       hoRecoMipMatch_key);
	 histogramBuilder.fillEnergyHistograms(bho_recoT->energy(),hoRecoMipMatch_key);
	
	 string l1MuonhoRecomipMatch_key = "L1MuonandHORecowithMipMatch";
	 histogramBuilder.fillDeltaEtaDeltaPhiHistograms(l1Muon_eta,horeco_eta,
							 l1Muon_phi, horeco_phi,
							 l1MuonhoRecomipMatch_key);
       }
     }
	 
     if(mipMatch){
       histogramBuilder.fillCountHistogram(l1MuonMipMatch_key);
       histogramBuilder.fillL1MuonPtHistograms(bl1Muon->pt(), l1MuonMipMatch_key);
       histogramBuilder.fillEtaPhiHistograms(bl1Muon->eta(), bl1Muon->phi(), l1MuonMipMatch_key);
     }
   }
}


// ------------ method called once each job just before starting event loop  ------------
void 
hoMuonAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
hoMuonAnalyzer::endJob() 
{
  //h1SAMuonPt->Write();
  //histogramBuilder.fillCountHistogram(eventCounter,"Events");
  //h1EventCounter->Write();
}

// ------------ method called when starting to processes a run  ------------

void 
hoMuonAnalyzer::beginRun(const edm::Run& iRun, 
			   const edm::EventSetup& evSetup)
{

  /*
   *Exploring Provenance
   */

  /*
  typedef std::vector<edm::Provenance const*> Provenances;
  Provenances provenances;
  iRun.getAllProvenance(provenances);
  for (Provenances::iterator itProv = provenances.begin(), itProvEnd =
         provenances.end(); itProv != itProvEnd; ++itProv) {
    cout << endl << (*itProv)->friendlyClassName() << ", " << (*itProv)->moduleLabel() << ", "
         << (*itProv)->productInstanceName() << ", " << (*itProv)->processName() << endl;
  }
  */

  /*
   * So that code runs faster. Trigger quantities that can be cached in
   * beginRun, are not cached again in analyze.  (Handled by Utility) 
   */

  // edm::InputTag l1GtTriggerMenuLiteInputTag;
  //m_l1GtUtils.getL1GtTriggerMenuLiteInputTag(iRun, l1GtTriggerMenuLiteInputTag);
  //cout << "Get L1 GT Trigger Menu Lite Input Tag" <<  l1GtTriggerMenuLiteInputTag << endl;
  
  
  bool useL1EventSetup = true;
  bool useL1GtTriggerMenuLite = true;
  cout << "getL1GtRunCache" << endl;
  cout << "UseL1EventSetup: " << useL1EventSetup << "UseL1GtTriggerMenuLite :" 
       << useL1GtTriggerMenuLite << endl;  
  m_l1GtUtils.getL1GtRunCache(iRun, evSetup, useL1EventSetup, useL1GtTriggerMenuLite);


  

  /*
  cout << "Trigger Menu From Provenance: " 
       << m_l1GtUtils.provL1GtTriggerMenuLiteInputTag() << endl;
  */
  
  //m_l1GtUtils.retrieveL1GtTriggerMenuLite(iRun, _l1GtTmLInputTag);

  /*
   * For the HLT
   */
  
}


// ------------ method called when ending the processing of a run  ------------

void 
hoMuonAnalyzer::endRun(const edm::Run& iRun, const edm::EventSetup& evSetup)
{
  
  //Only interested in unique values
  listL1MuonPt.sort();
  listL1MuonPt.unique();  //NB it is called after sort
  cout <<"The list contains " << listL1MuonPt.size() << "unique entries:";
  std::list<float>::iterator it; 
  for (it=listL1MuonPt.begin(); it!=listL1MuonPt.end(); ++it){
    cout << ' ' << *it;
  }
  cout << endl;
}


// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
hoMuonAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
hoMuonAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
hoMuonAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}



/*
 * Helper Functions for Filter
 */

bool hoBelowThreshold(HORecHit horeco){
  if(horeco.energy() < threshold) return true;
  return false;
}


bool isInsideRCut(float eta1, float eta2, float phi1, float phi2){
  
  float delta_eta, delta_phi;
  CommonFunctions commonFunctions;

  delta_eta = eta1 - eta2;
  delta_phi = commonFunctions.WrapCheck(phi1,phi2); //Finds difference in phi

  //The L1 Muon is compared with all HO Rec Hits above Threshold.
  if(pow(delta_eta,2)+pow(delta_phi,2) <= pow(deltaR_Max,2)) return true;
  return false;
}

void hoMuonAnalyzer::defineTriggersOfInterest(){
  
  /*
   * HLT Triggers
   */

  string hltIsoMu24_key = "hltIsoMu24";
  hltNamesOfInterest.insert(pair<string, string>(hltIsoMu24_key,"HLT_IsoMu24_v18"));
  hltFiltersOfInterest.insert(pair<string, edm::InputTag>(hltIsoMu24_key, 
							  edm::InputTag("hltL3crIsoL1sMu16L1f0L2f16QL3"
									"f24QL3crIsoRhoFiltered0p15",
									"","HLT")));

  string hltMu5_key = "hltMu5";
  hltNamesOfInterest.insert(pair<string, string>(hltMu5_key, "HLT_Mu5_v21"));
  hltFiltersOfInterest.insert(pair<string, edm::InputTag>(hltMu5_key, 
							  edm::InputTag("hltL3fL1sMu3L3Filtered5",
									"","HLT")));

}

//define this as a plug-in
DEFINE_FWK_MODULE(hoMuonAnalyzer);
