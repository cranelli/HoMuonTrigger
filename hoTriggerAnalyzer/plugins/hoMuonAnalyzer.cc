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
#include "HoMuonTrigger/hoTriggerAnalyzer/interface/hoMuonAnalyzer.h"

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

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
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

#include "HoMuonTrigger/hoTriggerAnalyzer/interface/HistogramBuilder.h"
#include "HoMuonTrigger/hoTriggerAnalyzer/interface/CommonFunctions.h"

#include <vector>
#include <iostream>
#include "math.h"

using namespace::std;

bool hoBelowThreshold(HORecHit horeco);
//bool isMipMatch(l1extra::L1MuonParticle l1muon, vector<HORecHit> & hoRecoHitsAboveThreshold);
bool isInsideRCut(float deltaR_Max, float eta1, float eta2, float phi1, float phi2);
//float WrapCheck(float phi1, float phi2);

hoMuonAnalyzer::hoMuonAnalyzer(const edm::ParameterSet& iConfig){
  
  //now do what ever initialization is needed
 
  //Get Input Tags from the Configuration

  _genInfoInput = iConfig.getParameter<edm::InputTag>("genInfoSrc");
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

   Handle<GenEventInfoProduct> genEventInfo;
   iEvent.getByLabel(_genInfoInput, genEventInfo);

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
   
   //List Triggers in Event
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

     
   /*
    * Generator Weight Information
    */

   string genEventInfo_key = "GenEventInfo";
   string genEventInfoNoWeight_key = "GenEventInfo_NoWeight";
   weight = genEventInfo->weight();
   //cout << "Event Generator Weight: " << weight << endl;
   histogramBuilder.fillWeightHistograms(weight, genEventInfo_key, weight);
   
   
   float qScale = genEventInfo->qScale();
   //cout << "Event Generator Qscale" << qScale << endl;
   histogramBuilder.fillPtHistograms(qScale,genEventInfo_key,weight);
   histogramBuilder.fillPtHistograms(qScale,genEventInfoNoWeight_key, weight);
   

   /*
    * Event Count
    */

   histogramBuilder.fillCountHistogram("Events",weight);

   /*
    * Level 1 Muons
    */

   string l1muon_key = "L1Muon";

   auto bl1Muon = l1Muons->cbegin();
   auto el1Muon = l1Muons->cend();

   for( ; bl1Muon != el1Muon; ++bl1Muon ) {
     histogramBuilder.fillCountHistogram(l1muon_key, weight);
     histogramBuilder.fillL1MuonPtHistograms(bl1Muon->pt(),
					     l1muon_key, weight);
     histogramBuilder.fillEtaPhiHistograms(bl1Muon->eta(), bl1Muon->phi(), 
					   l1muon_key, weight);
     //fillEtaPhiHistograms(bl1Muon->eta(), bl1Muon->phi(), l1muon_key);
     //For variable binning
     listL1MuonPt.push_back(bl1Muon->pt());
   }

   /*
    * Level 1 Muons in the Barrel
    */
   
   string l1muonB_key = "L1MuonBarrel";

   auto bl1MuonB = l1Muons->cbegin();

   for( ; bl1MuonB != l1Muons->cend(); ++bl1MuonB ) {
     if(abs(bl1MuonB->eta())>barrel_eta) continue;  // Only want L1 Muons in the Barrel
     histogramBuilder.fillCountHistogram(l1muonB_key,weight);
     histogramBuilder.fillL1MuonPtHistograms(bl1MuonB->pt(), l1muonB_key, weight);
     histogramBuilder.fillEtaPhiHistograms(bl1MuonB->eta(), bl1MuonB->phi(), 
					   l1muonB_key, weight);
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
	 histogramBuilder.fillTrigHistograms(hltTriggerResults->accept(hltNameIndex),
					     hlt_key, weight);
	 if(hltTriggerResults->accept(hltNameIndex)){
	   histogramBuilder.fillCountHistogram(hlt_key, weight);
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
    * HLT Filters and Trigger Objects
    */
   
   //Define TriggerObject Collection
   trigger::TriggerObjectCollection hltAllObjects = aodTriggerEvent->getObjects();
  
   //Store the Trigger Objects, from the Filters Of Interest
   auto bHLTFilters = hltFiltersOfInterest.begin();
   for(; bHLTFilters != hltFiltersOfInterest.end(); ++bHLTFilters){
     string hlt_key = bHLTFilters->first;
     InputTag hltLastFilterTag = bHLTFilters->second;
     hltTriggerObjects.clear();   //Each new filters, starts empty  
     
     //find the filter index
     size_t filterIndex = aodTriggerEvent->filterIndex(hltLastFilterTag);  
     if(filterIndex < aodTriggerEvent->sizeFilters()){  //Check that Filter is in Trigger Summary
       //cout << "Trigger Object should be Present" << endl;
       const trigger::Keys &keys = aodTriggerEvent->filterKeys(filterIndex);
       for(size_t j=0; j <keys.size(); j++){
	  //cout << "Size of Keys"<< keys.size()  << endl;
	 hltTriggerObjects.push_back(hltAllObjects[keys[j]]);  //Get Trigger Object using Keys
       }
     }
     //cout << hlt_key << "Trigger Objects Size" << hltTriggerObjects.size() << endl;
     hltTriggerObjectsOfInterest[hlt_key] = hltTriggerObjects;
     /*
      * Sub Code to see which triggers were in an event,
      * with a selected filter.
      for(int i =0; i< (int) hltTriggerResults->size(); i++) {
      if(!(hltTriggerResults->accept(i))) continue;
       std::string name(names.triggerName(i));
       //cout << name << endl;
       }
     */
     
     //Fill Histograms using selected Trigger Objects
     auto bTriggerObject = hltTriggerObjects.begin();
     for(; bTriggerObject != hltTriggerObjects.end(); ++bTriggerObject){
	
       histogramBuilder.fillEtaPhiHistograms(bTriggerObject->eta(),
					     bTriggerObject->phi(),
					     hlt_key, weight);
       histogramBuilder.fillPtHistograms(bTriggerObject->pt(),
					 hlt_key, weight);

       //Barrel HLT Trigger Objects
       if(abs(bTriggerObject->eta())<barrel_eta){
	 stringstream hltB_key_stream;
	 hltB_key_stream << hlt_key <<"_Barrel" << endl;
	 string hlt_Barrel_key = hltB_key_stream.str();
	 //cout << hlt_Barrel_key << endl;
	 histogramBuilder.fillEtaPhiHistograms(bTriggerObject->eta(),
					       bTriggerObject->phi(),
					       hlt_Barrel_key, weight);
	 histogramBuilder.fillPtHistograms(bTriggerObject->pt(),
					   hlt_Barrel_key, weight);
       }
     }
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
   histogramBuilder.fillCountHistogram(horeco_key, weight);
   for(; bho_reco != eho_reco; ++bho_reco){
     //cout << caloGeo->getPosition(bho_reco->id()).eta() << endl;
     
     //h1HORecoEnergy->HistogramBuilder.Fill( bho_reco->energy());

     histogramBuilder.fillEnergyHistograms(bho_reco->energy(), horeco_key, weight);

     float ho_eta, ho_phi;
     ho_eta = caloGeo->getPosition(bho_reco->id()).eta();
     ho_phi = caloGeo->getPosition(bho_reco->id()).phi();
     histogramBuilder.fillEtaPhiHistograms(ho_eta, ho_phi, horeco_key, weight);
   }


   /*
    * L1 Trigger Decisions
    */
   
   // Select on events that pass a specific L1Trigger Decision
   int iErrorCode = -1; 
   bool trigDecision = m_l1GtUtils.decision(iEvent, m_nameAlgTechTrig, iErrorCode);
   //cout << "Error Code: " << iErrorCode << " and Trigger Decision: " << trigDecision << endl;
   
   if(iErrorCode == 0){
     histogramBuilder.fillTrigHistograms(trigDecision,m_nameAlgTechTrig, weight);
     if(trigDecision){
       histogramBuilder.fillCountHistogram(m_nameAlgTechTrig, weight);
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
     histogramBuilder.fillCountHistogram(horecoT_key, weight);
     histogramBuilder.fillEnergyHistograms(bho_recoT->energy(), horecoT_key, weight);

     float hoT_eta, hoT_phi;
     hoT_eta = caloGeo->getPosition(bho_recoT->id()).eta();
     hoT_phi = caloGeo->getPosition(bho_recoT->id()).phi();
     histogramBuilder.fillEtaPhiHistograms(hoT_eta, hoT_phi, horecoT_key, weight);
   }

   /*
    * L1 Muons Matched to a MIP
    */

   string l1MuonBMipMatch_key = "L1MuonBarrelwithMipMatch";
   
   l1MuonsBarrelWithMip.clear();  //Initialize Vector to store l1muons with mip.

   for(int l1Muon_index =0; l1Muon_index < (int) l1Muons->size(); l1Muon_index++){
     l1extra::L1MuonParticle bl1MuonB = l1Muons->at(l1Muon_index);
     if(abs(bl1MuonB.eta())>barrel_eta) continue;  //Select Barrel L1 Muons
     
     
     float l1MuonB_phi, l1MuonB_eta;
     l1MuonB_eta = bl1MuonB.eta();
     l1MuonB_phi = bl1MuonB.phi();
     
     //bool isMipMatch=checkMipMatch();
     bho_recoT = hoRecoHitsAboveThreshold.begin();
     eho_recoT = hoRecoHitsAboveThreshold.end();
     bool mipMatch = false;
     for(; bho_recoT != eho_recoT; ++bho_recoT){
     
       float horeco_eta, horeco_phi;
       horeco_eta = caloGeo->getPosition(bho_recoT->id()).eta();
       horeco_phi = caloGeo->getPosition(bho_recoT->id()).phi();

       string l1MuonBhoReco_key = "L1MuonBarrelandHOReco";
       histogramBuilder.fillDeltaEtaDeltaPhiHistograms(l1MuonB_eta, horeco_eta,
						       l1MuonB_phi, horeco_phi,
						       l1MuonBhoReco_key, weight);

       //Make the Mip Match
       if(isInsideRCut(deltaRMip_Max, l1MuonB_eta, horeco_eta, l1MuonB_phi, horeco_phi)){
	 mipMatch=true; //Only need a single match
	 //NB It is possible for there to be more than one matched Mip.
	 string hoRecoMipMatch_key = "HORecowithMipMatch_All";
	 histogramBuilder.fillCountHistogram(hoRecoMipMatch_key, weight);
	 histogramBuilder.fillEtaPhiHistograms(horeco_eta, horeco_phi, 
					       hoRecoMipMatch_key,weight);
	 histogramBuilder.fillEnergyHistograms(bho_recoT->energy(),hoRecoMipMatch_key, weight);
	
	 string l1MuonBhoRecoMatch_key = "L1MuonBarrelandHORecowithMipMatch_All";
	 histogramBuilder.fillDeltaEtaDeltaPhiHistograms(l1MuonB_eta,horeco_eta,
							 l1MuonB_phi, horeco_phi,
							 l1MuonBhoRecoMatch_key,
							 weight);
       }
     }
	 
     if(mipMatch){
       //Store L1 Muons with Mips
       l1MuonsBarrelWithMip.push_back(bl1MuonB);
       //Make Histograms
       histogramBuilder.fillCountHistogram(l1MuonBMipMatch_key, weight);
       histogramBuilder.fillL1MuonPtHistograms(bl1MuonB.pt(), l1MuonBMipMatch_key,weight);
       histogramBuilder.fillEtaPhiHistograms(l1MuonB_eta, l1MuonB_phi,
					     l1MuonBMipMatch_key, weight);
     }
   }
   //cout << "Number of L1 Muons with Mips: " << l1MuonsWithMip.size() << endl;


   /*
    * L1 Muons Matched to HLT Muon Trigger Object
    */
   
   bl1MuonB = l1Muons->cbegin();
   for( ; bl1MuonB != l1Muons->cend(); ++bl1MuonB ){
     if(abs(bl1MuonB->eta())>barrel_eta) continue;  //Select Barrel L1 Muons
     hltTriggerObjects = hltTriggerObjectsOfInterest["hltMu5"];  //HLT 5 Gev Pt Muon
     //cout << "Size of HLT Muon Trigger Objects: " << hltTriggerObjects.size() << endl;

     auto bhltMuon = hltTriggerObjects.begin();
     
     bool hltMatch = false;
     for(; bhltMuon != hltTriggerObjects.end(); ++bhltMuon){
     
       float l1MuonB_eta, hltMuon_eta, l1MuonB_phi, hltMuon_phi;
       l1MuonB_eta = bl1MuonB->eta();
       l1MuonB_phi = bl1MuonB->phi();
       hltMuon_eta = bhltMuon->eta();
       hltMuon_phi = bhltMuon->phi();
       string l1MuonBhltMuon_key = "L1MuonBarrelandHLTMuon";
       histogramBuilder.fillDeltaEtaDeltaPhiHistograms(l1MuonB_eta, hltMuon_eta,
						       l1MuonB_phi, hltMuon_phi,
						       l1MuonBhltMuon_key, weight);
       if(isInsideRCut(deltaRHLT_Max,l1MuonB_eta, hltMuon_eta, l1MuonB_phi, hltMuon_phi)){
	 hltMatch=true; //Only need a single match
	 //NB It is possible for there to be more than one matched hlt.
	 string hltMuonhltMatch_key = "HltMuonwithHLTMatch_All";
	 histogramBuilder.fillCountHistogram(hltMuonhltMatch_key, weight);
	 histogramBuilder.fillEtaPhiHistograms(hltMuon_eta,
					       hltMuon_phi,
					       hltMuonhltMatch_key, weight);
	 histogramBuilder.fillPtHistograms(bhltMuon->pt(),hltMuonhltMatch_key, weight);
       
	 string l1MuonBhltMuonHltMatch_key = "L1MuonBarrelandHLTMuonwithHltMatch_All";
	 histogramBuilder.fillDeltaEtaDeltaPhiHistograms(l1MuonB_eta,hltMuon_eta,
							 l1MuonB_phi, hltMuon_phi,
							 l1MuonBhltMuonHltMatch_key,
							 weight);
       }
     }
     
     if(hltMatch){
       string l1MuonBHltMatch_key = "L1MuonBarrelwithHLTMatch";
       histogramBuilder.fillCountHistogram(l1MuonBHltMatch_key, weight);
       histogramBuilder.fillL1MuonPtHistograms(bl1MuonB->pt(), l1MuonBHltMatch_key, weight);
       histogramBuilder.fillEtaPhiHistograms(bl1MuonB->eta(), bl1MuonB->phi(), 
					     l1MuonBHltMatch_key, weight);
       
     }
   }


   /*
    * L1 Muons with Mip Matched to HLT Muon Trigger Object
    */
   
   auto bl1MuonBwMip = l1MuonsBarrelWithMip.cbegin();
   
   for( ; bl1MuonBwMip != l1MuonsBarrelWithMip.cend(); ++bl1MuonBwMip ){
     if(abs(bl1MuonBwMip->eta())>barrel_eta) continue;  //Select Barrel L1 Muons
     bool hltMatch = false;
     float l1MuonBwMip_eta, l1MuonBwMip_phi;
     l1MuonBwMip_eta = bl1MuonBwMip->eta();
     l1MuonBwMip_phi = bl1MuonBwMip->phi();
     hltTriggerObjects = hltTriggerObjectsOfInterest["hltMu5"];  //HLT 5 Gev Pt Muon
     //cout << "Size of HLT Muon Trigger Objects: " << hltTriggerObjects.size() << endl;

     auto bhltMuon = hltTriggerObjects.begin();
     for(; bhltMuon != hltTriggerObjects.end(); ++bhltMuon){
     
       float hltMuon_eta, hltMuon_phi;
       
       hltMuon_eta = bhltMuon->eta();
       hltMuon_phi = bhltMuon->phi();
       string l1MuonBwMiphltMuon_key = "L1MuonBarrelMipandHLTMuon";
       histogramBuilder.fillDeltaEtaDeltaPhiHistograms(l1MuonBwMip_eta, hltMuon_eta,
						       l1MuonBwMip_phi, hltMuon_phi,
						       l1MuonBwMiphltMuon_key, weight);
       if(isInsideRCut(deltaRHLT_Max,l1MuonBwMip_eta, hltMuon_eta, l1MuonBwMip_phi, hltMuon_phi)){
	 hltMatch=true; //Only need a single match
	 //NB It is possible for there to be more than one matched hlt.
	 string hltMuonhltMatchwMip_key = "HltMuonwithHLTMatchwithL1BarrelMip_All";
	 histogramBuilder.fillCountHistogram(hltMuonhltMatchwMip_key, weight);
	 histogramBuilder.fillEtaPhiHistograms(hltMuon_eta,
					       hltMuon_phi,
					       hltMuonhltMatchwMip_key, weight);
	 histogramBuilder.fillPtHistograms(bhltMuon->pt(),hltMuonhltMatchwMip_key, weight);
       
	 string l1MuonBwMiphltMuonMatch_key = "L1MuonBarrelMipandHLTMuonwithHltMatch_All";
	 histogramBuilder.fillDeltaEtaDeltaPhiHistograms(l1MuonBwMip_eta,hltMuon_eta,
							 l1MuonBwMip_phi, hltMuon_phi,
							 l1MuonBwMiphltMuonMatch_key,
							 weight);
       }
     }
     
     if(hltMatch){
       string l1MuonBwMipHltMatch_key = "L1MuonBarrelwithMipHLTMatch";
       histogramBuilder.fillCountHistogram(l1MuonBwMipHltMatch_key, weight);
       histogramBuilder.fillL1MuonPtHistograms(bl1MuonBwMip->pt(), l1MuonBwMipHltMatch_key,
					       weight);
       histogramBuilder.fillEtaPhiHistograms(l1MuonBwMip_eta, l1MuonBwMip_phi, 
					     l1MuonBwMipHltMatch_key, weight);
       
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


bool isInsideRCut(float deltaR_Max, float eta1, float eta2, float phi1, float phi2){
  
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

  /*
  string hltIsoMu24_key = "hltIsoMu24";
  hltNamesOfInterest.insert(pair<string, string>(hltIsoMu24_key,"HLT_IsoMu24_v18"));
  hltFiltersOfInterest.insert(pair<string, edm::InputTag>(hltIsoMu24_key, 
							  edm::InputTag("hltL3crIsoL1sMu16L1f0L2f16QL3"
									"f24QL3crIsoRhoFiltered0p15",
									"","HLT")));
  */

  string hltMu5_key = "hltMu5";
  hltNamesOfInterest.insert(pair<string, string>(hltMu5_key, "HLT_Mu5_v21"));
  hltFiltersOfInterest.insert(pair<string, edm::InputTag>(hltMu5_key, 
							  edm::InputTag("hltL3fL1sMu3L3Filtered5",
									"","HLT")));

}

//define this as a plug-in
DEFINE_FWK_MODULE(hoMuonAnalyzer);
