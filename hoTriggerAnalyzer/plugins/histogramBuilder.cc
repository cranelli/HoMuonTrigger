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

//histogramBuilder header file
#include "Analysis/hoTriggerAnalyzer/interface/histogramBuilder.h"

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

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
//#include "DataFormats/HcalDetId/interface/HcalDetId.h"

#include "TH1F.h"
#include "TH2.h"

#include <vector>
#include <iostream>
#include "math.h"

using namespace::std;

bool hoBelowThreshold(HORecHit horeco);
bool isInsideRCut(float delta_eta, float delta_phi);
double WrapCheck(float phi1, float phi2);

histogramBuilder::histogramBuilder(const edm::ParameterSet& iConfig){
  
  //now do what ever initialization is needed
 
  //Get Input Tags from the Configuration
  
  _genInput = iConfig.getParameter<edm::InputTag>("genSrc");
  _l1MuonInput = iConfig.getParameter<edm::InputTag>("l1MuonSrc");
  _horecoInput = iConfig.getParameter<edm::InputTag>("horecoSrc");
  //_stdMuInput = iConfig.getParameter<edm::InputTag>("stdMuSrc");
  
  // m_l1GtTmLInputTag = iConfig.getParameter<edm::InputTag> ("L1GtTmLInputTag");
  
  m_nameAlgTechTrig="L1_SingleMu7";
  //m_nameAlgTechTrig="L1_AlwaysTrue";
  
  //edm::Service<TFileService> _fileService;
    
  //eventCounter = 0;

  //initializeHistograms();
	
  /*
   * Vector to store hold unique L1MuonPt values, for variable binning.
   */

  //listL1MuonPt = new vector();
  
}


histogramBuilder::~histogramBuilder()
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
histogramBuilder::analyze(const edm::Event& iEvent, 
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

   fillCountHistogram("Events");
   
   

   /*
    * Level 1 Muons
    */

   string l1muon_key = "L1Muon";

   auto bl1Muon = l1Muons->cbegin();
   auto el1Muon = l1Muons->cend();

   for( ; bl1Muon != el1Muon; ++bl1Muon ) {
     fillL1MuonPtHistograms(bl1Muon->pt(), l1muon_key);
     fillEtaPhiHistograms(bl1Muon->eta(), bl1Muon->phi(), l1muon_key);
     //For variable binning
     listL1MuonPt.push_back(bl1Muon->pt());
   }
   

   /*
    * HO Reco Hits
    */

   string horeco_key = "horeco";

   //cout << hoRecoHits.size() << endl;
   auto bho_reco = hoRecoHits->begin();
   auto eho_reco = hoRecoHits->end();
   
   for(; bho_reco != eho_reco; ++bho_reco){
     //cout << caloGeo->getPosition(bho_reco->id()).eta() << endl;
     
     //h1HORecoEnergy->Fill( bho_reco->energy());

     fillEnergyHistograms(bho_reco->energy(), horeco_key);

     float ho_eta, ho_phi;
     ho_eta = caloGeo->getPosition(bho_reco->id()).eta();
     ho_phi = caloGeo->getPosition(bho_reco->id()).phi();
     fillEtaPhiHistograms(ho_eta, ho_phi, horeco_key);
   }


   /*
    * L1 Trigger Decisions
    */
   
   // Select on events that pass a specific L1Trigger Decision
   int iErrorCode = -1; 
   bool trigDecision = m_l1GtUtils.decision(iEvent, m_nameAlgTechTrig, iErrorCode);
   //cout << "Error Code: " << iErrorCode << " and Trigger Decision: " << trigDecision << endl;
   
   if(iErrorCode == 0){
     fillTrigHistograms(trigDecision,m_nameAlgTechTrig);
     if(trigDecision){
       fillCountHistogram(m_nameAlgTechTrig);
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
     h1SAMuonPt->Fill(bstandAloneMuon->pt());
      }
   */


   /*
    * HO Rec Hits Above Threshold
    */
   
   string horecoAT_key ="horecoAboveThreshold";

   //Filter out HO Rec Hits below Threshold.
    
   std::vector<HORecHit> hoRecoHitsAboveThreshold (hoRecoHits->size());
   
   auto it =std::remove_copy_if (hoRecoHits->begin(), hoRecoHits->end(),
				 hoRecoHitsAboveThreshold.begin(), hoBelowThreshold);
   hoRecoHitsAboveThreshold.resize(std::distance(hoRecoHitsAboveThreshold.begin(), it));
   
   //Handle<HORecHitCollection> hoRecHitsAboveThreshold;
   //HORecHitCollection* hoRecHitsAboveThreshold = new HORecHitCollection(hoRecoHits->size());
   
   auto bho_recoAT = hoRecoHitsAboveThreshold.begin();
   auto eho_recoAT = hoRecoHitsAboveThreshold.end();
   
   for(; bho_recoAT != eho_recoAT; ++bho_recoAT){
    
     fillEnergyHistograms(bho_recoAT->energy(), horecoAT_key);

     float hoAT_eta, hoAT_phi;
     hoAT_eta = caloGeo->getPosition(bho_recoAT->id()).eta();
     hoAT_phi = caloGeo->getPosition(bho_recoAT->id()).phi();
     fillEtaPhiHistograms(hoAT_eta, hoAT_phi, horecoAT_key);
   }

   /*
    * L1 Muons Matched to a MIP
    */

   string l1MuonMipMatch_key = "L1MuonwithMipMatch";

   bl1Muon = l1Muons->cbegin();
   el1Muon = l1Muons->cend();
   
   for( ; bl1Muon != el1Muon; ++bl1Muon ){
     bho_recoAT = hoRecoHitsAboveThreshold.begin();
     eho_recoAT = hoRecoHitsAboveThreshold.end();
   
     for(; bho_recoAT != eho_recoAT; ++bho_recoAT){

       float delta_eta, delta_phi;
       delta_eta = bl1Muon->eta() - caloGeo->getPosition(bho_recoAT->id()).eta();
       delta_phi = WrapCheck(bl1Muon->phi(), caloGeo->getPosition(bho_recoAT->id()).phi());
       
       if(isInsideRCut(delta_eta,delta_phi)){
	 fillL1MuonPtHistograms(bl1Muon->pt(), l1MuonMipMatch_key);
	 fillEtaPhiHistograms(bl1Muon->eta(), bl1Muon->phi(), l1MuonMipMatch_key);
       }
     }
   }
}


// ------------ method called once each job just before starting event loop  ------------
void 
histogramBuilder::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
histogramBuilder::endJob() 
{
  //h1SAMuonPt->Write();
  //fillCountHistogram(eventCounter,"Events");
  //h1EventCounter->Write();
}

// ------------ method called when starting to processes a run  ------------

void 
histogramBuilder::beginRun(const edm::Run& iRun, 
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
  cout << "UseL1EventSetup: " << useL1EventSetup << "UseL1GtTriggerMenuLite :"  << useL1GtTriggerMenuLite << endl;  
  m_l1GtUtils.getL1GtRunCache(iRun, evSetup, useL1EventSetup, useL1GtTriggerMenuLite);
  
  /*
  cout << "Trigger Menu From Provenance: " 
       << m_l1GtUtils.provL1GtTriggerMenuLiteInputTag() << endl;
  */
  
  //m_l1GtUtils.retrieveL1GtTriggerMenuLite(iRun, _l1GtTmLInputTag);
}


// ------------ method called when ending the processing of a run  ------------

void 
histogramBuilder::endRun(const edm::Run& iRun, const edm::EventSetup& evSetup)
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
histogramBuilder::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
histogramBuilder::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
histogramBuilder::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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


bool isInsideRCut(float delta_eta, float delta_phi){
  //The L1 Muon is compared with all HO Rec Hits above Threshold.
  if(pow(delta_eta,2)+pow(delta_phi,2) <= pow(deltaR_Max,2)) return true;
  return false;
}

/*
 * Wrap check calcuates the difference between two phi's,
 * making sure they are not more than 2 pi apart.
 */

 double WrapCheck(float phi1, float phi2){
   //double M_PI = (double) 3.14;                                                                     
   double delta_phi = phi1 - phi2;
   if(delta_phi < -M_PI){
     return (2*M_PI + delta_phi);
   }
   if(delta_phi > M_PI){
     return (delta_phi - 2*M_PI);
   }
   return delta_phi;
 }


/*
 * Helper Functions for Histograms
 */

/*
void histogramBuilder::initializeHistograms(){

   
  h1L1MuonPt = _fileService->make<TH1F>("h1L1MuonPt",
					"L1 Muon Pt distribution",
					2010, -5.0, 1000.0);
 
  h1L1Muon_Trig_Pt = _fileService->make<TH1F>("h1L1Muon_Trig_Pt",
					      "L1 Muon with" 
					      "L1 Trigger Pt distribution",
					      2000, 0.0, 1000.0);
 
  h1HORecoEnergy = _fileService->make<TH1F>("h1HORecoEnergy",
					"HO Reco Energy distribution",
					2100, -5.0, 100.0);
 
  h1Trig = _fileService->make<TH1F>("h1Trig", 
				    "Level 1 Trigger",2, 0,2);

  h1SAMuonPt = _fileService->make<TH1F>("h1SAMuonPt",
					"Stand Alone Muon Pt distribution",
					2000, 0.0, 1000.0);

  h1SAMuon_Trig_Pt = _fileService->make<TH1F>("h1SAMuon_Trig_Pt",
					      "Stand Alone Muon with" 
					      "L1 Trigger Pt distribution",
					      2000, 0.0, 1000.0);
}
*/

/*
 *Generic Form to build Counting Histograms 
 *Type of object (or filtering) specified by the key
 *Fills the 1 bin.
 */

void histogramBuilder::fillCountHistogram(std::string key){
  if(!_h1Counter.count(key)){
    _h1Counter[key] = _fileService->make<TH1F>(Form("%s_Count",key.c_str()), 
					       Form("%s Count",key.c_str()),
					       2, 0, 2);  
  }
  _h1Counter[key]->Fill(1);
}

/*
 *Generic Form to build Trigger Histograms 
 *Type of object (or filtering) specified by the key
 */

void histogramBuilder::fillTrigHistograms(bool trigDecision,std::string key){
  if(!_h1Trig.count(key)){
    _h1Trig[key] = _fileService->make<TH1F>(Form("%s_Trig",key.c_str()), 
					    Form("%s Trigger",key.c_str()),
					    2, 0, 2);  
  }
  _h1Trig[key]->Fill(trigDecision);
}

/*
 *Generic From to Build L1Muon Pt Histograms
 *Type of object (or filtering) specified by the key
 *has variable binning
*/

void histogramBuilder::fillL1MuonPtHistograms(float pt, std::string key){
  
  if(!_h1L1MuonPt.count(key)){
    float variableBinArray[] = {0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,12,14,16,18,20,25,30,35,40,45,50,60,70,80,100,120,140,180};
    
    _h1L1MuonPt[key] = _fileService->make<TH1F>(Form("%s_Pt",key.c_str()), 
						Form("%s Pt",key.c_str()),
						33,
						variableBinArray);
  }
  _h1L1MuonPt[key]->Fill(pt);
}


/*
 *Generic From to Build Energy Histograms
 *Type of object (or filtering) specified by the key
 */

void histogramBuilder::fillEnergyHistograms(float energy, std::string key){
  if(!_h1Energy.count(key)){
    _h1Energy[key] = _fileService->make<TH1F>(Form("%s_Energy",key.c_str()), 
					    Form("%s Energy",key.c_str()),
					    2100, -5.0, 100.0);  
  }
  _h1Energy[key]->Fill(energy);
}


/*
 *Generic Form to build Eta Phi Histograms 
 *Type of object (or filtering) specified by the key
 */

void histogramBuilder::fillEtaPhiHistograms(float eta, float phi, std::string key){
  if(!_h1Eta.count(key)){
    _h1Eta[key] = _fileService->make<TH1F>(Form("Eta_%s",key.c_str()), 
					   Form("Eta %s",key.c_str()),
					   500, -1.5, 1.5);  //HO has 72 iphis and 30 ietas   
  }
  _h1Eta[key]->Fill(eta);

  if(!_h1Phi.count(key)){
    _h1Phi[key] = _fileService->make<TH1F>(Form("Phi_%s",key.c_str()), 
					   Form("Phi %s",key.c_str()),
					   500, -3.14, 3.14);  //HO has 72 iphis and 30 ietas
  }
  _h1Phi[key]->Fill(phi);

  if(!_h2EtaPhiMap.count(key)){
    _h2EtaPhiMap[key] = _fileService->make<TH2F>(Form("EtaPhi_%s",key.c_str()), 
						 Form("EtaPhi %s",key.c_str()),
						 500, -1.5, 1.5, 500, -3.14, 3.14);  
  }
  _h2EtaPhiMap[key]->Fill(eta, phi);
  
}


//define this as a plug-in
DEFINE_FWK_MODULE(histogramBuilder);
