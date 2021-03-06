// -*- C++ -*-
//
// Class:      hoMuonTreeBuilder
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

//hoMuonTreeBuilder header file
#include "Analysis/hoTriggerAnalyzer/interface/hoMuonTreeBuilder.h"

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

//bool hoBelowThreshold(HORecHit horeco);
//bool isInsideRCut(float deltaR_Max, float eta1, float eta2, float phi1, float phi2);


hoMuonTreeBuilder::hoMuonTreeBuilder(const edm::ParameterSet& iConfig){
  
  //now do what ever initialization is needed
 


  //Get Input Tags from the Configuration
  
  _genInput = iConfig.getParameter<edm::InputTag>("genSrc");
  _l1MuonInput = iConfig.getParameter<edm::InputTag>("l1MuonSrc");
  _horecoInput = iConfig.getParameter<edm::InputTag>("horecoSrc");
  //_stdMuInput = iConfig.getParameter<edm::InputTag>("stdMuSrc");
  _hltSumAODInput = iConfig.getParameter<edm::InputTag>("hltSumAODSrc");
  //_hltTriggerResults = iConfig.getParameter<edm::InputTag>("hltResultsSrc");
  // m_l1GtTmLInputTag = iConfig.getParameter<edm::InputTag> ("L1GtTmLInputTag");
  
  //m_nameAlgTechTrig="L1_SingleMu7";


  defineTriggersOfInterest();

  
  
  // Tree
  
  ho_muon_tree = _fileService->make<TTree>("ho_muon_tree", "Generator, L1, L1Muon," 
					   "HLT, HLT Trigger Object, and HO Reco Data");

  initializeBranches();
  
}


hoMuonTreeBuilder::~hoMuonTreeBuilder()
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
hoMuonTreeBuilder::analyze(const edm::Event& iEvent, 
		       const edm::EventSetup& iSetup)
{
   using namespace edm;
   
   /*
    * Get Event Data and Event Setup
    */

   Handle<reco::GenParticleCollection> genParticles;
   iEvent.getByLabel(_genInput,genParticles);

   Handle<l1extra::L1MuonParticleCollection> l1Muons;
   iEvent.getByLabel(_l1MuonInput, l1Muons);

   //Handle<reco::TrackCollection> standAloneMuons;
   //iEvent.getByLabel(_stdMuInput,standAloneMuons);

   Handle<HORecHitCollection> hoRecoHits;
   iEvent.getByLabel(_horecoInput, hoRecoHits);

   Handle<trigger::TriggerEvent> triggerSummary;
   iEvent.getByLabel(_hltSumAODInput, triggerSummary);

   InputTag _hltTrigResultsLabel = edm::InputTag("TriggerResults", "", "HLT");
   Handle<TriggerResults> hltTriggerResults;
   iEvent.getByLabel(_hltTrigResultsLabel, hltTriggerResults);
   //const TriggerNames & names = iEvent.triggerNames(*hltTriggerResults);
   
   ESHandle<CaloGeometry> caloGeo;
   iSetup.get<CaloGeometryRecord>().get(caloGeo);


   /*
    * Set Up Level 1 Global Trigger Utility
    */
         
   bool useL1EventSetup = true;
   bool useL1GtTriggerMenuLite = true;

   m_l1GtUtils.getL1GtRunCache(iEvent, iSetup, useL1EventSetup, 
			       useL1GtTriggerMenuLite);
   


   fillGeneratorParticles(genParticles);
   fillL1Muons(l1Muons);
   fillHORecHits(hoRecoHits, caloGeo);
   
   auto bHLTFilters = hltFiltersOfInterest.begin();
   for(; bHLTFilters != hltFiltersOfInterest.end(); ++bHLTFilters){
     fillHLTTrigObjects(triggerSummary, bHLTFilters->first, bHLTFilters->second);
   }
   
   ho_muon_tree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
hoMuonTreeBuilder::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
hoMuonTreeBuilder::endJob() 
{
  //h1SAMuonPt->Write();
  //histogramBuilder.fillCountHistogram(eventCounter,"Events");
  //h1EventCounter->Write();
}

// ------------ method called when starting to processes a run  ------------

void 
hoMuonTreeBuilder::beginRun(const edm::Run& iRun, 
			   const edm::EventSetup& evSetup)
{
  
  bool useL1EventSetup = true;
  bool useL1GtTriggerMenuLite = true;
  cout << "getL1GtRunCache" << endl;
  cout << "UseL1EventSetup: " << useL1EventSetup << "UseL1GtTriggerMenuLite :" 
       << useL1GtTriggerMenuLite << endl;  
  m_l1GtUtils.getL1GtRunCache(iRun, evSetup, useL1EventSetup, useL1GtTriggerMenuLite);

  
}


// ------------ method called when ending the processing of a run  ------------

void 
hoMuonTreeBuilder::endRun(const edm::Run& iRun, const edm::EventSetup& evSetup)
{
  
  //Only interested in unique values
  /*
  listL1MuonPt.sort();
  listL1MuonPt.unique();  //NB it is called after sort
  cout <<"The list contains " << listL1MuonPt.size() << "unique entries:";
  std::list<float>::iterator it; 
  for (it=listL1MuonPt.begin(); it!=listL1MuonPt.end(); ++it){
    cout << ' ' << *it;
  }
  cout << endl;
  */
}


// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
hoMuonTreeBuilder::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
hoMuonTreeBuilder::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
hoMuonTreeBuilder::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

/*
 * Helper Functions for Filter
 */

/*
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
*/


void hoMuonTreeBuilder::defineTriggersOfInterest(){
  
  /*
   * HLT Triggers
   */
  
  basic_TrigObject hltIsoMu24;
  string hltIsoMu24_key = "hltIsoMu24";
  hltNamesOfInterest.insert(pair<string, string>(hltIsoMu24_key,"HLT_IsoMu24_v18"));
  hltFiltersOfInterest.insert(pair<string, edm::InputTag>(hltIsoMu24_key, 
							  edm::InputTag("hltL3crIsoL1sMu16L1f0L2f16QL3"
									"f24QL3crIsoRhoFiltered0p15",
									"","HLT")));

  mapHLTObjects.insert(pair<string,basic_TrigObject>(hltIsoMu24_key, hltIsoMu24));

  basic_TrigObject hltMu5;
  string hltMu5_key = "hltMu5";
  hltNamesOfInterest.insert(pair<string, string>(hltMu5_key, "HLT_Mu5_v21"));
  hltFiltersOfInterest.insert(pair<string, edm::InputTag>(hltMu5_key, 
							  edm::InputTag("hltL3fL1sMu3L3Filtered5",
									"","HLT")));
  mapHLTObjects.insert(pair<string,basic_TrigObject>(hltMu5_key, hltMu5));
}

void hoMuonTreeBuilder::fillGeneratorParticles(edm::Handle<reco::GenParticleCollection>
					       & genParticles){
  generator.pdgIds->clear();
  generator.etas->clear();
  generator.phis->clear();
  generator.pts->clear();
  
  auto bgen = genParticles->begin();
  
  for( ; bgen != genParticles->end(); ++bgen ) {
    generator.pdgIds->push_back(bgen->pdgId());
    generator.etas->push_back(bgen->eta());
    generator.phis->push_back(bgen->phi());
    generator.pts->push_back(bgen->pt());
  }
}

void hoMuonTreeBuilder::fillL1Muons(edm::Handle<l1extra::L1MuonParticleCollection> & l1Muons){
  l1muon.etas->clear();
  l1muon.phis->clear();
  l1muon.pts->clear();
  
  auto bL1muon = l1Muons->cbegin();
  
  for( ; bL1muon != l1Muons->cend(); ++bL1muon ) {
    l1muon.etas->push_back(bL1muon->eta());
    l1muon.phis->push_back(bL1muon->phi());
    l1muon.pts->push_back(bL1muon->pt());
  }
}

void hoMuonTreeBuilder::fillHORecHits(edm::Handle<HORecHitCollection> &hoRecoHits,
		   edm::ESHandle<CaloGeometry> & caloGeo){
  horeco.etas->clear();
  horeco.phis->clear();
  horeco.energies->clear();
  
  auto bho_reco = hoRecoHits->begin();
  
  for( ; bho_reco != hoRecoHits->end(); ++bho_reco ) {
    horeco.etas->push_back(caloGeo->getPosition(bho_reco->id()).eta());
    horeco.phis->push_back(caloGeo->getPosition(bho_reco->id()).phi());
    horeco.energies->push_back(bho_reco->energy());
  }
}

void hoMuonTreeBuilder::fillHLTTrigObjects(edm::Handle<trigger::TriggerEvent> triggerSummary,
					   string key, edm::InputTag lastFilterTag){
  mapHLTObjects[key].etas->clear();
  mapHLTObjects[key].phis->clear();
  mapHLTObjects[key].pts->clear();

  //Store the Trigger Objects, from the Filters Of Interest  
  //Define TriggerObject Collection                                                                  
                     
  trigger::TriggerObjectCollection hltObjectCollection = triggerSummary->getObjects();
  
  
  //find the filter index 
  size_t filterIndex = triggerSummary->filterIndex(lastFilterTag); 
  if(filterIndex < triggerSummary->sizeFilters()){  //Check that Filter is in Trigger Summary     
    //cout << "Trigger Object should be Present" << endl;                                          
    const trigger::Keys &filterKeys = triggerSummary->filterKeys(filterIndex);
    for(size_t j=0; j < filterKeys.size(); j++){
      //cout << "Size of Keys"<< keys.size()  << endl;                   
      //Get Trigger Object using Keys
      mapHLTObjects[key].etas->push_back(hltObjectCollection[filterKeys[j]].eta());
      mapHLTObjects[key].phis->push_back(hltObjectCollection[filterKeys[j]].phi());
      mapHLTObjects[key].pts->push_back(hltObjectCollection[filterKeys[j]].pt());
    }
  }

}


void hoMuonTreeBuilder::initializeBranches(){
  //Generator
  ho_muon_tree->Branch("Generator_pdgIds",
                       "std::vector<int>",&generator.pdgIds);
  ho_muon_tree->Branch("Generator_Etas",
                       "std::vector<Float_t>",&generator.etas);
  ho_muon_tree->Branch("Generator_Phis",
                       "std::vector<Float_t>",&generator.phis);
  ho_muon_tree->Branch("Generator_Pts",
                       "std::vector<Float_t>",&generator.pts);
  //L1Muon
  ho_muon_tree->Branch("L1Muon_Etas",
                       "std::vector<Float_t>",&l1muon.etas);
  ho_muon_tree->Branch("L1Muon_Phis",
                       "std::vector<Float_t>",&l1muon.phis);
  ho_muon_tree->Branch("L1Muon_Pts",
                       "std::vector<Float_t>",&l1muon.pts);
  //HOReco
  ho_muon_tree->Branch("HOReco_Etas",
                       "std::vector<Float_t>",&horeco.etas);
  ho_muon_tree->Branch("HOReco_Phis",
                       "std::vector<Float_t>",&horeco.phis);
  ho_muon_tree->Branch("HOReco_Pts",
                       "std::vector<Float_t>",&horeco.energies);

  // HLT Triggers Objects
  
  auto trigMapIt = mapHLTObjects.begin();
  for(; trigMapIt != mapHLTObjects.end(); ++trigMapIt){

    string hlt_key = trigMapIt->first;
    //basic_TrigObject hltTrigObject = trigMapIt->second; 
    std::stringstream etaBranchName, phiBranchName, ptBranchName;
 
    etaBranchName << hlt_key << "_Etas";
    ho_muon_tree->Branch(etaBranchName.str().c_str(),
			 "std::vector<Float_t>",&(trigMapIt->second.etas));

    phiBranchName << hlt_key << "_Phis";
    ho_muon_tree->Branch(phiBranchName.str().c_str(),
			 "std::vector<Float_t>",&(trigMapIt->second.phis));

    ptBranchName << hlt_key << "_Pts";
    ho_muon_tree->Branch(ptBranchName.str().c_str(),
			 "std::vector<Float_t>",&(trigMapIt->second.pts));
  }
    
}

//define this as a plug-in
DEFINE_FWK_MODULE(hoMuonTreeBuilder);
