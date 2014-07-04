#ifndef __HOMUON_HOMUONTREEBUILDER_H__
#define __HOMUON_HOMUONTREEBUILDER_H__
//
// Class:hoMuonTreeBuilder
// 
/*
 Description: Currently the header file for the entire HOMuon Trigger Analysis

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Christopher Anelli  

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

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"
#include "DataFormats/HcalRecHit/interface/HORecHit.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"

#include "HoMuonTrigger/hoTriggerAnalyzer/interface/HistogramBuilder.h"

#include <TTree.h>
#include <vector>
#include <iostream>
#include <map>
#include <list>

#include "math.h"

using namespace::std;


//Custom Structures

struct basic_Generator{
  vector<int>* pdgIds=0;
  vector<float>* etas=0;
  vector<float>* phis=0;
  vector<float>* pts=0;
};

struct basic_TrigObject{
  vector<float>* etas=0;
  vector<float>* phis=0;
  vector<float>* pts=0;
};

struct basic_HORecHit{
  vector<float>* etas=0;
  vector<float>* phis=0;
  vector<float>* energies=0;
};


//
// HO Muon Tree Builder Class Declaration
//

class hoMuonTreeBuilder : public edm::EDAnalyzer {
 public:
  explicit hoMuonTreeBuilder(const edm::ParameterSet&);
  ~hoMuonTreeBuilder();
  
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

  edm::InputTag _genInfoInput;
  edm::InputTag _genInput;
  edm::InputTag _l1MuonInput;
  //edm::InputTag _stdMuInput;
  edm::InputTag _horecoInput;
  //edm::InputTag _l1GtTmLInputTag;
  edm::InputTag _hltSumAODInput;
  HistogramBuilder histogramBuilder;
  
  //Utility for the Level 1 Information
  L1GtUtils m_l1GtUtils;
  //string m_nameAlgTechTrig;
  //bool trigDecision;

  //toFigureOutL1MuonVariableBinning
  //list<float> listL1MuonPt;
  
  //vector<l1extra::L1MuonParticle> l1MuonsBarrelWithMip;
  
  /*
   * Maps of selected hlt triggers to get the trigger decisions,
   * and hlt filters to get the trigger objects.
   */
  
  map<string, string> hltNamesOfInterest;
  map<string, edm::InputTag> hltFiltersOfInterest;
  void defineTriggersOfInterest();
  //vector<trigger::TriggerObject> hltTriggerObjects;
  //map<string, vector<trigger::TriggerObject>> hltTriggerObjectsOfInterest;
  
  
  //For Filtering
  // bool hoBelowThreshold(HORecHit horeco);
  
  
  //Tree Variables
  TTree *ho_muon_tree;
  void initializeBranches();
  void fillGeneratorWeights(edm::Handle<GenEventInfoProduct> & genEventInfo);
  void fillGeneratorParticles(edm::Handle<reco::GenParticleCollection> & genParticles);
  void fillL1Muons(edm::Handle<l1extra::L1MuonParticleCollection> & l1Muons);
  void fillHORecHits(edm::Handle<HORecHitCollection> &hoRecoHits,
		     edm::ESHandle<CaloGeometry> & caloGeo);
  void fillHLTTrigObjects(edm::Handle<trigger::TriggerEvent> triggerSummary,
			  string hlt_key, edm::InputTag lastFilterTag);
  void fillGenMuonPropParticles(edm::Handle<reco::GenParticleCollection> & genParticles,
				edm::ESHandle<MagneticField> & bField,
				edm::ESHandle<Propagator> & shProp);
  void fillGenMuonPropParticles(double radius, basic_Generator *genMuonProp, 
				edm::Handle<reco::GenParticleCollection> & genParticles,
				edm::ESHandle<MagneticField> & bField,
				edm::ESHandle<Propagator> & shProp);
  void fillHLTPropObjects(edm::Handle<trigger::TriggerEvent> & triggerSummary,
			  edm::ESHandle<MagneticField> & bField,
			  edm::ESHandle<Propagator> & shProp);
  void fillHLTPropObjects(double radius, basic_TrigObject *hltPropToRPC1, 
			   edm::Handle<trigger::TriggerEvent> & triggerSummary,
			   edm::ESHandle<MagneticField> & bField,
			   edm::ESHandle<Propagator> & shProp);
   TrajectoryStateOnSurface propagateToCylinder(FreeTrajectoryState & initial, double radius, 
						edm::ESHandle<Propagator> & shProp);
    
  //Particle Information
  double weight;
  basic_Generator generator;
  basic_TrigObject l1muon;
  basic_HORecHit horeco;
  basic_Generator genMuonPropToHO;
  basic_Generator genMuonPropToRPC1;
  basic_TrigObject hltPropToRPC1;
  map<string, basic_TrigObject> mapHLTObjects;

  // HLT Trigger Object Information (Possible to look at multiple filters)
  //map<string, vector<float>*> hltTrigObjects_etas;
  // map<string,vector<float>*> hltObjects_phis;
  //map<string,vector<float>*> hltObjects_pts;


  //Radiuses for Propagation
  
  //double ho_radius = 412.6; // In cm                                                  
  //From hcalouteralgo.xml Chose yLayer 1, same in all rings.                        
  
  double yLayer1 = 406.6;
  double xOff = -36.0674;
  double ho_radius = sqrt(pow(yLayer1,2)+pow(xOff,2));
  
  double rpc1_radius = 414.65;

  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
};

#endif
