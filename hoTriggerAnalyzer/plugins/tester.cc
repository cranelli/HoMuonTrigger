// -*- C++ -*-
//
// Package:    tester
// Class:      tester
// 
/**\class tester tester.cc L1TriggerDPGUpgrade/tester/plugins/tester.cc

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


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"

#include <vector>
#include <iostream>

#include "HoMuonTrigger/hoTriggerAnalyzer/interface/HistogramBuilder.h"


unsigned int NMaxL1AlgoBit = 128;

using namespace std;

//
// class declaration
//

class tester : public edm::EDAnalyzer {
public:
  explicit tester(const edm::ParameterSet&);
  ~tester();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;

  edm::Service<TFileService> _fileService
;

  edm::InputTag _genInput;
  edm::InputTag _l1GlobalTriggerReadoutInput;

  // TTree *ho_muon_tree;
  
  L1GtUtils m_l1GtUtils;

  HistogramBuilder histogramBuilder;
  
  //L1 Readout Record are put in a vector
  std::vector<bool> *vec_l1GlobalTrigger_Decision;
  
  //Generated particles are placed in a vector.
  /* 
  std::vector<int> *vec_generator_pdgId;
  std::vector<Float_t> *vec_generator_etas;
  std::vector<Float_t> * vec_generator_phis;
  std::vector<Float_t> * vec_generator_pts;
  */
      
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
tester::tester(const edm::ParameterSet& iConfig){
  
     //now do what ever initialization is needed
  
  
  _genInput = iConfig.getParameter<edm::InputTag>("genSrc");
  //_l1GlobalTriggerReadoutInput = iConfig.getParameter<edm::InputTag>("l1GlobalTriggerReadoutSrc");
  
  //ho_muon_tree = _fileService->make<TTree>("ho_muon_tree", 
  //				   "Generator, Propagator,and Reco Data");

    /*
     * Branches
     */
  /*
    ho_muon_tree->Branch("Generator_pdgId", "std::vector<int>",&vec_generator_pdgId);
    ho_muon_tree->Branch("Generator_Etas", "std::vector<Float_t>",&vec_generator_etas);
    ho_muon_tree->Branch("Generator_Phis", "std::vector<Float_t>",&vec_generator_phis);
    ho_muon_tree->Branch("Generator_Pts", "std::vector<Float_t>",&vec_generator_pts);
  */
  
    /*
     * Vectors
     */
  /*
    vec_l1GlobalTrigger_Decision =0;
    
    vec_generator_pdgId =0;
    vec_generator_etas =0;
    vec_generator_phis =0;
    vec_generator_pts =0;
  */
}


tester::~tester()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)


  
}


//
// member functions
//

// ------------ method called for each event  ------------
void
tester::analyze(const edm::Event& iEvent, 
		       const edm::EventSetup& iSetup)
{
   using namespace edm;



   //#ifdef THIS_IS_AN_EVENT_EXAMPLE
   //Handle<ExampleData> pIn;
   
   Handle<reco::GenParticleCollection> truthParticles;
   iEvent.getByLabel(_genInput,truthParticles);

   auto btruth = truthParticles->cbegin();
   auto etruth = truthParticles->cend();

   std::string test_gen_key = "test_gen";

   for(; btruth != etruth; ++btruth){
     //cout << btruth->eta() << endl;
     histogramBuilder.fillEtaPhiHistograms(btruth->eta(),btruth->phi(), test_gen_key);
   }

   /*
   delete vec_generator_pdgId; vec_generator_pdgId = new std::vector<int>();
   delete vec_generator_etas; vec_generator_etas = new std::vector<Float_t>();
   delete vec_generator_phis; vec_generator_phis = new std::vector<Float_t>();
   delete vec_generator_pts; vec_generator_pts = new std::vector<Float_t>();
  
   auto btruth = truthParticles->cbegin();
   auto etruth = truthParticles->cend();

   for( ; btruth != etruth; ++btruth ) {
     vec_generator_pdgId->push_back(btruth->pdgId());
     vec_generator_etas->push_back(btruth->eta());
     vec_generator_phis->push_back(btruth->phi());
     vec_generator_pts->push_back(btruth->pt());
   }
   */

   //Handle<L1GlobalTriggerReadoutRecord> l1GlobalTriggerReadouts;
   //iEvent.getByLabel(_l1GlobalTriggerReadoutInput,l1GlobalTriggerReadouts);

   /*
   for(unsigned int i=0; i< vec_generator_pdgId->size(); i++){
     std::cout << "Contains Particle Id: " << vec_generator_pdgId->at(i) << std::endl;
   }
   */
   
   //ho_muon_tree->Fill();
   
   //#endif

     /*   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
     */
}


// ------------ method called once each job just before starting event loop  ------------
void 
tester::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
tester::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------

void 
tester::beginRun(const edm::Run & iRun, const edm::EventSetup& evSetup)
{

  //cout << histogramBuilder.getMyNumber() << endl;

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
  
}


// ------------ method called when ending the processing of a run  ------------

void 
tester::endRun(const edm::Run& iRun, const edm::EventSetup& evSetup)
{

}


// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
tester::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
tester::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
tester::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(tester);
