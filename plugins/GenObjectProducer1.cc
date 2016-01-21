// -*- C++ -*-
//
// Package:    GenObjectProducer1
// Class:      GenObjectProducer1
// 
/**\class GenObjectProducer1 GenObjectProducer1.cc 
   BoostedTauAnalysis/GenMatchedRecoObjectProducer/src/GenObjectProducer1.cc

Description: produce a collection of gen objects from boosted di-tau objects

Implementation:

*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Thu Aug 23 11:23:58 CEST 2012
// $Id: GenObjectProducer1.cc,v 1.4 2012/10/10 09:11:43 yohay Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "Tools/Common/interface/GenTauDecayID.h"
#include "Tools/Common/interface/Common.h"

//code for any tau decay
#define TAU_ALL 3

//
// class declaration
//

class GenObjectProducer1 : public edm::EDFilter {
public:
  explicit GenObjectProducer1(const edm::ParameterSet&);
  ~GenObjectProducer1();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
      
  virtual bool beginRun(edm::Run&, edm::EventSetup const&);
  virtual bool endRun(edm::Run&, edm::EventSetup const&);
  virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);




  // ----------member data ---------------------------

  //input tag for gen particle collection
  edm::InputTag genParticleTag_;

  //minimum number of gen objects passing cuts that must be found for event to pass filter
  unsigned int minNumGenObjectsToPassFilter_;


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
GenObjectProducer1::GenObjectProducer1(const edm::ParameterSet& iConfig) :
  genParticleTag_(iConfig.getParameter<edm::InputTag>("genParticleTag")),
  minNumGenObjectsToPassFilter_
  (iConfig.getParameter<unsigned int>("minNumGenObjectsToPassFilter"))
{

  //now do what ever other initialization is needed
  produces<reco::GenParticleRefVector>();  
}


GenObjectProducer1::~GenObjectProducer1()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
bool GenObjectProducer1::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //get GEN particles
  edm::Handle<reco::GenParticleCollection> pGenParticles;
  iEvent.getByLabel(genParticleTag_, pGenParticles);

  //fill STL container of all decay products of pseudoscalar Higgses
  //std::vector<reco::GenParticleRef> aDecayProducts;
  std::auto_ptr<reco::GenParticleRefVector> aDecayProducts(new reco::GenParticleRefVector);
  for (reco::GenParticleCollection::const_iterator iGenParticle = pGenParticles->begin(); 
       iGenParticle != pGenParticles->end(); ++iGenParticle) {
    if((*iGenParticle).pdgId() == 35 && ((*iGenParticle).numberOfDaughters()==2) )
    {
      reco::GenParticleRef child0 = iGenParticle->daughterRef(0);
      reco::GenParticleRef child1 = iGenParticle->daughterRef(1);
      if(child0->pdgId() == 36 && child1->pdgId()==36 && (child0->numberOfDaughters()==2)&&(child1->numberOfDaughters()==2) )
      { 
        
        reco::GenParticleRef childchild00=child0->daughterRef(0);
        reco::GenParticleRef childchild01=child0->daughterRef(1);
        reco::GenParticleRef childchild10=child1->daughterRef(0);
        reco::GenParticleRef childchild11=child1->daughterRef(1);
        if(fabs(childchild00->pdgId())==13 &&( fabs(childchild10->pdgId())==15))
        {
          if(childchild00->pt() >childchild01->pt())
            aDecayProducts->push_back(childchild00);
          else
            aDecayProducts->push_back(childchild01);
        } 
        else if(fabs(childchild00->pdgId())==15 &&(fabs(childchild10->pdgId())==13))
        {
          if(childchild10->pt() > childchild11->pt())
          aDecayProducts->push_back(childchild10);
          else
          aDecayProducts->push_back(childchild11);
        }
      }
    }
  }

//flag indicating whether enough gen-matched reco objects were found
  bool foundGenObject = aDecayProducts->size() >= minNumGenObjectsToPassFilter_;


  //set the pT rank of the a decay product (highest rank is 0, next highest is 1, etc.)
  iEvent.put(aDecayProducts);
  return foundGenObject;
}

// ------------ method called once each job just before starting event loop  ------------
void GenObjectProducer1::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void GenObjectProducer1::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool GenObjectProducer1::beginRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool GenObjectProducer1::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool GenObjectProducer1::beginLuminosityBlock(edm::LuminosityBlock&, 
							   edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool GenObjectProducer1::endLuminosityBlock(edm::LuminosityBlock&, 
							 edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  -----------
void GenObjectProducer1::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(GenObjectProducer1);
