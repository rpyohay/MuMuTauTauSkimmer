// -*- C++ -*-
//
// Package:    GenMatchedRecoObjectProducer1
// Class:      GenMatchedRecoObjectProducer1
// 
/**\class GenMatchedRecoObjectProducer1 GenMatchedRecoObjectProducer1.cc 
   BoostedTauAnalysis/GenMatchedRecoObjectProducer1/src/GenMatchedRecoObjectProducer1.cc

Description: produce a collection of reco objects matched to gen boosted di-tau objects

Implementation:

*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Thu Aug 23 11:23:58 CEST 2012
// $Id: GenMatchedRecoObjectProducer1.cc,v 1.3 2012/09/25 11:44:34 yohay Exp $
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
#include "Tools/Common/interface/GenTauDecayID.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "Tools/Common/interface/Common.h"
#include "DataFormats/MuonReco/interface/Muon.h"

//
// class declaration
//

template<class T>
class GenMatchedRecoObjectProducer1 : public edm::EDFilter {
public:
  explicit GenMatchedRecoObjectProducer1(const edm::ParameterSet&);
  ~GenMatchedRecoObjectProducer1();

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


  /*input tag for gen particle collection to match
    count on the user to pass in a collection that will not lead to the same reco object being 
    matched to multiple different gen objects
    for example, if the input object is a boosted di-tau pair, only 1 member of the pair should be 
    in the input collection*/
  edm::InputTag selectedGenParticleTag_;

  //input tag for reco object collection
  edm::InputTag recoObjTag_;

  //input tag for base reco object collection
  edm::InputTag baseRecoObjTag_;
  
  //dR matching cut
  double dR_;

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
template<class T>
GenMatchedRecoObjectProducer1<T>::GenMatchedRecoObjectProducer1(const edm::ParameterSet& iConfig) :
  selectedGenParticleTag_(iConfig.getParameter<edm::InputTag>("selectedGenParticleTag")),
  recoObjTag_(iConfig.getParameter<edm::InputTag>("recoObjTag")),
  baseRecoObjTag_(iConfig.getParameter<edm::InputTag>("baseRecoObjTag")),
  dR_(iConfig.getParameter<double>("dR")),
  minNumGenObjectsToPassFilter_
  (iConfig.getParameter<unsigned int>("minNumGenObjectsToPassFilter"))
{
  //register your products
  produces<edm::RefVector<std::vector<T> > >();

  //now do what ever other initialization is needed
  
}


template<class T>
GenMatchedRecoObjectProducer1<T>::~GenMatchedRecoObjectProducer1()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
template<class T>
bool GenMatchedRecoObjectProducer1<T>::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //get selected gen particles
  edm::Handle<reco::GenParticleRefVector> pSelectedGenParticles;
  iEvent.getByLabel(selectedGenParticleTag_, pSelectedGenParticles);

  //get reco object collection
  edm::Handle<edm::RefVector<std::vector<T> > > pRecoObjs;
  iEvent.getByLabel(recoObjTag_, pRecoObjs);

  //get base reco object collection
  edm::Handle<std::vector<T> > pBaseRecoObjs;
  iEvent.getByLabel(baseRecoObjTag_, pBaseRecoObjs);

  //fill STL container of pointers to reco objects
  std::vector<T*> recoObjPtrs;
  for (typename edm::RefVector<std::vector<T> >::const_iterator iRecoObj = pRecoObjs->begin(); 
       iRecoObj != pRecoObjs->end(); ++iRecoObj) {
    recoObjPtrs.push_back(const_cast<T*>(iRecoObj->get()));
  }

  //make a copy of the reco object vector
  std::vector<T*> recoObjPtrsCopy = recoObjPtrs;


  //declare pointers to output collection to produce
  std::auto_ptr<edm::RefVector<std::vector<T> > >  genMatchedRecoObjs(new edm::RefVector<std::vector<T>>);


  //loop over selected gen particles
  for (reco::GenParticleRefVector::iterator iGenObj = pSelectedGenParticles->begin(); 
       iGenObj !=pSelectedGenParticles->end(); ++iGenObj) {
    int nearestRecoObjPTRank = 0; /*this is the index into recoObjPtrsCopy of the nearest object
				     since recoObjPtrsCopy is sorted in ascending order by pT, the 
				     pT rank of the nearest object is recoObjPtrsCopy.size() - 
				     nearestRecoObjPTRank - 1*/
    const T* nearestRecoObj = 
      Common::nearestObject(*iGenObj, recoObjPtrsCopy, nearestRecoObjPTRank);

    if ((nearestRecoObj != NULL) && (nearestRecoObjPTRank >= 0) && 
	(reco::deltaR(*nearestRecoObj, **iGenObj) < dR_)) {


	nearestRecoObjPTRank = 0; /*since only 1 output collection is produced in this case, the 
				    index into genMatchedRecoObjs should be 0*/
	genMatchedRecoObjs->
	  push_back(edm::Ref<std::vector<T> >(pBaseRecoObjs, 
					      pRecoObjs->at(nearestRecoObjPTRank).key()));
      
    }
  }

  //flag indicating whether right number of gen-matched reco objects were found
  bool foundGenMatchedRecoObject = genMatchedRecoObjs->size() >= minNumGenObjectsToPassFilter_;

  //put output collection into event
  iEvent.put(genMatchedRecoObjs);

  //stop processing if no gen-matched objects were found
  return foundGenMatchedRecoObject;
}

// ------------ method called once each job just before starting event loop  ------------
template<class T>
void GenMatchedRecoObjectProducer1<T>::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
template<class T>
void GenMatchedRecoObjectProducer1<T>::endJob() {
}

// ------------ method called when starting to processes a run  ------------
template<class T>
bool GenMatchedRecoObjectProducer1<T>::beginRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a run  ------------
template<class T>
bool GenMatchedRecoObjectProducer1<T>::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
template<class T>
bool GenMatchedRecoObjectProducer1<T>::beginLuminosityBlock(edm::LuminosityBlock&, 
							   edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
template<class T>
bool GenMatchedRecoObjectProducer1<T>::endLuminosityBlock(edm::LuminosityBlock&, 
							 edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  -----------
template<class T>
void 
GenMatchedRecoObjectProducer1<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
typedef GenMatchedRecoObjectProducer1<reco::Muon> GenMatchedMuonProducer1;
DEFINE_FWK_MODULE(GenMatchedMuonProducer1);
