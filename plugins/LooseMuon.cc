// -*- C++ -*-
//
// Package:    temp/LooseMuon
// Class:      LooseMuon
// 
/**\class LooseMuon LooseMuon.cc temp/LooseMuon/plugins/LooseMuon.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mengyao Shi
//         Created:  Wed, 25 Nov 2015 16:25:51 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//
//
// class declaration
//

class LooseMuon : public edm::EDFilter {
   public:
      explicit LooseMuon(const edm::ParameterSet&);
      ~LooseMuon();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
 edm::EDGetTokenT<edm::RefVector<std::vector<reco::Muon> > > muonTag_;
 unsigned int minNumObjsToPassFilter_;
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
LooseMuon::LooseMuon(const edm::ParameterSet& iConfig):
 muonTag_(consumes<edm::RefVector<std::vector<reco::Muon> > >(iConfig.getParameter<edm::InputTag>("muonTag"))),
 minNumObjsToPassFilter_(iConfig.getParameter<unsigned int>("minNumObjsToPassFilter"))
{

   //now do what ever initialization is needed
   produces<reco::MuonRefVector>();
}


LooseMuon::~LooseMuon()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
LooseMuon::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  bool isLooseMuon1=0;
  unsigned int nPassingMuons=0;
  edm::Handle<edm::RefVector<std::vector<reco::Muon> > > recoObjs;
  iEvent.getByToken(muonTag_, recoObjs);
  std::auto_ptr<reco::MuonRefVector> muonColl(new reco::MuonRefVector);
  for (typename edm::RefVector<std::vector<reco::Muon> >::const_iterator iRecoObj =
       recoObjs->begin(); iRecoObj != recoObjs->end();
       ++iRecoObj) 
  {
    isLooseMuon1=muon::isLooseMuon(**iRecoObj);
    muonColl->push_back(*iRecoObj);
    if(isLooseMuon1)
      nPassingMuons++;
  }
  
  iEvent.put(muonColl);

  return (nPassingMuons >= minNumObjsToPassFilter_);
}
   

// ------------ method called once each job just before starting event loop  ------------
void 
LooseMuon::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
LooseMuon::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
LooseMuon::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
LooseMuon::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
LooseMuon::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
LooseMuon::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
LooseMuon::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(LooseMuon);
