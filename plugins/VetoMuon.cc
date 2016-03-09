// -*- C++ -*-
//
// Package:    temp/VetoMuon
// Class:      VetoMuon
// 
/**\class VetoMuon VetoMuon.cc temp/VetoMuon/plugins/VetoMuon.cc

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

class VetoMuon : public edm::EDFilter {
   public:
      explicit VetoMuon(const edm::ParameterSet&);
      ~VetoMuon();

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
 edm::EDGetTokenT<reco::MuonRefVector> muonTag_;
 edm::EDGetTokenT<reco::MuonRefVector> vetoMuonTag_; 
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
VetoMuon::VetoMuon(const edm::ParameterSet& iConfig):
  muonTag_(consumes<reco::MuonRefVector>(iConfig.getParameter<edm::InputTag>("muonTag"))),
 vetoMuonTag_(consumes<reco::MuonRefVector>(iConfig.getParameter<edm::InputTag>("vetoMuonTag"))),
 minNumObjsToPassFilter_(iConfig.getParameter<unsigned int>("minNumObjsToPassFilter"))
{

   //now do what ever initialization is needed
   produces<reco::MuonRefVector>();
}


VetoMuon::~VetoMuon()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
VetoMuon::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  std::auto_ptr<reco::MuonRefVector> muonColl(new reco::MuonRefVector);

  edm::Handle<reco::MuonRefVector> pMuons;
  iEvent.getByToken(muonTag_, pMuons);

  edm::Handle<reco::MuonRefVector> pVetoMuons;
  iEvent.getByToken(vetoMuonTag_, pVetoMuons);


  std::vector<int> vetoMuonRefKeys;
  if (pVetoMuons.isValid()) {
    for (reco::MuonRefVector::const_iterator iVetoMuon = pVetoMuons->begin(); 
	 iVetoMuon != pVetoMuons->end(); ++iVetoMuon) {
      vetoMuonRefKeys.push_back(iVetoMuon->key());
    }
  }

  unsigned int nPassingMuons = 0;
  for(reco::MuonRefVector::const_iterator iMuon=pMuons->begin();
         iMuon!=pMuons->end();++iMuon){

    if (std::find(vetoMuonRefKeys.begin(), vetoMuonRefKeys.end(), 
		  iMuon->key()) == vetoMuonRefKeys.end()) {
      muonColl->push_back(*iMuon);
      ++nPassingMuons;

    }
  }
  iEvent.put(muonColl);

  return (nPassingMuons >= minNumObjsToPassFilter_);
}
   

// ------------ method called once each job just before starting event loop  ------------
void 
VetoMuon::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
VetoMuon::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
VetoMuon::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
VetoMuon::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
VetoMuon::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
VetoMuon::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
VetoMuon::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(VetoMuon);
