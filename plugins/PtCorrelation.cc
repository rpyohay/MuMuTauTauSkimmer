// -*- C++ -*-
//
// Package:    trigger_match/TriggerObjectFilter
// Class:      TriggerObjectFilter
// 
/**\class TriggerObjectFilter TriggerObjectFilter.cc trigger_match/TriggerObjectFilter/plugins/TriggerObjectFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mengyao Shi
//         Created:  Wed, 18 Nov 2015 13:18:21 GMT
//
//


// system include files
#include <memory>
#include <cmath>
// user include files
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TGraph.h"
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
//#include "BoostedTauAnalysis/Common/interface/Common.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "HLTrigger/HLTanalyzers/interface/HLTInfo.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTConfigData.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
using namespace std;
using namespace edm;
using namespace reco;
using namespace trigger;

//
// class declaration
//
template<class T>
class PtCorrelation : public edm::EDFilter {
   public:
      explicit PtCorrelation(const edm::ParameterSet&);
      ~PtCorrelation();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const& , edm::EventSetup const& )/* override*/;
      virtual void beginRun(const edm::Run& iRun, edm::EventSetup const& iSetup);
      virtual void beginJob();
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob();
      virtual void endRun(const edm::Run& iRun, edm::EventSetup const& iSetup);
//      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)/* override*/;
      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::RefVector<std::vector<T> > > recoObjTag17_;
      edm::EDGetTokenT<edm::RefVector<std::vector<T> > > recoObjTag8_;
      std::map<std::string, TH2D*> histos2D_;
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
PtCorrelation<T>::PtCorrelation(const edm::ParameterSet& iConfig):
  recoObjTag17_(consumes<edm::RefVector<std::vector<T> > >(iConfig.getParameter<edm::InputTag>("recoObjTag17"))),
  recoObjTag8_(consumes<edm::RefVector<std::vector<T> > >(iConfig.getParameter<edm::InputTag>("recoObjTag8"))),
  histos2D_()
{
}

template<class T>
PtCorrelation<T>::~PtCorrelation()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}
//template<class T>
//void PtCorrelation<T>::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
//{
//  std::cout<< "beginLuminosityBlock"<< std::endl;
//  return true;
//}


//
// member functions
//
template<class T>
void PtCorrelation<T>::beginRun(const edm::Run& iRun, edm::EventSetup const& iSetup)
{
}

template<class T>
void
PtCorrelation<T>::beginJob()
{
  //std::cout<< "beginJob" << std::endl;
  edm::Service<TFileService> fileService;
  histos2D_["ptTrigCand1"] =fileService->make< TH2D >("ptTrigCand1","Object 17 vs. Obj 8",150, 0., 150., 150, 0., 150.);
  histos2D_[ "ptTrigCand1" ]->SetXTitle( "obj17 p_{T} (GeV)" );
  histos2D_[ "ptTrigCand1" ]->SetYTitle( "obj8 p_{T} (GeV)" );


}

// ------------ method called on each new Event  ------------
template<class T>
bool
PtCorrelation<T>::filter( edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //get reco objects
  edm::Handle<edm::RefVector<std::vector<T> > > recoObjs17;
  iEvent.getByToken(recoObjTag17_, recoObjs17);
  edm::Handle<edm::RefVector<std::vector<T> > > recoObjs8;
  iEvent.getByToken(recoObjTag8_, recoObjs8);
  bool pass=0;
  if(recoObjs17->size() && recoObjs8->size()){
    pass=1;
    for (typename edm::RefVector<std::vector<T> >::const_iterator iRecoObj17 =
                recoObjs17->begin(); iRecoObj17 != recoObjs17->end();
              ++iRecoObj17){
      for (typename edm::RefVector<std::vector<T> >::const_iterator iRecoObj8 =
                recoObjs8->begin(); iRecoObj8 != recoObjs8->end();
              ++iRecoObj8){
         histos2D_[ "ptTrigCand1"]->Fill((*iRecoObj17)->pt(), (*iRecoObj8)->pt());      
 
      }
    }
  }
  return (pass);
}

// ------------ method called when starting to processes a run  ------------


// ------------ method called once each job just before starting event loop  ------------

// ------------ method called once each job just after ending the event loop  ------------
template<class T>
void 
PtCorrelation<T>::endJob() {
}


// ------------ method called when ending the processing of a run  ------------
template<class T>
void PtCorrelation<T>::endRun(const edm::Run& iRun, edm::EventSetup const& iSetup)
{
 // std::cout<< "endRun" << std::endl;
} 
// ------------ method called when starting to processes a luminosity block  ------------
 
// ------------ method called when ending the processing of a luminosity block  ------------
//template<class T>
//bool TriggerObjectFilter<T>::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
//{
//  std::cout<< "endLuminosityBlock"<< std::endl;
//  return true;
//}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
template<class T>
void
PtCorrelation<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
typedef PtCorrelation<reco::Muon> MuonPtCorrelation;
DEFINE_FWK_MODULE(MuonPtCorrelation);
