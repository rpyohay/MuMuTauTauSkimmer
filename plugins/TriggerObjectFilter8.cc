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
class TriggerObjectFilter8 : public edm::EDFilter {
   public:
      explicit TriggerObjectFilter8(const edm::ParameterSet&);
      ~TriggerObjectFilter8();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
//      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)/* override*/;
      virtual void beginRun(const edm::Run& iRun, edm::EventSetup const& iSetup);
      virtual void beginJob();
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob();
      //virtual bool endRun(edm::Run&, edm::EventSetup const&);
//      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)/* override*/;
      // ----------member data ---------------------------
      edm::InputTag recoObjTag_;
      edm::InputTag triggerEventTag_;
      edm::InputTag triggerResultsTag_;
      double delRMatchingCut_;
      std::vector<edm::InputTag> hltTags_;
      HLTConfigProvider hltConfig_;
      edm::InputTag theRightHLTTag_;
      edm::InputTag theRightHLTSubFilter0_;
      std::vector<edm::InputTag> HLTSubFilters_;
      unsigned int minNumObjsToPassFilter0_;
      std::map<std::string, TH1D*> histos1D_;
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
TriggerObjectFilter8<T>::TriggerObjectFilter8(const edm::ParameterSet& iConfig):
hltConfig_(),
histos1D_(), 
histos2D_()
{
  //now do what ever initialization is needed
   
  recoObjTag_ = iConfig.getParameter<edm::InputTag>("recoObjTag");
  const edm::InputTag dTriggerEventTag("hltTriggerSummaryAOD","","HLT");
  triggerEventTag_ = iConfig.getUntrackedParameter<edm::InputTag>("triggerEventTag",dTriggerEventTag);
  const edm::InputTag dTriggerResults("TriggerResults","","HLT");
  // By default, trigger results are labeled "TriggerResults" with process name "HLT" in the event.
  triggerResultsTag_ = iConfig.getUntrackedParameter<edm::InputTag>("triggerResultsTag",dTriggerResults);
  delRMatchingCut_ = iConfig.getUntrackedParameter<double>("triggerDelRMatch", 0.30);
  hltTags_ = iConfig.getParameter<std::vector<edm::InputTag> >("hltTags");
  //hltConfig_ = iConfig.getParameter<HLTConfigProvider>("hltConfig");
  theRightHLTTag_ = iConfig.getParameter<edm::InputTag>("theRightHLTTag");
  theRightHLTSubFilter0_ = iConfig.getParameter<edm::InputTag>("theRightHLTSubFilter0");
  
  //Whether using HLT trigger path name or the actual trigger filter name. Trigger path is default.
  HLTSubFilters_ = iConfig.getUntrackedParameter<std::vector<edm::InputTag> >("HLTSubFilters",std::vector<edm::InputTag>());
  minNumObjsToPassFilter0_ = iConfig.getParameter<unsigned int>("minNumObjsToPassFilter0");

  produces<edm::RefVector<std::vector<T> > >();
}

template<class T>
TriggerObjectFilter8<T>::~TriggerObjectFilter8()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}
//template<class T>
//bool TriggerObjectFilter<T>::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
//{
//  std::cout<< "beginLuminosityBlock"<< std::endl;
//  return true;
//}


//
// member functions
//
template<class T>
void TriggerObjectFilter8<T>::beginRun(const edm::Run& iRun, edm::EventSetup const& iSetup)
{
  bool changed_ = true;
  //std::cout << "beginRun"<< std::endl;
  
  if ( !hltConfig_.init(iRun,iSetup,hltTags_[0].process(),changed_) ){
    edm::LogError("TriggerObjectFilter") <<
                     "Error! Can't initialize HLTConfigProvider";
                         throw cms::Exception("HLTConfigProvider::init() returned non 0");
  }
}

template<class T>
void
TriggerObjectFilter8<T>::beginJob()
{
  //std::cout<< "beginJob" << std::endl;
  edm::Service<TFileService> fileService;
  //TCanvas *c1=new TCanvas("c1","graph with trigger efficiency",200,10,700,500);
  //c1->SetFillColor(42);
  //c1->SetGrid();
  //c1->GetFrame()->SetFillColor(21);
  //c1->GetFrame()->SetBorderSize(12);

  histos1D_["keysize0"]=fileService->make<TH1D>("keysize0","#of particles per event passing Mu7 leg statistics",10,0,10);
  histos2D_[ "ptTrigCand0"] =fileService->make< TH2D >("ptTrigCand0","Object vs. candidate_lower_p_{T} (GeV)",100, 0., 100., 100, 0., 100.);
  histos2D_[ "ptTrigCand0" ]->SetXTitle( "candidate p_{T} (GeV)" );
  histos2D_[ "ptTrigCand0" ]->SetYTitle( "object p_{T} (GeV)" );

}

// ------------ method called on each new Event  ------------
template<class T>
bool
TriggerObjectFilter8<T>::filter( edm::Event& iEvent, const edm::EventSetup& iSetup)
{
 // std::cout<< "enterfilter"<< std::endl;
  //create pointer to output collection
  std::auto_ptr<edm::RefVector<std::vector<T> > > recoObjColl(new edm::RefVector<std::vector<T> >);
  int index0 = 9999;
  //get reco objects
  edm::Handle<edm::RefVector<std::vector<T> > > recoObjs;
  iEvent.getByLabel(recoObjTag_, recoObjs);
  // Trigger Info
  edm::Handle<trigger::TriggerEvent> trgEvent;
  iEvent.getByLabel(triggerEventTag_,trgEvent);
  edm::Handle<edm::TriggerResults> pTrgResults;
  iEvent.getByLabel(triggerResultsTag_, pTrgResults);
  std::map<std::string, bool> triggerInMenu;
  std::string myHLTFilter = "";

  // get names of active HLT paths in this event
  std::vector<std::string> activeHLTPathsInThisEvent = hltConfig_.triggerNames();
  // loop over active HLT paths to search for desired path
  for (std::vector<std::string>::const_iterator iHLT = activeHLTPathsInThisEvent.begin(); iHLT != activeHLTPathsInThisEvent.end(); ++iHLT) { // active paths loop

    for (std::vector<edm::InputTag>::const_iterator iMyHLT = hltTags_.begin(); 
	 iMyHLT != hltTags_.end(); ++iMyHLT) {
      if ((*iMyHLT).label() == *iHLT) {
        //cout << "######## " << *iHLT << endl;
        myHLTFilter = (*iMyHLT).label();
	triggerInMenu[(*iMyHLT).label()] = true;
      //  std::cout << "(*iMyHLT).label() = " << (*iMyHLT).label() << std::endl;
 	//std::cout << "hltConfig_.prescaleValue(iEvent, iSetup, *iHLT) = ";
  	//std::cout << hltConfig_.prescaleValue(iEvent, iSetup, *iHLT) << std::endl;
      }
    }
  } // active paths loop
  edm::InputTag filterTag;
  // loop over these objects to see whether they match
  const trigger::TriggerObjectCollection& TOC( trgEvent->getObjects() );
  //choose the right sub-filter depending on the HLT path name
  std::vector<std::string> filters;
   try { filters = hltConfig_.moduleLabels( theRightHLTTag_.label() ); }
   catch (std::exception ex) { cout << "bad trigger\n"; }
   for(int i=0; i != trgEvent->sizeFilters(); ++i) {
     
     std::string label(trgEvent->filterTag(i).label());
     //std::cout << trgEvent->filterTag(i) << std::endl;
     if( label.find(theRightHLTSubFilter0_.label()) != std::string::npos )
       {
	 index0 = i;
       }
 
   }

   if (index0== 9999){
     index0 = 0;
     }
   const trigger::Keys& KEYS0(trgEvent->filterKeys(index0));
   const size_type nK0(KEYS0.size());
 //  int n=40;
  // double x[n];
   //double y[n];
   //for(int i=0; i<n;i++){
    // x[i]=-2.0+i*0.1;
    // y[i]=0.0;
  // }
   const edm::TriggerNames &trgNames = iEvent.triggerNames(*pTrgResults);
   const unsigned int trgIndex = trgNames.triggerIndex(myHLTFilter);
   bool firedHLT = (trgIndex < trgNames.size()) && (pTrgResults->accept(trgIndex));
   //for (unsigned int i = 0; i < trgNames.size(); i++)
   std::vector<unsigned int> passingRecoObjRefKeys0;
   std::vector<unsigned int> passingRecoObjRefKeys0_NoHLT;
 
   //bool isLooseMuon0=0;
   //bool isLooseMuon1=0;
   //not requiring HLT fired
   for(int ipart = 0; ipart != nK0; ++ipart) {
     const trigger::TriggerObject& TO0 = TOC[KEYS0[ipart]];
     for (typename edm::RefVector<std::vector<T> >::const_iterator iRecoObj =
                recoObjs->begin(); iRecoObj != recoObjs->end();
              ++iRecoObj) {
       if ((deltaR(**iRecoObj, TO0) < delRMatchingCut_) &&
                 (std::find(passingRecoObjRefKeys0_NoHLT.begin(), passingRecoObjRefKeys0_NoHLT.end(),
                          iRecoObj->key()) == passingRecoObjRefKeys0_NoHLT.end())) {
         passingRecoObjRefKeys0_NoHLT.push_back(iRecoObj->key());
         histos2D_[ "ptTrigCand0"]->Fill((*iRecoObj)->pt(),TO0.pt()); 
       }
     }
   }
   histos1D_["keysize0"]->Fill(passingRecoObjRefKeys0_NoHLT.size());
 

   if (firedHLT)
     { // firedHLT
       for(int ipart = 0; ipart != nK0; ++ipart) { 
	 const trigger::TriggerObject& TO0 = TOC[KEYS0[ipart]];
	 for (typename edm::RefVector<std::vector<T> >::const_iterator iRecoObj = 
	 	recoObjs->begin(); iRecoObj != recoObjs->end(); 
	      ++iRecoObj) {
    //         isLooseMuon0=muon::isLooseMuon(**iRecoObj);
         
	   if ((deltaR(**iRecoObj, TO0) < delRMatchingCut_) && 
	         (std::find(passingRecoObjRefKeys0.begin(), passingRecoObjRefKeys0.end(), 
	  		  iRecoObj->key()) == passingRecoObjRefKeys0.end())) {
	     recoObjColl->push_back(*iRecoObj);
        //     histos2D_[ "ptTrigCand0"]->Fill((*iRecoObj)->pt(),TO0.pt());
             passingRecoObjRefKeys0.push_back(iRecoObj->key());
	   }
	 }
       }
     }//firedH
   iEvent.put(recoObjColl);
   return (passingRecoObjRefKeys0.size() >= minNumObjsToPassFilter0_);
}

// ------------ method called when starting to processes a run  ------------


// ------------ method called once each job just before starting event loop  ------------

// ------------ method called once each job just after ending the event loop  ------------
template<class T>
void 
TriggerObjectFilter8<T>::endJob() {
}


// ------------ method called when ending the processing of a run  ------------
//template<class T>
//bool TriggerObjectFilter<T>::endRun(edm::Run&, edm::EventSetup const&)
//{
 // std::cout<< "endRun" << std::endl;
 // return true;
//} 
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
TriggerObjectFilter8<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in

typedef TriggerObjectFilter8<reco::Muon> MuonTriggerObjectFilter8;
DEFINE_FWK_MODULE(MuonTriggerObjectFilter8);
