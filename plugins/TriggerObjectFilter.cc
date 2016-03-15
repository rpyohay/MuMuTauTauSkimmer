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
class TriggerObjectFilter : public edm::EDFilter {
   public:
      explicit TriggerObjectFilter(const edm::ParameterSet&);
      ~TriggerObjectFilter();

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
      edm::EDGetTokenT<edm::RefVector<std::vector<T> > > recoObjTag_;
      edm::EDGetTokenT<trigger::TriggerEvent>  triggerEventTag_;
      edm::EDGetTokenT<edm::TriggerResults> triggerResultsTag_;
      double Cut_;
      std::vector<edm::InputTag> hltTags_;
      HLTConfigProvider hltConfig_;
      edm::InputTag theRightHLTTag_;
      edm::InputTag theRightHLTSubFilter1_;
      std::vector<edm::InputTag> HLTSubFilters_;
      unsigned int minNumObjsToPassFilter1_;
      std::map<std::string, TH1D*> histos1D_;
      std::map<std::string, TH2D*> histos2D_;
      int numberHLT=1;
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
TriggerObjectFilter<T>::TriggerObjectFilter(const edm::ParameterSet& iConfig):
  recoObjTag_(consumes<edm::RefVector<std::vector<T> > >(iConfig.getParameter<edm::InputTag>("recoObjTag"))),
////////////////  triggerEventTag_(consumes<trigger::TriggerEvent>(iConfig.getUntrackedParameter<edm::InputTag>("triggerEventTag"))),//,dTriggerEventTag))),
  // By default, trigger results are labeled "TriggerResults" with process name "HLT" in the event.
//////////  triggerResultsTag_(consumes<edm::TriggerResults>(iConfig.getUntrackedParameter<edm::InputTag>("triggerResultsTag"))),//,dTriggerResults)))
  Cut_(iConfig.getUntrackedParameter<double>("MatchCut")),
  hltTags_(iConfig.getParameter<std::vector<edm::InputTag> >("hltTags")),
  hltConfig_(),
  //hltConfig_ = iConfig.getParameter<HLTConfigProvider>("hltConfig");
  theRightHLTTag_(iConfig.getParameter<edm::InputTag>("theRightHLTTag")),
  theRightHLTSubFilter1_(iConfig.getParameter<edm::InputTag>("theRightHLTSubFilter1")),
  //Whether using HLT trigger path name or the actual trigger filter name. Trigger path is default.
  HLTSubFilters_(iConfig.getUntrackedParameter<std::vector<edm::InputTag> >("HLTSubFilters",std::vector<edm::InputTag>())),
  minNumObjsToPassFilter1_(iConfig.getParameter<unsigned int>("minNumObjsToPassFilter1")),
  histos1D_(),
  histos2D_()
{
  const edm::InputTag dTriggerEventTag("hltTriggerSummaryAOD","","HLT");
  triggerEventTag_ = (consumes<trigger::TriggerEvent>(iConfig.getUntrackedParameter<edm::InputTag>("triggerEventTag",dTriggerEventTag)));
  const edm::InputTag dTriggerResults("TriggerResults","","HLT");
  triggerResultsTag_ = (consumes<edm::TriggerResults>(iConfig.getUntrackedParameter<edm::InputTag>("triggerResultsTag",dTriggerResults)));

  produces<edm::RefVector<std::vector<T> > >();
}

template<class T>
TriggerObjectFilter<T>::~TriggerObjectFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}
//template<class T>
//void TriggerObjectFilter<T>::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
//{
//  std::cout<< "beginLuminosityBlock"<< std::endl;
//  return true;
//}


//
// member functions
//
template<class T>
void TriggerObjectFilter<T>::beginRun(const edm::Run& iRun, edm::EventSetup const& iSetup)
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
TriggerObjectFilter<T>::beginJob()
{
  //std::cout<< "beginJob" << std::endl;
  edm::Service<TFileService> fileService;
  histos1D_[ "etaDistri_num" ]=fileService->make<TH1D>("etaDistri_num","eta distribution of higher pt muon with fired HLT and trigger-reco match",60,-3.,3.);
  histos1D_["etaDistri_de"]=fileService->make<TH1D>("etaDistri_de","eta distribution of all reco muon without trigger fired or trigger-reco match",60,-3.,3.);
  histos1D_["num_divide_de"]=fileService->make<TH1D>("num_divide_de","eta distribution of Mu45 trigger+trigger matching efficiency",60,-3.,3.);
  histos1D_[ "etaDistri_num1" ]=fileService->make<TH1D>("etaDistri_num1","eta distribution of reco muon with fired HLT no trigger-reco match",60,-3.,3.);
  histos1D_["num_divide_de1"]=fileService->make<TH1D>("num_divide_de1","eta distribution of Mu45 trigger efficiency(no trigger-reco match done)",60,-3.,3.);

  histos1D_["etaDistri_de2"]=fileService->make<TH1D>("etaDistri_de","eta distribution of all reco muon with Mu45 trigger-reco match, no HLT fired",60,-3.,3.);
  histos1D_["num_divide_de2"]=fileService->make<TH1D>("num_divide_de","eta distribution of Mu45 trigger efficiency(with trigger-reco match)",60,-3.,3.);
  histos1D_["ptDistriTriggerObj"]=fileService->make<TH1D>("ptDistriTriggerObj", "DrellY pt distribution of trigger muon when trigger is fired",100,0,500);
  histos1D_["ptDistriTriggerObj"]->SetXTitle("P_{T}(GeV)");

  histos1D_["ptDistriRecoObj"]=fileService->make<TH1D>("ptDistriRecoObj","DrellY pt distribution of reco muon when trigger is fired", 100,0,500);
  histos1D_["ptDistriRecoObj"]->SetXTitle("P_{T}(GeV)");

  histos1D_["ptDistriTriggerObjNotFired"]=fileService->make<TH1D>("ptDistriTriggerObjNOT", "DrellY pt distribution of trigger muon when trigger is NOT fired",100,0,500);
  histos1D_["ptDistriTriggerObjNotFired"]->SetXTitle("P_{T}(GeV)");
  histos1D_["keysize1"]=fileService->make<TH1D>("keysize1","#of particles per event passing Mu45 leg statistics",10,0,10);
  histos2D_["ptTrigCand1"] =fileService->make< TH2D >("ptTrigCand1","DrellY Trigger vs. reco muon with highest p_{T} (GeV)",150, 0., 150., 150, 0., 150.);
  histos2D_[ "ptTrigCand1" ]->SetXTitle( "reco muon p_{T} (GeV)" );
  histos2D_[ "ptTrigCand1" ]->SetYTitle( "trigger muon p_{T} (GeV)" );
  histos1D_["HLTRate"]=fileService->make<TH1D>("HLTRate", "# of events fire HLT trigger", 10,0,10);
  histos1D_["HLTNotFireRate"]=fileService->make<TH1D>("HLTNotFireRate", "# of events not fire HLT trigger", 10,0,10);

}

// ------------ method called on each new Event  ------------
template<class T>
bool
TriggerObjectFilter<T>::filter( edm::Event& iEvent, const edm::EventSetup& iSetup)
{
 // std::cout<< "***********enterTriggerfilter, begin of a event *******************"<< std::endl;
  //create pointer to output collection
  std::auto_ptr<edm::RefVector<std::vector<T> > > recoObjColl(new edm::RefVector<std::vector<T> >);
  int index1 = 9999;

  //get reco objects
  edm::Handle<edm::RefVector<std::vector<T> > > recoObjs;
  iEvent.getByToken(recoObjTag_, recoObjs);

  double max=0.0;
  //double eta_of_max=0.0;
  for (typename edm::RefVector<std::vector<T> >::const_iterator iRecoObj =
                recoObjs->begin(); iRecoObj != recoObjs->end();
              ++iRecoObj){
    if(((*iRecoObj)->pt())>0.0)
      histos1D_["etaDistri_de"]->Fill((*iRecoObj)->eta());   
    if(max<((*iRecoObj)->pt()))
    {
      max=(*iRecoObj)->pt();
    }
    else
      continue;
    
  } 
  // Trigger Info
  edm::Handle<trigger::TriggerEvent> trgEvent;
  iEvent.getByToken(triggerEventTag_, trgEvent);

  edm::Handle<edm::TriggerResults> pTrgResults;
  iEvent.getByToken(triggerResultsTag_, pTrgResults);

  std::map<std::string, bool> triggerInMenu;
  std::string myHLTFilter = "";

  // get names of active HLT paths in this event
  std::vector<std::string> activeHLTPathsInThisEvent = hltConfig_.triggerNames();
  // loop over active HLT paths to search for desired path
  for (std::vector<std::string>::const_iterator iHLT = activeHLTPathsInThisEvent.begin(); iHLT != activeHLTPathsInThisEvent.end(); ++iHLT) { // active paths loop

    for (std::vector<edm::InputTag>::const_iterator iMyHLT = hltTags_.begin(); iMyHLT != hltTags_.end(); ++iMyHLT) {
      if ((*iMyHLT).label() == *iHLT) {
       // cout << "######## " << *iHLT << endl;
        myHLTFilter = (*iMyHLT).label();
	triggerInMenu[(*iMyHLT).label()] = true;
       // std::cout << "(*iMyHLT).label() = " << (*iMyHLT).label() << std::endl;
// 	std::cout << "hltConfig_.prescaleValue(iEvent, iSetup, *iHLT) = ";
 // 	std::cout << hltConfig_.prescaleValue(iEvent, iSetup, *iHLT) << std::endl;
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
//   std::cout<<"trgEvent->sizeFilters()"<<trgEvent->sizeFilters()<<std::endl;
   for(int i=0; i != trgEvent->sizeFilters(); ++i) {
     
     std::string label(trgEvent->filterTag(i).label());
     if( label.find(theRightHLTSubFilter1_.label()) != std::string::npos )
       {
        // std::cout<<"label.find(theRightHLTSubFilter1_label())="<<label.find(theRightHLTSubFilter1_.label())<<std::endl;
        // std::cout<<"std::string::npos="<<std::string::npos<<std::endl;
         index1 = i;
       }
 
   }

   if (index1== 9999){
     index1 = 0;
     }
  // std::cout<<"index1="<<index1<<std::endl;
   const trigger::Keys& KEYS1(trgEvent->filterKeys(index1));
   const size_type nK1(KEYS1.size());
   const edm::TriggerNames &trgNames = iEvent.triggerNames(*pTrgResults);
   const unsigned int trgIndex = trgNames.triggerIndex(myHLTFilter);
  // std::cout<<"trgIndex"<<trgIndex<<std::endl;
   bool firedHLT = (trgIndex < trgNames.size()) && (pTrgResults->accept(trgIndex));
   std::vector<unsigned int> passingRecoObjRefKeys1;
   std::vector<unsigned int> passingRecoObjRefKeys1_NoHLT;
 
   //bool isLooseMuon0=0;
   //bool isLooseMuon1=0;
   //not requiring HLT fired
   for(int ipart1 = 0; ipart1 != nK1; ++ipart1) {
     const trigger::TriggerObject& TO1 = TOC[KEYS1[ipart1]];
     for (typename edm::RefVector<std::vector<T> >::const_iterator iRecoObj =
                recoObjs->begin(); iRecoObj != recoObjs->end();
              ++iRecoObj) {
       if ((deltaR(**iRecoObj, TO1) < Cut_) &&
                 (std::find(passingRecoObjRefKeys1_NoHLT.begin(), passingRecoObjRefKeys1_NoHLT.end(),
                          iRecoObj->key()) == passingRecoObjRefKeys1_NoHLT.end())) {
         passingRecoObjRefKeys1_NoHLT.push_back(iRecoObj->key());
         histos1D_["etaDistri_de2"]->Fill((*iRecoObj)->eta());
         //histos2D_[ "ptTrigCand1"]->Fill((*iRecoObj)->pt(), TO1.pt());

       }
     }
   }
   histos1D_["keysize1"]->Fill(passingRecoObjRefKeys1_NoHLT.size());
 
   if(!firedHLT)
     {
       histos1D_["HLTNotFireRate"]->Fill(numberHLT);
       for(int ipart1=0; ipart1 != nK1; ++ipart1){
         const trigger::TriggerObject& TO1 =TOC[KEYS1[ipart1]];
      // std::cout<<"HLT not fired! TO1.pt()="<<TO1.pt()<<"\n TO1.eta()"<<TO1.eta()<<std::endl;
       
       histos1D_["ptDistriTriggerObjNotFired"]->Fill(TO1.pt());  
       for (typename edm::RefVector<std::vector<T> >::const_iterator iRecoObj =
                recoObjs->begin(); iRecoObj != recoObjs->end();
              ++iRecoObj) {
           if ((deltaR(**iRecoObj, TO1) < Cut_) &&
                 (std::find(passingRecoObjRefKeys1_NoHLT.begin(), passingRecoObjRefKeys1_NoHLT.end(),
                          iRecoObj->key()) == passingRecoObjRefKeys1_NoHLT.end())) {
           }
         }
       }
    }

   if (firedHLT)
     { // firedHLT
         histos1D_["HLTRate"]->Fill(numberHLT);
         for(int ipart1 = 0; ipart1 != nK1; ++ipart1) {

           const trigger::TriggerObject& TO1 = TOC[KEYS1[ipart1]];
           histos1D_["ptDistriTriggerObj"]->Fill(TO1.pt());
           for (typename edm::RefVector<std::vector<T> >::const_iterator iRecoObj =
                recoObjs->begin(); iRecoObj != recoObjs->end();
              ++iRecoObj) {
              histos1D_["ptDistriRecoObj"]->Fill((*iRecoObj)->pt());
              histos1D_["etaDistri_num1"]->Fill((*iRecoObj)->eta());
      //       isLooseMuon1=muon::isLooseMuon(**iRecoObj);
             if (/*(abs((*iRecoObj)->pt()- TO1.pt())/((*iRecoObj)->pt()) < Cut_)*/(deltaR(**iRecoObj, TO1) < Cut_) &&
                 (std::find(passingRecoObjRefKeys1.begin(), passingRecoObjRefKeys1.end(),
                            iRecoObj->key()) == passingRecoObjRefKeys1.end())) {

              
               histos1D_["etaDistri_num"]->Fill((*iRecoObj)->eta());
               //if((*iRecoObj)->eta()>(-2.0+i*0.1) && (*iRecoObj)->eta() < (-2.0+(i+1)*0.1)){
                // y[i]+=1;
              // }
               histos2D_[ "ptTrigCand1"]->Fill((*iRecoObj)->pt(), TO1.pt());
               passingRecoObjRefKeys1.push_back(iRecoObj->key());
             }
           }
         } 
     }//firedH
   iEvent.put(recoObjColl);
   return (passingRecoObjRefKeys1.size() >= minNumObjsToPassFilter1_);
}

// ------------ method called when starting to processes a run  ------------


// ------------ method called once each job just before starting event loop  ------------

// ------------ method called once each job just after ending the event loop  ------------
template<class T>
void 
TriggerObjectFilter<T>::endJob() {
   histos1D_["num_divide_de"]->Divide(  histos1D_[ "etaDistri_num" ],  histos1D_[ "etaDistri_de" ]);
   histos1D_["num_divide_de1"]->Divide(  histos1D_[ "etaDistri_num1" ],  histos1D_[ "etaDistri_de" ]);
   histos1D_["num_divide_de2"]->Divide(  histos1D_[ "etaDistri_num" ],  histos1D_[ "etaDistri_de2" ]);//look best
}


// ------------ method called when ending the processing of a run  ------------
template<class T>
void TriggerObjectFilter<T>::endRun(const edm::Run& iRun, edm::EventSetup const& iSetup)
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
TriggerObjectFilter<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
typedef TriggerObjectFilter<reco::Muon> MuonTriggerObjectFilter;
DEFINE_FWK_MODULE(MuonTriggerObjectFilter);
