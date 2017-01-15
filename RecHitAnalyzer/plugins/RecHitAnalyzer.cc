// -*- C++ -*-
//
// Package:    MLAnalyzer/RecHitAnalyzer
// Class:      RecHitAnalyzer
// 
/**\class RecHitAnalyzer RecHitAnalyzer.cc MLAnalyzer/RecHitAnalyzer/plugins/RecHitAnalyzer.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Michael Andrews
//         Created:  Sat, 14 Jan 2017 17:45:54 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "Geometry/CaloTopology/interface/EcalBarrelTopology.h"

#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDetId/interface/HcalElectronicsId.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TH2.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class RecHitAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
	public:
		explicit RecHitAnalyzer(const edm::ParameterSet&);
		~RecHitAnalyzer();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;

		// ----------member data ---------------------------
		edm::EDGetTokenT<EcalRecHitCollection> EBRecHitCollectionT_;
		edm::EDGetTokenT<HBHERecHitCollection> HBHERecHitCollectionT_;
		//edm::InputTag trackTags_; //used to select what tracks to read from configuration file
		//TH1D * histo; 
		TH2D * hEBEnergy; 
		TH2D * hEBTiming; 
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
RecHitAnalyzer::RecHitAnalyzer(const edm::ParameterSet& iConfig)
//	:
//		trackTags_(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))
{
	EBRecHitCollectionT_ = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedEBRecHitCollection"));
	HBHERecHitCollectionT_ = consumes<HBHERecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedHBHERecHitCollection"));
	//now do what ever initialization is needed
	usesResource("TFileService");
	edm::Service<TFileService> fs;
	//histo = fs->make<TH1D>("charge" , "Charges" , 200 , -2 , 2 );
	hEBEnergy = fs->make<TH2D>("E_rechit" , "E(iphi,ieta)" , EBDetId::MAX_IPHI , EBDetId::MIN_IPHI , EBDetId::MAX_IPHI, EBDetId::MAX_IETA , EBDetId::MIN_IETA , EBDetId::MAX_IETA );
	hEBTiming = fs->make<TH2D>("t_rechit" , "t(iphi,ieta)" , EBDetId::MAX_IPHI , EBDetId::MIN_IPHI , EBDetId::MAX_IPHI, EBDetId::MAX_IETA , EBDetId::MIN_IETA , EBDetId::MAX_IETA );

}


RecHitAnalyzer::~RecHitAnalyzer()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
	void
RecHitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;

	// ECAL Barrel
	edm::Handle<EcalRecHitCollection> EBRecHitsH;
	iEvent.getByToken(EBRecHitCollectionT_, EBRecHitsH);
	for(EcalRecHitCollection::const_iterator iRHit = EBRecHitsH->begin();
			iRHit != EBRecHitsH->end();                      
			++iRHit) {
		EBDetId ebId( iRHit->id() );
		hEBEnergy->Fill( ebId.iphi(),ebId.ieta(),iRHit->energy() );
		hEBTiming->Fill( ebId.iphi(),ebId.ieta(),iRHit->time() );
	}

	// HCAL Barrel
	edm::Handle<HBHERecHitCollection> HBHERecHitsH;
	iEvent.getByToken(HBHERecHitCollectionT_, HBHERecHitsH);
	for(HBHERecHitCollection::const_iterator iRHit = HBHERecHitsH->begin();
			iRHit != HBHERecHitsH->end();                      
			++iRHit) {
		HcalDetId hId( iRHit->id() );
	}

	/*
	using reco::TrackCollection;
	Handle<TrackCollection> tracks;
	iEvent.getByLabel(trackTags_,tracks);
	for(TrackCollection::const_iterator itTrack = tracks->begin();
			itTrack != tracks->end();                      
			++itTrack) {
		int charge = 0;
		charge = itTrack->charge();  
		histo->Fill( charge );
	}
	*/

#ifdef THIS_IS_AN_EVENT_EXAMPLE
	Handle<ExampleData> pIn;
	iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
	ESHandle<SetupData> pSetup;
	iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
	void 
RecHitAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
RecHitAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
RecHitAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);

	//Specify that only 'tracks' is allowed
	//To use, remove the default given above and uncomment below
	//ParameterSetDescription desc;
	//desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
	//descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RecHitAnalyzer);
