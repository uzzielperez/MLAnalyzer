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

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "Geometry/CaloTopology/interface/EcalBarrelTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDetId/interface/HcalElectronicsId.h"
#include "DQM/HcalCommon/interface/Constants.h"

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
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
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
    edm::EDGetTokenT<EcalRecHitCollection> redEBRecHitCollectionT_;
    edm::EDGetTokenT<EcalRecHitCollection> EBRecHitCollectionT_;
    edm::EDGetTokenT<EBDigiCollection>     EBDigiCollectionT_;
    edm::EDGetTokenT<HBHERecHitCollection> HBHERecHitCollectionT_;
    //edm::InputTag trackTags_; //used to select what tracks to read from configuration file

    //TH1D * histo; 
    TH2D * hEBEnergy; 
    TH2D * hEBEnergyRed; 
    TH2D * hEBTiming; 
    TH2D * hEB_adc[EcalDataFrame::MAXSAMPLES]; 
    TH2D * hHBHEEnergy; 

    TH1D * hHBHEDepth; 
    TH1D * hHBHER; 
    TTree* RHTree;

    std::vector<float> vEBEnergy_;
    std::vector<float> vEBTiming_;
    std::vector<float> vEBEnergyRed_;
    std::vector<float> vEBTimingRed_;

    TCanvas *cEB;
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
{
  redEBRecHitCollectionT_ = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedEBRecHitCollection"));
  EBRecHitCollectionT_ = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EBRecHitCollection"));
  //EBDigiCollectionT_ = consumes<EBDigiCollection>(iConfig.getParameter<edm::InputTag>("selectedEBDigiCollection"));
  EBDigiCollectionT_ = consumes<EBDigiCollection>(iConfig.getParameter<edm::InputTag>("EBDigiCollection"));
  HBHERecHitCollectionT_ = consumes<HBHERecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedHBHERecHitCollection"));

  //now do what ever initialization is needed
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  //histo = fs->make<TH1D>("charge" , "Charges" , 200 , -2 , 2 );

  // Histograms
  // ECAL
  // Rechits
  cEB = new TCanvas("cEB","cEB",600,300);
  hEBEnergy = fs->make<TH2D>("EB_rechitE", "E(i#phi,i#eta);i#phi;i#eta",
      EBDetId::MAX_IPHI  , EBDetId::MIN_IPHI-1, EBDetId::MAX_IPHI,
      2*EBDetId::MAX_IETA,-EBDetId::MAX_IETA,   EBDetId::MAX_IETA );
  hEBEnergyRed = fs->make<TH2D>("EB_rechitEred", "Ered(i#phi,i#eta);i#phi;i#eta",
      EBDetId::MAX_IPHI  , EBDetId::MIN_IPHI-1, EBDetId::MAX_IPHI,
      2*EBDetId::MAX_IETA,-EBDetId::MAX_IETA,   EBDetId::MAX_IETA );
  hEBTiming = fs->make<TH2D>("EB_rechitT", "t(i#phi,i#eta);i#phi;i#eta",
      EBDetId::MAX_IPHI  , EBDetId::MIN_IPHI-1, EBDetId::MAX_IPHI,
      2*EBDetId::MAX_IETA,-EBDetId::MAX_IETA,   EBDetId::MAX_IETA );
  // Digis
  char hname[50], htitle[50];
  for(int iS(0); iS < EcalDataFrame::MAXSAMPLES; iS++){
    sprintf(hname, "EB_adc%d",iS);
    sprintf(htitle,"adc%d(i#phi,i#eta);i#phi;i#eta",iS);
    hEB_adc[iS] = fs->make<TH2D>(hname, htitle,
        EBDetId::MAX_IPHI  , EBDetId::MIN_IPHI-1, EBDetId::MAX_IPHI,
        2*EBDetId::MAX_IETA,-EBDetId::MAX_IETA,   EBDetId::MAX_IETA );
  }
  // HCAL
  hHBHEEnergy = fs->make<TH2D>("HBHE_rechitE", "E(i#phi,i#eta);i#phi;i#eta",
      hcaldqm::constants::IPHI_NUM, hcaldqm::constants::IPHI_MIN-1,hcaldqm::constants::IPHI_MAX,
      hcaldqm::constants::IETA_NUM,-hcaldqm::constants::IETA_MAX,  hcaldqm::constants::IETA_MAX );
  hHBHEDepth = fs->make<TH1D>("HBHE_depth", "E(i#phi,i#eta);i#phi;i#eta",
      hcaldqm::constants::DEPTH_NUM, hcaldqm::constants::DEPTH_MIN, hcaldqm::constants::DEPTH_MAX+1);
  hHBHER = fs->make<TH1D>("HB_r" , "r" , 20 , 0. , 0.);

  // Output Tree
  RHTree    = fs->make<TTree>("RHTree", "RecHit tree");
  RHTree->Branch("EBenergy",	&vEBEnergy_);
  RHTree->Branch("EBtime",		&vEBTiming_);
  RHTree->Branch("EBenergyRed",	&vEBEnergyRed_);
  RHTree->Branch("EBtimeRed",		&vEBTimingRed_);
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

  // Get Calo Geometry
  edm::ESHandle<CaloGeometry> caloGeomH;
  iSetup.get<CaloGeometryRecord>().get(caloGeomH);
  const CaloGeometry* caloGeom = caloGeomH.product();
  //const CaloSubdetectorGeometry* towerGeometry = 
  //geo->getSubdetectorGeometry(DetId::Calo, CaloTowerDetId::SubdetId);

  //////////// EB //////////

  // EB rechit collection
  vEBEnergy_.assign(EBDetId::kSizeForDenseIndexing,0.);
  vEBTiming_.assign(EBDetId::kSizeForDenseIndexing,0);
  int iphi_,ieta_,idx;
  edm::Handle<EcalRecHitCollection> EBRecHitsH;
  iEvent.getByToken(EBRecHitCollectionT_, EBRecHitsH);
  for(EcalRecHitCollection::const_iterator iRHit = EBRecHitsH->begin();
      iRHit != EBRecHitsH->end();                      
      ++iRHit) {
    EBDetId ebId( iRHit->id() );
    iphi_ = ebId.iphi()-1;
    ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta();
    hEBEnergy->Fill( iphi_,ieta_,iRHit->energy() );
    //hEBEnergy->Fill( iphi_,ieta_ );
    hEBTiming->Fill( iphi_,ieta_,iRHit->time() );
    //std::cout << iRHit->time() << std::endl;
    idx   = ebId.hashedIndex(); // (ieta_+EBDetId::MAX_IETA)*EBDetId::MAX_IPHI + iphi_
    vEBEnergy_[idx] = iRHit->energy(); // c.f. [iphi][ieta]
    vEBTiming_[idx] = iRHit->time();   // c.f. [iphi][ieta]
  }
  cEB->cd();
  gPad->SetLogz(1);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  hEBEnergy->GetZaxis()->SetRangeUser(2.e-2, 9.e1);
  hEBEnergy->Draw("COL Z");
  char outFile[100];
  sprintf(outFile,"cEBEnergy_%llu.eps",iEvent.id().event());
  cEB->Print(outFile);

  // EB reduced rechit collection
  vEBEnergyRed_.assign(EBDetId::kSizeForDenseIndexing,0.);
  vEBTimingRed_.assign(EBDetId::kSizeForDenseIndexing,0);
  edm::Handle<EcalRecHitCollection> redEBRecHitsH;
  iEvent.getByToken(redEBRecHitCollectionT_, redEBRecHitsH);
  for(EcalRecHitCollection::const_iterator iRHit = redEBRecHitsH->begin();
      iRHit != redEBRecHitsH->end();                      
      ++iRHit) {
    EBDetId ebId( iRHit->id() );
    iphi_ = ebId.iphi()-1;
    ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta();
    hEBEnergyRed->Fill( iphi_,ieta_,iRHit->energy() );
    //hEBEnergy->Fill( iphi_,ieta_ );
    //hEBTimingRed->Fill( iphi_,ieta_,iRHit->time() );
    idx   = ebId.hashedIndex(); // (ieta_+EBDetId::MAX_IETA)*EBDetId::MAX_IPHI + iphi_
    vEBEnergyRed_[idx] = iRHit->energy(); // c.f. [iphi][ieta]
    vEBTimingRed_[idx] = iRHit->time();   // c.f. [iphi][ieta]
  }
  cEB->Clear();
  hEBEnergyRed->GetZaxis()->SetRangeUser(2.e-2, 9.e1);
  hEBEnergyRed->Draw("COL Z");
  sprintf(outFile,"cEBEnergyRed_%llu.eps",iEvent.id().event());
  cEB->Print(outFile);

  // EB digis
  edm::Handle<EBDigiCollection> EBDigisH;
  iEvent.getByToken(EBDigiCollectionT_, EBDigisH);
  for(EBDigiCollection::const_iterator iDigi = EBDigisH->begin();
      iDigi != EBDigisH->end();                      
      ++iDigi) {
    //DetId id( iDigi->id() );
    EBDetId ebId( iDigi->id() );
    iphi_ = ebId.iphi()-1;
    ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta();
    EcalDataFrame df(*iDigi);
    for(int iS(0); iS < EcalDataFrame::MAXSAMPLES; ++iS) {
      EcalMGPASample digiSample( df.sample(iS) );
      hEB_adc[iS]->Fill( iphi_, ieta_, digiSample.adc() );
      //std::cout << digiSample.adc() << std::endl;
    }
  }
  for(int iS(0); iS < EcalDataFrame::MAXSAMPLES; ++iS) {
    cEB->Clear();
    gPad->SetLogz(1);
    hEB_adc[iS]->GetZaxis()->SetRangeUser(190, 1.3e3);
    hEB_adc[iS]->Draw("COL Z");
    sprintf(outFile,"cEB_adc%d_%llu.eps",iS,iEvent.id().event());
    cEB->Print(outFile);
  }

  //////////// HBHE //////////

  // HCAL
  GlobalPoint pos;
  edm::Handle<HBHERecHitCollection> HBHERecHitsH;
  iEvent.getByToken(HBHERecHitCollectionT_, HBHERecHitsH);
  for(HBHERecHitCollection::const_iterator iRHit = HBHERecHitsH->begin();
      iRHit != HBHERecHitsH->end();                      
      ++iRHit) {
    HcalDetId hId( iRHit->id() );
    if (hId.subdet() != HcalSubdetector::HcalBarrel) continue;
    iphi_ = hId.iphi()-1;
    ieta_ = hId.ieta() > 0 ? hId.ieta()-1 : hId.ieta();
    hHBHEEnergy->Fill( iphi_,ieta_,iRHit->energy() );
    pos = caloGeom->getPosition(hId);
    hHBHEDepth->Fill( hId.depth() );
    float x = pos.x();
    //std::cout << x  << std::endl;
    //hHBHER->Fill( pos.x() );
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
  RHTree->Fill();

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
