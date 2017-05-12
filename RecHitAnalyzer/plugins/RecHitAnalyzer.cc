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
#include "DataFormats/EcalDetId/interface/EcalTrigTowerDetId.h"
#include "Geometry/CaloTopology/interface/EcalBarrelTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
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
#include "TMath.h"
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
    TH2D * hHBHEEnergy_EB; 

    TH1D * hHBHEDepth; 
    TH1D * hHBHER; 
    TTree* RHTree;

    std::vector<float> vEBEnergy_;
    std::vector<float> vEBTiming_;
    std::vector<float> vEBEnergyRed_;
    std::vector<float> vEBTimingRed_;
    std::vector<float> vEB_adc_[EcalDataFrame::MAXSAMPLES];

    //std::vector<float> vHBHEEnergy[hcaldqm::constants::DEPTH_NUM];
    std::vector<float> vHBHEEnergy_EB_;
  
    const unsigned HBHE_IETA_MAX = hcaldqm::constants::IETA_MAX_HB + 1;//17
    //const unsigned HBHE_IETA_MAX = 20;


    TCanvas *cEB, *cHBHE;
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
  //EBDigiCollectionT_ = consumes<EBDigiCollection>(iConfig.getParameter<edm::InputTag>("EBDigiCollection"));
  HBHERecHitCollectionT_ = consumes<HBHERecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedHBHERecHitCollection"));

  //now do what ever initialization is needed
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  //histo = fs->make<TH1D>("charge" , "Charges" , 200 , -2 , 2 );

  // Histograms
  // ECAL
  // Rechits
  cEB = new TCanvas("cEB","cEB",600,300);
  cHBHE = new TCanvas("cHBHE","cHBHE",600,300);
  hEBEnergy = fs->make<TH2D>("EB_rechitE", "E(i#phi,i#eta);i#phi;i#eta",
      EBDetId::MAX_IPHI  , EBDetId::MIN_IPHI-1, EBDetId::MAX_IPHI,
      2*EBDetId::MAX_IETA,-EBDetId::MAX_IETA,   EBDetId::MAX_IETA );
  hEBEnergyRed = fs->make<TH2D>("EB_rechitEred", "Ered(i#phi,i#eta);i#phi;i#eta",
      EBDetId::MAX_IPHI  , EBDetId::MIN_IPHI-1, EBDetId::MAX_IPHI,
      2*EBDetId::MAX_IETA,-EBDetId::MAX_IETA,   EBDetId::MAX_IETA );
  //hEBEnergyRed = fs->make<TH2D>("EB_rechitEred", "Ered(i#phi,i#eta);i#phi;i#eta",
  //    EcalTrigTowerDetId::kEBTowersInPhi*18  , EBDetId::MIN_IPHI-1, EBDetId::MAX_IPHI,
  //    EcalTrigTowerDetId::kEBTowersInEta*2,-EBDetId::MAX_IETA,   EBDetId::MAX_IETA );
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
  hHBHEEnergy_EB = fs->make<TH2D>("HBHE_rechitE_EB", "E(i#phi,i#eta);i#phi;i#eta",
      hcaldqm::constants::IPHI_NUM, hcaldqm::constants::IPHI_MIN-1,hcaldqm::constants::IPHI_MAX,
      2*HBHE_IETA_MAX,             -HBHE_IETA_MAX,                 HBHE_IETA_MAX );
  hHBHEDepth = fs->make<TH1D>("HBHE_depth", "E(i#phi,i#eta);i#phi;i#eta",
      hcaldqm::constants::DEPTH_NUM, hcaldqm::constants::DEPTH_MIN, hcaldqm::constants::DEPTH_MAX+1);
  hHBHER = fs->make<TH1D>("HB_r" , "r" , 20 , 0. , 0.);

  // Output Tree
  RHTree    = fs->make<TTree>("RHTree", "RecHit tree");
  RHTree->Branch("EBenergy",	&vEBEnergy_);
  RHTree->Branch("EBtime",		&vEBTiming_);
  RHTree->Branch("EBenergyRed",	&vEBEnergyRed_);
  RHTree->Branch("EBtimeRed",		&vEBTimingRed_);
  for(int iS(0); iS < EcalDataFrame::MAXSAMPLES; iS++){
    sprintf(hname, "EB_adc%d",iS);
    RHTree->Branch(hname, &vEB_adc_[iS]);
  }
  RHTree->Branch("HBHEenergy_EB", &vHBHEEnergy_EB_);
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
  GlobalPoint pos;
  const CaloCellGeometry *cell;

  //////////// EB //////////

  bool saveImgs = false;
  int iphi_,ieta_,idx;

  // EB rechit collection
  vEBEnergy_.assign(EBDetId::kSizeForDenseIndexing,0.);
  vEBTiming_.assign(EBDetId::kSizeForDenseIndexing,0);
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
    cell  = caloGeom->getGeometry(ebId);
    idx   = ebId.hashedIndex(); // (ieta_+EBDetId::MAX_IETA)*EBDetId::MAX_IPHI + iphi_
    //vEBEnergy_[idx] = iRHit->energy(); // c.f. [ieta][iphi]
    vEBEnergy_[idx] = iRHit->energy()/std::abs(cell->etaPos()); // c.f. [ieta][iphi]
    vEBTiming_[idx] = iRHit->time();   // c.f. [ieta][iphi]
  }
  cEB->cd();
  gPad->SetLogz(1);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  //hEBEnergy->GetZaxis()->SetRangeUser(2.e-2, 9.e1);
  //hEBEnergy->Draw("COL Z");
  char outFile[100];
  sprintf(outFile,"cEBEnergy_%llu.eps",iEvent.id().event());
  if (saveImgs) cEB->Print(outFile);

  // EB reduced rechit collection
  vEBEnergyRed_.assign(EBDetId::kSizeForDenseIndexing,0.);
  vEBTimingRed_.assign(EBDetId::kSizeForDenseIndexing,0);
  //vEBEnergyRed_.assign(EcalTrigTowerDetId::kEBTotalTowers,0.);
  //vEBTimingRed_.assign(EcalTrigTowerDetId::kEBTotalTowers,0);
  edm::Handle<EcalRecHitCollection> redEBRecHitsH;
  iEvent.getByToken(redEBRecHitCollectionT_, redEBRecHitsH);
  for(EcalRecHitCollection::const_iterator iRHit = redEBRecHitsH->begin();
      iRHit != redEBRecHitsH->end();                      
      ++iRHit) {
    EBDetId ebId( iRHit->id() );
    iphi_ = ebId.iphi()-1;
    ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta();
    hEBEnergyRed->Fill( iphi_,ieta_,iRHit->energy() );
    //hEBTimingRed->Fill( iphi_,ieta_,iRHit->time() );
    cell  = caloGeom->getGeometry(ebId);
    idx   = ebId.hashedIndex(); // (ieta_+EBDetId::MAX_IETA)*EBDetId::MAX_IPHI + iphi_
    //EcalTrigTowerDetId ttId( iRHit->id() );
    //idx = ttId.hashedIndex();
    //vEBEnergyRed_[idx] = iRHit->energy(); // c.f. [ieta][iphi]
    vEBEnergyRed_[idx] = iRHit->energy()/TMath::CosH(cell->etaPos()); // c.f. [ieta][iphi]
    vEBTimingRed_[idx] = iRHit->time();   // c.f. [ieta][iphi]
  }
  cEB->Clear();
  //hEBEnergyRed->GetZaxis()->SetRangeUser(2.e-2, 9.e1);
  //hEBEnergyRed->Draw("COL Z");
  //char outFile[100];
  sprintf(outFile,"cEBEnergyRed_%llu.eps",iEvent.id().event());
  if (saveImgs) cEB->Print(outFile);

  /*
  // EB digis
  for(int iS(0); iS < EcalDataFrame::MAXSAMPLES; ++iS)
    vEB_adc_[iS].assign(EBDetId::kSizeForDenseIndexing,0);
  edm::Handle<EBDigiCollection> EBDigisH;
  iEvent.getByToken(EBDigiCollectionT_, EBDigisH);
  for(EBDigiCollection::const_iterator iDigi = EBDigisH->begin();
      iDigi != EBDigisH->end();                      
      ++iDigi) {
    //DetId id( iDigi->id() );
    EBDetId ebId( iDigi->id() );
    iphi_ = ebId.iphi()-1;
    ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta();
    idx = ebId.hashedIndex(); // (ieta_+EBDetId::MAX_IETA)*EBDetId::MAX_IPHI + iphi_
    cell  = caloGeom->getGeometry(ebId);
    EcalDataFrame df(*iDigi);
    for(int iS(0); iS < EcalDataFrame::MAXSAMPLES; ++iS) {
      EcalMGPASample digiSample( df.sample(iS) );
      hEB_adc[iS]->Fill( iphi_, ieta_, digiSample.adc() );
      //std::cout << digiSample.adc() << std::endl;
      vEB_adc_[iS][idx] += digiSample.adc()/TMath::CosH(cell->etaPos()); // c.f. [ieta][iphi]
      //vEB_adc_[iS][idx] += digiSample.adc(); // c.f. [ieta][iphi]
    }
  }
  for(int iS(0); iS < EcalDataFrame::MAXSAMPLES; ++iS) {
    cEB->Clear();
    gPad->SetLogz(1);
    //hEB_adc[iS]->GetZaxis()->SetRangeUser(190, 1.3e3);
    //hEB_adc[iS]->Draw("COL Z");
    sprintf(outFile,"cEB_adc%d_%llu.eps",iS,iEvent.id().event());
    if (saveImgs) cEB->Print(outFile);
  }
  */

  //////////// HBHE //////////

  // HCAL
  cHBHE->cd();
  float maxEta = 0.;
  vHBHEEnergy_EB_.assign( hcaldqm::constants::IPHI_NUM*2*HBHE_IETA_MAX,0. );
  edm::Handle<HBHERecHitCollection> HBHERecHitsH;
  iEvent.getByToken(HBHERecHitCollectionT_, HBHERecHitsH);
  for(HBHERecHitCollection::const_iterator iRHit = HBHERecHitsH->begin();
      iRHit != HBHERecHitsH->end();                      
      ++iRHit) {
    HcalDetId hId( iRHit->id() );
    //if (hId.subdet() != HcalSubdetector::HcalBarrel) continue;
    iphi_ = hId.iphi()-1;
    ieta_ = hId.ieta() > 0 ? hId.ieta()-1 : hId.ieta();
    hHBHEEnergy->Fill( iphi_,ieta_,iRHit->energy() );
    hHBHEDepth->Fill( hId.depth() );

    if ( abs(hId.ieta()) > HBHE_IETA_MAX ) continue; 
    pos  = caloGeom->getPosition(hId);
    cell = caloGeom->getGeometry(hId);
    float x = pos.x();
    float rho = cell->rhoPos();
    float eta = cell->etaPos();
    float phi = cell->phiPos();
    hHBHEEnergy_EB->Fill( iphi_,ieta_,iRHit->energy() );

    idx = ( ieta_+(HBHE_IETA_MAX) )*hcaldqm::constants::IPHI_NUM + iphi_;
    vHBHEEnergy_EB_[idx] += iRHit->energy(); // c.f. [ieta][iphi]
    //vHBHEEnergy_EB_[idx] = iRHit->energy()/TMath::CosH(cell->etaPos()); // c.f. [ieta][iphi]

    //if (iRHit->energy() > 0.) {
      //std::cout << x  << std::endl;
      //std::cout << rho << ":" << eta << ":" << phi << std::endl;
    //}
    //if (eta > std::abs(maxEta)) maxEta = std::abs(eta);
    //hHBHER->Fill( pos.x() );
  }
  //std::cout << "maxEta: " << maxEta << std::endl;
  hHBHEEnergy->GetZaxis()->SetRangeUser(2.e-2, 9.e1);
  hHBHEEnergy->Draw("COL Z");
  sprintf(outFile,"cHBHEEnergy_%llu.eps",iEvent.id().event());
  if (saveImgs) cHBHE->Print(outFile);

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
