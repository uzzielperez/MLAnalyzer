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

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

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
    edm::EDGetTokenT<EBDigiCollection>     EBDigiCollectionT_;
    edm::EDGetTokenT<EcalRecHitCollection> EERecHitCollectionT_;
    edm::EDGetTokenT<HBHERecHitCollection> HBHERecHitCollectionT_;
    edm::EDGetTokenT<reco::GenParticleCollection> genParticleCollectionT_;
    //edm::InputTag trackTags_; //used to select what tracks to read from configuration file

    static const int EE_IZ_MAX = 2;
    //const unsigned HBHE_IETA_MAX = hcaldqm::constants::IETA_MAX_HB + 1;//17
    const unsigned HBHE_IETA_MAX = 20;
    const unsigned EE_NC_PER_ZSIDE = EEDetId::IX_MAX*EEDetId::IY_MAX; // 100*100

    // Initialize Calorimeter Geometry
    const CaloGeometry* caloGeom;
    //const CaloSubdetectorGeometry* towerGeom; 

    // Initializer global cell position
    //GlobalPoint pos;

    //TH1D * histo; 
    TH2D * hEB_energy; 
    TH2D * hEB_time; 
    TH2D * hEB_adc[EcalDataFrame::MAXSAMPLES]; 
    TH2D * hEE_energy[EE_IZ_MAX]; 
    TH2D * hEE_time[EE_IZ_MAX]; 

    TH2D * hHBHE_energy; 
    TH2D * hHBHE_energy_EB; 
    TH1D * hHBHE_depth; 

    TH1D * h_pT; 
    TH1D * h_E; 
    TH1D * h_eta; 
    TH1D * h_m0; 

    TTree* RHTree;

    float eventId_;
    float m0_;
    std::vector<float> vEB_energy_;
    std::vector<float> vEB_time_;
    std::vector<float> vEB_adc_[EcalDataFrame::MAXSAMPLES];
    std::vector<float> vEE_energy_[EE_IZ_MAX];
    std::vector<float> vEE_time_[EE_IZ_MAX];

    //std::vector<float> vHBHE_energy[hcaldqm::constants::DEPTH_NUM];
    std::vector<float> vHBHE_energy_EB_;
  
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
  //EBRecHitCollectionT_ = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EBRecHitCollection"));
  EBRecHitCollectionT_ = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedEBRecHitCollection"));
  //EBDigiCollectionT_ = consumes<EBDigiCollection>(iConfig.getParameter<edm::InputTag>("selectedEBDigiCollection"));
  //EBDigiCollectionT_ = consumes<EBDigiCollection>(iConfig.getParameter<edm::InputTag>("EBDigiCollection"));
  EERecHitCollectionT_ = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedEERecHitCollection"));
  //EERecHitCollectionT_ = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EERecHitCollection"));
  HBHERecHitCollectionT_ = consumes<HBHERecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedHBHERecHitCollection"));

  genParticleCollectionT_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticleCollection"));

  //now do what ever initialization is needed
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  //histo = fs->make<TH1D>("charge" , "Charges" , 200 , -2 , 2 );

  // Canvases
  cEB = new TCanvas("cEB","cEB",600,300);
  cHBHE = new TCanvas("cHBHE","cHBHE",600,300);

  // Histograms
  // EB rechits
  hEB_energy = fs->make<TH2D>("EB_energy", "E(i#phi,i#eta);i#phi;i#eta",
      EBDetId::MAX_IPHI  , EBDetId::MIN_IPHI-1, EBDetId::MAX_IPHI,
      2*EBDetId::MAX_IETA,-EBDetId::MAX_IETA,   EBDetId::MAX_IETA );
  //hEB_energy = fs->make<TH2D>("EB_energy", "E(i#phi,i#eta);i#phi;i#eta",
  //    EcalTrigTowerDetId::kEBTowersInPhi*18  , EBDetId::MIN_IPHI-1, EBDetId::MAX_IPHI,
  //    EcalTrigTowerDetId::kEBTowersInEta*2,-EBDetId::MAX_IETA,   EBDetId::MAX_IETA );
  hEB_time = fs->make<TH2D>("EB_time", "t(i#phi,i#eta);i#phi;i#eta",
      EBDetId::MAX_IPHI  , EBDetId::MIN_IPHI-1, EBDetId::MAX_IPHI,
      2*EBDetId::MAX_IETA,-EBDetId::MAX_IETA,   EBDetId::MAX_IETA );
  // EB Digis
  char hname[50], htitle[50];
  /*
  for(int iS(0); iS < EcalDataFrame::MAXSAMPLES; iS++){
    sprintf(hname, "EB_adc%d",iS);
    sprintf(htitle,"adc%d(i#phi,i#eta);i#phi;i#eta",iS);
    hEB_adc[iS] = fs->make<TH2D>(hname, htitle,
        EBDetId::MAX_IPHI  , EBDetId::MIN_IPHI-1, EBDetId::MAX_IPHI,
        2*EBDetId::MAX_IETA,-EBDetId::MAX_IETA,   EBDetId::MAX_IETA );
  }
  */

  // EE rechits
  for(int iz(0); iz < EE_IZ_MAX; iz++){
    const char *zside = (iz > 0) ? "p" : "m";
    sprintf(hname, "EE%s_energy",zside);
    sprintf(htitle,"E(ix,iy);ix;iy");
    hEE_energy[iz] = fs->make<TH2D>(hname, htitle,
      EEDetId::IX_MAX, EEDetId::IX_MIN-1, EEDetId::IX_MAX,
      EEDetId::IY_MAX, EEDetId::IY_MIN-1, EEDetId::IY_MAX );
    sprintf(hname, "EE%s_time",zside);
    sprintf(htitle,"t(ix,iy);ix;iy");
    hEE_time[iz] = fs->make<TH2D>(hname, htitle,
      EEDetId::IX_MAX, EEDetId::IX_MIN-1, EEDetId::IX_MAX,
      EEDetId::IY_MAX, EEDetId::IY_MIN-1, EEDetId::IY_MAX );
  } 

  // HBHE
  hHBHE_energy = fs->make<TH2D>("HBHE_energy", "E(i#phi,i#eta);i#phi;i#eta",
      hcaldqm::constants::IPHI_NUM,      hcaldqm::constants::IPHI_MIN-1, hcaldqm::constants::IPHI_MAX,
      2*hcaldqm::constants::IETA_MAX_HE,-hcaldqm::constants::IETA_MAX_HE,hcaldqm::constants::IETA_MAX_HE );
  hHBHE_energy_EB = fs->make<TH2D>("HBHE_energy_EB", "E(i#phi,i#eta);i#phi;i#eta",
      hcaldqm::constants::IPHI_NUM, hcaldqm::constants::IPHI_MIN-1,hcaldqm::constants::IPHI_MAX,
      2*HBHE_IETA_MAX,             -HBHE_IETA_MAX,                 HBHE_IETA_MAX );
  hHBHE_depth = fs->make<TH1D>("HBHE_depth", "Depth;depth;Hits",
      hcaldqm::constants::DEPTH_NUM, hcaldqm::constants::DEPTH_MIN, hcaldqm::constants::DEPTH_MAX+1);

  // Kinematics
  h_pT  = fs->make<TH1D>("h_pT" , "p_{T};p_{T};Particles", 100,  0., 500.);
  h_E   = fs->make<TH1D>("h_E"  , "E;E;Particles"        , 100,  0., 800.);
  h_eta = fs->make<TH1D>("h_eta", "#eta;#eta;Particles"  ,  50,-2.4, 2.4);
  h_m0  = fs->make<TH1D>("h_m0" , "m0;m0;Events"        ,   72, 50., 950.);

  // Output Tree
  RHTree = fs->make<TTree>("RHTree", "RecHit tree");
  RHTree->Branch("eventId",      &eventId_);
  RHTree->Branch("m0",           &m0_);
  RHTree->Branch("EB_energy",    &vEB_energy_);
  RHTree->Branch("EB_time",      &vEB_time_);
  /*
  for(int iS(0); iS < EcalDataFrame::MAXSAMPLES; iS++){
    sprintf(hname, "EB_adc%d",iS);
    RHTree->Branch(hname,       &vEB_adc_[iS]);
  }
  */
  for(int iz(0); iz < EE_IZ_MAX; iz++){
    const char *zside = (iz > 0) ? "p" : "m";
    sprintf(hname, "EE%s_energy",zside);
    RHTree->Branch(hname,       &vEE_energy_[iz]);
    sprintf(hname, "EE%s_time",zside);
    RHTree->Branch(hname,       &vEE_time_[iz]);
  }
  RHTree->Branch("HBHE_energy_EB", &vHBHE_energy_EB_);
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

  int nPho = 0;
  bool isDecayed = true;
  bool isHiggs = false;
  float etaCut = 1.4;
  float ptCut = 25.;
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleCollectionT_, genParticles);
  for (reco::GenParticleCollection::const_iterator iP = genParticles->begin();
       iP != genParticles->end();
       ++iP) {

    // PDG ID cut
    if ( std::abs(iP->pdgId()) != 22 ) continue;
    
    // Decay status
    //if ( iP->status() != 23 ) continue;
    if ( iP->status() != 1 ) continue; // NOT the same as Pythia status

    if ( isDecayed ) {
        // Check ancestry
        if ( !iP->mother() ) continue;
        if ( isHiggs ) {
            if ( std::abs(iP->mother()->pdgId()) != 25 && std::abs(iP->mother()->pdgId()) != 22 ) continue;
        } else {
            if ( std::abs(iP->mother()->status()) != 44 && std::abs(iP->mother()->status()) != 23 ) continue;
        }
        // Kinematic cuts
        if ( std::abs(iP->eta()) > etaCut ) continue;
        if ( std::abs(iP->pt()) < ptCut ) continue;
    } // apply cuts

    nPho++;

  } // genParticle loop: count good photons

  // Require exactly 2 gen-level photons
  // Indifferent about photons of status != 1
  if ( nPho != 2 ) return; 
  
  // Fill loop
  math::XYZTLorentzVector vDiPho;
  for (reco::GenParticleCollection::const_iterator iP = genParticles->begin();
       iP != genParticles->end();
       ++iP) {

    // PDG ID cut
    if ( std::abs(iP->pdgId()) != 22 ) continue;
    
    // Decay status
    //if ( iP->status() != 23 ) continue;
    if ( iP->status() != 1 ) continue;

    if ( isDecayed ) {
        // Check ancestry
        if ( isHiggs ) {
            if ( std::abs(iP->mother()->pdgId()) != 25 && std::abs(iP->mother()->pdgId()) != 22 ) continue;
        } else {
            if ( std::abs(iP->mother()->status()) != 44 && std::abs(iP->mother()->status()) != 23 ) continue;
        }
        // Kinematic cuts
        if ( std::abs(iP->eta()) > etaCut ) continue;
        if ( std::abs(iP->pt()) < ptCut ) continue;
    } // apply cuts

    if ( isDecayed ) { 
        std::cout << "status:" <<iP->status() << " pT:" << iP->pt() << " eta:" << iP->eta() << " E:" << iP->energy() << " mothId:" << iP->mother()->pdgId() << std::endl;
    } else {
        std::cout << "status:" <<iP->status() << " pT:" << iP->pt() << " eta:" << iP->eta() << " E:" << iP->energy() << std::endl;
    }
    vDiPho += iP->p4();
    // Fill histograms
    h_pT-> Fill( iP->pt()      );
    h_E->  Fill( iP->energy()  );
    h_eta->Fill( iP->eta()     );
  } // genParticle loop: fill hist
  h_m0->Fill( vDiPho.T() );
  std::cout << " m0: " << vDiPho.T() << std::endl;

  // ----- Get Calorimeter Geometry ----- //
  // Provides access to global cell position and coordinates below
  edm::ESHandle<CaloGeometry> caloGeomH;
  iSetup.get<CaloGeometryRecord>().get(caloGeomH);
  caloGeom = caloGeomH.product();
  //towerGeom = caloGeomH.product(); 
  //towerGeom->getSubdetectorGeometry(DetId::Calo, CaloTowerDetId::SubdetId);

  //////////// EB //////////

  bool saveImgs = false;
  int iphi_, ieta_, idx;

  // ----- EB reduced rechit collection ----- //
  // This contatins the reduced EB rechit collection after
  // the zero suppression and bad channel clean-up

  // Initialize arrays
  vEB_energy_.assign(EBDetId::kSizeForDenseIndexing,0.);
  vEB_time_.assign(EBDetId::kSizeForDenseIndexing,0.);
  //vEB_energy_.assign(EcalTrigTowerDetId::kEBTotalTowers,0.);
  //vEB_time_.assign(EcalTrigTowerDetId::kEBTotalTowers,0);

  // Record signal-full entries
  edm::Handle<EcalRecHitCollection> EBRecHitsH;
  iEvent.getByToken(EBRecHitCollectionT_, EBRecHitsH);
  for(EcalRecHitCollection::const_iterator iRHit = EBRecHitsH->begin();
      iRHit != EBRecHitsH->end();                      
      ++iRHit) {

    // Get detector id and convert to histogram-friendly coordinates
    EBDetId ebId( iRHit->id() );
    //EBDetId ebId( 1, 3 );
    //EcalTrigTowerDetId ttId( iRHit->id() );
    iphi_ = ebId.iphi()-1;
    ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta();
    //std::cout << "ECAL | (ieta,iphi): (" << ebId.ieta() << "," << ebId.iphi() << ")" <<std::endl;

    // Fill some histograms to monitor distributions
    // These will contain *cumulative* statistics and as such
    // such be used for monitoring purposes
    hEB_energy->Fill( iphi_,ieta_,iRHit->energy() );
    hEB_time->Fill( iphi_,ieta_,iRHit->time() );

    // Get Hashed Index
    // Hashed index provides a convenient index mapping
    // from [ieta][iphi] -> [idx]
    idx   = ebId.hashedIndex(); // (ieta_+EBDetId::MAX_IETA)*EBDetId::MAX_IPHI + iphi_
    //idx = ttId.hashedIndex();

    // Get global position of cell center
    //pos  = caloGeom->getPosition(ebId);
    //eta = pos.eta();
    //phi = pos.phi();
    //std::cout << "ECAL | (eta,phi,E): (" << eta << "," << phi << ","<< iRHit->energy()<<")" <<std::endl;

    // Fill event arrays
    // These are the actual inputs to the detector images
    vEB_energy_[idx] += iRHit->energy();
    vEB_time_[idx] += iRHit->time();
    //vEB_energy_[idx] += iRHit->energy()/TMath::CosH(cell->etaPos()); // pick out only transverse component

  } // EB reduced rechits

  // Write out state of histogram
  // For monitoring purposes only 
  char outFile[100];
  if (saveImgs) { 
    cEB->Clear();
    //hEB_energy->GetZaxis()->SetRangeUser(2.e-2, 9.e1);
    //hEB_energy->Draw("COL Z");
    sprintf(outFile,"cEB_energy_%llu.eps",iEvent.id().event());
    cEB->Print(outFile);
  }

  /*
  // ----- EB digis ----- //
  // This contains the raw EB digi collection:
  // Each digi is unpacked as a data frame containing info
  // from 10 time samples [iS] of the pulse shape
  // iS=[0-2]: Presample noise
  // iS=[3-9]: Nominal pulse shape
  // NOTE: This is the raw collection and includes
  // zero-supression and bad channel effects!

  // Initialize arrays
  for(int iS(0); iS < EcalDataFrame::MAXSAMPLES; ++iS)
    vEB_adc_[iS].assign(EBDetId::kSizeForDenseIndexing,0);

  // Record signal-full entries
  edm::Handle<EBDigiCollection> EBDigisH;
  iEvent.getByToken(EBDigiCollectionT_, EBDigisH);
  for(EBDigiCollection::const_iterator iDigi = EBDigisH->begin();
      iDigi != EBDigisH->end();                      
      ++iDigi) {

    // Get detector id and convert to histogram-friendly coordinates
    EBDetId ebId( iDigi->id() );
    //DetId id( iDigi->id() );
    iphi_ = ebId.iphi()-1;
    ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta();

    // Get Hashed Index & Cell Geometry
    // Hashed index provides a convenient index mapping
    // from [ieta][iphi] -> [idx]
    idx = ebId.hashedIndex(); // (ieta_+EBDetId::MAX_IETA)*EBDetId::MAX_IPHI + iphi_
    // Cell geometry provides access to (rho,eta,phi) coordinates of cell center
    //cell  = caloGeom->getGeometry(ebId);

    // Unpack the digi into a dataframe
    EcalDataFrame df(*iDigi);
    for(int iS(0); iS < EcalDataFrame::MAXSAMPLES; ++iS) {

      // Get the iS-th sample
      EcalMGPASample digiSample( df.sample(iS) );

      // Fill some histograms to monitor distributions
      // These will contain *cumulative* statistics and as such
      // such be used for monitoring purposes
      hEB_adc[iS]->Fill( iphi_, ieta_, digiSample.adc() );

      // Fill event arrays
      // These are the actual inputs to the detector images
      vEB_adc_[iS][idx] += digiSample.adc();
      //vEB_adc_[iS][idx] += digiSample.adc()/TMath::CosH(cell->etaPos()); // pick out only transverse component

    } // sample

  } // EB digi

  // Write out state of histogram
  // For monitoring purposes only 
  if (saveImgs) {
    for(int iS(0); iS < EcalDataFrame::MAXSAMPLES; ++iS) {
      cEB->Clear();
      gPad->SetLogz(1);
      //hEB_adc[iS]->GetZaxis()->SetRangeUser(190, 1.3e3);
      //hEB_adc[iS]->Draw("COL Z");
      sprintf(outFile,"cEB_adc%d_%llu.eps",iS,iEvent.id().event());
      cEB->Print(outFile);
    } // sample
  }
  */

  // ----- EE reduced rechit collection ----- //
  // This contatins the reduced EB rechit collection after
  // the zero suppression and bad channel clean-up

  int ix_, iy_, iz_; // NOTE: rows:iy, columns:ix

  // Initialize arrays
  for(int iz(0); iz < EE_IZ_MAX; ++iz) {
    vEE_energy_[iz].assign(EE_NC_PER_ZSIDE,0.);
    vEE_time_[iz].assign(EE_NC_PER_ZSIDE,0.);
  }

  // Record signal-full entries
  edm::Handle<EcalRecHitCollection> EERecHitsH;
  iEvent.getByToken(EERecHitCollectionT_, EERecHitsH);
  for(EcalRecHitCollection::const_iterator iRHit = EERecHitsH->begin();
      iRHit != EERecHitsH->end();                      
      ++iRHit) {

    // Get detector id and convert to histogram-friendly coordinates
    EEDetId eeId( iRHit->id() );
    ix_ = eeId.ix()-1;
    iy_ = eeId.iy()-1;
    iz_ = (eeId.zside() > 0) ? 1 : 0;
    //std::cout << "ECAL | (ix,iy): " << ix_ << "," << iy_ << "," << iz_ << "," << iRHit->energy() << std::endl;

    // Fill some histograms to monitor distributions
    // These will contain *cumulative* statistics and as such
    // such be used for monitoring purposes
    hEE_energy[iz_]->Fill( ix_,iy_,iRHit->energy() );
    hEE_time[iz_]->Fill( ix_,iy_,iRHit->time() );

    // Create hashed Index
    // Maps from [iy][ix] -> [idx]
    idx = iy_*EEDetId::IX_MAX + ix_; 

    // Fill event arrays
    // These are the actual inputs to the detector images
    vEE_energy_[iz_][idx] += iRHit->energy();
    vEE_time_[iz_][idx] += iRHit->time();
    //vEE_energy_[iz_][idx] += iRHit->energy()/TMath::CosH(cell->etaPos()); // pick out only transverse component

  } // EE reduced rechits

  //////////// HBHE //////////

  // ----- HBHE reduced rechit collection ----- //
  int depth_;

  // Initialize arrays
  vHBHE_energy_EB_.assign( hcaldqm::constants::IPHI_NUM*2*HBHE_IETA_MAX,0. );

  // Record signal-full entries
  edm::Handle<HBHERecHitCollection> HBHERecHitsH;
  iEvent.getByToken(HBHERecHitCollectionT_, HBHERecHitsH);
  for(HBHERecHitCollection::const_iterator iRHit = HBHERecHitsH->begin();
      iRHit != HBHERecHitsH->end();                      
      ++iRHit) {

    // Get detector id and convert to histogram-friendly coordinates
    // NOTE: HBHE detector ids are indexed by (ieta,iphi,depth)!
    HcalDetId hId( iRHit->id() );
    //HcalDetId hId( HcalSubdetector::HcalBarrel, 1, 71, 1 );
    //if (hId.subdet() != HcalSubdetector::HcalBarrel) continue;
    // WARNING: HBHE::iphi() is not aligned with EBRecHit::iphi()!
    // => Need to shift by 2 HBHE towers: HBHE::iphi: [1,...,71,72]->[3,4,...,71,72,1,2]
    iphi_  = hId.iphi()+2; // shift
    iphi_  = iphi_ > 72 ? iphi_-72 : iphi_; // wrap-around
    iphi_  = iphi_ -1; // make histogram-friendly
    //iphi_  = hId.iphi()-1;
    ieta_  = hId.ieta() > 0 ? hId.ieta()-1 : hId.ieta();
    depth_ = hId.depth();
    //std::cout << "HCAL | (ieta,iphi): (" << hId.ieta() << "," << iphi_+1 << ")" <<std::endl;

    // Fill some histograms to monitor distributions
    // These will contain *cumulative* statistics and as such
    // such be used for monitoring purposes
    hHBHE_depth->Fill( depth_ );

    // Restrict coverage of HBHE
    // HBHE_IETA_MAX == 17: match coverage of EB
    // HBHE_IETA_MAX == 20: match until granularity decreases
    if ( abs(hId.ieta()) > 20 ) {
      hHBHE_energy->Fill( iphi_  ,ieta_,iRHit->energy()*0.5 );
      hHBHE_energy->Fill( iphi_+1,ieta_,iRHit->energy()*0.5 );
      continue; 
    } else {
      hHBHE_energy->Fill( iphi_,ieta_,iRHit->energy() );
    }

    // Fill restricted coverage histograms
    hHBHE_energy_EB->Fill( iphi_,ieta_,iRHit->energy() );

    // Create hashed Index
    // Effectively sums energies over depth for a given (ieta,iphi)
    // Maps from [ieta][iphi] -> [idx]
    idx = ( ieta_+HBHE_IETA_MAX )*hcaldqm::constants::IPHI_NUM + iphi_;

    // Get global position of cell center
    //pos = caloGeom->getPosition(hId);
    //eta = pos.eta();
    //phi = pos.phi();
    //std::cout << "HCAL > (eta,phi,E): (" << eta << "," << phi << ","<< iRHit->energy()<<")" <<std::endl;

    // Fill event arrays
    // These are the actual inputs to the detector images
    vHBHE_energy_EB_[idx] += iRHit->energy();
    //vHBHE_energy_EB_[idx] += iRHit->energy()/TMath::CosH(cell->etaPos()); // pick out only transverse component

  }

  // Write out state of histogram
  // For monitoring purposes only 
  if (saveImgs) {
    cHBHE->cd();
    //std::cout << "maxEta: " << maxEta << std::endl;
    hHBHE_energy->GetZaxis()->SetRangeUser(2.e-2, 9.e1);
    hHBHE_energy->Draw("COL Z");
    sprintf(outFile,"cHBHE_energy_%llu.eps",iEvent.id().event());
    cHBHE->Print(outFile);
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
  // Write out event ID
  eventId_ = iEvent.id().event();
  m0_ = vDiPho.T();

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
