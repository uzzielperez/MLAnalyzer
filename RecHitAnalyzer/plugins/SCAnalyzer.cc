// -*- C++ -*-
//
// Package:    MLAnalyzer/SCAnalyzer
// Class:      SCAnalyzer
// 
/**\class SCAnalyzer SCAnalyzer.cc MLAnalyzer/SCAnalyzer/plugins/SCAnalyzer.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Michael Andrews
//         Created:  Mon, 17 Jul 2017 15:59:54 GMT
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

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h" // reco::PhotonCollection defined here

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
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

class SCAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
    explicit SCAnalyzer(const edm::ParameterSet&);
    ~SCAnalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    // ----------member data ---------------------------
    //edm::EDGetTokenT<edm::View<reco::GsfElectron>> electronCollectionT_;
    edm::EDGetTokenT<reco::GsfElectronCollection> electronCollectionT_;
    edm::EDGetTokenT<reco::PhotonCollection> photonCollectionT_;
    edm::EDGetTokenT<EcalRecHitCollection> EBRecHitCollectionT_;
    edm::EDGetTokenT<reco::GenParticleCollection> genParticleCollectionT_;

    static const int nPhotons = 1;
    static const int crop_size = 32;

    //TH1D * histo; 
    TH2D * hEB_energy; 
    TH2D * hSC_energy[nPhotons]; 
    TH2D * hSC_time[nPhotons]; 

    TTree* RHTree;

    float eventId_;
    std::vector<float> vEB_energy_;
    std::vector<float> vEB_SCenergy_;
    std::vector<float> vSC_energy_[nPhotons];
    std::vector<float> vSC_energyT_[nPhotons];
    std::vector<float> vSC_time_[nPhotons];
    std::vector<int> vGoodPhotonIdxs_;
    std::vector<float> vIphi_Emax;
    std::vector<float> vIeta_Emax;
    float vPho_pT_[nPhotons];
    float vPho_E_[nPhotons];
    float vPho_eta_[nPhotons];
    float vPho_phi_[nPhotons];

    // Initialize Calorimeter Geometry
    const CaloGeometry* caloGeom;
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
SCAnalyzer::SCAnalyzer(const edm::ParameterSet& iConfig)
{
  //EBRecHitCollectionT_ = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EBRecHitCollection"));
  //electronCollectionT_ = consumes<edm::View<reco::GsfElectron>>(iConfig.getParameter<edm::InputTag>("gsfElectronCollection"));
  electronCollectionT_ = consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("gsfElectronCollection"));
  photonCollectionT_ = consumes<reco::PhotonCollection>(iConfig.getParameter<edm::InputTag>("gedPhotonCollection"));
  EBRecHitCollectionT_ = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedEBRecHitCollection"));
  genParticleCollectionT_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticleCollection"));

  //now do what ever initialization is needed
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  // Histograms
  char hname[50], htitle[50];
  // EB rechits
  hEB_energy = fs->make<TH2D>("EB_energy", "E(i#phi,i#eta);i#phi;i#eta",
      EBDetId::MAX_IPHI  , EBDetId::MIN_IPHI-1, EBDetId::MAX_IPHI,
      2*EBDetId::MAX_IETA,-EBDetId::MAX_IETA,   EBDetId::MAX_IETA );
  // EB SCs
  for(int iPho (0); iPho < nPhotons; iPho++) {
    sprintf(hname, "SC_energy%d",iPho);
    sprintf(htitle,"E(i#phi,i#eta);iphi;ieta");
    hSC_energy[iPho] = fs->make<TH2D>(hname, htitle,
        crop_size, 0, crop_size,
        crop_size, 0, crop_size );
    sprintf(hname, "SC_time%d",iPho);
    sprintf(htitle,"t(i#phi,i#eta);iphi;ieta");
    hSC_time[iPho] = fs->make<TH2D>(hname, htitle,
        crop_size, 0, crop_size,
        crop_size, 0, crop_size );
  }
  // Output Tree
  RHTree = fs->make<TTree>("RHTree", "RecHit tree");
  RHTree->Branch("eventId",      &eventId_);
  RHTree->Branch("EB_energy",    &vEB_energy_);
  //RHTree->Branch("EB_SCenergy",  &vSC_energy_);
  for(int iPho (0); iPho < nPhotons; iPho++) {
    sprintf(hname, "SC_energy%d",iPho);
    RHTree->Branch(hname,          &vSC_energy_[iPho]);
    sprintf(hname, "SC_energyT%d",iPho);
    RHTree->Branch(hname,          &vSC_energyT_[iPho]);
    sprintf(hname, "SC_time%d",iPho);
    RHTree->Branch(hname,          &vSC_time_[iPho]);
    sprintf(hname, "pho_pT%d",iPho);
    RHTree->Branch(hname,      &vPho_pT_[iPho]);
    sprintf(hname, "pho_E%d",iPho);
    RHTree->Branch(hname,      &vPho_E_[iPho]);
    sprintf(hname, "pho_eta%d",iPho);
    RHTree->Branch(hname,      &vPho_eta_[iPho]);
    sprintf(hname, "pho_phi%d",iPho);
    RHTree->Branch(hname,      &vPho_phi_[iPho]);
  }
}


SCAnalyzer::~SCAnalyzer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
SCAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  edm::Handle<EcalRecHitCollection> EBRecHitsH;
  iEvent.getByToken(EBRecHitCollectionT_, EBRecHitsH);

  for(int iPho (0); iPho < nPhotons; iPho++) {
    vSC_energy_[iPho].assign(crop_size*crop_size,0.);
    vSC_energyT_[iPho].assign(crop_size*crop_size,0.);
    vSC_time_[iPho].assign(crop_size*crop_size,0.);
  }
  vEB_SCenergy_.assign(EBDetId::kSizeForDenseIndexing,0.);

  int iphi_, ieta_;

  vGoodPhotonIdxs_.clear();
  vIphi_Emax.clear();
  vIeta_Emax.clear();

  ////////// Apply selection and get coordinates of shower max //////////
  float Emax;
  int iphi_Emax, ieta_Emax;
  int nPho = 0;
  edm::Handle<reco::PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);
  std::cout << " >> PhoCol.size: " << photons->size() << std::endl;
  // Loop over photons
  for(reco::PhotonCollection::const_iterator iPho = photons->begin();
      iPho != photons->end();
      ++iPho) {

    // Apply reco selection
    std::cout << " >> pT: " << iPho->pt() << " eta: " << iPho->eta() << std::endl;
    if ( iPho->pt() < 10. ) continue;
    if ( iPho->eta() > 1.44 ) continue;

    // Get underlying super cluster
    reco::SuperClusterRef iSC = iPho->superCluster();
    //EcalRecHitCollection::const_iterator iRHit_( EBRecHitsH->find(iSC->seed()->seed()) );
    //std::cout << "Seed E: " << iRHit_->energy() << std::endl;
    std::vector<std::pair<DetId, float>> const& SCHits( iSC->hitsAndFractions() );
    std::cout << " >> SChits.size: " << SCHits.size() << std::endl;

    // Get Emax crystal
    Emax = 0.;
    iphi_Emax = -1;
    ieta_Emax = -1;

    // Loop over SC hits of photon
    for(unsigned iH(0); iH != SCHits.size(); ++iH) {

      // Get DetId
      if ( SCHits[iH].first.subdetId() != EcalBarrel ) continue;
      EcalRecHitCollection::const_iterator iRHit( EBRecHitsH->find(SCHits[iH].first) );
      if ( iRHit == EBRecHitsH->end() ) continue;

      // Convert coordinates to ordinals
      EBDetId ebId( iRHit->id() );
      ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta(); // [-85,...,-1,1,...,85]
      ieta_ += EBDetId::MAX_IETA; // [0,...,169]
      iphi_ = ebId.iphi()-1; // [0,...,359]

      // Keep coordinates of shower max
      if ( iRHit->energy() > Emax ) {
        Emax = iRHit->energy();
        iphi_Emax = iphi_;
        ieta_Emax = ieta_;
      }
      //std::cout << " >> " << iH << ": iphi_,ieta_,E: " << iphi_ << ", " << ieta_ << ", " << iRHit->energy() << std::endl; 
    } // SC hits

    // Apply selection on position of shower seed
    if ( ieta_Emax > 169 - 15 || ieta_Emax < 15 ) continue;
    vGoodPhotonIdxs_.push_back( nPho );
    vIphi_Emax.push_back( iphi_Emax );
    vIeta_Emax.push_back( ieta_Emax );
    std::cout << " >> Found: iphi_Emax,ieta_Emax: " << iphi_Emax << ", " << ieta_Emax << std::endl;
    nPho++;

  } // Photons

  // Enforce selection
  std::cout << " >> nPho: " << nPho << std::endl;
  if ( nPho != nPhotons ) return;
  std::cout << " >> Passed selection. " << std::endl;
    
  ////////// Store each shower crop //////////
  int i = 0;
  int iP = 0;
  for(reco::PhotonCollection::const_iterator iPho = photons->begin();
      iPho != photons->end();
      ++iPho) {

    // Skip failed photons
    if ( std::find(vGoodPhotonIdxs_.begin(), vGoodPhotonIdxs_.end(), i) == vGoodPhotonIdxs_.end() ) {
      i++;
      continue;
    }

    // Fill branch arrays
    vPho_pT_[iP] = iPho->pt();
    vPho_E_[iP] = iPho->energy();
    vPho_eta_[iP] = iPho->eta();
    vPho_phi_[iP] = iPho->phi();
    i++;
    iP++;

  } // photons

  int idx;
  int iphi_shift, ieta_shift;
  int iphi_crop, ieta_crop;
  // Provides access to global cell position and coordinates below
  edm::ESHandle<CaloGeometry> caloGeomH;
  iSetup.get<CaloGeometryRecord>().get(caloGeomH);
  caloGeom = caloGeomH.product();
  for(EcalRecHitCollection::const_iterator iRHit = EBRecHitsH->begin();
      iRHit != EBRecHitsH->end();
      ++iRHit) {

    // Convert detector coordinates to ordinals
    EBDetId ebId( iRHit->id() );
    iphi_ = ebId.iphi()-1; // [0,...,359]
    ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta(); // [-85,...,-1,0,...,84]
    ieta_ += EBDetId::MAX_IETA; // [0,...,169]

    for(unsigned iP(0); iP < nPhotons; iP++) {

      iphi_shift = vIphi_Emax[iP] - 15;
      ieta_shift = vIeta_Emax[iP] - 15;
      //std::cout << " >> Storing: iphi_Emax,ieta_Emax: " << vIphi_Emax[iP] << ", " << vIeta_Emax[iP] << std::endl;

      // Convert to [0,...,31][0,...,31]
      ieta_crop = ieta_ - ieta_shift;
      iphi_crop = iphi_ - iphi_shift;
      if ( iphi_crop >= EBDetId::MAX_IPHI ) iphi_crop = iphi_crop - EBDetId::MAX_IPHI; // get wrap-around hits

      if ( ieta_crop < 0 || ieta_crop > crop_size-1 ) continue;
      if ( iphi_crop < 0 || iphi_crop > crop_size-1 ) continue;
      
      // Convert to [0,...,32*32-1] 
      idx = ieta_crop*crop_size + iphi_crop;

      // Cell geometry provides access to (rho,eta,phi) coordinates of cell center
      auto cell = caloGeom->getGeometry(ebId);
      
      // Fill branch arrays 
      vSC_energy_[iP][idx] = iRHit->energy();
      vSC_energyT_[iP][idx] = iRHit->energy()/TMath::CosH(cell->etaPos());
      vSC_time_[iP][idx] = iRHit->time();
      vEB_SCenergy_[ebId.hashedIndex()] = iRHit->energy();
      //std::cout << " >> " << iP << ": iphi_,ieta_,E: " << iphi_crop << ", " << ieta_crop << ", " << iRHit->energy() << std::endl; 

      // Fill histograms to monitor cumulative distributions
      hSC_energy[iP]->Fill( iphi_crop,ieta_crop,iRHit->energy() );
      hSC_time[iP]->Fill( iphi_crop,ieta_crop,iRHit->time() );

    } // EB rechits
  } // photons

    /*
    // Get underlying super cluster
    reco::SuperClusterRef iSC = iPho->superCluster();
    //EcalRecHitCollection::const_iterator iRHit_( EBRecHitsH->find(iSC->seed()->seed()) );
    //std::cout << "Seed E: " << iRHit_->energy() << std::endl;
    std::vector<std::pair<DetId, float>> const& SCHits( iSC->hitsAndFractions() );

    // Fill array
    for(unsigned iH(0); iH != SCHits.size(); ++iH) {

      if ( SCHits[iH].first.subdetId() != EcalBarrel ) continue;
      EcalRecHitCollection::const_iterator iRHit( EBRecHitsH->find(SCHits[iH].first) );
      if ( iRHit == EBRecHitsH->end() ) continue;

      // Convert detector coordinates to ordinals
      EBDetId ebId( iRHit->id() );
      iphi_ = ebId.iphi()-1; // [0,...,359]
      ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta(); // [-85,...,-1,0,...,84]
      ieta_ += EBDetId::MAX_IETA; // [0,...,169]

      // Convert to [0,...,31][0,...,31]
      ieta_ = ieta_ - ieta_shift;
      iphi_ = iphi_ - iphi_shift;
      if ( iphi_ >= EBDetId::MAX_IPHI ) iphi_ = iphi_ - EBDetId::MAX_IPHI; // get wrap-around hits
      
      // Convert to [0,...,32*32-1] 
      idx = ieta_*crop_size + iphi_;

      // Fill branch arrays 
      vSC_energy_[iP][idx] = iRHit->energy();
      vSC_time_[iP][idx] = iRHit->time();
      vEB_SCenergy_[ebId.hashedIndex()] = iRHit->energy();
      //std::cout << " >> " << iH << ": iphi_,ieta_,E: " << iphi_ << ", " << ieta_ << ", " << iRHit->energy() << std::endl; 
    } // SC hits

    // Fill histograms to monitor cumulative distributions
    for(unsigned iD(0); iD < crop_size*crop_size; ++iD) {
      ieta_ = std::floor(1.*iD/crop_size);
      iphi_ = iD % crop_size;
      //std::cout << " >> " << iD << ": iphi_,ieta_,E: " << iphi_ << ", " << ieta_ << ", " << vSC_energy_[iP][iD] << std::endl; 
      hEB_SCenergy[iP]->Fill( iphi_,ieta_,vSC_energy_[iP][iD] );
      hEB_SCtime[iP]->Fill( iphi_,ieta_,vSC_time_[iP][iD] );
    }

    i++;
    iP++;
  } // photons
  */

  // Fill full EB for comparison
  vEB_energy_.assign(EBDetId::kSizeForDenseIndexing,0.);
  for(EcalRecHitCollection::const_iterator iRHit = EBRecHitsH->begin();
      iRHit != EBRecHitsH->end();
      ++iRHit) {

    // Get detector id and convert to histogram-friendly coordinates
    EBDetId ebId( iRHit->id() );
    iphi_ = ebId.iphi()-1;
    ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta();
    //std::cout << "ECAL | (ieta,iphi): (" << ebId.ieta() << "," << ebId.iphi() << ")" <<std::endl;

    // Fill some histograms to monitor distributions
    // These will contain *cumulative* statistics and as such
    // should be used for monitoring purposes only
    hEB_energy->Fill( iphi_,ieta_,iRHit->energy() );

    // Fill branch arrays
    idx = ebId.hashedIndex(); // (ieta_+EBDetId::MAX_IETA)*EBDetId::MAX_IPHI + iphi_
    vEB_energy_[idx] = iRHit->energy();
  } // EB rechits

  /*
  //edm::Handle<edm::View<reco::GsfElectron>> electrons;
  edm::Handle<reco::GsfElectronCollection> electrons;
  iEvent.getByToken(electronCollectionT_, electrons);
  //for(edm::View<reco::GsfElectron>::const_iterator iEle = electrons->begin();
  std::cout << "EleCol.size: " << electrons->size() << std::endl;
  int nEle = 0;
  for(reco::GsfElectronCollection::const_iterator iEle = electrons->begin();
      iEle != electrons->end();
      ++iEle) {

    reco::SuperClusterRef iSC = iEle->superCluster();
    std::vector<std::pair<DetId, float>> const& SCHits( iSC->hitsAndFractions() );
    std::cout << "SChits.size: " << SCHits.size() << std::endl;
    for(unsigned iH(0); iH != SCHits.size(); ++iH) {
      if ( SCHits[iH].first.subdetId() != EcalBarrel ) continue;
      EcalRecHitCollection::const_iterator iRHit( EBRecHitsH->find(SCHits[iH].first) );
      if ( iRHit == EBRecHitsH->end() ) continue;
      EBDetId ebId( iRHit->id() );
      iphi_ = ebId.iphi()-1;
      ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta();
      //std::cout << "iphi_,ieta_:" << iphi_ << "," << ieta_ << std::endl; 
      //if ( nEle == 0 ) hEB_energy->Fill( iphi_,ieta_,iRHit->energy() );
    } // SC hits
    nEle++;
  } // electrons
  */
  eventId_ = iEvent.id().event();

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
SCAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SCAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SCAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(SCAnalyzer);
