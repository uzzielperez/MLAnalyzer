// -*- C++ -*-
//
// Package:    MLAnalyzer/SCRegressor
// Class:      SCRegressor
// 
/**\class SCRegressor SCRegressor.cc MLAnalyzer/SCRegressor/plugins/SCRegressor.cc

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
#include "TH3.h"
#include "TTree.h"
#include "TStyle.h"
#include "TMath.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class SCRegressor : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
    explicit SCRegressor(const edm::ParameterSet&);
    ~SCRegressor();

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
    edm::EDGetTokenT<reco::GenJetCollection> genJetCollectionT_;

    static const int nPhotons = 2;
    //static const int nPhotons = 1;
    static const int crop_size = 32;
    //static const bool debug = true;
    static const bool debug = false;

    //TH1D * histo; 
    TH2D * hEB_energy; 
    TH2D * hSC_energy[nPhotons]; 
    TH2D * hSC_time[nPhotons]; 
    TH1F * hSC_mass; 
    TH1F * hDR; 
    TH1F * hdEta; 
    TH1F * hdPhi; 
    TH3F * hdPhidEta; 
    TH1F * hPt; 

    TTree* RHTree;

    float eventId_;
    std::vector<float> vEB_energy_;
    std::vector<float> vEB_SCenergy_;
    std::vector<float> vSC_energy_[nPhotons];
    std::vector<float> vSC_energyT_[nPhotons];
    std::vector<float> vSC_energyZ_[nPhotons];
    std::vector<float> vSC_time_[nPhotons];
    float vPho_pT_[nPhotons];
    float vPho_E_[nPhotons];
    float vPho_eta_[nPhotons];
    float vPho_phi_[nPhotons];
    float vPho_r9_[nPhotons];
    float vSC_mass_[nPhotons];
    float vSC_DR_[nPhotons];
    float vSC_pT_[nPhotons];

    float getEMJetMass( reco::GenJetRef, const GlobalPoint& );
    int nTotal, nPassed;

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
SCRegressor::SCRegressor(const edm::ParameterSet& iConfig)
{
  //EBRecHitCollectionT_ = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EBRecHitCollection"));
  //electronCollectionT_ = consumes<edm::View<reco::GsfElectron>>(iConfig.getParameter<edm::InputTag>("gsfElectronCollection"));
  electronCollectionT_ = consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("gsfElectronCollection"));
  photonCollectionT_ = consumes<reco::PhotonCollection>(iConfig.getParameter<edm::InputTag>("gedPhotonCollection"));
  EBRecHitCollectionT_ = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedEBRecHitCollection"));
  genParticleCollectionT_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticleCollection"));
  genJetCollectionT_ = consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJetCollection"));

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
  hSC_mass = fs->make<TH1F>("SC_mass", "m_{SC};m_{SC}",50, 0., 0.5);
  hDR = fs->make<TH1F>("DR_seed_subJet", "#DeltaR(seed,subJet);#DeltaR",50, 0., 50.*0.0174);
  hdEta = fs->make<TH1F>("dEta_seed_subJet", "#Delta#eta(seed,subJet);#Delta#eta",50, 0., 50.*0.0174);
  hdPhi = fs->make<TH1F>("dPhi_seed_subJet", "#Delta#phi(seed,subJet);#Delta#phi",50, 0., 50.*0.0174);
  hdPhidEta = fs->make<TH3F>("dPhidEta_GG", "#Delta(#phi,#eta,m);#Delta#phi(#gamma,#gamma);#Delta#eta(#gamma,#gamma);m",6, 0., 6.*0.0174, 6, 0., 6.*0.0174, 16., 0.,1.6);
  hPt = fs->make<TH1F>("Pt", "Pt", 65, 30., 160.);
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
    sprintf(hname, "SC_energyZ%d",iPho);
    RHTree->Branch(hname,          &vSC_energyZ_[iPho]);
    sprintf(hname, "SC_time%d",iPho);
    RHTree->Branch(hname,          &vSC_time_[iPho]);
    sprintf(hname, "SC_mass%d",iPho);
    RHTree->Branch(hname,          &vSC_mass_[iPho]);
    sprintf(hname, "SC_DR%d",iPho);
    RHTree->Branch(hname,          &vSC_DR_[iPho]);
    sprintf(hname, "SC_pT%d",iPho);
    RHTree->Branch(hname,          &vSC_pT_[iPho]);
    sprintf(hname, "pho_pT%d",iPho);
    RHTree->Branch(hname,      &vPho_pT_[iPho]);
    sprintf(hname, "pho_E%d",iPho);
    RHTree->Branch(hname,      &vPho_E_[iPho]);
    sprintf(hname, "pho_eta%d",iPho);
    RHTree->Branch(hname,      &vPho_eta_[iPho]);
    sprintf(hname, "pho_phi%d",iPho);
    RHTree->Branch(hname,      &vPho_phi_[iPho]);
    sprintf(hname, "pho_r9%d",iPho);
    RHTree->Branch(hname,      &vPho_r9_[iPho]);
  }
}


SCRegressor::~SCRegressor()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//
const math::XYZTLorentzVector getFinalP4( const reco::Candidate* phoCand ) {
  if ( phoCand->status() == 1 ) {
    return phoCand->p4();
  } else {
    if ( phoCand->numberOfDaughters() == 1 )
      return getFinalP4( phoCand->daughter(0) );
    else
      return phoCand->p4();
  }
}

float SCRegressor::getEMJetMass( reco::GenJetRef iJet, const GlobalPoint & seed )
{
    unsigned int nConstituents = iJet->getGenConstituents().size();
    //const std::vector <const reco::GenParticle*> jetConstituents = iJet->getGenConstituents();
    if ( debug ) std::cout << " >> nConst:" << nConstituents << " eta:" << iJet->eta() << " phi:" << iJet->phi() << std::endl;
    if ( debug ) std::cout << " >> seed eta:" << seed.eta() << " phi:" << seed.phi() << std::endl;

    int nSelConst = 0;
    float dR, dEta, dPhi;
    //math::PtEtaPhiELorentzVectorD jetP4;
    math::XYZTLorentzVector jetP4;
    for ( unsigned int j = 0; j < nConstituents; j++ ) {
      const reco::GenParticle* subJet = iJet->getGenConstituent( j );
      if ( debug ) {
        std::cout << " >> pdgId:" << subJet->pdgId() << " status:" << subJet->status()
        << " pT:" << subJet->pt() << " eta:" << subJet->eta() << " phi: " << subJet->phi() << " E:" << subJet->energy();
        std::cout << " moth:" << subJet->mother()->pdgId();
        std::cout << std::endl;
      }
      if ( std::abs(subJet->pdgId()) != 22 && std::abs(subJet->pdgId()) != 11 ) continue;
      if ( subJet->status() != 1 ) continue;
      if ( subJet->mother()->pdgId() != 111 && subJet->mother()->pdgId() != 35 ) continue;
      //dR = reco::deltaR( seed.eta(),seed.phi(), subJet->eta(),subJet->phi() );
      //dPhi = reco::deltaPhi( seed.phi(), subJet->phi() );
      //dEta = std::abs( seed.eta() - subJet->eta() );
      dR = reco::deltaR( iJet->eta(),iJet->phi(), subJet->eta(),subJet->phi() );
      dPhi = reco::deltaPhi( iJet->phi(), subJet->phi() );
      dEta = std::abs( iJet->eta() - subJet->eta() );
      //if ( dPhi > TMath::Pi() ) dPhi =- TMath::Pi();
      hDR->Fill( dR );
      hdEta->Fill( dEta );
      hdPhi->Fill( dPhi );
      if ( debug ) std::cout << " >> dR:" << dR << " dEta:" << dEta << " dPhi:" << dPhi << std::endl;
      if ( dR > (6*0.0174) ) continue;
      //if ( dR > (4*0.0174) ) continue;
      jetP4 += subJet->p4();
      nSelConst++;
    }
    if ( nSelConst == 0 ) return 0.;
    if ( jetP4.mass() < 1.e-3 ) return 0.;
    //std::cout << " EMJetMass:" << jetP4.mass() << std::endl;
    //std::cout << " >> pT:" << iJet->pt() << " eta:" << iJet->eta() << " phi: " << iJet->phi() << " E:" << iJet->energy() << std::endl;
    return jetP4.mass();
}

// ------------ method called for each event  ------------
void
SCRegressor::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  nTotal++;
  using namespace edm;

  edm::Handle<EcalRecHitCollection> EBRecHitsH;
  iEvent.getByToken(EBRecHitCollectionT_, EBRecHitsH);
  edm::Handle<reco::PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleCollectionT_, genParticles);
  // Provides access to global cell position and coordinates below
  edm::ESHandle<CaloGeometry> caloGeomH;
  iSetup.get<CaloGeometryRecord>().get(caloGeomH);
  const CaloGeometry* caloGeom = caloGeomH.product();
  if ( debug ) std::cout << " >> PhoCol.size: " << photons->size() << std::endl;

  for(int iPho (0); iPho < nPhotons; iPho++) {
    vSC_energy_[iPho].assign(crop_size*crop_size,0.);
    vSC_energyT_[iPho].assign(crop_size*crop_size,0.);
    vSC_energyZ_[iPho].assign(crop_size*crop_size,0.);
    vSC_time_[iPho].assign(crop_size*crop_size,0.);
  }
  vEB_SCenergy_.assign(EBDetId::kSizeForDenseIndexing,0.);

  int iphi_, ieta_;

  std::vector<int> vGenPi0Idxs_;
  std::vector<int> vGoodPhotonIdxs_;
  std::vector<float> vIphi_Emax;
  std::vector<float> vIeta_Emax;
  float dEta, dPhi, dR, mPi0;

  // ____________ Gen-level studies ________________ //
  int nPi0 = 0;
  int idx = -1;
  //std::vector<math::PtEtaPhiELorentzVectorD> vPhoPairs[nPhotons];
  std::vector<math::XYZTLorentzVector> vPhoPairs[nPhotons];
  for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin();
       iGen != genParticles->end();
       ++iGen) {

    idx++;
    // ID cuts
    if ( debug ) std::cout << " >> pdgId:"<< iGen->pdgId() << " status:" << iGen->status() << " nDaughters:" << iGen->numberOfDaughters() << std::endl;
    if ( std::abs(iGen->pdgId()) != 111 ) continue;
    if ( debug ) std::cout << " >> pdgId:111 nDaughters:" << iGen->numberOfDaughters() << std::endl;
    if ( iGen->numberOfDaughters() != 2 ) continue;
    //if ( iGen->mass() > 0.4 ) continue;
    vGenPi0Idxs_.push_back( idx );
    for ( unsigned int iD = 0; iD < iGen->numberOfDaughters(); iD++ ) {
      const reco::Candidate* phoCand = iGen->daughter(iD);
      vPhoPairs[nPi0].push_back( getFinalP4(phoCand) );
    }
    nPi0++;

  } // genParticle loop: count good photons
  if ( nPi0 != nPhotons ) return;

  ////////// Apply selection and get coordinates of shower max //////////
  float Emax;
  int iphi_Emax, ieta_Emax;
  int nPho = 0;
  GlobalPoint pos_Emax;
  float ptCut = 45., etaCut = 1.44;
  std::vector<GlobalPoint> vPos_Emax;
  bool isGenMatched;
  // Loop over photons
  idx = -1;
  for(reco::PhotonCollection::const_iterator iPho = photons->begin();
      iPho != photons->end();
      ++iPho) {

    idx++;

    // Apply reco selection
    if ( debug ) std::cout << " >> pT: " << iPho->pt() << " eta: " << iPho->eta() << std::endl;
    if ( iPho->pt() < ptCut ) continue;
    if ( std::abs(iPho->eta()) > etaCut ) continue;
    if ( debug ) std::cout << " >> pT: " << iPho->pt() << " eta: " << iPho->eta() << std::endl;
    if ( iPho->r9() < 0.7 ) continue;

    isGenMatched = false;
    for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin();
         iGen != genParticles->end();
         ++iGen) {

      // ID cuts
      if ( std::abs(iGen->pdgId()) != 22 ) continue;
      if ( iGen->status() != 1 ) continue; // NOT the same as Pythia status
      if ( !iGen->mother() ) continue;
      //if ( iGen->mother()->pdgId() != 35 && iGen->mother()->pdgId() != 22 ) continue;
      //if ( iGen->mother()->pdgId() != 35 ) continue;
      if ( iGen->mother()->pdgId() != 111 ) continue;
      // Kinematic cuts
      if ( std::abs(iGen->eta()) > etaCut ) continue;
      if ( std::abs(iGen->pt()) < ptCut ) continue;
      /*
      std::cout << " >> pdgId:" << iGen->pdgId() << " status:" << iGen->status()
        << " pT:" << iGen->pt() << " eta:" << iGen->eta() << " phi: " << iGen->phi() << " E:" << iGen->energy();
        std::cout << " moth:" << iGen->mother()->pdgId();
        std::cout << std::endl;
      */
      dR = reco::deltaR( iPho->eta(),iPho->phi(), iGen->eta(),iGen->phi() );
      if ( dR > 0.04 ) continue;
      isGenMatched = true;
      if ( debug ) std::cout << " >> GenPhoMatch | " << " pt:" << iGen->pt() << " dR: " << dR << std::endl;

    } // genParticle loop: count good photons

    //if ( !isGenMatched ) continue;
    if ( isGenMatched || !isGenMatched ) ;

    // Get underlying super cluster
    reco::SuperClusterRef const& iSC = iPho->superCluster();
    //EcalRecHitCollection::const_iterator iRHit_( EBRecHitsH->find(iSC->seed()->seed()) );
    //std::cout << "Seed E: " << iRHit_->energy() << std::endl;
    std::vector<std::pair<DetId, float>> const& SCHits( iSC->hitsAndFractions() );
    //std::cout << " >> SChits.size: " << SCHits.size() << std::endl;

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
        pos_Emax = caloGeom->getPosition(ebId);
      }
      //std::cout << " >> " << iH << ": iphi_,ieta_,E: " << iphi_ << ", " << ieta_ << ", " << iRHit->energy() << std::endl; 
    } // SC hits

    // Apply selection on position of shower seed
    //std::cout << " >> Found: iphi_Emax,ieta_Emax: " << iphi_Emax << ", " << ieta_Emax << std::endl;
    if ( Emax == 0. ) continue;
    if ( ieta_Emax > 169 - 15 || ieta_Emax < 15 ) continue;
    vGoodPhotonIdxs_.push_back( idx );
    vIphi_Emax.push_back( iphi_Emax );
    vIeta_Emax.push_back( ieta_Emax );
    vPos_Emax.push_back( pos_Emax );
    //std::cout << " >> Found: iphi_Emax,ieta_Emax: " << iphi_Emax << ", " << ieta_Emax << std::endl;
    nPho++;

  } // Photons

  // Enforce selection
  //std::cout << " >> nPho: " << nPho << std::endl;
  if ( nPho != nPhotons ) return;
  std::cout << " >> Passed selection. " << std::endl;

  /*
  // Get inv mass of SC
  // photons: massless
  // jet: sum all genJet constituents
  // Check if genJet
  edm::Handle<reco::GenJetCollection> genJets;
  iEvent.getByToken(genJetCollectionT_, genJets);
  std::vector<float> vJetFakePhoIdxs;
  //std::vector<float> vSC_mass_;
  //vSC_mass_.assign( nPhotons, 0. );
  float minDR;
  int minDR_idx;
  for ( unsigned int i = 0; i < nPhotons; i++ ) {
    reco::PhotonRef iPho( photons, vGoodPhotonIdxs_[i] );
    minDR = 100.;
    minDR_idx = -1;
    for ( unsigned int j = 0; j < genJets->size(); j++ ) {
      if ( std::find(vJetFakePhoIdxs.begin(),vJetFakePhoIdxs.end(),j) != vJetFakePhoIdxs.end() ) continue;
      reco::GenJetRef iJet( genJets, j );
      dR = reco::deltaR( iPho->eta(),iPho->phi(), iJet->eta(),iJet->phi() );
      if ( dR < minDR ) {
        minDR = dR;
        minDR_idx = j;
      }
    } // genJets
    if ( debug ) std::cout << " >> JetDR:" << minDR << std::endl;
    if ( minDR > 0.04 ) continue;

    reco::GenJetRef iJet_( genJets, minDR_idx );
    if ( debug ) std::cout << " >> JetDRmatched:" << minDR << " | iJet:" << minDR_idx << " pt:" << iJet_->pt() << std::endl;
    if ( debug ) std::cout << " >> EMjet mass:" << getEMJetMass( iJet_, vPos_Emax[i] ) << std::endl;
    //vSC_mass_[i] = getEMJetMass( iJet_, vPos_Emax[i] );
    //if ( vSC_mass_[i] < 0.1 ) std::cout << "zero mass " << vSC_mass_[i] << std::endl;
    //vSC_mass_[i] = 0.1;
    vJetFakePhoIdxs.push_back( minDR_idx );
  } // good photons
  if ( debug ) std::cout << "nDRmatchedJets: " << vJetFakePhoIdxs.size() << std::endl;
  if ( vJetFakePhoIdxs.size() != nPhotons ) return;
  //if ( vJetFakePhoIdxs.size() != 1 ) return;
  */

  ////////// Store each shower crop //////////
  for ( unsigned int i = 0; i < nPhotons; i++ ) {
    reco::PhotonRef iPho( photons, vGoodPhotonIdxs_[i] );
    // Fill branch arrays
    vPho_pT_[i] = iPho->pt();
    vPho_E_[i] = iPho->energy();
    vPho_eta_[i] = iPho->eta();
    vPho_phi_[i] = iPho->phi();
    //std::cout << "r9: " << iPho->r9() << std::endl;
    vPho_r9_[i] = iPho->r9();
  } // photons
  for ( int i = 0; i < nPi0; i++ ) {
    mPi0 = (vPhoPairs[i][0] + vPhoPairs[i][1]).mass();
    dR = reco::deltaR( vPhoPairs[i][0].eta(),vPhoPairs[i][0].phi(), vPhoPairs[i][1].eta(),vPhoPairs[i][1].phi() );
    dEta = std::abs( vPhoPairs[i][0].eta() - vPhoPairs[i][1].eta() );
    dPhi = reco::deltaPhi( vPhoPairs[i][0].phi(), vPhoPairs[i][1].phi() );
    if ( debug ) std::cout << " >> m0:" << mPi0 << " dR:" << dR << " dPhi:" << dPhi << std::endl;
    reco::GenParticleRef iGen( genParticles, vGenPi0Idxs_[i] );
    vSC_DR_[i] = dR;
    vSC_pT_[i] = iGen->pt(); 
    vSC_mass_[i] = mPi0;
    hPt->Fill( iGen->pt() );
    hdPhidEta->Fill( dPhi, dEta, mPi0 );
  }
  for ( int i = 0; i < nPhotons; i++ ) {
    //std::cout << "SC mass " << vSC_mass_[i] << std::endl;
    hSC_mass->Fill( vSC_mass_[i] );
  }

  //int idx;
  int iphi_shift, ieta_shift;
  int iphi_crop, ieta_crop;
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
      //auto cell = caloGeom->getGeometry(ebId);
      auto pos = caloGeom->getPosition(ebId);
      
      // Fill branch arrays 
      vSC_energy_[iP][idx] = iRHit->energy();
      //vSC_energyT_[iP][idx] = iRHit->energy()/TMath::CosH(vPho_eta_[iP]);
      vSC_energyT_[iP][idx] = iRHit->energy()/TMath::CosH(pos.eta());
      vSC_energyZ_[iP][idx] = iRHit->energy()*std::abs(TMath::TanH(pos.eta()));
      vSC_time_[iP][idx] = iRHit->time();
      vEB_SCenergy_[ebId.hashedIndex()] = iRHit->energy();
      //std::cout << " >> " << iP << ": iphi_,ieta_,E: " << iphi_crop << ", " << ieta_crop << ", " << iRHit->energy() << std::endl; 

      // Fill histograms to monitor cumulative distributions
      hSC_energy[iP]->Fill( iphi_crop,ieta_crop,iRHit->energy() );
      hSC_time[iP]->Fill( iphi_crop,ieta_crop,iRHit->time() );

    } // EB rechits
  } // photons

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

  eventId_ = iEvent.id().event();
  nPassed++;

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
SCRegressor::beginJob()
{
  nTotal = 0;
  nPassed = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SCRegressor::endJob() 
{
  std::cout << " selected: " << nPassed << "/" << nTotal << std::endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SCRegressor::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(SCRegressor);
