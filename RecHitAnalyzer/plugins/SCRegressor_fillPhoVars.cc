#include "MLAnalyzer/RecHitAnalyzer/interface/SCRegressor.h"

// Initialize branches _____________________________________________________//
void SCRegressor::branchesPhoVars ( TTree* tree, edm::Service<TFileService> &fs )
{

  tree->Branch("pho_pT",    &vPho_pT_);
  tree->Branch("pho_E",     &vPho_E_);
  tree->Branch("pho_eta",   &vPho_eta_);
  tree->Branch("pho_phi",   &vPho_phi_);

  tree->Branch("pho_r9",             &vPho_r9_);
  tree->Branch("pho_sieie",          &vPho_sieie_);
  tree->Branch("pho_phoIso",         &vPho_phoIso_);
  tree->Branch("pho_chgIso",         &vPho_chgIso_);
  tree->Branch("pho_chgIsoWrongVtx", &vPho_chgIsoWrongVtx_);
  tree->Branch("pho_Eraw",           &vPho_Eraw_);
  tree->Branch("pho_phiWidth",       &vPho_phiWidth_);
  tree->Branch("pho_etaWidth",       &vPho_etaWidth_);
  tree->Branch("pho_scEta",          &vPho_scEta_);
  tree->Branch("pho_sieip",          &vPho_sieip_);
  tree->Branch("pho_s4",             &vPho_s4_);

} // branchesPhoVars()

// Fill PhoVars rechits _________________________________________________________________//
void SCRegressor::fillPhoVars ( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{

  EcalClusterLazyTools clusterTools ( iEvent, iSetup, EBRecHitCollectionT_, EERecHitCollectionT_, ESRecHitCollectionT_);

  /*
  edm::Handle<double> rhoH;
  iEvent.getByToken( rhoLabel_, rhoH );
  //std::cout << "rho:" << *rhoH << std::endl;
  //std::cout << "rho:" << *(rhoH.product()) << std::endl;
  */

  edm::Handle<PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);
  /*
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleCollectionT_, genParticles);
  for ( auto const& mG : mGenPi0_RecoPho ) { 
    reco::GenParticleRef iGen( genParticles, mG.first );
    std::cout << "mass:" << iGen->mass() << std::endl;
    std::cout << "ptgen:" << iGen->pt() << std::endl;
  }
  */

  ////////// Store kinematics //////////

  vPho_pT_.clear();
  vPho_E_.clear();
  vPho_eta_.clear();
  vPho_phi_.clear();
  vPho_r9_.clear();
  vPho_sieie_.clear();
  for ( int iP : vRegressPhoIdxs_ ) {

    PhotonRef iPho( photons, iP );
    // Fill branch arrays
    vPho_pT_.push_back( iPho->pt() );
    vPho_E_.push_back( iPho->energy() );
    vPho_eta_.push_back( iPho->eta() );
    vPho_phi_.push_back( iPho->phi() );
  } // photons
  
  vPho_r9_.clear(); 
  vPho_sieie_.clear(); 
  vPho_phoIso_.clear(); 
  vPho_chgIso_.clear(); 
  vPho_chgIsoWrongVtx_.clear(); 
  vPho_Eraw_.clear(); 
  vPho_phiWidth_.clear(); 
  vPho_etaWidth_.clear(); 
  vPho_scEta_.clear(); 
  vPho_sieip_.clear(); 
  vPho_s4_.clear(); 

  //for ( auto const& mG : mGenPi0_RecoPho ) { 
  for ( int iP : vRegressPhoIdxs_ ) {

    //PhotonRef iPho( photons, mG.second[0] );
    PhotonRef iPho( photons, iP );
    reco::SuperClusterRef const& iSC = iPho->superCluster();
    std::vector<float> vCov = clusterTools.localCovariances( *(iSC->seed()) );

    vPho_r9_.push_back(             iPho->full5x5_r9() );
    vPho_sieie_.push_back(          iPho->full5x5_sigmaIetaIeta() );
    vPho_phoIso_.push_back(         iPho->photonIso() );
    vPho_chgIso_.push_back(         iPho->chargedHadronIso() );
    vPho_chgIsoWrongVtx_.push_back( iPho->chargedHadronIsoWrongVtx() );
    vPho_Eraw_.push_back(           iSC->rawEnergy() );
    vPho_phiWidth_.push_back(       iSC->phiWidth() );
    vPho_etaWidth_.push_back(       iSC->etaWidth() );
    vPho_scEta_.push_back(          iSC->eta() );
    vPho_sieip_.push_back(          vCov[1] );
    vPho_s4_.push_back(             clusterTools.e2x2( *(iSC->seed()) ) / clusterTools.e5x5( *(iSC->seed()) ) );
    /*
    std::cout << "pt:" << iPho->pt() << std::endl;
    std::cout << "E:" << iPho->energy() << std::endl;

    std::cout << "r9:"             << iPho->full5x5_r9() << std::endl;
    std::cout << "sieie:"          << iPho->full5x5_sigmaIetaIeta() << std::endl;
    std::cout << "phoIso:"         << iPho->photonIso() << std::endl;
    std::cout << "chgIso:"         << iPho->chargedHadronIso() << std::endl;
    std::cout << "chgIsoWrongVtx:" << iPho->chargedHadronIsoWrongVtx() << std::endl;

    std::cout << "Eraw:"           << iSC->rawEnergy() << std::endl;
    std::cout << "phiwidth:"       << iSC->phiWidth() << std::endl;
    std::cout << "etawidth:"       << iSC->etaWidth() << std::endl;
    std::cout << "eta:"            << iSC->eta() << std::endl;

    std::cout << "sieip:"          << vCov[1] << std::endl;
    std::cout << "s4:"             << clusterTools.e2x2( *(iSC->seed()) ) / clusterTools.e5x5( *(iSC->seed()) ) << std::endl;
    */
  } // photons

} // fillPhoVars()
