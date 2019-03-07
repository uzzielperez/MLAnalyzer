#include "MLAnalyzer/RecHitAnalyzer/interface/SCRegressor.h"

// Initialize branches _____________________________________________________//
void SCRegressor::branchesPhoSel_gamma ( TTree* tree, edm::Service<TFileService> &fs )
{
  hSC_pT = fs->make<TH1F>("SC_pT", "Pt", 27, 15., 150.);

  RHTree->Branch("SC_mass",   &vSC_mass_);
  RHTree->Branch("SC_DR",     &vSC_DR_);
  RHTree->Branch("SC_pT",     &vSC_pT_);

  RHTree->Branch("pho_pT",    &vPho_pT_);
  RHTree->Branch("pho_E",     &vPho_E_);
  RHTree->Branch("pho_eta",   &vPho_eta_);
  RHTree->Branch("pho_phi",   &vPho_phi_);
  //RHTree->Branch("pho_r9",    &vPho_r9_);
  //RHTree->Branch("pho_sieie", &vPho_sieie_);

}

// Run event selection ___________________________________________________________________//
bool SCRegressor::runPhoSel_gamma ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleCollectionT_, genParticles);

  edm::Handle<reco::PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);

  // Get reco (pT>5GeV) photons associated to gen photons
  float dR;
  float minDR = 100.;
  int minDR_idx = -1;
  mGenPi0_RecoPho.clear();
  for ( unsigned int iG = 0; iG < genParticles->size(); iG++ ) {

    reco::GenParticleRef iGen( genParticles, iG );

    if ( iGen->pdgId() != 22 ) continue;
    if ( iGen->status() != 1 ) continue;

    mGenPi0_RecoPho.insert( std::pair<unsigned int, std::vector<unsigned int>>(iG, std::vector<unsigned int>()) );

    // Loop over reco photon collection
    minDR = 100.;
    minDR_idx = -1;
    for ( unsigned int iP = 0; iP < photons->size(); iP++ ) {

      reco::PhotonRef iPho( photons, iP );
      if ( iPho->pt() < 5. ) continue;
      dR = reco::deltaR( iPho->eta(),iPho->phi(), iGen->eta(),iGen->phi() );
      if ( dR > minDR ) continue;

      minDR = dR;
      minDR_idx = iP;
      if ( debug ) std::cout << "   >> minDR_idx:" << minDR_idx << " " << minDR << std::endl;

    } // reco photons
    if ( minDR > 0.04 ) continue;
    mGenPi0_RecoPho[iG].push_back( minDR_idx );
    if ( debug ) std::cout << "   >> !minDR_idx:" << minDR_idx << std::endl;

  } // gen photons

  ////////// Apply selection //////////

  float ptCut = 15., etaCut = 1.44;

  if ( debug ) std::cout << " >> PhoCol.size: " << photons->size() << std::endl;
  // Ensure only 1 reco photon associated to each gen pho 
  std::vector<int> vRecoPhoIdxs_;
  for ( auto const& mG : mGenPi0_RecoPho ) {
    if ( mG.second.empty() ) continue;
    if ( mG.second.size() != 1 ) continue;
    vRecoPhoIdxs_.push_back( mG.second[0] );
  }
  if ( debug ) std::cout << " >> RecoPhos.size: " << vRecoPhoIdxs_.size() << std::endl;
  if ( vRecoPhoIdxs_.empty() ) return false;

  // Ensure each pre-selected photon is isolated
  bool isIso;
  vPreselPhoIdxs_.clear();
  for ( unsigned int i = 0; i < vRecoPhoIdxs_.size(); i++ ) {

    reco::PhotonRef iPho( photons, vRecoPhoIdxs_[i] );

    if ( iPho->pt() < ptCut ) continue;
    if ( std::abs(iPho->eta()) > etaCut ) continue;
    if ( debug ) std::cout << " >> pT: " << iPho->pt() << " eta: " << iPho->eta() << std::endl;
    if ( iPho->r9() < 0.5 ) continue;
    if ( iPho->hadTowOverEm() > 0.07 ) continue;
    if ( iPho->full5x5_sigmaIetaIeta() > 0.0105 ) continue;
    if ( iPho->hasPixelSeed() == true ) continue;

    ///*
    // Ensure pre-sel photons are isolated 
    isIso = true;
    for ( unsigned int j = 0; j < vRecoPhoIdxs_.size(); j++ ) {

      if ( i == j ) continue;
      reco::PhotonRef jPho( photons, vRecoPhoIdxs_[j] ); 
      dR = reco::deltaR( iPho->eta(),iPho->phi(), jPho->eta(),jPho->phi() );
      if ( debug ) std::cout << "   >> reco dR:" << dR << std::endl;
      if ( dR > 16*.0174 ) continue;
      isIso = false;
      break;

    } // reco photon j
    if ( !isIso ) continue;
    //*/

    vPreselPhoIdxs_.push_back( vRecoPhoIdxs_[i] );

  } // reco photon i
  if ( debug ) std::cout << " >> PreselPhos.size: " << vPreselPhoIdxs_.size() << std::endl;
  if ( vPreselPhoIdxs_.empty() ) return false;

  if ( debug ) std::cout << " >> Passed selection. " << std::endl;
  return true;
}

// Fill branches ___________________________________________________________________//
void SCRegressor::fillPhoSel_gamma ( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  edm::Handle<reco::PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleCollectionT_, genParticles);

  // Eliminate gen pi0s which are unmatched
  for ( auto const& mG : mGenPi0_RecoPho ) {
    if ( mG.second.empty() ) {
      mGenPi0_RecoPho.erase( mG.first );
      continue;
    }
    if ( std::find(vRegressPhoIdxs_.begin(), vRegressPhoIdxs_.end(), mG.second[0]) == vRegressPhoIdxs_.end() )
      mGenPi0_RecoPho.erase( mG.first );
  }
  if ( debug ) std::cout << " >> mGenPi0.size: " << mGenPi0_RecoPho.size() << std::endl;

  if ( debug ) {
    int iter_idx = 0;
    for ( auto const& mG : mGenPi0_RecoPho ) {
      std::cout << " >> mGenPi0_RecoPho[" << mG.first << "]: " << mG.second[0] << " vs. " << vRegressPhoIdxs_[iter_idx] << std::endl;
      iter_idx++;
    }
  }

  ////////// Store kinematics //////////
  vPho_pT_.clear();
  vPho_E_.clear();
  vPho_eta_.clear();
  vPho_phi_.clear();
  //vPho_r9_.clear();
  //vPho_sieie_.clear();
  for ( auto const& mG : mGenPi0_RecoPho ) {
    reco::PhotonRef iPho( photons, mG.second[0] );
    // Fill branch arrays
    vPho_pT_.push_back( iPho->pt() );
    vPho_E_.push_back( iPho->energy() );
    vPho_eta_.push_back( iPho->eta() );
    vPho_phi_.push_back( iPho->phi() );
    //vPho_r9_.push_back( iPho->full5x5_r9() );
    //vPho_sieie_.push_back( iPho->full5x5_sigmaIetaIeta() );
  } // photons

  vSC_DR_.clear();
  vSC_pT_.clear();
  vSC_mass_.clear();
  for ( auto const& mG : mGenPi0_RecoPho ) {
    reco::GenParticleRef iGen( genParticles, mG.first );
    vSC_DR_.push_back( 0. );
    vSC_pT_.push_back( iGen->pt() ); 
    vSC_mass_.push_back( 0. );

    hSC_pT->Fill( iGen->pt() );
  }

}
