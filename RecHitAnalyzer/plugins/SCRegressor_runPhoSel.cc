#include "MLAnalyzer/RecHitAnalyzer/interface/SCRegressor.h"

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

// Initialize branches _____________________________________________________//
void SCRegressor::branchesPhoSel ( TTree* tree, edm::Service<TFileService> &fs )
{
  hSC_mass = fs->make<TH1F>("SC_mass", "m_{SC};m_{SC}",50, 0., 0.5);
  hDR = fs->make<TH1F>("DR_seed_subJet", "#DeltaR(seed,subJet);#DeltaR",50, 0., 50.*0.0174);
  hdEta = fs->make<TH1F>("dEta_seed_subJet", "#Delta#eta(seed,subJet);#Delta#eta",50, 0., 50.*0.0174);
  hdPhi = fs->make<TH1F>("dPhi_seed_subJet", "#Delta#phi(seed,subJet);#Delta#phi",50, 0., 50.*0.0174);
  hdPhidEta = fs->make<TH3F>("dPhidEta_GG", "#Delta(#phi,#eta,m);#Delta#phi(#gamma,#gamma);#Delta#eta(#gamma,#gamma);m",6, 0., 6.*0.0174, 6, 0., 6.*0.0174, 16., 0.,1.6);
  hPt = fs->make<TH1F>("Pt", "Pt", 65, 30., 160.);

  char hname[50];
  for(int iPho (0); iPho < nPhotons; iPho++) {
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
    sprintf(hname, "pho_sieie%d",iPho);
    RHTree->Branch(hname,      &vPho_sieie_[iPho]);
  }

  hnPho = fs->make<TH2F>("nPho", "N(m_{#pi},p_{T,#pi})_{reco};m_{#pi^{0}};p_{T,#pi^0}",
      16, 0., 1.6, 17, 15., 100.);
}

// Run event selection ___________________________________________________________________//
bool SCRegressor::runPhoSel ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<reco::PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleCollectionT_, genParticles);

  //std::vector<int> vGenPi0Idxs_;
  vGenPi0Idxs_.clear();

  // ____________ Gen-level studies ________________ //
  int nPi0 = 0;
  float dR;
  //std::vector<math::PtEtaPhiELorentzVectorD> vPhoPairs[nPhotons];
  //std::vector<math::XYZTLorentzVector> vPhoPairs[nPhotons];
  std::map<unsigned int, std::vector<unsigned int>> mGenPi0_RecoPho;
  for ( unsigned int iG = 0; iG < genParticles->size(); iG++ ) {

    reco::GenParticleRef iGen( genParticles, iG );
    // ID cuts
    if ( debug ) std::cout << " >> pdgId:"<< iGen->pdgId() << " status:" << iGen->status() << " nDaughters:" << iGen->numberOfDaughters() << std::endl;
    if ( std::abs(iGen->pdgId()) != 111 ) continue;
    if ( debug ) std::cout << " >> pdgId:111 nDaughters:" << iGen->numberOfDaughters() << std::endl;
    if ( iGen->numberOfDaughters() != 2 ) continue;
    //if ( iGen->mass() > 0.4 ) continue;
    dR = reco::deltaR( iGen->daughter(0)->eta(),iGen->daughter(0)->phi(), iGen->daughter(1)->eta(),iGen->daughter(1)->phi() );
    if ( dR > 5*.0174 ) continue;

    vGenPi0Idxs_.push_back( iG );
    mGenPi0_RecoPho.insert( std::pair<unsigned int, std::vector<unsigned int>>(iG, std::vector<unsigned int>()) );
    nPi0++;

  } // genParticle loop: count good photons
  if ( debug ) std::cout << " >> mGenPi0.size: " << mGenPi0_RecoPho.size() << std::endl;

  if ( nPi0 != nPhotons ) return false;

  ////////// Apply selection //////////

  //float ptCut = 45., etaCut = 1.44;
  float ptCut = 15., etaCut = 1.44;

  // Get reco (pT>5GeV) photons associated to pi0
  float minDR = 100.;
  int minDR_idx = -1;
  // Loop over pi0s
  for ( auto& mG : mGenPi0_RecoPho ) {

    reco::GenParticleRef iGen( genParticles, mG.first );

    if ( debug ) std::cout << " >> pi0[" << mG.first << "]"<< std::endl;
    // Loop over photons from pi0
    for ( unsigned int iD = 0; iD < iGen->numberOfDaughters(); iD++ ) {

      const reco::Candidate* iGenPho = iGen->daughter(iD);

      // Loop over reco photon collection
      minDR = 100.;
      minDR_idx = -1;
      if ( debug ) std::cout << "  >> iD[" << iD << "]" << std::endl;
      for ( unsigned int iP = 0; iP < photons->size(); iP++ ) {

        reco::PhotonRef iPho( photons, iP );
        if ( iPho->pt() < 5. ) continue;
        dR = reco::deltaR( iPho->eta(),iPho->phi(), iGenPho->eta(),iGenPho->phi() );
        if ( dR > minDR ) continue;

        minDR = dR;
        minDR_idx = iP;
        if ( debug ) std::cout << "   >> minDR_idx:" << minDR_idx << " " << minDR << std::endl;

      } // reco photons
      if ( minDR > 0.04 ) continue;
      //mG[mG.first].push_back( minDR_idx );
      mG.second.push_back( minDR_idx );
      if ( debug ) std::cout << "   >> !minDR_idx:" << minDR_idx << std::endl;

    } // gen photons 
    dR = reco::deltaR( iGen->daughter(0)->eta(),iGen->daughter(0)->phi(), iGen->daughter(1)->eta(),iGen->daughter(1)->phi() );
    if ( debug ) std::cout << "   >> gen dR:" << dR << std::endl;

  } // gen pi0s

  // Ensure only 1 reco photon associated to each gen pi0
  std::vector<int> vRecoPhoIdxs_;
  for ( auto const& mG : mGenPi0_RecoPho ) {
    if ( debug ) std::cout << " >> pi0[" << mG.first << "] size:" << mG.second.size() << std::endl; 
    // Possibilities:
    // 1) pho_gen1: unmatched, pho_gen2: unmatched => both reco pho failed pT cut or dR matching, reject
    // 2) pho_gen1: reco-matched1, pho_gen1: unmatched => one reco pho failed pT cut or dR matching, can accept
    // 3) pho_gen1: reco-matched1, pho_gen2: reco-matched2
    //    3a) reco-matched1 == reco-matched2 => merged, accept
    //    3b) reco-matched1 != reco-matched2 => resolved, reject
    if ( mG.second.size() == 0 ) {
      //mGenPi0_RecoPho.erase( mG.first );
      continue;
    }
    if ( mG.second.size() == 2 && mG.second[0] != mG.second[1] ) {
      //mGenPi0_RecoPho.erase( mG.first );
      continue; // 2 resolved reco photons 
    }
    vRecoPhoIdxs_.push_back( mG.second[0] );
  } 
  if ( debug ) std::cout << " >> RecoPhos.size: " << vRecoPhoIdxs_.size() << std::endl;
  if ( debug ) std::cout << " >> mGenPi0.size: " << mGenPi0_RecoPho.size() << std::endl;

  // Ensure each of pi0-matched reco photons passes pre-selection
  // NOTE: at this point there is 1-1 correspondence between pi0 and reco pho,
  // so suffices to check there are as many pre-selected phos as pi0s
  bool isIso;
  vPreselPhoIdxs_.clear();
  if ( debug ) std::cout << " >> PhoCol.size: " << photons->size() << std::endl;
  for ( unsigned int i = 0; i < vRecoPhoIdxs_.size(); i++ ) {

    reco::PhotonRef iPho( photons, vRecoPhoIdxs_[i] );

    if ( iPho->pt() < ptCut ) continue;
    if ( std::abs(iPho->eta()) > etaCut ) continue;
    if ( debug ) std::cout << " >> pT: " << iPho->pt() << " eta: " << iPho->eta() << std::endl;
    if ( iPho->r9() < 0.5 ) continue;
    if ( iPho->hadTowOverEm() > 0.07 ) continue;
    if ( iPho->full5x5_sigmaIetaIeta() > 0.0105 ) continue;
    if ( iPho->hasPixelSeed() == true ) continue;

    // Ensure pre-sel photons are isolated 
    isIso = true;
    for ( unsigned int j = 0; j < vRecoPhoIdxs_.size(); j++ ) {

      if ( i == j ) continue;
      reco::PhotonRef jPho( photons, vRecoPhoIdxs_[j] ); 
      dR = reco::deltaR( iPho->eta(),iPho->phi(), jPho->eta(),jPho->phi() );
      if ( debug ) std::cout << "   >> reco dR:" << dR << std::endl;
      if ( dR > 12*.0174 ) continue;
      isIso = false;
      break;

    } // reco photon j
    if ( !isIso ) continue;

    vPreselPhoIdxs_.push_back( vRecoPhoIdxs_[i] );

  } // reco photon i
  if ( debug ) std::cout << " >> PreSelPhos.size: " << vRecoPhoIdxs_.size() << std::endl;

  if ( vPreselPhoIdxs_.size() != nPhotons ) return false;
  //TODO: allow variable number of pre-selected photons to be saved per event
  //TODO: get rid of gen pi0s which are not regressed

  if ( debug ) std::cout << " >> Passed selection. " << std::endl;
  return true;
}

// Fill branches ___________________________________________________________________//
void SCRegressor::fillPhoSel ( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  edm::Handle<reco::PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);

  ////////// Store each shower crop //////////
  for ( unsigned int i = 0; i < vRegressPhoIdxs_.size(); i++ ) {
    reco::PhotonRef iPho( photons, vRegressPhoIdxs_[i] );
    // Fill branch arrays
    vPho_pT_[i] = iPho->pt();
    vPho_E_[i] = iPho->energy();
    vPho_eta_[i] = iPho->eta();
    vPho_phi_[i] = iPho->phi();
    //std::cout << "r9: " << iPho->r9() << std::endl;
    vPho_r9_[i] = iPho->r9();
    vPho_sieie_[i] = iPho->full5x5_sigmaIetaIeta(); 
  } // photons

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleCollectionT_, genParticles);

  float dEta, dPhi, dR, mPi0, ptPi0;
  //for ( int i = 0; i < nPi0; i++ ) {
  for ( unsigned int i = 0; i < vGenPi0Idxs_.size(); i++ ) {
    reco::GenParticleRef iGen( genParticles, vGenPi0Idxs_[i] );
    //mPi0 = (vPhoPairs[i][0] + vPhoPairs[i][1]).mass();
    //dR = reco::deltaR( vPhoPairs[i][0].eta(),vPhoPairs[i][0].phi(), vPhoPairs[i][1].eta(),vPhoPairs[i][1].phi() );
    //dEta = std::abs( vPhoPairs[i][0].eta() - vPhoPairs[i][1].eta() );
    //dPhi = reco::deltaPhi( vPhoPairs[i][0].phi(), vPhoPairs[i][1].phi() );
    mPi0 = iGen->mass();
    ptPi0 = iGen->pt();

    dR = reco::deltaR( iGen->daughter(0)->eta(),iGen->daughter(0)->phi(), iGen->daughter(1)->eta(),iGen->daughter(1)->phi() );
    dEta = std::abs( iGen->daughter(0)->eta() - iGen->daughter(1)->eta() );
    dPhi = reco::deltaPhi( iGen->daughter(0)->phi(), iGen->daughter(1)->phi() );

    if ( debug ) std::cout << " >> m0:" << mPi0 << " dR:" << dR << " dPhi:" << dPhi << std::endl;
    vSC_DR_[i] = dR;
    vSC_pT_[i] = iGen->pt(); 
    vSC_mass_[i] = mPi0;

    hPt->Fill( ptPi0 );
    hdPhidEta->Fill( dPhi, dEta, mPi0 );
    //hnPho->Fill( mPi0, iGen->pt() );
    hnPho->Fill( mPi0, ptPi0 );
    hSC_mass->Fill( mPi0 );
  }

}
