#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"

// Fill EB rec hits ////////////////////////////////
// Store event rechits in a vector of length equal
// to number of crystals in EB (ieta:170 x iphi:360)

TH1D *h_pT;
TH1D *h_E;
TH1D *h_eta;
TH1D *h_m0;
TH1D *h_leadJetPt;
float eventId_;
float m0_;
//float diPhoE_;
//float diPhoPt_;

// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel ( TTree* tree, edm::Service<TFileService> &fs ) {

  h_pT  = fs->make<TH1D>("h_pT" , "p_{T};p_{T};Particles", 100,  0., 500.);
  h_E   = fs->make<TH1D>("h_E"  , "E;E;Particles"        , 100,  0., 800.);
  h_eta = fs->make<TH1D>("h_eta", "#eta;#eta;Particles"  , 100, -5., 5.);
  h_m0  = fs->make<TH1D>("h_m0" , "m0;m0;Events"        ,   72, 50., 950.);
  //h_m0  = fs->make<TH1D>("h_m0" , "m0;m0;Events"        ,   50, 90., 240.);
  //h_m0  = fs->make<TH1D>("h_m0" , "m0;m0;Events"        ,   50, 0., 2.);
  h_leadJetPt  = fs->make<TH1D>("h_leadJetPt" , "p_{T};p_{T};Events", 100,  0., 500.);

  RHTree->Branch("eventId",        &eventId_);
  RHTree->Branch("m0",             &m0_);
  //RHTree->Branch("diPhoE",         &diPhoE_);
  //RHTree->Branch("diPhoPt",        &diPhoPt_);

} // branchesEB()

// Fill EB rechits _________________________________________________________________//
bool RecHitAnalyzer::runEvtSel ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<reco::PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleCollectionT_, genParticles);

  int nPho = 0;
  std::cout << "PhoCol.size: " << photons->size() << std::endl;
  math::XYZTLorentzVector vDiPho;

  for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin();
       iGen != genParticles->end();
       ++iGen) {

    // ID cuts
    if ( std::abs(iGen->pdgId()) != 22 ) continue;
    if ( iGen->status() != 1 ) continue; // NOT the same as Pythia status
    if ( !iGen->mother() ) continue;
    if ( iGen->mother()->pdgId() != 35 && iGen->mother()->pdgId() != 22 ) continue;
    //std::cout << "status:" <<iGen->status() << " pT:" << iGen->pt() << " eta:" << iGen->eta() << " E:" << iGen->energy() << " mothId:" << iGen->mother()->pdgId() << std::endl;
    nPho++;
    vDiPho += iGen->p4();

  } // genParticle loop: count good photons

  // Require exactly 2 gen-level photons
  // Indifferent about photons of status != 1
  std::cout << "nPho:" << nPho << std::endl;
  if ( nPho != 4 ) return false;
  //if ( vDiPho.mass() < 80. ) return false;

  // Fill loop
  for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin();
       iGen != genParticles->end();
       ++iGen) {

    // PDG ID cut
    if ( std::abs(iGen->pdgId()) != 22 ) continue;
    if ( iGen->status() != 1 ) continue; // NOT the same as Pythia status
    if ( !iGen->mother() ) continue;
    if ( iGen->mother()->pdgId() != 35 && iGen->mother()->pdgId() != 22 ) continue;

    // Fill histograms
    h_pT-> Fill( iGen->pt()      );
    h_E->  Fill( iGen->energy()  );
    h_eta->Fill( iGen->eta()     );
  } // genParticle loop: fill hist
  h_m0->Fill( vDiPho.mass() );
  std::cout << " m0: " << vDiPho.mass() <<" (" << vDiPho.T() << ")" << std::endl;

  m0_ = vDiPho.mass();

  // Write out event ID
  eventId_ = iEvent.id().event();

  return true;

} // fillEB()

/*
//____ Apply event selection cuts _____//
bool RecHitAnalyzer::runSelections_H24G ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<reco::PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleCollectionT_, genParticles);

  int nPho = 0;
  //bool isHiggs = true;
  //bool isDecayed = true;
  float etaCut = 1.44;
  //float etaCut = 2.3;
  //float etaCut = 2.5;
  float ptCut = 18.;
  //float dRCut = 0.4;
  //float dR, dEta, dPhi;
  std::cout << " >> recoPhoCol.size: " << photons->size() << std::endl;
  math::PtEtaPhiELorentzVectorD vDiPho;
  //math::XYZTLorentzVector vDiPho;
  std::vector<float> vE, vPt, vEta, vPhi;
  float leadPhoPt = 0;

  // Apply diphoton trigger-like selection
  for(reco::PhotonCollection::const_iterator iPho = photons->begin();
      iPho != photons->end();
      ++iPho) {

    // Kinematic cuts
    if ( std::abs(iPho->eta()) > etaCut ) continue;
    if ( std::abs(iPho->pt()) < ptCut ) continue;

    nPho++;

    // Record kinematics
    vDiPho += iPho->p4();
    vE.push_back(   iPho->energy() );
    vPt.push_back(  iPho->pt()     );
    vEta.push_back( iPho->eta()    );
    vPhi.push_back( iPho->phi()    );
    if ( std::abs(iPho->pt()) > leadPhoPt ) leadPhoPt = std::abs(iPho->pt());

  } // recoPhotons

  // Apply diphoton trigger-like selection
  if ( nPho != 2 ) return false;
  if ( leadPhoPt < 30. ) return false;
  m0_ = vDiPho.mass();
  if ( m0_ < 90. ) return false;
  //// Check dR
  //dEta = std::abs( vEta[0] - vEta[1] );
  //dPhi = std::abs( vPhi[0] - vPhi[1] );
  //dR = TMath::Power(dEta,2.) + TMath::Power(dPhi,2.);
  //dR = TMath::Sqrt(dR);
  //if ( dR < dRCut ) return false;
  std::cout << " >> passed trigger" << std::endl;

  // Apply good photon selection
  int i = 0;
  nPho = 0;
  leadPhoPt = 0.;
  for(reco::PhotonCollection::const_iterator iPho = photons->begin();
      iPho != photons->end();
      ++iPho) {

    // Kinematic cuts
    if ( std::abs(iPho->eta()) > etaCut ) continue;
    //if ( std::abs(iPho->eta()) > 1.44 && std::abs(iPho->eta()) < 1.57 ) continue;
    if ( std::abs(iPho->pt()) < m0_/4. ) continue;

    if ( std::abs(iPho->pt()) > leadPhoPt ) leadPhoPt = std::abs(iPho->pt());
    vPho_[i] = iPho->p4();
    nPho++;
    i++;

  } // recoPhotons
  if ( nPho != 2 ) return false;
  if ( leadPhoPt < m0_/3. ) return false;

  // Fill histograms
  diPhoE_  = 0.;
  diPhoPt_ = 0.;
  for(int i = 0; i < 2; i++) {
    std::cout << " >> pT:" << vPt[i] << " eta:" << vEta[i] << " phi: " << vPhi[i] << " E:" << vE[i] << std::endl;
    h_pT-> Fill( vPt[i]  );
    h_E->  Fill( vE[i]   );
    h_eta->Fill( vEta[i] );
    diPhoE_  += vE[i];
    diPhoPt_ += vPt[i];
  }
  //vDiPho = vDiPho/m0_;
  //vDiPho = vDiPho/diPhoE_;
  //vDiPho = vDiPho/diPhoPt_;
  //h_m0->Fill( vDiPho.mass() );
  //std::cout << " >> m0: " << vDiPho.mass() << " diPhoPt: " << diPhoPt_ << " diPhoE: " << diPhoE_ << std::endl;
  h_m0->Fill( m0_ );
  std::cout << " >> m0: " << m0_ << " diPhoPt: " << diPhoPt_ << " diPhoE: " << diPhoE_ << std::endl;

  //for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin();
  //     iGen != genParticles->end();
  //     ++iGen) {

  //  // ID cuts
  //  if ( std::abs(iGen->pdgId()) != 22 ) continue;
  //  if ( iGen->status() != 1 ) continue; // NOT the same as Pythia status
  //  if ( !iGen->mother() ) continue;
  //  if ( iGen->mother()->pdgId() != 35 && iGen->mother()->pdgId() != 22 ) continue;
  //  //std::cout << "status:" <<iGen->status() << " pT:" << iGen->pt() << " eta:" << iGen->eta() << " E:" << iGen->energy() << " mothId:" << iGen->mother()->pdgId() << std::endl;
  //  // Kinematic cuts
  //  if ( std::abs(iGen->eta()) > etaCut ) continue;
  //  if ( std::abs(iGen->pt()) < ptCut ) continue;
  //  nPho++;
  //  vDiPho += iGen->p4();
  //  if ( std::abs(iGen->pt()) > leadPt ) leadPt = std::abs(iGen->pt());

  //} // genParticle loop: count good photons

  //// Require exactly 2 gen-level photons
  //// Indifferent about photons of status != 1
  //std::cout << "GenCollection: " << nPho << std::endl;
  //if ( nPho != 4 ) return false;
  ////if ( vDiPho.mass() < 80. ) return false;

  //// Fill loop
  //for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin();
  //     iGen != genParticles->end();
  //     ++iGen) {

  //  // PDG ID cut
  //  if ( std::abs(iGen->pdgId()) != 22 ) continue;
  //  if ( iGen->status() != 1 ) continue; // NOT the same as Pythia status
  //  if ( !iGen->mother() ) continue;
  //  if ( iGen->mother()->pdgId() != 35 && iGen->mother()->pdgId() != 22 ) continue;
  //  // Kinematic cuts
  //  if ( std::abs(iGen->eta()) > etaCut ) continue;
  //  if ( std::abs(iGen->pt()) < ptCut ) continue;
  //  std::cout << " pT:" << iGen->pt() << " eta:" << iGen->eta() << " phi: " << iGen->phi() << " E:" << iGen->energy() << std::endl;

  //  // Fill histograms
  //  //h_pT-> Fill( iGen->pt()      );
  //  h_E->  Fill( iGen->energy()  );
  //  h_eta->Fill( iGen->eta()     );
  //} // genParticle loop: fill hist
  //h_pT-> Fill( leadPt );
  //h_m0->Fill( vDiPho.mass() );
  //std::cout << "leadPt: " << leadPt << std::endl;
  //std::cout << " m0: " << vDiPho.mass() <<" (" << vDiPho.T() << ")" << std::endl;

  //m0_ = vDiPho.mass();
  //std::cout << "PhoCol.size: " << photons->size() << std::endl;
  //for(reco::PhotonCollection::const_iterator iPho = photons->begin();
  //    iPho != photons->end();
  //    ++iPho) {
  //  if ( std::abs(iPho->eta()) > etaCut ) continue;
  //  //if ( std::abs(iPho->pt()) < ptCut-2. ) continue;
  //  std::cout << " pT:" << iPho->pt() << " eta:" << iPho->eta() << " phi: " << iPho->phi() << " E:" << iPho->energy() << std::endl;
  //}

  // Check leading jet in reco jet collection
  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);
  std::cout << " >> PFJetCol.size: " << jets->size() << std::endl;
  for(reco::PFJetCollection::const_iterator iJet = jets->begin();
      iJet != jets->end();
      ++iJet) {
    //std::cout << " pT:" << iJet->pt() << " eta:" << iJet->eta() << " phi: " << iJet->phi() << " E:" << iJet->energy() << std::endl;
  }

  // Check leading jet in gen jet collection
  float leadJetPt = 0.;
  edm::Handle<reco::GenJetCollection> genJets;
  iEvent.getByToken(genJetCollectionT_, genJets);
  std::cout << " >> GenJetCol.size: " << jets->size() << std::endl;
  for(reco::GenJetCollection::const_iterator iJet = genJets->begin();
      iJet != genJets->end();
      ++iJet) {
    if ( std::abs(iJet->pt()) > leadJetPt ) leadJetPt = std::abs(iJet->pt());
    //std::cout << " >> pT:" << iJet->pt() << " eta:" << iJet->eta() << " phi: " << iJet->phi() << " E:" << iJet->energy() << std::endl;
  }
  std::cout << " >> leadJetPt: " << leadJetPt << std::endl;
  h_leadJetPt->Fill( leadJetPt );

  return true;

} // runSelections_H24G

//____ Apply event selection cuts _____//
bool RecHitAnalyzer::runSelections_H2GG ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<reco::PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleCollectionT_, genParticles);

  int nPho = 0;
  //bool isHiggs = true;
  //bool isDecayed = true;
  //float etaCut = 1.4;
  //float etaCut = 2.3;
  //float etaCut = 5.;
  //float ptCut = 0.;
  //float dRCut = 0.4;
  //float dR, dEta, dPhi;
  std::cout << "PhoCol.size: " << photons->size() << std::endl;
  math::XYZTLorentzVector vDiPho;

  for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin();
       iGen != genParticles->end();
       ++iGen) {

    // ID cuts
    if ( std::abs(iGen->pdgId()) != 22 ) continue;
    if ( iGen->status() != 1 ) continue; // NOT the same as Pythia status
    if ( !iGen->mother() ) continue;
    if ( iGen->mother()->pdgId() != 25 && iGen->mother()->pdgId() != 22 ) continue;
    //std::cout << "status:" <<iGen->status() << " pT:" << iGen->pt() << " eta:" << iGen->eta() << " E:" << iGen->energy() << " mothId:" << iGen->mother()->pdgId() << std::endl;
    nPho++;
    vDiPho += iGen->p4();

  } // genParticle loop: count good photons

  // Require exactly 2 gen-level photons
  // Indifferent about photons of status != 1
  std::cout << nPho << std::endl;
  if ( nPho != 2 ) return false;
  //if ( vDiPho.mass() < 80. ) return false;

  // Fill loop
  for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin();
       iGen != genParticles->end();
       ++iGen) {

    // PDG ID cut
    if ( std::abs(iGen->pdgId()) != 22 ) continue;
    if ( iGen->status() != 1 ) continue; // NOT the same as Pythia status
    if ( !iGen->mother() ) continue;
    if ( iGen->mother()->pdgId() != 25 && iGen->mother()->pdgId() != 22 ) continue;

    // Fill histograms
    h_pT-> Fill( iGen->pt()      );
    h_E->  Fill( iGen->energy()  );
    h_eta->Fill( iGen->eta()     );
  } // genParticle loop: fill hist
  h_m0->Fill( vDiPho.mass() );
  std::cout << " m0: " << vDiPho.mass() <<" (" << vDiPho.T() << ")" << std::endl;

  m0_ = vDiPho.mass();

  return true;

} // runSelections_H2GG()

//____ Apply event selection cuts _____//
bool RecHitAnalyzer::runSelections ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  // Initialize data collection pointers
  edm::Handle<reco::PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleCollectionT_, genParticles);

  int nPho = 0;
  //float etaCut = 1.4;
  float etaCut = 2.3;
  float ptCut = 25.;
  std::cout << "PhoCol.size: " << photons->size() << std::endl;
  for(reco::PhotonCollection::const_iterator iPho = photons->begin();
      iPho != photons->end();
      ++iPho) {

    // Kinematic cuts
    if ( std::abs(iPho->eta()) > etaCut ) continue;
    if ( std::abs(iPho->pt()) < ptCut ) continue;

    nPho++;

  } // recoPhotons

  // Require at least 2 passed reco photons
  // Will also include PU photons
  if ( nPho < 2 ) return false;

  float dRCut = 0.4;
  float dEta, dPhi, dR;
  for(reco::PhotonCollection::const_iterator iPho = photons->begin();
      iPho != photons->end();
      ++iPho) {

    // Kinematic cuts
    if ( std::abs(iPho->eta()) > etaCut ) continue;
    if ( std::abs(iPho->pt()) < ptCut ) continue;

    std::cout << "nPho:" << nPho << " pT:" << iPho->pt() << " eta:" << iPho->eta() << " E:" << iPho->energy() << std::endl;

    for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin();
        iGen != genParticles->end();
        ++iGen) {

      // ID cuts
      if ( std::abs(iGen->pdgId()) != 22 ) continue;
      if ( iGen->status() != 1 ) continue; // NOT the same as Pythia status
      if ( !iGen->mother() ) continue;
      if ( std::abs(iGen->mother()->status()) != 44 && std::abs(iGen->mother()->status()) != 23 ) continue;

      // Match by dR
      dEta = std::abs( iPho->eta() - iGen->eta() );
      dPhi = std::abs( iPho->phi() - iGen->phi() );
      dR = TMath::Power(dEta,2.) + TMath::Power(dPhi,2.);
      dR = TMath::Sqrt(dR);
      if ( dR < dRCut ) {
        h_pT-> Fill( iGen->pt()      );
        h_E->  Fill( iGen->energy()  );
        h_eta->Fill( iGen->eta()     );
        break;
      }

    } // genParticle loop: count good photons

  } // recoPhotons

  return true;

} // runSelections()

*/
