#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"
#include "MLAnalyzer/RecHitAnalyzer/interface/HighMassPhotonID.h"

// Run event selection ////////////////////////////////
//const std::string jetSelection = "dijet_gg_qq"; // TODO: put switch at cfg level

TH1F *h_m0;
TH1F *h_nJet;
TH1F *h_phoPt;
TH1F *h_phoE;
TH1F *h_phoEta;
TH1F *h_phoR9;
TH1F *h_phoSieie;
TH1F *h_phoMva;
TH1F *h_jetPt;
TH1F *h_jetE;
TH1F *h_jetEta;

unsigned long long eventId_;
unsigned int runId_;
unsigned int lumiId_;
float m0_;
//float nJet_;
float diPhoE_;
float diPhoPt_;
std::vector<float> vFC_inputs_;


float m0cut = 500.;
//float m0cut = 80.;
float minPt = 75.;
float hOverEcut = 0.05;
float chIsocut = 5.;
float sIeIeCut = -1.;

// passedElectronVeto
// passedHadTowerOverEmCut(photon)
// passedChargedHadronCut(photon)

// passedSigmaIetaIetaCut(photon,isSat)
// passedCorPhoIsoHighPtID(photon,rho)


// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel_highMass ( TTree* tree, edm::Service<TFileService> &fs ) {

    // Default mbandrews selection

    h_m0     = fs->make<TH1F>("h_m0"    , "m0;m0;Events"         ,  100, m0cut, 13000.);
    h_phoPt  = fs->make<TH1F>("h_phoPt" , "p_{T};p_{T};Particles", 100,  0., 4000.);
    h_phoE   = fs->make<TH1F>("h_phoE"  , "E;E;Particles"        , 100,  0., 13000.);
    h_phoEta = fs->make<TH1F>("h_phoEta", "#eta;#eta;Particles"  , 100, -5., 5.);
    h_phoR9  = fs->make<TH1F>("h_phoR9" , "R9;R9;Particles"  , 50, 0., 1.);
    h_phoSieie  = fs->make<TH1F>("h_phoSieie" , "Sieie;Sieie;Particles"  , 50, 0., 0.1);

  // ExoDiPhoton
  tree->Branch("eventId",        &eventId_);
  tree->Branch("runId",          &runId_);
  tree->Branch("lumiId",         &lumiId_);
  tree->Branch("m0",             &m0_);
  tree->Branch("FC_inputs",      &vFC_inputs_);
  tree->Branch("diPhoE",         &diPhoE_);
  tree->Branch("diPhoPt",        &diPhoPt_);

} // branchesEvtSel()

// Run event selection _______________________________________________________________//
bool RecHitAnalyzer::runEvtSel_highMass ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<reco::PhotonCollection> photons;
  //edm::Handle<pat::PhotonCollection> photons;
  iEvent.getByToken( photonCollectionT_, photons );

  int nPhoTrg = 0;
  float dR, m0;
  float leadPhoPt = 0.;
  math::PtEtaPhiELorentzVectorD vDiPho;
  std::vector<int> vPhoIdxs;
  // edm::Handle<edm::View<pat::Photon> > photons;
  // std::vector<edm::Ptr<pat::Photon>> goodPhotons;

  // Photon pre-selection
  // High Mass Photon ID
  // Store Photons in vector

  for ( unsigned int iP = 0; iP < photons->size(); iP++ ) {
    reco::PhotonRef iPho( photons, iP );
    //pat::PhotonRef iPho( photons, iP );
    if ( std::abs(iPho->pt()) < minPt ) continue;
    //std::cout << iPho->full5x5_sigmaIetaIeta() << std::endl;
    if ( std::abs(iPho->eta()) > 1.44 ) continue; // barrel photons only
    if ( iPho->r9() < 0.5 ) continue;
    if ( iPho->hadTowOverEm() > hOverEcut ) continue;
    if ( iPho->full5x5_sigmaIetaIeta() > 0.0105 ) continue; // re-implement with passSigmaIetaIetaCut
    if ( iPho->hasPixelSeed() == true ) continue;
    //if ( std::abs(iPho->eta()) > 2.1 ) continue;
    //if ( std::abs(iPho->eta()) > 1.44 && std::abs(iPho->eta()) < 1.57 ) continue;
    if ( iPho->hasPixelSeed() == true ) continue;
    // if ( iPho->passElectronVeto() ) continue; // isElectron
    // Implement passCorPhoIso (photon, rho) - understand rho variable
    if (debug) std::cout << " >> pT:" << iPho->pt() << " eta:" << iPho->eta() << " phi: " << iPho->phi() << " E:" << iPho->energy() << std::endl;
    nPhoTrg++;
    if ( std::abs(iPho->pt()) > leadPhoPt ) leadPhoPt = std::abs(iPho->pt());
    vDiPho += iPho->p4();
    vPhoIdxs.push_back( iP );
  } // reco photons

  m0 = vDiPho.mass();
  if ( m0 < m0cut ) return false;
  if ( nPhoTrg != 2 ) return false;
  //if ( leadPhoPt < 30. ) return false;
  if ( leadPhoPt < minPt ) return false;

  // Apply selection
  int nPho = 0;
  int leadPho = -1;
  leadPhoPt = 0.;
  for ( int iP = 0; iP < nPhoTrg; iP++ ) {
    reco::PhotonRef iPho( photons, vPhoIdxs[iP] );
    //pat::PhotonRef iPho( photons, vPhoIdxs[iP] );
    // Get leading photon pt
    if ( std::abs(iPho->pt()) > leadPhoPt ) {
      leadPhoPt = std::abs(iPho->pt());
      leadPho = iP;
    }
    // Minimum pt/m0 cut
    if ( std::abs(iPho->pt()) < m0/4. ) continue;
    nPho++;
  }
  if ( nPho != 2 ) return false;
  //if ( leadPhoPt < m0/3 ) return false;
  nPho = nPhoTrg;
  //std::cout << " n:" << nPho << " m0:" << m0 << std::endl;

  // Get photon order
  int ptOrder[2] = {0, 1};
  if ( leadPho == 1 ) {
    ptOrder[0] = 1;
    ptOrder[1] = 0;
  }

  // Fill kinematic variables
  h_m0->Fill( m0 );
  diPhoE_  = 0.;
  diPhoPt_ = 0.;
  float dphi[2] = {0., 0.};
  vFC_inputs_.clear();
  for ( int iP = 0; iP < nPho; iP++ ) {
    reco::PhotonRef iPho( photons, vPhoIdxs[ptOrder[iP]] );
    //pat::PhotonRef iPho( photons, vPhoIdxs[ptOrder[iP]] );
    h_phoPt->Fill( iPho->pt() );
    h_phoE->Fill( iPho->energy() );
    h_phoEta->Fill( iPho->eta() );
    h_phoR9->Fill( iPho->r9() );
    h_phoSieie->Fill( iPho->full5x5_sigmaIetaIeta() );

    diPhoE_  += std::abs( iPho->energy() );
    diPhoPt_ += std::abs( iPho->pt() );
    vFC_inputs_.push_back( iPho->pt()/m0 );
    vFC_inputs_.push_back( iPho->eta() );
    dphi[iP] = iPho->phi();
  }
  vFC_inputs_.push_back( TMath::Cos(reco::deltaPhi(dphi[0], dphi[1])) );

  // Write out event
  m0_ = m0;
  //nJet_ = nJet;
  eventId_ = iEvent.id().event();
  runId_ = iEvent.id().run();
  lumiId_ = iEvent.id().luminosityBlock();

  return true;

} // runEvtSel()
