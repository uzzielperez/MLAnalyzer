#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"
#include "Calibration/IsolatedParticles/interface/DetIdFromEtaPhi.h"
/*
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
*/

// Run jet event selection ////////////////////////////////

//const bool debug = true;
const bool debug = false;
const unsigned nJets = 2;
const int search_window = 7;
//const int image_padding = 14;
const int image_padding = 12;
TH1D *h_jet_pT;
TH1D *h_jet_E;
TH1D *h_jet_eta;
TH1D *h_jet_m0;
TH1D *h_jet_nJet;
TH1D *h_nGG;
TH1D *h_nQQ;
unsigned int jet_eventId_;
unsigned int jet_runId_;
unsigned int jet_lumiId_;
float vJet_m0_[nJets];
float vJet_pt_[nJets];
float vJetSeed_iphi_[nJets];
float vJetSeed_ieta_[nJets];
float vJetIsQuark_[nJets];
float vJetIds_[nJets];
std::vector<float> vSubJetE_[nJets];
std::vector<float> vSubJetPx_[nJets];
std::vector<float> vSubJetPy_[nJets];
std::vector<float> vSubJetPz_[nJets];
std::vector<int> vGoodJetIdxs;
std::vector<int> vJetIds;

// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel_jet ( TTree* tree, edm::Service<TFileService> &fs ) {

  h_jet_pT    = fs->make<TH1D>("h_jet_pT"  , "p_{T};p_{T};Particles", 100,  0., 500.);
  h_jet_E     = fs->make<TH1D>("h_jet_E"   , "E;E;Particles"        , 100,  0., 800.);
  h_jet_eta   = fs->make<TH1D>("h_jet_eta" , "#eta;#eta;Particles"  , 100, -5., 5.);
  h_jet_nJet  = fs->make<TH1D>("h_jet_nJet", "nJet;nJet;Events"     ,  10,  0., 10.);
  h_jet_m0    = fs->make<TH1D>("h_jet_m0"  , "m0;m0;Events"         ,  60, 50., 110.);
  h_nGG       = fs->make<TH1D>("h_nGG"     , "nGG;nGG;Events"       ,   3,  0.,   3.);
  h_nQQ       = fs->make<TH1D>("h_nQQ"     , "nQQ;nQQ;Events"       ,   3,  0.,   3.);

  char hname[50];
  RHTree->Branch("eventId",        &jet_eventId_);
  RHTree->Branch("runId",          &jet_runId_);
  RHTree->Branch("lumiId",         &jet_lumiId_);
  for ( unsigned iJ = 0; iJ != nJets; iJ++ ) {
    sprintf(hname, "jetM%d",  iJ);
    RHTree->Branch(hname,            &vJet_m0_[iJ]);
    sprintf(hname, "jetPt%d", iJ);
    RHTree->Branch(hname,            &vJet_pt_[iJ]);
    sprintf(hname, "jetSeed_iphi%d", iJ);
    RHTree->Branch(hname,            &vJetSeed_iphi_[iJ]);
    sprintf(hname, "jetSeed_ieta%d", iJ);
    RHTree->Branch(hname,            &vJetSeed_ieta_[iJ]);
    sprintf(hname, "jetIsQuark%d", iJ);
    RHTree->Branch(hname,            &vJetIsQuark_[iJ]);
    sprintf(hname, "jetIds%d", iJ);
    RHTree->Branch(hname,            &vJetIds_[iJ]);


    sprintf(hname, "subJet%d_E", iJ);
    RHTree->Branch(hname,            &vSubJetE_[iJ]);
    sprintf(hname, "subJet%d_Px", iJ);
    RHTree->Branch(hname,            &vSubJetPx_[iJ]);
    sprintf(hname, "subJet%d_Py", iJ);
    RHTree->Branch(hname,            &vSubJetPy_[iJ]);
    sprintf(hname, "subJet%d_Pz", iJ);
    RHTree->Branch(hname,            &vSubJetPz_[iJ]);
  }

} // branchesEvtSel_jet()

// Run event selection ___________________________________________________________________//
bool RecHitAnalyzer::runEvtSel_jet ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<HBHERecHitCollection> HBHERecHitsH_;
  iEvent.getByLabel( HBHERecHitCollectionT_, HBHERecHitsH_ );

  edm::ESHandle<CaloGeometry> caloGeomH_;
  iSetup.get<CaloGeometryRecord>().get( caloGeomH_ );
  const CaloGeometry* caloGeom = caloGeomH_.product();

  /*
  edm::ESHandle<CaloTopology> caloTopoH_;
  iSetup.get<CaloTopologyRecord>().get( caloTopoH_ );
  const CaloTopology *caloTopo = caloTopoH_.product();
  */
  bool genPassed = has_dijet( iEvent, iSetup );
  //bool genPassed = has_w2jet_z2invisible( iEvent, iSetup );
  if ( !genPassed ) return false; 
  //if ( genPassed ) return true; 

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByLabel(jetCollectionT_, jets);
  if ( debug ) std::cout << " >> PFJetCol.size: " << jets->size() << std::endl;

  HcalSubdetector subdet_;
  HcalDetId seedId;
  float seedE;
  //int nJet = 0;
  //int jetIdx = -1;
  int iphi_, ieta_, ietaAbs_;

  // Loop over jets
  for ( unsigned iJ(0); iJ != nJets; ++iJ ) {

    reco::PFJetRef iJet( jets, vGoodJetIdxs[iJ] );

    // Get closest HBHE tower to jet position
    // This will not always be the most energetic deposit
    HcalDetId hId( spr::findDetIdHCAL( caloGeom, iJet->eta(), iJet->phi(), false ) );
    if ( hId.subdet() != HcalBarrel && hId.subdet() != HcalEndcap ) continue;
    HBHERecHitCollection::const_iterator iRHit( HBHERecHitsH_->find(hId) );
    seedE = ( iRHit == HBHERecHitsH_->end() ) ? 0. : iRHit->energy();
    seedId = hId;
    if ( debug ) std::cout << " >> hId.ieta:" << hId.ieta() << " hId.iphi:" << hId.iphi() << " E:" << seedE << std::endl;

    // Look for the highest HBHE tower deposit within a search window
    for ( int ieta = 0; ieta < search_window; ieta++ ) {
      ieta_   = hId.ieta() - (search_window/2)+ieta;
      if ( std::abs(ieta_) > HBHE_IETA_MAX_HE-1 ) continue;
      if ( std::abs(ieta_) < HBHE_IETA_MIN_HB ) continue;
      subdet_ = std::abs(ieta_) > HBHE_IETA_MAX_HB ? HcalEndcap : HcalBarrel;
      for ( int iphi = 0; iphi < search_window; iphi++ ) {
        iphi_   = hId.iphi() - (search_window/2)+iphi;
        if ( iphi_ > HBHE_IPHI_MAX ) {
          iphi_ = iphi_-HBHE_IPHI_MAX;
        } else if ( iphi_ < HBHE_IPHI_MIN ) {
          iphi_ = HBHE_IPHI_MAX-abs(iphi_); 
        }

        //if ( debug ) std::cout << " !! ieta_:" << ieta_ << " iphi_:" << iphi_ << std::endl;
        HcalDetId hId_( subdet_, ieta_, iphi_, 1 );
        HBHERecHitCollection::const_iterator iRHit( HBHERecHitsH_->find(hId_) );
        if ( iRHit == HBHERecHitsH_->end() ) continue;
        if ( iRHit->energy() <= seedE ) continue;
        //if ( iRHit->energy() <= seedE ) ;
        if ( debug ) std::cout << " !! hId.ieta:" << hId_.ieta() << " hId.iphi:" << hId_.iphi() << " E:" << iRHit->energy() << std::endl;
        seedE = iRHit->energy();
        seedId = hId_;

      } // iphi 
    } // ieta

    // NOTE: HBHE iphi = 1 does not correspond to EB iphi = 1!
    // => Need to shift by 2 HBHE towers: HBHE::iphi: [1,...,71,72]->[3,4,...,71,72,1,2]
    //float eta, phi;
    //GlobalPoint pos;
    //EBDetId ebId( 1, 1 );
    //pos = caloGeom->getPosition( ebId );
    //eta = pos.eta();
    //phi = pos.phi();
    //HcalDetId hcalebId( spr::findDetIdHCAL( caloGeom, eta, phi, false ) );
    //seedId = hcalebId; 
    iphi_  = seedId.iphi() + 2; // shift
    iphi_  = iphi_ > HBHE_IPHI_MAX ? iphi_-HBHE_IPHI_MAX : iphi_; // wrap-around
    iphi_  = iphi_ - 1; // make histogram-friendly
    ietaAbs_  = seedId.ietaAbs() == HBHE_IETA_MAX_HE ? HBHE_IETA_MAX_HE-1 : seedId.ietaAbs(); // last HBHE ieta embedded
    ieta_  = seedId.zside() > 0 ? ietaAbs_-1 : -ietaAbs_;
    ieta_  = ieta_+HBHE_IETA_MAX_HE-1;

    // If the seed is too close to the edge of HE, discard event
    // Required to keep the seed at the image center
    if ( HBHE_IETA_MAX_HE-1 - ietaAbs_ < image_padding ) return false;

    // Save position of highest HBHE tower
    // in EB-aligned coordinates
    if ( debug ) std::cout << " !! ieta_:" << ieta_ << " iphi_:" << iphi_ << " ietaAbs_:" << ietaAbs_ << " E:" << seedE << std::endl;
    vJetSeed_iphi_[iJ] = iphi_;
    vJetSeed_ieta_[iJ] = ieta_;

    // Plot kinematics
    h_jet_pT->Fill( std::abs(iJet->pt()) );
    h_jet_eta->Fill( iJet->eta() );
    h_jet_E->Fill( iJet->energy() );
    h_jet_m0->Fill( iJet->mass() );

    // Write out event ID
    vJet_pt_[iJ] = iJet->pt();
    vJet_m0_[iJ] = iJet->mass();
    vJetIsQuark_[iJ] = vJetIds[iJ] < 4 ? true : false;
    vJetIds_[iJ] = vJetIds[iJ];

    vSubJetE_[iJ].clear();
    vSubJetPx_[iJ].clear();
    vSubJetPy_[iJ].clear();
    vSubJetPz_[iJ].clear();
    //std::vector<reco::PFCandidatePtr> jetConstituents = iJet->getPFConstituents();
    unsigned int nConstituents = iJet->getPFConstituents().size();
    for ( unsigned int j = 0; j < nConstituents; j++ ) {
      const reco::PFCandidatePtr subJet = iJet->getPFConstituent( j );
      if ( debug ) std::cout << " >> " << j << ": E:" << subJet->energy() << " px:" << subJet->px() << " py:" << subJet->py() << " pz:" << subJet->pz() << std::endl;
      vSubJetE_[iJ].push_back( subJet->energy() );
      vSubJetPx_[iJ].push_back( subJet->px() );
      vSubJetPy_[iJ].push_back( subJet->py() );
      vSubJetPz_[iJ].push_back( subJet->pz() );
    }

  } // good jets 

  jet_eventId_ = iEvent.id().event();
  jet_runId_ = iEvent.id().run();
  jet_lumiId_ = iEvent.id().luminosityBlock();
  if ( debug ) std::cout << " >> analyze: passed" << std::endl;
  return true;

} // runEvtSel_jet()

bool RecHitAnalyzer::has_dijet( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel( genParticleCollectionT_, genParticles );

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByLabel(jetCollectionT_, jets);
  float dR;

  vGoodJetIdxs.clear();
  vJetIds.clear();
  /*
  edm::Handle<reco::GenJetCollection> genJets;
  iEvent.getByLabel(genJetCollectionT_, genJets);
  std::vector<float> vJetFakePhoIdxs;
  */

  unsigned nGG = 0;
  unsigned nQQ = 0;
  //int nQcdIn = 0;
  //int nQcdOut = 0;

  if ( debug ) std::cout << " >>>>>>>>>>>>>>>>>>>> evt:" << std::endl;
  for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin();
      iGen != genParticles->end();
      ++iGen) {

    if ( iGen->numberOfMothers() != 2 ) continue;
    if ( iGen->status() != 3 ) continue;
    if ( debug ) std::cout << " >> id:" << iGen->pdgId() << " status:" << iGen->status() << " nDaught:" << iGen->numberOfDaughters() << " pt:"<< iGen->pt() << " eta:" <<iGen->eta() << " phi:" <<iGen->phi() << " nMoms:" <<iGen->numberOfMothers()<< std::endl;

    // Loop over jets
    for ( unsigned iJ(0); iJ != jets->size(); ++iJ ) {
      //if ( debug ) std::cout << " >>>>>> jet[" << iJ << "]" << std::endl;
      reco::PFJetRef iJet( jets, iJ );
      if ( std::abs(iJet->pt())  < 70. ) continue;
      if ( std::abs(iJet->eta()) > 2. ) continue;
      dR = reco::deltaR( iJet->eta(),iJet->phi(), iGen->eta(),iGen->phi() );
      if ( debug ) {
        std::cout << " >>>>>> jet[" << iJ << "] Pt:" << iJet->pt() << " jetEta:" << iJet->eta() << " jetPhi:" << iJet->phi() 
        << " dR:" << dR << std::endl;
      }
      if ( dR > 0.4 ) continue; 
      vGoodJetIdxs.push_back( iJ );
      vJetIds.push_back( std::abs(iGen->pdgId()) );
      if ( debug ) {
        std::cout << " >>>>>> DR matched: jet[" << iJ << "] pdgId:" << std::abs(iGen->pdgId()) << std::endl;
      }
      break;
    } // reco jets

  } // gen particles 
  h_jet_nJet->Fill( vGoodJetIdxs.size() );
  if ( vGoodJetIdxs.size() != nJets ) return false;

  for ( unsigned iJ(0); iJ != nJets; ++iJ ) {
    if ( vJetIds[iJ] < 4 ) nQQ++;
    else if ( vJetIds[iJ] == 21 ) nGG++;
  }
  if ( vGoodJetIdxs[0] == vGoodJetIdxs[1] ) return false; // protect against double counting: only valid for nJets==2
  if ( nQQ+nGG != nJets ) return false; // require dijet
  if ( nQQ != nJets && nGG != nJets ) return false; // require gg or qq final state
  //if ( nQQ != 1 && nGG != 1 ) return false; // require qg final state

  h_nGG->Fill( nGG );
  h_nQQ->Fill( nQQ );
    /*
    for ( unsigned int k = 0; k < genJets->size(); k++ ) {
      reco::GenJetRef iJet( genJets, k );
      if ( std::abs(iJet->pt()) < 50. ) continue;
      dR = reco::deltaR( iJet->eta(),iJet->phi(), iGen->eta(),iGen->phi() );
      if ( debug ) std::cout << " >>>>>> genJet[" << k << "]Pt:" << iJet->pt() << " jetEta:" << iJet->eta() << " jetPhi:" << iJet->phi() << " dR:" << dR << std::endl; 
    } // genJets
    */
    /*
    // Get incoming QCD partons
    for ( unsigned int i = 0; i < iGen->numberOfDaughters(); i++ ) {
      const reco::Candidate* ppParton = iGen->daughter(i);
      if ( ppParton->status() != 3 ) continue;
      if ( debug ) std::cout << " >>>> id:" << ppParton->pdgId() << " status:" << ppParton->status() << " nDaught:" << ppParton->numberOfDaughters() << " pt:" << ppParton->pt() << " eta:" <<ppParton->eta() << " phi:" <<ppParton->phi() <<std::endl;
      // Get outgoing QCD partons
      for ( unsigned int j = 0; j < ppParton->numberOfDaughters(); j++ ) {
        const reco::Candidate* qcdIn = ppParton->daughter(j);
        if ( qcdIn->status() != 3 ) continue;
        if ( debug ) std::cout << " >>>>>> id:" << qcdIn->pdgId() << " status:" << qcdIn->status() << " nDaught:" << qcdIn->numberOfDaughters() << " pt:" << qcdIn->pt() << " eta:" <<qcdIn->eta() << " phi:" <<qcdIn->phi() << std::endl;
        nQcdIn++;
        for ( unsigned int k = 0; k < qcdIn->numberOfDaughters(); k++ ) {
          const reco::Candidate* qcdOut = qcdIn->daughter(k);
          if ( qcdOut->status() != 3 ) continue;
          if ( debug ) std::cout << " >>>>>>>> id:" << qcdOut->pdgId() << " status:" << qcdOut->status() << " nDaught:" << qcdOut->numberOfDaughters() << " pt:" << qcdOut->pt() << " eta:" <<qcdOut->eta() << " phi:" <<qcdOut->phi() << std::endl;
          //nQcdOut++;
          if ( std::abs(qcdOut->pdgId()) != 21 ) continue; 
          nGG++;

        // Loop over jets
        for ( unsigned iJ(0); iJ != jets->size(); ++iJ ) {

          //if ( debug ) std::cout << " >>>>>> jet[" << iJ << "]" << std::endl;
          reco::PFJetRef iJet( jets, iJ );
          if ( std::abs(iJet->pt()) < 50. ) continue;
          //if ( std::abs(iJet->eta()) > 2. ) continue;
          dR = reco::deltaR( iJet->eta(),iJet->phi(), qcdOut->eta(),qcdOut->phi() );
          if ( debug ) {
            std::cout << " >>>>>> jet[" << iJ << "]Pt:" << iJet->pt() << " jetEta:" << iJet->eta() << " jetPhi:" << iJet->phi() 
            << " dR:" << dR << std::endl;
          }
        } // PFjets
        for ( unsigned int k = 0; k < genJets->size(); k++ ) {
          reco::GenJetRef iJet( genJets, k );
          if ( std::abs(iJet->pt()) < 50. ) continue;
          dR = reco::deltaR( iJet->eta(),iJet->phi(), qcdOut->eta(),qcdOut->phi() );
          if ( debug ) std::cout << " >>>>>> genJet[" << k << "]Pt:" << iJet->pt() << " jetEta:" << iJet->eta() << " jetPhi:" << iJet->phi() << " dR:" << dR << std::endl; 
        } // genJets

        } // final
      } // outgoing partons
    } // incoming partons

  } // gen particles
    */

  /*
  if ( debug ) std::cout << " >> genJets.size: " << genJets->size() << std::endl;
  for ( unsigned int i = 0; i < genJets->size(); i++ ) {
    reco::GenJetRef iJet( genJets, i );
    if ( debug ) std::cout << " >> genJet[" << i << "]Pt:" << iJet->pt() << " jetEta:" << iJet->eta() << " jetPhi:" << iJet->phi() << std::endl; 
    unsigned int nConstituents = iJet->getGenConstituents().size();
    for ( unsigned int j = 0; j < nConstituents; j++ ) {
      const reco::GenParticle* subJet = iJet->getGenConstituent( j );
      //if ( std::abs(subJet->pdgId()) != 21 ) continue;
      dR = reco::deltaR( iJet->eta(),iJet->phi(), subJet->eta(),subJet->phi() );
      if ( debug ) {
        std::cout << " >> pdgId:" << subJet->pdgId() << " status:" << subJet->status()
        << " pT:" << subJet->pt() << " eta:" << subJet->eta() << " phi: " << subJet->phi() << " E:" << subJet->energy();
        std::cout << " moth:" << subJet->mother()->pdgId();
        std::cout << " dR: " << dR << std::endl;
      }
    }
  }
  if ( debug ) std::cout << " >> nQcdIn:" << nQcdIn << " nQcdOut:" << nQcdOut << " nGG:" << nGG << std::endl;
  if ( nGG != 2 ) return false;
  if ( nQcdIn !=2 || nQcdOut !=2 ) return false;
  */

  if ( debug ) std::cout << " >> has_dijet: passed" << std::endl;
  return true;

} // has_dijet
