#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"
/*
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
*/
using std::vector;
// Run jet event selection ////////////////////////////////

TH1D *h_jet_pT;
TH1D *h_jet_E;
TH1D *h_jet_eta;
TH1D *h_jet_m0;
TH1D *h_jet_nJet;
float jet_eventId_;
float jet_m0_;
//const int nJets = 2;
const int search_window = 7;
//const int image_padding = 14;
const int image_padding = 12;
vector<float> vJetSeed_iphi_;
vector<float> vJetSeed_ieta_;

// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel_jet ( TTree* tree, edm::Service<TFileService> &fs ) {

  h_jet_pT    = fs->make<TH1D>("h_jet_pT"  , "p_{T};p_{T};Particles", 100,  0., 500.);
  h_jet_E     = fs->make<TH1D>("h_jet_E"   , "E;E;Particles"        , 100,  0., 800.);
  h_jet_eta   = fs->make<TH1D>("h_jet_eta" , "#eta;#eta;Particles"  , 100, -5., 5.);
  h_jet_nJet  = fs->make<TH1D>("h_jet_nJet", "nJet;nJet;Events"     ,  10,  0., 10.);
  h_jet_m0    = fs->make<TH1D>("h_jet_m0"  , "m0;m0;Events"         ,  60, 50., 110.);

  RHTree->Branch("eventId",        &jet_eventId_);
  RHTree->Branch("m0",             &jet_m0_);
  RHTree->Branch("jetSeed_iphi",   &vJetSeed_iphi_);
  RHTree->Branch("jetSeed_ieta",   &vJetSeed_ieta_);

} // branchesEvtSel_jet()

// Run event selection ___________________________________________________________________//
bool RecHitAnalyzer::runEvtSel_jet ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<HBHERecHitCollection> HBHERecHitsH_;
  iEvent.getByToken( HBHERecHitCollectionT_, HBHERecHitsH_ );

  edm::ESHandle<CaloGeometry> caloGeomH_;
  iSetup.get<CaloGeometryRecord>().get( caloGeomH_ );
  const CaloGeometry* caloGeom = caloGeomH_.product();

  /*
  edm::ESHandle<CaloTopology> caloTopoH_;
  iSetup.get<CaloTopologyRecord>().get( caloTopoH_ );
  const CaloTopology *caloTopo = caloTopoH_.product();
  */

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);
  if ( debug ) std::cout << " >> PFJetCol.size: " << jets->size() << std::endl;


  int nJet = 0;
  vector<int>       jetIdx;
  vJetSeed_iphi_.clear();
  vJetSeed_ieta_.clear();

  // Loop over jets
  for ( unsigned iJ(0); iJ != jets->size(); ++iJ ) {

    reco::PFJetRef iJet( jets, iJ );
    if ( std::abs(iJet->pt()) < minJetPt_ ) continue;
    if ( std::abs(iJet->eta()) > maxJetEta_) continue;
    //if ( iJet->mass() < 50. || iJet->mass() > 110. ) continue;
    if ( debug ) std::cout << " >> jet[" << iJ << "]Pt:" << iJet->pt() << " jetE:" << iJet->energy() << " jetM:" << iJet->mass() << std::endl;

    // Get closest HBHE tower to jet position
    // This will not always be the most energetic deposit
    HcalDetId hId( spr::findDetIdHCAL( caloGeom, iJet->eta(), iJet->phi(), false ) );
    if ( hId.subdet() != HcalBarrel && hId.subdet() != HcalEndcap ) continue;
    HBHERecHitCollection::const_iterator iRHit( HBHERecHitsH_->find(hId) );
    float seedE = ( iRHit == HBHERecHitsH_->end() ) ? 0. : iRHit->energy() ;
    HcalDetId seedId = hId;
    int iphi_, ieta_, ietaAbs_;
    if ( debug ) std::cout << " >> hId.ieta:" << hId.ieta() << " hId.iphi:" << hId.iphi() << " E:" << seedE << std::endl;

    // Look for the highest HBHE tower deposit within a search window
    for ( int ieta = 0; ieta < search_window; ieta++ ) {
      ieta_   = hId.ieta() - (search_window/2)+ieta;
      if ( std::abs(ieta_) > HBHE_IETA_MAX_HE-1 ) continue;
      if ( std::abs(ieta_) < HBHE_IETA_MIN_HB ) continue;
      HcalSubdetector subdet_ = std::abs(ieta_) > HBHE_IETA_MAX_HB ? HcalEndcap : HcalBarrel;
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
    if ( HBHE_IETA_MAX_HE-1 - ietaAbs_ < image_padding ) {
      if ( debug ) std::cout << " Fail HE edge cut " << std::endl;
      continue;
    }


    if ( debug ) std::cout << " !! ieta_:" << ieta_ << " iphi_:" << iphi_ << " ietaAbs_:" << ietaAbs_ << " E:" << seedE << std::endl;
    
    jetIdx        .push_back(iJ);
    vJetSeed_iphi_.push_back(iphi_);
    vJetSeed_ieta_.push_back(ieta_);

    nJet++;
    if ( (nJets_ > 0) && (nJet >= nJets_) ) break;

    // Doesnt seem to work:
    //const CaloSubdetectorTopology* topo = caloTopo->getSubdetectorTopology( hId );
    //const CaloSubdetectorTopology* topo = caloTopo->getSubdetectorTopology( DetId::Hcal, hId.subdet() );
    //std::vector<DetId> hId_5x5 = topo->getWindow( hId, 5, 5 );
    //std::cout << "window size:" << hId_5x5.size() << std::endl;
    //std::cout << " pT:" << iJet->pt() << " eta:" << iJet->eta() << " phi: " << iJet->phi() << " E:" << iJet->energy() << std::endl;
  } // jets

  h_jet_nJet->Fill( nJet );
  if ( debug ) {
    for(int thisJetIdx : jetIdx)
      std::cout << " >> jetIdx:" << thisJetIdx << std::endl;
  }

  if ( (nJets_ > 0) && (nJet != nJets_) ){
    if ( debug ) std::cout << " Fail jet multiplicity:  " << nJet << " != " << nJets_ << std::endl;
    return false;
  }

  if ( jetIdx.size() == 0){
    if ( debug ) std::cout << " No passing jets...  " << std::endl;
    return false;
  }


  // Plot kinematics
  for(int thisJetIdx : jetIdx){
    reco::PFJetRef thisJet( jets, thisJetIdx );
    if ( debug ) std::cout << " >> Jet[" << thisJetIdx << "] Pt:" << thisJet->pt() << std::endl;
    h_jet_pT->Fill( std::abs(thisJet->pt()) );
    h_jet_eta->Fill( thisJet->eta() );
    h_jet_E->Fill( thisJet->energy() );
    h_jet_m0->Fill( thisJet->mass() );

  }

  // Write out event ID
  reco::PFJetRef leadJet( jets, jetIdx.at(0) );
  jet_m0_ = leadJet->mass();
  jet_eventId_ = iEvent.id().event();

  return true;

} // runEvtSel_jet()
