#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"

// Run jet event selection ////////////////////////////////

const int search_window = 7;
//const int image_padding = 14;
const int image_padding = 12;
unsigned int jet_runId_;
unsigned int jet_lumiId_;
unsigned long long jet_eventId_;
std::vector<float> vJetSeed_iphi_;
std::vector<float> vJetSeed_ieta_;

// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel_jet ( TTree* tree, edm::Service<TFileService> &fs ) {

  tree->Branch("eventId",        &jet_eventId_);
  tree->Branch("runId",          &jet_runId_);
  tree->Branch("lumiId",         &jet_lumiId_);
  tree->Branch("jetSeed_iphi",   &vJetSeed_iphi_);
  tree->Branch("jetSeed_ieta",   &vJetSeed_ieta_);

  branchesEvtSel_dijet_gg_qq( tree, fs );

} // branchesEvtSel_jet()

// Run event selection ___________________________________________________________________//
bool RecHitAnalyzer::runEvtSel_jet ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<HBHERecHitCollection> HBHERecHitsH_;
  iEvent.getByToken( HBHERecHitCollectionT_, HBHERecHitsH_ );

  edm::ESHandle<CaloGeometry> caloGeomH_;
  iSetup.get<CaloGeometryRecord>().get( caloGeomH_ );
  const CaloGeometry* caloGeom = caloGeomH_.product();

  vGoodJetIdxs.clear();
  vJetSeed_iphi_.clear();
  vJetSeed_ieta_.clear();

  // Run individual selection
  bool hasPassed;
  hasPassed = runEvtSel_dijet_gg_qq( iEvent, iSetup );

  if ( !hasPassed ) return false; 

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);
  if ( debug ) std::cout << " >> PFJetCol.size: " << jets->size() << std::endl;

  HcalSubdetector subdet_;
  HcalDetId seedId;
  float seedE;
  int iphi_, ieta_, ietaAbs_;

  // Loop over jets
  for ( int thisJetIdx : vGoodJetIdxs ) {

    reco::PFJetRef iJet( jets, thisJetIdx );

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
    vJetSeed_iphi_.push_back( iphi_ );
    vJetSeed_ieta_.push_back( ieta_ );

  } // good jets 

  jet_eventId_ = iEvent.id().event();
  jet_runId_ = iEvent.id().run();
  jet_lumiId_ = iEvent.id().luminosityBlock();
  if ( debug ) std::cout << " >> analyze: passed" << std::endl;
  return true;

} // runEvtSel_jet()
