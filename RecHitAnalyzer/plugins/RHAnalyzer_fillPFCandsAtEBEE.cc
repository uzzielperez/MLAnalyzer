#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"

// Fill PFCands in EB+EE ////////////////////////////////
// Store PFCands in EB+EE projection

std::vector<float> vEndTracksPt_EE_[nEE];
std::vector<float> vMuonsPt_EE_[nEE];
std::vector<float> vEndTracksPt_EB_;
std::vector<float> vMuonsPt_EB_;

// Initialize branches ____________________________________________________________//
void RecHitAnalyzer::branchesPFCandsAtEBEE ( TTree* tree, edm::Service<TFileService> &fs ) {

  // Branches for images
  tree->Branch("EndTracksPt_EB", &vEndTracksPt_EB_);
  tree->Branch("MuonsPt_EB", &vMuonsPt_EB_);

  char hname[50];
  for ( int iz(0); iz < nEE; iz++ ) {
    // Branches for images
    const char *zside = (iz > 0) ? "p" : "m";
    sprintf(hname, "EndTracksPt_EE%s",zside);
    tree->Branch(hname,        &vEndTracksPt_EE_[iz]);
    sprintf(hname, "MuonsPt_EE%s",zside);
    tree->Branch(hname,        &vMuonsPt_EE_[iz]);
  } // iz

} // branchesEB()

// Fill TRK rechits at EB/EE ______________________________________________________________//
void RecHitAnalyzer::fillPFCandsAtEBEE ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  int ix_, iy_, iz_;
  int idx_; // rows:ieta, cols:iphi
  float eta, phi;
  GlobalPoint pos;

  vEndTracksPt_EB_.assign( EBDetId::kSizeForDenseIndexing, 0. );
  vMuonsPt_EB_.assign( EBDetId::kSizeForDenseIndexing, 0. );
  for ( int iz(0); iz < nEE; iz++ ) {
    vEndTracksPt_EE_[iz].assign( EE_NC_PER_ZSIDE, 0. );
    vMuonsPt_EE_[iz].assign( EE_NC_PER_ZSIDE, 0. );
  }

  edm::Handle<PFCollection> pfCandsH_;
  iEvent.getByToken( pfCollectionT_, pfCandsH_ );

  // Provides access to global cell position
  edm::ESHandle<CaloGeometry> caloGeomH_;
  iSetup.get<CaloGeometryRecord>().get( caloGeomH_ );
  const CaloGeometry* caloGeom = caloGeomH_.product();


  for ( PFCollection::const_iterator iPFC = pfCandsH_->begin();
        iPFC != pfCandsH_->end(); ++iPFC ) {
    const reco::Track* thisTrk = iPFC->bestTrack();
    if(!thisTrk) continue;

    const math::XYZPointF& ecalPos = iPFC->positionAtECALEntrance();
    eta = ecalPos.eta();
    phi = ecalPos.phi();
    
    if ( std::abs(eta) > 3. ) continue;
    DetId id( spr::findDetIdECAL( caloGeom, eta, phi, false ) );
    if ( id.subdetId() == EcalBarrel ) {
      EBDetId ebId( id );

      idx_ = ebId.hashedIndex(); // (ieta_+EB_IETA_MAX)*EB_IPHI_MAX + iphi_
      // Fill vectors for images
      
      vEndTracksPt_EB_[idx_] += thisTrk->pt();
      if(iPFC->particleId() == 3)
	vMuonsPt_EB_[idx_] += thisTrk->pt();

    } else if ( id.subdetId() == EcalEndcap ) {
      EEDetId eeId( id );
      ix_ = eeId.ix() - 1;
      iy_ = eeId.iy() - 1;
      iz_ = (eeId.zside() > 0) ? 1 : 0;
        
      // Create hashed Index: maps from [iy][ix] -> [idx_]
      idx_ = iy_*EE_MAX_IX + ix_;
      // Fill vectors for images
      vEndTracksPt_EE_[iz_][idx_] += thisTrk->pt();
      if(iPFC->particleId() == 3)
	vMuonsPt_EE_[iz_][idx_] += thisTrk->pt();
    } 
  }//PF Candidates

} // fillEB()
