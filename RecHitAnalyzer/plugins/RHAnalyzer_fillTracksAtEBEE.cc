#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"

// Fill Tracks in EB+EE ////////////////////////////////
// Store tracks in EB+EE projection

TH2F *hTracks_EE[nEE];
TH2F *hTracks_EB;
TH2F *hTracksPt_EE[nEE];
TH2F *hTracksPt_EB;
std::vector<float> vTracksPt_EE_[nEE];
std::vector<float> vEndTracksPt_EE_[nEE];
std::vector<float> vTracksQPt_EE_[nEE];
std::vector<float> vTracks_EE_[nEE];
std::vector<float> vTracksPt_EB_;
std::vector<float> vEndTracksPt_EB_;
std::vector<float> vTracksQPt_EB_;
std::vector<float> vTracks_EB_;
TH1F *hTracks_minDr;
TH1F *hTracks_matchPt;
TH1F *hTracks_matchEta;
TH1F *hTracks_matchPhi;

TH1F *hTracks_NomatchPt;
TH1F *hTracks_NomatchEta;
TH1F *hTracks_NomatchPhi;

TH1F *hPFlow_minDr;
TH1F *hPFlow_matchPt;
TH1F *hPFlow_matchEta;
TH1F *hPFlow_matchPhi;
TH1F *hPFlow_matchID;

TH1F *hPFlow_NomatchPt;
TH1F *hPFlow_NomatchEta;
TH1F *hPFlow_NomatchPhi;
TH1F *hPFlow_NomatchID;


// Initialize branches ____________________________________________________________//
void RecHitAnalyzer::branchesTracksAtEBEE ( TTree* tree, edm::Service<TFileService> &fs ) {

  // Branches for images
  tree->Branch("Tracks_EB",   &vTracks_EB_);
  tree->Branch("TracksPt_EB", &vTracksPt_EB_);
  tree->Branch("EndTracksPt_EB", &vEndTracksPt_EB_);
  tree->Branch("TracksQPt_EB", &vTracksQPt_EB_);

  // Histograms for monitoring
  hTracks_EB = fs->make<TH2F>("Tracks_EB", "N(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX  , EB_IPHI_MIN-1, EB_IPHI_MAX,
      2*EB_IETA_MAX,-EB_IETA_MAX,   EB_IETA_MAX );
  hTracksPt_EB = fs->make<TH2F>("TracksPt_EB", "pT(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX  , EB_IPHI_MIN-1, EB_IPHI_MAX,
      2*EB_IETA_MAX,-EB_IETA_MAX,   EB_IETA_MAX );


  hTracks_minDr = fs->make<TH1F>("Tracks_minDr", "MinDr;dR;Tracks",100, -0.01,0.1);

  hTracks_matchPt  = fs->make<TH1F>("Tracks_matchPt",  "MatchPt;PT;Tracks",  100, -0.01,100);
  hTracks_matchEta = fs->make<TH1F>("Tracks_matchEta", "MatchEta;Eta;Tracks",100, -2.6,2.6);
  hTracks_matchPhi = fs->make<TH1F>("Tracks_matchPhi", "MatchPhi;Phi;Tracks",100, -3.2,3.2);

  hTracks_NomatchPt  = fs->make<TH1F>("Tracks_NomatchPt",  "NoMatchPt;PT;Tracks",  100, -0.01,100);
  hTracks_NomatchEta = fs->make<TH1F>("Tracks_NomatchEta", "NoMatchEta;Eta;Tracks",100, -2.6,2.6);
  hTracks_NomatchPhi = fs->make<TH1F>("Tracks_NomatchPhi", "NoMatchPhi;Phi;Tracks",100, -3.2,3.2);

  hPFlow_minDr = fs->make<TH1F>("PFlow_minDr", "MinDr;dR;PFlow",100, -0.01,0.1);

  hPFlow_matchPt  = fs->make<TH1F>("PFlow_matchPt",  "MatchPt;PT;PFlow",  100, -0.01,100);
  hPFlow_matchEta = fs->make<TH1F>("PFlow_matchEta", "MatchEta;Eta;PFlow",100, -2.6,2.6);
  hPFlow_matchPhi = fs->make<TH1F>("PFlow_matchPhi", "MatchPhi;Phi;PFlow",100, -3.2,3.2);
  hPFlow_matchID  = fs->make<TH1F>("PFlow_matchID",  "MatchID;ID;PFlow",  10,  -0.5,9.5);

  hPFlow_NomatchPt  = fs->make<TH1F>("PFlow_NomatchPt",  "NoMatchPt;PT;PFlow",  100, -0.01,100);
  hPFlow_NomatchEta = fs->make<TH1F>("PFlow_NomatchEta", "NoMatchEta;Eta;PFlow",100, -2.6,2.6);
  hPFlow_NomatchPhi = fs->make<TH1F>("PFlow_NomatchPhi", "NoMatchPhi;Phi;PFlow",100, -3.2,3.2);
  hPFlow_NomatchID  = fs->make<TH1F>("PFlow_NomatchID",  "NoMatchID;ID;PFlow",  10,  -0.5,9.5);


  char hname[50], htitle[50];
  for ( int iz(0); iz < nEE; iz++ ) {
    // Branches for images
    const char *zside = (iz > 0) ? "p" : "m";
    sprintf(hname, "Tracks_EE%s",zside);
    tree->Branch(hname,        &vTracks_EE_[iz]);
    sprintf(hname, "TracksPt_EE%s",zside);
    tree->Branch(hname,        &vTracksPt_EE_[iz]);
    sprintf(hname, "EndTracksPt_EE%s",zside);
    tree->Branch(hname,        &vEndTracksPt_EE_[iz]);
    sprintf(hname, "TracksQPt_EE%s",zside);
    tree->Branch(hname,        &vTracksQPt_EE_[iz]);

    // Histograms for monitoring
    sprintf(hname, "Tracks_EE%s",zside);
    sprintf(htitle,"N(ix,iy);ix;iy");
    hTracks_EE[iz] = fs->make<TH2F>(hname, htitle,
        EE_MAX_IX, EE_MIN_IX-1, EE_MAX_IX,
        EE_MAX_IY, EE_MIN_IY-1, EE_MAX_IY );
    sprintf(hname, "TracksPt_EE%s",zside);
    sprintf(htitle,"pT(ix,iy);ix;iy");
    hTracksPt_EE[iz] = fs->make<TH2F>(hname, htitle,
        EE_MAX_IX, EE_MIN_IX-1, EE_MAX_IX,
        EE_MAX_IY, EE_MIN_IY-1, EE_MAX_IY );
  } // iz

} // branchesEB()

// Fill TRK rechits at EB/EE ______________________________________________________________//
void RecHitAnalyzer::fillTracksAtEBEE ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  int ix_, iy_, iz_;
  int iphi_, ieta_, idx_; // rows:ieta, cols:iphi
  float eta, phi;
  GlobalPoint pos;

  vTracks_EB_.assign( EBDetId::kSizeForDenseIndexing, 0. );
  vTracksPt_EB_.assign( EBDetId::kSizeForDenseIndexing, 0. );
  vEndTracksPt_EB_.assign( EBDetId::kSizeForDenseIndexing, 0. );
  vTracksQPt_EB_.assign( EBDetId::kSizeForDenseIndexing, 0. );
  for ( int iz(0); iz < nEE; iz++ ) {
    vTracks_EE_[iz].assign( EE_NC_PER_ZSIDE, 0. );
    vTracksPt_EE_[iz].assign( EE_NC_PER_ZSIDE, 0. );
    vEndTracksPt_EE_[iz].assign( EE_NC_PER_ZSIDE, 0. );
    vTracksQPt_EE_[iz].assign( EE_NC_PER_ZSIDE, 0. );
  }

  edm::Handle<reco::TrackCollection> tracksH_;
  iEvent.getByToken( trackCollectionT_, tracksH_ );

  edm::Handle<PFCollection> pfCandsH_;
  iEvent.getByToken( pfCollectionT_, pfCandsH_ );

  // Provides access to global cell position
  edm::ESHandle<CaloGeometry> caloGeomH_;
  iSetup.get<CaloGeometryRecord>().get( caloGeomH_ );
  const CaloGeometry* caloGeom = caloGeomH_.product();

  reco::Track::TrackQuality tkQt_ = reco::Track::qualityByName("highPurity");

  for ( reco::TrackCollection::const_iterator iTk = tracksH_->begin();
        iTk != tracksH_->end(); ++iTk ) {
    if ( !(iTk->quality(tkQt_)) ) continue;

    eta = iTk->eta();
    phi = iTk->phi();
    if ( std::abs(eta) > 3. ) continue;

    float minDr = 10;
    const reco::PFCandidate* pfCand = getPFCand(pfCandsH_,eta,phi,minDr);
    hTracks_minDr->Fill(minDr);
    if(!pfCand){
      std::cout << "No PFCand: " << iTk->pt() << " " << eta << " " << phi << std::endl;
      getPFCand(pfCandsH_, eta, phi, minDr, true);
      hTracks_NomatchPt  -> Fill(iTk->pt());
      hTracks_NomatchEta -> Fill(iTk->eta());
      hTracks_NomatchPhi -> Fill(iTk->phi());
    }else{
      hTracks_matchPt  -> Fill(iTk->pt());
      hTracks_matchEta -> Fill(iTk->eta());
      hTracks_matchPhi -> Fill(iTk->phi());
    }


    DetId id( spr::findDetIdECAL( caloGeom, eta, phi, false ) );
    if ( id.subdetId() == EcalBarrel ) {
      EBDetId ebId( id );
      iphi_ = ebId.iphi() - 1;
      ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta();
      // Fill histograms for monitoring
      hTracks_EB->Fill( iphi_, ieta_ );
      hTracksPt_EB->Fill( iphi_, ieta_, iTk->pt() );
      idx_ = ebId.hashedIndex(); // (ieta_+EB_IETA_MAX)*EB_IPHI_MAX + iphi_
      // Fill vectors for images
      vTracks_EB_[idx_] += 1.;
      vTracksPt_EB_[idx_] += iTk->pt();
      vTracksQPt_EB_[idx_] += (iTk->charge()*iTk->pt());
    } else if ( id.subdetId() == EcalEndcap ) {
      EEDetId eeId( id );
      ix_ = eeId.ix() - 1;
      iy_ = eeId.iy() - 1;
      iz_ = (eeId.zside() > 0) ? 1 : 0;
      // Fill histograms for monitoring
      hTracks_EE[iz_]->Fill( ix_, iy_ );
      hTracksPt_EE[iz_]->Fill( ix_, iy_, iTk->pt() );
      // Create hashed Index: maps from [iy][ix] -> [idx_]
      idx_ = iy_*EE_MAX_IX + ix_;
      // Fill vectors for images
      vTracks_EE_[iz_][idx_] += 1.;
      vTracksPt_EE_[iz_][idx_] += iTk->pt();
      vTracksQPt_EE_[iz_][idx_] += (iTk->charge()*iTk->pt());
    } 
  } // tracks


  for ( PFCollection::const_iterator iPFC = pfCandsH_->begin();
        iPFC != pfCandsH_->end(); ++iPFC ) {
    const reco::Track* thisTrk = iPFC->bestTrack();
    if(!thisTrk) continue;

    //if ( !(thisTrk->quality(tkQt_)) ) continue;

    float minDr = 10;
    const reco::Track* trackCand = getTrackCand(tracksH_,thisTrk->eta(),thisTrk->phi(),minDr);
    hPFlow_minDr->Fill(minDr);
    if(!trackCand){
      std::cout << "No trackCand: " << thisTrk->pt() << " " << thisTrk->eta() << " " << thisTrk->phi() << std::endl;
      getTrackCand(tracksH_, thisTrk->eta(), thisTrk->phi(), minDr, true);
      hPFlow_NomatchPt  -> Fill(thisTrk->pt());
      hPFlow_NomatchEta -> Fill(thisTrk->eta());
      hPFlow_NomatchPhi -> Fill(thisTrk->phi());
      hPFlow_NomatchID -> Fill(iPFC->particleId());
    }else{
      hPFlow_matchPt  -> Fill(thisTrk->pt());
      hPFlow_matchEta -> Fill(thisTrk->eta());
      hPFlow_matchPhi -> Fill(thisTrk->phi());
      hPFlow_matchID -> Fill(iPFC->particleId());
    }

    const math::XYZPointF& ecalPos = iPFC->positionAtECALEntrance();
    eta = ecalPos.eta();
    phi = ecalPos.phi();
    
    if ( std::abs(eta) > 3. ) continue;
    DetId id( spr::findDetIdECAL( caloGeom, eta, phi, false ) );
    if ( id.subdetId() == EcalBarrel ) {
      EBDetId ebId( id );
      iphi_ = ebId.iphi() - 1;
      ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta();

      idx_ = ebId.hashedIndex(); // (ieta_+EB_IETA_MAX)*EB_IPHI_MAX + iphi_
      // Fill vectors for images
      vEndTracksPt_EB_[idx_] += thisTrk->pt();

    } else if ( id.subdetId() == EcalEndcap ) {
      EEDetId eeId( id );
      ix_ = eeId.ix() - 1;
      iy_ = eeId.iy() - 1;
      iz_ = (eeId.zside() > 0) ? 1 : 0;
        
      // Create hashed Index: maps from [iy][ix] -> [idx_]
      idx_ = iy_*EE_MAX_IX + ix_;
      // Fill vectors for images
      vEndTracksPt_EE_[iz_][idx_] += thisTrk->pt();
    } 
  }//PF Candidates

} // fillEB()
