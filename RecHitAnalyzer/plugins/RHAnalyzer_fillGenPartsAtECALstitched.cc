#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"
#include "Calibration/IsolatedParticles/interface/DetIdFromEtaPhi.h"

// Fill stitched EEm_EB_EEp image /////////////////////............/
// Store all ECAL event rechits into a stitched EEm_EB_EEp image 
// segmented in iphi,ieta spannning the full -3 < eta < 3. 
// Use EB-like granularity giving an extended range ieta=[-140,140].
//
// For endcaps, project EE hits into a helper histogram binned by 
// phi,eta before filling the full extended ECAL(iphi,eta) image.
// For barrel, fill EB hits directly since geometries are 1:1. 
//
// 'ieta_global' keeps track of the global ieta index count used
// for filling the extended image vector vECAL_genPartPt.
// 'ieta_signed' keeps track of the position along [-140,140] used
// for filling the monitoring histogram hECAL_genPartPt.

TH2F *hEvt_EE_genEMPartPt[nEE];
TH2F *hEvt_EE_genHadPartPt[nEE];
TProfile2D *hECAL_genEMPartPt;
TProfile2D *hECAL_genHadPartPt;
std::vector<float> vECAL_genEMPartPt_;
std::vector<float> vECAL_genHadPartPt_;

// Initialize branches _______________________________________________________________//
void RecHitAnalyzer::branchesGenPartsAtECALstitched ( TTree* tree, edm::Service<TFileService> &fs ) {

  // Branches for images
  tree->Branch("ECAL_genEMPartPt",     &vECAL_genEMPartPt_);
  tree->Branch("ECAL_genHadPartPt",    &vECAL_genHadPartPt_);
  // Intermediate helper histogram (single event only)
  hEvt_EE_genEMPartPt[0] = new TH2F("evt_EEm_genEMPartPt", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEm );
  hEvt_EE_genEMPartPt[1] = new TH2F("evt_EEp_genEMPartPt", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEp );
  hEvt_EE_genHadPartPt[0] = new TH2F("evt_EEm_genHadPartPt", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEm );
  hEvt_EE_genHadPartPt[1] = new TH2F("evt_EEp_genHadPartPt", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEp );

  // Histograms for monitoring
  hECAL_genEMPartPt = fs->make<TProfile2D>("ECAL_genEMPartPt", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX,    EB_IPHI_MIN-1, EB_IPHI_MAX,
      2*ECAL_IETA_MAX_EXT, -ECAL_IETA_MAX_EXT,   ECAL_IETA_MAX_EXT );
  hECAL_genHadPartPt = fs->make<TProfile2D>("ECAL_genHadPartPt", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX,    EB_IPHI_MIN-1, EB_IPHI_MAX,
      2*ECAL_IETA_MAX_EXT, -ECAL_IETA_MAX_EXT,   ECAL_IETA_MAX_EXT );

} // branchesGenPartsAtECALstitched()

// Function to map EE(phi,eta) histograms to ECAL(iphi,ieta) vector _______________________________//
void fillGenPartsAtECAL_with_EEproj ( TH2F *hEvt_EE_genPartPt_, int ieta_global_offset, int ieta_signed_offset, int genType ) {

  int ieta_global_, ieta_signed_;
  int ieta_, iphi_, idx_;
  float genPt_;

  for (int ieta = 1; ieta < hEvt_EE_genPartPt_->GetNbinsY()+1; ieta++) {
    ieta_ = ieta - 1;
    ieta_global_ = ieta_ + ieta_global_offset;
    ieta_signed_ = ieta_ + ieta_signed_offset;
    for (int iphi = 1; iphi < hEvt_EE_genPartPt_->GetNbinsX()+1; iphi++) {

      genPt_ = hEvt_EE_genPartPt_->GetBinContent( iphi, ieta );
      if ( genPt_ <= zs ) continue;
      // NOTE: EB iphi = 1 does not correspond to physical phi = -pi so need to shift!
      iphi_ = iphi  + 5*38; // shift
      iphi_ = iphi_ > EB_IPHI_MAX ? iphi_-EB_IPHI_MAX : iphi_; // wrap-around
      iphi_ = iphi_ - 1;
      idx_  = ieta_global_*EB_IPHI_MAX + iphi_;
      // EM
      if ( genType == 0 ) {
        // Fill vector for image
        vECAL_genEMPartPt_[idx_] = genPt_;
        // Fill histogram for monitoring
        hECAL_genEMPartPt->Fill( iphi_, ieta_signed_, genPt_ );
      } else {
      // Had
        // Fill vector for image
        vECAL_genHadPartPt_[idx_] = genPt_;
        // Fill histogram for monitoring
        hECAL_genHadPartPt->Fill( iphi_, ieta_signed_, genPt_ );
      }

    } // iphi_
  } // ieta_

} // fillGenPartsAtECAL_with_EEproj

// Fill stitched EE-, EB, EE+ rechits ________________________________________________________//
void RecHitAnalyzer::fillGenPartsAtECALstitched ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  int iphi_, ieta_, iz_, idx_;
  int ieta_global, ieta_signed;
  int ieta_global_offset, ieta_signed_offset;
  float eta, phi, genPt_;
  GlobalPoint pos;

  vECAL_genEMPartPt_.assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
  vECAL_genHadPartPt_.assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
  for ( int iz(0); iz < nEE; ++iz ) {
    hEvt_EE_genEMPartPt[iz]->Reset();
    hEvt_EE_genHadPartPt[iz]->Reset();
  }

  edm::Handle<EcalRecHitCollection> EBRecHitsH_;
  iEvent.getByLabel( EBRecHitCollectionT_, EBRecHitsH_ );
  edm::Handle<EcalRecHitCollection> EERecHitsH_;
  iEvent.getByLabel( EERecHitCollectionT_, EERecHitsH_ );
  edm::ESHandle<CaloGeometry> caloGeomH_;
  iSetup.get<CaloGeometryRecord>().get( caloGeomH_ );
  const CaloGeometry* caloGeom = caloGeomH_.product();

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel( genParticleCollectionT_, genParticles );

  // Fill endcap histogram binned in phi, eta
  for ( reco::GenParticleCollection::const_iterator iGen = genParticles->begin();
      iGen != genParticles->end();
      ++iGen ) {

    if ( iGen->status() != 1 ) continue;
    eta = iGen->eta();
    phi = iGen->phi();
    if ( std::abs(eta) > 3. ) continue;
    DetId id( spr::findDetIdECAL( caloGeom, eta, phi, false ) );
    if ( id.subdetId() == EcalBarrel ) continue;
    if ( id.subdetId() == EcalEndcap ) {
      iz_ = (eta > 0.) ? 1 : 0;
      // Fill intermediate helper histogram by eta,phi
      if ( abs(iGen->pdgId()) == 11 || abs(iGen->pdgId()) == 22 ) {
        hEvt_EE_genEMPartPt[iz_]->Fill( phi, eta, iGen->pt() );
      } else {
        hEvt_EE_genHadPartPt[iz_]->Fill( phi, eta, iGen->pt() );
      }
    }
  } // genParts 

  // Map EE-(phi,eta) to bottom part of ECAL(iphi,ieta)
  ieta_global_offset = 0;
  ieta_signed_offset = -ECAL_IETA_MAX_EXT;
  fillGenPartsAtECAL_with_EEproj( hEvt_EE_genEMPartPt[0], ieta_global_offset, ieta_signed_offset, 0 );
  fillGenPartsAtECAL_with_EEproj( hEvt_EE_genHadPartPt[0], ieta_global_offset, ieta_signed_offset, 1 );

  // Fill middle part of ECAL(iphi,ieta) with the EB rechits.
  ieta_global_offset = 55;

  // Fill endcap histogram binned in phi, eta
  for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin();
      iGen != genParticles->end();
      ++iGen) {

    if ( iGen->status() != 1 ) continue;
    eta = iGen->eta();
    phi = iGen->phi();
    genPt_ = iGen->pt();
    if ( std::abs(eta) > 3. ) continue;
    DetId id( spr::findDetIdECAL( caloGeom, eta, phi, false ) );
    if ( id.subdetId() == EcalEndcap ) continue;
    if ( id.subdetId() == EcalBarrel ) { 
      EBDetId ebId( id );
      iphi_ = ebId.iphi() - 1;
      ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta();
      if ( genPt_ <= zs ) continue;
      // Fill vector for image
      ieta_signed = ieta_;
      ieta_global = ieta_ + EB_IETA_MAX + ieta_global_offset;
      idx_ = ieta_global*EB_IPHI_MAX + iphi_; 
      if ( abs(iGen->pdgId()) == 11 || abs(iGen->pdgId()) == 22 ) {
        vECAL_genEMPartPt_[idx_] += genPt_;
        // Fill histogram for monitoring
        hECAL_genEMPartPt->Fill( iphi_, ieta_signed, genPt_ );
      } else {
        vECAL_genHadPartPt_[idx_] += genPt_;
        // Fill histogram for monitoring
        hECAL_genHadPartPt->Fill( iphi_, ieta_signed, genPt_ );
      }
    }

  } // EB

  // Map EE+(phi,eta) to upper part of ECAL(iphi,ieta)
  ieta_global_offset = ECAL_IETA_MAX_EXT + EB_IETA_MAX;
  ieta_signed_offset = EB_IETA_MAX;
  fillGenPartsAtECAL_with_EEproj( hEvt_EE_genEMPartPt[1], ieta_global_offset, ieta_signed_offset, 0 );
  fillGenPartsAtECAL_with_EEproj( hEvt_EE_genHadPartPt[1], ieta_global_offset, ieta_signed_offset, 1 );

} // fillGenPartsAtECALstitched()
