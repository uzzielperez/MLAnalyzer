#include "MLAnalyzer/RecHitAnalyzer/interface/SCRegressor.h"

// Initialize branches _____________________________________________________//
void SCRegressor::branchesSC ( TTree* tree, edm::Service<TFileService> &fs )
{
  char hname[50], htitle[50];
  for(int iPho (0); iPho < nPhotons; iPho++) {
    sprintf(hname, "SC_energy%d",iPho);
    sprintf(htitle,"E(i#phi,i#eta);iphi;ieta");
    hSC_energy[iPho] = fs->make<TH2D>(hname, htitle,
        crop_size, 0, crop_size,
        crop_size, 0, crop_size );
    sprintf(hname, "SC_time%d",iPho);
    sprintf(htitle,"t(i#phi,i#eta);iphi;ieta");
    hSC_time[iPho] = fs->make<TH2D>(hname, htitle,
        crop_size, 0, crop_size,
        crop_size, 0, crop_size );
  }
  for(int iPho (0); iPho < nPhotons; iPho++) {
    sprintf(hname, "SC_energy%d",iPho);
    RHTree->Branch(hname,          &vSC_energy_[iPho]);
    sprintf(hname, "SC_energyT%d",iPho);
    RHTree->Branch(hname,          &vSC_energyT_[iPho]);
    sprintf(hname, "SC_energyZ%d",iPho);
    RHTree->Branch(hname,          &vSC_energyZ_[iPho]);
    sprintf(hname, "SC_time%d",iPho);
    RHTree->Branch(hname,          &vSC_time_[iPho]);
  }
} // branchesSC()

// Fill SC rechits _________________________________________________________________//
void SCRegressor::fillSC ( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{

  edm::Handle<EcalRecHitCollection> EBRecHitsH;
  iEvent.getByToken(EBRecHitCollectionT_, EBRecHitsH);

  edm::ESHandle<CaloGeometry> caloGeomH;
  iSetup.get<CaloGeometryRecord>().get(caloGeomH);
  const CaloGeometry* caloGeom = caloGeomH.product();

  for(int iPho (0); iPho < nPhotons; iPho++) {
    vSC_energy_[iPho].assign(crop_size*crop_size,0.);
    vSC_energyT_[iPho].assign(crop_size*crop_size,0.);
    vSC_energyZ_[iPho].assign(crop_size*crop_size,0.);
    vSC_time_[iPho].assign(crop_size*crop_size,0.);
  }

  int iphi_, ieta_, idx_; // rows:ieta, cols:iphi
  int iphi_shift, ieta_shift;
  int iphi_crop, ieta_crop;
  for(EcalRecHitCollection::const_iterator iRHit = EBRecHitsH->begin();
      iRHit != EBRecHitsH->end();
      ++iRHit) {

    if ( iRHit->energy() < zs ) continue;

    // Convert detector coordinates to ordinals
    EBDetId ebId( iRHit->id() );
    iphi_ = ebId.iphi()-1; // [0,...,359]
    ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta(); // [-85,...,-1,0,...,84]
    ieta_ += EBDetId::MAX_IETA; // [0,...,169]

    for ( unsigned iP(0); iP < nPhotons; iP++ ) {

      iphi_shift = vIphi_Emax[iP] - 15;
      ieta_shift = vIeta_Emax[iP] - 15;
      //std::cout << " >> Storing: iphi_Emax,ieta_Emax: " << vIphi_Emax[iP] << ", " << vIeta_Emax[iP] << std::endl;

      // Convert to [0,...,31][0,...,31]
      ieta_crop = ieta_ - ieta_shift;
      iphi_crop = iphi_ - iphi_shift;
      if ( iphi_crop >= EBDetId::MAX_IPHI ) iphi_crop = iphi_crop - EBDetId::MAX_IPHI; // get wrap-around hits

      if ( ieta_crop < 0 || ieta_crop > crop_size-1 ) continue;
      if ( iphi_crop < 0 || iphi_crop > crop_size-1 ) continue;
      
      // Convert to [0,...,32*32-1] 
      idx_ = ieta_crop*crop_size + iphi_crop;

      // Cell geometry provides access to (rho,eta,phi) coordinates of cell center
      //auto cell = caloGeom->getGeometry(ebId);
      auto pos = caloGeom->getPosition(ebId);
      
      // Fill branch arrays 
      vSC_energy_[iP][idx_] = iRHit->energy();
      //vSC_energyT_[iP][idx_] = iRHit->energy()/TMath::CosH(vPho_eta_[iP]);
      vSC_energyT_[iP][idx_] = iRHit->energy()/TMath::CosH(pos.eta());
      vSC_energyZ_[iP][idx_] = iRHit->energy()*std::abs(TMath::TanH(pos.eta()));
      vSC_time_[iP][idx_] = iRHit->time();
      /*
      vEB_SCenergy_[ebId.hashedIndex()] = iRHit->energy();
      */
      //std::cout << " >> " << iP << ": iphi_,ieta_,E: " << iphi_crop << ", " << ieta_crop << ", " << iRHit->energy() << std::endl; 

      // Fill histograms to monitor cumulative distributions
      hSC_energy[iP]->Fill( iphi_crop,ieta_crop,iRHit->energy() );
      hSC_time[iP]->Fill( iphi_crop,ieta_crop,iRHit->time() );

    } // photons
  } // EB rechits

} // fillSC()
