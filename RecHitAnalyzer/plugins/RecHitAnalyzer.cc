// -*- C++ -*-
//
// Package:    MLAnalyzer/RecHitAnalyzer
// Class:      RecHitAnalyzer
// 
//
// Original Author:  Michael Andrews
//         Created:  Sat, 14 Jan 2017 17:45:54 GMT
//
//
#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"

//
// constructors and destructor
//
RecHitAnalyzer::RecHitAnalyzer(const edm::ParameterSet& iConfig)
{
  //EBRecHitCollectionT_ = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EBRecHitCollection"));
  EBRecHitCollectionT_ = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedEBRecHitCollection"));
  //EBDigiCollectionT_ = consumes<EBDigiCollection>(iConfig.getParameter<edm::InputTag>("selectedEBDigiCollection"));
  //EBDigiCollectionT_ = consumes<EBDigiCollection>(iConfig.getParameter<edm::InputTag>("EBDigiCollection"));
  EERecHitCollectionT_ = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedEERecHitCollection"));
  //EERecHitCollectionT_ = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EERecHitCollection"));
  HBHERecHitCollectionT_ = consumes<HBHERecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedHBHERecHitCollection"));

  genParticleCollectionT_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticleCollection"));
  photonCollectionT_ = consumes<reco::PhotonCollection>(iConfig.getParameter<edm::InputTag>("gedPhotonCollection"));

  //now do what ever initialization is needed
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  // Histograms
  // EB rechits
  hEB_energy = fs->make<TH2D>("EB_energy", "E(i#phi,i#eta);i#phi;i#eta",
      EBDetId::MAX_IPHI  , EBDetId::MIN_IPHI-1, EBDetId::MAX_IPHI,
      2*EBDetId::MAX_IETA,-EBDetId::MAX_IETA,   EBDetId::MAX_IETA );
  //hEB_energy = fs->make<TH2D>("EB_energy", "E(i#phi,i#eta);i#phi;i#eta",
  //    EcalTrigTowerDetId::kEBTowersInPhi*18  , EBDetId::MIN_IPHI-1, EBDetId::MAX_IPHI,
  //    EcalTrigTowerDetId::kEBTowersInEta*2,-EBDetId::MAX_IETA,   EBDetId::MAX_IETA );
  hEB_time = fs->make<TH2D>("EB_time", "t(i#phi,i#eta);i#phi;i#eta",
      EBDetId::MAX_IPHI  , EBDetId::MIN_IPHI-1, EBDetId::MAX_IPHI,
      2*EBDetId::MAX_IETA,-EBDetId::MAX_IETA,   EBDetId::MAX_IETA );
  // EB Digis
  char hname[50], htitle[50];
  /*
  for(int iS(0); iS < EcalDataFrame::MAXSAMPLES; iS++){
    sprintf(hname, "EB_adc%d",iS);
    sprintf(htitle,"adc%d(i#phi,i#eta);i#phi;i#eta",iS);
    hEB_adc[iS] = fs->make<TH2D>(hname, htitle,
        EBDetId::MAX_IPHI  , EBDetId::MIN_IPHI-1, EBDetId::MAX_IPHI,
        2*EBDetId::MAX_IETA,-EBDetId::MAX_IETA,   EBDetId::MAX_IETA );
  }
  */

  // EE rechits
  for(int iz(0); iz < EE_IZ_MAX; iz++){
    const char *zside = (iz > 0) ? "p" : "m";
    sprintf(hname, "EE%s_energy",zside);
    sprintf(htitle,"E(ix,iy);ix;iy");
    hEE_energy[iz] = fs->make<TH2D>(hname, htitle,
        EEDetId::IX_MAX, EEDetId::IX_MIN-1, EEDetId::IX_MAX,
        EEDetId::IY_MAX, EEDetId::IY_MIN-1, EEDetId::IY_MAX );
    sprintf(hname, "EE%s_time",zside);
    sprintf(htitle,"t(ix,iy);ix;iy");
    hEE_time[iz] = fs->make<TH2D>(hname, htitle,
        EEDetId::IX_MAX, EEDetId::IX_MIN-1, EEDetId::IX_MAX,
        EEDetId::IY_MAX, EEDetId::IY_MIN-1, EEDetId::IY_MAX );
  } 

  // HBHE
  hHBHE_energy = fs->make<TH2D>("HBHE_energy", "E(i#phi,i#eta);i#phi;i#eta",
      hcaldqm::constants::IPHI_NUM,      hcaldqm::constants::IPHI_MIN-1, hcaldqm::constants::IPHI_MAX,
      2*hcaldqm::constants::IETA_MAX_HE,-hcaldqm::constants::IETA_MAX_HE,hcaldqm::constants::IETA_MAX_HE );
  hHBHE_energy_EB = fs->make<TH2D>("HBHE_energy_EB", "E(i#phi,i#eta);i#phi;i#eta",
      hcaldqm::constants::IPHI_NUM, hcaldqm::constants::IPHI_MIN-1,hcaldqm::constants::IPHI_MAX,
      2*HBHE_IETA_MAX_EB,          -HBHE_IETA_MAX_EB,              HBHE_IETA_MAX_EB );
  //hHBHE_depth = fs->make<TH1D>("HBHE_depth", "Depth;depth;Hits",
  //    hcaldqm::constants::DEPTH_NUM, hcaldqm::constants::DEPTH_MIN, hcaldqm::constants::DEPTH_MAX+1);

  /*
  hEvt_HBHE_energy = fs->make<TH2D>("evt_HBHE_energy", "E(i#phi,i#eta);i#phi;i#eta",
      hcaldqm::constants::IPHI_NUM,      hcaldqm::constants::IPHI_MIN-1, hcaldqm::constants::IPHI_MAX,
      2*(hcaldqm::constants::IETA_MAX_HE-1),-(hcaldqm::constants::IETA_MAX_HE-1),hcaldqm::constants::IETA_MAX_HE-1 );
  hEvt_HBHE_EMenergy = fs->make<TH2D>("evt_HBHE_EMenergy", "E(i#phi,i#eta);i#phi;i#eta",
      hcaldqm::constants::IPHI_NUM, -TMath::Pi(), TMath::Pi(), 2*(hcaldqm::constants::IETA_MAX_HE-1), eta_bins_HBHE );
  */

  // Kinematics
  h_pT  = fs->make<TH1D>("h_pT" , "p_{T};p_{T};Particles", 100,  0., 500.);
  h_E   = fs->make<TH1D>("h_E"  , "E;E;Particles"        , 100,  0., 800.);
  h_eta = fs->make<TH1D>("h_eta", "#eta;#eta;Particles"  ,  50,-2.4, 2.4);
  h_m0  = fs->make<TH1D>("h_m0" , "m0;m0;Events"        ,   72, 50., 950.);

  // Output Tree
  // These will be use to create the actual images later on
  RHTree = fs->make<TTree>("RHTree", "RecHit tree");
  RHTree->Branch("eventId",      &eventId_);
  RHTree->Branch("m0",           &m0_);
  RHTree->Branch("ECAL_energy",  &vECAL_energy_);
  RHTree->Branch("EB_energy",    &vEB_energy_);
  RHTree->Branch("EB_time",      &vEB_time_);
  /*
     for(int iS(0); iS < EcalDataFrame::MAXSAMPLES; iS++){
     sprintf(hname, "EB_adc%d",iS);
     RHTree->Branch(hname,       &vEB_adc_[iS]);
     }
     */
  for(int iz(0); iz < EE_IZ_MAX; iz++){
    const char *zside = (iz > 0) ? "p" : "m";
    sprintf(hname, "EE%s_energy",zside);
    RHTree->Branch(hname,       &vEE_energy_[iz]);
    sprintf(hname, "EE%s_time",zside);
    RHTree->Branch(hname,       &vEE_time_[iz]);
  }
  RHTree->Branch("HBHE_energy_EB",     &vHBHE_energy_EB_);
  RHTree->Branch("HBHE_energy",        &vHBHE_energy_);
  RHTree->Branch("HBHE_EMenergy",      &vHBHE_EMenergy_);
}

//
// member functions
//
// ------------ method called for each event  ------------
void
RecHitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;

  // ----- Apply event selection cuts ----- //

  edm::Handle<reco::PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleCollectionT_, genParticles);

  bool passedSelection = false;
  //passedSelection = runSelections( photons, genParticles );
  passedSelection = runSelections_H2GG( photons, genParticles, true );

  if ( !passedSelection ) return;

  // ----- Get Calorimeter Geometry ----- //

  // Provides access to global cell position and coordinates below
  edm::ESHandle<CaloGeometry> caloGeomH;
  iSetup.get<CaloGeometryRecord>().get(caloGeomH);
  caloGeom = caloGeomH.product();
  //towerGeom = caloGeomH.product(); 
  //towerGeom->getSubdetectorGeometry(DetId::Calo, CaloTowerDetId::SubdetId);

  //////////// EB //////////

  // EB reduced rechit collection //
  // This contatins the reduced EB rechit collection after
  // the zero suppression and bad channel clean-up
  edm::Handle<EcalRecHitCollection> EBRecHitsH;
  iEvent.getByToken(EBRecHitCollectionT_, EBRecHitsH);
  TH2D *hEvt_HBHE_EMenergy = new TH2D("evt_HBHE_EMenergy", "E(i#phi,i#eta);i#phi;i#eta",
      hcaldqm::constants::IPHI_NUM, -TMath::Pi(), TMath::Pi(),
      2*(hcaldqm::constants::IETA_MAX_HE-1), eta_bins_HBHE );

  fillEBrechits( EBRecHitsH, caloGeom, hEvt_HBHE_EMenergy );

  // EB digis //
  // This contains the raw EB digi collection:
  // Each digi is unpacked as a data frame containing info
  // from 10 time samples [iS] of the pulse shape
  // iS=[0-2]: Presample noise
  // iS=[3-9]: Nominal pulse shape
  // NOTE: This is the raw collection and includes
  // selective-readout and bad channel effects!
  /*
  edm::Handle<EBDigiCollection> EBDigisH;
  iEvent.getByToken(EBDigiCollectionT_, EBDigisH);
  fillEBdigis( EBDigisH, caloGeom );
  */

  //////////// EE //////////

  // EE reduced rechit collection //
  // This contatins the reduced EE rechit collection after
  // the zero suppression and bad channel clean-up
  edm::Handle<EcalRecHitCollection> EERecHitsH;
  iEvent.getByToken(EERecHitCollectionT_, EERecHitsH);
  TH2D *hEvt_EE_energy[EE_IZ_MAX];
  hEvt_EE_energy[0] = new TH2D("evt_EEm_energy", "E(i#phi,i#eta);i#phi;i#eta",
      EBDetId::MAX_IPHI, -TMath::Pi(), TMath::Pi(),
      5*(hcaldqm::constants::IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEm );
  hEvt_EE_energy[1] = new TH2D("evt_EEp_energy", "E(i#phi,i#eta);i#phi;i#eta",
      EBDetId::MAX_IPHI, -TMath::Pi(), TMath::Pi(),
      5*(hcaldqm::constants::IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEp );

  fillEErechits( EERecHitsH, caloGeom, hEvt_HBHE_EMenergy, hEvt_EE_energy );

  //////////// ECAL //////////

  // ECAL @HCAL granularity
  fillECALatHCAL( hEvt_HBHE_EMenergy );

  // EB stitched to EEs(ieta,iphi)
  fillECALstitched( hEvt_EE_energy );

  //////////// HBHE //////////

  // HBHE reduced rechit collection //
  edm::Handle<HBHERecHitCollection> HBHERecHitsH;
  iEvent.getByToken( HBHERecHitCollectionT_, HBHERecHitsH );

  fillHBHErechits( HBHERecHitsH, caloGeom );

  //////////// Bookkeeping //////////

  // Write out event ID
  eventId_ = iEvent.id().event();

  // Fill RHTree
  RHTree->Fill();

} // analyze()


// ------------ method called once each job just before starting event loop  ------------
void 
RecHitAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
RecHitAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
RecHitAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//____ Apply event selection cuts _____//
bool RecHitAnalyzer::runSelections_H2GG ( edm::Handle<reco::PhotonCollection> &photons, edm::Handle<reco::GenParticleCollection> &genParticles, bool isHiggs ) {

  int nPho = 0;
  bool isDecayed = true;
  float etaCut = 1.4;
  //float etaCut = 2.3;
  float ptCut = 25.;
  float dRCut = 0.4;
  float dR, dEta, dPhi;
  std::cout << "PhoCol.size: " << photons->size() << std::endl;
  math::XYZTLorentzVector vDiPho;
  for(reco::PhotonCollection::const_iterator iPho = photons->begin();
      iPho != photons->end();
      ++iPho) {

      // Kinematic cuts
      if ( std::abs(iPho->eta()) > etaCut ) continue;
      if ( std::abs(iPho->pt()) < ptCut ) continue;

      for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin();
           iGen != genParticles->end();
           ++iGen) {

          // ID cuts
          if ( std::abs(iGen->pdgId()) != 22 ) continue;
          if ( iGen->status() != 1 ) continue; // NOT the same as Pythia status
          if ( !iGen->mother() ) continue;
          if ( isHiggs ) {
              if ( std::abs(iGen->mother()->pdgId()) != 25 && std::abs(iGen->mother()->pdgId()) != 22 ) continue;
          } else {
              if ( std::abs(iGen->mother()->status()) != 44 && std::abs(iGen->mother()->status()) != 23 ) continue;
          }

          // Kinematic cuts
          if ( std::abs(iGen->eta()) > etaCut ) continue;
          if ( std::abs(iGen->pt()) < ptCut ) continue;

          // Match by dR
          dEta = std::abs( iPho->eta() - iGen->eta() );
          dPhi = std::abs( iPho->phi() - iGen->phi() );
          dR = TMath::Power(dEta,2.) + TMath::Power(dPhi,2.);
          dR = TMath::Sqrt(dR);
          if ( dR < dRCut ) {
             nPho++;
             vDiPho += iPho->p4();
             break;
          }

      } // genParticle loop: count good photons

  } // recoPhotons

  // Require exactly 2 gen-level photons
  // Indifferent about photons of status != 1
  if ( nPho != 2 ) return false;
  //if ( vDiPho.mass() < 80. ) return false;

  // Fill loop
  for (reco::GenParticleCollection::const_iterator iP = genParticles->begin();
       iP != genParticles->end();
       ++iP) {

    // PDG ID cut
    if ( std::abs(iP->pdgId()) != 22 ) continue;

    // Decay status
    //if ( iP->status() != 23 ) continue;
    if ( iP->status() != 1 ) continue;

    if ( isDecayed ) {
        // Check ancestry
        if ( isHiggs ) {
            if ( std::abs(iP->mother()->pdgId()) != 25 && std::abs(iP->mother()->pdgId()) != 22 ) continue;
        } else {
            if ( std::abs(iP->mother()->status()) != 44 && std::abs(iP->mother()->status()) != 23 ) continue;
        }
        // Kinematic cuts
        if ( std::abs(iP->eta()) > etaCut ) continue;
        if ( std::abs(iP->pt()) < ptCut ) continue;
    } // apply cuts

    if ( isDecayed ) {
        std::cout << "status:" <<iP->status() << " pT:" << iP->pt() << " eta:" << iP->eta() << " E:" << iP->energy() << " mothId:" << iP->mother()->pdgId() << std::endl;
    } else {
        std::cout << "status:" <<iP->status() << " pT:" << iP->pt() << " eta:" << iP->eta() << " E:" << iP->energy() << std::endl;
    }
    // Fill histograms
    h_pT-> Fill( iP->pt()      );
    h_E->  Fill( iP->energy()  );
    h_eta->Fill( iP->eta()     );
  } // genParticle loop: fill hist
  h_m0->Fill( vDiPho.mass() );
  std::cout << " m0: " << vDiPho.mass() <<" (" << vDiPho.T() << ")" << std::endl;

  m0_ = vDiPho.mass();

  return true;

} // runSelections_H2GG()

//____ Apply event selection cuts _____//
bool RecHitAnalyzer::runSelections ( edm::Handle<reco::PhotonCollection> &photons, edm::Handle<reco::GenParticleCollection> &genParticles ) {

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

  /* 
  for(reco::PhotonCollection::const_iterator iPho = photons->begin();
      iPho != photons->end();
      ++iPho) {

    // Kinematic cuts
    if ( std::abs(iPho->eta()) > etaCut ) continue;
    if ( std::abs(iPho->pt()) < ptCut ) continue;

    h_pT-> Fill( iPho->pt()      );
    h_E->  Fill( iPho->energy()  );
    h_eta->Fill( iPho->eta()     );
    std::cout << "nPho:" << nPho << " pT:" << iPho->pt() << " eta:" << iPho->eta() << " E:" << iPho->energy() << std::endl;

  } // recoPhotons
  */

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

//____ Fill EB rechits _____//
void RecHitAnalyzer::fillEBrechits ( edm::Handle<EcalRecHitCollection> &EBRecHitsH,
                                     const CaloGeometry* &caloGeom,
                                     TH2D* &hEvt_HBHE_EMenergy )
{

  int iphi_, ieta_, idx;
  float eta, phi;

  // Initialize arrays
  vEB_energy_.assign(EBDetId::kSizeForDenseIndexing,0.);
  vEB_time_.assign(EBDetId::kSizeForDenseIndexing,0.);
  //vEB_energy_.assign(EcalTrigTowerDetId::kEBTotalTowers,0.);
  //vEB_time_.assign(EcalTrigTowerDetId::kEBTotalTowers,0);

  // Record signal-full entries
  for(EcalRecHitCollection::const_iterator iRHit = EBRecHitsH->begin();
      iRHit != EBRecHitsH->end();                      
      ++iRHit) {

    // Get detector id and convert to histogram-friendly coordinates
    EBDetId ebId( iRHit->id() );
    //EBDetId ebId( 1, 3 );
    //EcalTrigTowerDetId ttId( iRHit->id() );
    iphi_ = ebId.iphi()-1;
    ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta();
    //std::cout << "ECAL | (ieta,iphi): (" << ebId.ieta() << "," << ebId.iphi() << ")" <<std::endl;

    // Fill some histograms to monitor distributions
    // These will contain *cumulative* statistics and as such
    // should be used for monitoring purposes only
    hEB_energy->Fill( iphi_,ieta_,iRHit->energy() );
    hEB_time->Fill( iphi_,ieta_,iRHit->time() );

    // Get Hashed Index
    // Hashed index provides a convenient index mapping
    // from [ieta][iphi] -> [idx]
    idx   = ebId.hashedIndex(); // (ieta_+EBDetId::MAX_IETA)*EBDetId::MAX_IPHI + iphi_
    //idx = ttId.hashedIndex();

    // Get global position of cell center
    pos  = caloGeom->getPosition(ebId);
    eta = pos.eta();
    phi = pos.phi();
    //std::cout << "ECAL | (eta,phi,E): (" << eta << "," << phi << ","<< iRHit->energy()<<")" <<std::endl;

    // Fill event arrays
    // These are the actual inputs to the detector images
    vEB_energy_[idx] += iRHit->energy();
    vEB_time_[idx] += iRHit->time();
    //vEB_energy_[idx] += iRHit->energy()/TMath::CosH(eta); // pick out only transverse component

    hEvt_HBHE_EMenergy->Fill( phi, eta, iRHit->energy() );

  } // EB rechits

  /*
  // FOR GEOMETRY DEBUGGING
  float phi_edge;
  for(int i = 1; i < 360+1; i++) {

    // Get detector id and convert to histogram-friendly coordinates
    //EBDetId ebId( iRHit->id() );
    EBDetId ebId( 85, i );
    //EcalTrigTowerDetId ttId( iRHit->id() );
    iphi_ = ebId.iphi()-1;
    ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta();
    //std::cout << "ECAL | (ieta,iphi): (" << ebId.ieta() << "," << ebId.iphi() << ")" <<std::endl;

    // Get global position of cell center
    pos  = caloGeom->getPosition(ebId);
    eta = pos.eta();
    phi = pos.phi();

    cell = caloGeom->getGeometry(ebId);
    phi_edge = cell->getCorners()[3].phi();

    //std::cout << "ECAL | (eta,phi,E): (" << eta << "," << phi << ","<< iRHit->energy()<<")" <<std::endl;
    std::cout << iphi_+1 << ":" << phi << ":" << phi_edge << std::endl;
    //std::cout << iphi_+1 << ":" << phi << ":" << cell->getCorners()[0].phi() << "," << cell->getCorners()[1].phi() 
    //    << phi_edge << "," << cell->getCorners()[3].phi() << std::endl;
    //std::cout << phi_edge << ","<<std::endl;
  }
  */

} // fillEBrechits()

//____ Fill EE rechits _____//
void RecHitAnalyzer::fillEErechits ( edm::Handle<EcalRecHitCollection> &EERecHitsH,
                                     const CaloGeometry* &caloGeom,
                                     TH2D* &hEvt_HBHE_EMenergy,
                                     TH2D* (&hEvt_EE_energy)[EE_IZ_MAX] )
{

  int idx;
  int ix_, iy_, iz_; // NOTE: rows:iy, columns:ix
  float eta, phi;

  // Initialize arrays
  for(int iz(0); iz < EE_IZ_MAX; ++iz) {
    vEE_energy_[iz].assign(EE_NC_PER_ZSIDE,0.);
    vEE_time_[iz].assign(EE_NC_PER_ZSIDE,0.);
  }

  // Record signal-full entries
  for(EcalRecHitCollection::const_iterator iRHit = EERecHitsH->begin();
      iRHit != EERecHitsH->end();                      
      ++iRHit) {

    // Get detector id and convert to histogram-friendly coordinates
    EEDetId eeId( iRHit->id() );
    ix_ = eeId.ix()-1;
    iy_ = eeId.iy()-1;
    iz_ = (eeId.zside() > 0) ? 1 : 0;
    //std::cout << "ECAL | (ix,iy): " << ix_ << "," << iy_ << "," << iz_ << "," << iRHit->energy() << std::endl;

    // Fill some histograms to monitor distributions
    // These will contain *cumulative* statistics and as such
    // should be used for monitoring purposes only
    hEE_energy[iz_]->Fill( ix_,iy_,iRHit->energy() );
    hEE_time[iz_]->Fill( ix_,iy_,iRHit->time() );

    // Create hashed Index
    // Maps from [iy][ix] -> [idx]
    idx = iy_*EEDetId::IX_MAX + ix_; 

    // Fill event arrays
    // These are the actual inputs to the detector images
    vEE_energy_[iz_][idx] += iRHit->energy();
    vEE_time_[iz_][idx] += iRHit->time();
    //vEE_energy_[iz_][idx] += iRHit->energy()/TMath::CosH(cell->etaPos()); // pick out only transverse component

    // Get global position of cell center
    pos  = caloGeom->getPosition(eeId);
    eta = pos.eta();
    phi = pos.phi();

    hEvt_HBHE_EMenergy->Fill( phi, eta, iRHit->energy() );

    // Create fake ieta, iphi geometry for EEs and fill with EE rechits
    hEvt_EE_energy[iz_]->Fill( phi, eta, iRHit->energy() );

  } // EE reduced rechits

} // fillEErechits()

//____ Fill HBHE rechits _____//
void RecHitAnalyzer::fillHBHErechits ( edm::Handle<HBHERecHitCollection> &HBHERecHitsH, const CaloGeometry* &caloGeom ) {

  int iphi_, ieta_, ietaAbs_, idx;
  //int depth_;
  //float eta, phi;

  // Initialize arrays
  vHBHE_energy_EB_.assign( hcaldqm::constants::IPHI_NUM*2*HBHE_IETA_MAX_EB,0. );
  vHBHE_energy_.assign( hcaldqm::constants::IPHI_NUM*2*(hcaldqm::constants::IETA_MAX_HE-1),0. );
  TH2D *hEvt_HBHE_energy = new TH2D("evt_HBHE_energy", "E(i#phi,i#eta);i#phi;i#eta",
      hcaldqm::constants::IPHI_NUM,           hcaldqm::constants::IPHI_MIN-1,    hcaldqm::constants::IPHI_MAX,
      2*(hcaldqm::constants::IETA_MAX_HE-1),-(hcaldqm::constants::IETA_MAX_HE-1),hcaldqm::constants::IETA_MAX_HE-1 );

  // Record signal-full entries
  for(HBHERecHitCollection::const_iterator iRHit = HBHERecHitsH->begin();
      iRHit != HBHERecHitsH->end();                      
      ++iRHit) {

    // Get detector id and convert to histogram-friendly coordinates
    // NOTE: HBHE detector ids are indexed by (ieta,iphi,depth)!
    HcalDetId hId( iRHit->id() );
    //HcalDetId hId( HcalSubdetector::HcalBarrel, 1, 71, 1 );
    //if (hId.subdet() != HcalSubdetector::HcalBarrel) continue;
    // WARNING: HBHE::iphi() is not aligned with EBRecHit::iphi()!
    // => Need to shift by 2 HBHE towers: HBHE::iphi: [1,...,71,72]->[3,4,...,71,72,1,2]
    iphi_  = hId.iphi()+2; // shift
    iphi_  = iphi_ > hcaldqm::constants::IPHI_MAX ? iphi_-hcaldqm::constants::IPHI_MAX : iphi_; // wrap-around
    iphi_  = iphi_ - 1; // make histogram-friendly
    ietaAbs_  = hId.ietaAbs() == hcaldqm::constants::IETA_MAX_HE ? hcaldqm::constants::IETA_MAX_HE-1 : hId.ietaAbs();
    ieta_  = hId.zside() > 0 ? ietaAbs_-1 : -ietaAbs_;
    //depth_ = hId.depth();
    //std::cout << "HCAL | (ieta,iphi): (" << hId.ieta() << "," << iphi_+1 << ")" <<std::endl;

    // Fill some histograms to monitor distributions
    // These will contain *cumulative* statistics and as such
    // should be used for monitoring purposes only
    //hHBHE_depth->Fill( depth_ );

    // Get global position of cell center
    //pos = caloGeom->getPosition(hId);
    //eta = pos.eta();
    //phi = pos.phi();
    //std::cout << "HCAL > (eta,phi,E): (" << eta << "," << phi << ","<< iRHit->energy()<<")" <<std::endl;

    // Fill histos/arrays by ieta coverage of HBHE
    // hId.ieta() <= 17: match coverage of EB
    // hId.ieta() <= 20: match until iphi granularity decreases
    // hId.ieta()  > 20: match till end of HE

    // Full HBHE coverage: hId.ieta()  > 20
    // After iphi granularity drop, fill adjacent iphi and split energy evenly
    if ( hId.ietaAbs() > HBHE_IETA_MAX_iEta20 ) {
      hHBHE_energy->Fill( iphi_  ,ieta_,iRHit->energy()*0.5 );
      hHBHE_energy->Fill( iphi_+1,ieta_,iRHit->energy()*0.5 );
      hEvt_HBHE_energy->Fill( iphi_  ,ieta_,iRHit->energy()*0.5 );
      hEvt_HBHE_energy->Fill( iphi_+1,ieta_,iRHit->energy()*0.5 );
      continue;
    } else {
      hHBHE_energy->Fill( iphi_,ieta_,iRHit->energy() );
      hEvt_HBHE_energy->Fill( iphi_,ieta_,iRHit->energy() );
    }

    // Fill HBHE coverage with EB overlap: hId.ieta() <= 17
    if ( hId.ietaAbs() > HBHE_IETA_MAX_EB ) continue;
    hHBHE_energy_EB->Fill( iphi_,ieta_,iRHit->energy() );

    // Create hashed Index
    // Effectively sums energies over depth for a given (ieta,iphi)
    // Maps from [ieta][iphi] -> [idx]
    idx = ( ieta_+HBHE_IETA_MAX_EB )*hcaldqm::constants::IPHI_NUM + iphi_;

    // Fill event arrays
    // These are the actual inputs to the detector images
    vHBHE_energy_EB_[idx] += iRHit->energy();
    //vHBHE_energy_EB_[idx] += iRHit->energy()/TMath::CosH(eta); // pick out only transverse component
  }

  // Fill HBHE energy branch
  // This stores the actual HBHE image
  for (int ieta = 1; ieta < hEvt_HBHE_energy->GetNbinsY()+1; ieta++) {
    for (int iphi = 1; iphi < hEvt_HBHE_energy->GetNbinsX()+1; iphi++) {
      idx = (ieta-1)*hcaldqm::constants::IPHI_NUM + (iphi-1); 
      vHBHE_energy_[idx] += hEvt_HBHE_energy->GetBinContent( iphi, ieta );
    } 
  }  

  /*
  // FOR GEOMETRY DEBUGGING
  float phi;
  for(int i=1; i<72+1;i++) {
  //for(int i=1; i<16+1;i++) {
  //for(int i=16; i<29+1;i++) {

    // Get detector id and convert to histogram-friendly coordinates
    // NOTE: HBHE detector ids are indexed by (ieta,iphi,depth)!
    //HcalDetId hId( iRHit->id() );
    HcalDetId hId( HcalSubdetector::HcalBarrel, 1, i, 1 );
    //HcalDetId hId( HcalSubdetector::HcalBarrel, i, 1, 1 );
    //HcalDetId hId( HcalSubdetector::HcalEndcap, i, 1, 1 );
    //if (hId.subdet() != HcalSubdetector::HcalBarrel) continue;
    // WARNING: HBHE::iphi() is not aligned with EBRecHit::iphi()!
    // => Need to shift by 2 HBHE towers: HBHE::iphi: [1,...,71,72]->[3,4,...,71,72,1,2]
    iphi_  = hId.iphi()+2; // shift
    iphi_  = iphi_ > 72 ? iphi_-72 : iphi_; // wrap-around
    iphi_  = iphi_ -1; // make histogram-friendly
    //iphi_  = hId.iphi()-1;
    ieta_  = hId.ieta() > 0 ? hId.ieta()-1 : hId.ieta();
    //std::cout << "HCAL | (ieta,iphi): (" << hId.ieta() << "," << iphi_+1 << ")" <<std::endl;

    // Get global position of cell center
    pos = caloGeom->getPosition(hId);
    //eta = pos.eta();
    phi = pos.phi();
    //std::cout << "HCAL > (eta,phi): (" << eta << "," << phi << ","<< std::endl;
    std::cout <<  hId.iphi() << ":" << iphi_+1 << ":" << phi << std::endl;
    //std::cout <<  hId.ieta() << ":" << eta << std::endl;
  }
  */

} // fillHBHErechits()

//____ Fill ECAL @HCAL granularity _____//
void RecHitAnalyzer::fillECALatHCAL( TH2D* &hEvt_HBHE_EMenergy ) {

  float ieta, iphi, idx;

  // Convert EM energy HCAL histogram to vector for conversion to image later
  vHBHE_EMenergy_.assign( hcaldqm::constants::IPHI_NUM*2*(hcaldqm::constants::IETA_MAX_HE-1),0. );
  for (int ieta_ = 1; ieta_ < hEvt_HBHE_EMenergy->GetNbinsY()+1; ieta_++) {
    ieta = ieta_ - 1;
    //std::cout << ieta+1 << ":" << hEvt_HBHE_EMenergy->GetYaxis()->GetBinCenter( ieta+1 ) << std::endl;
    for (int iphi_ = 1; iphi_ < hEvt_HBHE_EMenergy->GetNbinsX()+1; iphi_++) {
      // WARNING: EB detector iphi=1 does not correspond to physical phi=-pi so need to shift
      iphi = iphi_ + 38; // shift
      iphi = iphi > hcaldqm::constants::IPHI_MAX ? iphi-hcaldqm::constants::IPHI_MAX : iphi; // wrap-around
      iphi = iphi - 1;
      idx = ieta*hcaldqm::constants::IPHI_NUM + iphi; 
      vHBHE_EMenergy_[idx] += hEvt_HBHE_EMenergy->GetBinContent( iphi_, ieta_ );
      //std::cout << iphi_ << ":" << iphi+1 << ":" << hEvt_HBHE_EMenergy->GetXaxis()->GetBinCenter( iphi_ ) << std::endl;
    } 
    //break;
  }  

} // fillECALatHCAL()

//____ Fill stitched EE-, EB, EE+ geometry _____//
void RecHitAnalyzer::fillECALstitched ( TH2D* (&hEvt_EE_energy)[EE_IZ_MAX] ) {

  int iphi, ieta, idx;

  vECAL_energy_.assign( 2*(EBDetId::MAX_IETA+55)*EBDetId::MAX_IPHI,0. );

  // Fill fake EE-(ieta,iphi) geometry with EE- hits
  int offset = 0;
  for (int ieta_ = 1; ieta_ < hEvt_EE_energy[0]->GetNbinsY()+1; ieta_++) {
    ieta = ieta_ - 1;
    //std::cout << ieta+1 << ":" << hEvt_EE_energy[0]->GetYaxis()->GetBinCenter( ieta+1 ) << std::endl;
    for (int iphi_ = 1; iphi_ < hEvt_EE_energy[0]->GetNbinsX()+1; iphi_++) {
      // WARNING: HBHE::iphi() is not aligned with EBRecHit::iphi()!
      // => Need to shift by 2 HBHE towers: HBHE::iphi: [1,...,71,72]->[3,4,...,71,72,1,2]
      iphi = iphi_ + 190; // shift
      iphi = iphi > EBDetId::MAX_IPHI ? iphi-EBDetId::MAX_IPHI : iphi; // wrap-around
      iphi = iphi - 1;
      idx = ieta*EBDetId::MAX_IPHI + iphi; 
      vECAL_energy_[idx] += hEvt_EE_energy[0]->GetBinContent( iphi_, ieta_ );
      //if (ieta_ == 1) std::cout << iphi_ << ":" << iphi+1 << ":" << hEvt_EE_energy[0]->GetXaxis()->GetBinCenter( iphi_ ) << std::endl;
      offset++;
    } 
  }  

  // Append true EB hits
  float offset_ = offset;
  for (int idx_ = 0; idx_ < 2*EBDetId::MAX_IETA*EBDetId::MAX_IPHI; idx_++) {
    idx = idx_ + offset_;
    vECAL_energy_[idx] += vEB_energy_[idx_];
    offset++;
  }
  offset_ = offset;

  // Append fake EE+(ieta,iphi) geometry with EE+ hits
  for (int ieta_ = 1; ieta_ < hEvt_EE_energy[1]->GetNbinsY()+1; ieta_++) {
    ieta = ieta_ - 1;
    //std::cout << ieta+1 << ":" << hEvt_EE_energy[1]->GetYaxis()->GetBinCenter( ieta+1 ) << std::endl;
    for (int iphi_ = 1; iphi_ < hEvt_EE_energy[1]->GetNbinsX()+1; iphi_++) {
      // WARNING: HBHE::iphi() is not aligned with EBRecHit::iphi()!
      // => Need to shift by 2 HBHE towers: HBHE::iphi: [1,...,71,72]->[3,4,...,71,72,1,2]
      iphi = iphi_ + 190; // shift
      iphi = iphi > EBDetId::MAX_IPHI ? iphi-EBDetId::MAX_IPHI : iphi; // wrap-around
      iphi = iphi - 1;
      idx = ieta*EBDetId::MAX_IPHI + iphi + offset_; 
      vECAL_energy_[idx] += hEvt_EE_energy[1]->GetBinContent( iphi_, ieta_ );
      //std::cout << iphi_ << ":" << iphi+1 << ":" << hEvt_EE_energy[1]->GetXaxis()->GetBinCenter( iphi_ ) << std::endl;
      offset++;
    } 
    //break;
  }  

  if (offset != 280*EBDetId::MAX_IPHI) std::cout << "!!! full extended ECAL not traversed !!!: " << offset << std::endl;

} // fillECALstitched()

//____ Fill EB digis _____//
void RecHitAnalyzer::fillEBdigis ( edm::Handle<EBDigiCollection> &EBDigisH, const CaloGeometry* &caloGeom ) {

  int iphi_, ieta_, idx;

  // Initialize arrays
  for(int iS(0); iS < EcalDataFrame::MAXSAMPLES; ++iS)
    vEB_adc_[iS].assign(EBDetId::kSizeForDenseIndexing,0);

  // Record signal-full entries
  for(EBDigiCollection::const_iterator iDigi = EBDigisH->begin();
      iDigi != EBDigisH->end();                      
      ++iDigi) {

    // Get detector id and convert to histogram-friendly coordinates
    EBDetId ebId( iDigi->id() );
    //DetId id( iDigi->id() );
    iphi_ = ebId.iphi()-1;
    ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta();

    // Get Hashed Index & Cell Geometry
    // Hashed index provides a convenient index mapping
    // from [ieta][iphi] -> [idx]
    idx = ebId.hashedIndex(); // (ieta_+EBDetId::MAX_IETA)*EBDetId::MAX_IPHI + iphi_
    // Cell geometry provides access to (rho,eta,phi) coordinates of cell center
    //cell  = caloGeom->getGeometry(ebId);

    // Unpack the digi into a dataframe
    EcalDataFrame df(*iDigi);
    for(int iS(0); iS < EcalDataFrame::MAXSAMPLES; ++iS) {

      // Get the iS-th sample
      EcalMGPASample digiSample( df.sample(iS) );

      // Fill some histograms to monitor distributions
      // These will contain *cumulative* statistics and as such
      // should be used for monitoring purposes only
      hEB_adc[iS]->Fill( iphi_, ieta_, digiSample.adc() );

      // Fill event arrays
      // These are the actual inputs to the detector images
      vEB_adc_[iS][idx] += digiSample.adc();
      //vEB_adc_[iS][idx] += digiSample.adc()/TMath::CosH(cell->etaPos()); // pick out only transverse component

    } // sample

  } // EB digi

} // fillEBdigis()


//define this as a plug-in
DEFINE_FWK_MODULE(RecHitAnalyzer);
