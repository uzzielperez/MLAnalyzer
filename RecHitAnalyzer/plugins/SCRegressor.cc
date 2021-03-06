// -*- C++ -*-
//
// Package:    MLAnalyzer/SCRegressor
// Class:      SCRegressor
// 
//
// Original Author:  Michael Andrews
//         Created:  Mon, 17 Jul 2017 15:59:54 GMT
//
//

#include "MLAnalyzer/RecHitAnalyzer/interface/SCRegressor.h"
//#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"

//
// constructors and destructor
//
SCRegressor::SCRegressor(const edm::ParameterSet& iConfig)
{
  //EBRecHitCollectionT_ = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EBRecHitCollection"));
  //electronCollectionT_ = consumes<edm::View<reco::GsfElectron>>(iConfig.getParameter<edm::InputTag>("gsfElectronCollection"));
  electronCollectionT_ = consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("gsfElectronCollection"));
  photonCollectionT_ = consumes<PhotonCollection>(iConfig.getParameter<edm::InputTag>("photonCollection"));
  EBRecHitCollectionT_ = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedEBRecHitCollection"));
  EERecHitCollectionT_    = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedEERecHitCollection"));
  ESRecHitCollectionT_    = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedESRecHitCollection"));
  genParticleCollectionT_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticleCollection"));
  genJetCollectionT_ = consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJetCollection"));
  trackCollectionT_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("trackCollection"));
  rhoLabel_ = consumes<double>(iConfig.getParameter<edm::InputTag>("rhoLabel"));
  trgResultsT_ = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("trgResults"));

  //now do what ever initialization is needed
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  // Output Tree
  RHTree = fs->make<TTree>("RHTree", "RecHit tree");
  RHTree->Branch("eventId", &eventId_);
  RHTree->Branch("runId",   &runId_);
  RHTree->Branch("lumiId",  &lumiId_);

  RHTree->Branch("SC_iphi", &vIphi_Emax_);
  RHTree->Branch("SC_ieta", &vIeta_Emax_);

  //branchesPiSel ( RHTree, fs );
  //branchesPhotonSel ( RHTree, fs );
  branchesDiPhotonSel ( RHTree, fs );
  branchesH2aaSel ( RHTree, fs );
  branchesSC     ( RHTree, fs );
  //branchesEB     ( RHTree, fs );
  //branchesTracksAtEBEE     ( RHTree, fs );
  branchesPhoVars     ( RHTree, fs );

}

SCRegressor::~SCRegressor()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//
//
// ------------ method called for each event  ------------
void
SCRegressor::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  edm::Handle<EcalRecHitCollection> EBRecHitsH;
  iEvent.getByToken(EBRecHitCollectionT_, EBRecHitsH);

  edm::Handle<PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);

  // Provides access to global cell position and coordinates below
  edm::ESHandle<CaloGeometry> caloGeomH;
  iSetup.get<CaloGeometryRecord>().get(caloGeomH);
  const CaloGeometry* caloGeom = caloGeomH.product();

  // Run explicit jet selection
  bool hasPassed;
  vPreselPhoIdxs_.clear();
  nTotal += nPhotons;
  //hasPassed = runPiSel ( iEvent, iSetup ); //TODO: add config-level switch
  //hasPassed = runPhotonSel ( iEvent, iSetup );
  hasPassed = runDiPhotonSel ( iEvent, iSetup );
  //hasPassed = runH2aaSel ( iEvent, iSetup );
  if ( !hasPassed ) return;
  bool runGen = runH2aaSel ( iEvent, iSetup );

  // Get coordinates of photon supercluster seed
  nPho = 0;
  int iphi_Emax, ieta_Emax;
  float Emax;
  GlobalPoint pos_Emax;
  std::vector<GlobalPoint> vPos_Emax;
  vIphi_Emax_.clear();
  vIeta_Emax_.clear();
  vRegressPhoIdxs_.clear();
  int iphi_, ieta_; // rows:ieta, cols:iphi
  for ( unsigned int i = 0; i < vPreselPhoIdxs_.size(); i++ ) {

    ///*
    PhotonRef iPho( photons, vPreselPhoIdxs_[i] );

    // Get underlying super cluster
    reco::SuperClusterRef const& iSC = iPho->superCluster();
    //EcalRecHitCollection::const_iterator iRHit_( EBRecHitsH->find(iSC->seed()->seed()) );
    //std::cout << "Seed E: " << iRHit_->energy() << std::endl;
    std::vector<std::pair<DetId, float>> const& SCHits( iSC->hitsAndFractions() );
    //std::cout << " >> SChits.size: " << SCHits.size() << std::endl;

    // Get Emax crystal
    Emax = 0.;
    iphi_Emax = -1;
    ieta_Emax = -1;

    // Loop over SC hits of photon
    for(unsigned iH(0); iH != SCHits.size(); ++iH) {

      // Get DetId
      if ( SCHits[iH].first.subdetId() != EcalBarrel ) continue;
      EcalRecHitCollection::const_iterator iRHit( EBRecHitsH->find(SCHits[iH].first) );
      if ( iRHit == EBRecHitsH->end() ) continue;

      // Convert coordinates to ordinals
      EBDetId ebId( iRHit->id() );
      //EBDetId ebId( iSC->seed()->seed() );
      ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta(); // [-85,...,-1,1,...,85]
      ieta_ += EBDetId::MAX_IETA; // [0,...,169]
      iphi_ = ebId.iphi()-1; // [0,...,359]

      // Keep coordinates of shower max
      if ( iRHit->energy() > Emax ) {
        Emax = iRHit->energy();
        iphi_Emax = iphi_;
        ieta_Emax = ieta_;
        pos_Emax = caloGeom->getPosition(ebId);
      }
      //std::cout << " >> " << iH << ": iphi_,ieta_,E: " << iphi_ << ", " << ieta_ << ", " << iRHit->energy() << std::endl;
    } // SC hits

    // Apply selection on position of shower seed
    //std::cout << " >> Found: iphi_Emax,ieta_Emax: " << iphi_Emax << ", " << ieta_Emax << std::endl;
    if ( Emax <= zs ) continue;
    if ( ieta_Emax > 169 - 15 || ieta_Emax < 15 ) continue;
    vIphi_Emax_.push_back( iphi_Emax );
    vIeta_Emax_.push_back( ieta_Emax );
    vPos_Emax.push_back( pos_Emax );
    vRegressPhoIdxs_.push_back( vPreselPhoIdxs_[i] );
    //std::cout << " >> Found: iphi_Emax,ieta_Emax: " << iphi_Emax << ", " << ieta_Emax << std::endl;
    //*/
    nPho++;

  } // Photons

  // Enforce selection
  if ( debug ) std::cout << " >> nPho: " << nPho << std::endl;
  //if ( nPho == 0 ) return;
  if ( nPho != 2 ) return;
  if ( debug ) std::cout << " >> Passed cropping. " << std::endl;

  //fillPiSel ( iEvent, iSetup );
  //fillPhotonSel ( iEvent, iSetup );
  fillDiPhotonSel ( iEvent, iSetup );
  fillH2aaSel ( iEvent, iSetup );
  fillSC     ( iEvent, iSetup );
  //fillEB     ( iEvent, iSetup );
  //fillTracksAtEBEE     ( iEvent, iSetup );
  fillPhoVars     ( iEvent, iSetup );

  eventId_ = iEvent.id().event();
  runId_ = iEvent.id().run();
  lumiId_ = iEvent.id().luminosityBlock();
  //nPassed++;
  nPassed += nPho;

  RHTree->Fill();

#ifdef THIS_IS_AN_EVENT_EXAMPLE
  Handle<ExampleData> pIn;
  iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void
SCRegressor::beginJob()
{
  nTotal = 0;
  nPassed = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void
SCRegressor::endJob()
{
  std::cout << " selected: " << nPassed << "/" << nTotal << std::endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SCRegressor::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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

//define this as a plug-in
DEFINE_FWK_MODULE(SCRegressor);
