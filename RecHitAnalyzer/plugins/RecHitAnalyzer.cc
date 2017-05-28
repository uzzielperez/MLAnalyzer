// -*- C++ -*-
//
// Package:    MLAnalyzer/RecHitAnalyzer
// Class:      RecHitAnalyzer
// 
/**\class RecHitAnalyzer RecHitAnalyzer.cc MLAnalyzer/RecHitAnalyzer/plugins/RecHitAnalyzer.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Michael Andrews
//         Created:  Sat, 14 Jan 2017 17:45:54 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalDetId/interface/EcalTrigTowerDetId.h"
#include "Geometry/CaloTopology/interface/EcalBarrelTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDetId/interface/HcalElectronicsId.h"
#include "DQM/HcalCommon/interface/Constants.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class RecHitAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
    explicit RecHitAnalyzer(const edm::ParameterSet&);
    ~RecHitAnalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    // ----------member data ---------------------------
    edm::EDGetTokenT<EcalRecHitCollection> redEBRecHitCollectionT_;
    edm::EDGetTokenT<EcalRecHitCollection> EBRecHitCollectionT_;
    edm::EDGetTokenT<EBDigiCollection>     EBDigiCollectionT_;
    edm::EDGetTokenT<HBHERecHitCollection> HBHERecHitCollectionT_;
    //edm::InputTag trackTags_; //used to select what tracks to read from configuration file

    //TH1D * histo; 
    TH2D * hEBenergy; 
    TH2D * hEBenergyRed; 
    TH2D * hEBtiming; 
    TH2D * hEB_adc[EcalDataFrame::MAXSAMPLES]; 

    TH2D * hHBHEenergy; 
    TH2D * hHBHEenergy_EB; 
    TH1D * hHBHEdepth; 

    TTree* RHTree;

    std::vector<float> vEBenergy_;
    std::vector<float> vEBtiming_;
    std::vector<float> vEBenergyRed_;
    std::vector<float> vEBtimingRed_;
    std::vector<float> vEB_adc_[EcalDataFrame::MAXSAMPLES];

    //std::vector<float> vHBHEenergy[hcaldqm::constants::DEPTH_NUM];
    std::vector<float> vHBHEenergy_EB_;
  
    //const unsigned HBHE_IETA_MAX = hcaldqm::constants::IETA_MAX_HB + 1;//17
    const unsigned HBHE_IETA_MAX = 20;


    TCanvas *cEB, *cHBHE;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
RecHitAnalyzer::RecHitAnalyzer(const edm::ParameterSet& iConfig)
{
  redEBRecHitCollectionT_ = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedEBRecHitCollection"));
  EBRecHitCollectionT_ = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EBRecHitCollection"));
  //EBDigiCollectionT_ = consumes<EBDigiCollection>(iConfig.getParameter<edm::InputTag>("selectedEBDigiCollection"));
  EBDigiCollectionT_ = consumes<EBDigiCollection>(iConfig.getParameter<edm::InputTag>("EBDigiCollection"));
  HBHERecHitCollectionT_ = consumes<HBHERecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedHBHERecHitCollection"));

  //now do what ever initialization is needed
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  //histo = fs->make<TH1D>("charge" , "Charges" , 200 , -2 , 2 );

  // Canvases
  cEB = new TCanvas("cEB","cEB",600,300);
  cHBHE = new TCanvas("cHBHE","cHBHE",600,300);

  // Histograms
  // EB rechits
  hEBenergy = fs->make<TH2D>("EB_rechitE", "E(i#phi,i#eta);i#phi;i#eta",
      EBDetId::MAX_IPHI  , EBDetId::MIN_IPHI-1, EBDetId::MAX_IPHI,
      2*EBDetId::MAX_IETA,-EBDetId::MAX_IETA,   EBDetId::MAX_IETA );
  hEBenergyRed = fs->make<TH2D>("EB_rechitEred", "Ered(i#phi,i#eta);i#phi;i#eta",
      EBDetId::MAX_IPHI  , EBDetId::MIN_IPHI-1, EBDetId::MAX_IPHI,
      2*EBDetId::MAX_IETA,-EBDetId::MAX_IETA,   EBDetId::MAX_IETA );
  //hEBenergyRed = fs->make<TH2D>("EB_rechitEred", "Ered(i#phi,i#eta);i#phi;i#eta",
  //    EcalTrigTowerDetId::kEBTowersInPhi*18  , EBDetId::MIN_IPHI-1, EBDetId::MAX_IPHI,
  //    EcalTrigTowerDetId::kEBTowersInEta*2,-EBDetId::MAX_IETA,   EBDetId::MAX_IETA );
  hEBtiming = fs->make<TH2D>("EB_rechitT", "t(i#phi,i#eta);i#phi;i#eta",
      EBDetId::MAX_IPHI  , EBDetId::MIN_IPHI-1, EBDetId::MAX_IPHI,
      2*EBDetId::MAX_IETA,-EBDetId::MAX_IETA,   EBDetId::MAX_IETA );
  // EB Digis
  char hname[50], htitle[50];
  for(int iS(0); iS < EcalDataFrame::MAXSAMPLES; iS++){
    sprintf(hname, "EB_adc%d",iS);
    sprintf(htitle,"adc%d(i#phi,i#eta);i#phi;i#eta",iS);
    hEB_adc[iS] = fs->make<TH2D>(hname, htitle,
        EBDetId::MAX_IPHI  , EBDetId::MIN_IPHI-1, EBDetId::MAX_IPHI,
        2*EBDetId::MAX_IETA,-EBDetId::MAX_IETA,   EBDetId::MAX_IETA );
  }
  // HBHE
  hHBHEenergy = fs->make<TH2D>("HBHE_rechitE", "E(i#phi,i#eta);i#phi;i#eta",
      hcaldqm::constants::IPHI_NUM,      hcaldqm::constants::IPHI_MIN-1, hcaldqm::constants::IPHI_MAX,
      2*hcaldqm::constants::IETA_MAX_HE,-hcaldqm::constants::IETA_MAX_HE,hcaldqm::constants::IETA_MAX_HE );
  hHBHEenergy_EB = fs->make<TH2D>("HBHE_rechitE_EB", "E(i#phi,i#eta);i#phi;i#eta",
      hcaldqm::constants::IPHI_NUM, hcaldqm::constants::IPHI_MIN-1,hcaldqm::constants::IPHI_MAX,
      2*HBHE_IETA_MAX,             -HBHE_IETA_MAX,                 HBHE_IETA_MAX );
  hHBHEdepth = fs->make<TH1D>("HBHE_depth", "E(i#phi,i#eta);i#phi;i#eta",
      hcaldqm::constants::DEPTH_NUM, hcaldqm::constants::DEPTH_MIN, hcaldqm::constants::DEPTH_MAX+1);

  // Output Tree
  RHTree = fs->make<TTree>("RHTree", "RecHit tree");
  RHTree->Branch("EBenergy",    &vEBenergy_);
  RHTree->Branch("EBtime",      &vEBtiming_);
  RHTree->Branch("EBenergyRed", &vEBenergyRed_);
  RHTree->Branch("EBtimeRed",   &vEBtimingRed_);
  for(int iS(0); iS < EcalDataFrame::MAXSAMPLES; iS++){
    sprintf(hname, "EB_adc%d",iS);
    RHTree->Branch(hname,       &vEB_adc_[iS]);
  }
  RHTree->Branch("HBHEenergy_EB", &vHBHEenergy_EB_);
}


RecHitAnalyzer::~RecHitAnalyzer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
RecHitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  // ----- Get Calorimeter Geometry ----- //
  // Provides access to global cell position and coordinates below
  edm::ESHandle<CaloGeometry> caloGeomH;
  iSetup.get<CaloGeometryRecord>().get(caloGeomH);
  const CaloGeometry* caloGeom = caloGeomH.product();
  //const CaloSubdetectorGeometry* towerGeom = caloGeomH.product(); 
  //towerGeom->getSubdetectorGeometry(DetId::Calo, CaloTowerDetId::SubdetId);

  // Initializer for cell position
  // e.g. provides access to x, y, z coordinates of cell center
  GlobalPoint pos;

  // Initializer for cell geometry
  // e.g. provides access, to rho, eta, phi coordinates of cell center
  const CaloCellGeometry *cell;

  //////////// EB //////////

  bool saveImgs = true;
  int iphi_, ieta_, idx;

  // ----- EB rechit collection ----- //
  // This contatins the raw EB rechit collection before
  // the zero suppression and bad channel clean-up

  // Initialize arrays
  vEBenergy_.assign(EBDetId::kSizeForDenseIndexing,0.);
  vEBtiming_.assign(EBDetId::kSizeForDenseIndexing,0);

  // Record signal-full entries
  edm::Handle<EcalRecHitCollection> EBRecHitsH;
  iEvent.getByToken(EBRecHitCollectionT_, EBRecHitsH);
  for(EcalRecHitCollection::const_iterator iRHit = EBRecHitsH->begin();
      iRHit != EBRecHitsH->end();                      
      ++iRHit) {

    // Get detector id and convert to histogram-friendly coordinates
    EBDetId ebId( iRHit->id() );
    iphi_ = ebId.iphi()-1;
    ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta();

    // Fill some histograms to monitor distributions
    // These will contain *cumulative* statistics and as such
    // such be used for monitoring purposes
    hEBenergy->Fill( iphi_,ieta_,iRHit->energy() );
    //hEBenergy->Fill( iphi_,ieta_ ); // if occupancy is preferred
    hEBtiming->Fill( iphi_,ieta_,iRHit->time() );

    // Get Hashed Index & Cell Geometry
    // Hashed index provides a convenient index mapping
    // from [ieta][iphi] -> [idx]
    idx   = ebId.hashedIndex(); // (ieta_+EBDetId::MAX_IETA)*EBDetId::MAX_IPHI + iphi_
    // Cell geometry provides access to (rho,eta,phi) coordinates of cell center
    //cell  = caloGeom->getGeometry(ebId);

    // Fill event arrays
    // These are the actual inputs to the detector images
    vEBenergy_[idx] += iRHit->energy();
    //vEBenergy_[idx] = iRHit->energy()/std::abs(cell->etaPos()); // pick out only transverse component
    vEBtiming_[idx] += iRHit->time();

  } // EB rechits

  // Write out state of histogram
  // For monitoring purposes only 
  char outFile[100];
  if (saveImgs) {
    cEB->cd();
    gPad->SetLogz(1);
    gStyle->SetPalette(1);
    gStyle->SetOptStat(0);
    //hEBenergy->GetZaxis()->SetRangeUser(2.e-2, 9.e1);
    //hEBenergy->Draw("COL Z");
    sprintf(outFile,"cEBenergy_%llu.eps",iEvent.id().event());
    cEB->Print(outFile);
  }

  // ----- EB reduced rechit collection ----- //
  // This contatins the reduced EB rechit collection after
  // the zero suppression and bad channel clean-up

  // Initialize arrays
  vEBenergyRed_.assign(EBDetId::kSizeForDenseIndexing,0.);
  vEBtimingRed_.assign(EBDetId::kSizeForDenseIndexing,0);
  //vEBenergyRed_.assign(EcalTrigTowerDetId::kEBTotalTowers,0.);
  //vEBtimingRed_.assign(EcalTrigTowerDetId::kEBTotalTowers,0);

  // Record signal-full entries
  edm::Handle<EcalRecHitCollection> redEBRecHitsH;
  iEvent.getByToken(redEBRecHitCollectionT_, redEBRecHitsH);
  for(EcalRecHitCollection::const_iterator iRHit = redEBRecHitsH->begin();
      iRHit != redEBRecHitsH->end();                      
      ++iRHit) {

    // Get detector id and convert to histogram-friendly coordinates
    EBDetId ebId( iRHit->id() );
    //EBDetId ebId( 1, 3 );
    //EcalTrigTowerDetId ttId( iRHit->id() );
    iphi_ = ebId.iphi()-1;
    ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta();
    std::cout << "ECAL | (ieta,iphi): (" << ebId.ieta() << "," << ebId.iphi() << ")" <<std::endl;

    // Fill some histograms to monitor distributions
    // These will contain *cumulative* statistics and as such
    // such be used for monitoring purposes
    hEBenergyRed->Fill( iphi_,ieta_,iRHit->energy() );
    //hEBtimingRed->Fill( iphi_,ieta_,iRHit->time() );

    // Get Hashed Index & Cell Geometry
    // Hashed index provides a convenient index mapping
    // from [ieta][iphi] -> [idx]
    idx   = ebId.hashedIndex(); // (ieta_+EBDetId::MAX_IETA)*EBDetId::MAX_IPHI + iphi_
    //idx = ttId.hashedIndex();
    // Cell geometry provides access to (rho,eta,phi) coordinates of cell center
    cell  = caloGeom->getGeometry(ebId);
    float eta = cell->etaPos();
    float phi = cell->phiPos();
    std::cout << "ECAL | (eta,phi,E): (" << eta << "," << phi << ","<< iRHit->energy()<<")" <<std::endl;
    pos  = caloGeom->getPosition(ebId);
    eta = pos.eta();
    phi = pos.phi();
    std::cout << "ECAL | (eta,phi,E): (" << eta << "," << phi << ","<< iRHit->energy()<<")" <<std::endl;

    // Fill event arrays
    // These are the actual inputs to the detector images
    //vEBenergyRed_[idx] = iRHit->energy();
    vEBtimingRed_[idx] += iRHit->time();
    //vEBenergyRed_[idx] += iRHit->energy()/TMath::CosH(cell->etaPos()); // pick out only transverse component

  } // EB reduced rechits

  // Write out state of histogram
  // For monitoring purposes only 
  if (saveImgs) { 
    cEB->Clear();
    //hEBenergyRed->GetZaxis()->SetRangeUser(2.e-2, 9.e1);
    //hEBenergyRed->Draw("COL Z");
    sprintf(outFile,"cEBenergyRed_%llu.eps",iEvent.id().event());
    cEB->Print(outFile);
  }


  // ----- EB digis ----- //
  // This contains the raw EB digi collection:
  // Each digi is unpacked as a data frame containing info
  // from 10 time samples [iS] of the pulse shape
  // iS=[0-2]: Presample noise
  // iS=[3-9]: Nominal pulse shape
  // NOTE: This is the raw collection and includes
  // zero-supression and bad channel effects!

  // Initialize arrays
  for(int iS(0); iS < EcalDataFrame::MAXSAMPLES; ++iS)
    vEB_adc_[iS].assign(EBDetId::kSizeForDenseIndexing,0);

  // Record signal-full entries
  edm::Handle<EBDigiCollection> EBDigisH;
  iEvent.getByToken(EBDigiCollectionT_, EBDigisH);
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
      // such be used for monitoring purposes
      hEB_adc[iS]->Fill( iphi_, ieta_, digiSample.adc() );

      // Fill event arrays
      // These are the actual inputs to the detector images
      vEB_adc_[iS][idx] += digiSample.adc();
      //vEB_adc_[iS][idx] += digiSample.adc()/TMath::CosH(cell->etaPos()); // pick out only transverse component

    } // sample

  } // EB digi

  // Write out state of histogram
  // For monitoring purposes only 
  if (saveImgs) {
    for(int iS(0); iS < EcalDataFrame::MAXSAMPLES; ++iS) {
      cEB->Clear();
      gPad->SetLogz(1);
      //hEB_adc[iS]->GetZaxis()->SetRangeUser(190, 1.3e3);
      //hEB_adc[iS]->Draw("COL Z");
      sprintf(outFile,"cEB_adc%d_%llu.eps",iS,iEvent.id().event());
      cEB->Print(outFile);
    } // sample
  }


  //////////// HBHE //////////

  // ----- HBHE reduced rechit collection ----- //
  int depth_;

  // Initialize arrays
  vHBHEenergy_EB_.assign( hcaldqm::constants::IPHI_NUM*2*HBHE_IETA_MAX,0. );

  // Record signal-full entries
  edm::Handle<HBHERecHitCollection> HBHERecHitsH;
  iEvent.getByToken(HBHERecHitCollectionT_, HBHERecHitsH);
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
    iphi_  = iphi_ > 72 ? iphi_-72 : iphi_; // wrap-around
    iphi_  = iphi_ -1; // make histogram-friendly
    //iphi_  = hId.iphi()-1;
    ieta_  = hId.ieta() > 0 ? hId.ieta()-1 : hId.ieta();
    depth_ = hId.depth();
    //std::cout << "HCAL | (ieta,iphi): (" << hId.ieta() << "," << hId.iphi() << ")" <<std::endl;
    std::cout << "HCAL | (ieta,iphi): (" << hId.ieta() << "," << iphi_+1 << ")" <<std::endl;

    // Fill some histograms to monitor distributions
    // These will contain *cumulative* statistics and as such
    // such be used for monitoring purposes
    hHBHEdepth->Fill( depth_ );

    // Restrict coverage of HBHE
    // HBHE_IETA_MAX == 17: match coverage of EB
    // HBHE_IETA_MAX == 20: match until granularity decreases
    if ( abs(hId.ieta()) > 20 ) {
      hHBHEenergy->Fill( iphi_  ,ieta_,iRHit->energy()*0.5 );
      hHBHEenergy->Fill( iphi_+1,ieta_,iRHit->energy()*0.5 );
      continue; 
    } else {
      hHBHEenergy->Fill( iphi_,ieta_,iRHit->energy() );
    }

    // Fill restricted coverage histograms
    hHBHEenergy_EB->Fill( iphi_,ieta_,iRHit->energy() );

    // Get array index by hand:
    // Effectively sums energies over depth for a given (ieta,iphi)
    // Maps from [ieta][iphi] -> [idx]
    idx = ( ieta_+(HBHE_IETA_MAX) )*hcaldqm::constants::IPHI_NUM + iphi_;

    // Get cell geometry
    // Cell geometry provides access to (rho,eta,phi) coordinates of cell center
    cell = caloGeom->getGeometry(hId);
    //float rho = cell->rhoPos();
    float eta = cell->etaPos();
    float phi = cell->phiPos();
    //std::cout << "HCAL | (eta,phi): (" << eta << "," << phi << ")" <<std::endl;
    std::cout << "HCAL | (eta,phi,E): (" << eta << "," << phi << ","<< iRHit->energy()<<")" <<std::endl;
    // Cell position provides access to global (x,y,z) coordinates of cell center
    pos  = caloGeom->getPosition(hId);
    eta = pos.eta();
    phi = pos.phi();
    std::cout << "HCAL > (eta,phi,E): (" << eta << "," << phi << ","<< iRHit->energy()<<")" <<std::endl;

    // Fill event arrays
    // These are the actual inputs to the detector images
    vHBHEenergy_EB_[idx] += iRHit->energy();
    //vHBHEenergy_EB_[idx] += iRHit->energy()/TMath::CosH(cell->etaPos()); // pick out only transverse component

  }

  // Write out state of histogram
  // For monitoring purposes only 
  if (saveImgs) {
    cHBHE->cd();
    //std::cout << "maxEta: " << maxEta << std::endl;
    hHBHEenergy->GetZaxis()->SetRangeUser(2.e-2, 9.e1);
    hHBHEenergy->Draw("COL Z");
    sprintf(outFile,"cHBHEenergy_%llu.eps",iEvent.id().event());
    cHBHE->Print(outFile);
  }

  /*
     using reco::TrackCollection;
     Handle<TrackCollection> tracks;
     iEvent.getByLabel(trackTags_,tracks);
     for(TrackCollection::const_iterator itTrack = tracks->begin();
     itTrack != tracks->end();                      
     ++itTrack) {
     int charge = 0;
     charge = itTrack->charge();  
     histo->Fill( charge );
     }
     */
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

//define this as a plug-in
DEFINE_FWK_MODULE(RecHitAnalyzer);
