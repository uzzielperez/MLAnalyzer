#ifndef SCRegressor_h
#define SCRegressor_h
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

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h" // reco::PhotonCollection defined here

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TStyle.h"
#include "TMath.h"
#include "TProfile2D.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class SCRegressor : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
    explicit SCRegressor(const edm::ParameterSet&);
    ~SCRegressor();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    // ----------member data ---------------------------
    //edm::EDGetTokenT<edm::View<reco::GsfElectron>> electronCollectionT_;
    edm::EDGetTokenT<reco::GsfElectronCollection> electronCollectionT_;
    edm::EDGetTokenT<reco::PhotonCollection> photonCollectionT_;
    edm::EDGetTokenT<EcalRecHitCollection> EBRecHitCollectionT_;
    edm::EDGetTokenT<reco::GenParticleCollection> genParticleCollectionT_;
    edm::EDGetTokenT<reco::GenJetCollection> genJetCollectionT_;


    static const int nPhotons = 2;
    //static const int nPhotons = 1;

    TProfile2D *hEB_energy;
    TProfile2D *hEB_time;
    std::vector<float> vEB_energy_;
    std::vector<float> vEB_time_;

    //TH1D * histo;
    TH2D * hSC_energy[nPhotons];
    TH2D * hSC_time[nPhotons];
    TH1F * hSC_mass;
    TH1F * hDR;
    TH1F * hdEta;
    TH1F * hdPhi;
    TH3F * hdPhidEta;
    TH1F * hPt;

    TTree* RHTree;

    float eventId_;

    void branchesSC ( TTree*, edm::Service<TFileService>& );
    void branchesEB ( TTree*, edm::Service<TFileService>& );

    void fillSC     ( const edm::Event&, const edm::EventSetup& );
    void fillEB     ( const edm::Event&, const edm::EventSetup& );

    void branchesPhoSel ( TTree*, edm::Service<TFileService>& );
    bool runPhoSel ( const edm::Event&, const edm::EventSetup& );
    void fillPhoSel     ( const edm::Event&, const edm::EventSetup& );

    std::vector<int> vGenPi0Idxs_;
    std::vector<int> vPreselPhoIdxs_;
    std::vector<int> vRegressPhoIdxs_;
    std::vector<float> vIphi_Emax;
    std::vector<float> vIeta_Emax;

    std::vector<float> vEB_SCenergy_;
    std::vector<float> vSC_energy_[nPhotons];
    std::vector<float> vSC_energyT_[nPhotons];
    std::vector<float> vSC_energyZ_[nPhotons];
    std::vector<float> vSC_time_[nPhotons];
    float vPho_pT_[nPhotons];
    float vPho_E_[nPhotons];
    float vPho_eta_[nPhotons];
    float vPho_phi_[nPhotons];
    float vPho_r9_[nPhotons];
    float vPho_sieie_[nPhotons];
    float vSC_mass_[nPhotons];
    float vSC_DR_[nPhotons];
    float vSC_pT_[nPhotons];

    int nTotal, nPassed;

    //TProfile2D * hnPho;
    TH2F * hnPho;
    TH2F * hnPhoGt2;
    TH1F * hdR_nPhoGt2;
    TH2F * hdPhidEta_nPhoGt2;
    TProfile2D * hdPhidEta_jphoPt_o_iphoPt;
    TH1F * hjphoPt_o_iphoPt;

};

//
// constants, enums and typedefs
//
static const float zs = 0.;

static const int crop_size = 32;
static const bool debug = true;
//static const bool debug = false;

static const int EB_IPHI_MIN = EBDetId::MIN_IPHI;//1;
static const int EB_IPHI_MAX = EBDetId::MAX_IPHI;//360;
static const int EB_IETA_MIN = EBDetId::MIN_IETA;//1;
static const int EB_IETA_MAX = EBDetId::MAX_IETA;//85;
//
// static data member definitions
//

#endif
